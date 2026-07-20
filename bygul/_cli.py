import os
from Bio import SeqIO
import click
from tqdm import tqdm
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict


@click.group(context_settings={"show_default": True})
@click.version_option("3.2.0")
def cli():
    pass


@cli.command(context_settings=dict(
    ignore_unknown_options=True,
    allow_extra_args=True,
))
@click.pass_context
@click.option("--genomes",
              default="NA",
              help="Comma-separated list of"
              "genome file paths in fasta format")
@click.option(
    "--reference",
    default="NA",
    help=(
        "Reference used to generate the primer file"
    )
)
@click.option(
    "--primers",
    default="NA",
    help=(
        "primers to use for amplicon sequencing in bed format"
    ),
)
@click.option(
    "--proportions",
    default="NA",
    type=str,
    help=(
        "Read proportions for each sample, e.g.(0.8,0.2) must sum to 1.0. "
        "If not provided, the program will randomly assign proportions"
    ),
)
@click.option(
    "--csv",
    default="NA",
    type=str,
    help=(
        "Path to a CSV file containing sample paths "
        "and proportions. The CSV should have two columns: "
        "'sample_path' and 'proportion'. "
        "Cannot be used together with --proportions --genomes"
    ),
)
@click.option(
    "--multifasta",
    default="NA",
    type=str,
    help=(
        "Path to a multifasta file containing all sample sequences"
        "Cannot be used together with --genomes and --proportions."
        "Must be used together with --csv"
    ),
)
@click.option(
    "--outdir",
    default="results",
    type=click.Path(exists=False),
    help="Output directory",
    show_default=True,
)
@click.option(
    "--simulator",
    default="wgsim",
    type=click.Choice(["wgsim", "mason"], case_sensitive=False),
    help="Select the simulator to use (wgsim or mason)",
)
@click.option(
    "--simulation_mode",
    default="amplicon",
    type=click.Choice(["amplicon", "metagenomics"], case_sensitive=False),
    help="Select type of simulation",
)
@click.option(
    "--maxmismatch",
    default=1,
    show_default=True,
    help="Maximum number of mismatches allowed in primer region",
)
@click.option(
    "--wgsim_insert_size", default=150,
    help="Outer distance for simulation using wgsim in amplicon"
    "simulation mode."
)
@click.option(
    "--wgsim_read_length", default=150,
    help="Read length for simulation using wgsim in amplicon"
    "simulation mode."
)
@click.option(
    "--wgsim_error_rate", default=0.0001,
    help="Error rate for simulation using wgsim in amplicon"
    "simulation mode."
)
@click.option("--readcnt", default=500, help="Number of reads per amplicon")
@click.option(
    "--redo",
    is_flag=True,
    default=False,
    help="Overwrite the output directory if it already exists.",
)
def simulate_proportions(
    ctx,
    genomes,
    proportions,
    reference,
    primers,
    wgsim_insert_size,
    wgsim_read_length,
    wgsim_error_rate,
    outdir,
    readcnt,
    maxmismatch,
    simulator,
    redo,
    simulation_mode,
    csv,
    multifasta
):
    from bygul.utils import (
        preprocess_primers,
        merge_fastq_files,
        validate_simulation_args,
        check_dir,
        process_sample_proportions,
        process_amplicon_worker,
        process_genome_worker
    )
    # validare simulation arugments
    validate_simulation_args(simulation_mode, primers,
                             reference, proportions,
                             csv, genomes, multifasta)
    if reference != "NA":
        # read the reference sequence
        reference = next(SeqIO.parse(reference, "fasta"))
    # needed to pass simulation specific flags
    extra_simulator_flags = ctx.args
    ctx = click.get_current_context()
    genome_map = defaultdict(list)
    if csv != "NA" and multifasta != "NA":
        df = pd.read_csv(csv)
        sample_names = df["sample_name"].tolist()
        proportions = df["proportion"].tolist()
        for record in SeqIO.parse(multifasta, "fasta"):
            sample = record.id.split("|")[0]
            genome_map[sample].append(record)
        # Ensure every sample in the CSV exists
        missing = set(sample_names) - set(genome_map)
        if missing:
            raise click.ClickException(
                "Samples missing from multifasta:"
                f" {', '.join(sorted(missing))}"
            )
    else:
        sample_paths = str(genomes).split(",")
        for path in sample_paths:
            for record in SeqIO.parse(path, "fasta"):
                sample = record.id.split("|")[0]
                genome_map[sample].append(record)
        sample_names = list(genome_map.keys())
    # check directory exists- if redo specified make again
    check_dir(outdir, redo, sample_names)
    # process the proportions and give warnings if necessary
    # assign proportions randomly if not provided
    proportions = process_sample_proportions(proportions,
                                             sample_names,
                                             outdir,
                                             csv)
    # read counts defined pased on proportions
    rc = int(readcnt)
    read_cnts = np.random.multinomial(rc, proportions)

    # Global outputs for merging
    output_path1 = os.path.join(os.path.abspath(outdir), "reads_1.fastq")
    output_path2 = os.path.join(os.path.abspath(outdir), "reads_2.fastq")

    if simulation_mode == "amplicon":
        print("Reading and preprocessing the primer file...")
        df_primers_template = preprocess_primers(primers, reference)

        task_args = [
            (name, genome_map[name], cnt,
             df_primers_template, maxmismatch, outdir,
             simulator, wgsim_insert_size,
             wgsim_read_length, wgsim_error_rate,
             extra_simulator_flags)
            for name, cnt in zip(sample_names, read_cnts)
        ]

        worker_func = process_amplicon_worker

    else:
        # Setup for standard/else clause simulation mapping
        task_args = [
            (name, path, cnt, outdir, simulator,
             wgsim_insert_size,
             wgsim_read_length, wgsim_error_rate,
             extra_simulator_flags)
            for name, path, cnt in zip(sample_names, sample_paths, read_cnts)
        ]
        worker_func = process_genome_worker
    # Run Parallel Pool Execution
    print(f"Spinning up parallel execution for {len(sample_names)} samples...")
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(worker_func,
                                   arg): arg[0] for arg in task_args}

        with tqdm(total=len(sample_names),
                  desc="Simulation progress...") as pbar:
            for future in as_completed(futures):
                sample_name = futures[future]
                try:
                    result = future.result()
                    if result[0] == "success":
                        _, name, r_path1, r_path2 = result
                        print(f"\nMerging reads for {name}...")
                        merge_fastq_files(r_path1, output_path1)
                        merge_fastq_files(r_path2, output_path2)
                    elif result[0] == "warning":
                        print(f"\n[{sample_name}] {result[1]}")
                except Exception as e:
                    print(f"\nError processing sample {sample_name}: {e}")
                finally:
                    pbar.update(1)

    print("Finished all samples!")


@cli.command()
@click.argument("genomes", type=str)
@click.argument(
    "primers", type=str
)
@click.argument("reference", type=str)
@click.option(
    "--maxmismatch",
    default=1,
    show_default=True,
    help="Maximum number of mismatches allowed in primer region",
)
@click.option(
    "--outdir",
    default=".",
    type=click.Path(exists=False),
    help="Output directory",
    show_default=True,
)
def check_primers(genomes, primers,
                  reference, maxmismatch,
                  outdir):
    from bygul.utils import (
        preprocess_primers,
        process_primer_check_worker,
    )
    # read the reference sequence
    reference = next(SeqIO.parse(reference, "fasta"))
    genome_map = defaultdict(list)
    for record in SeqIO.parse(genomes, "fasta"):
        sample = record.id.split("_")[0]
        genome_map[sample].append(record)
    df_primers_template = preprocess_primers(primers, reference)
    task_args = [
        (name, genome_map[name],
            df_primers_template, maxmismatch, outdir)
        for name in genome_map.keys()
    ]
    worker_func = process_primer_check_worker

    # Run Parallel Pool Execution
    print("Spinning up parallel execution for "
          f"{len(genome_map.keys())} samples...")

    dfs = []
    with ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(worker_func, arg): arg[0]
            for arg in task_args
        }

        with tqdm(total=len(genome_map.keys()),
                  desc="Simulation progress...") as pbar:
            for future in as_completed(futures):
                sample_name = futures[future]
                try:
                    status, result = future.result()

                    if status == "success":
                        print(f"Successfully processed sample {sample_name}")
                        dfs.append(result)   # result is the DataFrame

                    elif status == "warning":
                        print(f"\n[{sample_name}] {result}")

                except Exception as e:
                    print(f"\nError processing sample {sample_name}: {e}")

                finally:
                    pbar.update(1)

    # Concatenate all successful DataFrames
    if dfs:
        final_df = pd.concat(dfs, ignore_index=True)
    else:
        final_df = pd.DataFrame()
    final_df.to_csv(os.path.join(outdir, "amplicon_stats.csv"), index=False)


if __name__ == "__main__":
    cli()
