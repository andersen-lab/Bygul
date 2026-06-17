import os
from Bio import SeqIO
import click
from tqdm import tqdm
import numpy as np
import pandas as pd


@click.group(context_settings={"show_default": True})
@click.version_option("3.2.0")
def cli():
    pass


@cli.command(context_settings=dict(
    ignore_unknown_options=True,
    allow_extra_args=True,
))
@click.pass_context
@click.argument("genomes", type=str)
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
    simulation_mode
):
    from bygul.utils import (
        preprocess_primers,
        create_valid_primer_combinations,
        make_amplicon,
        write_fasta_group,
        run_simulation_on_fasta,
        run_simulation_on_fasta_single_genome,
        merge_fastq_files,
        find_closest_primer_match,
        assess_genome_quality_from_fasta,
        validate_simulation_args,
        check_dir,
        process_sample_proportions
    )
    # validare simulation arugments
    validate_simulation_args(simulation_mode, primers, reference)
    # needed to pass simulation specific flags
    extra_simulator_flags = ctx.args
    ctx = click.get_current_context()
    # split the sample names and paths into a list
    sample_names = [fp.split("/")[-1].split(".")[0]
                    for fp in str(genomes).split(",")]
    # check directory exists- if redo specified make again
    check_dir(outdir, redo, sample_names)
    sample_paths = str(genomes).split(",")
    # Print information about the quality of the provided file
    for genome in sample_paths:
        assess_genome_quality_from_fasta(genome)
    # process the proportions and give warnings if necessary
    # assign proportions randomly if not provided
    proportions = process_sample_proportions(proportions,
                                             sample_names,
                                             sample_names,
                                             outdir)
    # read counts defined pased on proportions
    read_cnts = [i * int(readcnt) for i in proportions]
    if simulation_mode == "amplicon":
        df_primers_template = preprocess_primers(primers, reference)
        print("Reading and preprocessing the primer file...")

        with tqdm(total=len(sample_names),
                  desc="Simulation progress...") as pbar:
            for name, path, cnt in zip(sample_names, sample_paths, read_cnts):
                print(f"Processing sample {name}...")
                # List to store amplicon
                # DataFrames from all contigs in this sample
                sample_amplicons_list = []
                # Iterate through EVERY contig/sequence in the FASTA file
                for genome_seq in SeqIO.parse(path, "fasta"):
                    print("  Extracting amplicons"
                          f"from contig: {genome_seq.id}...")
                    # Find matches on the current contig
                    contig_df = find_closest_primer_match(
                        df_primers_template.copy(),
                        str(genome_seq.seq),
                        maxmismatch)
                    all_amplicons = create_valid_primer_combinations(contig_df)
                    all_amplicons = all_amplicons.fillna(0)
                    # Calculate lengths
                    all_amplicons["amplicon_length"] = np.where(
                        (all_amplicons["primer_start"] != 0) &
                        (all_amplicons["primer_end"] != 0),
                        all_amplicons["primer_end"] -
                        all_amplicons["primer_start"] +
                        all_amplicons["primer_seq_y"].str.len(),
                        0,
                    )

                    # Generate sequences for this contig
                    all_amplicons["amplicon_sequence"] = all_amplicons.apply(
                        lambda row: make_amplicon(
                            row["primer_start"],
                            row["primer_end"],
                            row["primer_seq_y"],
                            genome_seq.seq,
                        ),
                        axis=1,
                    )
                    # Tag with contig ID to avoid
                    # confusion if IDs overlap between contigs
                    all_amplicons["contig_id"] = genome_seq.id
                    sample_amplicons_list.append(all_amplicons)

                # Combine results from all contigs for this sample
                if not sample_amplicons_list:
                    print(f"Warning: No sequences found in {path}")
                    continue
                full_sample_df = pd.concat(sample_amplicons_list,
                                           ignore_index=True)
                # Output directory setup
                amp_out_dir = os.path.join(outdir, name, "amplicons")
                os.makedirs(amp_out_dir, exist_ok=True)
                full_sample_df.to_csv(os.path.join(amp_out_dir,
                                                   "amplicon_stats.csv"),
                                      index=False)
                # Write FASTA files for simulation
                full_sample_df["amplicon_suffix"] =\
                    full_sample_df["amplicon_number"].apply(
                    lambda x: x.split("_")[0] if "_" in x else x
                )
                for n, g in full_sample_df.groupby("amplicon_suffix"):
                    write_fasta_group(g, n, amp_out_dir)

                print("Starting read simulation...")
                read_dir = os.path.join(outdir, name, "reads")
                os.makedirs(read_dir, exist_ok=True)

                fasta_files = [
                    os.path.join(amp_out_dir, f)
                    for f in os.listdir(amp_out_dir)
                    if f.endswith((".fasta", ".fa"))
                ]

                for fasta_file in fasta_files:
                    run_simulation_on_fasta(
                        fasta_file,
                        read_dir,
                        cnt,
                        simulator,
                        wgsim_insert_size,
                        wgsim_read_length,
                        wgsim_error_rate,
                        extra_flags=extra_simulator_flags
                    )
                # Merging logic remains the same
                read_path1 = os.path.join(os.path.abspath(outdir),
                                          name, "reads/merged_reads_1.fastq")
                read_path2 = os.path.join(os.path.abspath(outdir),
                                          name, "reads/merged_reads_2.fastq")
                output_path1 = os.path.join(os.path.abspath(outdir),
                                            "reads_1.fastq")
                output_path2 = os.path.join(os.path.abspath(outdir),
                                            "reads_2.fastq")
                print(f"Merging reads for {name}...")
                merge_fastq_files(read_path1, output_path1)
                merge_fastq_files(read_path2, output_path2)
                pbar.update(1)
            print("Finished all samples!")
    else:
        with tqdm(total=len(sample_names),
                  desc="Simulation progress...") as pbar:
            for name, path, cnt in zip(sample_names, sample_paths, read_cnts):
                if not os.path.exists(os.path.join(outdir, name, "reads")):
                    os.makedirs(os.path.join(outdir, name, "reads"))
                run_simulation_on_fasta_single_genome(
                        path,
                        os.path.join(outdir, name, "reads"),
                        cnt,
                        simulator,
                        extra_flags=extra_simulator_flags
                    )
                read_path1 = os.path.join(
                    os.path.abspath(outdir), name, "reads/reads_1.fastq"
                )
                read_path2 = os.path.join(
                    os.path.abspath(outdir), name, "reads/reads_2.fastq"
                )
                output_path1 = os.path.join(
                    os.path.abspath(outdir), "reads_1.fastq")
                output_path2 = os.path.join(
                    os.path.abspath(outdir), "reads_2.fastq")
                print("Merging all reads...")
                merge_fastq_files(read_path1, output_path1)
                merge_fastq_files(read_path2, output_path2)
                print("Finished!")
            pbar.update(1)


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
def check_primers(genomes, primers, reference, maxmismatch):
    from bygul.utils import (
        assess_genome_quality_from_fasta,
        preprocess_primers,
        find_closest_primer_match,
        create_valid_primer_combinations,
    )

    assess_genome_quality_from_fasta(genomes)
    primer_df = preprocess_primers(primers, reference)
    print("Reading and preprocessing the primer file...")
    all_results = []

    # Parse ALL sequences in the multifasta
    for genome_record in SeqIO.parse(genomes, "fasta"):
        genome_id = genome_record.id
        genome_seq = str(genome_record.seq)
        df = find_closest_primer_match(primer_df, genome_seq,
                                       maxmismatch)
        all_amplicons = create_valid_primer_combinations(df)
        all_amplicons = all_amplicons.fillna(0)
        all_amplicons["amplicon_length"] = np.where(
                    (all_amplicons["primer_start"] != 0)
                    & (all_amplicons["primer_end"] != 0),
                    all_amplicons["primer_end"]
                    - all_amplicons["primer_start"]
                    + all_amplicons["primer_seq_y"].str.len(),
                    0,
                )
        all_amplicons["genome_id"] = genome_id

        all_results.append(all_amplicons)
        # Combine results from all fasta entries
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
    else:
        final_df = pd.DataFrame()
    final_df.to_csv(
        os.path.join("amplicon_stats.csv"),
        index=False,
    )


if __name__ == "__main__":
    cli()
