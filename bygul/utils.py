import math
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
import pandas as pd
import os
import itertools
import numpy as np
import regex as re
import warnings
import sys
import shutil


def process_sample_proportions(
    proportions,
    sample_names,
    sample_paths,
    outdir,
    csv,
):
    """
    Processes sample proportions for read simulation.

    Parameters:
        proportions (str or list): Either "NA", a comma-separated string,
            or a list of proportions from a CSV.
    """

    if csv:
        proportions = [float(p) for p in proportions]

    elif proportions == "NA":
        if len(sample_names) == 1:
            print(
                "Only one sample provided. "
                "Using 1.0 as the sample proportion."
            )
            proportions = [1.0]

        else:
            print(
                "Read simulation proportions not provided. "
                "Generating proportions randomly..."
            )

            proportions = generate_random_values(len(sample_names))

            with open(
                os.path.join(outdir, "sample_proportions.txt"),
                "w",
            ) as file:
                for name, proportion in zip(sample_names, proportions):
                    file.write(f"{name}: {proportion}\n")

    else:
        proportions = [float(x) for x in proportions.split(",")]

    # Validate lengths
    if not (
        len(sample_names)
        == len(proportions)
        == len(sample_paths)
    ):
        raise Exception(
            "Number of samples, proportions, and sample paths should match!"
        )

    # Validate proportions sum
    if round(sum(proportions), 6) != 1.0:
        raise Exception(
            "Sum of all proportions should equal 1.0!"
        )

    return proportions


def check_dir(outdir, redo, sample_names):
    # Define the specific targets we want to protect/overwrite
    targets = [
        os.path.join(outdir, "reads_1.fastq"),
        os.path.join(outdir, "reads_2.fastq")
    ]
    for sample in sample_names:
        targets.append(os.path.join(outdir, sample))

    # Check if any of these specific targets already exist
    existing_targets = [t for t in targets if os.path.exists(t)]

    if existing_targets:
        if not redo:
            print("Error: The following output files/"
                  f"directories already exist in '{outdir}':")
            for t in existing_targets:
                print(f"  - {os.path.basename(t)}")
            print("Use --redo to overwrite them.")
            sys.exit(1)
        else:
            print("Notice: Overwriting existing outputs "
                  f"in '{outdir}' because --redo was set.")
            for t in existing_targets:
                if os.path.isdir(t):
                    shutil.rmtree(t)
                else:
                    os.remove(t)
            os.makedirs(outdir, exist_ok=True)
    else:
        os.makedirs(outdir, exist_ok=True)


def validate_simulation_args(simulation_mode, primers, reference,
                             proportions, proportions_csv, genomes):
    if simulation_mode == "amplicon" and primers == "NA":
        print("Primer file is required for simulation mode amplicon")
        sys.exit(1)
    if simulation_mode == "metagenomics" and primers != "NA":
        print("Primer file not needed for metagenomics simulation")
        sys.exit(1)
    if simulation_mode == "metagenomics" and reference != "NA":
        print("Reference file not needed for metagenomics simulation")
        sys.exit(1)
    if simulation_mode == "amplicon" and reference == "NA":
        print("Reference file is required for simulation mode amplicon")
        sys.exit(1)
    if proportions_csv != "NA":
        if proportions != "NA":
            print("Cannot use --proportions with --csv")
            sys.exit(1)

        if genomes != "NA":
            print("Cannot use --genomes with --csv")
            sys.exit(1)


def assess_genome_quality_from_fasta(fasta_path):
    """
    Parses a FASTA genome file and assesses quality by:
    - Counting ambiguous (non-ACGT) bases.
    - Reporting the length of each contig.
    - Emitting warnings for ambiguous bases or multiple contigs.
    - Printing a genome quality summary.

    Parameters:
        fasta_path (str): Path to the FASTA file.

    Returns:
        dict: {
            'total_ambiguous_bases': int,
            'contig_lengths': dict of {contig_id: length}
        }
    """

    ambiguous_bases = {'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N'}
    total_ambiguous = 0
    contig_lengths = {}

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        contig_lengths[record.id] = len(seq)
        total_ambiguous += sum(1 for base in seq if base in ambiguous_bases)

    report = {
        'total_ambiguous_bases': total_ambiguous,
        'contig_lengths': contig_lengths
    }

    num_contigs = len(report['contig_lengths'])
    num_ambiguous = report['total_ambiguous_bases']

    # Warnings
    if num_ambiguous > 0:
        warnings.warn(
            f"{fasta_path}: Contains {num_ambiguous} ambiguous base(s). "
            "Please choose a better quality genome."
        )

    if num_contigs > 1:
        warnings.warn(
            f"{fasta_path}: Contains {num_contigs} contigs. "
            "Does your organism have more than one chromosome? "
            "Are you providing high quality assemblies?"
        )

    # Print results
    print(f"\nGenome: {fasta_path}")
    print(f"  Total ambiguous bases: {num_ambiguous}")
    print("  Contig lengths:")

    for contig, length in report['contig_lengths'].items():
        print(f"    {contig}: {length}")

    return report


def extract_sequence(reference, chrom, start, end):
    """Extracts sequence from reference based on coordinates and strand."""
    if reference.id not in chrom:
        raise ValueError(f"Chromosome {chrom} not found in reference.")
    seq = reference.seq[start:end]  # Extract sequence
    return str(seq)


def validate_primer_bed(df):
    """
    Validates a DataFrame representing a primer BED file.

    Args:
        df (pd.DataFrame): DataFrame to validate.

    Returns:
        pd.DataFrame: The original DataFrame if valid.

    Raises:
        ValueError: If the DataFrame does not conform to expected format.
    """

    # Check column count
    if df.shape[1] < 6 or df.shape[1] > 7:
        raise ValueError(
            "DataFrame must have 6 or 7"
            " columns. Please refer to example primer file."
        )

    # Check column types
    if not all(isinstance(df.iloc[i, 0], str) for i in range(len(df))):
        raise ValueError(
            "First column must contain only strings.(Chromosome name)")

    if not pd.api.types.is_numeric_dtype(df.iloc[:, 1]):
        raise ValueError("Second column must be numeric.(Primer start)")

    if not pd.api.types.is_numeric_dtype(df.iloc[:, 2]):
        raise ValueError("Third column must be numeric.(Primer end)")
    # Check that third column is greater than second column
    if not (df.iloc[:, 2] > df.iloc[:, 1]).all():
        raise ValueError(
            "Third column values must be greater than second column values."
            "Primer start coordinates cannot be greater than primer end."
        )

    # Check fourth column format
    # This allows optional alternative
    # suffixes like _alt or _1 after LEFT/RIGHT,
    # but prevents the double-number format in the middle.
    pattern = re.compile(r'^[A-Za-z0-9-]+_\d+_(LEFT|RIGHT)(_[A-Za-z0-9]+)?$')

    # Strip any leading/trailing spaces and ensure the values are strings
    if not all(
        isinstance(val, str) and
            pattern.match(val.strip()) for val in df.iloc[:, 3]
    ):
        raise ValueError(
            "Fourth column format is incorrect."
            " Expected 'string_number_LEFT/RIGHT'"
            " with possible extra characters"
            " at the end indicating whether the primer is alternative."
        )

    if not pd.api.types.is_numeric_dtype(df.iloc[:, 4]):
        raise ValueError("Fifth column must be numeric.(strand +/-)")

    if not all(isinstance(df.iloc[i, 5], str) for i in range(len(df))):
        raise ValueError(
            "Sixth column must contain only strings.(Primer sequence)")
    return df


def generate_random_values(N):
    # Generate N random values
    random_values = np.random.rand(N)

    # Normalize the values so that they sum to 1
    normalized_values = random_values / np.sum(random_values)

    # Convert the numpy array to a list
    return normalized_values.tolist()


def merge_fastq_files(fastq_file, output_file):
    """
    Merges a FASTQ file into an output FASTQ file using subprocess.call.

    Parameters:
    fastq_file (str): Path to the input FASTQ file.
    output_file (str): Path to the output FASTQ file.
    """
    with open(fastq_file, "rb") as src, open(output_file, "ab") as dst:
        shutil.copyfileobj(src, dst)


def create_valid_primer_combinations(df):
    # Ensure the column exists before assignment
    if "valid_combinations" not in df.columns:
        df["valid_combinations"] = None

    valid_primers = []  # Use a list instead of concatenating DataFrames

    for i in range(len(df)):
        # Pair coordinates with their mismatch maps
        left_coords = zip(df.at[i, "left_primer_loc"],
                          df.at[i, "left_match"])
        right_coords = zip(df.at[i, "right_primer_loc"],
                           df.at[i, "right_match"])
        # Safe assignment using .at[]
        df.at[i, "valid_combinations"] = evaluate_matches(left_coords,
                                                          right_coords)
        for (
            primer_start,
            primer_end,
            left_match,
            right_match
        ) in df.at[i, "valid_combinations"]:
            valid_primers.append(
                {
                    "amplicon_number": df.at[i, "amplicon_number"],
                    "primer_start": primer_start,
                    "primer_end": primer_end,
                    "left_match": left_match,
                    "right_match": right_match
                }
            )

    # Check if we found any valid primers
    if not valid_primers:
        raise ValueError(
            "No primer matches found, please check your primer file.")

    # Convert collected data to DataFrame efficiently
    valid_primers_df = pd.DataFrame.from_records(valid_primers)
    # Merge with original DataFrame to include additional columns
    df = df[["amplicon_number", "primer_seq_x",
             "primer_seq_y", "ambiguous_bases"]]
    all_amplicons = df.merge(
        valid_primers_df,
        how="left",
        on="amplicon_number",
    )

    return all_amplicons


def preprocess_primers(primer_file, reference):
    # define column names to read primers bed file
    col_names = [
        "ref",
        "start",
        "end",
        "left_right",
        "primer_pool",
        "strand",
        "primer_seq",
    ]

    primer_bed = pd.read_csv(
        primer_file,
        sep="\t",
        names=col_names,
        comment="#",
    )

    primer_bed = validate_primer_bed(primer_bed)

    if (
        "primer_seq" not in primer_bed.columns
        or primer_bed["primer_seq"].isna().all()
        or (primer_bed["primer_seq"].astype(str).str.strip() == "").all()
    ):
        warnings.warn(
            "primer_seq column missing or empty; "
            "extracting sequences from reference.. "
            "This is not recommended.",
            UserWarning,
        )
        primer_bed["primer_seq"] = [
            extract_sequence(reference, row.ref, row.start, row.end)
            for row in primer_bed.itertuples(index=False)
            ]
    else:
        neg_strand = primer_bed["strand"] == "-"
        _rc = str.maketrans("ACGTacgt", "TGCAtgca")
        primer_bed.loc[neg_strand, "primer_seq"] = [
            seq.translate(_rc)[::-1]
            for seq in primer_bed.loc[neg_strand, "primer_seq"]
            ]
    # split the amplicon name into number and left/right
    primer_bed["amplicon_number"] = primer_bed["left_right"].str.split(
        "_").str[1]
    # merge the df with itself to have right and left primer on one row
    df = pd.merge(
        primer_bed.loc[primer_bed["left_right"].str.contains("LEFT")],
        primer_bed.loc[primer_bed["left_right"].str.contains("RIGHT")],
        on=["amplicon_number", "primer_pool"],
    )
    # select needed columns
    df = df[["amplicon_number", "primer_seq_x", "primer_seq_y"]]
    mask = df.duplicated(subset=["amplicon_number"], keep=False)
    # Apply a function to append an index for duplicated amplicon_number values
    df.loc[mask, "amplicon_number"] = (
        df.loc[mask, "amplicon_number"].astype(str)
        + "_"
        + df.loc[mask].groupby("amplicon_number").cumcount().astype(str)
    )
    return df


def count_contigs(fasta_file):
    """Counts the number of contigs in the FASTA file."""
    contig_count = 0
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                contig_count += 1
    return contig_count


def run_simulation_on_fasta(
    fasta_file,
    output_dir,
    read_cnt,
    simulator,
    wgsim_insert_size,
    wgsim_read_length,
    wgsim_error_rate,
    extra_flags=None
):
    """Runs simulator on a single FASTA file with the given parameters."""
    # Count the number of contigs in the FASTA file
    num_contigs = count_contigs(fasta_file)

    if num_contigs == 0:
        raise ValueError("No contigs found in the FASTA file.")

    # Calculate the number of reads per contig
    reads_per_contig = read_cnt // num_contigs  # Integer division

    # Prepare output filenames
    output_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
    merged_output1 = os.path.join(output_dir, "merged_reads_1.fastq")
    merged_output2 = os.path.join(output_dir, "merged_reads_2.fastq")

    # Initialize merged output files (if not already created)
    open(merged_output1, "a").close()  # Create or append empty file
    open(merged_output2, "a").close()  # Create or append empty file

    # Loop through the contigs and run wgsim for each
    for contig_idx in range(num_contigs):
        # Temporary output files for each contig
        output1 = os.path.join(
            output_dir, f"{output_prefix}_contig{contig_idx + 1}_1.fastq"
        )
        output2 = os.path.join(
            output_dir, f"{output_prefix}_contig{contig_idx + 1}_2.fastq"
        )

        if simulator == "wgsim":
            command = [
                "wgsim",
                "-N",
                str(reads_per_contig),
                "-d",
                str(wgsim_insert_size),
                "-e",
                str(wgsim_error_rate),
                "-1",
                str(wgsim_read_length),
                "-2",
                str(wgsim_read_length),
                fasta_file,
                output1,
                output2,
            ]
            if extra_flags:
                command.extend(extra_flags)

        else:
            # Adjust Mason command
            command = [
                "mason_simulator",
                "-ir",
                fasta_file,
                "-n",
                str(int(reads_per_contig)),
                "-o",
                output1,
                "-or",
                output2,
            ]
            if extra_flags:
                command.extend(extra_flags)
        # Run the simulator command and capture any errors
        try:
            subprocess.run(
                command, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running the command: {e}")
        # Merge the contig-specific output into the final merged output files
        merge_fastq_files(output1, merged_output1)
        merge_fastq_files(output2, merged_output2)


def run_simulation_on_fasta_single_genome(
    fasta_file,
    output_dir,
    read_cnt,
    simulator,
    extra_flags=None
):
    """Runs simulator on a single FASTA file with the given parameters."""

    output1 = os.path.join(
            output_dir, "reads_1.fastq"
    )
    output2 = os.path.join(
            output_dir, "reads_2.fastq"
    )

    if simulator == "wgsim":
        command = [
            "wgsim",
            "-N",
            str(read_cnt),
            fasta_file,
            output1,
            output2,
        ]
        if extra_flags:
            command.extend(extra_flags)
    else:
        # Adjust Mason command
        command = [
            "mason_simulator",
            "-ir",
            fasta_file,
            "-n",
            str(read_cnt),
            "-o",
            output1,
            "-or",
            output2,
        ]
        if extra_flags:
            command.extend(extra_flags)
    # Run the simulator command and capture any errors
    try:
        subprocess.run(
            command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running the command: {e}")


# Extended ambiguity-aware mismatch display
def mismatch_alignment(primer, matched_seq):
    """
    Returns matched sequence with mismatches shown in parentheses.
    Also returns a flag if primer contains any ambiguous base.
    """
    ambiguous_bases = {'R', 'Y', 'S', 'W',
                       'K', 'M', 'B', 'D',
                       'H', 'V', 'N'}
    has_ambiguity = any(base in ambiguous_bases
                        for base in matched_seq.upper())
    aligned = []
    for p, m in zip(primer.upper(), matched_seq.upper()):
        aligned.append(m if p == m else f"({m})")
    return "".join(aligned), has_ambiguity


def find_closest_primer_match(df, reference_seq, maxmismatch):
    """
    For each row in df, find all left/right primer match positions (as lists),
    allowing up to `maxmismatch` mismatches. Ensures both primers are found
    on the same strand. Returns original df columns + matches, mismatch maps,
    strand, and whether primers contain ambiguous bases (IUPAC codes).
    """
    results = []
    warned = False
    for row in df.itertuples(index=False):
        primer_left = row.primer_seq_x
        primer_right = row.primer_seq_y

        pattern_left = f"({primer_left}){{s<={maxmismatch}}}"
        pattern_right = f"({primer_right}){{s<={maxmismatch}}}"

        # --- 1. INITIALIZE REVERSE VARIABLES UP FRONT ---
        left_rev, right_rev = [], []
        left_rev_actual, right_rev_actual = [], []
        left_rev_mismatch_map, right_rev_mismatch_map = [], []
        left_rev_has_ambig, right_rev_has_ambig = [], []

        # --- 2. FORWARD STRAND SEARCH ---
        left_fwd_data = [
            (m.start(), m.group())
            for m in re.finditer(pattern_left, reference_seq,
                                 flags=re.IGNORECASE, overlapped=True)
        ]
        right_fwd_data = [
            (m.start(), m.group())
            for m in re.finditer(pattern_right, reference_seq,
                                 flags=re.IGNORECASE, overlapped=True)
        ]
        left_fwd = [pos for pos, _ in left_fwd_data]
        right_fwd = [pos for pos, _ in right_fwd_data]

        left_fwd_actual = [seq for _, seq in left_fwd_data]
        right_fwd_actual = [seq for _, seq in right_fwd_data]

        left_fwd_mismatch_map = []
        left_fwd_has_ambig = []
        for seq in left_fwd_actual:
            aligned, has_ambig = mismatch_alignment(primer_left, seq)
            left_fwd_mismatch_map.append(aligned)
            left_fwd_has_ambig.append(has_ambig)
        right_fwd_mismatch_map = []
        right_fwd_has_ambig = []
        for seq in right_fwd_actual:
            aligned, has_ambig = mismatch_alignment(primer_right, seq)
            right_fwd_mismatch_map.append(aligned)
            right_fwd_has_ambig.append(has_ambig)

        # --- 3. REVERSE STRAND SEARCH (IF FORWARD FAILS) ---
        if not left_fwd or not right_fwd:
            complement_table = str.maketrans("ATCGatcg", "TAGCtagc")
            # Reverse complement without Bio Python
            right_rev_seq = primer_left.translate(complement_table)[::-1]
            left_rev_seq = primer_right.translate(complement_table)[::-1]
            pattern_left_rev = f"({left_rev_seq}){{s<={maxmismatch}}}"
            pattern_right_rev = f"({right_rev_seq}){{s<={maxmismatch}}}"
            left_rev_data = [
                (m.start(), m.group())
                for m in re.finditer(pattern_left_rev, reference_seq,
                                     flags=re.IGNORECASE, overlapped=True)
            ]
            right_rev_data = [
                (m.start(), m.group())
                for m in re.finditer(pattern_right_rev, reference_seq,
                                     flags=re.IGNORECASE, overlapped=True)
            ]
            left_rev = [pos for pos, _ in left_rev_data]
            right_rev = [pos for pos, _ in right_rev_data]

            left_rev_actual = [seq for _, seq in left_rev_data]
            right_rev_actual = [seq for _, seq in right_rev_data]
            for seq in left_rev_actual:
                aligned, has_ambig = mismatch_alignment(left_rev_seq, seq)
                left_rev_mismatch_map.append(aligned)
                left_rev_has_ambig.append(has_ambig)
            for seq in right_rev_actual:
                aligned, has_ambig = mismatch_alignment(right_rev_seq, seq)
                right_rev_mismatch_map.append(aligned)
                right_rev_has_ambig.append(has_ambig)

        # --- 4. AMBIGUOUS BASE CHECK ---
        has_ambiguous_base = any([
            any(b in primer_left.upper() for b in "RYSWKMBDHVN"),
            any(b in primer_right.upper() for b in "RYSWKMBDHVN"),
            any(left_fwd_has_ambig),
            any(right_fwd_has_ambig),
            any(left_rev_has_ambig),
            any(right_rev_has_ambig),
        ])
        if has_ambiguous_base and not warned:
            warnings.warn("One or more primers contain ambiguous "
                          "bases (e.g., N, R, Y, etc). "
                          "Matches may be unreliable.")
            warned = True

        # --- 5. BUILD RESULT DICTIONARY ---
        result_row = row._asdict()
        result_row["ambiguous_bases"] = has_ambiguous_base

        if left_fwd and right_fwd:
            result_row.update({
                "left_primer_loc": left_fwd,
                "right_primer_loc": right_fwd,
                "left_seq_actual": left_fwd_actual,
                "right_seq_actual": right_fwd_actual,
                "left_match": left_fwd_mismatch_map,
                "right_match": right_fwd_mismatch_map,
            })
        elif left_rev and right_rev:  # Fixed variable naming conflict here
            result_row.update({
                "left_primer_loc": left_rev,
                "right_primer_loc": right_rev,
                "left_seq_actual": left_rev_actual,
                "right_seq_actual": right_rev_actual,
                "left_match": left_rev_mismatch_map,
                "right_match": right_rev_mismatch_map,
            })
        else:
            result_row.update({
                "left_primer_loc": [],
                "right_primer_loc": [],
                "left_seq_actual": [],
                "right_seq_actual": [],
                "left_match": [],
                "right_match": []
            })

        results.append(result_row)
    return pd.DataFrame(results)


def make_amplicon(left_primer_loc, right_primer_loc, primer_seq_y, reference):
    """function to create amplicons based on the string match location"""
    if math.isnan(left_primer_loc) or math.isnan(right_primer_loc):
        # if there is either no left and right primer match
        amplicon = ""
        # length of the right primer is added to include the
        # right primer in the amplicon
    else:
        amplicon = str(
            reference[
                int(left_primer_loc - 1):
                int(right_primer_loc + len(primer_seq_y))
            ]
        )
    return amplicon


def evaluate_matches(left_primer_coordinates, right_primer_coordinates):
    """Find valid primer pairs that can produce an amplicon with mismatches."""
    if left_primer_coordinates and right_primer_coordinates:
        valid_combinations = []
        combinations = itertools.product(left_primer_coordinates,
                                         right_primer_coordinates)
        for left, right in combinations:
            left_pos, left_match = left
            right_pos, right_match = right
            amplicon_length = right_pos - left_pos

            if 0 < amplicon_length <= 2000:
                valid_combinations.append(
                    (left_pos, right_pos,
                     left_match,
                     right_match)
                )

        return valid_combinations
    else:
        return []


def write_fasta_group(group, amplicon_number, output_dir):
    fasta_filename = os.path.join(
        output_dir, f"amplicon_{amplicon_number}.fasta")

    filtered_records = [
        SeqRecord(Seq(seq), id=f"{amplicon_number}_{i}", description="")
        for i, seq in enumerate(group["amplicon_sequence"])
        if 1 < len(seq) < 10000
    ]

    if filtered_records:
        SeqIO.write(filtered_records, fasta_filename, "fasta")
    else:
        print(
            f"No valid sequences for amplicon {amplicon_number},"
            "no file written.please checkout the amplicon_stats.csv "
            "file for more information."
        )


def process_amplicon_worker(args):
    """Worker for the 'amplicon' simulation mode."""
    (name, path, cnt, df_primers_template, maxmismatch, outdir,
     simulator, wgsim_insert_size, wgsim_read_length, wgsim_error_rate,
     extra_simulator_flags) = args

    sample_amplicons_list = []

    for genome_seq in SeqIO.parse(path, "fasta"):
        contig_df = find_closest_primer_match(
            df_primers_template.copy(),
            str(genome_seq.seq),
            maxmismatch
        )
        all_amplicons = create_valid_primer_combinations(contig_df)
        all_amplicons = all_amplicons.fillna(0)

        all_amplicons["amplicon_length"] = np.where(
            (all_amplicons["primer_start"] != 0) &
            (all_amplicons["primer_end"] != 0),
            all_amplicons["primer_end"] -
            all_amplicons["primer_start"] +
            all_amplicons["primer_seq_y"].str.len(),
            0,
        )

        all_amplicons["amplicon_sequence"] = all_amplicons.apply(
            lambda row: make_amplicon(
                row["primer_start"],
                row["primer_end"],
                row["primer_seq_y"],
                genome_seq.seq,
            ), axis=1,
        )
        all_amplicons["contig_id"] = genome_seq.id
        sample_amplicons_list.append(all_amplicons)

    if not sample_amplicons_list:
        return ("warning", f"Warning: No sequences found in {path}")

    full_sample_df = pd.concat(sample_amplicons_list, ignore_index=True)
    amp_out_dir = os.path.join(outdir, name, "amplicons")
    os.makedirs(amp_out_dir, exist_ok=True)
    full_sample_df.to_csv(os.path.join(amp_out_dir,
                                       "amplicon_stats.csv"),
                          index=False)

    full_sample_df["amplicon_suffix"] = full_sample_df[
        "amplicon_number"].apply(
        lambda x: x.split("_")[0] if "_" in x else x
    )
    for n, g in full_sample_df.groupby("amplicon_suffix"):
        write_fasta_group(g, n, amp_out_dir)

    read_dir = os.path.join(outdir, name, "reads")
    os.makedirs(read_dir, exist_ok=True)

    fasta_files = [
        os.path.join(amp_out_dir, f)
        for f in os.listdir(amp_out_dir)
        if f.endswith((".fasta", ".fa"))
    ]

    for fasta_file in fasta_files:
        run_simulation_on_fasta(
            fasta_file, read_dir, cnt, simulator,
            wgsim_insert_size, wgsim_read_length, wgsim_error_rate,
            extra_flags=extra_simulator_flags
        )

    # Expected paths for merging step in main thread
    read_path1 = os.path.join(os.path.abspath(outdir),
                              name, "reads/merged_reads_1.fastq")
    read_path2 = os.path.join(os.path.abspath(outdir),
                              name, "reads/merged_reads_2.fastq")
    return ("success", name, read_path1, read_path2)


def process_genome_worker(args):
    """Worker for the default/standard genome simulation mode (else clause)."""
    name, path, cnt, outdir, simulator, extra_simulator_flags = args
    read_dir = os.path.join(outdir, name, "reads")
    os.makedirs(read_dir, exist_ok=True)

    run_simulation_on_fasta_single_genome(
        path,
        read_dir,
        cnt,
        simulator,
        extra_flags=extra_simulator_flags
    )

    # Expected paths for merging step in main thread
    read_path1 = os.path.join(os.path.abspath(outdir),
                              name, "reads/reads_1.fastq")
    read_path2 = os.path.join(os.path.abspath(outdir),
                              name, "reads/reads_2.fastq")
    return ("success", name, read_path1, read_path2)
