import os
import unittest

import pandas as pd


class ProportionTests(unittest.TestCase):

    def _run_freyja_and_check(self, expected):
        """Run Freyja demix and verify expected abundances."""
        ref = "bygul/tests/data/reference.fasta"

        os.system(
            f"minimap2 -ax sr {ref} "
            "results/reads_1.fastq results/reads_2.fastq | "
            "samtools sort -o results/merged.bam"
        )
        os.system("samtools index results/merged.bam")

        os.system(
            "freyja variants results/merged.bam "
            "--variants results/variants.tsv "
            f"--depths results/depths.tsv --ref {ref}"
        )
        os.system("freyja update --outdir .")
        os.system(
            "freyja demix results/variants.tsv "
            "results/depths.tsv "
            "--output results/demix.tsv "
            "--depthcutoff 20 "
            "--lineageyml lineages.yml"
        )

        self.assertTrue(os.path.exists("results/demix.tsv"))

        df = pd.read_csv("results/demix.tsv", sep="\t", index_col=0)
        abundances = [
            float(x) for x in df.loc["abundances"].values[0].split()
        ]

        for observed, exp in zip(abundances, expected):
            self.assertAlmostEqual(observed, exp, delta=0.1)

    def test_prop_with_freyja(self):
        os.system(
            "bygul simulate-proportions "
            "--genomes bygul/tests/data/BCN-SEARCH-105346.fasta,"
            "bygul/tests/data/CA-SEARCH-43254.fasta "
            "--primers bygul/tests/data/ARTIC_V4-1.bed "
            "--reference bygul/tests/data/reference.fasta "
            "--proportions 0.8,0.2 "
            "--redo"
        )

        self.assertTrue(os.path.exists("results/reads_1.fastq"))
        self._run_freyja_and_check([0.8, 0.2])

    def test_prop_with_metagenomics(self):
        os.system(
            "bygul simulate-proportions "
            "--genomes bygul/tests/data/BCN-SEARCH-105346.fasta,"
            "bygul/tests/data/CA-SEARCH-43254.fasta "
            "--simulation_mode metagenomics "
            "--proportions 0.8,0.2 --redo"
        )

        self.assertTrue(os.path.exists("results/reads_1.fastq"))
        self._run_freyja_and_check([0.8, 0.2])

    def test_prop_with_metagenomics_csv(self):
        os.system(
            "bygul simulate-proportions "
            "--csv bygul/tests/data/sample_proportions.csv "
            "--multifasta bygul/tests/data/sample_genomes.fasta "
            "--redo --simulation_mode metagenomics "
            "--illumina-read-length 200"
        )

        self.assertTrue(os.path.exists("results/reads_1.fastq"))
        self._run_freyja_and_check([0.8, 0.2])

    def test_simulation_with_csv(self):
        os.system(
            "bygul simulate-proportions "
            "--primers bygul/tests/data/ARTIC_V4-1.bed "
            "--csv bygul/tests/data/sample_proportions.csv "
            "--multifasta bygul/tests/data/sample_genomes.fasta "
            "--redo "
            "--simulator mason "
            "--reference bygul/tests/data/reference.fasta "
            "--illumina-read-length 200"
        )

        self.assertTrue(os.path.exists("results/reads_1.fastq"))

        # Expected values should match sample_proportions.csv
        self._run_freyja_and_check([0.8, 0.2])

    def test_simulation_with_multiple_contig_samples(self):
        os.system(
            "bygul simulate-proportions "
            "--primers bygul/tests/data/ARTIC_V4-1.bed "
            "--csv bygul/tests/data/sample_proportions.csv "
            "--multifasta bygul/tests/data/sample_genomes_multi_contig.fasta "
            "--redo "
            "--simulator mason "
            "--reference bygul/tests/data/reference.fasta "
            "--illumina-read-length 200"
        )

        self.assertTrue(os.path.exists("results/reads_1.fastq"))

        # Expected values should match sample_proportions.csv
        self._run_freyja_and_check([0.8, 0.2])


if __name__ == "__main__":
    unittest.main()
