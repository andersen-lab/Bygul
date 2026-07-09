import pandas as pd
import os
import unittest


class ProportionTests(unittest.TestCase):

    def run_freyja(self):
        ref = "bygul/tests/data/reference.fasta"

        self.assertEqual(
            os.system(
                f"minimap2 -ax sr {ref} "
                "results/reads_1.fastq results/reads_2.fastq | "
                "samtools sort -o results/merged.bam"
            ),
            0,
        )

        self.assertEqual(
            os.system("samtools index results/merged.bam"),
            0,
        )

        self.assertEqual(
            os.system(
                "freyja variants results/merged.bam "
                "--variants results/variants.tsv "
                f"--depths results/depths.tsv --ref {ref}"
            ),
            0,
        )

        self.assertEqual(
            os.system("freyja update --outdir ."),
            0,
        )

        self.assertEqual(
            os.system(
                "freyja demix results/variants.tsv "
                "results/depths.tsv "
                "--output results/demix.tsv "
                "--depthcutoff 20 "
                "--lineageyml lineages.yml"
            ),
            0,
        )

        self.assertTrue(os.path.exists("results/demix.tsv"))

        df = pd.read_csv(
            "results/demix.tsv",
            sep="\t",
            index_col=0,
        )

        abundances = [
            float(x)
            for x in df.loc["abundances"].values[0].split()
        ]

        return sorted(abundances, reverse=True)

    def test_prop_with_freyja_manual_proportions(self):
        # Test explicit --proportions workflow
        self.assertEqual(
            os.system(
                "bygul simulate-proportions "
                "--genomes bygul/tests/data/BCN-SEARCH-105346.fasta,"
                "bygul/tests/data/CA-SEARCH-43254.fasta "
                "--primers bygul/tests/data/ARTIC_V4-1.bed "
                "--reference bygul/tests/data/reference.fasta "
                "--proportions 0.8,0.2 "
                "--redo"
            ),
            0,
        )

        self.assertTrue(os.path.exists("results/reads_1.fastq"))

        abundances = self.run_freyja()

        self.assertAlmostEqual(abundances[0], 0.8, delta=0.1)
        self.assertAlmostEqual(abundances[1], 0.2, delta=0.1)

    def test_prop_with_freyja_csv_mason(self):
        # Test CSV + multifasta workflow
        self.assertEqual(
            os.system(
                "bygul simulate-proportions "
                "--csv bygul/tests/data/sample_proportions.csv "
                "--multifasta bygul/tests/data/sample_genomes.fasta "
                "--redo "
                "--simulator mason "
                "--illumina-read-length 200"
            ),
            0,
        )

        self.assertTrue(os.path.exists("results/reads_1.fastq"))

        abundances = self.run_freyja()

        expected = pd.read_csv(
            "bygul/tests/data/sample_proportions.csv"
        )

        expected_props = sorted(
            expected["proportion"].tolist(),
            reverse=True,
        )

        self.assertEqual(
            len(abundances),
            len(expected_props),
        )

        for observed, expected in zip(abundances, expected_props):
            self.assertAlmostEqual(
                observed,
                expected,
                delta=0.1,
            )


if __name__ == "__main__":
    unittest.main()
