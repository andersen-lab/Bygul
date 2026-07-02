import pandas as pd
import os
import unittest


class ProportionTests(unittest.TestCase):

    def test_prop_with_freyja(self):
        # 1. Run the simulation
        os.system(
            "bygul simulate-proportions "
            "--genomes bygul/tests/data/BCN-SEARCH-105346.fasta,"
            "bygul/tests/data/CA-SEARCH-43254.fasta "
            "--primers bygul/tests/data/ARTIC_V4-1.bed "
            "--reference bygul/tests/data/reference.fasta "
            "--proportions 0.8,0.2 --redo")
        self.assertTrue(os.path.exists("results/reads_1.fastq"))

        # 2. Run Freyja Pipeline (Align -> Variants -> Demix)
        # Note: Replace 'minimap2' with your preferred aligner
        ref = "bygul/tests/data/reference.fasta"
        os.system(f"minimap2 -ax sr {ref} "
                  "results/reads_1.fastq results/reads_2.fastq | "
                  "samtools sort -o results/merged.bam")
        os.system("samtools index results/merged.bam")
        # Freyja commands
        os.system("freyja variants results/merged.bam "
                  "--variants results/variants.tsv "
                  f"--depths results/depths.tsv --ref {ref}")
        os.system("freyja update --outdir .")
        os.system("freyja demix results/variants.tsv "
                  "results/depths.tsv --output "
                  "results/demix.tsv --depthcutoff 20 "
                  "--lineageyml lineages.yml")

        # 3. Validate Proportions
        self.assertTrue(os.path.exists("results/demix.tsv"))

        # Read the Freyja output
        df = pd.read_csv("results/demix.tsv", sep='\t', index_col=0)

        # Access the string of numbers in the 'abundances' row
        abundances_str = df.loc['abundances'].values[0]

        abundances = [float(x) for x in abundances_str.split()]

        top_abundance = abundances[0]
        second_abundance = abundances[1]

        self.assertAlmostEqual(top_abundance, 0.8, delta=0.1)
        self.assertAlmostEqual(second_abundance, 0.2, delta=0.1)


if __name__ == '__main__':
    unittest.main()
