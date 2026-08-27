import unittest
import os


def file_exists(directory, filename):
    file_path = os.path.join(directory, filename)
    return os.path.exists(file_path)


class CommandLineTests(unittest.TestCase):
    def test_version(self):
        os.system("bygul --version")

    def test_simulation(self):
        os.system(
            "bygul simulate-proportions \
            --genomes bygul/tests/data/ATM-2FFMD73N3.fasta \
            --primers bygul/tests/data/ARTIC_V4-1.bed --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_simulator(self):
        os.system(
            "bygul simulate-proportions \
            --genomes bygul/tests/data/ATM-2FFMD73N3.fasta \
            --primers bygul/tests/data/ARTIC_V4-1.bed \
            --simulator wgsim --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_max_mismatch(self):
        os.system(
            "bygul simulate-proportions \
            --genomes bygul/tests/data/ATM-2FFMD73N3.fasta \
            --primers bygul/tests/data/ARTIC_V4-1.bed \
            --maxmismatch 2 --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_readcnt(self):
        os.system(
            "bygul simulate-proportions \
            --genomes bygul/tests/data/ATM-2FFMD73N3.fasta \
            --primers bygul/tests/data/ARTIC_V4-1.bed \
            --readcnt 200 --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_readlength(self):
        os.system(
            "bygul simulate-proportions \
            --genomes bygul/tests/data/ATM-2FFMD73N3.fasta \
            --primers bygul/tests/data/ARTIC_V4-1.bed \
            --readlength 150 --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_2genomes(self):
        os.system(
            "bygul simulate-proportions "
            "--genomes bygul/tests/data/ATM-2FFMD73N3.fasta,"
            "bygul/tests/data/KR-SEARCH-120354.fasta "
            "--primers bygul/tests/data/ARTIC_V4-1.bed --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_2genomes_proportions(self):
        os.system(
            "bygul simulate-proportions "
            "--genomes bygul/tests/data/ATM-2FFMD73N3.fasta,"
            "bygul/tests/data/KR-SEARCH-120354.fasta "
            "--primers bygul/tests/data/ARTIC_V4-1.bed "
            "--proportions 0.8,0.2 --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_metagenomics(self):
        os.system(
            "bygul simulate-proportions "
            "--genomes bygul/tests/data/ATM-2FFMD73N3.fasta,"
            "bygul/tests/data/KR-SEARCH-120354.fasta "
            "--simulation_mode metagenomics "
            "--proportions 0.8,0.2 --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_mason(self):
        os.system(
            "bygul simulate-proportions "
            "--primers bygul/tests/data/ARTIC_V4-1.bed "
            "--genomes bygul/tests/data/ATM-2FFMD73N3.fasta,"
            "bygul/tests/data/KR-SEARCH-120354.fasta "
            "--proportions 0.8,0.2 --redo --simulator mason "
            "--read_length 200 --insert_size_sd 50"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_art(self):
        os.system(
            "bygul simulate-proportions "
            "--primers bygul/tests/data/ARTIC_V4-1.bed "
            "--genomes bygul/tests/data/ATM-2FFMD73N3.fasta,"
            "bygul/tests/data/KR-SEARCH-120354.fasta "
            "--proportions 0.8,0.2 --redo --simulator art "
            "--read_length 200 --insert_size_sd 50"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_art_metagenomics(self):
        os.system(
            "bygul simulate-proportions "
            "--genomes bygul/tests/data/ATM-2FFMD73N3.fasta,"
            "bygul/tests/data/KR-SEARCH-120354.fasta "
            "--simulation_mode metagenomics "
            "--proportions 0.8,0.2 --redo --simulator art "
            "--read_length 200 --insert_size 200"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_csv(self):
        os.system(
            "bygul simulate-proportions "
            "--primers bygul/tests/data/ARTIC_V4-1.bed "
            "--csv bygul/tests/data/sample_proportions.csv "
            "--multifasta bygul/tests/data/sample_genomes.fasta "
            "--redo --simulator mason "
            "--reference bygul/tests/data/reference.fasta "
            "--read_length 200 --insert_size_sd 50"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_csv_metagenomics(self):
        os.system(
            "bygul simulate-proportions "
            "--csv bygul/tests/data/sample_proportions.csv "
            "--multifasta bygul/tests/data/sample_genomes.fasta "
            "--redo --simulation_mode metagenomics "
            "--read_length 200 --insert_size 200"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_check_primers(self):
        os.system(
            "bygul check-primers "
            "bygul/tests/data/ATM-2FFMD73N3.fasta "
            "bygul/tests/data/ARTIC_V4-1.bed "
        )
        self.assertTrue(file_exists(".", "amplicon_stats.csv"))


if __name__ == "__main__":
    unittest.main()
