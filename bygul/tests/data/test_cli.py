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
            bygul/tests/data/ATM-2FFMD73N3.fasta \
            bygul/tests/data/ARTIC_V4-1.bed \
            bygul/tests/data/reference.fasta --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_simulator(self):
        os.system(
            "bygul simulate-proportions \
            bygul/tests/data/ATM-2FFMD73N3.fasta \
            bygul/tests/data/ARTIC_V4-1.bed \
            bygul/tests/data/reference.fasta \
            --simulator wgsim --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_max_mismatch(self):
        os.system(
            "bygul simulate-proportions \
            bygul/tests/data/ATM-2FFMD73N3.fasta \
            bygul/tests/data/ARTIC_V4-1.bed \
            bygul/tests/data/reference.fasta \
            --maxmismatch 2 --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_readcnt(self):
        os.system(
            "bygul simulate-proportions \
            bygul/tests/data/ATM-2FFMD73N3.fasta \
            bygul/tests/data/ARTIC_V4-1.bed --readcnt 200 \
            bygul/tests/data/reference.fasta --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_readlength(self):
        os.system(
            "bygul simulate-proportions \
            bygul/tests/data/ATM-2FFMD73N3.fasta \
            bygul/tests/data/ARTIC_V4-1.bed --readlength 130 \
            bygul/tests/data/reference.fasta --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_2genomes(self):
        os.system(
            "bygul simulate-proportions \
            bygul/tests/data/ATM-2FFMD73N3.fasta,\
            bygul/tests/data/KR-SEARCH-120354.fasta \
            bygul/tests/data/ARTIC_V4-1.bed \
            bygul/tests/data/reference.fasta --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))

    def test_simulation_with_2genomes_proportions(self):
        os.system(
            "bygul simulate-proportions \
            bygul/tests/data/ATM-2FFMD73N3.fasta,\
            bygul/tests/data/KR-SEARCH-120354.fasta \
            bygul/tests/data/ARTIC_V4-1.bed \
            bygul/tests/data/reference.fasta \
            --proportions 0.8,0.2 --redo"
        )
        self.assertTrue(file_exists(".", "results/reads_1.fastq"))


if __name__ == "__main__":
    unittest.main()
