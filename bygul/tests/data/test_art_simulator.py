import os
import tempfile
import unittest
from unittest.mock import patch

from bygul.utils import (
    run_simulation_on_fasta,
    run_simulation_on_fasta_single_genome,
)


def write_art_outputs(command):
    output_prefix = command[command.index("-o") + 1]
    with open(f"{output_prefix}1.fq", "w") as read1:
        read1.write("@read/1\nACGT\n+\nIIII\n")
    with open(f"{output_prefix}2.fq", "w") as read2:
        read2.write("@read/2\nTGCA\n+\nIIII\n")


class ArtSimulatorTests(unittest.TestCase):
    def test_amplicon_art_command_and_merge(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = os.path.join(tmpdir, "amplicon.fasta")
            with open(fasta, "w") as handle:
                handle.write(">contig1\nACGTACGT\n>contig2\nTGCATGCA\n")

            commands = []

            def fake_run(command, **kwargs):
                commands.append(command)
                write_art_outputs(command)

            with patch("bygul.utils.subprocess.run", side_effect=fake_run):
                run_simulation_on_fasta(
                    fasta,
                    tmpdir,
                    10,
                    "art",
                    150,
                    150,
                    0.0001,
                    125,
                    300,
                    25,
                    "HS25",
                )

            self.assertEqual(len(commands), 2)
            self.assertEqual(commands[0][0], "art_illumina")
            self.assertNotIn("-amp", commands[0])
            self.assertIn("-p", commands[0])
            self.assertIn("-na", commands[0])
            self.assertEqual(commands[0][commands[0].index("-c") + 1], "5")
            self.assertEqual(commands[0][commands[0].index("-l") + 1], "125")
            self.assertEqual(commands[0][commands[0].index("-m") + 1], "300")
            self.assertEqual(commands[0][commands[0].index("-s") + 1], "25")
            self.assertTrue(
                os.path.exists(os.path.join(tmpdir, "merged_reads_1.fastq"))
            )
            self.assertTrue(
                os.path.exists(os.path.join(tmpdir, "merged_reads_2.fastq"))
            )

    def test_metagenomics_art_outputs_are_fastq(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = os.path.join(tmpdir, "genome.fasta")
            with open(fasta, "w") as handle:
                handle.write(">genome\nACGTACGT\n")

            with patch(
                "bygul.utils.subprocess.run",
                side_effect=lambda command, **kwargs: write_art_outputs(
                    command
                ),
            ):
                run_simulation_on_fasta_single_genome(
                    fasta,
                    tmpdir,
                    12,
                    "art",
                    150,
                    200,
                    10,
                    "HS25",
                )

            self.assertTrue(os.path.exists(os.path.join(tmpdir,
                                                        "reads_1.fastq")))
            self.assertTrue(os.path.exists(os.path.join(tmpdir,
                                                        "reads_2.fastq")))


if __name__ == "__main__":
    unittest.main()
