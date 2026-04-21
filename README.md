# Bygul: Amplicon & Metagenomics Read Simulator

**Bygul** is a Python 3 tool designed for simulating sequencing reads in wastewater surveillance and other metagenomic applications. It allows users to simulate complex multi-sample datasets with customizable proportions using industry-standard backends like `wgsim` and `mason`.

---

## 🏗 Installation

Bygul requires **Python 3**. Since it relies on external simulators (`wgsim` and `mason`), we recommend using Conda to manage dependencies.

### Option 1: Via Conda (Recommended)
```bash
conda create -n bygul bioconda::bygul
```

### Option 2: Via PyPI
```bash
pip install bygul
```
*Note: Some binary dependencies (wgsim/mason) may need to be installed manually or built from source if using this method.*

### Option 3: Local Build from Source
```bash
git clone [https://github.com/andersen-lab/Bygul](https://github.com/andersen-lab/Bygul)
cd Bygul
pip install -e .
```

---

## 🧬 Usage: Amplicon Sequencing Mode
Use this mode when simulating specific genomic regions defined by a primer set.

### Basic Command
```bash
bygul simulate-proportions [SAMPLES.fasta] --primers [primer.bed] --reference [ref.fasta] --proportions [0.8,0.2] --outdir [output_dir]
```

### Advanced Examples
* **Random Proportions & Mismatches:**
    Simulate with random proportions and allow up to 2 SNPs in primer regions.
    ```bash
    bygul simulate-proportions s1.fasta,s2.fasta --primers p.bed --reference r.fasta --outdir results/ --maxmismatch 2
    ```
* **Switching Simulators:**
    Use `mason` instead of the default `wgsim`.
    ```bash
    bygul simulate-proportions s1.fasta,s2.fasta --primers p.bed --reference r.fasta --simulator mason
    ```
* **Custom Error Rates & Lengths:**
    Pass simulator-specific parameters (e.g., error rate `-e`, read lengths `-1` and `-2`, indel fraction `-R`) directly.
    ```bash
    bygul simulate-proportions s1.fasta,s2.fasta --primers p.bed --reference r.fasta -e 0.001 -1 400 -2 400 -R 0.01
    ```

---

## 🌍 Usage: Metagenomics Mode
Simulate reads from entire samples without requiring a primer BED file or a reference sequence.

### Basic Metagenomics Simulation
```bash
bygul simulate-proportions s1.fasta,s2.fasta --outdir results/ --simulation_mode metagenomics
```

### Metagenomics with Specific Parameters
```bash
bygul simulate-proportions s1.fasta,s2.fasta --proportions 0.5,0.5 --outdir results/ --simulation_mode metagenomics --simulator mason --illumina-read-length 200
```

---

## 📝 Technical Notes

### Parameter Handling
Bygul acts as a wrapper. While most flags are passed directly to the underlying simulators, the following are managed directly by Bygul for better realism:
- `--readcnt`: Number of reads per amplicon.
- `--wgsim_insert_size`: Insert size for wgsim.
- `--wgsim_read_length` / `--wgsim_error_rate`.

To see all available backend flags, run:
```bash
wgsim --help
mason_simulator --help
```

### Best Practices
* **Read Counts:** Set `--readcnt` higher than the number of contigs in your amplicon file. Too few reads can result in empty files for certain amplicons.
* **Primer Files:** The BED file **must** include a column with the primer sequence. Bygul allows 1 SNP mismatch by default; use `--maxmismatch` to change this.

### Output Files
* **Consolidated Reads:** Simulated reads from all samples are at `outdir/reads.fastq`.
* **Proportions:** Assigned proportions are recorded in `results/sample_proportions.txt`.
* **Quality Metrics:** Check `outdir/[sample_name]/amplicon_stats.csv` for information on **amplicon dropouts**, mismatches, and ambiguous bases.

---

## 🎓 Citation
If you use this workflow in a paper, please cite the original repository:
[https://github.com/andersen-lab/Bygul](https://github.com/andersen-lab/Bygul)