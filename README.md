# Bygul: Amplicon & Metagenomics Read Simulator

**Bygul** is a Python 3 tool designed for simulating sequencing reads in wastewater surveillance and other metagenomic applications. It allows users to simulate complex multi-sample datasets with customizable proportions using industry-standard backends like `wgsim`, `mason`, and ART Illumina.

---

## 🏗 Installation

Bygul requires **Python 3**. Since it relies on external simulators (`wgsim`, `mason`, and ART), we recommend using Conda to manage dependencies.For more info on <a href="https://github.com/lh3/wgsim">wgsim</a>, <a href="https://github.com/seqan/seqan/blob/main/apps/mason2/README.mason_simulator">mason simulator</a>, and <a href="https://www.niehs.nih.gov/research/resources/software/biostatistics/art">ART</a> please check their documentations.

### Option 1: Via Conda (Recommended)
```bash
conda create -n bygul bioconda::bygul
```

### Option 2: Via PyPI
```bash
pip install bygul
```
*Note: Some binary dependencies (wgsim/mason/ART) may need to be installed manually or built from source if using this method.*

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
bygul simulate-proportions [SAMPLE1.fasta,SAMPLE2.fasta] --primers [primer.bed] --reference [reference.fasta] --proportions [0.8,0.2] --outdir [output_dir]
```

### Advanced Examples
* **Random Proportions & Mismatches:**
    Simulate with random proportions and allow up to 2 SNPs in primer regions.
    ```bash
    bygul simulate-proportions sample1.fasta,sample2.fasta --primers primer.bed --reference reference.fasta --outdir results/ --maxmismatch 2
    ```
* **Switching Simulators:**
    Use `mason` or `art` instead of the default `wgsim`.
    ```bash
    bygul simulate-proportions sample1.fasta,sample2.fasta --primers primer.bed --reference reference.fasta --simulator mason
    bygul simulate-proportions sample1.fasta,sample2.fasta --primers primer.bed --reference reference.fasta --simulator art
    ```
* **Custom Error Rates & Lengths:**
    Pass simulator-specific parameters (e.g. indel fraction `-R`) directly.
    ```bash
    bygul simulate-proportions sample1.fasta,sample2.fasta --primers primer.bed --reference reference.fasta -R 0.01
    ```
---

## 🌍 Usage: Metagenomics Mode
Simulate reads from entire samples without requiring a primer BED file or a reference sequence.

### Basic Metagenomics Simulation
```bash
bygul simulate-proportions sample1.fasta,sample2.fasta --outdir results/ --simulation_mode metagenomics
```

### Metagenomics with Specific Parameters
```bash
bygul simulate-proportions sample1.fasta,sample2.fasta --proportions 0.5,0.5 --outdir results/ --simulation_mode metagenomics --simulator mason --illumina-read-length 200
```

---

## 📝 Technical Notes

### Parameter Handling
Bygul acts as a wrapper. While most flags are passed directly to the underlying simulators, the following are managed directly by Bygul for more realistic simulations(amplicon simulation mode only):
- `--readcnt`: Number of reads per amplicon.
- `--wgsim_insert_size`: Insert size for wgsim.
- `--wgsim_read_length` / `--wgsim_error_rate`.
- `--art_read_length` / `--art_insert_size` / `--art_insert_sd` / `--art_seq_system`: paired-end ART Illumina defaults.

To see all available backend flags, run:
```bash
wgsim --help
mason_simulator --help
art_illumina --help
```
Please note that some dependencies are not available through pypi.
You need to install them using conda or build from source.

# Amplicon sequencing mode
## Example commands

Run the tool using the following command.

If you are just checking your primers without wanting to run the simulation, you can provide a multi-fasta file including all the sequences along with the primer file and the reference used to generate the primer file. This will create amplicon_stats.csv file as described below.

 ```
bygul check-primers sequences.fasta primer.bed reference.fasta
 ```

Example simualtion command
 ```
bygul simulate-proportions [SAMPLE1.fasta,SAMPLE2.fasta,..] --primers [primer.bed] --reference [reference.fasta] --proportions [0.8,0.2,..] --outdir [output_directory]
 ```

Simulate reads from different samples without defining proportions (will be assigned randomly, proportions can be found in `results/sample_proportions.txt`) and allowing upto 2 SNPs mistmatches in the primer regions.
 ```
bygul simulate-proportions sample.fasta,sample2.fasta --primers primer.bed --reference reference.fasta --outdir results/ --maxmismatch 2
 ```
Simulate reads with user-defined proportions and specifing read simulator.
bygul uses wgsim as a simulator but you can change it to mason or ART.
 ```
bygul simulate-proportions sample.fasta,sample2.fasta --primers primer.bed --reference reference.fasta --proportions 0.2,0.8 --simulator mason
bygul simulate-proportions sample.fasta,sample2.fasta --primers primer.bed --reference reference.fasta --proportions 0.2,0.8 --simulator art --art_read_length 150
 ```
Simulate reads with user-defined proportions and number of reads per amplicon.
 ```
bygul simulate-proportions sample.fasta,sample2.fasta --primers primer.bed --reference reference.fasta --proportions 0.2,0.8 --readcnt 1000
 ```

Simulate reads with additional parameters such as base error rate, read length and indels fraction
 ```
bygul simulate-proportions sample.fasta,sample2.fasta --primers primer.bed --reference reference.fasta --proportions 0.2,0.8 --readcnt 1000 -e 0.001 -1 400 -2 400 -R 0.01
 ```
## Notes
#### Number of reads per amplicon
It is recommended to define the number of reads per amplicon to be greater than the number of contigs in your amplicon file. This is particularly important when your primers are designed for whole genome sequencing, where each amplicon may contain a substantial number of contigs. Setting too few reads per amplicon may result in empty read files for certain amplicons, leading to incomplete simulated reads.
#### Primer bed file
### 🧬 Input BED File Format
The pipeline expects a tab-delimited BED file where the first six columns represent standard genomic coordinates (`chrom`, `chromStart`, `chromEnd`, `name`, `poolName`, `strand`). Crucially, the fourth column (`name`) must follow a strict naming convention to prevent downstream parsing failures in variant-calling tools: **`[Scheme-Name]_[AmpliconNumber]_[Direction]_[OptionalSuffix]`** (e.g., `SARS-CoV-2_3_LEFT` or `SARS-CoV-2_3_LEFT_alt`). To ensure structural boundaries are parsed correctly, the prefix must not contain underscores, and any optional trailing modifiers must be restricted to standard alternative tags (`_alt`, `_ALT1`) or tracking indexes (`_0`, `_1`). Multi-level pool formatting, such as `SARS-CoV-2_400_1_LEFT_1`, is malformed and will fail validation.. The maximum number of mismatches allowed for each primer sequence is 1 SNP. To change this number, you may use the `--maxmismatches` flag.
#### Complete set of available parameters
To learn more about how to adjust other parameters for the simulator please read the documentation for wgsim, mason simulator, and ART. Users can pass any simulator parameter directly in their command. The only parameters set through bygul are `--readcnt`, `--wgsim_insert_size`, `--wgsim_read_length`, `--wgsim_error_rate`, `--art_read_length`, `--art_insert_size`, `--art_insert_sd`, and `--art_seq_system`.
#### Simulated reads output
Simulated reads from all samples are located in `provided_output_path/reads.fastq`
#### Information about amplicon dropouts
In order to find more about amplicon dropouts, please refer to `provided_output_path/sample_name/amplicon_stats.csv` file. Please note that primer_seq_x and primer_seq_y define the left and right primer sequence whereas left_match and right_match shows the actual sequence found in the sample for a better comparison of mismatching bases in the primer sequence. Additionally, if there are any ambiguous bases present in the matching sequence, the ambiguous_bases value returns true. 

# Metagenomics mode
Users can now simulate reads from different samples in a metagenomics setting without specifying a primer bed file. Providing a reference sequence is not required for this setting.

## Example commands

Simulate reads from different samples without defining proportions (will be assigned randomly, proportions can be found in `results/sample_proportions.txt`).
 ```
bygul simulate-proportions sample.fasta,sample2.fasta --outdir results/ --simulation_mode metagenomics
 ```

Specify proportions for each sample and add other simulator specific parameters. To access simulator parameters, please read wgsim, mason, and ART documentation.
 ```
bygul simulate-proportions sample.fasta,sample2.fasta --proportions 0.5,0.5 --outdir results/ --simulation_mode metagenomics --simulator mason --illumina-read-length 200
bygul simulate-proportions sample.fasta,sample2.fasta --proportions 0.5,0.5 --outdir results/ --simulation_mode metagenomics --simulator art --art_read_length 150
 ```
