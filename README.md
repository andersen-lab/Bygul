# Bygul: Amplicon / Metagenomics Read Simulator

A tool for Amplicon/Metagenomics read simulation for waste water sequencing or other aplications. Users can easily simulate reads from mutiple samples with different proportions using the tool.

## Usage
If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <https://github.com/andersen-lab/Bygul> repository.

## Installation
Bygul is written in python 3 but it requires <a href="https://github.com/lh3/wgsim">wgsim</a> and 
<a href="https://github.com/seqan/seqan/blob/main/apps/mason2/README.mason_simulator">mason simulator</a>
to simulate reads.

### Local build from source
```
git clone https://github.com/andersen-lab/Bygul
cd Bygul
pip install -e .
```
Please note that pip does not install all the requirements,
some packages need to be installed via Conda or be built from source.

### Installing via Conda
```
conda create -n bygul bioconda::bygul
```

### Installing via Pypi
```
pip install bygul
```
Please note that some dependencies are not available through pypi.
You need to install them using conda or build from source.

# Amplicon sequencing mode
## Example commands

Run the tool using the following command.
 ```
bygul simulate-proportions [SAMPLE1.fasta,SAMPLE2.fasta,..] --primers [primer.bed] --reference [reference.fasta] --proportions [0.8,0.2,..] --outdir [output_directory]
 ```

Simulate reads from different samples without defining proportions (will be assigned randomly, proportions can be found in `results/sample_proportions.txt`) and allowing upto 2 SNPs mistmatches in the primer regions.
 ```
bygul simulate-proportions sample.fasta,sample2.fasta --primers primer.bed --reference reference.fasta --outdir results/ --maxmismatch 2
 ```
Simulate reads with user-defined proportions and specifing read simulator.
bygul uses wgsim as a simulator but you can change it to mason.
 ```
bygul simulate-proportions sample.fasta,sample2.fasta --primers primer.bed --reference reference.fasta --proportions 0.2,0.8 --simulator mason
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
Please remember that the primer file must contain a column containing primer sequence. The maximum number of mismatches allowed for each primer sequence is 1 SNP. To change this number, you may use the `--maxmismatches` flag.
#### Complete set of available parameters
To learn more about how to adjust other parameters for the simulator please read the documentation for wgsim and mason simulator. Users can pass any simulator parameter directly in their command. The only parameters set through bygul are `--readcnt` and `--wgsim_insert_size` for amplicon sequencing mode.
#### Simulated reads output
Simulated reads from all samples are located in `provided_output_path/reads.fastq`
#### Information about amplicon dropouts
In order to find more about amplicon dropouts, please refer to `provided_output_path/sample_name/amplicon_stats.csv` file. Please note that primer_seq_x and primer_seq_y define the left and right primer sequence whereas left_mismatch_map and right_mismatch_map shows the actual sequence found in the sample for a better comparison of mismatching bases in the primer sequence. Additionally, if there are any ambiguous bases present in the matching sequence, the ambiguous_bases value returns true. 

# Metagenomics mode
Users can now simulate reads from different samples in a metagenomics setting without specifying a primer bed file. Providing a reference sequence is not required for this setting.

## Example commands

Simulate reads from different samples without defining proportions (will be assigned randomly, proportions can be found in `results/sample_proportions.txt`).
 ```
bygul simulate-proportions sample.fasta,sample2.fasta --outdir results/ --simulation_mode metagenomics
 ```

Specify proportions for each sample and add other simulator specific parameters. To access simulator parameters, please read wgsim and mason documentation.
 ```
bygul simulate-proportions sample.fasta,sample2.fasta --proportions 0.5,0.5 --outdir results/ --simulation_mode metagenomics --simulator mason --illumina-read-length 200
 ```