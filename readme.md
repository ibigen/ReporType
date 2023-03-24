# Detect_type - RAW VERSION


Detect_type is a user-friendly [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline that identifies the variant of a given species (virus or bacteria) under study.

This pipeline uses [Abricate](https://github.com/tseemann/abricate) for variant identification and other independent software that prepares samples priviously analysed with illumina, nanopore and sanger technologies, to be analysed by Abricate.
The other software used are specified in the illustrative scheme.

Detect_type accepts as samples input .fasta, .fastq and .ab1 files, zipped or not, and will also require a database name or database fasta file input. 
The main output consists of a .csv file with the sample name, the genes found, coverage and identity percentage, the database used and the accession.
It will also produce the detailed Abricate output files and the intermidiate files that are produced by other software (trimmed samples, etc...).


![alt text](https://github.com/ibigen/loci_screening_typing/blob/main/images/Illustrative%20scheme.png.png)


## Instalation
You need to have  [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed.
All the other dependencies will be automatically installed with detect_type.
For installation, you need to:


1. Download this git repository:<br>
`$git clone https://github.com/ibigen/loci_screening_typing/`<br>
`$cd loci_screening_typing`

2. Install running:<br>
`$chmod +x install.sh`<br>
`$./install.sh`<br>


## Usage

First of all, you need to activate the detect_type environment:<br>
`$alias detect_type="conda activate detect_type && snakemake"`
`$conda activate detect_type`

Now open de "config.yaml" file and fill it with your options.

Execute the command line:
`detect_type --cores all `

Note that `--cores all` is need for a snakefile to run, but theres other Snakemake exucutable options  [here](https://snakemake.readthedocs.io/en/v5.1.4/executable.html).

<br>
<br>


When you are donne using detect_type you can deactivate the environment with:<br>
`conda deactivate detect_type`


## Uninstall
## Citation

