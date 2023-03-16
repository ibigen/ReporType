# Detect_type - RAW VERSION


Detect_type is a user-friendly [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline that identifies the variant of a given species (virus or bacteria) under study.

This pipeline uses [Abricate](https://github.com/tseemann/abricate) for variant identification and other independent software that prepares samples priviously analysed with illumina, nanopore and sanger technologies, to be analysed by Abricate.
The other software used are specified in the illustrative scheme.

Detect_type accepts as samples input .fasta, .fastq and .ab1 files, zipped or not, and will also require a database name or database fasta file input. 
The main output consists of a .csv file with the sample name, the genes found, coverage and identity percentage, the database used and the accession.
It will also produce the detailed Abricate output files and the intermidiate files that are produced by other software (trimmed samples, etc...).


![alt text](https://github.com/ibigen/loci_screening_typing/blob/main/detect_type/images/Illustrative%20scheme.png.png)


## Instalation
You need to have  [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed.
All the other dependencies will be automatically installed with detect_type.
For installation, you need to:


1. Download this git repository:<br>
`git clone https://github.com/ibigen/loci_screening_typing/`

2. Install running:<br>
`chmod +x install.sh`<br>
`./install.sh`


## Usage

First of all, you need to activate the detect_type environment:<br>
`conda activate detect_type`


Because detect_type is based on a Snakemake pipeline, all the commands available for Sankemake are also 
available for detec_type. You can check the Snakemake original commands [here](https://snakemake.readthedocs.io/en/v5.1.4/executable.html).

In addition to all the options that Snakemake offers, there are some that are unique to detect_type:

   - **-d, --database** &rarr;  (Mandatory) Insert the name of the database for the samples you want to analyse (standard: influenza). If is the first time using this databse you need to add the path to the fasta file contining the database, a new database will be created with the name of the given fasta file (standard: /path_to_database_file/influenza.fasta)


   - **-i, --input**  &rarr; (Mandatory) Insert the path to the folder with the samples you wish to analyse (standard:path_to/my_samples). This folder can contain samples from different technologies, as long as they are all analyzed according to the same database.


   - **-o, --output**  &rarr; Chose the name of your final csv output file (default="all_samples").
 

   - **-t, --tecnology**  &rarr;  Opcionally, especify the tecnologies you are going to analyse (default="any"). If you leave it with the default, all samples of the given folder will be analysed. Choices: 'fasta', 'nanopore', 'illumina_single', 'illumina_paired', 'sanger' or 'any'.


   - **-th, --threads** &rarr; Choose how many threads you which to use (default=2).

**Abricate params**

   - **-minid, --abricate_minid** &rarr; Abricate params: Minimum DNA %identity (default=1).


   - **-mincov, --abricate_mincov** &rarr; Abricate params: Minimum DNA %coverage (default=1).

**Sanger params**

   - **-startbase, --abiview_startbase** &rarr; Abiview params: First base to report or display (default=20)


   - **-endbase, --abiview_endbase** &rarr; Abiview params: Last sequence base to report or display (default=800).


### Command line

A simple exemple of a executable command line woud be:<br>
`detect_type --cores all -d /path_to_database_file/influenza.fasta -i path_to/my_samples`

Note that `--cores all` is need for a snakefile to run, but theres other Snakemake exucutable options  [here](https://snakemake.readthedocs.io/en/v5.1.4/executable.html).

<br>
<br>


When you are donne using detect_type you can deactivate the environment with:<br>
`conda deactivate detect_type`


## Uninstall
## Citation

