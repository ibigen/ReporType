# Detect_type - RAW VERSION


Detect_type is a user-friendly [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline that identifies the variant of a given species (virus or bacteria) under study.

This pipeline uses [Abricate](https://github.com/tseemann/abricate) for variant identification and other independent software that prepares samples priviously analysed with illumina, nanopore and sanger technologies, to be analysed by Abricate.
The other software used are specified in the illustrative scheme.

Detect_type accepts as samples input .fasta, .fastq and .ab1 files, zipped or not, and will also require a database name or database fasta file input. 
The main output consists of a .csv file with the sample name, the genes found, coverage and identity percentage, the database used and the accession.
It will also produce the detailed Abricate output files and the intermidiate files that are produced by other software (trimmed samples, etc...).


![alt text](https://github.com/ibigen/loci_screening_typing/blob/main/detect_type_work_flow.png)


## Instalation
You need to have  [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed.
All the other dependencies will be automatically installed with detect_type.
For installation, you need to:


1. Download this git repository:<br>
`$ git clone https://github.com/ibigen/loci_screening_typing/`<br>
`$ cd loci_screening_typing`

2. Install running:<br>
`$ chmod +x install.sh`<br>
`$ ./install.sh`<br>


## Usage

First of all, you need to activate the detect_type environment with the commands:<br>
`$ alias detect_type="snakemake"`<br>
`$ conda activate detect_type`<br>

Now you must configure your entery params. You have to options, you can open de "config.yaml" file and fill it with your options you configurate them through the command line.<br>

There are some mandatory params for configuration listed below. <br>

**Database input params:** <br>

If you have already created a database or have a formatated fasta file with your database:<br>
> **database**: name of the database you wish to use. If is the first time using this database you need to add the path to the fasta file contaning the database, a new database will be created with the name of the given fasta file (example: database=path/to/my_database.fasta or database=my_database).<br>

If you don't have a database file already formatated for abricate, you can provide two files and the final location of your database fasta file:<br>
> **fasta_db**: fasta file with the sequences for your database (example: fasta_db=path/to/sequences.fasta).<br>
> **table_db**:  table (tsv) with sequece, id and acession for each sequence (example table_db=path/to/table.tsv)
> **final_fasta_db**: path to your final database fasta file, new database will be created with the name of the given fasta file (example: final_fasta_bd=path/to/my_database.fasta) <br>


**Samples params:** <br>

> **sample_directory**: path to the folder with the samples you wish to analyse. This folder can contain samples from different technologies, as long as they are all analyzed according to the same database (example: sample_directory=path/to/my_samples_folder/).<br>

Detect_type optional configuration params includes: <br>

> **output_name**: name of your final csv output file (default: output_name=all_samples)<br>
> **output_directory**: directory for your results (default: output_directory=results/)<br>
> **tecnology**: especify the tecnologies you are going to analyse. If you leave it with the default, all samples of the given folder will be analysed. Your opcions are: fasta,nanopore,illumina_single,illumina_paired,sanger, or any. You must separete them with a coma (default: tecnology=any)<br>
> **multi_fasta**: if you are going to analyse any multi-fasta files, give the name of each multi-fasta file. You can chosse "all" if all of your fasta files are multi-fasta(default: multi_fasta=none).<br>
> **threads**: threads you which to use (default: threads=2).<br>

You can also specify some software params.<br>

**Abricate params:**<br>

> **minid**: minimum DNA %identity (default: minid=1).<br>
> **mincov**: minimum DNA %coverage (default: mincov=1).<br>

**Illumina params:** (for single and paired reads)<br>
> **illuminaclip_single** and **illuminaclip_paired**: Trimmomatic Illuminaclip, directory of your illumina adapters, as well as specific cleaning informations for your file (default: illuminaclip=ILLUMINACLIP:primers/adapters.fasta:3:30:10:6:true)<br>
> **slidingwindow_single** and **slidingwindow_paired**: Trimmomatic Slidingwindow, minimum average quality established for each sequence according to a certain number of bases (default: slidingwindow=SLIDINGWINDOW:4:15).<br>
> **minlen_single** and **minlen_paired**: Trimmomatic Minlen, minimum read size (default: minlen=MINLEN:36).<br>
> **leading_single** and **leading_paired**: Trimmomatic Leading, bases to remove at the beginning of the read (default: leading=LEADING:3)<br>
> **trailing_single** and **trailing_paired**: Trimmomatic Trailing, bases to remove at the end of the read (default: trailing=TRAILING:3)<br>

**Nanopore params:**<br>
> **quality**: Nanofilt minimum quality mean for read (default: quality=8).<br>
> **length**: Nanofilt minimun length per read (default: length=50).<br>
> **maxlength**: Nanofilt maximum length per read (default: maxlength=50000)<br>
> **headcrop**: Nanofilt headcrop, bases to remove at the beginning of the read (default: headcrop=30).<br>
> **tailcrop**: Nanofilt tailcrop, bases to remove at the end of the read (default: tailcrop=30).<br>
> **kmer**: Raven k-mer, length of minimizers used to find overlaps (default: kmer=15). <br>
> **polishing**: Raven polishing-rounds, number of times racon is invoked (default: polishing=2).<br>

**Sanger params:**
> **startbase**: Abiview first sequence base to report or display (default: startbase=20).<br>
> **endbase**: Abiview last sequence base to report or display (default: endbase=800). <br>


The optional configuration params also include all the configuration params for Snakemake, that can be consulted [here](https://snakemake.readthedocs.io/en/v5.1.4/executable.html). The most relevant Snakemake executable params are: <br>
> **--cores**: number of CPU to be used, it is mandatory (example: --cores all).<br>
> **-np**: dry-run to verify the jobs you are submiting. <br>
> **--config**: you must use this command before configurate the detect_type params previosly refered (example: --config database=my_database).<br>
> **--snakefile**: you can execute detect_type in any directory using this command to specify the directory for the snakefile of detect_type (example: --snakefile path/to/loci_screening_typing/snakefile).<br> 
> **--configfile**: you can execute detect_type in any directory using this command to specify the directory for the config file of detect_type (example: --configfile path/to/loci_screening_typing/config.yaml).<br>


### Command line examples:

If you configurate the config.yaml file, you can only run:<br>
`$ detect_type --cores all `<br>

If you wish to configurate throug the command line, here are some examples:<br>

>To run out of loci_screening_typing directory:<br>
`$ detect_type --cores all --snakefile path/to/loci_screening_typing/snakefile  --configfile path/to/loci_screening_typing/config.yaml --config database=my_database sample_directory=path/to/my_samples_folder/`<br>
>To run in loci_screening_typing directory:<br>
`$ detect_type --cores all --config database=my_database sample_directory=path/to/my_samples_folder/`<br>
>Now with some opcional configuration params:<br>  
`$ detect_type --cores all --config database=my_database sample_directory=path/to/my_samples_folder/ output_name=all_samples_database tecnology=fasta,illumina_single multi_fasta=multi_fasta_file_1,multi_fasta_file_2 `<br>


<br>
<br>


When you are donne using detect_type you can deactivate the environment with:<br>
`$ conda deactivate detect_type`<br>


## Uninstall

To uninstall detect_type, you need to delete the conda environment with: <br>
`$ conda env remove --name detect_type`<br>

## Citation

