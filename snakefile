import pandas as pd
import os
import re
from statistics import mean
import subprocess
from Bio import SeqIO
import gzip
import shutil
import tempfile
import argparse

######################  INPUTS  ###############################

# Define as opções de linha de comando
parser = argparse.ArgumentParser(description='Detect_type extra opcions (using snakemake)')

parser.add_argument('-d', '--database', required=True, type=str,
                     help='Insert the name of the database for the samples you want to analyse (standard: influenza).If is the first time using this databse you need to add the path to the fasta file contining the database, a new database will be created with the name of the given fasta file (standard: /path_to_database_file/influenza.fasta).')

parser.add_argument('-i', '--input', required=True, type=str,
                     help='Insert the path to the folder with the samples you wish to analyse (standard:path_to/my_samples). This folder can contain samples from different technologies, as long as they are all analyzed according to the same database.')

parser.add_argument('-o', '--output', default="all_samples", type=str,
                     help='Chose the name of your final csv output file (default="all_samples").')

parser.add_argument('-t', '--tecnology', nargs='+', type=str, default=['any'], choices=['fasta','nanopore','illumina_single','illumina_paired','sanger','any'],
                     help='Opcionally, especify the tecnologies you are going to analyse (default="any"). If you leave it with the default, all samples of the given folder will be analysed.')

parser.add_argument('-th', '--threads', default=2, type=int,
                     help='Choose how many threads you which to use (default=2).')

parser.add_argument('-minid', '--abricate_minid', default=1, type=int,
                     help='Abricate params: Minimum DNA %identity (default=1).')

parser.add_argument('-mincov', '--abricate_mincov', default=1, type=int,
                     help='Abricate params: Minimum DNA %coverage (default=1).')

parser.add_argument('-startbase', '--abiview_startbase', default=20, type=int,
                     help='Abiview params: First base to report or display (default=20).')

parser.add_argument('-endbase', '--abiview_endbase', default=800, type=int,
                     help='Abiview params: Last sequence base to report or display (default=800).')


args, _ = parser.parse_known_args()

db = args.database
sample_path = args.input
table = "table_configuration.py"   #########CHECK!!!!!!!!!!!!!!!
threads = args.threads
minid=args.abricate_minid
mincov=args.abricate_mincov
output=args.output
startbase=args.abiview_startbase
endbase=args.abiview_endbase
tec_input=args.tecnology

###################### DB INPUT ###############################

def path_check(path):
    exist = os.path.exists(path)
    return (exist)


def db_check(db):
    #shell("abricate --setupdb")
    for i in shell("abricate --list", iterable=True):
        if db in i:
            check = True
            break
        else:
            check = False
    return (check)


def fasta_import_to_abricate(path):
    path = path.replace('"', "")
    path = str(path)
    return (path)


def abri_default_db():
    x = subprocess.check_output("command abricate --help | grep datadir", shell=True)
    x = str(x)
    x = x.strip()
    x = x.split("]")
    x = x[-2]
    x=x.split("[")
    path=x[-1]
    return (path)


def name_db(path):
    name = path.strip()
    name = name.split("/")
    name = name[-1]
    name = name.split(".")
    name = name[0]
    return (name)


fasta_check = path_check(db)

if fasta_check == False:
    db_check = db_check(db)
    if db_check == False:
        print("Input error: The database name is incorrect or the given path to fasta file does not exist!")
    else:
        print('Database available! Starting analysis...')
else:
    fasta_path = fasta_import_to_abricate(db)
    join = ['"', '"']
    fasta_path = fasta_path.join(join)
    db = name_db(fasta_path)
    check = db_check(db)
    if check == False:
        abricate_db = abri_default_db()
        shell("mkdir {abricate_db}/{db}")
        shell("cp {fasta_path} {abricate_db}/{db}/sequences")
        shell("cd {abricate_db}/{db} ; makeblastdb -in sequences -title {db} -dbtype nucl -hash_index")
        check=db_check(db)
        if check == False:
            print("Input error: Failed creating new database! Check the fasta file and respective path.")
        else:
            print('Database available! Starting analysis...')
    else:
        print('Database available! Starting analysis...')




############## SAMPLES INPUT ##############################


FASTA_EXTENSION = [".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fna.gz"]
FASTQ_EXTENSION = [".fastq", ".fastq.gz", ".fq", ".fq.gz"]
AB1_EXTENSION = [".ab1"]

### size limit of the sequences
MAX_LENGHT_ILLUMINA_FASQC_SEQ = 500     ## because of IONtorrent
MIN_LENGHT_MINION_FASQC_SEQ = 100

################ SEPARATION #######################

def __test_ends_with(file_name, list_ends):
    """ check if the file as this etention """
    for end_str in list_ends:
        if file_name.endswith(end_str): return True
    return False

def __is_nanopore(file_name):
    """
    return True for nanopore, default illumina or IONtorrent
    raise exception if can not detected
    """

    vect_length = []
    with (gzip.open(file_name, mode='rt') if file_name.endswith(".gz") \
        else open(file_name, mode='r')) as handle_read:
        ### read 100 lines
        try:
            count = 0
            for record in SeqIO.parse(handle_read, "fastq"):
                vect_length.append(len(str(record.seq)))
                count += 1
                if (count > 100): break
                #print("mean(vect_length): {} ".format(mean(vect_length)))
        except:
            pass

    ### if read something in last SeqIO.parse
    if (len(vect_length) > 1):
        if (mean(sorted(vect_length, reverse=True)[:5]) <= MAX_LENGHT_ILLUMINA_FASQC_SEQ): return False
        if (mean(vect_length) > MIN_LENGHT_MINION_FASQC_SEQ): return True
        raise Exception("Can not detect file format. Ensure Illumina fastq file.")
    raise Exception("File is not in fastq.gz format.")


def __is_fasta_file(file_name):
    """ test fasta format """
    with (gzip.open(file_name, mode='rt') if file_name.endswith(".gz") \
        else open(file_name, mode='r')) as handle_read:
        try:
            for record in SeqIO.parse(handle_read, "fasta"):
                handle_read.close()
                return True
        except:
            pass
    return False

def remove_extensions_file_name(file_name, vect_to_remove=FASTQ_EXTENSION):
    if os.path.exists (file_name):
        file_name=file_name.strip()
        file_name = file_name.split("/")
        file_name = file_name[-1]
    for to_search in vect_to_remove:
        if (file_name.endswith(to_search)):
            return file_name[:len(file_name) - len(to_search)]
    return file_name


def __get_prefix_file_name(file_name):
    """ returnprefix file name based on patterns """
    m = re.search('[a-zA-Z0-9_\.]+(_[lL]\d+_[rR]\d+_\d+)_[a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+(_[lL]\d+_[rR]\d+)[a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+(_[rR]\d+_[lL]\d+)[a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+(_[lL]\d+)[a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+(_[rR]\d+)[a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+([_]reverse)[a-zA-Z0-9_\.]+', file_name.lower())
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+([_]forward)[a-zA-Z0-9_\.]+', file_name.lower())
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+(_\d+)[\.][a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+(_\d+)[_\.][a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+(_\d+)[_][a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    return remove_extensions_file_name(file_name)


def get_number_file(file_name):
    """ return the number of file
    :out None if not found  """

    m = re.search('[a-zA-Z0-9_\.]+(_[rR]\d+)[a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[m.regs[1][0]:m.regs[1][1]].lower().replace('r', '').replace('_', '')

    m = re.search('[a-zA-Z0-9_\.]+(_\d+)[_\.][a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[m.regs[1][0]:m.regs[1][1]].lower().replace('_', '')

    m = re.search('[a-zA-Z0-9_\.]+([_]reverse)[a-zA-Z0-9_\.]+', file_name.lower())
    if (not m is None): return file_name[m.regs[1][0]:m.regs[1][1]].lower().replace('_', '').replace('reverse', '2')

    m = re.search('[a-zA-Z0-9_\.]+([_]forward)[a-zA-Z0-9_\.]+', file_name.lower())
    if (not m is None): return file_name[m.regs[1][0]:m.regs[1][1]].lower().replace('_', '').replace('forward', '1')
    return None


def get_order_pair(b_test_first_in_pair, vect_names):
    """ from a list retruen the first that belongs to first or second pair """
    for path_and_file in vect_names:
        file_name = os.path.basename(path_and_file)
        if b_test_first_in_pair and get_number_file(file_name) == "1": return path_and_file
        if not b_test_first_in_pair and get_number_file(file_name) == "2": return path_and_file
    return None


def collect_files(*args):
    ### list all possible files
    list_fasta = []
    list_fastq_illumina_single = []
    ## paired end files { 'r1' : [ 'file_1.fastq.gz', 'file_r1.fastq.gz', .. ]
    ##                    'r2' : [ 'file_2.fastq.gz', 'file_r2.fastq.gz', .. ] }
    dt_fastq_illumina_pair = { 'r1' : [], 'r2' : [] }
    list_fastq_nanopore = []
    list_ab1 = []

    ## list all files names
    dt_list_all_files_names = {}

    for item in args:
        for p, _, files in os.walk(os.path.abspath(item)):
            for file in files:

                ## file repeated
                path_file = os.path.join(p, file)
                if file in dt_list_all_files_names:
                    print("File name repeated: " + file)
                    continue

                ## testing file
                if __test_ends_with(path_file, FASTA_EXTENSION):
                    if __is_fasta_file(path_file):
                        dt_list_all_files_names[file] = 1
                        list_fasta.append(path_file)
                    else:
                        print("Error: some problem with this fasta file: " + file)
                    continue
                if __test_ends_with(path_file, FASTQ_EXTENSION):
                    try:
                        if __is_nanopore(path_file):
                            list_fastq_nanopore.append(path_file)
                        else: ### first, all of them are single end
                            list_fastq_illumina_single.append(path_file)
                        dt_list_all_files_names[file] = 1
                    except Exception as e:
                        print("Error: " + str(e) + "\nFile: " + path_file)
                    continue
                if __test_ends_with(file, AB1_EXTENSION):
                    dt_list_all_files_names[file] = 1
                    list_ab1.append(path_file)
                    continue

                ## didn't found any technology
                print("File not recognize: " + path_file)

    ### try to find close
    dt_fastq_illumina_pair_temp = {}
    for path_and_file_name in list_fastq_illumina_single:
        sample_name = __get_prefix_file_name(os.path.basename(path_and_file_name))

        if sample_name in dt_fastq_illumina_pair_temp: dt_fastq_illumina_pair_temp[sample_name].append(path_and_file_name)
        else: dt_fastq_illumina_pair_temp[sample_name] = [path_and_file_name]

    ### get the candidates and join them
    dt_to_remove = {}   ## files to remove because are in pair
    for sample_name in dt_fastq_illumina_pair_temp:
        first_in_pair = get_order_pair(True, dt_fastq_illumina_pair_temp[sample_name])
        second_in_pair = get_order_pair(False, dt_fastq_illumina_pair_temp[sample_name])

        if not first_in_pair is None and not second_in_pair is None and first_in_pair != second_in_pair:
            dt_to_remove[first_in_pair] = 1
            dt_to_remove[second_in_pair] = 1
            dt_fastq_illumina_pair['r1'].append(first_in_pair)
            dt_fastq_illumina_pair['r2'].append(second_in_pair)

    ###  remove the ones in the paired files
    for key in dt_to_remove:
        list_fastq_illumina_single.remove(key)

    ## test files
    if len(list_fasta) + len(list_fastq_illumina_single) +\
           (len(dt_fastq_illumina_pair['r1']) * 2) + len(list_fastq_nanopore) == 0:
        sys.exit("Error, no file founded at '{}' path".format(args[0]))

    return list_fasta, list_fastq_illumina_single, dt_fastq_illumina_pair,\
           list_fastq_nanopore, list_ab1




########################PRE ANALISIS###############################

## test exists
exist = path_check(sample_path)

if exist == False:
    print("Input error: The sample path given is not correct!")
else:
    ## test if windows or not and add the back slash to the end
    if os.name == 'nt':
        sample_path += '\'
    elif not sample_path.endswith('/'):
        sample_path += '/'

    list_fasta, list_illumina_fastq_single, dt_fastq_illumina_pair,\
        list_fastq_nanopore, list_ab1 = collect_files(sample_path)


### get all samples name, till the end
SAMPLES_NAME_global = []


##############TECNOLOGY CONFIRMATION#################

### get all the detected technologies
tec_app=[]
if len(dt_fastq_illumina_pair['r1']) > 0:
    tec_app.append("illumina_paired")

if len (list_illumina_fastq_single) > 0:
    tec_app.append("illumina_single")

if len (list_fasta)>0:
    tec_app.append("fasta")

if len (list_fastq_nanopore)>0:
    tec_app.append("nanopore")

if len (list_ab1)>0:
    tec_app.append("sanger")

### check useres input
for tecnology_input in tec_input:
    if tecnology_input == "any":
        tec_input.remove("any")

###remove non listed tecnoloogies

if len(tec_input)>0:
    if "fasta" not in tec_input:
        list_fasta=[]
    if "nanopore" not in tec_input:
        list_fastq_nanopore=[]
    if "sanger" not in tec_input:
        list_ab1=[]
    if "illumina_single" not in tec_input:
        list_illumina_fastq_single=[]
    if "illumina_paired" not in tec_input:
        dt_fastq_illumina_pair = { 'r1' : [], 'r2' : [] }


########################ANALISIS###############################
rule all:
    input:
        expand("results/{output}.csv",output=output)
    params:
        tec_input=tec_input,
        tec_app=tec_app
    run:
        for tecnology_input in params.tec_input:
            if tecnology_input != "any":
                if tecnology_input not in tec_app:
                    print("Warning: no files detect for", tecnology_input, "tecnology. Verify if you enter the rigth tecnology.")
        if len (params.tec_input) > 0:
            if params.tec_app not in params.tec_input:
                print("Warning: there are samples that are not", ', '.join(tec_input), "files, those samples will not be analysed. Verify if you enter the rigth tecnologies.")



if len(dt_fastq_illumina_pair['r1']) > 0:

    ## list of the samples
    SAMPLES_NAME_global.extend([__get_prefix_file_name(os.path.basename(file_name)) for file_name in dt_fastq_illumina_pair['r1']])

    dt_extention_r1 = { __get_prefix_file_name(os.path.basename(file_name)) : dt_fastq_illumina_pair['r1'][index] for index, file_name in enumerate(dt_fastq_illumina_pair['r1']) }
    dt_extention_r2 = { __get_prefix_file_name(os.path.basename(file_name)) : dt_fastq_illumina_pair['r2'][index] for index, file_name in enumerate(dt_fastq_illumina_pair['r2']) }

    rule pre_illumina_paired:
        input:
            r1 = lambda wildcards: dt_extention_r1[wildcards.sample],
            r2 = lambda wildcards: dt_extention_r2[wildcards.sample]
        output:
            s1="results/intermidiate/trimm_paired_sur_1/{sample}.fastq.gz",
            d1="results/intermidiate/trimm_paired_rem_1/{sample}.fastq.gz",
            s2="results/intermidiate/trimm_paired_sur_2/{sample}.fastq.gz",
            d2="results/intermidiate/trimm_paired_rem_2/{sample}.fastq.gz"
        params:
            threads = threads
        shell:
            """trimmomatic PE -threads {params.threads} {input.r1} {input.r2} {output.s1} {output.d1} {output.s2} {output.d2} ILLUMINACLIP:primers/adapters.fasta:3:30:10:6:true SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36 || (touch {output.s1} {output.s2} && echo Warning: trimmomatic failed, were created empty files)"""

    rule illumina_paired:
        input:
            r1 ="results/intermidiate/trimm_paired_sur_1/{sample}.fastq.gz",
            r2 ="results/intermidiate/trimm_paired_sur_2/{sample}.fastq.gz"
        output:
             "results/intermidiate/fasta_files/{sample}.fasta"
        params:
            threads = threads
        shell:
            """spades --pe1-1 {input.r1} --pe1-2 {input.r2} --only-assembler -t {params.threads} -o results/intermidiate/spades/{wildcards.sample} ||  echo Warning: spades failed, were created empty files;"""
            """cp results/intermidiate/spades/{wildcards.sample}/contigs.fasta {output} || touch {output}"""#; rm -rf results/intermidiate/spades/{wildcards.sample}"""


if len(list_illumina_fastq_single) > 0:

    ### add sample names for global
    SAMPLES_NAME_global.extend(__get_prefix_file_name(os.path.basename(file_name)) for file_name in list_illumina_fastq_single)

    dt_extention_single = { __get_prefix_file_name(os.path.basename(file_name)) : list_illumina_fastq_single[index] for index, file_name in enumerate(list_illumina_fastq_single) }

    rule pre_illumina_single:
        input:
            r1 = lambda wildcards: dt_extention_single[wildcards.sample]
        output:
            s1="results/intermidiate/trimm_single_sur_1/{sample}.fastq.gz"
        params:
            threads = threads
        shell:
            """trimmomatic SE -threads {params.threads} {input.r1} {output.s1} ILLUMINACLIP:primers/adapters.fasta:3:30:10:6:true SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36 || (touch {output.s1} && echo Warning: trimmomatic failed, were created empty files)"""

    rule illumina_single:
        input:
            r1 = "results/intermidiate/trimm_single_sur_1/{sample}.fastq.gz"
        output:
            "results/intermidiate/fasta_files/{sample}.fasta"
        params:
            threads = threads
        shell:
            """spades -s {input.r1} --only-assembler -t {params.threads} -o results/intermidiate/spades/{wildcards.sample} || echo Warning: spades failed, were created empty files;"""
            """cp results/intermidiate/spades/{wildcards.sample}/contigs.fasta {output} || touch {output}"""#; rm -rf results/intermidiate/spades/{wildcards.sample}"""


if len(list_fasta) > 0:

    ### sample name for fasta

    SAMPLES_NAME_global.extend(remove_extensions_file_name(file_name,FASTA_EXTENSION) for file_name in list_fasta)

    dt_extention_fasta_file = {remove_extensions_file_name(os.path.basename(file_name),FASTA_EXTENSION) : list_fasta[index] for index, file_name in enumerate(list_fasta)}

   
    rule fasta:
        input:
            lambda wildcards: dt_extention_fasta_file[wildcards.sample]
        output:
            "results/intermidiate/fasta_files/{sample}.fasta"
        shell:
            "cp {input} {output}"


if len(list_fastq_nanopore) > 0:

    ### sample name for nanopore
    SAMPLES_NAME_global.extend(remove_extensions_file_name(file_name) for file_name in list_fastq_nanopore)

    dt_extention_nano = { remove_extensions_file_name(os.path.basename(file_name)) : list_fastq_nanopore[index] for index, file_name in enumerate( list_fastq_nanopore)}
    
    rule pre_nanopore:
        input:
           lambda wildcards: dt_extention_nano[wildcards.sample]
        output:
            "results/intermidiate/nanofilt_filtred_files/{sample}.fastq.gz"
        shell:
            """(gunzip -c {input} | NanoFilt -q 8 -l 50 --headcrop 30 --tailcrop 30 --maxlength 50000 | gzip > {output}) || (touch {output} && echo Warning: nanofilt failed, were created empty files)"""
 

    rule nanopore:
        input:
            "results/intermidiate/nanofilt_filtred_files/{sample}.fastq.gz"
        output:
            gfa = "results/intermidiate/raven_gfa_files/{sample}.fasta",
            fasta = "results/intermidiate/fasta_files/{sample}.fasta"
        params:
            threads=threads
        shell:
            """raven -t {params.threads} {input} --graphical-fragment-assembly {output.gfa} --disable-checkpoints || (touch {output.gfa} && echo Warning: raven failed, were created empty files) ; awk \'/^S/{{ printf(">%s\\n%s\\n", $2, $3) }}\' {output.gfa} > {output.fasta} || touch {output.fasta} """


if len(list_ab1) > 0:

    SAMPLES_NAME_global.extend(remove_extensions_file_name(file_name,AB1_EXTENSION) for file_name in list_ab1)

    dt_extention_sanger = { remove_extensions_file_name(os.path.basename(file_name),AB1_EXTENSION) : list_ab1[index] for index, file_name in enumerate(list_ab1)}

    rule sanger:
        input:
            lambda wildcards: dt_extention_sanger[wildcards.sample]
        output:
            "results/intermidiate/fasta_files/{sample}.fasta"
        params:
            startbase=startbase,
            endbase=endbase
        shell:
            "abiview {input} -graph none -startbase {params.startbase} -endbase {params.endbase} -osformat2 fasta -outseq {output}  || (touch {output} && echo Warning: abiview failed, were created empty files)"



rule abricate:
    input:
        "results/intermidiate/fasta_files/{sample}.fasta"
    output:
        "results/detailed/{sample}.tab"
    params:
        db = db,
        minid=minid,
        mincov=mincov
    shell:
        "abricate --db {params.db} --nopath {input} --minid {params.minid} --mincov {params.mincov} > {output}"


rule table:
    input:
        expand("results/detailed/{sample}.tab", sample=SAMPLES_NAME_global)
    output:
        expand("results/{output}.csv",output=output)
    script:
        table

