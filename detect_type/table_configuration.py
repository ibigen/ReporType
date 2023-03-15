import pandas as pd
import os
import csv
import numpy as np

pd.options.mode.chained_assignment = None  

def results_path(path):
    FILES=[]
    for p,_,files in os.walk(os.path.abspath(path)):
        for file in files:
            file=os.path.join(p, file)
            FILES.append (file)
    return (FILES)


def types (list,coverage,id,gene,acce,seq):
    SUB_I=[]
    id_type,cov_type,gene_type,acce_type,seq_type,id_sub,cov_sub,gene_sub,acce_sub,seq_sub=(None,None,None,None,None,None,None,None,None,None)
    for i in range(len(list)):
        if list[i] == 'species' or list[i] == 'lineage' or list[i] == 'type' :#######acrescentar outros
            id_type=gene[i]+'-'+str(id[i])
            cov_type=gene[i]+'-'+str(coverage[i])
            acce_type=gene[i]+'-'+str(acce[i])
            gene_type=gene[i]
            seq_type=gene[i]+'-'+seq[i]
        if list[i] == 'subtype' or list[i] == 'genus':
            SUB_I.append(int(i))
            ID=[]
            COV=[]
            GENE=[]
            ACCE=[]
            SEQ=[]
            for x in SUB_I:
                id_sueb=str(gene[x])+'-'+str(id[x])
                cov_sueb=str(gene[x])+'-'+str(coverage[x])
                acce_sueb=str(gene[x])+'-'+str(acce[x])
                gene_sueb=str(gene[x])
                seq_sueb=str(gene[x])+'-'+str(seq[x])
                ID.append(id_sueb)
                COV.append(cov_sueb)
                GENE.append(gene_sueb)
                ACCE.append(acce_sueb)
                SEQ.append(seq_sueb)
            GENE=sorted(GENE)
            ID=sorted(ID)
            COV=sorted(COV)
            ACCE=sorted(ACCE)
            SEQ=sorted(SEQ)
            join=' : '
            join1=''
            id_sub=join.join(ID)
            cov_sub=join.join(COV)
            gene_sub=join1.join(GENE)
            acce_sub=join.join(ACCE)
            seq_sub=join.join(SEQ)       
    return id_type,cov_type,gene_type,acce_type,seq_type,id_sub,cov_sub,gene_sub,acce_sub,seq_sub

def check_only_type(database,TYPES):
    for i in range(len(TYPES)):
        if database == TYPES[i]:
            unique=True
        else:
            unique=False
    if database=="ompA":
        unique=True
    if database=="influenza" or database=="coronavirus":
        unique=False
    return(unique)

FILES = snakemake.input
if len(FILES)==0:
    print("  ERROR: Some error occured during the analysis, verify your samples  ")
    
    
file_name=snakemake.output[0]


def transform_in_list(obj):
   LIST=[]
   for line in obj:
       line=str(line)
       LIST.append(line)
   return(LIST)

for file in FILES:
    file_in_name=str(file)
    file_in_name=file_in_name.strip()
    file_in_name=file_in_name.split("/")[-1]
    file_in_name=file_in_name.split(".")[0]
    file = pd.read_table(file,header=0,usecols=['#FILE','SEQUENCE','GENE','%COVERAGE','%IDENTITY','DATABASE','ACCESSION'])
    file.sort_values(by=['GENE'])
    n=len(file)
    if n == 0:
        file= pd.DataFrame({'#FILE':[file_in_name],'SEQUENCE':["---"],'GENE':["0 genes found (check if there is any warning during the analysis)"],'%COVERAGE':["---"],'%IDENTITY':["---"],'DATABASE':["---"],'ACCESSION':["---"],})
        file.to_csv(file_name,sep="\t",header=False, index=False,mode='a')
    else:
        file=file.sort_values(by=['%COVERAGE'], ascending=False)
        file=file.drop_duplicates(subset="GENE")
        n=len(file)
        type=file["DATABASE"]
        TYPES=[]
        DBS=[]
        for line in type:
            line=line.strip()
            database=line.split('_')[0]
            type=line.split('_')[-1]
            if type !='DATABASE':
                TYPES.append(type)
            if database!='DATABASE':
                DBS.append(database)
        database=DBS[0]
        unique= check_only_type(database,TYPES)
        if unique == True:
            file.to_csv(file_name,sep="\t",header=False, index=False,mode='a')
        else:
            file["DATABASE"]=[database]*n
            coverage=file["%COVERAGE"]
            coverage=transform_in_list(coverage)
            id=file['%IDENTITY']
            id=transform_in_list(id)
            gene=file['GENE']
            gene=transform_in_list(gene)
            acce=file['ACCESSION']
            acce=transform_in_list(acce)
            seq=file['SEQUENCE']
            seq=transform_in_list(seq)
            id_type,cov_type,gene_type,acce_type,seq_type,id_sub,cov_sub,gene_sub,acce_sub,seq_sub=types(TYPES,coverage,id,gene,acce,seq)
            if gene_type is None:
                file['%IDENTITY']=[str(id_sub)]*n
                file['GENE']=[str(gene_sub)]*n
                file['%COVERAGE']=[str(cov_sub)]*n
                file['ACCESSION']=[str(acce_sub)]*n
                file['SEQUENCE']=[str(seq_sub)]*n
            elif gene_sub is None:
                file['%IDENTITY']=[str(id_type)]*n
                file['GENE']=[str(gene_type)]*n
                file['%COVERAGE']=[str(cov_type)]*n
                file['ACCESSION']=[str(acce_type)]*n
                file['SEQUENCE']=[str(seq_type)]*n        
            else:
                join=" : "
                file['%IDENTITY']=[str(id_type)+str(join)+str(id_sub)]*n
                file['GENE']=[str(gene_type)+'-'+str(gene_sub)]*n
                file['%COVERAGE']=[str(cov_type)+str(join)+str(cov_sub)]*n
                file['ACCESSION']=[str(acce_type)+str(join)+str(acce_sub)]*n
                file['SEQUENCE']=[str(seq_type)+str(join)+str(seq_sub)]*n
            file=file[:1]
            file.to_csv(file_name,sep="\t",header=False, index=False,mode='a')
        


final_file=pd.read_csv(file_name, sep='\t',names=['SAMPLE','SEQUENCE','GENE','COVERAGE(%)','IDENTITY(%)','DATABASE','ACCESSION'])

final_file.to_csv(file_name,header=True, index=False,mode='w')


f_file=pd.read_csv(snakemake.output[0])

duplicados=f_file.duplicated(subset=["SAMPLE"],keep=False)



for x in range(len(duplicados)):
    if duplicados[x]==True:
        f_file.loc[[x],["SAMPLE"]]=np.nan
        
        
f_file["SAMPLE"]=f_file["SAMPLE"].fillna(f_file["SEQUENCE"])

for x in f_file.index:
    f_sample=str(f_file.loc[[x],["SAMPLE"]])
    f_sample=f_sample.strip()
    f_sample=f_sample.split(".")[0]
    f_sample=f_sample.split()[2]
    f_file["SAMPLE"][x]=f_sample


f_file=f_file.drop(columns=["SEQUENCE"])


f_file.to_csv(file_name,sep='\t',header=True, index=False,mode='w')





