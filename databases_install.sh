#!/bin/bash
echo "Starting ReporType installation!"
output=$(command abricate --help | sed -n 's/.*--datadir \[\(.*\)\].*/\1/p')
path=$(echo "$output" | cut -d '[' -f 2 | cut -d ']' -f 1)
inicial_path=$(pwd)

##Creating measles db
mkdir "$path"/measles
cp "$inicial_path"/databases/measles.fasta "$path"/measles/sequences
cd "$path"/measles && makeblastdb -in sequences -title measles -dbtype nucl -hash_index
abricate --setupdb


##Creating c_trachomatis db
mkdir "$path"/c_trachomatis
cp "$inicial_path"/databases/c_trachomatis.fasta "$path"/c_trachomatis/sequences
cd "$path"/c_trachomatis && makeblastdb -in sequences -title c_trachomatis -dbtype nucl -hash_index
abricate --setupdb

##Creating influenza db
mkdir "$path"/influenza
cp "$inicial_path"/databases/influenza.fasta "$path"/influenza/sequences
cd "$path"/influenza && makeblastdb -in sequences -title influenza -dbtype nucl -hash_index
abricate --setupdb

##Creating newcastle db
mkdir "$path"/newcastle
cp "$inicial_path"/databases/newcastle.fasta "$path"/newcastle/sequences
cd "$path"/newcastle ; makeblastdb -in sequences -title newcastle -dbtype nucl -hash_index
abricate --setupdb    

##Creating dengue db
mkdir "$path"/dengue
cp "$inicial_path"/databases/dengue.fasta "$path"/dengue/sequences
cd "$path"/dengue ; makeblastdb -in sequences -title dengue -dbtype nucl -hash_index
abricate --setupdb

##Creating lp_subspecies_prediction db
mkdir "$path"/lp_subspecies_prediction
cp "$inicial_path"/databases/lp_subspecies_prediction.fasta "$path"/lp_subspecies_prediction/sequences
cd "$path"/lp_subspecies_prediction ; makeblastdb -in sequences -title lp_subspecies_prediction -dbtype nucl -hash_index
abricate --setupdb 


##Creating lp_serogroup_typing db
mkdir "$path"/lp_serogroup_typing
cp "$inicial_path"/databases/lp_serogroup_typing.fasta "$path"/lp_serogroup_typing/sequences
cd "$path"/lp_serogroup_typing ; makeblastdb -in sequences -title lp_serogroup_typing -dbtype nucl -hash_index
abricate --setupdb   

##Creating lp_dot_icm db
mkdir "$path"/lp_dot_icm
cp "$inicial_path"/databases/lp_dot_icm.fasta "$path"/lp_dot_icm/sequences
cd "$path"/lp_dot_icm ; makeblastdb -in sequences -title lp_dot_icm -dbtype nucl -hash_index
abricate --setupdb 


echo""
echo""
echo "Databases available:"
abricate_list=$(abricate --list)
echo "$abricate_list"   