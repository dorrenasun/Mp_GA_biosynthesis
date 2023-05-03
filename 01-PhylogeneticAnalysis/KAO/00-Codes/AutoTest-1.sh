#!/bin/bash

#Sources
PATH_ONEKP="../../../../00-OneKP/Extracted/OneKP-Nonseed-merged.fasta"
PATH_GENOME="../../../../00-Genomes/Combined/AllGenomes-merged.fasta"
PATH_DONG="../../../../00-Dong42/Combined/Dong42-merged.fasta"
REPEAT_BLAST=True #True or False


#Create folder for output if not exits
mkdir -p ../02-Output

echo "Path_OneKP, $PATH_ONEKP" > ../02-Output/timestamp.txt
echo "Path_Genome, $PATH_GENOME" >> ../02-Output/timestamp.txt
echo "Path_Dong, $PATH_DONG" >> ../02-Output/timestamp.txt


echo "==============================================="
echo "REPEAT_BLAST: $REPEAT_BLAST"
BLAST_OUT_ONEKP="01-OneKP-Nonseed-blast.tsv"
BLAST_OUT_GENOME="01-AnnotatedGenomes-blast.tsv"
BLAST_OUT_DONG="01-Dong42-blast.tsv"

if [[ $REPEAT_BLAST == True ]] || [[ ! -f "../02-Output/$BLAST_OUT_GENOME" ]] || [[ ! -f "../02-Output/$BLAST_OUT_DONG" ]] || [[ ! -f "../02-Output/$BLAST_OUT_ONEKP" ]]
then
	# #Blast against known genomes
	echo "==============================================="
	echo "Running blast against annotated genome!"
	cd ../02-Output/
	makeblastdb -in $PATH_GENOME -parse_seqids -dbtype prot -out 00-Genomes_db
	blastp -db 00-Genomes_db -query ../01-Input/Blast_Input.fasta -out 01-AnnotatedGenomes-blast.tsv -evalue '0.001' -max_target_seqs 50000 -outfmt "6 std ppos" -num_threads 8
	rm 00-Genomes_db.*

	#Blast against OneKP
	echo "Running blast against OneKP!"
	makeblastdb -in $PATH_ONEKP -parse_seqids -dbtype prot -out 00-OneKP_db
	blastp -db 00-OneKP_db -query ../01-Input/Blast_Input.fasta -out 01-OneKP-Nonseed-blast.tsv -evalue '0.001' -max_target_seqs 50000 -outfmt "6 std ppos" -num_threads 8
	rm 00-OneKP_db.*

	#Blast against Dong42
	echo "Running blast against Dong42!"
	makeblastdb -in $PATH_DONG -parse_seqids -dbtype prot -out 00-Dong_db
	blastp -db 00-Dong_db -query ../01-Input/Blast_Input.fasta -out 01-Dong-blast.tsv -evalue '0.001' -max_target_seqs 50000 -outfmt "6 std ppos" -num_threads 8
	rm 00-Dong_db.*
	cd ../00-Codes


else
	echo "REPEAT_BLAST is False; Skipping blast."
fi

#Create a timestamp for a particular run
runID=$(date +"%Y%m%d_%H%M%S")
echo "runID, $runID" >> ../02-Output/timestamp.txt

#Create folder for output
outDir="../02-Output/$runID"
mkdir -p $outDir

#Create log file for each run
LOG_FILE="${outDir}/AutoTest-1.log"
{
echo "Copying scripts.."
cd $outDir
mkdir -p 00-Codes
cp ../../00-Codes/AutoTest-1.sh ./00-Codes/
cp ../../00-Codes/00-Functions.R ./00-Codes/
cp ../../00-Codes/01-RetrieveFasta.R ./00-Codes/
cp ../../00-Codes/02-MergeFasta.R ./00-Codes/
cp ../../00-Codes/03-MaskAlignment.R ./00-Codes/
cp ../../00-Codes/04-TranslateTree.R ./00-Codes/

mkdir -p 01-Input
cp ../../01-Input/Blast_Input.fasta ./01-Input
cp ../../01-Input/p450.hmm ./01-Input



echo "==============================================="
echo "Retrieving seqs!"
cd ../../00-Codes

Rscript 01-RetrieveFasta.R
echo "==============================================="
Rscript 02-MergeFasta.R



# # echo "Checking with self-blast!"
# cd $outDir
# makeblastdb -in 02-Candidates-all.fasta -parse_seqids -dbtype prot -out Candidates-all_db
# blastp -db Candidates-all_db -query 02-Candidates-all.fasta -out 03-Candidates-all-self-blast.tsv -evalue 0.001 -outfmt "6 std ppos" -num_threads 8
# rm Candidates-all_db.*

echo "==============================================="
echo "Checking with hmmsearch!"
cd $outDir
hmmsearch --tblout 03-Candidates-all_P450.tsv ../../01-Input/p450.hmm 02-Candidates-all.fasta > 03-Candidates-all_P450.out


# align with MAFFT
echo "==============================================="
echo "Making alignment with mafft!"
mafft --retree 2 02-Candidates-all-fl.fasta > 03-Candidates-all_FFTNS2.fasta

cd ../../00-Codes
Rscript 03-MaskAlignment.R

# Test Evolutionary models & Build tree with iqtree
echo "==============================================="
echo "Initial tree with iqtree2!"
cd $outDir
iqtree -nt AUTO -ntmax 8 -s ./03-Candidates-all_FFTNS2.80.phy -m TEST -mset phyml -fast -alrt 1000 --prefix ./04-Candidates-all_FFTNS2_80
# iqtree -nt AUTO -ntmax 8 -s ./03-Candidates-all_EINSI.80.phy -m TEST -mset phyml -B 1000 -alrt 1000 --prefix ./04-Candidates-all_EINSI_80

cd ../../00-Codes
Rscript 04-TranslateTree.R
say "I finished!"
# echo "Mafft alignment!"
# cd 03-Retrieved



echo "Test finished!"

}&>$LOG_FILE

