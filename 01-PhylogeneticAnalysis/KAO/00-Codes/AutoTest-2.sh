#!/bin/bash
#Save current directory
codeDir=$PWD

#Read output directory
timestamp="../02-Output/timestamp.txt"

#Loop through timestamp file to find the correct output folder
while read line; do
#Reading each line, separate by comma
entry=(${line//,/ })

if [[ ${entry[0]} == "runID" ]]; then
	echo "run ID found"
	break
fi
done < $timestamp

outDir="../02-Output/${entry[1]}"
LOG_FILE="${outDir}/AutoTest-2.log"
{
echo "Copying scripts.."
cd $outDir
mkdir -p 00-Codes
cp ../../00-Codes/AutoTest-2.sh ./00-Codes/
cp ../../00-Codes/05-ExtractFromSubtree.R ./00-Codes/
cp ../../00-Codes/06-BlastMatrix.R ./00-Codes/

mkdir -p 01-Input
cp ../../01-Input/Outgroup.txt ./01-Input
cp ../../01-Input/Exclude.txt ./01-Input

echo "OutDir located at ${outDir}"

#Retrieve sequences from the manual subtree
cd $codeDir
Rscript 05-ExtractFromSubtree.R


#Run alignment
echo "==============================================="
echo "Aligning with MAFFT!"
cd $outDir
mafft_base="./06-Subtree_LINSI"
mafft_out="${mafft_base}.fasta"
if [[ -f $mafft_out ]]; then
	echo "Previous ${mafft_out} exists! Removed!"
	rm $mafft_out
fi
# mafft --thread 8 --genafpair --maxiterate 1000 05-Subtree.fasta > $mafft_out
mafft --thread 8 --localpair --maxiterate 1000 05-Subtree.fasta > $mafft_out

#Mask alignments again
echo "==============================================="
mask_base="${mafft_base}.98"
if [[ -f "${mask_base}.fasta" ]]; then
	echo "Removing previously masked alignments ..."
	ls |grep "${mask_base}.*"
	rm "${mask_base}.*"
fi

cd $codeDir
Rscript 03-MaskAlignment.R

# Test Evolutionary models & Build tree with iqtree
echo "==============================================="
cd $outDir
tree_base=${mask_base/"06"/"07"}

if [[ -f "${tree_base}.log" ]]; then
	echo "Removing previous iq-tree outputs ..."
	ls |grep "${tree_base}.*"
	rm "${tree_base}.*"
fi

iqtree -nt AUTO -ntmax 8 -s "${mask_base}.phy" -m TEST -mset phyml -B 1000 -alrt 1000 --prefix $tree_base

#Translate the tree
echo "==============================================="
cd $codeDir
Rscript 04-TranslateTree.R

#Run Self blast
echo "==============================================="
echo "Running subtree self blast!"
cd $outDir

makeblastdb -in 05-Subtree.fasta -parse_seqids -dbtype prot -out 05-Subtree_db
	blastp -db 05-Subtree_db -query 05-Subtree.fasta -out 08-Subtree-self_blast.tsv -max_target_seqs 50000 -outfmt "6 std ppos" -num_threads 8
	rm 05-Subtree_db.*

echo "==============================================="
cd $codeDir
#Rscript 06-BlastMatrix.R


say "I finished!"

echo "Test finished!"
}&>$LOG_FILE
