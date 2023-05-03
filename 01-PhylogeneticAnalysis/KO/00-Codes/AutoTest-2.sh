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

runID=$(date +"%Y%m%d_%H%M%S")


#Create folder for output
outDir="${codeDir}/../02-Output/${entry[1]}"
echo "runID, $runID" > ${outDir}/timestamp.txt

#mkdir -p $outDir


LOG_FILE="${outDir}/AutoTest-2.log"
{

echo "OutDir located at ${outDir}"
echo "Copying scripts while running.."
cd $outDir
mkdir -p $runID
cd $runID

mkdir -p 00-Codes
mkdir -p 01-Input

cp ${codeDir}/AutoTest-2.sh ./00-Codes/
cp ${codeDir}/00-Functions.R ./00-Codes/
cp ${codeDir}/05-ExtractFromSubtree.R ./00-Codes/

cp ${codeDir}/../01-Input/Outgroup.txt ./01-Input
cp ${codeDir}/../01-Input/Exclude.txt ./01-Input


#Retrieve sequences from the manual subtree
cd $codeDir
Rscript 05-ExtractFromSubtree.R
cp ${outDir}/05-* ${outDir}/${runID}/


#Run alignment
echo "==============================================="
echo "Aligning with MAFFT!"
cd $outDir

algo="LINSI"
mafft_base="06-Subtree_${algo}"
mafft_out="${mafft_base}.fasta"
if [[ -f $mafft_out ]]; then
	echo "Previous ${mafft_out} exists! Removed!"
	rm $mafft_out
fi
if [ $algo = "EINSI" ]; then
	echo "Running EINSI algorithm"
	mafft --thread 8 --genafpair --maxiterate 1000 05-Subtree.fasta > $mafft_out
fi

if [ $algo = "LINSI" ]; then
	echo "Running LINSI algorithm"
	mafft --thread 8 --localpair --maxiterate 1000 05-Subtree.fasta > $mafft_out
fi

cp $mafft_out ./${runID}/

#Mask alignments again
echo "==============================================="
ratio=98


mask_base="${mafft_base}.$ratio"
if [[ -f "${mask_base}.fasta" ]]; then
	echo "Removing previously masked alignments ..."
	files_oldmask=`(ls | grep "${mask_base}*")`
	for file in $files_oldmask; do
	  # delete each file found
	  rm "$file"
	  echo "Deleted $file"
	done
fi

cd $codeDir
cp 03-MaskAlignment.R ${outDir}/${runID}/00-Codes/
Rscript 03-MaskAlignment.R $ratio

cd $outDir
cp ${mask_base}* ./${runID}/


# Test Evolutionary models & Build tree with iqtree
echo "==============================================="

tree_base=${mask_base/"06"/"07"}

if [[ -f "${tree_base}.log" ]]; then
	echo "Removing previous iq-tree outputs ..."
	files_oldtree=`(ls | grep "${tree_base}*")`
	for file in $files_oldtree; do
	  # delete each file found
	  rm "$file"
	  echo "Deleted $file"
	done
fi
iqtree -nt AUTO -ntmax 8 -s "${mask_base}.phy" -m TEST -mset phyml -B 1000 -alrt 1000 --prefix $tree_base -redo

#Translate the tree
echo "==============================================="
cd $codeDir
cp 04-TranslateTree.R ${outDir}/${runID}/00-Codes/
Rscript 04-TranslateTree.R
cd $outDir
cp ${tree_base}* ./${runID}/


#Run Self blast
echo "==============================================="
echo "Running subtree self blast!"

makeblastdb -in 05-Subtree.fasta -parse_seqids -dbtype prot -out 05-Subtree_db
	blastp -db 05-Subtree_db -query 05-Subtree.fasta -out 08-Subtree-self_blast.tsv -max_target_seqs 50000 -outfmt "6 std ppos" -num_threads 8
	rm 05-Subtree_db.*

cp 08-Subtree-self_blast.tsv ./${runID}/

echo "==============================================="
#Rscript 06-BlastMatrix.R

echo "Copying log.."
cd ${outDir}

cp AutoTest-2.log ./${runID}/


}&>$LOG_FILE


say "I finished!"
echo "Test finished!"
