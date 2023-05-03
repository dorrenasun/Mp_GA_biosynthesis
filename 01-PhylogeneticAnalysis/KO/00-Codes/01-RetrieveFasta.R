rm(list = ls())
message("Running 01-RetrieveFasta.R !")
if (!require(readxl)) install.packages("readxl")
library(readxl)

if (!exists("read_fasta", mode="function")) source("00-Functions.R")
cutoff=50

#Retrieve blast sequences from annotated genomes
{
  message("\nRetrieving fasta for genome hits!")
  inFile.fasta.Genome<-list.files(path = Genomedir,pattern = "merged.fasta",full.names = T)
  inFile.ID2Src.Genome<-list.files(path = Genomedir,pattern = "ID2Source",full.names = T)
  inFile.tsv.Genome<-"../02-Output/01-AnnotatedGenomes-blast.tsv"
  inFile.list.Genome<-file.path(Genomedir,"..","List_of_sources.xlsx")
  
  #Copy the input files to output directory
  InputCopydir<-file.path(outdir,"01-Input")
  if (!dir.exists(InputCopydir)) dir.create(InputCopydir,recursive = T)
  file.copy(inFile.ID2Src.Genome,InputCopydir)
  file.copy(inFile.list.Genome,InputCopydir)
  file.copy(inFile.tsv.Genome,outdir)
  
  
  #Load blast results
  tsv.Genome<-fread(inFile.tsv.Genome,header = F)
  names(tsv.Genome)<-unlist(strsplit("qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos",split=" "))

  #Correct "SM_"issues caused by blast
  list.SM<-unique(tsv.Genome[grep(saccver,pattern = "SM_",fixed = T),saccver])
  if (length(list.SM)>0){
    message(paste0("Correct names for:\n",paste0(list.SM,collapse = "\n" )))
    tsv.Genome[saccver %in% list.SM,saccver:=sub("SM_","Sm_",saccver)]
  }
  #Correct "ADC"issues caused by blast
  list.ADC<-unique(tsv.Genome[grep(saccver,pattern = "ADC",fixed = T),saccver])
  if (length(list.ADC)>0){
    message(paste0("Correct names for:\n",paste0(list.ADC,collapse = "\n" )))
    tsv.Genome[saccver %in% list.ADC,saccver:=sub("ADC","Adc",saccver)]
  }
  
  ID2Src.Genome<-fread(inFile.ID2Src.Genome)
  #merge with ID2Src
  tsv.Genome<-merge(tsv.Genome,ID2Src.Genome,by.x = "saccver",by.y = "pID",all.x = T)
  #Check any IDs not found in ID2Src
  if (sum(is.na(tsv.Genome$Source))>0){
    message("Some blast hits not found in Genome ID2Src! Problem with parsing")
  }
   
  #Load list of Genomes to include  
  list.Genome<-as.data.table(read_excel(inFile.list.Genome,sheet = 1,range = cell_cols("A:J") ))
  #Check if there are any unmatches in Source names
  list.unmatch<-ID2Src.Genome[!Source %in% list.Genome$`File name`,Source]
  if (length(list.unmatch)>0) {
    message("Unmatching sources! Check your genome parsing!")
    message(paste(list.unmatch,sep = "\n"))
  }  
  list.Src<-list.Genome[(!is.na(Include))&(Include=="T"),`File name`]

  #Filter ID by selected genome
  tsv.Genome.sele<-tsv.Genome[Source %in% list.Src]

  #Select the top N genes from each species
  setkeyv(tsv.Genome.sele,c("qaccver","Source","evalue"))
  tsv.Genome.sele[,No:=match(Gene,unique(Gene)),by=c("qaccver","Source")]
  tsv.Genome.cut<-tsv.Genome.sele[No<cutoff,]
  
  #Keep only the best hit from the primary transcript for each gene
  setkeyv(tsv.Genome.cut,c("Gene","var","evalue"))
  tsv.Genome.cut[,primary:=(var==min(var) & !duplicated(var)),by=.(Gene)]
  tsv.Genome.final<-tsv.Genome.cut[primary==T,]
  
  

  #Load fasta sequences
  fasta.Genome<-read_fasta(inFile.fasta.Genome)
  #Check if there are any unmatches in blast hits
  list.missing<-setdiff(tsv.Genome.final$saccver, fasta.Genome$ID)
  if (length(list.missing)>0) {
    message("Unmatching IDs! Some Genome blast hits not found!") 
    message(paste(list.missing,sep = "\n"))
  } 
  
  fasta.Genome.sele<-fasta.Genome[ID %in% tsv.Genome.final$saccver,]
 
  write(fasta.Genome.sele$Text,paste0(outdir,"01-AnnotatedGenomes-blast_hits.fasta"))
  message("Fasta retrieved for genome hits!")
  
}

# Retrieve OneKP results
{
  message("\nRetrieving fasta for onekp hits!")
  inFile.fasta.OneKP<-list.files(path = OneKPdir,pattern = "merged.fasta",full.names = T)
  inFile.tsv.OneKP<-"../02-Output/01-OneKP-Nonseed-blast.tsv"
  inFile.list.OneKP<-"../01-Input/1kP-Sample-List-2022.tsv"
  
  #Copy the input files to output directory
  file.copy(inFile.tsv.OneKP,outdir)
  file.copy(inFile.list.OneKP,InputCopydir)
  
  tsv.OneKP<-fread(inFile.tsv.OneKP,header = F)
  names(tsv.OneKP)<-unlist(strsplit("qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos",split=" "))
 
  
  #Load list of OneKP items to include 
  list.OneKP<-fread(inFile.list.OneKP,header = T)
  #Check if there are any unmatches in OneKP codes
  tsv.OneKP[,Code:=sub(".*_([A-Z]+)-[0-9].*","\\1",saccver)]
  list.unmatch<-tsv.OneKP[!Code %in% list.OneKP$Code,saccver]
  if (length(list.unmatch)>0) warning("Unmatching Codes! Check your OneKP parsing!")  
  
  
  #Filter by OneKP sources and ranking in the same species
  list.OneKP.sele<-list.OneKP[(!is.na(Include))&(Include=="T"),Code]
  tsv.OneKP.sele<-tsv.OneKP[Code %in% list.OneKP.sele,]
  setkeyv(tsv.OneKP.sele,c("qaccver","Code","evalue"))
  tsv.OneKP.sele[,No:=seq(1:.N),by=c("qaccver","Code")]
  tsv.OneKP.cut<-tsv.OneKP.sele[No<cutoff,]
  
  #Keep the best hit for each gene
  setkeyv(tsv.OneKP.cut,c("saccver","evalue"))
  tsv.OneKP.final<-tsv.OneKP.cut[!duplicated(saccver)]
  
  #Load OneKP fasta
  fasta.OneKP<-read_fasta(inFile.fasta.OneKP)
  #Check if there are any missing IDs
  list.missing<-setdiff(tsv.OneKP.final$saccver, fasta.OneKP$ID)
  if (length(list.missing)>0) {
    message("Unmatching IDs! Some OneKP blast hits not found!")  
    message(paste0(list.missing,sep = "\n"))
    }
  
  fasta.OneKP.sele<-fasta.OneKP[ID %in% tsv.OneKP.final$saccver,]
  write(fasta.OneKP.sele$Text,paste0(outdir,"01-OneKP-blast_hits.fasta"))
  message("Fasta retrieved for onekp hits!")
  
    
}

# Retrieve Dong42 results
{
  message("\nRetrieving fasta for Dong42 hits!")
  inFile.fasta.Dong<-list.files(path = Dongdir,pattern = "merged.fasta",full.names = T)
  inFile.tsv.Dong<-"../02-Output/01-Dong-blast.tsv"
  
  file.copy(inFile.tsv.Dong,outdir)
  # inFile.list.Dong<-"../01-Input/1kP-Sample-List-2022.tsv"
  
  tsv.Dong<-fread(inFile.tsv.Dong,header = F)
  names(tsv.Dong)<-unlist(strsplit("qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos",split=" "))
  
  #Exclude seqs from Marchantia, P.patens or Smoellendorffi
  tsv.Dong.sele<-tsv.Dong[!grepl("Mpolymorpha",saccver) & !grepl("Ppatens",saccver) & !grepl("Smoellendorffii",saccver),]
  
  #Retrieve species from each seq
  tsv.Dong.sele[,Source:=sub("_.*","",saccver)]
  if (nrow(tsv.Dong.sele[grepl("[0-9]+",Source),])>0) message("Problem parsing saccver sources!")
  
  #Recongize gene, isoforms and different proteins from TRINITY output
  tsv.Dong.sele[,`:=`(Gene=sub("_i.*\\.p.*","",saccver),
                 var=sub(".*_(i.*\\.p.*)","\\1",saccver)
  )]
  
  
  #Select the top N genes from each species
  setkeyv(tsv.Dong.sele,c("qaccver","Source","evalue"))
  tsv.Dong.sele[,No:=match(Gene,unique(Gene)),by=c("qaccver","Source")]
  tsv.Dong.cut<-tsv.Dong.sele[No<cutoff,]
  
  
  #Keep only the best hit from each gene
  setkeyv(tsv.Dong.cut,c("Gene","evalue"))
  tsv.Dong.cut[,primary:=(var==min(var) & !duplicated(var)),by=.(Gene)]
  tsv.Dong.final<-tsv.Dong.cut[primary==T,]
  
  
  #Load Dong fasta
  fasta.Dong<-read_fasta(inFile.fasta.Dong)
  #Check if there are any missing IDs
  list.missing<-setdiff(tsv.Dong.final$saccver, fasta.Dong$ID)
  if (length(list.missing)>0) {
    message("Unmatching IDs! Some Dong blast hits not found!")  
    message(paste(list.missing,sep = "\n"))
  }
  
  
  fasta.Dong.sele<-fasta.Dong[ID %in% tsv.Dong.final$saccver,]
  write(fasta.Dong.sele$Text,paste0(outdir,"01-Dong-blast_hits.fasta"))
  message("Fasta retrieved for Dong hits!")
  
  
}
