rm(list = ls())
if (!require(data.table)) install.packages("data.table")
if (!require(R.utils)) install.packages("R.utils")
library(data.table)
library(R.utils)
library(stringr)

#Function to read all fasta files
read_fasta<-function(file){
  if (!file.exists(file)) stop("Fasta file not found!")
  Fasta<-fread(file,header = F,sep="\t",col.names = "Full")
  ID.positions<-grep(">",Fasta$Full,fixed = T)
  ID.lengths<-c(ID.positions[-1],nrow(Fasta)+1)-ID.positions
  
  Fasta[,`:=`(Num=rep(1:length(ID.positions),ID.lengths),
              ID=rep(sub(" .*","",Full[ID.positions]),ID.lengths))]
  Fasta[,`:=`(ID=sub(">","",ID,fixed = T),
              Sequence=Full)]
  Fasta[ID.positions,Sequence:=""]
  Fasta.trans<-Fasta[,.(Sequence=paste(Sequence,collapse = "")),by=.(Num,ID)]
  # Fasta.trans[,`:=`(ID=gsub("|","+",ID,fixed = T),
  #                   Sequence=gsub("*","",Sequence,fixed = T))]
  Fasta.trans[,`:=`(ID=gsub("|","+",ID,fixed = T))]
  Fasta.trans[,ID:=sub("[ |\t].*","",ID)]
  Fasta.trans[,ID:=sub("scaffold-","",ID,fixed = T)]
  Fasta.trans[,ID:=gsub("\\(","{",ID)]
  Fasta.trans[,ID:=gsub("\\)","}",ID)]
  Fasta.trans[grepl("PIUF",ID,fixed=T),ID:=sub("Pellia_sp._{cf_epiphylla_{L.}_Corda}","Pellia_cf._epiphylla",ID,fixed = T)]
  
  #Only keep the first two or three words in taxon name
  Fasta.trans[grepl("cf",ID,fixed = T),newID:=sub("^([A-Z]+\\-[0-9]+\\-[[:alpha:]]+_cf\\._[[:alpha:]]+).*","\\1",ID)]
  Fasta.trans[!grepl("cf",ID,fixed = T),newID:=sub("^([A-Z]+\\-[0-9]+\\-[[:alpha:]]+_[[:alpha:]]+).*","\\1",ID)]
  
  #Fasta.trans[,Text:=paste(paste0(">",newID),Sequence,sep = "\n")]
  
  return(Fasta.trans)
}

# inDir<-c("Download/Mosses","Download/Liverworts","Download/Hornworts")
inDir<-"./Download"
list.file<-list.files(inDir,pattern = "translated-protein",full.names = T,recursive = T)

fasta.all<-data.table()
for (f in list.file){
  print(paste0("Processing: ",f))
  temp<-read_fasta(f)
  temp$Source<-f
  #Remove "scaffold" from gene names
  fasta.all<-rbind(fasta.all,temp)
}
# stop()
#Check if there are redundant IDs
message("Checking redundancy in IDs..")
if (sum(duplicated(fasta.all$ID))>0) warning("Duplicated IDs! Wrong parsing!")

#Add marks to IDs, and truncate IDs longer than 50
prefixes<-fread("Prefixes.csv")
fasta.all[,Clade:=sub(".*\\/","",dirname(Source))]
fasta.all<-merge(fasta.all,prefixes,by="Clade",all.x = T)
if (sum(is.na(fasta.all$Prefix)>0)) message("Sequences of unknown source!")
fasta.all[,newID1:=paste0(Prefix,"_",newID)]

#Confirm ID format with sampling
fasta.sample<-fasta.all[!duplicated(Source),]

strlim<-50
fasta.all[nchar(newID1)>strlim,newID1:=str_trunc(newID1,strlim,ellipsis = "..")]
#Check if there are still IDs longer than 50
if (max(nchar(fasta.all$newID1))>strlim) warning("IDs longer than 50!")

fasta.all[,Text:=paste(paste0(">",newID1),Sequence,sep = "\n")]
sources<-unique(fasta.all[,.(Clade,Source)])
sources[,Code:=sub(".*\\/([A-Z]+)-.*","\\1",Source)]

#Retrieve Species Names
Info<-fread(file.path(inDir,"Sample-List-with-Taxonomy.tsv.csv"))
fasta.all[,Code:=sub("([A-Z]+)-.*","\\1",ID)]
fasta.all[,Species:=Info$Species[match(Code,Info$`1kP_Sample`)]]
fasta.all[,Clade:=Info$Clade[match(Code,Info$`1kP_Sample`)]]
#Substitute some Clade names
# fasta.all[grepl("Monilophytes", Clade), Clade:="Fern"]

{
  outDir<-"Extracted/"
  outbase<-unlist(strsplit(inDir,"[\\.\\/]+"))
  outbase<-"OneKP-Nonseed"
  
  list.old<-list.files(outDir,full.names = T,recursive = T)
  if (sum(grep(outbase,list.old))>0) file.remove(list.old[grep(outbase,list.old)])
  
  write(fasta.all$Text,paste0(outDir,outbase,"-merged.fasta"))
  write.table(fasta.all[,.(newID1,Source,Code,Clade,Species)],paste0(outDir,outbase,"-ID2Source.txt"),sep = "\t",quote = F,row.names = F,col.names = T)

}
message("OneKP Genomes merged!")