rm(list = ls())

if (!exists("read_fasta", mode="function")) source("00-Functions.R")
inDir<-"./Resources/"


{ #Subsetting
  if(!require(readxl)) install.packages("readxl")
  library(readxl)

  fSource<-"./List_of_sources.xlsx"
  list_src<-read_excel(fSource,sheet = 1,range = cell_cols("A:J"),col_names = T)
  subset<-list_src$`File name`[!is.na(list_src$Include) & list_src$Include=="T"]
  list.file<-paste0(inDir,unique(subset))
}

# list.file<-list.files(inDir,full.names = T)
fasta.all<-data.table()
nonexist<-c()
for (f in list.file){
  if (!file.exists(f)){
    warning(paste0(f,"does not exist!"))
    nonexist<-c(nonexist,f)
    next
  } 
  temp<-read_fasta(f)
  temp$Source<-basename(f)
  fasta.all<-rbind(fasta.all,temp)
}

#Write non-exist files to a file
if (length(nonexist)>0) write(nonexist,file="01-Genomes/Non_exist.txt")

#Check if there are redundant IDs
if (sum(duplicated(fasta.all$ID))>0) {
  warning("Duplicated IDs! Wrong parsing!")
  list.dup<-fasta.all[duplicated(ID),ID]
  fasta.dup<-fasta.all[ID %in% list.dup,]
}
##Parsing IDs to genes and variants
fasta.all[,pID:=ID]

#data.table for checking file names
fasta.all.test<-fasta.all[!duplicated(Source),]

#Remove the .p suffix from Phytozome datasets
list.endp<-unique(fasta.all[grep(".p$",ID),Source])
if (length(list.endp)>0) {
  fasta.all[Source %in% list.endp,pID:=sub("(.*).p$","\\1",ID)]
  message(paste0("Removed .p for sequences from:\n",paste(list.endp,collapse = "\n")))
  message("\n")
  }


#Change -XX style of transcript variants to .XX
list.dash<-unique(fasta.all[grep("-[0-9]+$",pID),Source])
if (length(list.dash)>0) {
  fasta.all[Source %in% list.dash,pID:=sub("-([0-9]+)$",".\\1",pID)]
  message(paste0("Transformed from:\n",paste(list.dash,collapse = "\n")))
  message("\n")
}

#Remove the "evm.model." prefix if any
list.evm<-unique(fasta.all[grep("evm.*model.",ID),Source])
if (length(list.evm)>0) {
  fasta.all[Source %in% list.evm,pID:=sub("evm.*model.","",ID)]
  message(paste0("Removed evm.model. for sequences from:\n",paste(list.evm,collapse = "\n"))) 
  message("\n")
}

#Remove the "jgi+" prefix if any
list.jgi<-unique(fasta.all[grep("jgi",ID),Source])
if (length(list.evm)>0) {
  fasta.all[Source %in% list.jgi,pID:=sub("jgi+","",ID,fixed = T)]
  message(paste0("Removed jgi+ for sequences from:\n",paste(list.jgi,collapse = "\n"))) 
  message("\n")
}


#Add labels for certain genomes
#Ginkgo
fasta.all[Source=="Gbi_new_chr.pep.fa",pID:=paste0("Gbi_",pID)]
#C.papaya
fasta.all[Source=="Cpapaya_113_ASGPBv0.4.protein.fa.gz",pID:=paste0("Cpa_",pID)]
#Selaginella
fasta.all[Source=="Smoellendorffii_91_v1.0.protein.fa.gz",pID:=paste0("Sm_",pID)]



#Sort sequence by transcript variant
fasta.all[,`:=`(Gene=sub("(.*)\\.[0-9]+$","\\1",pID),
             var=sub(".*\\.([0-9]+)$","\\1",pID))]
fasta.all[grep("[A-z]",var),var:=0]
fasta.all[,var:=as.integer(var)]

#Check if there are non-variant IDs
fasta.vars<-fasta.all[,.(nvar=length(unique(var))),by=Source]
list.nonvar<-unique(c(fasta.all[var>100,Source],fasta.vars[nvar==1,Source]))

if (length(list.nonvar)>0) {
  fasta.all[Source %in% list.nonvar,`:=`(Gene=pID,
                                         var=0)]
  message(paste0("Genomes with no transcript variant:\n",paste(list.nonvar,collapse = "\n"))) 
  message("\n")
}

#Retrieve Species and taxon names
fasta.all[,Species:=list_src$Species[match(Source,list_src$`File name`)]]
fasta.all[,Clade:=list_src$Clade[match(Source,list_src$`File name`)]]

#Check Sequence Integrity
fasta.all[,`:=`(Start=substr(Sequence,1,1),
                End=substr(Sequence,nchar(Sequence),nchar(Sequence))
)]

fasta.all[,`:=`(N_comp=(Start=="M"),
                C_comp=(End=="*"))]
#Check if the integrity info was provided in the original file
Integrity<-fasta.all[,.(Species=unique(Species),Total=.N,Sum_N=sum(N_comp),Sum_C=sum(C_comp)),by=.(Source)]
Integrity[,NoStop:=(Sum_C/Total)<0.1]
list.NoStop<-unlist(Integrity[!(is.na(NoStop)) & (NoStop==T),.(Source)])
fasta.all[Source %in% list.NoStop,C_comp:=NA]

fasta.all[!is.na(C_comp) & C_comp==T, Sequence:=substr(Sequence,1,nchar(Sequence)-1)]

#Combine final IDs and sequences
fasta.all[,Text:=paste0(">",pID,"\n",Sequence)]

fasta.example<-fasta.all[!duplicated(Source),.(ID,pID,Gene,var,Species),by=Source]

#Write files
outdir<-"./Combined/"
if (!dir.exists(outdir)) dir.create(outdir)
outbase<-"AllGenomes"
list.old<-list.files(outdir,full.names = T,recursive = T)
if (sum(grep(outbase,list.old))>0) file.remove(list.old[grep(outbase,list.old)])

write(fasta.all$Text,paste0(outdir,outbase,"-merged.fasta"))
write.table(fasta.all[,.(pID,Gene,var,Source,Clade,Species)],paste0(outdir,outbase,"-ID2Source.txt"),col.names = T,row.names = F,quote = F,sep = "\t")
write.table(fasta.example,paste0(outdir,outbase,"-example.txt"),col.names = T,row.names = F,quote = F,sep = "\t")
write.table(Integrity[!(is.na(NoStop)) & (NoStop==T),.(Species,File=Source)],paste0(outdir,"Genome_of_No_Integrity.tsv"),col.names = T,row.names = F,quote = F,sep = "\t")

message("Genomes merged!")
