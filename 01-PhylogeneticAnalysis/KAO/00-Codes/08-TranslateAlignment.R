rm(list=ls())
message("Running 08-TranslateAlignment.R !")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",update = T)
if (!exists("read_fasta", mode="function")|!exists("read_hmm", mode="function")) source("00-Functions.R")
# outdir<-file.path(outdir,"..","20230421_011847")

#Load information for species etc.
file.ids<-list.files(outdir,"02-Candidates.*id$",full.names = T)
ids<-fread(file.ids)
ids_long<-ids[,.(Name_sep=unlist(strsplit(Name,"[;\\|]"))),by=.(Num,Name,Clade,Species,Database)]

#read fasta alignment file
OutDir<-outdir
list.phy<-c(file.path(outdir,"05-Subtree.fasta"),
            file.path(outdir,"06-Subtree_LINSI.fasta"),
            file.path(outdir,"06-Subtree_LINSI.98.fasta")
            )

for (f in list.phy){
  if (!file.exists(f)) stop("Target file not found!")
  if (grepl(".phy$",f)){
    phy<-fread(f,header=T)
    names(phy)<-c("Code","Seq_seg")
    phy[,No:=1:.N]
    fasta<-phy[,.(Sequence=paste0(Seq_seg,collapse = "")),by=.(Code)]
    fasta[,newName:=ids_long$Name[match(Code,ids_long$Num)]]
    outname<-file.path(OutDir,sub(".phy",".forhuman.fasta",basename(f)))
    
  }
  if (grepl(".fasta",f)){
    fasta<-read_fasta(f)
    fasta[,newName:=ids_long$Name[match(ID,ids_long$Num)]]
    outname<-file.path(OutDir,sub(".fasta",".forhuman.fasta",basename(f)))
  }
  fasta[,newText:=paste0(">",newName,"\n",Sequence)]
  write(fasta$newText,outname)
  message(paste0(f, " is translated! Outputfile is ",outname))
}
