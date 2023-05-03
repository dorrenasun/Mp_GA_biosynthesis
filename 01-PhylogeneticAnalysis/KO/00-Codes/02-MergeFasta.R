rm(list = ls())
message("Running 02-MergeFasta.R !")

if (!exists("read_fasta", mode="function")) source("00-Functions.R")
if (!require(stringr)) install.packages("stringr")
library(stringr)

# outdir<-"../02-Output/20220606_114002/"

#Retrieve Species information
inFile.ID2Src.Genome<-list.files(path = Genomedir,pattern = "ID2Source",full.names = T)
inFile.ID2Src.OneKP<-list.files(path = OneKPdir,pattern = "ID2Source",full.names = T)
inFile.ID2Src.Dong<-list.files(path = Dongdir,pattern = "ID2Source",full.names = T)

InputCopydir<-file.path(outdir,"01-Input")
if (!dir.exists(InputCopydir)) dir.create(InputCopydir,recursive = T)
file.copy(inFile.ID2Src.Genome,InputCopydir)
file.copy(inFile.ID2Src.OneKP,InputCopydir)
file.copy(inFile.ID2Src.Dong,InputCopydir)
ID2Src.Genome<-fread(inFile.ID2Src.Genome,sep = "\t")
ID2Src.OneKP<-fread(inFile.ID2Src.OneKP,sep = "\t")
ID2Src.Dong<-fread(inFile.ID2Src.Dong,sep = "\t")
ID2Src.all<-rbind(ID2Src.Genome[,.(ID=pID,Clade,Species,Database="Genome")],
                  ID2Src.OneKP[,.(ID=newID1,Clade,Species,Database="OneKP")],
                  ID2Src.Dong[,.(ID=pID,Clade="Liverworts",Species,Database="Dong")]
                  )


#Read fasta files

# list.man<-c("../01-Input/Subtree_seqs-manual.fasta")
list.man<-c()
OneKP.hit<-file.path(outdir,"01-OneKP-blast_hits.fasta")
Genome.hit<-file.path(outdir,"01-AnnotatedGenomes-blast_hits.fasta")
Dong.hit<-file.path(outdir,"01-Dong-blast_hits.fasta")

list.hit<-c(Genome.hit,OneKP.hit,Dong.hit)
fasta<-data.table()
for (f in c(list.man,list.hit)){
  message(paste0("Reading ",f))
  temp<-read_fasta(f)
  temp$Source<-f
  fasta<-rbind(fasta,temp)
}
fasta[,Species:=ID2Src.all$Species[match(ID,ID2Src.all$ID)]]
fasta[,Clade:=ID2Src.all$Clade[match(ID,ID2Src.all$ID)]]
fasta[,Database:=ID2Src.all$Database[match(ID,ID2Src.all$ID)]]

#Substitute some Clade names
fasta[Clade %in% c("Liverwort","Liverworts"),Clade := "Liverwort"]
fasta[Clade %in% c("Moss","Mosses"),Clade := "Moss"]
fasta[Clade %in% c("Hornwort","Hornworts"),Clade := "Hornwort"]
fasta[Clade %in% c("Lycophytes","Lycophyte"),Clade := "Lycophyte"]
fasta[grepl("Monilophytes", Clade), Clade:="Fern"]
fasta[Clade %in% c("Gymnosperms"),Clade := "Gymnosperm"]
fasta[Clade %in% c("Basal Angiosperm","Magnoliids","Eudicots","Monocot"),Clade := "Angiosperm"]


manuals<-unique(fasta[Source %in% list.man,ID])
list.longID<-data.table()
merge_ID<-function(names,by){
  # names<-fasta.uni[1,Name]
  # by=";"
  namelist<-data.table(Name=unique(unlist(strsplit(names,by,fixed = T))))
  namelist[,`:=`(Man=(Name %in% manuals),
                 IsCYP=grepl("CYP",Name,fixed = T),
                 Len=nchar(Name)
                 ),]
  namelist<-namelist[order(-Man,IsCYP,Len,Name)]
  newname<-paste(namelist$Name,collapse = by)
  return(newname)
}
#Merge duplicated sequence
# fasta.uni<-fasta[,.(Name=),by=Sequence]
fasta.uni<-fasta[,.(Name=merge_ID(paste0(ID,collapse = ";"),by=";"),
                    Clade=merge_ID(paste0(Clade,collapse = ";"),by=";"),
                    Species=merge_ID(paste0(Species,collapse = ";"),by=";"),
                    Database=merge_ID(paste0(Database,collapse = ";"),by=";")
                        ),by=Sequence]

#Mark the fragmental sequences (length<300)
fasta.uni[nchar(Sequence)<300,`:=`(Frag=TRUE,
                                   Name=paste0(Name,".f"))]
fasta.uni[is.na(Frag),Frag:=FALSE]

#Assign numbers to replace the IDs
fasta.uni[,`:=`(Num=str_pad(1:nrow(fasta.uni), 5, pad = "0"))]
#Mark full-length and fragmental sequences
fasta.uni[Frag==TRUE,Num:=paste0(Num,"S")] #S for short
fasta.uni[Frag==FALSE,Num:=paste0(Num,"L")] #L for long

#Combine numbers with the sequence
fasta.uni[,Text:=paste0(">",Num,"\n",Sequence)]

# if (nchar(newname)>50) {
#   newname1<-namelist$Name[1]
#   list.longID<<-rbind(list.longID,data.table(ID_long=newname,ID_short=newname1))
#   newname<-newname1
# }

# #Retrieve sources of sequences
# IDSrc<-fread(paste0(outDir,"01-Genomes-ID2Source.txt"))
# IDcompare<-merge(fasta.pri[,.(ID,Source)],IDSrc,by="ID",all.x = T)
# IDcompare[is.na(Source.y),Source.y:=Source.x]
# IDcompare[is.na(Gene),Gene:=ID]
# # IDcompare[,Gene:=gsub("evm.model.","",Gene,fixed = T)]
# setnames(IDcompare, "Source.y", "Source");

# stop()
write(fasta.uni$Text,paste0(outdir,"02-Candidates-all.fasta"))
write(fasta.uni[!Frag==TRUE,Text],paste0(outdir,"02-Candidates-all-fl.fasta"))
write(fasta.uni[Frag==TRUE,Text],paste0(outdir,"02-Candidates-all-frag.fasta"))
write(fasta.uni[Database=="Genome",Text],paste0(outdir,"02-Candidates-Anotated.fasta"))

write.table(fasta.uni[,.(Num,Name,Clade,Species,Database)],paste0(outdir,"02-Candidates-all.id"),col.names = T,row.names = F,quote = F,sep = "\t")

message("Fastas merged!")