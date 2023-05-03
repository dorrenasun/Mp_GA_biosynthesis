rm(list = ls())

if (!exists("read_fasta", mode="function")) source("00-Functions.R")

#Download sequence files
inDir<-"./pep2/"
if (!dir.exists(inDir)){
  tarfile<-"pep2.tar.gz"
  if (!file.exists(tarfile)) {
    old_timeout<-getOption("timeout")
    options(timeout = 600)
    tryCatch(
    download.file("https://figshare.com/ndownloader/files/27933150",tarfile,mode = "wb", quiet = FALSE),
    error = function(e){
      message(paste(tarfile, 'download failed!'))
      file.remove(tarfile)
    } 
    
    )
    options(timeout = old_timeout)
    }
  untar(tarfile)
}

list.file<-list.files(inDir,full.names = T)

{ #Subsetting
  if(!require(readxl)) install.packages("readxl")
  library(readxl)
  
  fSource<-"./List_of_Dong42.xlsx"
  list_src<-read_excel(fSource,sheet = 1,range = cell_cols("A:I"),col_names = T)
  subset<-list_src$`File Name`[!is.na(list_src$Include) & list_src$Include=="T"]
  list.file<-paste0(inDir,unique(subset))
}

refs<-c("../02-OneKP/Extracted/OneKP-Nonseed-merged.fasta",
        "../01-Genomes/Combined/AllGenomes-merged.fasta"
        )
fasta.ref<-data.table()
for (f in refs){
  if (file.exists(f)) {
    ref.temp<-read_fasta(f)
    fasta.ref<-rbind(fasta.ref, ref.temp)
  }
}



fasta.all<-data.table()
nonexist<-c()
for (f in list.file){
  if (!file.exists(f)){
    warning(paste0(f,"does not exist!"))
    nonexist<-c(nonexist,f)
    next
  } 
  temp<-read_fasta(f)
  temp[,`:=`(Source=basename(f),
             Cover=(Sequence %in% fasta.ref$Sequence)
             )]
  fasta.all<-rbind(fasta.all,temp)
}

#Check if these files are actually from OneKP
Cover<-fasta.all[,.(Covered=sum(Cover),Total=.N),by=Source]
Cover[Covered==Total,IsOneKP:="Yes"]
Cover[!Covered==Total,IsOneKP:="No"]



#Merge the rest
list.ori<-Cover[IsOneKP=="No",Source]
fasta.ori<-fasta.all[Source %in% list.ori,]

#Check if there are redundant IDs
if (sum(duplicated(fasta.ori$ID))>0) warning("Duplicated IDs! Wrong parsing!")

#Remove the "TRINITY" prefix if any
list.tri<-unique(fasta.all[grep("\\+TRINITY",ID),Source])
if (length(list.tri)>0) {
  fasta.ori[Source %in% list.tri,pID:=sub("\\+TRINITY","",ID)]
  message(paste0("Removed +TRINITY for sequences from:\n",paste(list.tri,collapse = "\n"))) 
  message("\n")
}

#Attach Species name
fasta.ori[,Species:=list_src$Taxa[match(Source,list_src$`File Name`)]]

#Check sequence integrity
fasta.ori[,`:=`(Start=substr(Sequence,1,1),
                End=substr(Sequence,nchar(Sequence),nchar(Sequence))
                )]

fasta.ori[,`:=`(N_comp=(Start=="M"),
                C_comp=(End=="*"))]
#Remove ending "*"
fasta.ori[C_comp==T,Sequence:=substr(Sequence,1,nchar(Sequence)-1)]

#Replace with new IDs
fasta.ori[is.na(pID),pID:=ID]
fasta.ori[,Text:=paste0(">",pID,"\n",Sequence)]


#Write files
outdir<-"./Combined/"
if (!dir.exists(outdir)) dir.create(outdir)
outbase<-"Dong42"
list.old<-list.files(outdir,full.names = T,recursive = T)
if (sum(grep(outbase,list.old))>0) file.remove(list.old[grep(outbase,list.old)])

#Write the table for OneKP checking
write.table(Cover,paste0(outdir,outbase,"-IsOneKP.tsv"),col.names = T,row.names = F,quote = F,sep = "\t")
#Write fasta
write(fasta.ori$Text,paste0(outdir,outbase,"-merged.fasta"))
#Write ID Sources
write.table(fasta.ori[,.(pID,Source,Species,N_comp,C_comp)],paste0(outdir,outbase,"-ID2Source.txt"),col.names = T,row.names = F,quote = F,sep = "\t")

message("Transcriptomes merged!")
