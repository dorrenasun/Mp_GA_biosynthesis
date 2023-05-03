rm(list = ls())
message("Running 03-MaskAlignment.R !")
if (!exists("read_fasta", mode="function")) source("00-Functions.R")
library(stringr)

Repeat=TRUE
ratio=0.98


mask_alignment<-function(file,ratio){
  message(paste0("Processing ",file, " ..."))
  alignment<-read_fasta(file)
  #Check if the original alignment have equal length
  alignment[,len:=nchar(Sequence)]
  if (length(unique(alignment$len))>1) warning("Alignment sequence lengths not equal!")
  seqlen<-max(alignment$len)
  
  aln.seq<-alignment[,.(ID,Sequence)]
  for (i in 1:seqlen){
    lab<-paste0(i,"/",seqlen)
    col<-substr(aln.seq$Sequence,i,i)
    rat<-sum(col=="-")/length(col)
    if (rat<ratio){
      lab<-paste0(lab,"-con")
      aln.seq[,new:=col]
      setnames(aln.seq,"new",paste0("X",i))
      # print(lab)
    }

  }
  
  keepcol<-setdiff(names(aln.seq),c("ID","Sequence"))
  aln.seq.keep<-aln.seq[,.SD,.SDcols=c("ID",keepcol)]
  aln.seq.keep[,newSeq:=paste(.SD,collapse = ""),by="ID",.SDcols=keepcol]
  
  aln.seq.new<-aln.seq.keep[,.(ID,newSeq)]
  #Shorten "_seed_" from seq names
  aln.seq.new[grepl("_seed_",ID,fixed = T),`:=`(FL=T,ID=gsub("_seed_","",ID,fixed = T))]
  aln.seq.new[is.na(FL),FL:=F]
  
  #Check if there are empty sequences
  aln.seq.new[,nGap:=lengths(regmatches(newSeq, gregexpr("-", newSeq)))]
  list.puregaps<-aln.seq.new[nGap==nchar(newSeq),ID]
  if (length(list.puregaps)>0){
    message("The following seqs contain only gaps after masking:")
    print(list.puregaps)
  }
  
  #Merge identical sequences after masking
  aln.seq.uniq<-aln.seq.new[!ID %in% list.puregaps,.(ID=paste0(ID,collapse = "|"),
                               FL=sum(FL)
                              ),by=.(newSeq)]
  
  #Check for empty sequences
  aln.seq.uniq[,`:=`(Text=paste0(">",ID,"\n",newSeq),
                     FL=NULL)]
  return(aln.seq.uniq)
}

list.aln<-list.files(outdir,pattern = "(EINSI|FFTNS2|LINSI)\\.fasta$",full.names = T)
list.files(outdir)

for (f in list.aln){
  #output filebase
  file.base<-paste0(outdir,"/",sub(".fasta$","",basename(f)),".",ratio*100)
  file.fasta<-paste0(file.base,".fasta")
  file.phy<-paste0(file.base,".phy")
  
  if ((!Repeat==T) &file.exists(file.phy) & file.exists(file.fasta)){
    message(paste0(file.base,".* files exist. Skipping!"))
    next
  }
  #Mask the alignments
  keep50<-mask_alignment(f,ratio)
  

  #write fasta file
  write(keep50$Text,file.fasta)
  
  #write to phylip file
  header<-paste0(" ",nrow(keep50)," ",unique(nchar(keep50$newSeq)))
  write(header,file.phy)
  
  keep50[,`:=`(Num=str_pad(1:nrow(keep50), 6, pad = "0"))]
  keep50[,phy:=paste0(ID,"    ",newSeq)]
  write(keep50$phy,file.phy,append = T)
  # write.table(keep50[,.(Num,ID)],paste0(file.base,".id"),quote=F,sep="\t",row.names = F)
  
}
message("Alignments masked!")
