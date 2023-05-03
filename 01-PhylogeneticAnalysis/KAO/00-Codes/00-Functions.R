if (!require(data.table)) install.packages("data.table")
#if (!require(R.utils)) install.packages("R.utils")
library(data.table)
#library(R.utils)

para<-fread("../02-Output/timestamp.txt",header = F)
outdir<-paste0("../02-Output/",para[V1=="runID",V2],"/")
if (!dir.exists(outdir)) dir.create(outdir)
OneKPdir<-dirname(para[V1=="Path_OneKP",V2])
Genomedir<-dirname(para[V1=="Path_Genome",V2])
Dongdir<-dirname(para[V1=="Path_Dong",V2])


#Read all fasta files
read_fasta<-function(file){
  if (!file.exists(file)) stop("Fasta file not found!")
  Fasta<-fread(file,header = F,sep="\t",col.names = "Full")
  ID.positions<-grep(">",Fasta$Full,fixed = T)
  ID.lengths<-c(ID.positions[-1],nrow(Fasta)+1)-ID.positions
  IsIsoetes<-sum(grepl("Isoetes",Fasta$Full))
  if (IsIsoetes>0){
    Fasta[,`:=`(Num=rep(1:length(ID.positions),ID.lengths),
                       ID=rep(sub(".*(Itaiw_v1_.*_[0-9]+_.*R.).*","\\1",Full[ID.positions]),ID.lengths))]
  }else{
    Fasta[,`:=`(Num=rep(1:length(ID.positions),ID.lengths),
                ID=rep(sub(" .*","",Full[ID.positions]),ID.lengths))]
  }
  
  Fasta[,`:=`(ID=sub(">","",ID,fixed = T),
              Sequence=Full)]
  Fasta[ID.positions,Sequence:=""]
  Fasta.trans<-Fasta[,.(Sequence=paste(Sequence,collapse = "")),by=.(Num,ID)]
  Fasta.trans[,`:=`(ID=gsub("|","+",ID,fixed = T))]
  Fasta.trans[,ID:=sub("[ |\t].*","",ID)]
  Fasta.trans[,Text:=paste(paste0(">",ID),Sequence,sep = "\n")]
  
  return(Fasta.trans)
}

get_primary<-function(fasta){
  
  list.IDSrc<-list.files(Genomedir,"ID2Source.txt",full.names = T)
  IDSrc<-data.table()
  for (l in list.IDSrc){
    temp<-fread(l)
    IDSrc<-rbind(IDSrc,temp)
  }
  names(IDSrc)<-c("ID","Gene","var","IDSrc")
  fasta1<-copy(fasta)

  fasta1<-merge(fasta1,IDSrc,by="ID",all.x = T)
  #Make sure the variant is integer
  fasta1[,var:=as.integer(var)]
  fasta1[is.na(Gene),`:=`(Gene=ID,var=0)]
  #Convert variants to integers
  fasta1<-fasta1[order(Gene, var)]
  
  #Extract primary transcripts
  fasta.pri<-fasta1[!duplicated(Gene),]
  fasta.pri[,c("var","IDSrc"):=NULL]
  
  return(fasta.pri)
}

read_hmm<-function(hmmfile){
  hmm<-read.table(hmmfile,header = F,comment.char = "#")
  names(hmm)<-c("Target_name","Target_Acc","Query_name","Query_acc",
                "Evalue.Full","Score.Full","Bias.Full",
                "Evalue.Best1","Score.Best1","Bias.Best1",
                "Num_exp","Num_reg","Num_clu","Num_ov","Num_env","Num_dom","Num_rep","Num_inc",
                "Description"
  )
  return(as.data.table(hmm))
  
}
read_hmm_dom<-function(hmmdomfile){
  hmmdom<-read.table(hmmdomfile,header = F,comment.char = "#")
  names(hmmdom)<-c("Target_name","Target_Acc","Tlen","Query_name","Query_acc","Qlen",
                "Evalue.Full","Score.Full","Bias.Full",
                "No_Dom","Num_Dom","cEvalue.Dom","iEvalue.Dom","Score.Dom","Bias.Dom",
                "Hmm.From","Hmm.To", "Ali.From","Ali.To","Env.From","Env.To",
                "Accuracy","Description"
  )
  return(as.data.table(hmmdom))
  
}
