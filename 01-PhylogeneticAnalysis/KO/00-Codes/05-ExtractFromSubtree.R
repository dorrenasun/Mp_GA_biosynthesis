rm(list=ls())
message("Running 05-ExtractFromSubtree.R !")

if (!exists("read_fasta", mode="function")) source("00-Functions.R")

# outdir<-"../02-Output/20230402_121150/"

subtree<-file.path(outdir,"05-Subtree.nex")
if (!file.exists(subtree)) stop("Subtree file not found! Exit!")
data<-fread(subtree,header=F,quote="",sep = "\t",fill = T)
Line_begin<-grep("begin trees;",data$V1)
Line_end<-grep("end;",data$V1)
trees<-data$V1[(Line_begin+1):(Line_end-1)]
all_seqs<-data.table()
for (entry in trees){
  #Remove all annotations from FigTree
  clean1<-gsub("\\[.*?\\]","",entry)
  clean2<-sub(".*=","",clean1)
  #Remove the scores and quotes
  clean3<-gsub(":[0-9\\.E\\-]+","",clean2)
  #Remove ending semicolon
  if (substr(clean3, nchar(clean3), nchar(clean3))==";"){
    clean3<-substr(clean3, 1, nchar(clean3)-1)
  }
  # clean4<-gsub(" |'|\\(|\\)|","",clean3)
  seqs<-data.table(ID_comb=unlist(strsplit(clean3,",")))
  #Remove leading and trailing spaces
  seqs[,ID_comb:=trimws(ID_comb,"both")]
  #Remove leading and trailing brackets
  seqs[,ID_comb:=trimws(ID_comb,"both",whitespace = "[\\(\\)]")]
  #Remove leading and trailing quotes
  seqs[,ID_comb:=trimws(ID_comb,"both",whitespace = "'")]
  
  seqs[,No:=1:.N]
  #Split the IDs by semicolons
  seqs_sep<-seqs[,.(No,ID_sep=unlist(strsplit(ID_comb,"[;\\|]"))),by=ID_comb]
  #Remove quotation marks
  seqs_sep[,ID_sep:=gsub("\\'","",ID_sep)]
  setkeyv(seqs,"No")
  all_seqs<-rbind(all_seqs,seqs_sep)
}



#Add outgroup sequences from file
outgroup<-data.frame()
try(
  outgroup<-read.table("../01-Input/Outgroup.txt",header = F,sep = "\t",strip.white = T,comment.char = "#"),
  silent = TRUE
)
if (nrow(outgroup)>0){
  outgroup<-as.data.table(outgroup)
  setnames(outgroup,"V1","ID_comb")
  
# outgroup<-rbind(outgroup,data.table(ID_comb=tsv.out))
outgroup_sep<-outgroup[!grepl("^#",ID_comb),.(ID_sep=unlist(strsplit(ID_comb,"[;\\|]"))),by=ID_comb]
outgroup_sep[,No:=0]

all_seqs<-rbind(all_seqs,outgroup_sep)
}

#Remove outliers
exlude<-data.frame()
try(
  exclude<-read.table("../01-Input/Exclude.txt",header = F,sep = "\t",comment.char = "#")
  
)
if (exists("exclude") & nrow(exclude)>0){
  exclude<-as.data.table(exclude)
  if (sum(!exclude$V1 %in% all_seqs$ID_comb)) {
    message("Sequences to exclude not found:")
    message(paste0(exclude[!V1 %in% all_seqs$ID_comb,V1],"\n"))
    }
  all_seqs<-all_seqs[!ID_comb %in% exclude$V1,]
  
}


#Read the source fasta and id files
fasta.file<-paste0(outdir,"02-Candidates-all.fasta")
id.file<-sub(".fasta",".id",fasta.file,fixed = T)
for (f in c(fasta.file,id.file)){
  if (!file.exists(f)) stop("Source fasta file not found! Exit!")
}
fasta<-read_fasta(fasta.file)
id<-fread(id.file)

id.sub<-id[,.(IDnew=unlist(strsplit(Name,"[;\\|]")),Num),by=Name]
id.sele<-merge(id.sub,all_seqs,by.x="IDnew",by.y="ID_sep")

#Check if all wanted sequences were covered
coverage<-all_seqs[,.(cover=sum(ID_sep %in% id.sele$IDnew)),by=ID_comb]
if (sum(coverage$cover<1)>0){
  message("Some sequences not found!")
  message(paste0(coverage[cover==0,ID_comb],"\n"))
}


#Select the needed fasta
fasta.sub<-fasta[ID %in% id.sele$Num,]
#Translate the name to human-readable form
fasta.sub[,Name:=id.sele[match(fasta.sub$ID,Num),Name]]
fasta.sub[,No:=id.sele[match(fasta.sub$ID,Num),No]]
fasta.sub[,Textnew:=paste0(">",Name,"\n",Sequence)]

setkeyv(fasta.sub,"No")


write(fasta.sub$Text,paste0(outdir,"05-Subtree.fasta"))
write(fasta.sub$Textnew,paste0(outdir,"05-Subtree_for_human.fasta"))
write(fasta.sub$Name,paste0(outdir,"05-Subtree.ids"))

#Split files for InterproScan
# MaxSeg<-100
# fasta.sub[,Count:=1:.N]
# for (i in 1:ceiling(nrow(fasta.sub)/MaxSeg)){
#   start<-(i-1)*MaxSeg+1
#   end<-min(max(fasta.sub$Count),MaxSeg*i)
#   fasta.seg<-fasta.sub[Count>=start & Count<=end,]
#   outfile.seg<-file.path(outdir,paste0("05-Subtree-seg",i,".fasta"))
#     write(fasta.seg$Text,outfile.seg)
# }

message("Sequences extracted from subtree!")
