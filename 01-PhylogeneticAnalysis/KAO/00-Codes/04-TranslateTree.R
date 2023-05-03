rm(list=ls())
message("Running 04-TranslateTree.R !")

if (!exists("read_fasta", mode="function")) source("00-Functions.R")

# msg<-"
# This is a script to parse the output files from http://gs.bs.s.u-tokyo.ac.jp/ to standard newick format.
# Usage:
#   parseGS(newick,OTUList,output=NULL)
#   
# Newick and OTUList files could be downloaded from the results page. If you wish to change the name of proteins, modify the second column in OTUList.
# If no OTUList file was provided, the input files will still be parsed, but without changing the IDs.
# Specify the name of parsed file with the \"output\" argument. If not provided, the input newick file name will be used with a \"Parsed-\" prefix.  
# Last modified: 2020-06-12 by Rui Sun
# "
# message(msg)
parseGS<-function(newick,OTUList=NULL,output=NULL){
  # newick<-list
  # OTUList<-paste0(outdir,"02-Candidates-all.id")
  
  if (!is.null(OTUList)){
    dict<-fread(OTUList,header=T,colClasses = "character")
    setnames(dict,"Name","ID")
    testlen<-Inf
    # suffix<-paste0(rep("T",testlen-max(nchar(dict$Num))),collapse = "")
    dict[nchar(ID)>testlen,ID:=substr(ID,1,testlen)]
    dict[duplicated(ID),ID:=Num]
    dict[,ID:=gsub(";","|",ID)]
    dict[,ID:=gsub("[\\(\\)]","",ID)]
    # dict[,ID:=paste0(suffix,Num)]
    #dict[!substr(ID,1,1)=="'",ID:=paste0("'",ID,"'")]
  }
  
  data<-fread(newick,header=F,quote="",sep = "\t",stringsAsFactors = F,colClasses = "character")
  names(data)<-"Original"
  data[,No:=1:.N]
  
  # data$Split<-lapply(data$Original,function(x)(strsplit(x,",")))
  myReplace<-function(entry){
    # entry<-data$Original
    temp1<-gsub("([\\(:\\);\\|/,]+)","@\\1@",entry)
    temp2<-data.table(Split1=unlist(strsplit(temp1,"@",fixed = T)))
    temp2[Split1 %in% dict$Num,Replace:=dict[match(Split1,Num),ID]]
    temp2[is.na(Replace),Replace:=Split1]
    #Check if any problem with splitting
    mimic<-paste0(temp2$Split1,collapse = "")
    if (!mimic==entry) message("Something wrong with splitting!")
    out<-paste0(temp2$Replace,collapse = "")
    return(out)
  }
  data[,Replace:=myReplace(Original),by=No]
  if (is.null(output)) {
    output<-paste0(newick,".newick")
  }
  for (i in data$No){
    if (i>1) output<-paste0(output,".",i+1)
    write.table(data[No==i,Replace],output,quote = F,row.names = F,col.names = F) 
  }
  return(data)
}
# outdir<-"../02-Output/20220610_120945/"
list<-list.files(outdir,pattern="\\.(contree|treefile|support)$",full.names = T,recursive = F)

OTUList<-paste0(outdir,"02-Candidates-all.id")
for (file in list){
  
  outfile<-paste0(file,".newick")
  # if (file.exists(outfile)){
  #   message(paste0(outfile, " exists. Skipping!"))
  # }else{
    message(paste0("Processing ",file))
    result<-parseGS(file,OTUList,outfile)
  # }
  
  }

