message("Loading some common functions and parameters")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",update = T)

check_package<-function(pkg){
  if (!require(pkg,character.only = T,quietly = T)){
    # message(paste(pkg,"is not here. Install."))
    BiocManager::install(pkg,ask=F)
    library(pkg,character.only = T)
  }
}

pkg_list<-c("R.utils","data.table","ggplot2","scales")
for (p in pkg_list) {check_package(p)}


#Threshold for DEG
pthre<-0.01
FCthre<-0.585
path.out<-paste0("./Output_a",pthre,"FC",FCthre)
if (!dir.exists(path.out)) dir.create(path.out,recursive = T)

# FCthre<-1


#Color Pallette for ggplot2
Lwd=0.3514
fontsize<-6
#theme
theme_pub<-theme_bw()+theme(
  aspect.ratio = 1.5,
  # plot.margin = unit(c(2,1,1,1), "cm"),
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour="black",linewidth =1*Lwd),
  text = element_text(size = fontsize, colour = "black"),
  strip.background = element_blank(),
  strip.text=element_text(size = fontsize, colour = "black"),
  # axis.title.x = element_blank(),
  axis.text  = element_text(size = fontsize, colour = "black"),
  # axis.text.x = element_text(family = "sans",size = 6, angle=45,hjust = 1,colour = "black"),
  axis.ticks = element_line(colour="black",linewidth=0.5*Lwd),
  # axis.text.x = element_blank(),
  legend.background = element_blank(),
  legend.key.size=unit(0.5,"cm"),
  legend.text.align = 0#,
  # legend.position = "none"
)



pal <- c("grey30","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d","white")
pie(rep(1,16),col = pal)
#Common functions

read_gff<-function(gff_file,select="mRNA"){
  #Read as whole lines
  gff<-fread(gff_file,sep = "\n",header=F,strip.white = T) 
  #Remove comment lines
  gff.clean<-gff[!substr(V1,1,1)=="#",tstrsplit(V1,"\t")]
  
  gff.sel<-gff.clean[V3 %in% select,.(type=V3,V9,Num=1:.N)]
  gff.sel.l<-gff.sel[,.(item=unlist(strsplit(V9,";"))),by=.(type,Num)]
  gff.sel.l[,`:=`(Tag=sub("=.*","",item),
                  Value=sub(".*=","",item)
  )]
  gff.t<-dcast(gff.sel.l,type+Num~Tag,value.var = "Value")
  return(gff.t[,Num:=NULL])
}

call_prep<-function(){
  path.curr<-getwd()
  source("../Prepare_Annotations.R",chdir = T)
  setwd(path.curr)
}
#Load annotations
path.annot<-"../RefGenome"
if (!dir.exists(path.annot)) call_prep()

#Convert v6 transcript to gene
if (!exists("ids.v6")){
  file.gff.v6<-file.path(path.annot,"MpTak_v6.1r1.gff.gz")
  if (!file.exists(file.gff.v6)) {
    download.file(url = "https://marchantia.info/download/MpTak_v6.1/MpTak_v6.1r1.gff.gz",
                  destfile = file.gff.v6)}
  mRNA.v6<-read_gff(file.gff.v6,select = "mRNA")
  ids.v6<-mRNA.v6[,.(transcript=ID,gene=Parent)]
}

#Retrieve nomenclatures
file.Mpnames<-file.path(path.annot,"MpAnnot_clean.tsv")
if (!file.exists(file.Mpnames)) {
  call_prep()
  }
Mpnames<-fread(file.Mpnames)
#   
file.MpGO<-file.path(path.annot,"MpGenetoGO_clean.tsv")
if (!file.exists(file.MpGO)) {
  call_prep()
}
MpGO<-fread(file.MpGO)

file.MpSyn<-file.path(path.annot,"MpGeneSyn_clean.tsv")
if (!file.exists(file.MpSyn)) {
  call_prep()
}
MpSyn<-fread(file.MpSyn)


#Get Salmon quant files
path.salmon<-file.path("..","Salmon_Tak1")
file.quant<-list.files(path=path.salmon,pattern="\\.sf$",recursive = T,full.names = T)

{#Creat a sampleTable
  sampleTable<-data.table(files=file.quant)
  sampleTable[,names:=sub(".*\\/((Mp|Tak).*-(cW|FR)-[0-9]+).*","\\1",files)]
  sampleTable[,`:=`(Genotype=sub("-(cW|FR).*","",names),
                    Treatment=sub(".*-(cW|FR).*","\\1",names)
  )]
  sampleTable[,Genotype:=sub("-4","",Genotype)]
  sampleTable[,Genotype:=sub("-1","1",Genotype)]
  sampleTable[,Genotype:=sub("-","_",Genotype)]
  
  
  
  #Factorize
  order.gt<-c("Tak1","Mpcps_Mock","Mpcps_KA","Mpdella")
  order.trt<-c("cW","FR")
  
  #define colors by gt
  col.gt<-pal[c(1,7,8,13)]
  names(col.gt)<-order.gt
  
  sampleTable[,`:=`(Genotype=factor(Genotype,levels = order.gt),
                    Treatment=factor(Treatment,levels = order.trt)
  )]
  
  #sort by Treatment and Genotype, then combine them into Cond
  setorder(sampleTable,Treatment,Genotype)
  sampleTable[,Cond:=paste0(Genotype,".",Treatment)]
  sampleTable[,Cond:=factor(Cond,levels=unique(Cond))]
  
}
