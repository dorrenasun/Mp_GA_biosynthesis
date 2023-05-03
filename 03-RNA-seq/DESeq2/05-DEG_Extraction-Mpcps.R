rm(list = ls())
#################### Environment ####################
source("00-Common.R")

pkg_list<-c("data.table","topGO","GO.db","ggh4x")
for (p in pkg_list) {check_package(p)}

#################Parameters#######################
inDir<-file.path(path.out,"01-DESeq2-Mpcps")
# topN<-1000
# Ont<-"BP"

outPath<-file.path(path.out,"05-Heatmaps-Mpcps")
if(dir.exists(outPath)) unlink(sub("/$","",outPath),recursive = T)
if(!dir.exists(outPath))dir.create(outPath)

theme_hm<-theme_pub+theme(
  #Remove border
  panel.border = element_blank(),
  #Remove titles
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  #Remove x&y ticks
  axis.ticks = element_blank(),
  axis.ticks.length.x =unit(0, "pt"),
  axis.ticks.length.y =unit(5, "pt"),
  #Lengend position
  legend.position = "right",
  legend.key.size = unit(0.2, "cm"),
  #Format x and y texts
  axis.text.x.bottom =element_text(size = fontsize, angle = 45, hjust = 1,color = "grey50"),
  axis.text.x.top =element_text(size = fontsize, angle = 90, hjust = 0,vjust=0.5, colour = "black"),
  axis.text.y.left=element_text(face = "italic",size = fontsize, colour = "black"),
  axis.text.y.right=element_text(size = fontsize)
)
#################### Other Functions ####################
toVer6<-function(dt,to_translate="tracking_id"){
  # dt<-dt.temp
  # to_translate<-"tracking_id"
  ori_col<-names(dt)
  dt$myquery<-dt[,..to_translate]
  #Translate to ver 6
  dt[,`:=`(myrank=1:.N,gene=sub("\\..*","",myquery))]
  dt[gene=="Mp4g0958", gene:="Mp4g09580"]
  dt.v6<-merge(dt,MpSyn,by.x="gene",by.y="query",all.x=T)
  dt.v6[MpTak_v6.1=="-",MpTak_v6.1:=NA]
  setorder(dt.v6,myrank)
  list.nonver6<-dt.v6[is.na(MpTak_v6.1),..to_translate]
  if (length(list.nonver6)>0){
    message("The following GOI(s) does not exist in ver6:")
    message(paste(list.nonver6,collapse = "\n"))
  }
  out_cols<-c(ori_col,"MpTak_v6.1")
  dt.out<-dt.v6[!is.na(MpTak_v6.1)&!duplicated(MpTak_v6.1),..out_cols]
  return(dt.out)
}
myHM<-function(data){
  #Theme for barplot

  # data<-data.sub[gene %in% gene[sig]]
  # data<-data.sub
  setorder(data,Gene,Sample)
  fir.labels<-unique(data$Gene)
  sec.labels<-data[match(fir.labels,Gene),name]
  sec.labels.dt<-data.table(Pos=1:length(sec.labels),Label=sec.labels)
  sec.labels.dt1<-sec.labels.dt[,.(Pos=mean(Pos),
                                   Start=min(Pos)-0.35,
                                   End=max(Pos)+0.35,
                                   n=length(Pos)
  ),by=Label]
  sec.labels.dt1[n>10,Label:=paste0(Label,"(",n,")")]
  sec.labels.dt1[,Height:=length(unique(data$Sample))+0.6]
  
  data[,`:=`(y=as.integer(factor(Gene,levels = fir.labels)),
             isSig=factor(isSig)
             )]
  # data.FRres<-unique(data[FR.Res>0,c("Gene","FR.Res","y")])
  p<-ggplot(data, aes(x=y,y=Sample)) +
    geom_tile(aes(fill = log2FoldChange),color="black",width=0.9, height=0.9, linewidth=0.3*Lwd) +
    geom_point(data = data[!isSig==0,],aes(x=y,y=Sample,color=isSig),shape=20,size=Lwd)+
    geom_segment(data=sec.labels.dt1,aes(x=Start,y = Height,xend=End,yend=Height),color="black",linewidth=0.5*Lwd)+
    geom_text(data=sec.labels.dt1,aes(x=Pos,y = Height+0.1,label=Label),color="black",size=fontsize*Lwd,angle=90,hjust=0)+

    scale_fill_distiller(breaks = c(-4,0,4), limits=c(-4.2,4.2),palette = "RdBu",oob=squish)+
    theme_hm+theme(aspect.ratio = (length(unique(data$Sample))+0.2)/length(unique(data$Gene)))
  
  p<-p+scale_y_discrete(expand = expansion(add = c(0, 0.7)))+
    scale_x_continuous(limits=c(0.5, length(fir.labels)+0.5),
                       expand = expansion(),
                       breaks = 1:length(fir.labels),labels = fir.labels#,
                       # sec.axis = sec_axis(~.,breaks = sec.labels.dt1$Pos,labels = sec.labels.dt1$Label)
                       )+
    scale_color_manual(breaks=c(0,1,2),values = c(NA,"grey70","black"))+
    force_panelsizes(rows = unit(1.25, "cm"),
                     cols = unit(0.25*length(fir.labels), "cm"))
  
  return(p)
}



################# Main #######################
#Load data from DESeq2
list.files<-c(
  file.path(inDir, "DEG-FR_vs_cW_Tak1.tsv"),
  file.path(inDir, "DEG-FR_vs_cW_MpcpsMock.tsv"),
  file.path(inDir, "DEG-MpcpsMock_vs_Tak1_FR.tsv"),
  file.path(inDir, "DEG-MpcpsKA_vs_Tak1_FR.tsv"),
  file.path(inDir, "DEG-Mpcps_KAvsMock_FR.tsv")
)
tag_all<-sub(".*DEG-(.*)\\.tsv","\\1",list.files)
  
DEG.all<-data.table()
for (f in list.files){
  DEG.temp<-fread(f)
  DEG.temp[,Sample:=sub(".*DEG-(.*)\\.tsv","\\1",f)]
  DEG.all<-rbind(DEG.all,DEG.temp[,.(Gene,log2FoldChange,padj,Sample)])
  
}
DEG.all[,Sample:=factor(Sample,levels = rev(unique(Sample)))]
DEG.all[,isSig:=as.integer(padj<pthre)]
DEG.all[abs(log2FoldChange)>FCthre,isSig:=2*isSig]
DEG.all[,ofInt:=sum(isSig[Sample %in% tag_all[c(3)]]>1)>0,by=.(Gene)]

DEG.all[,cps:=(log2FoldChange[Sample==tag_all[3]]>0),by=.(Gene)]

#Some manual GOIs
list.GOI<-fread("./list-GOI.csv",header = T)

order.pathway<-c("Flavonoid","ABA","GA","Cytokinin")

list.GOI.v6<-toVer6(list.GOI[pathway %in% order.pathway],to_translate = "ver6")
#filter out genes not exist in ver6
list.GOI.v6[,rank:=1:.N]
setorder(list.GOI.v6,rank,MpTak_v6.1)
order.name<-unique(list.GOI.v6$name)
order.gene<-unique(list.GOI.v6$MpTak_v6.1)


#Merge list of GOI with expression data
#Warning: genes that has no expression data will be omitted
DEG.sel<-merge(DEG.all[Gene %in% order.gene],
               list.GOI.v6[,.(pathway,name,MP_GENES,Gene=MpTak_v6.1)],
               by="Gene")
list.missing<-list.GOI.v6[!MpTak_v6.1 %in% DEG.sel$Gene,paste0(name,"(",MpTak_v6.1,")")]
if (length(list.missing)>0){
  message("The following gene(s) has no expression data and is discarded:")
  message(paste(list.missing,collapse = "\n"))
}

#Fix order
DEG.sel[,`:=`(Gene=factor(Gene,levels = order.gene),
              name=factor(name,levels = order.name)
)]

setorder(DEG.sel,name,Gene)

#Attach gene names
# # DEG.sel[!MP_GENES=="",Gene:=paste0(Gene,"(",MP_GENES,")")]
# DEG.sel[,`:=`(Gene=factor(Gene,levels = unique(Gene))
#               )]
for (p in order.pathway){
  DEG.sub<-DEG.sel[pathway==p,]
  #Check if significant genes are more thant one
  DEG.sig<-DEG.sub[ofInt==T,]
  if (nrow(DEG.sig)>0){
    p.hm.sig<-myHM(DEG.sub[ofInt==T,])
    print(p.hm.sig)
    
    ggsave(filename = file.path(outPath,paste0("Heatmap_Mpcps-",p,"-Sig.pdf")),
           plot = p.hm.sig,
           height = 5,
           width=30,
           units = "cm"
    )
    
  }else{
    message(paste0("No significant DEGs in pathway: ",p))
  }
  #If it is a small list, generate the plot as it is
  if (length(unique(DEG.sub$Gene))<30){
    p.hm<-myHM(DEG.sub)
    ggsave(filename = file.path(outPath,paste0("Heatmap_Mpcps-",p,".pdf")),
           plot = p.hm,
           height = 10,
           width=50,
           units = "cm"
    )
  }
    
}

#Large lists
file.lists<-list(
  CYP="./list-CYP.csv",
  DOXC="./list-DOXC.csv",
  UGT="./list-UGT.csv",
  CellWall="./list-CellWall.csv",
  POD="./list-POD.csv",
  TF="./list-TF.tsv"
)
for (ll in names(file.lists)){
  list.sub<-fread(file.lists[[ll]],header = T)
  DEG.sub<-DEG.all[Gene %in% list.sub[ver6.1!="-",ver6.1] & ofInt==T,]
  DEG.sub[,name:=list.sub[match(Gene,ver6.1),Clade]]
  DEG.sub[,cpsFC:=log2FoldChange[Sample=="MpcpsMock_vs_Tak1_FR"],by=.(Gene)]

  DEG.Down<-DEG.sub[cpsFC<0,]
  if (nrow(DEG.Down)>0){
    DEG.Down[,minFC:=min(log2FoldChange[Sample=="MpcpsMock_vs_Tak1_FR"]),by=.(name)]
    setorder(DEG.Down,minFC,cpsFC)
    # if (ll=="POD") setorder(DEG.Down,Gene)
    DEG.Down[,Gene:=factor(Gene,levels = unique(Gene))]
    p.hm.Down<-myHM(DEG.Down)
    print(p.hm.Down)
    ggsave(filename = file.path(outPath,paste0("Heatmap_Mpcps-",ll,"-Down.pdf")),
           plot = p.hm.Down,
           height = 5,
           width=30,
           units = "cm"
    )
  }else{
    message(paste0("No down-regulated DEGs in pathway: ",ll))
  }
  DEG.Up<-DEG.sub[cpsFC>0,]
  if (nrow(DEG.Up)>0){
    DEG.Up[,maxFC:=max(log2FoldChange[Sample=="MpcpsMock_vs_Tak1_FR"]),by=.(name)]
  setorder(DEG.Up,maxFC,cpsFC)
  # if (ll=="POD") setorder(DEG.Up,Gene)
  DEG.Up[,Gene:=factor(Gene,levels = unique(Gene))]
  p.hm.Up<-myHM(DEG.Up)
  print(p.hm.Up)
  ggsave(filename = file.path(outPath,paste0("Heatmap_Mpcps-",ll,"-Up.pdf")),
         plot = p.hm.Up,
         height = 5,
         width=30,
         units = "cm"
  )
  }else{
    message(paste0("No up-regulated DEGs in pathway: ",ll))
  }

}


