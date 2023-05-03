rm(list = ls())
#################### Environment ####################
source("00-Common.R")

pkg_list<-c("data.table","DESeq2","tximeta","ComplexHeatmap","ggrastr",#"readxl",
            "eulerr","ggplot2","scales","ggrepel")
for (p in pkg_list) {check_package(p)}
# if(!require(org.Mpolymorpha.eg.db)) source("00-GO_clustering.R")


if (!require(lasso2)) {
  packageurl <-"https://cran.r-project.org/src/contrib/Archive/lasso2/lasso2_1.2-22.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
}

#################Parameters#######################

outDir<-file.path(path.out,"01-DESeq2-Mpcps")
if (dir.exists(outDir)) unlink(sub("\\/$","",outDir), recursive=TRUE)
if (!dir.exists(outDir)) dir.create(outDir,recursive = T)
SpecialGenes<-c(MpCPS="Mp2g07200",
                MpKS="Mp6g05950",
                MpKOL1="Mp3g18320",
                MpKAOL1="Mp4g23680",
               # MpKAOL2="Mp1g25410",
                MpBNB="Mp3g23300"
                )



################# Functions #######################
make_dds<-function(coldata,design){
  se <- tximeta(coldata,type="salmon",skipMeta=T,txOut=F,tx2gene = tx2gene.Tak1)
  
  # Full interaction model
  dds <- DESeqDataSet(se, design = design)
  dds.keep<-rowSums(counts(dds)) >= 10
  dds.kept <- dds[dds.keep,]
  
  return(dds.kept)
}

myOutput<-function(DESeqDataSet,contrast,tag){
  res.DES<-results(DESeqDataSet, contrast = contrast,alpha = pthre,independentFiltering = F)
  
  res <- suppressWarnings(as.data.table(res.DES,keep.rownames="Gene"))
  res[,Name:=Mpnames[match(Gene,MpTak_v6.1),MPGENES]]
  res[is.na(Name),Name:=""]
  
  cols<-c("Gene","Name",setdiff(names(res),c("Gene","Name")))
  
  res.out<-res[,..cols]
  file<-file.path(outDir,paste0("DEG-",tag,".tsv"))
  write.table(res.out,
              file=file,
              row.names = F,
              col.names = T,
              quote = F,
              sep = "\t")
  res[,Sig:=((!is.na(padj)) & (padj<pthre) & (abs(log2FoldChange)>FCthre))]
  
  res.sig<-res[Sig==TRUE,..cols]
  # res.sig[,Sig:=NULL]
  setkey(res.sig,log2FoldChange)
  file.sig<-file.path(outDir,paste0("DEG-",tag,"-sig.tsv"))
  
  write.table(res.sig,
              file=file.sig,
              row.names = F,
              col.names = T,
              quote = F,
              sep = "\t")
  res.out[,Tag:=tag]
  return(res.out)
}
  
myVolcano<-function(res.dt,tag,spec_genes=SpecialGenes,xlim=12,ylim=200){
  res<-copy(res.dt)
  res[,Tag:=factor(Tag,levels = unique(Tag))]
 # res.dt<-res.MpcpsMock.cW
  # tag<-unique(res.dt$Tag)[1]
  file.volcano<-file.path(outDir,paste0("Volcano-",tag,".pdf"))
  
  #Set categories
  res[,Cat:=as.integer(padj<pthre)]
  res[padj<pthre & abs(log2FoldChange)>FCthre,Cat:=2]
  res[padj<pthre &  log2FoldChange<(-FCthre) & Gene%in%spec_genes,Cat:=3]
  res[padj<pthre &  log2FoldChange>FCthre & Gene%in%spec_genes,Cat:=4]
  
  res.sum<-res[Cat>1,.(Up=paste0("Up:\n",sum(log2FoldChange>FCthre)," genes"),
                       Down=paste0("Down:\n",sum(log2FoldChange<(-FCthre))," genes")
                       ),by=(Tag)]
  
  #labels
  
  
  labels<-melt(res.sum,id.vars = "Tag",variable.name = "Side",value.name = "Label")
  labels[Side=="Up",`:=`(x=0.5*xlim,y=0.7*ylim)]
  labels[Side=="Down",`:=`(x=-0.5*xlim,y=0.7*ylim)]

  #Fix all p>ylim to ylim
  # res[-log10(padj)>ylim,padj:=10^(-ylim-1)]
  
  if (!is.null(names(spec_genes))){
    res[Gene %in% spec_genes,Name:=names(spec_genes)[match(Gene,spec_genes)]]

  }

  
  p.vol<-ggplot(res[!is.na(padj),],aes(x=log2FoldChange,y=-log10(padj)))+
    rasterize(geom_point(data=res[Cat<=2,],aes(fill=as.factor(Cat)),shape=21,size=.8*Lwd,alpha=0.6,stroke=0),dpi=1200)+
    geom_point(data=res[Cat>2,],aes(fill=as.factor(Cat)),shape=21,size=2*Lwd,alpha=1,color="white",stroke=0.013)+
    geom_text(data=labels,aes(x,y,label=Label),hjust=0.5,size=6*Lwd)+
    
    geom_text_repel(data=res[Cat>2,],
                    aes(x=log2FoldChange,
                        y=-log10(padj),
                        label=sub("\\[(.*)\\].*","\\1",Name)#,
                        # color=as.factor(Cat)
                        ),
                    alpha=0.8,
                    color="#ff6600",
                    min.segment.length =0,
                    segment.size=0.2*Lwd,
                    max.overlaps=500,
                    # hjust=0.5,nudge_y = 1,
                    size=5*Lwd)+
    facet_wrap(~Tag,nrow=1,strip.position = "top")+
    theme_pub+theme(aspect.ratio = 1.2,legend.position = "none")+
    scale_fill_manual(values = c(`0`="yellow",`1`="yellow",`2`="black",`3`="blue",`4`="#ff6600"))+
    # scale_color_manual(values = c(`0`="yellow",`1`="yellow",`2`="black",`3`="blue",`4`="#ff6600"))+
    scale_x_continuous(limits = c(-xlim,xlim),oob = squish)+
    scale_y_continuous(limits = c(-.5,ylim+2),oob=squish)
  print(p.vol)
  ggsave(file.volcano,plot = p.vol,width = 20,height = 4.5,units = "cm")
  return(p.vol)
}

myHist<-function(res.dt,tag,spec_genes=SpecialGenes,xlim=10,ylim=800){
  res<-copy(res.dt)
  res[,Tag:=factor(Tag,levels = unique(Tag))]
  # res.dt<-res.MpcpsMock.cW
  # tag<-unique(res.dt$Tag)[1]
  file.hist<-file.path(outDir,paste0("Hist-",tag,".pdf"))
  
  #Set categories
  res.sig<-res[(padj<pthre)&(abs(log2FoldChange)>FCthre),]
  # res[,Cat:=as.integer(padj<pthre)]
  # res[padj<pthre & abs(log2FoldChange)>FCthre,Cat:=2]
  # res[padj<pthre &  log2FoldChange<(-FCthre) & Gene%in%spec_genes,Cat:=3]
  # res[padj<pthre &  log2FoldChange>FCthre & Gene%in%spec_genes,Cat:=4]
  
  res.sum<-res.sig[,.(Up=paste0("Up:\n",sum(log2FoldChange>FCthre)," genes"),
                       Down=paste0("Down:\n",sum(log2FoldChange<(-FCthre))," genes")
  ),by=(Tag)]
  
  
  #labels
  labels<-melt(res.sum,id.vars = "Tag",variable.name = "Side",value.name = "Label")
  labels[Side=="Up",`:=`(x=0.5*xlim,y=0.8*ylim)]
  labels[Side=="Down",`:=`(x=-0.5*xlim,y=0.8*ylim)]
  
  {
    p.hist<-ggplot(res.sig, aes(x=log2FoldChange)) + 
      geom_histogram(color="white", fill="black",bins = 60,linewidth=0.2*Lwd)+
      geom_text(data=labels,aes(x,y,label=Label),hjust=0.5,size=6*Lwd)+
      
      facet_wrap(~Tag,nrow=1)+
      ylab("Number of genes")+
      theme_pub+theme(aspect.ratio = 1.0,legend.position = "none")
    print(p.hist)
    
  }

  ggsave(file.hist,plot = p.hist,width = 20,height = 4,units = "cm")
  return(p.hist)
}


myFisher<-function(A,B,v,total){
  # A<-"MpDELLAox.Up"
  # B<-"pif_up"
  x<-v$original.values[paste0(A,"&",B)]
  m<-sum(v$original.values[grep(A,names(v$original.values))])
  n<-total-m
  k<-sum(v$original.values[grep(B,names(v$original.values))])
  return(phyper(x-1,m,n,k,lower.tail = FALSE))
}


res2mat<-function(res.dt){
  res.sigmat<-res.dt[,.(Gene=Gene,
                     Up=(padj<pthre&(log2FoldChange >FCthre)),
                     Down=(padj<pthre&(log2FoldChange < -FCthre))
  )
  ]
  # setnames(res.sigmat,"Up",)
  if ("Tag" %in% names(res.dt)){
    tag<-unique(res.dt$Tag)[1]
    if (!is.na(tag)){
      names(res.sigmat)[2:3]<-paste(tag,names(res.sigmat)[2:3],sep=".")
    }
  }
  
  return(res.sigmat)
}

res2list<-function(inFile,Tag="Data",rm=NULL){
  DEGdata<-fread(inFile,header = T,sep = "\t")
  Sig.up<-DEGdata[!is.na(padj) & padj<pthre & log2FoldChange>FCthre & (!Gene %in% rm),Gene]
  Sig.down<-DEGdata[!is.na(padj) & padj<pthre & log2FoldChange<(-FCthre) & (!Gene %in% rm),Gene]
  list.out<-list(Sig.up,Sig.down)
  names(list.out)<-paste0(Tag,c(".Up",".Down"))
  return(list.out)
}

################# Main #######################
message("Running DEG analysis...")


#Generate tx2gene list
tx2gene.Tak1<-ids.v6[!grepl("MpUg",gene),]

# data table for making UpSet plot
UpSet.all<-data.table(Gene=unique(tx2gene.Tak1$gene))

# MpGO.ori[,gene:=sub("\\..*","",query)]
# UpSet.all[,Annot:=Gene %in% MpGO$MpTak_v6.1]


{#Build DESeqDataSet & filter by total counts
  dds.cps<-make_dds(sampleTable[!Genotype=="Mpdella",],design = ~ Genotype+Treatment+Genotype:Treatment)

  #log-scale transformations
  rlg.cps<-rlog(dds.cps,blind = F)
  rlg_mat<-assay(rlg.cps)
  rlg_mat.dt<-as.data.table(rlg_mat,keep.rownames = "Gene")

  vsd.cps <- vst(dds.cps, blind=FALSE)
  vsd_mat<-assay(vsd.cps)
  vsd_mat.dt<-as.data.table(vsd_mat,keep.rownames = "Gene")
  write.table(vsd_mat.dt,
              file=file.path(outDir,"VST_cps.tsv"),
              row.names = F,
              col.names = T,
              quote = F,
              sep = "\t")
  
} 

{  #Visualization of sample distances
  pcaData <- plotPCA(vsd.cps, intgroup=c("Genotype","Treatment","names"),ntop=length(vsd.cps), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  {
    p.pca<-ggplot(pcaData, topn=length(vsd.cps),aes(x=PC1,y=PC2, color=Genotype)) +
      geom_point(aes(shape=Treatment),size=4*Lwd,alpha=0.8) +
      geom_text(aes(label=names),size=2,show.legend = F,nudge_x = 1,hjust=0,vjust=0.5)+
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      scale_color_manual(values = col.gt)+
      scale_alpha_manual(values =c(0.6,0.8,1)) +
      scale_shape_manual(values =c(16,17)) +
      coord_fixed() + theme_pub + 
      theme(aspect.ratio = 1,
            legend.position = "right")
    print(p.pca)
    out.pdf<-file.path(outDir,"Sample_distance_Mpcps_VST.pdf")
    ggsave(file=out.pdf, plot=p.pca, width=10, height=10,units = "cm")
  }
}


{  #Pair-wise comparisons
  # Ref:https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html
  test.cps <- DESeq(dds.cps)
  #cW comparisons
  res.MpcpsMock.cW<-myOutput(test.cps,contrast = c("Genotype","Mpcps_Mock","Tak1"),tag="MpcpsMock_vs_Tak1_cW")
  res.MpcpsKA.cW<-myOutput(test.cps, contrast = c("Genotype","Mpcps_KA","Tak1"),tag = "MpcpsKA_vs_Tak1_cW")
  res.KA.cW<-myOutput(test.cps, contrast = c("Genotype","Mpcps_KA","Mpcps_Mock"),tag = "Mpcps_KAvsMock_cW")

  p.cW<-myVolcano(rbind(res.MpcpsMock.cW,res.MpcpsKA.cW,res.KA.cW),tag="cW")
  #FR comparisons
  res.MpcpsMock.FR<-myOutput(test.cps,
                             contrast = list( c("Genotype_Mpcps_Mock_vs_Tak1",
                                                "GenotypeMpcps_Mock.TreatmentFR") ),
                             tag="MpcpsMock_vs_Tak1_FR")
  res.MpcpsKA.FR<-myOutput(test.cps,
                           contrast = list( c("Genotype_Mpcps_KA_vs_Tak1",
                                              "GenotypeMpcps_KA.TreatmentFR") ),
                           tag="MpcpsKA_vs_Tak1_FR")
  
  res.KA.FR<-myOutput(test.cps, contrast = c(0,-1,1,0,-1,1),tag="Mpcps_KAvsMock_FR")
  p.FR<-myVolcano(rbind(res.MpcpsMock.FR,res.MpcpsKA.FR,res.KA.FR),tag="FR")
  
  #FR response comparisons
  res.FR.Tak1<-myOutput(test.cps, contrast = c("Treatment","FR","cW"),tag = "FR_vs_cW_Tak1")
  res.FR.MpcpsMock<-myOutput(test.cps, list( c("Treatment_FR_vs_cW","GenotypeMpcps_Mock.TreatmentFR")),tag = "FR_vs_cW_MpcpsMock")
  res.FR.MpcpsKA<-myOutput(test.cps, list( c("Treatment_FR_vs_cW","GenotypeMpcps_KA.TreatmentFR")),tag = "FR_vs_cW_MpcpsKA")

  res.FRresp<-rbind(res.FR.Tak1,res.FR.MpcpsMock)
  p.FR.resp<-myVolcano(res.FRresp,tag="FRresp")#+theme(aspect.ratio = 0.4)
  p.hist.FRresp<-myHist(res.FRresp,tag="FRresp")
  
  #Check by pairwise loading
  {
    dds.p1<-make_dds(sampleTable[Cond %in% c("Tak1.cW","Mpcps_Mock.cW"),],design = ~ Genotype)
    res.p1<-myOutput(DESeq(dds.p1), contrast = c("Genotype","Mpcps_Mock","Tak1"),tag = "MpcpsMock_vs_Tak1_cW.pw")
    
    dds.p2<-make_dds(sampleTable[Cond %in% c("Tak1.cW","Mpcps_KA.cW"),],design = ~ Genotype)
    res.p2<-myOutput(DESeq(dds.p2), contrast = c("Genotype","Mpcps_KA","Tak1"),tag = "MpcpsKA_vs_Tak1_cW.pw")
    
    dds.p3<-make_dds(sampleTable[Cond %in% c("Mpcps_KA.cW","Mpcps_Mock.cW"),],design = ~ Genotype)
    res.p3<-myOutput(DESeq(dds.p3), contrast = c("Genotype","Mpcps_KA","Mpcps_Mock"),tag = "Mpcps_KAvsMock_cW.pw")
    
    dds.p4<-make_dds(sampleTable[Cond %in% c("Tak1.FR","Mpcps_Mock.FR"),],design = ~ Genotype)
    res.p4<-myOutput(DESeq(dds.p4), contrast = c("Genotype","Mpcps_Mock","Tak1"),tag = "MpcpsMock_vs_Tak1_FR.pw")
    
    dds.p5<-make_dds(sampleTable[Cond %in% c("Tak1.FR","Mpcps_KA.FR"),],design = ~ Genotype)
    res.p5<-myOutput(DESeq(dds.p5), contrast = c("Genotype","Mpcps_KA","Tak1"),tag = "MpcpsKA_vs_Tak1_FR.pw")
    
    dds.p6<-make_dds(sampleTable[Cond %in% c("Mpcps_KA.FR","Mpcps_Mock.FR"),],design = ~ Genotype)
    res.p6<-myOutput(DESeq(dds.p6), contrast = c("Genotype","Mpcps_KA","Mpcps_Mock"),tag = "Mpcps_KAvsMock_FR.pw")
    
    dds.p7<-make_dds(sampleTable[Cond %in% c("Tak1.FR","Tak1.cW"),],design = ~ Treatment)
    res.p7<-myOutput(DESeq(dds.p7), c("Treatment","FR","cW"),tag = "FR_vs_cW_Tak1.pw")
    
    dds.p8<-make_dds(sampleTable[Cond %in% c("Mpcps_Mock.FR","Mpcps_Mock.cW"),],design = ~ Treatment)
    res.p8<-myOutput(DESeq(dds.p8), c("Treatment","FR","cW"),tag = "FR_vs_cW_MpcpsMock.pw")
    
    dds.p9<-make_dds(sampleTable[Cond %in% c("Mpcps_KA.FR","Mpcps_KA.cW"),],design = ~ Treatment)
    res.p9<-myOutput(DESeq(dds.p9), c("Treatment","FR","cW"),tag = "FR_vs_cW_MpcpsKA.pw")
    
  }
  
  {#Compare contrast-design results with pairwise comparison
    list.pw<-list.files(outDir,pattern = "DEG-.*pw.tsv",full.names = T)
  
    dt.compare<-data.table()
    for (file.pw in list.pw){
      file.con<-sub(".pw","",file.pw)
      tag<-sub(".*DEG-(.*).pw.tsv","\\1",file.pw)
      res.pw<-res2list(file.pw)
      res.con<-res2list(file.con)
      entry<-data.table(
        Tag=tag,
        Up.coverage=length(intersect(res.pw$Data.Up,res.con$Data.Up))/length(res.pw$Data.Up),
        Down.coverage=length(intersect(res.pw$Data.Down,res.con$Data.Down))/length(res.pw$Data.Down)
      )
      dt.compare<-rbind(dt.compare,entry)
    }
    file.compare<-file.path(outDir,"Contrast_validation.tsv")
    write.table(dt.compare,file.compare,row.names = F,col.names = T)
  }
  

  #Make matrix for all comparisons
  res.all<-rbind(res.MpcpsMock.cW,
                 res.MpcpsKA.cW,
                 res.KA.cW,
                 
                 res.MpcpsMock.FR,
                 res.MpcpsKA.FR,
                 res.KA.FR,
                 
                 res.FR.Tak1,
                 res.FR.MpcpsMock,
                 res.FR.MpcpsKA
                 )
  UpSet.res.l<-res.all[,.(Up=(padj<pthre)&(log2FoldChange>FCthre),
                       Down=(padj<pthre)&(log2FoldChange< -FCthre)
                       ),by=.(Gene,Tag)]
  UpSet.cps<-dcast(UpSet.res.l,Gene~Tag, value.var = c("Up","Down"),sep = ".")
  
}
#Venndiagram
# library(eulerr)
{#Count sums
cols.all<-setdiff(names(UpSet.cps),c("Gene","Sum"))
UpSet.cps[is.na(UpSet.cps)] <- F
UpSet.cps.sum<-data.table(cond=cols.all,
                           sum=as.numeric(UpSet.cps[, lapply(.SD, sum),.SDcols=cols.all])
)


# UpSet.cps[is.na(UpSet.cps)] <- 0
}
{ # euler plots
  
  cols.cps.cW<-c("Up.MpcpsMock_vs_Tak1_cW","Down.MpcpsMock_vs_Tak1_cW",
                 "Up.Mpcps_KAvsMock_cW","Down.Mpcps_KAvsMock_cW"
  )
  cols.cps.FR<-c("Up.MpcpsMock_vs_Tak1_FR","Down.MpcpsMock_vs_Tak1_FR",
                 "Up.Mpcps_KAvsMock_FR","Down.Mpcps_KAvsMock_FR"
  )
  cols.FR<-c("Up.FR_vs_cW_Tak1","Down.FR_vs_cW_Tak1",
             "Up.FR_vs_cW_MpcpsMock","Down.FR_vs_cW_MpcpsMock"
  )
  cols.cps.down<-c("Down.MpcpsMock_vs_Tak1_FR",
                   "Up.Mpcps_KAvsMock_FR",
                   "Up.FR_vs_cW_Tak1"#,
                   # "FR_vs_cW_Tak1.Down"
  )
  cols.cps.up<-c("Up.MpcpsMock_vs_Tak1_FR",
                 "Down.Mpcps_KAvsMock_FR",
                 "Down.FR_vs_cW_Tak1"
  )
  list.cps<-list(cps.cW=cols.cps.cW,
                 cps.FR=cols.cps.FR,
                 FR=cols.FR,
                 cps.down=cols.cps.down,
                 cps.up=cols.cps.up
  )
  
  for (l in names(list.cps)){ 
    cols.sel<-unlist(list.cps[l])
    data.venn<-as.data.frame(UpSet.cps[,..cols.sel])
    rownames(data.venn)<-UpSet.cps$Gene
    set.seed(123)
    v <- euler(data.venn,shape="ellipse")
    labels<-UpSet.cps.sum[match(cols.sel, cond),paste0(cond,"\n(Total:",sum,")")]
    # error_plot(v)
    eulerr_options(edges=list(lwd=1),
                   labels=list(fontsize = 4,fontfamily="sans",font=3),
                   quantities=list(fontsize = 8),
                   padding=unit(2,"mm")
                   )
    p<-plot(v, 
            #fills = c("#ffd5d5","#d5e5ff","#ff8080","#87aade"),
            fills = F,
            labels=labels,
            
            quantities = TRUE)
    print(p)
    ggsave(file.path(outDir,paste0("Euler-",l,".pdf")),plot = p,width = 6,height = 6,units = "cm")
    
    #Try to find genes in each sector
    list.sect<-as.data.table(v$original.values,keep.rownames = T)
    names(list.sect)<-c("Sect","Num")
    setorder(list.sect,-Num)
    cts<-UpSet.cps[,c("Gene",cols.sel),with=F]
    
    for (s in list.sect[Num>0,Sect]){
      s.conds<-unlist(strsplit(s,"&"))
      s.otherconds<-setdiff(cols.sel,s.conds)
      cts[,ct1:=rowSums(.SD,na.rm=T),.SDcols=s.conds]
      if (length(s.otherconds)>0){
        cts[,ct2:=rowSums(.SD,na.rm=T),.SDcols=s.otherconds]
      }else{
        cts[,ct2:=0]
      }
      
      
      
      cts[ct1==length(s.conds)&ct2==0,Tag:=s]
      
      if (nrow(cts[ct1==length(s.conds)&ct2==0,])!=list.sect[Sect==s,Num]) warning("Something wrong!")
    }
    
    # list.euler<-
    list.euler<-merge(cts[!is.na(Tag),c("Tag","Gene")],
                      Mpnames[,.(MpTak_v6.1,Name,annots)],
                      by.x="Gene",by.y="MpTak_v6.1",all.x=T
    )
    list.euler[!is.na(Name)&!Name=="",GeneName:=paste0(Gene," [",Name,"]")]
    # list.euler[!is.na(GeneName),GeneName:=Gene]
    list.euler[is.na(list.euler)]<-""
    
    list.euler[,Tag:=factor(Tag,levels = list.sect$Sect)]
    setorder(list.euler,Tag,Gene)
    list.euler[,No.:=1:.N,by=.(Tag)]
    write.table(list.euler[,.(Tag,No.,Gene,Name,GeneName,annots)],
                file=file.path(outDir,paste0("Euler-",l,"_GeneList.tsv")),
                row.names = F,col.names = T,
                sep = "\t",quote = T)
    
    
  }
}




message("DESeq analysis finished!")



