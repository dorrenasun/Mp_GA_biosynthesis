rm(list = ls())
#################### Environment ####################
source("00-Common.R")

pkg_list<-c("data.table","tximeta","DESeq2",
            "multcomp","multcompView","rstatix",
            "ggplot2","scales","ggrepel")
for (p in pkg_list) {check_package(p)}
# if(!require(org.Mpolymorpha.eg.db)) source("00-GO_clustering.R")


if (!require(lasso2)) {
  packageurl <-"https://cran.r-project.org/src/contrib/Archive/lasso2/lasso2_1.2-22.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
}

make_cld<-function(res_test,order,threshold=0.05){
  res.test<-copy(res_test)
  res.test[,`:=`(group1=factor(group1,levels = order),
                 group2=factor(group2,levels = order)
  )]
  setorder(res.test,group1,group2)
  res.sym<-rbind(res.test[,.(group1,group2,p.adj)],res.test[,.(group1=group2,group2=group1,p.adj)])
  res.test.w<-dcast(res.sym,group1~group2, value.var = "p.adj")
  res.test.m<-as.matrix(res.test.w,rownames = "group1")
  Letters.test<-multcompLetters(res.test.m,threshold=threshold)$Letters
  CLD.test<-data.table(Cond=names(Letters.test),`.group`=Letters.test)
  # CLD.test<-data.table(Group=names(Letters.test),`.group`=Letters.test)
  return(CLD.test)
}

#################Parameters#######################

outDir<-file.path(path.out,"06-Barplot-Mpcps")
if (dir.exists(outDir)) unlink(sub("\\/$","",outDir), recursive=TRUE)
if (!dir.exists(outDir)) dir.create(outDir,recursive = T)

sampleTable.cps<-sampleTable[!Genotype=="Mpdella",]
tx2gene.Tak1<-ids.v6[!grepl("MpUg",gene),]

se <- tximeta(sampleTable.cps,type="salmon",skipMeta=T,txOut=F,tx2gene = tx2gene.Tak1)
dds <- DESeqDataSet(se, design =  ~ Genotype+Treatment+Genotype:Treatment)
dds <- estimateSizeFactors(dds)

list.goi<-list(
  ABA=c(
    MpABA4="Mp5g19490",
    MpNCED="Mp2g07800",
    MpCYP707A="Mp2g25940",
    MpHYP="Mp5g18910",
   # MpLEA1="Mp5g23710",
   # MpLEA2="Mp6g13390",
    MpLEAL1="Mp4g09300",
    MpLEAL5="Mp4g05760",
   # MpLEAL7="Mp2g26410",
    MpDHN1="Mp6g07540"
  ),
  GA=c(
    MpCPS="Mp2g07200",
    MpKS="Mp6g05950",
    MpKOL1="Mp3g18320",
    MpKAOL1="Mp4g23680",
    MpKAOL2="Mp1g25410",
    MpBNB="Mp3g23300"
  )
  
)



# tpm<-as.data.table(assay(se,"abundance"),keep.rownames = "Gene")

{
  Lwd=0.3514
  #theme
  theme_pub<-theme_bw()+theme(
    aspect.ratio = 1.3,
    # plot.margin = unit(c(2,1,1,1), "cm"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour="black",linewidth=1*Lwd),
    text = element_text(family = "sans",size = 6, colour = "black"),
    strip.background = element_blank(),
    strip.text=element_text(family = "sans",size = 6, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y =  element_text(family = "sans",size = 6, colour = "black"),
    axis.text  = element_text(family = "sans",size = 6, colour = "black"),
    # axis.text.x = element_text(family = "sans",size = 6, angle=45,hjust = 1,colour = "black"),
    axis.ticks = element_line(colour="black",linewidth=0.5*Lwd),
    # axis.text.x = element_blank(),
    legend.background = element_blank(),
    legend.key.size=unit(0.5,"cm"),
    legend.text.align = 0#,
    # legend.position = "none"
  )

}

counts.dt<-as.data.table(counts(dds,normalized=T, replaced=FALSE),keep.rownames = "Gene")
counts.sel<-melt(counts.dt[Gene %in% unlist(list.goi),],id.vars = "Gene",variable.name = "Sample",value.name = "Counts")

for (l in names(list.goi)){
  goi<-data.table(Name=names(list.goi[[l]]),
                  Gene=unlist(list.goi[[l]])
                  )

  counts.l<-merge(goi,counts.sel,by="Gene",all.x=T)
  
  counts.l[,Gene:=factor(Gene,levels=goi$Gene)]
  setorder(counts.l,Gene)
  counts.l[,Name:=paste0(Name,"(",Gene,")")]
  counts.l[,Name:=factor(Name,levels=unique(Name))]
  counts.l<-merge(counts.l,sampleTable.cps[,.(names,Genotype,Treatment,Cond)],by.x="Sample",by.y="names",all.x=T)
  counts.l.sum<-counts.l[,.(Ave=mean(Counts),SD=sd(Counts)),by=.(Genotype,Treatment,Cond,Name,Gene)]
  
  
  CLD.all<-data.table()
  test.all<-list()
  for (g in goi$Gene){
    data.sub<-counts.l[Gene==g,]
    
    # Levene's test for equal variances (The data must be normally distributed.)
    # Reference: http://www.sthda.com/english/wiki/compare-multiple-sample-variances-in-r
    # the null hypothesis is that all populations variances are equal; (rejected if p<0.05)
    # the alternative hypothesis is that at least two of them differ.
    # 
    message(paste0("Levene's test for Gene: ", g))
    res.Levene<-levene_test(data=data.sub,Counts~Cond)
    print(res.Levene)
    
    #Two-way ANNOVA
    # message(paste0("ANNOVA for Gene: ", g)
    res.aov<-aov(Counts~Cond, data = data.sub)
    # print(summary(res.aov))
    
    #Tukey's HSD test
    res.tukey<-as.data.table(tukey_hsd(res.aov))
    CLD.test<-make_cld(res.tukey,levels(sampleTable.cps$Cond))
    CLD.test$Gene<-g
    CLD.all<-rbind(CLD.all,CLD.test)
    
    
  }

  # CLD.all$.group<-gsub(" ","",CLD.all$.group)
  
  tpm.sum1<-merge(counts.l.sum,CLD.all,all.x=T,by=c("Gene","Cond"))
  
  p<-ggplot(counts.l.sum)+
    geom_bar(aes(x=Genotype,y=Ave,fill=Treatment),
             stat ="identity",
             position=position_dodge(width=0.9),
             colour="black",
             linewidth=Lwd*0.5,
             # alpha=0.6,
             width = 0.8)+
    geom_point(data=counts.l,aes(x=Genotype,y=Counts,fill=Treatment),
               position = position_dodge2(width = 0.9,padding = 0.02),
               size=0.4,
               shape=16,
               color="magenta",
               alpha=0.6
    )+
    geom_errorbar(data=counts.l.sum,
                  aes(x=Genotype,ymin=Ave-SD,ymax=Ave+SD,group=Treatment),
                  position=position_dodge(0.9),
                  width=.2,linewidth=Lwd*0.5)+
    geom_text(data=tpm.sum1,
              aes(x=Genotype,y=Ave+SD*1.5,label=.group,group=Treatment),
              position=position_dodge(0.9),
              size=6*Lwd,
              color="black",
              vjust=-0.5,
              hjust=0.5
    )+
    scale_fill_manual(values=c("#cccccc","#ffd5d5"))+
    facet_wrap(~Name,scales = "free",nrow=1)+
    labs(y="Normalized counts")+
    # scale_y_continuous(limits = c(-0.5,max(dt.sel$Level)+5))+
    theme_pub
  print(p)
  outfile<-file.path(outDir,paste0("Barplot-Mpcps-",l,".pdf"))
  ggsave(plot = p,filename = outfile,width = 25,height = 3.2,units = "cm")
  
}


