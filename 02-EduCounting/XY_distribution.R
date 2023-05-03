check_package<-function(pkg){
  if (!require(pkg,character.only = T,quietly = T)){
    # message(paste(pkg,"is not here. Install."))
    BiocManager::install(pkg,ask=F)
    library(pkg,character.only = T)
  }
}
pkg_list<-c("data.table","ggplot2",
            "multcompView","rstatix",
            "showtext","svglite")
for (p in pkg_list) {check_package(p)}

shapiro.test1<-function(x){
  if (sd(x)>0 & length(x)>=3)
    return(ifelse(shapiro.test(x)$p.value<0.05,"Not Normal","Normal"))
  else return("NA")
}

make_cld<-function(res_test,threshold=0.05){
  res.test<-copy(res_test)
  setorder(res.test,group1,group2)
  res.sym<-rbind(res.test[,.(group1,group2,p.adj)],res.test[,.(group1=group2,group2=group1,p.adj)])
  res.test.w<-dcast(res.sym,group1~group2, value.var = "p.adj")
  res.test.m<-as.matrix(res.test.w,rownames = "group1")
  Letters.test<-multcompLetters(res.test.m, threshold=threshold)$Letters
  CLD.test<-data.table(Group=names(Letters.test),`.group`=Letters.test)
  return(CLD.test)
}
############# Plot settings###########################
#Color-blind friendly pallette
pal <- c("grey40","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d","white")
pie(rep(1,15), col=pal)
colors.fill<-pal[c(1,16,6)]
colors<-pal[c(1,9,6)]
Lwd=0.3514
#theme
theme_barplot<-theme_bw()+theme(
  aspect.ratio = 0.8,
  #plot.margin = unit(c(2,1,1,1), "cm"),
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour="black",linewidth=1*Lwd),
  panel.spacing = unit(1, "mm"),
  text = element_text(family = "sans",size = 6, colour = "black"),
  strip.background = element_blank(),
  strip.text=element_text(family = "sans",size = 6, colour = "black"),
  axis.title = element_text(family = "sans",size = 6, colour = "black"),
  axis.title.x = element_blank(),
  axis.text = element_text(family = "sans",size = 6, colour = "black"),
  axis.ticks = element_line(colour="black",linewidth = 0.5*Lwd),
  legend.text.align = 0,
  legend.position = "none"
)
#############Export settings###########################
Fig.wd<-15
Fig.ht<-4


HomeDir="../Original_stacks-Results-preserved/Coordinates"
OutDir="./"

list=list.files(pattern = ".csv",path=HomeDir,full.names = T)
# list=list[!grepl("Statistics",list)]


# Read all data points. Values in micron units.
data=data.table()
data.sum=data.table()
data.hist=data.table()
#histogram parameters
bw=50
hlim=800

for (f in list){
  temp=fread(f)
  #Parse for genotype
  #Correct y values from imageJ
  temp[,`:=`(Genotype=sub("-DAPI.*","",Label),
    YM=min(YM)+max(YM)-YM)]
  
  #Correct coordinates by linear regression
  xc=mean(temp$XM)
  yc=mean(temp$YM)
  fit=lm(temp$XM~temp$YM)
  k1=coef(fit)[["temp$YM"]]
  b1=coef(fit)[["(Intercept)"]]
  #The other axis
  k2=-1/k1
  b2=xc-k2*yc
  
  #Reference point for calculating the angle
  yr=yc+100
  xr=k1*yr+b1
  
  #New coordinates
  temp[,`:=`(X1=-(k1*YM-XM+b1)/sqrt(1+k1^2),
             Y1=(k2*YM-XM+b2)/sqrt(1+k2^2),
             X2=XM-xc,
             Y2=YM-yc
  )]
  
  #Export fitting parameters and total sum
  sum.temp=data.table(
    Label=unique(temp$Label),
    Genotype=unique(temp$Genotype),
    Count=nrow(temp),
    xc=xc,
    yc=yc,
    k1=k1,
    b1=b1,
    k2=k2,
    b2=b2,
    xr=xr,
    yr=yr,
    angle=atan2(yr-yc,xr-xc)/pi*180
    
  )
  #Histogram of X distribution
  htemp=hist(temp$X2,breaks = seq(-hlim,hlim,by = bw),plot = F)
  htemp.dt=data.table(Label=unique(temp$Label),
                    Genotype=unique(temp$Genotype),
                    Mids=htemp$mids,
                    Ct=htemp$counts,
                    Dens=htemp$density)
  data=rbind(data,temp[,.(Label,Genotype,XM,YM,X1,Y1,X2,Y2)])
  data.sum=rbind(data.sum,sum.temp)
  data.hist=rbind(data.hist,htemp.dt)
}
#Reorder by genotype
order.Genotype=c("Tak1",
                 "Mpcps-LD_L4","Mpcps-LD_L27",
                 "MpCPS-comp_L6","MpCPS-comp_L7"
                 )
data[, Genotype:=factor(Genotype,levels = order.Genotype)]
data.sum[, Genotype:=factor(Genotype,levels = order.Genotype)]
data.hist[, Genotype:=factor(Genotype,levels = order.Genotype)]
order.Type=unique(sub("_.*","",order.Genotype))
data.sum[,Type:=sub("_.*","",Genotype)]
data.sum[,Type:=factor(Type,levels=order.Type)]
data.hist[,Type:=sub("_.*","",Genotype)]
data.hist[,Type:=factor(Type,levels=order.Type)]

#Figure for Total counts
{
  #Shapiro test for normality
  message("Shapiro test for normal distribution")
  res.shapiro<-data.sum[,.(isNorm=shapiro.test1(Count)),by=.(Genotype)]
  print(res.shapiro)
  
  # Levene's test for equal variances (The data must be normally distributed.)
  # Reference: http://www.sthda.com/english/wiki/compare-multiple-sample-variances-in-r
  # the null hypothesis is that all populations variances are equal; (rejected if p<0.05)
  # the alternative hypothesis is that at least two of them differ.
  # 
  message("Levene's test for variance homogeneity")
  res.Levene<-levene_test(data=data.sum,Count~Genotype)
  print(res.Levene)
  

  #One-way annova
  res.aov <- aov(Count ~ Genotype, data = data.sum)
  # Summary of the analysis
  if (summary(res.aov)[[1]]$`Pr(>F)`[1]<0.05) message("One-way anova passed!")
  
  
  #post-hoc analysis
  res.test<-as.data.table(tukey_hsd(res.aov))
  
res.test[,`:=`(group1=factor(group1,levels = order.Genotype),
               group2=factor(group2,levels = order.Genotype)
)]

CLD<-merge(make_cld(res.test),data.sum[,.(Group=Genotype,max=max(Count)),by=.(Genotype)],by="Group")
# setnames(CLD,"Group","Genotype")

#Make plot for counts
  p1<-ggplot(data=data.sum,aes(x=Genotype,y=Count,color=Type))+
  stat_summary(color="black",
               width=0.5,
               fatten=1.5,
               size=1*Lwd,
               fun = mean,
               geom = "crossbar")+
    geom_point(position = position_jitter(seed = 1,width = 0.1),
               alpha=0.6,
               stroke=0,
               size=1
    )+
    geom_text(data=CLD,aes(x=Genotype,y=max,label=.group),
            size=6*Lwd,
            color="black",vjust=-1,hjust=0.5
            )+
  scale_y_continuous(limits = c(0,2000))+
  scale_color_manual(values = colors)+
  theme_barplot+theme(aspect.ratio = 0.8)
  ggsave(file=paste0(OutDir,"/TotalCounts.pdf"), plot=p1, width=Fig.wd, height=Fig.ht,units = "cm")
  
  print(p1)
}
# Figure for histogram
{
  # hist.sum<-data.hist[,.(Ave=mean(Ct),SD=sd(Ct)),by=.(Genotype,Mids)]
  hist.sum<-data.hist[,.(Ave=mean(Dens*bw),SD=sd(Dens*bw)),by=.(Genotype,Type,Mids)]
  hist.sum<-hist.sum[!Ave==0,]
  p2<-ggplot(hist.sum) +
    geom_bar(aes(x=Mids,y=Ave,fill=Type),
             # color="black",
             alpha=0.6,size=0.5*Lwd,width=45, stat="identity")+
    geom_errorbar(aes(x=Mids,y=Ave,ymin = Ave - SD,ymax = Ave + SD),
                  size=0.5*Lwd,width=10,position ="identity" ) + 
    facet_wrap(~Genotype,ncol=1,strip.position = "right")+
    ylab("Relative Frequency")+
    xlab("Distance to the middle (micron)")+
    scale_y_continuous(limits = c(-0.05,0.35),breaks = c(0,0.3))+
    scale_fill_manual(values = colors)+
    theme_barplot+theme(aspect.ratio = .1)
  print(p2)
  ggsave(file=paste0(OutDir,"/X_Frequency.pdf"), plot=p2, width=Fig.wd, height=Fig.ht,units = "cm")
  
}
#library(gridExtra)
#grid.arrange(p0, p1, nrow = 1)

