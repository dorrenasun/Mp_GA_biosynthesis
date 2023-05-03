rm(list = ls())
#################### Environment ####################
source("00-Common.R")
pkg_list<-c("data.table","GO.db","ggforce","ggrepel",
            "dendextend",
            "igraph","ggraph",
            "GOSemSim","Rtsne","viridisLite")
for (p in pkg_list) {check_package(p)}

if (!dir.exists(path.annot)) stop("Annotation directory does not exist")
# if(!require(org.Mpolymorpha.eg.db,lib.loc = path.annot)) source("00-Prepare_Annotations.R")

#################Parameters#######################
#Directory for reading inputs
inDir<-file.path(path.out,"03-TopGO")
if(!dir.exists(inDir)) stop("Input directory does not exist!")

#Directory for writing outputs
outDir<-file.path(path.out,"04-TopGOFigs")
if(dir.exists(outDir)) unlink(outDir,recursive = T)
if(!dir.exists(outDir)) dir.create(outDir)

#Method of calculation Semantic similarity
methods<-c("Resnik","Lin","Jiang","Rel","Wang")
method<-methods[2] # one of methods

#Targeting ontology
Ont<-"BP"

{#Retrieve offspring and ancester relationships for all BP terms
  term.all<-as.data.table(GOTERM)
  MpGO[,Ontology:=term.all[match(GO,go_id),Ontology]]
  MpGO.ont<-MpGO[Ontology==Ont,.(MpTak_v6.1,GO)]
  anc.all<-switch(Ont,
                  MF = as.data.table(GOMFANCESTOR),
                  BP = as.data.table(GOBPANCESTOR),
                  CC = as.data.table(GOCCANCESTOR))
  names(anc.all)<-c("go","anc")
  # anc.all[,tag:=paste(go,anc,sep="|")]
  
  off.all<-switch(Ont,
                  MF = as.data.table(GOMFOFFSPRING),
                  BP = as.data.table(GOBPOFFSPRING),
                  CC = as.data.table(GOCCOFFSPRING))
  names(off.all)<-c("off","go")
  # off.all[,tag:=paste(off,go,sep="|")]
  
}
# Number of clusters 
Nclust=12

#Limits for visualization
maxOverlap<-100
maxLog10p<-25

#Theme for ggplot2
Lwd=0.3514
fontsize<-6

theme_go<-theme_bw()+theme(
  aspect.ratio = 0.8,
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour="black",linewidth =1*Lwd),
  text = element_text(size = fontsize, colour = "black"),
  strip.background = element_blank(),
  strip.text=element_text(size = fontsize, colour = "black"),
  
  axis.ticks = element_blank(),
  axis.text =  element_blank(),
  axis.title = element_blank(),
  legend.position = "bottom",
  legend.key.size = unit(0.2,"cm"),
  # legend.key.width = unit(0.2,"cm"),
  legend.direction = "horizontal",
  legend.box = "vertical",
  legend.title = element_text(size=fontsize),
  legend.text = element_text(size=fontsize)
)

################# Functions #######################
#Mostly copied from https://rdrr.io/bioc/GOSemSim/src/R/computeIC.R, but 
# removed the dependence on pre-calculated gotbl

computeIC_local <- function(goAnno, ont=Ont) {

  goids<-unique(term.all[Ontology==ont,go_id])
  ## all GO terms appearing in an given ontology ###########
  goterms<-goAnno$GO
  gocount <- table(goterms)
  ## goid of specific organism and selected category.
  goname  <- names(gocount) 
  
  #retrieve child terms
  off.sub<-unique(off.all[off %in% goname,])
  off.sub[,cnt:=gocount[off]]
  off.sum<-off.sub[,.(sum=sum(cnt)),by=go]
  off.sum.all<-data.table(go=goids)
  off.sum.all[go %in% goids,sum:=off.sum[match(goids,go),sum]]
  off.sum.all[is.na(sum),sum:=0]
  
  ## ensure goterms not appearing in the specific annotation have 0 frequency..
  go.diff        <- setdiff(goids, goname)
  m              <- double(length(go.diff))
  names(m)       <- go.diff
  gocount        <- as.vector(gocount)
  names(gocount) <- goname
  gocount        <- c(gocount, m)
  
  cnt <- gocount[goids] + as.vector(off.sum.all[match(goids,go),sum])
  # cnt<-gocount[goids] + sapply(goids, function(i) sum(gocount[off.sub[go==i,off]], na.rm=TRUE))
  names(cnt) <- goids

  
  ## the probabilities of occurrence of GO terms in a specific corpus.
  p <- cnt/sum(gocount)
  ## IC of GO terms was quantified as the negative log likelihood.
  IC <- -log(p)
  return(IC)
}

#Mostly copied from https://rdrr.io/bioc/GOSemSim/src/R/godata.R, but 
# 1) removed the dependence on species database
# 2) retrieve all ancestral terms for each GO annotation

godata_local <- function(ont=Ont, computeIC = TRUE) {
  ont <- toupper(ont)
  ont <- match.arg(ont, c("BP", "CC", "MF"))
  
  message('preparing gene to GO mapping data...')

  goAnno.all<-rbind(MpGO.ont[,.(MpTak_v6.1,GO,anc=GO)],
    merge(MpGO.ont,anc.all[!is.na(anc) & !anc=="all",],by.x="GO",by.y="go",allow.cartesian=T))
  goAnno<-unique(goAnno.all[,.(GID=MpTak_v6.1,GO=anc,ONTOLOGY=ont)])
  # goAnno<-unique(MpGO.ont[,.(GID=MpTak_v6.1,GO=GO,ONTOLOGY=ont)])
  
  res <- new("GOSemSimDATA",
             keys = unique(goAnno$GID),
             ont = ont,
             geneAnno = goAnno#,
             # metadata = metadata(OrgDb)
             )
  
  if (computeIC) {
    message('preparing IC data...')
    IC <- computeIC_local(goAnno, ont)
    res@IC <- IC
  }
  return(res)
}

#Mostly copied from https://rdrr.io/bioc/rrvgo/src/R/rrvgo.R, but removing the dependence on species database
calculateSimMatrix_local <- function(x,
                               # orgdb,
                               semdata=GOSemSim::godata(orgdb, ont=ont),
                               ont=c("BP", "MF", "CC"),
                               method=c("Resnik", "Lin", "Rel", "Jiang", "Wang")) {
  
  # check function args
  ont <- match.arg(ont) 
  method <- match.arg(method) 
  
  # # load orgdb object
  # if(all(is(orgdb) != "OrgDb")) {
  #   orgdb <- loadOrgdb(orgdb)
  # }
  # 
  # filter GO terms not in orgdb
  x <- unique(x)
  found <- x %in% names(semdata@IC)
  hasAncestor <- !is.na(sapply(x, function(x) tryCatch(GOSemSim:::getAncestors(ont)[x], error=function(e) NA)))
  if(all(!found)) {
    warning("No terms were found in orgdb for ", ont,
            "\nCheck that the organism and ontology match the ones provided by orgdb")
    return(NA)
  } else if(!all(found)) {
    warning("Removed ", length(x) - sum(found), " terms that were not found in orgdb for ", ont)
  }
  x <- x[found & hasAncestor]
  
  
  # return the similarity matrix
  m <- matrix(GOSemSim::goSim(x, x, semData=semdata, measure=method),
              ncol=length(x), dimnames=list(x, x))
  
  # removing terms which the similarity couldn't be calculated
  out <- apply(m, 2, function(x) all(is.na(x)))
  return(m[!out, !out])
}

################# Main #######################
{#Load TopGO results to be visualized
  files.topGO<-list.files(inDir,"Fisher-.*",full.names = T)
  resFisher.all<-data.table()
  for (f in files.topGO[grepl("FR",files.topGO)]){
    resFisher<-fread(f)
    resFisher[,Tag:=sub(".*Fisher-BP_(.*)\\.tsv","\\1",f)]
    setorder(resFisher,classicFisher)
    resFisher[,rank:=1:.N]
    resFisher.all<-rbind(resFisher.all,resFisher)
  }
}
list.sets<-list(
  All=c("MpcpsMock_vs_Tak1_FR.Down","Mpcps_KAvsMock_FR.Up",
        "MpcpsMock_vs_Tak1_FR.Up","Mpcps_KAvsMock_FR.Down")
)

resFisher.sel<-resFisher.all[Tag %in% unlist(list.sets)]

#Determine the position and clustering using all enriched GO terms
message(paste0("Detecting clusters ..."))
{ #Calculate semantic similarity
  go.sel<-unique(resFisher.sel$GO.ID)
  anc.sel<-anc.all[go %in% go.sel,]
  off.sel<-off.all[go %in% go.sel,]
  
  node.Mp<-unique(term.all[go_id %in% go.sel,.(go=go_id,term=Term)])
  
  #Calculate semantic similarity among GO terms

  MpGOdata<-godata_local(ont = "BP")
  simMatrix <- calculateSimMatrix_local(node.Mp$go,
                                  # orgdb = "test",
                                  ont = "BP",
                                  semdata = MpGOdata,
                                  method=method
  )

}

{ #Visualiztion with tsne algorithm
  dist<-as.dist(1-simMatrix)
  set.seed(2) #For replicable figures
  #For perplexity and theta, see https://distill.pub/2016/misread-tsne/
  tsne<-Rtsne(dist,is_distance = T,perplexity =7,max_iter=1000,theta=0.1,check_duplicates = FALSE)
  
  rownames(tsne$Y)<-rownames(simMatrix)
  go.pos<-data.table(GO.ID=rownames(tsne$Y),X=tsne$Y[,1],Y=tsne$Y[,2])
  node.Mp[,`:=`(x=go.pos[match(go,GO.ID),X],
                y=go.pos[match(go,GO.ID),Y]
  )]
  
  #Clustering based on tsne distance
  hct<-as.dendrogram(hclust(dist(tsne$Y)))
  hct.cut<-dendextend::cutree(hct,k=Nclust, order_clusters_as_data = F)
  
  go.tsne.cl<-data.table(GO.ID=names(hct.cut),Clust=factor(hct.cut))
  node.Mp[,clust:=go.tsne.cl[match(go,GO.ID),Clust]]
  
  #Save the figure for clustering
  pal<-viridis(n=length(unique(hct.cut)),begin = 0,end=1)
  names(pal)<-as.character(unique(hct.cut))
  hct<-dendextend::color_branches(hct,Nclust,col=pal,
                                  groupLabels = T)
  newlabel<-paste0("[",labels(hct),"] ",term.all[match(labels(hct),go_id),Term])
  labels(hct) <- newlabel
  pdf(file.path(outDir,paste0("GO_clustering.pdf")),width = 7,height = 25)
  par(cex=0.3,mar=c(3,3,3,50))
  plot(hct,horiz=T)
  # axis(1,)
  dev.off()
  
  #A function to decide representative terms
  find_label<-function(GOs){
    #Extract the needed columns
    dt.sub<-data.table(go_col=unique(GOs),clust_col=node.Mp[match(unique(GOs),go),clust])
    names(dt.sub)<-c("go_col","clust_col")
    
    #How many offsprings does a GO term have in its own cluster
    dt.sub[,off_own_cluster:=sapply(go_col,function(x){
      cl<-dt.sub[go_col==x,clust_col]
      off<-off.all[go==x,off]
      own_cluster<-dt.sub[clust_col==cl,go_col]
      return(length(intersect(off,own_cluster)))
    },simplify=T)
    ]
    #Sort the table by cluster and offsprings in its own cluster
    setorder(dt.sub,clust_col,-off_own_cluster)
    
    #How many offsprings does a GO term have all clusters
    dt.sub[,off_all_cluster:=sapply(go_col,function(x){
      off<-off.all[go==x,off]
      return(length(intersect(off,dt.sub$go_col)))
    },simplify=T)
    ]
    #If a GO term's offsprings are mostly in its own cluster, define it as "full-covered"
    dt.sub[,full_cover:=ifelse(off_all_cluster==0|off_own_cluster/off_all_cluster>0.7,1,0)]
    go.fullcover<-dt.sub[full_cover==1,go_col]
    
    #How many full-covered ancestors does a GO term have in its own cluster
    dt.sub[,anc_own_cluster:=sapply(go_col,function(x){
      cl<-dt.sub[go_col==x,clust_col]
      anc<-anc.all[go==x,anc]
      own_cluster<-dt.sub[clust_col==cl & full_cover==1,go_col]
      return(length(intersect(anc,own_cluster)))
    },simplify = T)
    ]
    
    #Show the GO terms which 1) is full-covered AND 2) has no full-covered ancesters in its own cluster
    # dt.sub[,show_tag:=ifelse(full_cover==1&anc_own_cluster==0,1,0)]
    # Alternatively, you can show the most specific terms with no offsprings in the same cluster
    dt.sub[,show_tag:=ifelse(off_own_cluster==0,1,0)]
    
    return(dt.sub)
  }
  
  #Determine which tags to show for all terms
  showTag<-find_label(node.Mp$go)
  node.Mp[,show_tag:=showTag[match(go,go_col),show_tag]]
  setorder(node.Mp,clust)
  outfile.cluster<-file.path(outDir,"All_Clusters.tsv")
  write.table(node.Mp,outfile.cluster,sep="\t",row.names = F,col.names = T)
  

  #Generate a figure for all GO terms
  p.all<-ggplot(node.Mp) +
    geom_point(data=node.Mp[show_tag==0,],
                    aes(x=x,y=y,fill=clust),
                    stroke=0,
                    shape=21,
                    size=2,
                    alpha=0.6,
                    color="grey50")+
    
    geom_point(data=node.Mp[show_tag==1,],
                    aes(x=x,y=y,fill=clust,
                        stroke=0.5*Lwd*show_tag),
                    shape=21,
                    size=2,
                    alpha=0.7,
                    color="black")+
    scale_fill_manual(values = pal,
                      breaks=levels(node.Mp$clust),
                      guide = guide_legend(title="Functional cluster",
                                           direction = "horizontal",
                                           ncol=3,#order=3,
                                           title.position = "top")
                      )+
    geom_text_repel(data=node.Mp[show_tag==1,],aes(x=x,y=y,label = term),
                    color="black",
                    segment.size= 0.2*Lwd,
                    min.segment.length = 0,
                    max.overlaps=50,
                    box.padding = 0.2*Lwd,
                    seed=123,
                    size=2*Lwd
    )+
    theme_go
  print(p.all)
  outfile.plot.all<-file.path(outDir,paste0("GOplot-",Ont,"_everthing.pdf"))
  ggsave(filename = outfile.plot.all,plot = p.all,height = 10,width = 10,units = "cm")
}

{ #Generate figure for each set of enriched terms
  file.label<-file.path(".","All_Clusters_manual.tsv")
  if (file.exists(file.label)){
    man_lab<-fread(file.label)
    
  }
  for (s in names(list.sets)){
    tags<-unlist(list.sets[[s]])
    message(paste0("Processing ",s," ..."))
    #Extract for each tag
    resFisher.sub<-resFisher.all[Tag%in%tags,]
    resFisher.sub[,Tag:=factor(Tag,levels = tags)]
    #Determine if GO terms are connected by their shared genes
    #Build the edges by gene overlaps between GO terms
    sigGO.l1<-resFisher.sub[,.(Gene=unlist(strsplit(Genes,";"))),by=.(Tag,GO.ID)]
    sigGO.w1<-dcast.data.table(sigGO.l1,Gene~GO.ID+Tag,fun.aggregate = function(x){1},value.var= "GO.ID",fill = 0,drop = F)
    mat <- crossprod(as.matrix(sigGO.w1,rownames="Gene"))
    mat.dt<-as.data.table(mat,keep.rownames = "from")
    edges.ori<-melt.data.table(mat.dt,id.vars = "from",variable.name = "to",value.name = "overlap",variable.factor = F)
    edges.ori[,`:=`(go.from=sub("_.*","",from),
                    tag.from=sub(".*?_","",from),
                    go.to=sub("_.*","",to),
                    tag.to=sub(".*?_","",to)
                    )
      
    ]
    #Calculate overlap fraction
    edges.reform<-edges.ori[(!go.from==go.to) & (tag.from==tag.to) & (overlap>0),.(tag=tag.from,from=go.from,to=go.to,overlap)]
    edges.reform[,Nfrom:=sigGO.l1[GO.ID==from & Tag==tag,length(Gene)],by=.(tag,from)]
    edges.reform[,Nto:=sigGO.l1[GO.ID==to & Tag==tag,length(Gene)],by=.(tag,to)]
    
    edges.reform[,`:=`(Per=overlap/max(Nfrom,Nto),
                       check_redun=paste0(tag,sort(c(from,to)),collapse  ="|")),by=.(tag,from,to)]
    edges.filtered<-edges.reform[Per>.3,]
    #Remove symmetric edges
    edges.uniq<-edges.filtered[!duplicated(check_redun),]
    
    #Keep the needed columns & #Merge XY coordinates
    edges<-edges.uniq[,.(Tag=tag,from,to,overlap)]
    edges[,`:=`(x=node.Mp[match(from,go),x],
                y=node.Mp[match(from,go),y],
                xend=node.Mp[match(to,go),x],
                yend=node.Mp[match(to,go),y],
                Tag=factor(Tag,levels = tags)
                )]
    # #Fix the super-large values to limits of visualization
    # edges[overlap>maxOverlap,overlap:=maxOverlap]
    
    
    #Retrieve cluster and XY postions
    resFisher.sub.1<-merge(resFisher.sub,node.Mp[,.(go,clust,x,y)],by.x="GO.ID",by.y="go",all.x=T)
    if (exists("man_lab")){
      resFisher.sub.1[,show_tag:=man_lab[match(GO.ID,go),show_tag]]
      
    }else{
      showTag.sub<-find_label(resFisher.sub.1$GO.ID)
      
      resFisher.sub.1[,show_tag:=showTag.sub[match(GO.ID,go_col),show_tag]]
      
    }
    
    setorder(resFisher.sub.1,Tag,clust)
    # #Fix the super-large values to limits of visualization
    # resFisher.sub.1[-log10(classicFisher)>maxLog10p,classicFisher:=10^(-maxLog10p)]
    
    
    #Make plot
    #Refer to https://mran.microsoft.com/snapshot/2017-08-20/web/packages/ggrepel/vignettes/ggrepel.html for ggrepel
    p.sub<-ggplot(resFisher.sub.1) +
      # geom_segment(data=edges,
      #              aes(x=x,y=y,xend=xend,yend=yend),
      #              alpha=0.4,
      #              linewidth=0.2*Lwd,
      #              color="grey70") +
      geom_point(data=resFisher.sub.1[is.na(show_tag)|show_tag==0,],
                      aes(x=x,y=y,fill=clust,size=-log10(classicFisher)),
                      stroke=0.2*Lwd,
                      shape=21,
                      alpha=0.5,
                      color="black")+
      
      geom_point(data=resFisher.sub.1[!is.na(show_tag) & show_tag==1,],
                      aes(x=x,y=y,fill=clust,size=-log10(classicFisher)),
                      stroke=0,
                      shape=21,
                      alpha=0.5,
                      color="black")+
      geom_point(data=resFisher.sub.1[!is.na(show_tag) & show_tag==1,],
                 aes(x=x,y=y,fill=clust,size=-log10(classicFisher)),
                 stroke=0.5*Lwd,
                 shape=1,
                 alpha=1,
                 color="black")+
      scale_size(range = c(0.5,4),
                 breaks = seq(0,maxLog10p,by=5),
                 guide = guide_legend(title="-log10(p)",
                                      direction = "horizontal",
                                      nrow=1,order=1,
                                      label.position ="right",
                                      title.position = "left"))+
      scale_fill_manual(values = pal,
                        breaks = levels(resFisher.sub.1$clust),
                        guide = guide_legend(title="Functional cluster",
                                             direction = "horizontal",
                                             nrow=1,order=2,
                                             label.position ="top",
                                             title.position = "left"))+
      geom_text_repel(data=resFisher.sub.1[show_tag==1,],aes(x=x,y=y,label = Term,color=clust),
                     
                     segment.size= 0.2*Lwd,
                     min.segment.length = 0,
                     max.overlaps=150,
                     box.padding = 0.1*Lwd,
                     seed=123,
                     size=6*Lwd
      )+
      scale_color_manual(values = pal,
                        breaks = levels(resFisher.sub.1$clust))+
      facet_wrap(~Tag,scales = "fixed",nrow = 2)+
      theme_go
      
    print(p.sub)
    outfile.plot<-file.path(outDir,paste0("GOplot-",Ont,"_",s,".pdf"))
    ggsave(filename = outfile.plot,plot = p.sub,height = 15,width = 11,units = "cm")
  
    p.sub1<-ggplot(resFisher.sub.1) +
      geom_segment(data=edges,
                   aes(x=x,y=y,xend=xend,yend=yend),
                   alpha=0.4,
                   linewidth=0.1*Lwd,
                   color="grey80") +

      geom_point(data=resFisher.sub.1,
                 aes(x=x,y=y,fill=clust,size=-log10(classicFisher),
                     stroke=0.5*Lwd),
                 shape=21,
                 alpha=0.7,
                 color="black")+
      scale_size(range = c(0.5,4),
                 breaks = seq(0,maxLog10p,by=5),
                 guide = guide_legend(title="-log10(p)",
                                      direction = "horizontal",
                                      ncol=1,order=1,
                                      title.position = "top"))+
      scale_fill_manual(values = pal,
                        breaks = levels(resFisher.sub.1$clust),
                        guide = guide_legend(title="Functional cluster",
                                             direction = "horizontal",
                                             ncol=2,order=1,
                                             title.position = "top"))+
      geom_text_repel(data=resFisher.sub.1,aes(x=x,y=y,label = Term),
                      color="black",
                      segment.size= 0.2*Lwd,
                      min.segment.length = 0,
                      max.overlaps=50,
                      box.padding = 0.1*Lwd,
                      seed=123,
                      size=1*Lwd
      )+
      facet_wrap(~Tag,scales = "fixed",nrow = 2)+
      theme_go
    outfile.plot1<-file.path(outDir,paste0("GOplot-",Ont,"_",s,".all.pdf"))
    ggsave(filename = outfile.plot1,plot = p.sub1,height = 15,width = 11,units = "cm")
    
    
# {     
#   setorder(resFisher.sub.1,clust,Tag,rank)
#   resFisher.sub.1[,GO.ID:=factor(GO.ID,levels = unique(GO.ID))]
#   resFisher.sub.1[,Tag:=factor(Tag,levels = rev(levels(Tag)))]
# 
#   p.sub2<-ggplot(resFisher.sub.1) +
#       geom_tile(data=resFisher.sub.1,
#                 aes(x=GO.ID,y=Tag,fill=clust,alpha=-log10(classicFisher)),
#                 height=0.95,
#                 width=0.95,
#                 linewidth=0)+
# 
#       scale_fill_manual(values = pal,
#                         breaks = levels(resFisher.sub.1$clust),
#                         guide = guide_legend(title="Functional cluster",
#                                              direction = "horizontal",
#                                              ncol=2,order=3,
#                                              title.position = "top"))+
#     theme_pub+theme(aspect.ratio = 0.2,panel.border = element_blank())
#     print(p.sub2)}
    
  }
}

