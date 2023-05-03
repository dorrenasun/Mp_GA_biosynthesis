rm(list=ls())
message("Running 07-Visualization.R !")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",update = T)

check_package<-function(pkg){
  if (!require(pkg,character.only = T,quietly = T)){
    # message(paste(pkg,"is not here. Install."))
    BiocManager::install(pkg,ask=F)
    library(pkg,character.only = T)
  }
}
pkg_list<-c("data.table","ggplot2","ggtree",
            # "phytools",
            "treeio","dplyr","tidytree","ggnewscale",
            # "aplot",
            "ggtext")
for (p in pkg_list) {check_package(p)}

if (!exists("read_fasta", mode="function")|!exists("read_hmm", mode="function")) source("00-Functions.R")
#================================================================================================
outdir<-file.path(outdir,"..","20230421_011847")
OutDir<-file.path(outdir,"Non-parametric-noLBA")
# OutDir<-outdir
#Inputs & outputs
tree.file<-list.files(path = OutDir,pattern = ".*07.*ordered\\.newick",full.names = T)
subname.file<-list.files(path = "../01-Input/",pattern = "Name_sub",full.names = T)

# align.file<-file.path(OutDir,"06-Subtree_EINSI.fasta")
# # dom.file<-file.path(OutDir,"05-Subtree_dom.tsv")
# # tm.file<-file.path(OutDir,"06-Subtree_TMHMM.txt")
# ipr.file<-file.path(outdir,"05-Subtree-IPR.tsv")
#================================================================================================
#Parameters
#Create plot
lwd <- 0.3514
fontsize <- 5
colors_Clade <- c(
  "Charophyte" = "tan4",
  "Moss" = "#8DD35F",
  "Liverwort" = "#2CA02C",
  "Hornwort" = "#88AA00",
  "Lycophyte" = "#AB37C8",
  "Fern" = "#5F5FD3",
  "Gymnosperm" = "#5599FF",
  "Angiosperm" = "#0044AA"
)
colors_CKI1 <- c("CKI1" = "grey20",
                 "Others" = "grey50")
clustal.c1 <- "orange"
clustal.c2 <- "red3"
clustal.c3 <- "royalblue"
clustal.c4 <- "springgreen4"
        
fill.clustal <- c(
  "G" = clustal.c1,
  "P" = clustal.c1,
  "S" = clustal.c1,
  "T" = clustal.c1,
  "H" = clustal.c2,
  "K" = clustal.c2,
  "R" = clustal.c2,
  "F" = clustal.c3,
  "W" = clustal.c3,
  "Y" = clustal.c3,
  "I" = clustal.c4,
  "L" = clustal.c4,
  "M" = clustal.c4,
  "V" = clustal.c4
)

fill.domain<-c(
  # "black",
  "#0072B2",
  # "#56B4E9",
  "#CC79A7",
  
  "orange",
  "#009E73"
)
#================================================================================================
      
      

tree1<-read.newick(tree.file)
tree1$tip.label<-gsub("[ ']","",tree1$tip.label)


#Load information for species etc.
file.ids<-list.files(outdir,"02-Candidates.*id$",full.names = T)
ids<-fread(file.ids)
ids_long<-ids[,.(Name_sep=unlist(strsplit(Name,"[;\\|]"))),by=.(Num,Name,Clade,Species,Database)]

# #Read alternative names & species info from file
# names.file<-"../01-Input/Name_substitute.csv"
# names<-as.data.table(read.csv(names.file))

#Extract tip names
annot<-data.table(No=1:length(tree1$tip.label),Name=tree1$tip.label)
annot_long<-annot[,.(Name_sep=unlist(strsplit(Name,"[;\\|]"))),by=.(No,Name)]

annot_long1<-merge(annot_long,ids_long[,.(Name_sep,Num,Clade,Species,Database)],by = "Name_sep",all.x=T)
if (sum(is.na(annot_long1$Num))>0) stop("Something wrong with ID parsing!")

#Truncate some names:
annot_long1[Database=="OneKP",newName_sep:=sub(".*_([A-Z]+-[0-9]+).*","\\1",Name_sep)]
annot_long1[grepl("AagrBONN",Name_sep),newName_sep:=sub("_Sc2ySwM_","_",Name_sep)]
annot_long1[is.na(newName_sep),newName_sep:=Name_sep]

#Merge back the names
annot_sub<-annot_long1[!newName_sep=="Lc_composite_25679.2_fr0",.(Num=paste0(unique(Num),collapse = "|"),
                           label=paste0(unique(newName_sep),collapse = "|"),
                           Clade=paste0(unique(Clade),collapse = "|"),
                           Species=paste0(unique(Species),collapse = "|"),
                           Database=paste0(unique(Database),collapse = "|")
                           ),by=.(No,Name)]
#Manual name info
subname<-fread(subname.file)
subname<-subname[!Short_name=="",]
#Check if all entries could be found in annot_sub
list.mismatch<-setdiff(subname$Long_Name,annot_sub$Name)
if (length(list.mismatch)>0) warning("Some manual names not found in the tree!")
annot_sub[Name %in% subname$Long_Name, label:=paste(label,subname$Short_name[match(Name,subname$Long_Name)],sep = "_")]

#Create label with species:
annot_sub[,Full_Name:=paste0("plain('[')*italic('",Species,"')*plain('] ",label,"')")]

#Fix order of legend:
annot_sub[,Clade:=factor(Clade,levels = names(colors_Clade))]

#Substitue tree labels
tree1$tip.label<-annot_sub[match(tree1$tip.label,Name),label]

#Convert to tibble
tree1.tbl<-as_tibble(tree1)
#remove "'"
tree1.tbl$label<-gsub("'","",tree1.tbl$label,fixed = T)

# stop()
{ #Extract various nodes
  node.og.seq1<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("Mp6g18420.1",tree1.tbl$label)]
  node.og.seq2<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("AT5G04660.1",tree1.tbl$label)]
  node.outgroup<-MRCA(tree1.tbl,node.og.seq1,node.og.seq2)$node

  node.LwKO1.seq1<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("Mcrispata_DN24541_c0_g1_i2.p1",tree1.tbl$label)]
  node.LwKO1.seq2<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("Mp2g01940.1",tree1.tbl$label)]
  node.LwKO1<-MRCA(tree1.tbl,node.LwKO1.seq1,node.LwKO1.seq2)$node

  node.LwKO2.seq1<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("Tlacunosa_DN40656_c0_g1_i7.p1",tree1.tbl$label)]
  node.LwKO2.seq2<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("Mp3g18320.1",tree1.tbl$label)]
  node.LwKO2<-MRCA(tree1.tbl,node.LwKO2.seq1,node.LwKO2.seq2)$node
  
  nodes<-data.table(node=c(node.outgroup,node.LwKO1,node.LwKO2),
                    label=c("Outgroup","Liverwort\nClade 1","Liverwort\nClade 2")
                    )
  nodes[,label:=factor(label,levels = unique(label))]
}


#Join info
tree2<-full_join(as.treedata(tree1.tbl),annot_sub,by="label")
tree2.tbl<-as_tibble(tree2)
tree2.tbl.sub<-as.data.table(tree2.tbl[!is.na(tree2.tbl$Full_Name),c("label","branch.length","Clade","Full_Name")])
# annot_sub[,IsCKI:=ifelse(label %in% sub.CKI1$label,yes = "CKI1",no="Others")]
xlim=6
p1 <- ggtree(tree2, ladderize = F, color="black",size=0.5*lwd) +
  # ggplot(tree2,ladderize=F)+
  # geom_tree(aes(color=IsCKI,size=IsCKI),size=1*lwd)+
  geom_nodelab(
    mapping = aes(label = label),
    color="grey30",
    size = 3 * lwd,
    nudge_x = -0.02,
    nudge_y = 0.2,
    vjust=0,
    hjust = 1
  ) +
  geom_tippoint(aes(color = Clade),size=2,alpha=0,shape=15)+
  scale_color_manual(values = colors_Clade,guide=guide_legend(override.aes = list(alpha = 1))) +
  ggnewscale::new_scale_color()+
  geom_tiplab(
    mapping = aes(color = Clade, label = Full_Name),
    align = F,
    alpha=1,
    geom = "text",
    parse=T,
    offset = 0.02,
    size = fontsize * lwd,
    linetype = "dotted",
    linesize = 0.5 * lwd
  ) +
  scale_color_manual(values = colors_Clade,guide="none") +
  #Scale bar
  geom_treescale(x=0,
                 y=-nrow(annot_sub),
                 linesize=lwd,
                 fontsize = fontsize*lwd,
                 width = 0.4,
                 color='black'
                   )+
  # #Highlights
  geom_hilight(data = nodes,
               mapping=aes(node=node,fill=label),
               type="gradient",
               gradient.direction = 'tr',
               align="right",
               alpha=.6,
               extendto=xlim-0.05,
               to.bottom = T) +
  scale_fill_manual(values = c("grey70","#e9f6c5","#caf1ed"),guide="none")+
  geom_cladelab(data=nodes,
                mapping=aes(node=node,label=label),
                offset=3.2,
                align=T,
                offset.text = -0.05,
                # hjust=0.5,
                hjust=1,
                vjust=0.5,
                # angle=270,
                fontsize=fontsize*2*lwd,
                alpha=1,
                color="black"#,
                # color=c(NULL,"black")#,
                # barsize=0
                )+
  scale_x_continuous(limits = c(0,xlim)) +
  scale_y_reverse(
    expand = expansion(0, 1)
    # name = "IDs",
    # breaks = unique(msa$No),
    # labels = unique(msa$label),
    # position = "right"
  ) +
  
  # theme_grey() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.title = element_blank(),
    axis.text.y=element_blank(),
    # axis.text.y = element_text(
    #   family = "sans",
    #   size = fontsize,
    #   color = "black",
    #   face = "plain"
    # ),
    legend.text = element_text(size=fontsize,color = "black"),
    legend.title = element_text(size=fontsize,color = "black"),
    legend.key.size = unit(0.2,"cm"),
    legend.position = "left"
  )

plot(p1)
outfile1<-file.path(OutDir,sub("newick","pdf",basename(tree.file)))
ggsave(outfile1,p1,width = 20, height = 24,units = "cm")


