rm(list = ls())
message("Running 06-BlastMatrix.R !")

if (!exists("read_fasta", mode="function")) source("00-Functions.R")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",update = T)

check_package<-function(pkg){
  if (!require(pkg,character.only = T,quietly = T)){
    # message(paste(pkg,"is not here. Install."))
    BiocManager::install(pkg,ask=F)
    library(pkg,character.only = T)
  }
}
pkg_list<-c("data.table","ggplot2","treeio","ggtext","pheatmap")
for (p in pkg_list) {check_package(p)}

# outdir<-"../02-Output/20220607_211847-Archived/"
#load gene ids from phylogentic tree
tree.file<-list.files(path = outdir,pattern = ".*07.*ordered\\.newick",full.names = T)
tree1<-read.newick(tree.file)
tree1$tip.label<-gsub("[ ']","",tree1$tip.label)

#Extract tip names
annot<-data.table(No=1:length(tree1$tip.label),Name=tree1$tip.label)
annot_long<-annot[,.(Name_sep=unlist(strsplit(Name,"[;\\|]"))),by=.(No,Name)]

#Load information for species etc.
file.ids<-list.files(outdir,"02-Candidates.*id$",full.names = T)
ids<-fread(file.ids)
ids_long<-ids[,.(Name_sep=unlist(strsplit(Name,"[;\\|]"))),by=.(Num,Name,Clade,Species,Database)]

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


#Find some nodes:
{ #Extract various nodes
  tree1.tbl<-as_tibble(tree1)
  node.og.seq1<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("Mp7g14540.1",tree1.tbl$label)]
  node.og.seq2<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("AT5G38970.1",tree1.tbl$label)]
  node.outgroup<-MRCA(tree1.tbl,node.og.seq1,node.og.seq2)$node
  
  node.g2.seq1<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("YFGP-2106464",tree1.tbl$label)]
  node.g2.seq2<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("Psubtropica_DN49014_c0_g1_i3.p1",tree1.tbl$label)]
  node.g2<-MRCA(tree1.tbl,node.g2.seq1,node.g2.seq2)$node
  
  node.MpKAOL2<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("Mp1g25410.1",tree1.tbl$label)]
  node.someCYP817<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("YFGP-2107049",tree1.tbl$label)]
  node.CYP817<-MRCA(tree1.tbl,node.MpKAOL2,node.someCYP817)$node
  
  node.g3.seq1<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("AKXB-2011361",tree1.tbl$label)]
  node.g3.seq2<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("Apun_utg000128l.60.1",tree1.tbl$label)]
  node.g3<-MRCA(tree1.tbl,node.g3.seq1,node.g3.seq2)$node
  
  node.OsCYP729<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("LOC_Os10g23160.1",tree1.tbl$label)]
  node.distCYP729<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("Adc01116",tree1.tbl$label)]
  node.CYP729<-MRCA(tree1.tbl,node.OsCYP729,node.distCYP729)$node
  
  
  node.AtKAO1<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("AT1G05160.1",tree1.tbl$label)]
  node.distKAO<-tree1.tbl$node[!is.na(tree1.tbl$label)& grepl("AANG009306",tree1.tbl$label)]
  node.KAO<-MRCA(tree1.tbl,node.AtKAO1,node.distKAO)$node
  
  list.nodes<-c(node.outgroup,node.g2,node.CYP817,node.g3,node.CYP729,node.KAO)
  for (n in 1:length(list.nodes)){
    offspr<-offspring(tree1.tbl,list.nodes[n],tiponly = T)$node
    annot_sub[No %in% offspr, Group:=n]
  }
  
}


# # ids.sele<-ids[Name %in% groups[Clade=="CYP729",ID],]
# ids.sele<-ids[ID_sep %in% groups$ID_sep,]
# ids.m<-merge(ids.sele,groups,by="ID_sep")
# # ids.new<-unique(ids.m[!grep(".f$",ID),.(Num,ID,No,Clade)])
# ids.new<-unique(ids.m[!grep("S$",Num),.(Num,ID,No,Clade)])
# setkeyv(ids.new,c("No"))

# stop()

#load self-blast results
blast.file<-list.files(path = outdir,pattern="08-.*blast\\.tsv$",full.names = T)
if (!file.exists(blast.file)) stop("Self-blast file not found! Stopped!")
blast<-fread(blast.file,header = F)
names(blast)<-unlist(strsplit("qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos",split=" "))


#Extract blast terms
gap=0.5
blast.sele<-blast[(qaccver %in% annot_sub$Num) & (saccver %in% annot_sub$Num) & (evalue<0.001),]
blast.sele[,`:=`(qname=annot_sub[match(qaccver,Num),label],
                 x=annot_sub[match(qaccver,Num),No] + (annot_sub[match(qaccver,Num),Group]-1)*gap,
                 # xg=,
                 sname=annot_sub[match(saccver,Num),label],
                 y=annot_sub[match(saccver,Num),No] + (annot_sub[match(saccver,Num),Group]-1)*gap
                 # yg=annot_sub[match(qaccver,Num),Group]
)]

setorder(annot_sub,"No")
blast.sele[,`:=`(qname=factor(qname,levels = unique(annot_sub$label)),
                 sname=factor(sname,levels = rev(unique(annot_sub$label)))
)]

blast.sele.m<-blast.sele[,.(qname,sname,x,y,pident)]
blast.sele.m[pident<40,pident:=NA]


#load MSA and calculate pid (pid(AB)=# matches/total length of B)
msa.file<-list.files(path = outdir,pattern="06.*LINSI.fasta$",full.names = T)
if (file.exists(msa.file)){
  msa<-read_fasta(msa.file)
  pIdent<-function(x,y){
    if (!nchar(x)==nchar(y)) message("Warning! sequences of different length!")
    temp<-data.table(x=unlist(strsplit(x,"")),y=unlist(strsplit(y,"")))
    temp[,match:=(x==y)]
    pid<-sum(temp[!y=='-',match])/nrow(temp[!y=='-',])
    return(pid)
  }
  msa.pid<-data.table()
  for (i in 1:nrow(msa)){
    for (j in 1:nrow(msa)){
      idI<-msa[i,ID]
      idJ<-msa[j,ID]
      if ((idI %in% annot_sub$Num)&(idJ %in% annot_sub$Num)){
        seqI<-msa[i,Sequence]
        seqJ<-msa[j,Sequence]
        PID12=pIdent(seqI,seqJ)
        
      }
      blast.sele[saccver==idI & qaccver == idJ, msa_pid:=PID12]
      
    }
  }
  msa.sele<-blast.sele[,.(qname,sname,x,y,pident=msa_pid)]
  msa.sele[pident<.4,pident:=NA]
  

}


#Plot the matrix
if (!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
fontsize=2.5
lwd=0.3514
theme_hm<-theme_bw()+theme(aspect.ratio = 1,
                           plot.background = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.text.x.top = element_text(family = "sans",angle=90,hjust=0,size=fontsize),
                           axis.text.y = element_text(family = "sans",size=fontsize),
                           axis.title = element_blank(),
                           axis.ticks.x = element_line(linewidth = lwd*0.5,color="black"),
                           axis.ticks.y = element_line(linewidth = lwd*0.5,color="black"),
                           axis.ticks.length = unit(0.1,"cm"),
                           legend.key.height = unit(0.3,"cm"),
                           legend.position = "bottom"
                
                           )
p1<-ggplot(msa.sele,aes(x=x,y=y,fill=pident))+
  #ggplot(blast.sele.m,aes(x=qaccver,y=saccver,fill=pident))+
  geom_tile()+
  scale_x_continuous(position = "top",
                     expand = expansion(mult=0,add=0.5),
                     # limits = c(-0.1,max(blast.sele$x+1)),
                     # breaks=c(0,50,100),
                     # labels=c("A","B","c"),
                     breaks=sort(unique(blast.sele$x)),
                     # labels=sort(unique(blast.sele$x))
                     labels=unique(annot_sub$label)
  )+
  scale_y_reverse(expand = expansion(mult=0,add=0.5),
                  breaks=sort(unique(blast.sele$x)),
                  labels=unique(annot_sub$label)
  )+
  theme_hm#+
ggsave(paste0(outdir,"08-PerIdent_matrix-msa.pdf"),p1,width = 16,height = 20,units = "cm")

  # scale_fill_gradient2(midpoint = .4, low = "blue", mid = "white",
  #                      high = "red", space = "Lab" )
  #scale_y_discrete(position = "")
{p2<-#ggplot(msa.sele,aes(x=Name_1,y=Name_2,fill=pident))+
  ggplot(blast.sele.m,aes(x=x,y=y,fill=pident))+
  geom_tile()+
  scale_x_continuous(position = "top",
                     expand = expansion(mult=0,add=0.5),
                    # limits = c(-0.1,max(blast.sele$x+1)),
                   # breaks=c(0,50,100),
                   # labels=c("A","B","c"),
                   breaks=sort(unique(blast.sele$x)),
                   # labels=sort(unique(blast.sele$x))
                   labels=unique(annot_sub$label)
                   )+
  scale_y_reverse(expand = expansion(mult=0,add=0.5),
                  breaks=sort(unique(blast.sele$x)),
                  labels=unique(annot_sub$label)
                  )+
  # theme_grey()
  theme_hm#+

print(p2)
ggsave(paste0(outdir,"08-PerIdent_matrix-blast.pdf"),p2,width = 16.5,height = 20,units = "cm")

}


