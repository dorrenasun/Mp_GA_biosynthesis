# rm(list = ls())
message("Now preparing some annotation files...")
#################### Functions ################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",update = T)

check_package<-function(pkg_list){
  pkg_new<-setdiff(pkg_list, installed.packages())
   
  if (length(pkg_new)>0){
    message(paste(pkg_new,"is not here. Install.\n"))
    BiocManager::install(pkg_new,update=T,ask=F)
  }
  for (pkg in pkg_list) library(pkg,character.only = T)
}

read_fasta<-function(file){
  
  if (!file.exists(file)) stop("Fasta file not found!")
  Fasta<-fread(file,header = F,sep="\t",col.names = "Full")
  ID.positions<-grep(">",Fasta$Full,fixed = T)
  ID.lengths<-c(ID.positions[-1],nrow(Fasta)+1)-ID.positions
  
  Fasta[,`:=`(Num=rep(1:length(ID.positions),ID.lengths),
              ID=rep(sub(" .*","",Full[ID.positions]),ID.lengths))]
  Fasta[,`:=`(ID=sub(">","",ID,fixed = T),
              Sequence=Full)]
  Fasta[ID.positions,Sequence:=""]
  Fasta.trans<-Fasta[,.(Text=paste(Full,collapse = "\n"),
                        Sequence=paste(Sequence,collapse = "")
  ),by=.(Num,ID)]
  return(Fasta.trans)
}
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

#################### Environment ################################

pkg_list<-c("GO.db","R.utils","data.table")#,"AnnotationDbi","AnnotationForge")
check_package(pkg_list)

#################### Main ################################
main<-function(){
  path_annot<-"./RefGenome"
  if (!dir.exists(path_annot)) dir.create(path_annot)
  
  #Load ver 6.1 and ver 3.1 gene ids
  message("Loading ids from v6.1 and ver3.1...")
  
  file.gff.v6<-file.path(path_annot,"MpTak_v6.1r1.gff.gz")
  if (!file.exists(file.gff.v6)) {
    download.file(url = "https://marchantia.info/download/MpTak_v6.1/MpTak_v6.1r1.gff.gz",
                  destfile = file.gff.v6)}
  mRNA.v6<-read_gff(file.gff.v6,select = "mRNA")
  ids.v6<-mRNA.v6[,.(transcript=ID,gene=Parent)]
  
  
  file.gff.v3<-file.path(path_annot,"Mpolymorphav3.1.allTrs.gff3.gz")
  if (!file.exists(file.gff.v3)) {
    download.file(url = "https://marchantia.info/download/download/Mpolymorphav3.1.allTrs.gff3.gz",
                  destfile = file.gff.v3)}
  gene.v3<-read_gff(file.gff.v3,select = "gene")
  
  
  
  #Gene ID correspondance between versions
  message("Sorting gene ids of different version...")
  
  file.corr<-file.path(path_annot,"gene_correspondence_MpTak_v6.1.tsv")
  if (!file.exists(file.corr)) download.file(url ="https://marchantia.info/download/MpTak_v6.1/gene_correspondence_MpTak_v6.1.tsv",destfile = file.corr)
  MpGeneSyn<-fread(file.corr)
  #Fix "Mapoly0018s0018,Mapoly0018s0019" issue
  MpGeneSyn.fixed<-MpGeneSyn[, .(JGI3.1 = unlist(strsplit(JGI3.1, ","))), by = .(gene_id,MpTak1_v5.1,MpTak1_v5.1r2,MpTak_v6.1)]
  #Transform to long format for easy query
  MpGeneSyn.long<-melt(MpGeneSyn.fixed,id.vars = c("MpTak_v6.1"),variable.name = "version",value.name = "query")
  #Add missing genes from ver 3
  MpGeneSyn.uniq<-unique(rbind(MpGeneSyn.long[!query=="-",.(query,MpTak_v6.1)],
                               data.table(query=setdiff(gene.v3$Name,MpGeneSyn.long$query),
                                          MpTak_v6.1="-")))
  #Check if any query have redundant v6.1 hit with "-"
  MpGeneSyn.check<-MpGeneSyn.uniq[,.(nHit=.N,
                                     Hits=paste0(sort(MpTak_v6.1),collapse = ";")
  ),by=.(query)]
  MpGeneSyn.check[nHit>1,Hits:=gsub("-;","",Hits,fixed = T)]
  MpGeneSyn.clean<-MpGeneSyn.check[,.(MpTak_v6.1=unlist(strsplit(Hits, ";"))),by=.(query)]
  write.table(MpGeneSyn.clean,file=file.path(path_annot,"MpGeneSyn_clean.tsv"),sep = "\t",row.names = F)
  
  #Check if there are any missing ver 6.1 genes
  list.miss<-setdiff(ids.v6$gene,MpGeneSyn.clean$MpTak_v6.1)
  if (length(list.miss)>0) message("Some ver6 genes not covered!")
  
  list.miss<-setdiff(gene.v3$Name,MpGeneSyn.clean$query)
  if (length(list.miss)>0) message("Some ver3 genes not covered!")
  
  
  
  #Process GO annotations
  message("Preparing GO annotation...")
  
  file.GO<-file.path(path_annot,"Blast2GO_clean_RuiSun202011.tsv")
  if (!file.exists(file.GO)) download.file(url="https://raw.githubusercontent.com/dorrenasun/MpDELLA/0d44e587b230d4df11847745eb25da0ac41c2f34/03-RNA-seq_Analysis/Inputs/Blast2GO_clean_RuiSun202011.tsv",
                                           destfile = file.GO)
  MpGO.ori<-fread(file.GO)
  MpGO.ori[,gene:=sub("\\..*","",query)]
  #Just in case there is any ver 6 bug
  MpGO.ori[gene=="Mp4g0958", gene:="Mp4g09580"]
  
  #Check if there are any genes not covered
  list.miss<-setdiff(MpGO.ori$gene,MpGeneSyn.clean$query)
  if (length(list.miss)>0) message("Some Blast2GO genes not found in ver6!")
  
  #Convert to ver6
  MpGO.fixed<-merge(MpGO.ori,MpGeneSyn.clean,all.x=T,by.x="gene", by.y="query")
  
  
  # Fix inconsistencies between the go.obo used during annotation and GO.db
  # A list of all annotations that is not in current GO.db
  goids.NA<-data.table(go=unique(MpGO.ori[!GO %in% keys(GO.db),GO]))
  
  #Obsolete terms should be removed.
  goids.NA[go %in% keys(GOOBSOLETE), `:=`(Tag="obs",Ont=GOOBSOLETE[[go]]@Ontology),by=go]
  
  #Secondary IDs should be replaced
  ListSec<-data.table(GOid=keys(GO.db),Sec=Secondary(keys(GO.db)))
  goids.NA[!go %in% keys(GOOBSOLETE) & go %in% unlist(ListSec$Sec), `:=`(Tag="sec")]
  if (sum(goids.NA$Tag=="sec",na.rm = T)>0) goids.NA[Tag=="sec",Ori:=ListSec[grep(pattern=go,Sec),GOid],by=go]
  MpGO.fixed[GO %in% goids.NA[Tag=="sec",go],GO:=goids.NA[go==GO,Ori],by=GO]
  #Any other problem?
  if (goids.NA[is.na(Tag),length(go)]>0) {
    warning("Some GO IDs still not fixed!")
    print(goids.NA[is.na(Tag),go])
  }
  
  MpGO.sel<-unique(MpGO.fixed[!MpTak_v6.1=="-",.(MpTak_v6.1,GO)])
  write.table(MpGO.sel,file=file.path(path_annot,"MpGenetoGO_clean.tsv"),sep = "\t",row.names = F)
  
  #Integrate annotations
  
  {#Retrieve nomenclatures
    message("Retrieving gene names...")
    
    file.Mpnames<-file.path(path_annot,"./nomenlatures.txt")
    if (!file.exists(file.Mpnames)) download.file(url ="https://marchantia.info/nomenclature/nomenlatures.txt",destfile = file.Mpnames)
    Mpnames<-fread(file.Mpnames,check.names = T)
    
    #ver5 & 6 names
    Mpnames[grepl("Mp.g[0-9]+",GeneID.Location),gene:=sub(".*(Mp.g[0-9]+).*","\\1",GeneID.Location)]
    Mpnames[gene=="Mp4g0958", gene:="Mp4g09580"]
    
    #ver3 names
    Mpnames[is.na(gene) & grepl("Mapoly",GeneID.Location), gene:=sub(".*(Mapoly.*s[0-9]+).*","\\1",GeneID.Location)]
    
    #convert to ver 6
    Mpnames.valid<-merge(Mpnames[!is.na(gene),],MpGeneSyn.clean,all.x=T,by.x="gene",by.y="query")
    #Combine synonyms
    Mpnames.valid[!synonym=="",Name:=paste0(gene_symbol,"|",synonym)]
    Mpnames.valid[synonym=="",Name:=gene_symbol]
    #Remove space from Names
    Mpnames.valid[,Name:=gsub(" ","",Name)]
    
    #check if there are any missing terms
    list.miss<-Mpnames.valid[is.na(MpTak_v6.1),gene]
    if (length(list.miss)>0) message("Some named genes not found in ver6!")
    

    #Sort out duplicates
    Mpnames.uniq<-Mpnames.valid[!MpTak_v6.1=="-",.(Name=paste0(unique(unlist(strsplit(Name,"[\\|;]"))),collapse = "|"),
                                                   product=paste0(unique(product),collapse = ";")),by=.(MpTak_v6.1)]
    Mpnames.uniq[,MPGENES:=paste0("[",Name,"]: ",product)]
  }
  
  { #Integrate annotations
    message("Integrating gene annotations...")
    file.Mpannot<-file.path(path_annot,"MpTak_v6.1_func_annotation.tsv")
    if (!file.exists(file.Mpannot)) download.file(url ="https://marchantia.info/download/MpTak_v6.1/MpTak_v6.1_func_annotation.tsv",destfile = file.Mpannot)
    Mpannot<-fread(file.Mpannot,sep = "\t",quote = "")
    
    annot.sel<-c("KEGG","KOG","Pfam","PANTHER")
    Mpannot.sel<-Mpannot[db_type %in% annot.sel,]
    
    #convert transcript id to gene id
    Mpannot.sel[,gene:=ids.v6[match(gene_id,transcript),gene]]
    Mpannot.sel[!description=="-",annot:=paste0("[",ref_id,"]: ",description)]
    Mpannot.sel[description=="-",annot:=paste0("[",ref_id,"]")]
    
    Mpannot.sel[,`:=`(dup=duplicated(ref_id)),by=.(gene)]
    setorder(Mpannot.sel,gene,db_type)
    Mpannot.uniq<-Mpannot.sel[!dup==T & db_type %in% annot.sel,.(annots=paste0(annot,collapse = ";\n")),by=.(gene)]
    
    list.problem<-setdiff(Mpannot.uniq$gene,ids.v6$gene)
    if (length(list.problem)>0) {
      message("Problem in IDs of Gene annotations!")
      message(paste0(list.problem, " not a ver6.1 gene!\n"))
    }
    MpNameAnnot<-merge(Mpnames.uniq,Mpannot.uniq,all=T,by.x="MpTak_v6.1",by.y="gene")
    MpNameAnnot[is.na(MpNameAnnot)]<-""
    write.table(MpNameAnnot,file=file.path(path_annot,"MpAnnot_clean.tsv"),sep = "\t",row.names = F)
    
    
  }
  
  message("Annotation preparation finished!")
}

main()