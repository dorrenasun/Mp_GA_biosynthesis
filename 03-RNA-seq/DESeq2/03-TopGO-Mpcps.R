rm(list = ls())
#################### Environment ####################
source("00-Common.R")

pkg_list<-c("data.table","topGO","GO.db")
for (p in pkg_list) {check_package(p)}

#################### Other Functions ####################

#Convert DESeq2 results to matrix
Res2List<-function(inFile,Tag="Data",rm=NULL){
  DEGdata<-fread(inFile,header = T,sep = "\t")
  Sig.up<-DEGdata[!is.na(padj) & padj<pthre & log2FoldChange>FCthre & (!Gene %in% rm),Gene]
  Sig.down<-DEGdata[!is.na(padj) & padj<pthre & log2FoldChange<(-FCthre) & (!Gene %in% rm),Gene]
  list.out<-list(Sig.up,Sig.down)
  names(list.out)<-paste0(Tag,c(".Up",".Down"))
  return(list.out)
}
Euler2List<-function(inFile,rm=NULL){
  # inFile<-"Output_a0.01FC0.585/Mpcps/Euler-cps.cW_GeneList.tsv"
  EulerData<-fread(inFile,header = T,sep = "\t")
  EulerData.list<-EulerData[,.(list=.(unique(Gene))),by="Tag"]
  Eulerlist<-EulerData.list$list
  names(Eulerlist)<-EulerData.list$Tag
  return(Eulerlist)
  
}

#################Parameters#######################
inDir<-file.path(path.out,"01-DESeq2-Mpcps")
topN<-1000
Ont<-"BP"

outDir<-file.path(path.out,"03-TopGO")
if(dir.exists(outDir)) unlink(outDir,recursive = T)
if(!dir.exists(outDir))dir.create(outDir)

################# Main #######################
#Convert GO terms to a list
MpGOlist<-MpGO[,.(list=.(unique(GO))),by="MpTak_v6.1"]
MpgeneID2GO<-MpGOlist$list
names(MpgeneID2GO)<-MpGOlist$MpTak_v6.1

# # Merge gene names with ID
# Mpnames[!(Name==""),Syn:=Name]
# Mpnames[Name=="",Syn:=MpTak_v6.1]

#Extract unique gene IDs
geneID.Tak1<-ids.v6[!grepl("MpUg",gene),unique(gene)]

#Retrieve data and collect them in a list
# files.DEG<-list.files(inDir,pattern="^DEG-.*tsv",full.names = T)
# files.sel<-files.DEG[!grepl("sig|pw",files.DEG)]

files.sel<-c(file.path(inDir,"DEG-FR_vs_cW_Tak1.tsv"),
             file.path(inDir,"DEG-MpcpsMock_vs_Tak1_FR.tsv"),
             file.path(inDir,"DEG-Mpcps_KAvsMock_FR.tsv")#,
             )

list.all<-list()
for (f in files.sel){
  tag<-sub(".*DEG-(.*)\\.tsv","\\1",f)
  DEG.temp<-Res2List(f,tag,rm="Mp2g07200")
  list.all<-c(list.all,DEG.temp)
}

#Classic Fisher enrichment analysis
for (l in names(list.all)){
  
  ListGOI<-unlist(list.all[l])

  MpgeneList<-factor(as.integer(geneID.Tak1 %in% ListGOI))
  names(MpgeneList)<-geneID.Tak1
  
  MpGOdata<-new("topGOdata",description="test",ontology=Ont,
                allGenes=MpgeneList,#geneSel = topDiffGenes,
                annot = annFUN.gene2GO,gene2GO=MpgeneID2GO,nodeSize = 5)
  resultFisher <- runTest(MpGOdata, algorithm = "classic", statistic = "fisher")
  
  resFisher<-as.data.table(GenTable(MpGOdata,classicFisher=resultFisher,orderBy="classicFisher",
                                    topNodes=topN,numChar=100))
  resFisher[,classicFisher:=as.numeric(classicFisher)]
  resFisher<-resFisher[classicFisher<pthre,]
  
  #Retrieve significant genes for the enriched GOs
  allGO <- genesInTerm(MpGOdata)[resFisher$GO.ID]
  allGO.sig<-lapply(allGO,function(x)x[x %in% ListGOI])
  sigGO<-data.table(GO.ID=rep(names(allGO.sig), lengths(allGO.sig)),
                   Gene=unlist(allGO.sig))
  sigGO[,Syn:=Mpnames[match(sigGO$Gene,MpTak_v6.1),Name]]
  sigGO[!is.na(Syn)&(!Syn==""),Syn:=paste0(Syn,"(",Gene,")")]
  
  sigGO[is.na(Syn)|Syn=="",Syn:=Gene]
  
  #Add Gene names to resFisher
  sigGO.genes<-sigGO[,.(Genes=paste0(Syn,collapse = ";")),by=.(GO.ID)]
  resFisher[,Genes:=sigGO.genes$Genes[match(resFisher$GO.ID,sigGO.genes$GO.ID)]]
  
  outFile<-file.path(outDir,paste0("Fisher-",Ont,"_",l,".tsv"))
  write.table(resFisher,file=outFile,quote=F,row.names = F,sep="\t")
}

write.table(GO_dbInfo(),file=file.path(outDir,"GO_version_timestamp.txt"),row.names = F,sep = "\t",quote = F)



