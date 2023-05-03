rm(list = ls())
#Load required packages
if (!require(data.table)) install.packages("data.table")
library(data.table)
if (!require(rvest)) install.packages("rvest")
library(rvest)

if (!require(tools)) install.packages("tools")
library(tools)

#Folder for downloaded files
outdir<-"./Download"
if (!dir.exists(outdir)) dir.create(outdir)
file.error<-file.path(outdir,"Download_failed.txt")
if (file.exists(file.error)) file.remove(file.error)

#Load all sample information for onekp
file.list.all<-file.path(outdir,"Sample-List-with-Taxonomy.tsv.csv")
if (!file.exists(file.list.all)) {
  download.file("https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100627/Sample-List-with-Taxonomy.tsv.csv",
                file.list.all)
}
list.all<-fread(file.list.all)
setnames(list.all,"1kP_Sample","Code")

#Download file with md5 info
file.md5<-file.path(outdir,"100627.md5")
if (!file.exists(file.md5)){
  download.file("https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100627/100627.md5",
                file.md5)
}
md5<-fread(file.md5,header = F)
names(md5)<-c("md5","file")
md5[,basename:=basename(file)]

#Load selected clades
file.prefix<-"Prefixes.csv"
if (!file.exists(file.prefix)) stop(paste(file.prefix," not found! Stop!"))
prefixes<-read.csv(file.prefix,comment.char = "#")

#Select clades based on Prefix
list.sele<-list.all[Clade %in% prefixes$Clade,Code]

homepage<-"https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100627/assemblies/"
home<-read_html(homepage)
species<-home %>% html_nodes("a") %>% html_text2()
species<-setdiff(species,c(""))
list.ftp<-data.table(ID=species,
                     Code=substr(species,1,4))
list.ftp[,Url_pep:=paste0(homepage,species,Code,"-translated-protein.fa.gz")]
list.ftp[,Url_nuc:=paste0(homepage,species,Code,"-translated-nucleotides.fa.gz")]
list.ftp[,Url_assembly:=paste0(homepage,species,Code,"-SOAPdenovo-Trans-assembly.fa.gz")]

#Compare local list with scraping
list.target<-list.ftp[Code %in% list.sele,]
if (sum(!list.sele %in% list.target$Code)){
  message("Some species not available! Check consistency!")
  print(list.all[!Code %in% list.target$Code,])
} else (message("All species are in the server!"))

#Define subfolder by clade
list.target<-merge(list.target,list.all[,.(Code,Clade)],by="Code",all.x=T)
list.target[,Dest:=file.path(outdir,Clade,"")]

#Check if all md5 values were provided by server
list.allfile<-melt(list.target[,.(Dest,Url_pep,Url_nuc,Url_assembly)],id.vars = "Dest",variable.name = "Cat",value.name = "Url")
list.allfile[,basename:=basename(Url)]
list.allfile[,Filename:=file.path(Dest,basename),by=.(basename)]

list.no_md5<-list.allfile[!basename %in% md5$basename,]

if (nrow(list.no_md5)>0){
  setorder(list.no_md5,"basename")
  file.testgz<-"./Testgz.sh"
  write("#!/bin/bash",file.testgz)
  for (ff in 1:nrow(list.no_md5)){
    command<-paste0("gunzip -v -t '",list.no_md5$Filename[ff],"'")
    write(command,file.testgz,append = T)
  }
  
}

#Function to download files
my_download<-function(url,dir,force=F){
  # url<-list.target$Url_pep[62]
  file<-file.path(dir,basename(url))
  md5_server<-unlist(md5[basename==basename(url),md5])
  
  #Decide whether to download the file
  to_download<-F
  if (!file.exists(file)||force==T) {
    to_download<-T
  }else{
      #check md5 of the existing file
    md5_current<-md5sum(file)
    if (length(md5_server)==1 & !(md5_current %in% md5_server)) to_download<-T
    }
  
  if (to_download==T){
    message(paste0("Downloading ",basename(url)))
    tryCatch(download.file(url,file,mode = "wb", quiet = FALSE),
             error = function(e){
               message(paste(file, 'download failed!'))
               write(url,file=file.error,append = T)
               file.remove(file)
               } 
    )
    if (!file.exists(file)){
      message(paste(file, 'Failed to download'))
      write(url,file=file.error,append = T)
    }else if (file.exists(file)) {
      md5_check<-md5sum(file)
      if (length(md5_server)==1){
        if (!(md5_check %in% md5_server)){
          message(paste(file, 'has wrong md5! Removed!'))
          write(url,file=file.error,append = T)
          file.remove(file)
          
        }
      }
    }
  }else{
    message(paste0(basename(url)," already exists!Skipping!"))
  }
}
oldw <- getOption("warn")
options(warn = -1)
old_timeout<-getOption("timeout")
options(timeout = 300)

#Download to subfolders
for (i in 1:nrow(list.target)){
  dir<-list.target$Dest[i]
  if (!dir.exists(dir)) dir.create(dir)
  my_download(list.target$Url_pep[i],dir,force = F)
  my_download(list.target$Url_nuc[i],dir,force = F)
  my_download(list.target$Url_assembly[i],dir,force = F)
}
options(warn = oldw)
options(timeout = old_timeout)
message("Finished downloading OneKP!")