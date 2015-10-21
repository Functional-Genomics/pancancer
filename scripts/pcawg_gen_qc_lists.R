#!/usr/bin/env Rscript
# =========================================================
# Copyright 2015,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this code.  If not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================

# main QC table
qc.syn.id <- 'syn4991537'
# folder 
qc.folder.syn.id <- 'syn3107128'

#metadata.file <- "http://www.ebi.ac.uk/~nf/pcawg/metadata_freeze3_v4.tsv"
#
metadata.syn.id <- "syn5002504"


args <- commandArgs(trailingOnly=TRUE)

if ( length(args) != 3 ) {
  cat("ERROR: usage: pcawg_gen_qc_lists.R synapse_login synapse_password upload@{yes,no}\n")
  cat("This script will download the main QC table from synapse and generate the files with the white and black list of libraries, samples and donors.\n")
  q(status=1)
}

login <- args[1]
pass <- args[2]
upload <- args[3]
  


###########
if (!require(synapseClient)) {
  source('http://depot.sagebase.org/CRAN.R')
  pkgInstall("synapseClient")
  pkgInstall("Rsftp")
  pkgInstall("RCurl")
  pkgInstall("rjson")
  pkgInstall("digest")
  pkgInstall("RUnit")
  require(synapseClient)
}
if (! require("devtools")) {
  install.packages("devtools")
  require(devtools)
  install_github("brian-bot/rGithubClient")
}
require(rGithubClient)
require(RCurl)


cat("Connecting to Synapse...\n")
synapseLogin(login,pass)
cat("Connecting to Synapse...done.\n")


cat("Downloading data from Synapse...\n") 
qc.link <- NULL

try(qc.link<-synGet(qc.syn.id,downloadLocation="./",ifcollision = "overwrite.local"))
if( is.null(qc.link) ) {
  cat("ERROR: Unable to download QC table from Synapse\n")
  cat("For information on about to setup the credentials check https://annaisystems.zendesk.com/hc/en-us/article_attachments/201188358/Pan-Cancer_Researcher_Guide.pdf")
  q(status=1)
}
cat("Downloaded:",getFileLocation(qc.link), "\n")
qc.file <- basename(getFileLocation(qc.link))
cat("Downloading data from Synapse...done.\n")
  
###########
# Metadata
cat("Downloading and reading metadata...\n")
try(md.link<-synGet(metadata.syn.id,downloadLocation="./",ifcollision = "overwrite.local"))
if( is.null(md.link) ) {
  cat("ERROR: Unable to download metadata table from Synapse\n")
  cat("For information on about to setup the credentials check https://annaisystems.zendesk.com/hc/en-us/article_attachments/201188358/Pan-Cancer_Researcher_Guide.pdf")
  q(status=1)
}
cat("Downloaded:",getFileLocation(md.link), "\n")
md.file <- basename(getFileLocation(md.link))
cat("Downloading data from Synapse...done.\n")

metadata <- NULL
try(metadata <- read.table(md.file,header=T,check.names=F,sep="\t",quote="\"",comment.char=""))
if(is.null(metadata)) {
  cat("Unable to load ",md.file)
  q(status=1)
}
metadata$library_type[is.na(metadata$library_type)] <- "fr-unstranded"
rownames(metadata) <- metadata$analysis_id
cat("Downloading and reading metadata...done.\n")
cat(paste("Read ",nrow(metadata)," rows in ",metadata.file,"\n",sep=""))


# donors per project
ss <- unique(metadata$submitted_donor_id)
cat("Donors: "); cat(length(ss)); cat ("\n");


#####
# QC
qc <- read.table(qc.file,sep="\t",header=T,check.names=F)
cat("Read ",nrow(qc)," rows in ",qc.file,"\n")
#dim(qc)
rownames(qc) <- qc[,1]
qc <- qc[,-1]

samples.in.qc <- metadata[unique(gsub("/.*","",colnames(qc))),]
cat("samples in QC ",nrow(samples.in.qc),"\n")
# 2202
normal <- sum(samples.in.qc$is_tumour=="no",na.rm=T)
#normal
tumour <- sum(samples.in.qc$is_tumour=="yes",na.rm=T)
#tumour

####################################
# what was the correlation threshold 0.95
# what was the criteria 3/bias
# some analysis id appear multiple times with the  suffix {1,2,3...}
# to which library does the id correspond to?
vars <- grep("_fail",rownames(qc),value=T)
print(vars)
flagged.libs <- list()
flagged.sum <- list()
for ( crit in vars ) {
  x <- as.character(unlist(qc[crit,]))
  names(x) <- colnames(qc)
  flagged.sum[[crit]] <- unlist(table(x))
  flagged.libs[[crit]] <- names(x[x=="True"])
}

flagged.libs.v <- unique(unlist(flagged.libs))
cat("Libraries flagged:",length(flagged.libs.v),"\n")
# 267
perc.libs.flagged <- length(flagged.libs.v)/ncol(qc)
cat("Libraries flagged (%):",perc.libs.flagged,"\n")

# create a subset of the initial matrix
qc.sel <- qc[vars,flagged.libs.v]
qc.bool <- qc.sel=="True"
qc.vot <- colSums(qc.bool)

####################################################
# Blacklists - libraries, samples and donors
# libs excluded
black.list.libs <- names(qc.vot[qc.vot>=2])
num.libraries.black.listed <- length(black.list.libs)
cat("Libraries blacklisted:",num.libraries.black.listed,"\n")

# analysis ids
black.list.libs2 <- gsub("/.*","",black.list.libs)
x <- table(black.list.libs2)
cat("Samples with blacklisted libraries:",length(x),"\n")

# expected number of libs per sample
v <- c()
for ( i in seq(1,nrow(metadata))) {
  a <- as.character(metadata$rg_label[i])
  v <- append(v,length(unlist(strsplit(a,split=","))))
}
names(v) <- metadata$analysis_id

analysis.ids <- v[names(x)]-x
samples.black.listed <- names(analysis.ids[analysis.ids==0])
#length(samples.black.listed)

a <- table(as.character(metadata[samples.black.listed,"submitted_donor_id"]))
b <- table(as.character(metadata[metadata$submitted_donor_id%in% names(a),"submitted_donor_id"]))
ab <- a/b
cat("Donors blacklisted:",length(ab[ab==1]),"\n")
black.list.donors <- unique(names(ab[ab==1]))

###################################################
# White list - libraries, samples, and donors
qc.white <- qc[vars,!colnames(qc) %in% flagged.libs.v]
white.list.libs <- names(qc.white)
num.libraries.white.listed <- length(white.list.libs)
cat("Libraries whitelisted:",num.libraries.white.listed,"\n")

# analysis ids
white.list.libs2 <- gsub("/.*","",white.list.libs)
x <- table(white.list.libs2)
cat("Samples with whitelisted libraries:",length(x),"\n")

# expected number of libs per sample
v <- c()
for ( i in seq(1,nrow(metadata))) {
  a <- as.character(metadata$rg_label[i])
  v <- append(v,length(unlist(strsplit(a,split=","))))
}
names(v) <- metadata$analysis_id

analysis.ids <- v[names(x)]-x
samples.white.listed <- names(analysis.ids[analysis.ids==0])

a <- table(as.character(metadata[samples.white.listed,"submitted_donor_id"]))
b <- table(as.character(metadata[metadata$submitted_donor_id%in% names(a),"submitted_donor_id"]))
ab <- a/b
cat("Donors whitelisted:",length(ab[ab==1]),"\n")
white.list.donors <- unique(names(ab[ab==1]))


# Donors
write.table(white.list.donors,file=paste("donors_white_list.tsv",sep=""),quote=F,row.names=F,col.names=F)
write.table(black.list.donors,file=paste("donors_black_list.tsv",sep=""),quote=F,row.names=F,col.names=F)

# Samples
write.table(metadata[samples.white.listed,c("analysis_id","submitted_donor_id")],file=paste("samples_white_list.tsv",sep=""),quote=F,row.names=F,col.names=F)
write.table(metadata[samples.black.listed,c("analysis_id","submitted_donor_id")],file=paste("samples_black_list.tsv",sep=""),quote=F,row.names=F,col.names=F)


# Libraries
df <- list()
df$libs <- colnames(qc)
df$samples <- gsub("/.*","",colnames(qc))
df$donors <- metadata[df$samples,"donor_id"]
df <- data.frame(df)
rownames(df) <- df$libs
write.table(df[white.list.libs,],file=paste("libs_white_list.tsv",sep=""),quote=F,row.names=F,col.names=F)
write.table(df[black.list.libs,],file=paste("libs_black_list.tsv",sep=""),quote=F,row.names=F,col.names=F)

cat("All files created.\n")

if ( upload == "yes" ) {
  thisCode <- getPermlink(getRepo("Functional-Genomics/pancancer"), "scripts/pcawg_gen_qc_lists.R")

  cat("Uploading files to synapse...")
  qc.folder <- synGet(qc.folder.syn.id)
  # samples
  f1 <- synStore(File(path="samples_white_list.tsv", parentId=qc.folder$properties$id),
                 activityName="Sample - whitelist",
                 used=list(
                   list(url=thisCode, name=basename(thisCode), wasExecuted=T),
                   list(entity=qc.link, wasExecuted=F),
                   list(url=metadata.file, wasExecuted=F)))

  f2 <- synStore(File(path="samples_black_list.tsv", parentId=qc.folder$properties$id),
                activityName="Sample - blacklist",
                used=list(
                  list(url=thisCode, name=basename(thisCode), wasExecuted=T),
                  list(entity=qc.link, wasExecuted=F),
                  list(url=metadata.file, wasExecuted=F)))

  f3 <- synStore(File(path="donors_white_list.tsv", parentId=qc.folder$properties$id),
                activityName="Donor - whitelist",
                used=list(
                  list(url=thisCode, name=basename(thisCode), wasExecuted=T),
                  list(entity=qc.link, wasExecuted=F),
                  list(url=metadata.file, wasExecuted=F)))

  f4 <- synStore(File(path="donors_black_list.tsv", parentId=qc.folder$properties$id),
                activityName="Donor - blacklist",
                used=list(
                  list(url=thisCode, name=basename(thisCode), wasExecuted=T),
                  list(entity=qc.link, wasExecuted=F),
                  list(url=metadata.file, wasExecuted=F)))

}
cat("All done,bye.\n")
q(status=0)




