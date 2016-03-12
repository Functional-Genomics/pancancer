#!/usr/bin/env Rscript
# =========================================================
# Copyright 2015-2016,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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

# samples excluded
excluded.syn.id <- "syn5648208"
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
library(rGithubClient)
library(RCurl)


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
cat(paste("Read ",nrow(metadata)," rows in ",md.file,"\n",sep=""))


# donors
ss <- unique(metadata$submitted_donor_id)
aids <- unique(metadata$aliquot_id)
cat("Donors: "); cat(length(ss)); cat ("\n");
cat("Aliquots: "); cat(length(aids)); cat ("\n");


############################
# Excluded samples/libraries
excluded.link <- NULL

try(excluded.link<-synGet(excluded.syn.id,downloadLocation="./",ifcollision = "overwrite.local"))
if( is.null(excluded.link) ) {
  cat("ERROR: Unable to download excluded samples table from Synapse\n")
  cat("For information on about to setup the credentials check https://annaisystems.zendesk.com/hc/en-us/article_attachments/201188358/Pan-Cancer_Researcher_Guide.pdf")
  q(status=1)
}
cat("Downloaded:",getFileLocation(excluded.link), "\n")
excluded.file <- basename(getFileLocation(excluded.link))
cat("Downloading data from Synapse...done.\n")

excluded.libs.table <- NULL
try(excluded.libs.table <- read.table(excluded.file,header=T,check.names=F,sep="\t",quote="\"",comment.char=""))
if(is.null(excluded.libs.table)) {
  cat("Unable to load ",excluded.file)
  q(status=1)
}
cat("Found ",nrow(excluded.libs.table)," in ",excluded.file,"\n")

excluded.samples <- as.character(excluded.libs.table$analysis_id[as.character(excluded.libs.table$fastq_file)=="all"])
cat("Samples excluded ",length(excluded.samples),"\n")
#excluded.samples
cat("Downloading and reading excluded samples...done.\n")


#####
# QC
qc <- read.table(qc.file,sep="\t",header=T,check.names=F)
cat("Read ",nrow(qc)," rows in ",qc.file,"\n")
#dim(qc)
rownames(qc) <- qc[,1]
qc <- qc[,-1]
colnames(qc) <- unlist(qc["orig_id_tophat",])

samples.in.qc <- metadata[unique(gsub("/.*","",colnames(qc))),]
cat("samples in QC ",nrow(samples.in.qc),"\n")

#rownames(qc)
# transpose the qc matrix (libs by feature)
t.qc <- t(qc)

####################################
# what was the correlation threshold 0.95
# what was the criteria 3/bias
vars <- grep("_fail",colnames(t.qc),value=T)
print(vars)
flagged.libs <- list()
flagged.sum <- list()
#crit <- vars[1]
for ( crit in vars ) {
  x <- as.character(unlist(t.qc[,crit]))
  names(x) <- rownames(t.qc)
  flagged.sum[[crit]] <- unlist(table(x))
  flagged.libs[[crit]] <- names(x[x=="True"])
}

flagged.libs.v <- unique(unlist(flagged.libs))
cat("Libraries flagged:",length(flagged.libs.v),"\n")
# 
perc.libs.flagged <- round(length(flagged.libs.v)/ncol(qc)*100,2)
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

# check if the results are consistent
previously.blacklisted <- rownames(t.qc)[as.character(t.qc[,"blacklist"])=="True"]
if ( length(black.list.libs)==length(previously.blacklisted)) {
  if ( sum(!black.list.libs%in%previously.blacklisted)!=0 ) {
    cat("ERROR: There is a mismatch between the number of blacklisted libraries")
    q(status=1)
  }
} else {
  cat("ERROR: There is a mismatch between the number of blacklisted libraries")
  q(status=1)
}

# expected number of libs per sample (based on the metadata)
# single read group label (rg_label) for "057da4ba-421e-4f39-afa8-c7de2ca665e2"
# go back to rg_label
v <- c()
for ( i in seq(1,nrow(metadata))) {
  a <- as.character(metadata$rg_label[i])
  v <- append(v,length(unlist(strsplit(a,split=" "))))
}
names(v) <- as.character(metadata$analysis_id)
expLibsPerSample <- v

analysis.ids <- expLibsPerSample[names(x)]-x
samples.black.listed <- names(analysis.ids[analysis.ids==0])
cat("Samples blacklisted:",length(samples.black.listed),"\n")
cat("Samples excluded:",length(excluded.samples),"\n")

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

# check if the results are consistent
previously.whitelisted <- rownames(t.qc)[as.character(t.qc[,"whitelist"])=="True"]
if ( length(white.list.libs)==length(previously.whitelisted)) {
  if ( sum(!white.list.libs%in%previously.whitelisted)!=0 ) {
    cat("ERROR: There is a mismatch between the number of white listed libraries")
    q(status=1)
  }
} else {
  cat("ERROR: There is a mismatch between the number of white listed libraries")
  q(status=1)
}

# analysis ids
white.list.libs2 <- gsub("/.*","",white.list.libs)
x <- table(white.list.libs2)
white.list.libs3 <- x[!names(x)%in%excluded.samples]
cat("Samples with whitelisted libraries:",length(x),"\n")
cat("Samples with whitelisted libraries:",length(white.list.libs3),"\n")

#
analysis.ids <- expLibsPerSample[names(x)]-x

samples.white.listed <- names(analysis.ids[analysis.ids==0])
samples.white.listed <- samples.white.listed[!samples.white.listed%in%excluded.samples]
samples.white.listed <- samples.white.listed[!samples.white.listed%in%samples.black.listed]
cat("Samples whitelisted:",length(samples.white.listed),"\n")

#
a <- table(as.character(metadata[samples.white.listed,"submitted_donor_id"]))
b <- table(as.character(metadata[metadata$submitted_donor_id%in% names(a),"submitted_donor_id"]))
ab <- a/b
cat("Donors whitelisted:",length(ab[ab==1]),"\n")
white.list.donors <- unique(names(ab[ab==1]))

###################################################
# Donors
write.table(white.list.donors,file=paste("donors_white_list.tsv",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
write.table(black.list.donors,file=paste("donors_black_list.tsv",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

# Samples
write.table(metadata[samples.white.listed,c("analysis_id","submitted_donor_id")],file=paste("samples_white_list.tsv",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
write.table(metadata[samples.black.listed,c("analysis_id","submitted_donor_id")],file=paste("samples_black_list.tsv",sep=""),quote=F,row.names=F,col.names=F,sep="\t")


# Libraries
df <- list()
df$libs <- colnames(qc)
df$samples <- gsub("/.*","",colnames(qc))
df$donors <- metadata[df$samples,"donor_id"]
df <- data.frame(df)
rownames(df) <- df$libs
# remove the excluded samples
df <- df[!df$samples %in% excluded.samples,]

write.table(df[rownames(df)%in%white.list.libs,],file=paste("libs_white_list.tsv",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
write.table(df[rownames(df)%in%append(black.list.libs,excluded.samples),],file=paste("libs_black_list.tsv",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

cat("All files created.\n")

if ( upload == "yes" ) {
  thisCode <- getPermlink(getRepo("Functional-Genomics/pancancer"), "scripts/pcawg_gen_qc_lists.R")
  cat("Uploading files to synapse...")
  qc.folder <- synGet(qc.folder.syn.id)
  # samples
  f1 <- synStore(File(path="samples_white_list.tsv", parentId=qc.folder$properties$id,name="Sample whitelist"),
                 
                 activityName="Sample - whitelist",
                 used=list(
                   list(url=thisCode, name=basename(thisCode), wasExecuted=T),
                   list(entity=qc.link, wasExecuted=F),
                   list(entity=excluded.link,wasExecuted=F),
                   list(entity=md.link, wasExecuted=F)))

  f2 <- synStore(File(path="samples_black_list.tsv", parentId=qc.folder$properties$id,name="Sample blacklist"),
                activityName="Sample - blacklist",
                used=list(
                  list(url=thisCode, name=basename(thisCode), wasExecuted=T),
                  list(entity=qc.link, wasExecuted=F),
                  list(entity=excluded.link,wasExecuted=F),
                  list(entity=md.link, wasExecuted=F)))

  f3 <- synStore(File(path="donors_white_list.tsv", parentId=qc.folder$properties$id,name="Donor whitelist"),
                activityName="Donor - whitelist",
                used=list(
                  list(url=thisCode, name=basename(thisCode), wasExecuted=T),
                  list(entity=qc.link, wasExecuted=F),
                  list(entity=excluded.link,wasExecuted=F),
                  list(entity=md.link, wasExecuted=F)))

  f4 <- synStore(File(path="donors_black_list.tsv", parentId=qc.folder$properties$id,name="Donor blacklist"),
                activityName="Donor - blacklist",
                used=list(
                  list(url=thisCode, name=basename(thisCode), wasExecuted=T),
                  list(entity=qc.link, wasExecuted=F),
                  list(entity=excluded.link,wasExecuted=F),
                  list(entity=md.link, wasExecuted=F)))

}
cat("All done,bye.\n")
q(status=0)




