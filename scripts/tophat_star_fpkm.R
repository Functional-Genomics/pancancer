#!/usr/bin/env Rscript
# =========================================================
# Copyright 2016,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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


# gene expression
data.in.synapse <- list()
data.in.synapse$tophat2.fpkm.syn <- list(syn.id="syn5553983",
                                 file="tophat2.pcawg.fpkm.tsv.gz")
  
data.in.synapse$star.fpkm.syn<- list(syn.id="syn5553981",
                                 file="star.pcawg.fpkm.tsv.gz")

data.in.synapse$tophat2.raw.syn <- list(syn.id="syn5522827",
                                 file="tophat2.raw.tsv.gz")

data.in.synapse$star.raw.syn <- list(syn.id="syn3159796",
                                      file="star.raw.tsv.gz")


########################################################
args <- commandArgs(trailingOnly=TRUE)

if ( length(args) != 4 ) {
  cat("This script will download data from synapse.\n")
  cat("ERROR: usage: tophat_star_fpkm.R synapse_login synapse_password gene_lengths.tsv protein_coding.tsv\n")
  q(status=1)
}

login <- args[1]
pass <- args[2]
glength.file <- args[3]
pcgenes.file <- args[4]

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
if (!require("data.table")) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("data.table")
}
library(data.table)

cat("Connecting to Synapse...\n")
synapseLogin(login,pass)
cat("Connecting to Synapse...done.\n")


cat("Downloading data from Synapse...\n")
for (n in names(data.in.synapse) ) {
  cat("<<<< ",n,"\n")
  try(syn.file<-synGet(data.in.synapse[[n]]$syn.id,downloadFile=T,downloadLocation="./",ifcollision = "overwrite.local",))
  if( is.null(syn.file) ) {
    cat("For information on about to setup the credentials check https://annaisystems.zendesk.com/hc/en-us/article_attachments/201188358/Pan-Cancer_Researcher_Guide.pdf")
    cat("ERROR: Unable to download ",data.in.synapse[[n]]$syn.id,"\n")
    q(status=1)
  }  
  system(paste("mv ",getFileLocation(syn.file)," ",data.in.synapse[[n]]$file,"\n",sep=""))
}
cat("Downloading data from Synapse...done\n")
# All files downloaded

#####################################
# Open and load the gene lengths file
tsv.data <- fread(input=glength.file,sep = "\t", header=T,check.names=FALSE,data.table=FALSE)
gene.lengths <- tsv.data[,2]
names(gene.lengths) <- as.character(tsv.data[,1])
head(gene.lengths)


#####################################
# Open the file with the list of protein coding genes
tsv.data <- fread(input=pcgenes.file,sep = "\t", header=T,check.names=FALSE,data.table=FALSE)
pc.genes <- as.character(tsv.data[,1])
head(pc.genes)

######################################
# Code from iRAP
countstable2rpkms <- function(table,lens,mass.labels=NULL,exitonerror=TRUE) {
  # check if there missing features
  missing.feat <- (!rownames(table) %in% names(lens))

  if ( sum(missing.feat) ) {
    message("ERROR: Length of ",paste(rownames(table)[missing.feat],sep=",")," not found.")
    if (exitonerror) { q(status=1) }
    return(NULL)
  }
  v.compute.rpkm <-  function(l,lens,mass.labels=NULL) {
    #(l*1e6)/(sum(l)*lens[names(l)]/1000)

    if ( is.null(mass.labels) ) {
      tot.mass <- sum(l)   
    } else {
      tot.mass <- sum(l[mass.labels])
      if ( tot.mass == 0 ) {
        print("Tot.mass==0!!!!")
        tot.mass <- 1
      }
      message("Tot.mass:",tot.mass)
    }
    return(10^9*l/(tot.mass*lens[names(l)]))
  }
  if ( is.vector(table) ) {
    return(round(v.compute.rpkm(table,lens,mass.labels),2))
  } else {
    return(round(apply(table,2,v.compute.rpkm,lens,mass.labels),2))
  }
}
quant.load <- function(f,header=T) {
  tsv.data <- NULL
  if ( ! file.exists(f) ) {
    stop("File does not exist:",f,"\n")
  }
  if ( sum(grep(".gz$",f)) ) {
    f <- paste("zcat ",f,sep="")
  }
  tsv.data <- fread(input=f,sep = "\t", header=header,check.names=FALSE,data.table=FALSE)
  if ( !is.null(tsv.data) ) {
    rownames(tsv.data) <- as.character(tsv.data[,1])
    tsv.data <- tsv.data[,-1,drop=FALSE]
    return(tsv.data)
  }
  return(NULL)
}

#####################
quant.file <- "tophat2.raw.tsv.gz"
th2.quant <- quant.load(quant.file,header=T)
# % of reads mapping PC genes
tot.expr <- colSums(th2.quant)
pc.expr <- colSums(th2.quant[pc.genes,])
sum(pc.genes %in% rownames(th2.quant))
th.pc.expr <- round(pc.expr/tot.expr*100,2)

#quant.file <- "star.raw.tsv.gz"
#star.quant <- quant.load(quant.file,header=T)
#star.quant <- star.quant[grepl("ENS*",rownames(star.quant)),]

# remove non-genes from the table
# % of reads mapping PC genes
#tot.expr <- colSums(star.quant)
#pc.expr <- colSums(star.quant[pc.genes,])
#sum(pc.genes %in% rownames(star.quant))
#star.pc.expr <- round(pc.expr/tot.expr*100,2)

#png("th2_star_pc_expr_bp.png")
#boxplot(list("Topha2"=th.pc.expr,"Star"=star.pc.expr))
#dev.off()

png("th2_pc_expr_bp.png")
par(bty="n")
boxplot(list("Topha2"=th.pc.expr),main="Reads mapping protein coding genes per aliquot_id",ylab="%")
dev.off()

###################
# (Re) Compute the fpkms
th2.fpkm.repl <- countstable2rpkms(th2.quant,lens=gene.lengths,mass.labels=pc.genes)
th2.fpkm.orig <- countstable2rpkms(th2.quant,lens=gene.lengths)

# Save the FPKMs to a file
Gene <- rownames(th2.fpkm.repl)
write.table(cbind(Gene,th2.fpkm.repl),file="th2_fpkm.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

Gene <- rownames(th2.fpkm.orig)
write.table(cbind(Gene,th2.fpkm.orig),file="th2_fpkm_orig.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
################
# Official FPKMS

# FPKMs
th2.fpkms.pub <- quant.load("tophat2.pcawg.fpkm.tsv.gz",header=T)
# same order of columns
th2.fpkms.pub <- th2.fpkms.pub[,colnames(th2.fpkm.orig)]
# FPKMs_UQ
th2.fpkms_uq.pub <- quant.load("tophat2.pcawg.fpkm_uq.tsv.gz",header=T)
# same order of columns
th2.fpkms_uq.pub <- th2.fpkms_uq.pub[,colnames(th2.fpkm.orig)]

png("FPKMS.png",width=1000,height=1000,res=150)
par(bty="n",mfrow=c(2,2))
for (sample in c(1,10,200,500) ) {
  sample.name <- colnames(th2.fpkm.orig)[sample]
  cat("----------------------------------------------------\n")
  cat(sample.name)
  boxplot(list("FPKM\nPC"=th2.fpkm.repl[,sample]+1,
               "FPKM\nVanilla"=th2.fpkm.orig[,sample]+1,
               "Official\nFPKM"=th2.fpkms.pub[,sample]+1,
               "Official\nFPKM-UQ"=th2.fpkms_uq.pub[,sample]+1),
          log="y",main=sample.name,cex=0.7,las=2)
  cat("FPKM PC\n")
  print(summary(th2.fpkm.repl[,sample]))
  cat("FPKM\n")
  print(summary(th2.fpkm.orig[,sample]))
  cat("Official FPKM\n")
  print(summary(th2.fpkms.pub[,sample]))
  cat("Official FPKM-UQ\n")
  print(summary(th2.fpkms_uq.pub[,sample]))
}
dev.off()

q(status=0)



