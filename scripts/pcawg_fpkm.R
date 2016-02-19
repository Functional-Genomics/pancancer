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

args <- commandArgs(trailingOnly=TRUE)

if ( length(args) != 4 ) {
  cat(" Script to compute FPKMs as defined by the PCAWG group 3.\n")
  cat("ERROR: usage: pcawg_fpkm.R quant_file gene_lengths.tsv protein_coding.tsv out_file\n")
  q(status=1)
}

quant.file <- args[1]
glength.file <- args[2]
pcgenes.file <- args[3]
out.file <- args[4]

for ( f in c(quant.file,glength.file,pcgenes.file) ) {
  if ( ! file.exists(f) ) {
    cat("ERROR: unable to open/find ",f,"\n")
    q(status=1)
  }
}

if (!require("data.table")) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("data.table")
}
library(data.table)

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

#####################################
# Open and load the gene lengths file
tsv.data <- fread(input=glength.file,sep = "\t", header=T,check.names=FALSE,data.table=FALSE)
if ( is.null(tsv.data) ) {
  cat("Error while loading ",glength.file,"\n")
  q(status=1)
}
gene.lengths <- tsv.data[,2]
names(gene.lengths) <- as.character(tsv.data[,1])
cat("Found ",length(gene.lengths)," entries in the file ",glength.file,"\n")


#####################################
# Open the file with the list of protein coding genes
tsv.data <- fread(input=pcgenes.file,sep = "\t", header=T,check.names=FALSE,data.table=FALSE)
if ( is.null(tsv.data) ) {
  cat("Error while loading ",pcgenes.file,"\n")
  q(status=1)
}
pc.genes <- as.character(tsv.data[,1])
cat("Found ",length(pc.genes)," entries in the file ",pcgenes.file,"\n")

######################################
quant <- quant.load(quant.file,header=T)

###############################
# check if the gene lists match
if ( sum(!pc.genes %in% rownames(quant))!=0 ) {
  cat("ERROR: some of the genes in ",quant.file," do not appear in ",pcgenes.file,"\n")
  q(status=1)
}
if ( sum(!  rownames(quant) %in% names(gene.lengths) )!=0 ) {
  cat("ERROR: some of the genes in ",glength.file," do not appear in ",quant.file,"\n")
  q(status=1)
}

###############################################
# (Re) Compute the fpkms
fpkm.repl <- countstable2rpkms(quant,lens=gene.lengths,mass.labels=pc.genes)

# Save the FPKMs to a file
Gene <- rownames(fpkm.repl)
cat("FPKMs computed! Saving matrix to ",out.file,"\n")
write.table(cbind(Gene,fpkm.repl),file=out.file,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

q(status=0)



