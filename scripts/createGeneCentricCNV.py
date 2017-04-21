#!/usr/bin/env python 
import scipy as sp

import pdb
import h5py
import sys
import os
import glob
import fnmatch
import pandas as pd

from intervaltree import Interval, IntervalTree

def readcnv(fn, intervalltree):
    cnvdata = sp.loadtxt(fn, delimiter = '\t', dtype = 'string')
    cnvdata = cnvdata[1:,:]
    if 'X' in cnvdata[:,0]:
        bv = cnvdata[:,0] == 'X'
        cnvdata = cnvdata[~bv] #exclude CNV on chrX
    if 'Y' in cnvdata[:,0]:
        bv = cnvdata[:,0] == 'Y'
        cnvdata = cnvdata[~bv] #exclude CNV on chrY
    dict_results = {}
    for i,c in enumerate(cnvdata):
        #if int(c[-1]) <= 1: ###
            #continue
	if c[3] == 'NA':
            continue # skip if total count are NA
        if int(c[3]) == 2: ### normal CNV
            continue
        gene_region =  intervalltree[int(eval(c[0]))][ int(eval(c[1])): int(eval(c[2]))]
        if len(gene_region) == 0:
            continue
        cnvgenes = sp.array([x.data[1].split(' ')[-1] for x in gene_region])
        cnvgenes = sp.array([x.strip('\'').strip('\"') for x in cnvgenes])
        for g in cnvgenes:
            if g in dict_results:
                dict_results[g].append(int(c[3]))
            else:
                dict_results[g] = [int(c[3])]
    gid = []
    cnv = []
    for g in dict_results.keys():
        gid.append(g)
        cnv.append(sp.mean(dict_results[g]))
        #cnv.append(dict_results[g][0])
    return gid, cnv
#            if r.data[1] in dict_results:
                

if __name__ == "__main__":
    if len(sys.argv[1:]) < 5:
        sys.stderr.write('ERROR: missing argument\n')
        sys.stderr.write('Usage: script.py <gencode.gtf> <mapping file> <CNVs dir> <outfilename> <samples not kept.tsv>\n')
        sys.exit(1)    
    ### gencode file in gtf format
    fn_anno = sys.argv[1]
    ### aliquot list in tsv (claudia format)
    fn_aliquotlist =sys.argv[2] 
    ### directory of cnv calls (one aliquot per file)
    dir_cnv = sys.argv[3] 
    ### out file
    fn_out = sys.argv[4] 
    fn_out_not_kept = sys.argv[5]


    ### read annotation
    anno_data  = sp.loadtxt(fn_anno, delimiter = '\t', dtype = 'string', usecols=[0,2,3,4,8])
    # filter to gene annotation only
    anno_data = anno_data[anno_data[:,1] == 'gene',:]
    # remove X/Y/MT (since we do not care in QTL)
    #iAutosome = sp.array([x[0].isdigit() for x in anno_data])
    i_notAutosome= map(lambda x:x in ['chrM','chrX','chrY'], anno_data[:,0])
    anno_data = anno_data[~sp.array(i_notAutosome),:]
    print '{0} genes kept from the annotation file\n'.format(anno_data.shape[0])
    #anno_data = anno_data[iAutosome,:]
    
    ## create dictionary of interval trees
    intervalltree={}
    for i in xrange(1,23):
        anno_data_chr = anno_data[anno_data[:,0] == 'chr'+str(i),:]
        ttree = IntervalTree()
        for gn in anno_data_chr:
            ttree[int(gn[2]):int(gn[3])] = [gn[0], gn[4].split(';')[0]]
        intervalltree[i] = ttree



    ### load aliquot data
    aq_data = sp.loadtxt(fn_aliquotlist, delimiter = '\t', dtype = 'string', usecols = [2])
    aq_data = aq_data[1:]
    

    ### get aliquot file list
    cnv_files = sp.array(os.listdir(dir_cnv))

    allgids = sp.array([x.split(';')[0].split(' ')[1].strip('\"') for x in anno_data[:,4]])
    allgids = sp.sort(allgids)
    cnvdata = sp.zeros(( aq_data.shape[0], anno_data.shape[0] ))
    #cnvdata[:] = sp.nan
    cnvdata[:] = 2 #set all the rest to euploidy
    counter = 0
    aliquot_not_kept=[]
    for i,aq in enumerate(aq_data): ### read each aliquot
        #print aq
        #ix_f = sp.array([x.startswith(aq) for x in cnv_files])
        filename = glob.glob(dir_cnv+aq+'*.gz')
        if len(filename)==1:
	    print '{0} found.'.format(aq)
            counter+=1
        else:
            aliquot_not_kept.append(aq)
            continue
#        if not ix_f.any():
 #           continue

        #gid, cnv = readcnv(os.path.join(dir_cnv, cnv_files[ix_f][0]), intervalltree)
        gid,cnv = readcnv(filename[0],intervalltree)
        gid = sp.array(gid)
        cnv = sp.array(cnv)
    
        sidx = sp.argsort(gid)
        gid  = gid[sidx]
        cnv = cnv[sidx]
        
        midx = sp.in1d(allgids, gid)
        cnvdata[i,midx] = cnv 

    #print number of aliquot IDs analysed
    print 'number of aliquot ids analysed: {0}\n'.format(counter)
    ### OUTPUT statement
    OUT = h5py.File(fn_out, 'w')
    #OUT.create_dataset(name = 'cnv', data = cnvdata.T)
    OUT.create_dataset(name = 'cnv',data=cnvdata[:])
    OUT.create_dataset(name = 'gene_name', data = allgids)
    OUT.create_dataset(name = 'wgs_aliquot_id', data = aq_data)
    OUT.close()
    out = open(fn_out_not_kept,'w')
    for a in aliquot_not_kept:
        out.write(a+'\n')
    out.close()
    print 'Done.'
    sys.exit(0)
