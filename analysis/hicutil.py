import argparse
import numpy as np
import csv
import math
import glob
import sys
import scipy.stats

def readLocalDiffFiles(filelist,celltypelist):
    compdata = {}
    for idx,filename in enumerate(filelist):
        chrloc = filename.find('_chr')
        chrnum = int(filename[chrloc+4:-4])
        fileparts = filename.split('_')
        celltype = fileparts[1]
        for i in xrange(3,len(fileparts)):
            ct = '_'.join(fileparts[1:i])
            if ct in celltypelist:
                celltype = ct
        for fsplit in fileparts:
            if 'replicates' not in fsplit and 'rep' in fsplit and fsplit[-1] != 's':
                celltype = celltype + fsplit[-1]
        #print filename, chrnum, celltype
        dictkey = (celltype, chrnum)
        dictdata = []
        with open(filename, 'rb') as f:
            freader = csv.reader(f, delimiter='\t')
            freader.next() #skip header line
            for line in freader:
                if not math.isnan(float(line[2])):
                    dictdata.append([float(x) for x in line])
            dictdata = np.array(dictdata)
        compdata[dictkey] = dictdata
    #sys.exit()
    return compdata

def readNonRepLocalDiffFiles(filelist,celltypelist):
    compdata = {}
    for idx,filename in enumerate(filelist):
        chrloc = filename.find('_chr')
        chrnum = int(filename[chrloc+4:-4])
        fileparts = filename.split('_')
        celltype1 = fileparts[1]
        ct1end = 2
        for i in xrange(3,len(fileparts)):
            ct = '_'.join(fileparts[1:i])
            if ct in celltypelist:
                celltype1 = ct
                ct1end = i
        celltype2 = fileparts[ct1end]
        for j in xrange(ct1end+2,len(fileparts)):
            ct = '_'.join(fileparts[ct1end:j])
            if ct in celltypelist:
                celltype2 = ct
        #print filename, chrnum, celltype
        dictkey = ((celltype1,celltype2), chrnum)
        dictdata = []
        with open(filename, 'rb') as f:
            freader = csv.reader(f, delimiter='\t')
            freader.next() #skip header line
            for line in freader:
                if not math.isnan(float(line[2])):
                    dictdata.append([float(x) for x in line])
            dictdata = np.array(dictdata)
        compdata[dictkey] = dictdata
    return compdata

def calcPercSimilarity(compdata, chrlengths):
    # compdata: dictionary with ((celltype1,celltype2),chrnum) where each entry is the list of significant, dominating intervals [[start,end,VI,pval],...]

    totallen = np.sum(chrlengths)
    totalsim = {} #will be dictionary w/ key = celltypepair, value = %sim
    for key,data in compdata.iteritems():
        chrnum = key[1]
        if len(data) > 0:
            sigints = np.zeros(int(max(data[:,1])))
            for row in data:
                sigints[int(row[0]):int(row[1])] = 1
            chrsum = np.sum(sigints)
        else:
            chrsum = 0
        celltype = key[0]
        if celltype in totalsim:
            totalsim[celltype] += float(chrsum)/totallen
        else:
            totalsim[celltype] = float(chrsum)/totallen
        #idx = -1
        #for i,sublist in enumerate(totalsim):
        #    if sublist[0] == celltype:
        #        idx = i
        #if idx > -1:
        #    totalsim[idx][1] += float(chrsum)/totallen
        #else:
        #    totalsim.append([celltype, float(chrsum)/totallen])
    return totalsim

def writePercSimToFile(totalsim, filename):

    with open(filename,'wb') as f:
        fwriter = csv.writer(f, delimiter='\t')
        for row in totalsim:
            fwriter.writerow([row[0][0], row[0][1], row[1]])
    print 'wrote percent similarity values to', filename

def calcBinCons(compdata, chrlengths, res = 100000):
    # format data for Figure 4, calculate average conservation per bin

    chrdict = {} # keys will be chromosome numbers
    for key,data in compdata.iteritems(): # compdata[(celltype1,celltype2,chrnum)] = [[start,end,VI,qval]...]
        chrnum = key[1]
        # for each chromosome, make a matrix with a row for each key, and columns for each bin, indicator if bin is included in overlaps
        if chrnum not in chrdict:
            chrmat = np.zeros((len(compdata),chrlengths[chrnum-1]))
        else:
            chrmat = chrdict[chrnum]
        rowidx = np.where(~chrmat.any(axis=1))[0][0]
        for interval in data:
            chrmat[rowidx,int(interval[0]):int(interval[1])] = 1
        chrdict[chrnum] = chrmat

    bincons = [[] for i in range(22)] # initialize empty list of length 22 for each chromosome #
    #avgperbin = 0
    for chrnum,chrmat in chrdict.iteritems():
        bincons[chrnum-1] = np.sum(chrmat, axis=0)/float(len(compdata)/float(22))
        #avgperbin += np.sum(bincons[chrnum-1])
    #print 'Average conservation per bin = ', avgperbin/np.sum(chrlengths)
    return bincons
