import argparse
import numpy as np
import csv
import math
import glob
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.stats
from sklearn import manifold, cluster, metrics
import umap
from scipy.spatial.distance import pdist, squareform
from sklearn import datasets
from scipy.cluster.hierarchy import linkage
import hicutil
import tadanalysisutil as tad

newpalette = ['#8dd3c7','#58508d','#bc5090','#003f5c','#ff6361']
sns.set(style='whitegrid',palette=sns.color_palette(newpalette),color_codes=False,font='Ubuntu')

def readAllFiles(bgfileseed,repfileseed,donorfileseed,celltypelist):

    with open(celltypelist,'rb') as f:
        freader = csv.reader(f)
        ctlist = []
        for line in freader:
            ctlist.extend(line)

    bgfilelist = []
    for fileseed in bgfileseed:
        bgfilelist += glob.glob(fileseed)
    bgdata = hicutil.readNonRepLocalDiffFiles(bgfilelist, ctlist)

    repfilelist = glob.glob(repfileseed)
    repdata = hicutil.readRepLocalDiffFiles(repfilelist, ctlist)

    donorfilelist = glob.glob(donorfileseed)
    donordata = hicutil.readNonRepLocalDiffFiles(donorfilelist, ctlist)

    return bgdata,repdata,donordata,ctlist

def resFragAnalysis(alltadsim,measure,outputloc):

    resfragdata = [['hESC_DpnII_rep', 'hESC_HindIII', 'hESC_NcoI'], ['HFF-hTERT_DpnII', 'HFF-hTERT_HindIII_beads', 'HFF-hTERT_HindIII_plate', 'HFF-hTERT_MboI', 'HFF-hTERT_NcoI']]
    resfragscores = [[],[]]
    background = []
    for key, score in alltadsim.iteritems():
        if (key[0] in resfragdata[0] and key[1] in resfragdata[0]):
            resfragscores[0].extend([score])
            #print key, score
        elif (key[0] in resfragdata[1] and key[1] in resfragdata[1]):
            resfragscores[1].extend([score])
            #print key, score
        else:
            background.extend([score])

    tad.plotResFragResults(resfragscores, background, measure, outputloc+'tadsim_resfragboxplot_violin.pdf')

def acrossLabAnalysis(alltadsim,reptadsim,measure,outputloc):

    acrosslabdata = [['hESC_Dixon2015','hESC_FP','hESC_HindIII','hESC_Jin','hESC'],['IMR90dixon', 'IMR90rao', 'IMR90_normal', 'IMR90_RI']]
    acrosslabscores = [[],[]]
    background = []
    for key,score in alltadsim.iteritems():
        if (key[0] in acrosslabdata[0] and key[1] in acrosslabdata[0]):
            acrosslabscores[0].extend([score])
            #print key, acrosslabscores[0][-1]
        elif (key[0] in acrosslabdata[1] and key[1] in acrosslabdata[1]):
            acrosslabscores[1].extend([score])
            #print key,acrosslabscores[1][-1]
        else:
            background.extend([score])

    tad.plotAcrossLabResults(acrosslabscores, background, reptadsim.values(), measure, outputloc+'tadsim_acrosslabboxplot.pdf')

def insituDilutionAnalysis(alltadsim,reptadsim,measure,outputloc):

    dil = ['A549','Caki2','G401','LNCaP-FGC','NCI-H460','Panc1','RPMI-7951','SKMEL5','SKNDZ','SKNMC','T47D','IMR90dixon','hESC','hESC_HindIII','HFF-hTERT_HindIII_techreps','HG00733','HG00732','HG00731','HG00514','HG00513','HG00512','GM19238','GM19239','GM19240','hESC_Jin','IMR90_flav','IMR90_normal','IMR90_tnfa','hESC_Dixon2015','GM20431','BrainMicroEndo','AstrocyteCerebellum','AstrocyteSpinalCord','DLD1','BrainPericyte','EndomMicroEndoth','HepSin','ACHN','IMR90_RI','hESC_FP','Adrenal','Bladder','DPC','Hippocampus','Lung','Ovary','Pancreas','Psoas','RightVentricle','SmallBowel','Spleen']
    insitu = ['HFF-hTERT_HindIII_beads','HFF-hTERT_HindIII_plate', 'IMR90rao','GM12878','HMEC','HUVEC','K562rao','KBM7','NHEK','HFF-hTERT_MboI','SkelMuscle','TransColon','hESC_NcoI','HFF-hTERT_NcoI','hESC_DpnII', 'hESC_DpnII_rep','HFF-hTERT_DpnII','HFF-c6']

    ispairs = []
    dilpairs = []
    isdil = []
    samecelltypediffexp = []
    samecelltypesameexp = [[],[]] # 1 list for in situ, 1 for dilution
    for key,score in alltadsim.iteritems():
        if (key[0] in dil and key[1] in dil):
            dilpairs.extend([score])
            if '_flav' in key[0] or '_flav' in key[1] or '_tnfa' in key[0] or '_tnfa' in key[1]: continue
            if ('IMR90' in key[0] and 'IMR90' in key[1]) or ('hESC' in key[0] and 'hESC' in key[1]) or ('HFF-hTERT' in key[0] and 'HFF-hTERT' in key[1]):
                samecelltypesameexp[1].extend([score])
        elif key[0] in insitu and key[1] in insitu:
            ispairs.extend([score])
            if '_flav' in key[0] or '_flav' in key[1] or '_tnfa' in key[0] or '_tnfa' in key[1]: continue
            if ('IMR90' in key[0] and 'IMR90' in key[1]) or ('hESC' in key[0] and 'hESC' in key[1]) or ('HFF-hTERT' in key[0] and 'HFF-hTERT' in key[1]):
                samecelltypesameexp[0].extend([score])
        else:
            isdil.extend([score])
            if '_flav' in key[0] or '_flav' in key[1] or '_tnfa' in key[0] or '_tnfa' in key[1]: continue
            if ('IMR90' in key[0] and 'IMR90' in key[1]) or ('hESC' in key[0] and 'hESC' in key[1]) or ('HFF-hTERT' in key[0] and 'HFF-hTERT' in key[1]):
                samecelltypediffexp.extend([score])

    diffcelltype = []
    for key,score in alltadsim.iteritems():
        ct1split = key[0].split('_')
        ct2split = key[1].split('_')
        if ct1split[0] != ct2split[0]:
            diffcelltype.extend([score])

    tad.plotISDilPairs(ispairs,dilpairs,isdil,measure,outputloc+'tadsim_protocolpairboxplot.pdf')
    tad.plotSameCellTypeBoxplots(samecelltypediffexp, samecelltypesameexp, diffcelltype, measure, outputloc+'tadsim_samecelltype_protocol_comparison.pdf')

    dilreps = []
    isreps = []
    for key,score in reptadsim.iteritems():
        if key in dil or key[:-2] in dil:
            dilreps.extend([score])
        elif key in insitu or key[:-2] in insitu:
            isreps.extend([score])
    print 'total number of dilution replicate pairs =', len(dilreps)
    print 'total number of in situ replicate pairs =', len(isreps)
    tad.protocolRepBoxplots(isreps,dilreps,measure,outputloc+'tadsim_protocolrepboxplot.pdf')

def tissueAnalysis(alltadsim,tadsim_donor,reptadsim,measure,outputloc):
    
    # where do donor values fit in here??
    withintissue = tadsim_donor.values()
    acrosstissue = []
    background = []
    replicates = []
    tissuesamps = ['SkelMuscle','TransColon','Adrenal','Bladder','DPC','Hippocampus','Lung','Ovary','Pancreas','Psoas','Righ\
tVentricle','SmallBowel','Spleen']
    for key,score in alltadsim.iteritems():
        if key[0] in tissuesamps and key[1] in tissuesamps:
            acrosstissue.extend([score])
        else:
            background.extend([score])
    tad.plotTissueBoxplots(withintissue, acrosstissue, reptadsim.values(), background, measure, outputloc+'tadsim_tissueboxplots.pdf')

def trioAnalysis(alltadsim,reptadsim,measure,outputloc):

    withintrio = []
    acrosstrio = []
    background = []
    bloodlymph = []
    parentchild = []
    parentparent = []
    children = ['HG00733', 'HG00514', 'GM19240']
    bloodlymphtypes = ['GM12878','GM20431','GM19240','GM19239','GM19238','HG00733','HG00732','HG00731','HG00514','HG00513','\
HG00512']
    for key,score in alltadsim.iteritems():
        if key[0][:5] == key[1][:5] and (key[0][:4] == 'HG00' or key[0][:4] == 'GM19'):
            #print 'within trio,',key
            withintrio.extend([score])
            if key[0] in children or key[1] in children:
                parentchild.extend([score])
                #print 'parentchild',key,score
            elif key[0] != key[1]:
                parentparent.extend([score])
                #print 'parentparent',key,score
        elif key[0][:5] != key[1][:5] and (key[0][:4] == 'HG00' or key[0][:4] == 'GM19') and (key[1][:4] == 'HG00' or key[1]\
[:4] == 'GM19'):
            acrosstrio.extend([score])
            bloodlymph.extend([score])
            #print 'across trio,', key
        elif key[0] in bloodlymphtypes and key[1] in bloodlymphtypes:
            #print key
            bloodlymph.extend([score])
        else:
            background.extend([score])
    trioreps = []
    for key,score in reptadsim.iteritems():
        if key[0][:4] == 'GM19' or key[0][:4] == 'HG00':
            trioreps.extend([score])

    tad.plotTrioBoxplots(parentchild,parentparent,trioreps,bloodlymph,background,measure,outputloc+'tadsim_trioboxplots_byfamily.pdf')

def main(repfileseed, bgfileseed, donorfileseed, res, rawcountfilename, celltypelist, outputloc):

    measure = 'TADsim'

    bgdata,repdata,donordata,celltypelist = readAllFiles(bgfileseed,repfileseed,donorfileseed,celltypelist)
    print 'done reading all data files'

    chrlengths = [248956422, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566] # hg19 chr lengths from https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
    for i,chrlen in enumerate(chrlengths):
        chrlengths[i] = int(np.ceil(float(chrlen)/res))

    alltadsim = hicutil.calcPercSimilarity(bgdata, chrlengths)
    reptadsim = hicutil.calcPercSimilarity(repdata, chrlengths)
    donortadsim = hicutil.calcPercSimilarity(donordata, chrlengths)

    tad.plotReplicateVNonReplicate(reptadsim.values(), alltadsim.values(), [], measure, outputloc+'tadsim_replicateVnonreplicate.pdf')
    if len(rawcountfilename) > 0:
        rawcontactcounts = tad.readRawCountsFile(rawcountfilename)
        tad.plotSimVContactCounts(rawcontactcounts, reptadsim, measure, outputloc+'tadsim_rawcountsVrep.pdf')

    simmat,labels = tad.generateHeatMap(alltadsim, measure, outputloc+'tadsim_fullheatmap.png')

    resFragAnalysis(alltadsim,measure,outputloc)
    acrossLabAnalysis(alltadsim,reptadsim,measure,outputloc)
    insituDilutionAnalysis(alltadsim,reptadsim,measure,outputloc)
    tissueAnalysis(alltadsim,donortadsim,reptadsim,measure,outputloc)
    trioAnalysis(alltadsim,reptadsim,measure,outputloc)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-rf', type=str, help='fileseed for comparisons of replicate TAD sets')
    parser.add_argument('-nrf', type=str, nargs='+', help='fileseed for comparisons of non-replicate TAD sets')
    parser.add_argument('-td',type=str, help='fileseed for comparisons of TAD sets from tissue samples of different donors')
    parser.add_argument('-rcf', type=str, default='', help='Raw count file')
    parser.add_argument('-res', type=int, help='Resolution of Hi-C data')
    parser.add_argument('-o', type=str, help='Location for output files')
    parser.add_argument('-c', type=str, help='file containing list of cell types')
    args = parser.parse_args()

    main(args.rf, args.nrf, args.td, args.res, args.rcf, args.c, args.o)
