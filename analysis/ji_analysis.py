import csv
import numpy as np
import scipy.stats
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import glob
import sys
from sklearn import manifold, cluster, metrics
import umap
from scipy.spatial.distance import pdist, squareform
from sklearn import datasets
from scipy.cluster.hierarchy import linkage
import tadanalysisutil as tad

newpalette = ['#8dd3c7','#58508d','#bc5090','#003f5c','#ff6361']
sns.set(style='whitegrid',palette=sns.color_palette(newpalette),color_codes=False,font='Ubuntu')

def readFilenames(filename):
    filelist = []
    with open(filename,'rb') as f:
        freader = csv.reader(f)
        for line in freader:
            filelist.extend(line)
    return filelist

def readTADfile(filename, res):
    # written for Armatus files
    tadbdys = []
    with open(filename, 'rb') as tf:
        freader = csv.reader(tf, delimiter = '\t')
        for line in freader:
            tadbdys.append([int(line[1]), int(line[2])])
    if len(tadbdys) == 0: return np.array(tadbdys)
    if tadbdys[0][0] > tadbdys[-1][1]:
        tadbdys = np.flipud(np.array(tadbdys))
    tadbdys = np.array(tadbdys)
    tadbdys[:,0] = np.floor(np.divide(tadbdys[:,0],res))
    tadbdys[:,1] = np.floor(np.divide((tadbdys[:,1]+1),res)-1)
    return tadbdys

def readTADs(filelist,res):
    # this function is dependent on my file structure, where all TAD files are in an 'armatusresults' directory, by folder named by cell type. If your directory structure is different, please adapt the first few lines to accommodate your naming system.
    tads = {}
    for filename in filelist:
        #need to extract celltype and chrnumber from filename
        fsplit = filename.split('/')
        idx = fsplit.index('armatusresults')
        celltype = fsplit[idx+1]
        if celltype == 'SkelMuscle' or celltype == 'TransColon':
            celltype = filename.split('_')[1]
        chrloc = filename.find('Chr')
        if chrloc == -1:
            chrloc = filename.find('chr')
        try:
            chrnum = int(filename[chrloc+3:chrloc+5])
        except:
            chrnum = int(filename[chrloc+3])
        if (celltype,chrnum) in tads:
            tads[(celltype,chrnum)].append(readTADfile(filename,res))
        else:
            tads[(celltype,chrnum)] = [readTADfile(filename,res)]
    return tads

def calcJI(tadset1,tadset2):
    # have to return number of overlaps and total number by chromosome, then later add them up for genome-wide value 
    tadbdys1 = tadset1.flatten()
    tadbdys2 = tadset2.flatten()
    intersectsize = len(set(tadbdys1) & set(tadbdys2))
    unionsize = len(set(tadbdys2) | set(tadbdys2))
    return intersectsize,unionsize

def plot2Dembed(simvals,labels,filename):

    distvals = 1-simvals
    umapvar = umap.UMAP(metric='precomputed')
    coords = umapvar.fit_transform(distvals)

    listforfile = [[x, 0] for x in labels]
    listforfile = [['cell type', 'K1']]+listforfile

    print 'computing spectral clusters'
    for nclus in range(2,11):
        specclus = cluster.SpectralClustering(n_clusters=nclus,eigen_solver='arpack',affinity='precomputed')
        clusters = specclus.fit(simvals)
        cluslabels = clusters.labels_

        listforfile[0].extend(['K'+str(nclus)])
        for i,sample in enumerate(listforfile[1:]):
            sample.extend([cluslabels[i]])

        silscore = metrics.silhouette_score(distvals, cluslabels, metric='precomputed')
        print '# of clusters =',nclus,', silhouette score =',silscore

        fig,ax1 = plt.subplots(figsize=(10,10))
        col = ['b', 'g', 'r', 'm', 'c', 'k','y']
        for clus in range(nclus):
            #print coords[cluslabels==clus,0], coords[cluslabels==clus,1]
            plt.scatter(coords[cluslabels == clus,0],coords[cluslabels == clus,1])

        for label,x,y in zip(labels,coords[:,0],coords[:,1]):
            plt.annotate(label, xy=(x,y))

        plt.tight_layout()
        newfilename = filename[:-4]+'_nclus'+str(nclus)+'.pdf'
        plt.savefig(newfilename,bbox_inches='tight')
        plt.close(fig)
        print '2D embedding figure saved as', newfilename

def readPercSimFile(filename):
    
    with open(filename,'rb') as f:
        freader = csv.reader(f, delimiter='\t')
        percsimvals = {}
        for line in freader:
            percsimvals[(line[0], line[1])] = float(line[2])
    return percsimvals

def computeAllRepJI(repfilelist,res):

    reptads = readTADs(repfilelist,res)
    repji = {}
    for key,tadsets in reptads.iteritems():
        numreps = len(tadsets)
        for i in xrange(numreps-1):
            for j in xrange(i+1,numreps):
                jinum, jidenom = calcJI(tadsets[i], tadsets[j])
                if (key[0],i,j) in repji:
                    repji[(key[0],i,j)].append([jinum, jidenom])
                else:
                    repji[(key[0],i,j)] = [[jinum,jidenom]]
    repjivals = []
    for key,jinums in repji.iteritems():
        jinums = np.array(jinums)
        repjivals.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
    return repji, repjivals

def computeAllNonRepJI(nonrepfilelist,res):

    nonreptads = readTADs(nonrepfilelist, res)
    nonrepji = {}
    pairsdone = []
    for key1,tadset1 in nonreptads.iteritems():
        for key2,tadset2 in nonreptads.iteritems():
            if key1[1] != key2[1] or key1[0] == key2[0]: continue
            if [key1[0], key2[0], key1[1]] in pairsdone or [key2[0], key1[0], key1[1]] in pairsdone: continue
            jinum, jidenom = calcJI(tadset1[0], tadset2[0])
            if (key1[0], key2[0]) in nonrepji:
                nonrepji[(key1[0],key2[0])].append([jinum, jidenom])
            elif (key2[0], key1[0]) in nonrepji:
                nonrepji[(key2[0],key1[0])].append([jinum,jidenom])
            else:
                nonrepji[(key1[0],key2[0])] = [[jinum,jidenom]]
            pairsdone.append([key1[0],key2[0],key1[1]])
    nonrepjivals = []
    for key,jinums in nonrepji.iteritems():
        jinums = np.array(jinums)
        nonrepjivals.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
    return nonrepji, nonrepjivals

def computeAllDonorJI(donorfilelist,res):

    tissuedonortads = readTADs(donorfilelist, res)
    donorji = {}
    pairsdone = []
    for key1,tadset1 in tissuedonortads.iteritems():
        if len(tadset1) == 2:
            jinum, jidenom = calcJI(tadset1[0], tadset1[1])
            if key1[0] in donorji:
                donorji[key1[0]].append([jinum,jidenom])
            else:
                donorji[key1[0]] = [[jinum,jidenom]]
        else:
            for key2,tadset2 in tissuedonortads.iteritems():
                if key1[1] != key2[1] or key1[0] == key2[0]: continue
                if key1[0][:6] == key2[0][:6] and [key1[0],key2[0]] not in pairsdone and [key2[0], key1[0]] not in pairsdone:
                    jinum,jidenom = calcJI(tadset1[0], tadset2[0])
                    if (key1[0],key2[0]) in donorji:
                        donorji[(key1[0],key2[0])].append([jinum,jidenom])
                    else:
                        donorji[(key1[0],key2[0])] = [[jinum,jidenom]]
                    pairsdone.append([key1[0],key2[0]])
    donorjivals = []
    for key,jinums in donorji.iteritems():
        jinums = np.array(jinums)
        donorjivals.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
    return donorjivals

def resFragAnalysis(nonrepji, measure, outputloc):

    resfragdata = [['hESC_DpnII_rep', 'hESC_HindIII', 'hESC_NcoI'], ['HFF-hTERT_DpnII', 'HFF-hTERT_HindIII_beads', 'HFF-hTERT_HindIII_plate', 'HFF-hTERT_MboI', 'HFF-hTERT_NcoI']]
    resfragjis = [[],[]]
    nonrepjivals = []
    for key, jinums in nonrepji.iteritems():
        jinums = np.array(jinums)
        if (key[0] in resfragdata[0] and key[1] in resfragdata[0]):
            resfragjis[0].extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
        elif (key[0] in resfragdata[1] and key[1] in resfragdata[1]):
            resfragjis[1].extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
        else:
            nonrepjivals.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
    
    tad.plotResFragResults(resfragjis, nonrepjivals, measure, outputloc+'ji_resfragboxplot_withnew_violin.pdf')

def acrossLabAnalysis(nonrepji,repjivals,measure,outputloc):

    acrosslabdata = [['hESC_Dixon2015','hESC_FP','hESC_HindIII','hESC_Jin','hESC'],['IMR90dixon', 'IMR90rao', 'IMR90_normal', 'IMR90_RI']]
    acrosslabjis = [[],[]]
    nonrepjivals = []
    for key,jinums in nonrepji.iteritems():
        jinums = np.array(jinums)
        if (key[0] in acrosslabdata[0] and key[1] in acrosslabdata[0]):
            acrosslabjis[0].extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
        elif (key[0] in acrosslabdata[1] and key[1] in acrosslabdata[1]):
            acrosslabjis[1].extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
        else:
            nonrepjivals.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
    
    tad.plotAcrossLabResults(acrosslabjis, nonrepjivals, repjivals, measure, outputloc+'ji_acrosslabboxplot_withnew_violin.pdf')

def insituDilutionAnalysis(nonrepji, repjidict, measure, outputloc):
    
    dil = ['A549','Caki2','G401','LNCaP-FGC','NCI-H460','Panc1','RPMI-7951','SKMEL5','SKNDZ','SKNMC','T47D','IMR90dixon','hESC','hESC_HindIII','HFF-hTERT_HindIII_techreps','HG00733','HG00732','HG00731','HG00514','HG00513','HG00512','GM19238','GM19239','GM19240','hESC_Jin','IMR90_flav','IMR90_normal','IMR90_tnfa','hESC_Dixon2015','GM20431','BrainMicroEndo','AstrocyteCerebellum','AstrocyteSpinalCord','DLD1','BrainPericyte','EndomMicroEndoth','HepSin','ACHN','IMR90_RI','hESC_FP','Adrenal','Bladder','DPC','Hippocampus','Lung','Ovary','Pancreas','Psoas','RightVentricle','SmallBowel','Spleen']
    insitu = ['HFF-hTERT_HindIII_beads','HFF-hTERT_HindIII_plate', 'IMR90rao','GM12878','HMEC','HUVEC','K562rao','KBM7','NHEK','HFF-hTERT_MboI','SkelMuscle','TransColon','hESC_NcoI','HFF-hTERT_NcoI','hESC_DpnII', 'hESC_DpnII_rep','HFF-hTERT_DpnII','HFF-c6']    

    ispairs = []
    dilpairs = []
    isdil = []
    samecelltypediffexp = []
    samecelltypesameexp = [[],[]] # 1 list for in situ, 1 for dilution
    for key,jinums in nonrepji.iteritems():
        jinums = np.array(jinums)
        if (key[0] in dil and key[1] in dil):
            dilpairs.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
            if '_flav' in key[0] or '_flav' in key[1] or '_tnfa' in key[0] or '_tnfa' in key[1]: continue
            if ('IMR90' in key[0] and 'IMR90' in key[1]) or ('hESC' in key[0] and 'hESC' in key[1]) or ('HFF-hTERT' in key[0] and 'HFF-hTERT' in key[1]):
                samecelltypesameexp[1].extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
        elif key[0] in insitu and key[1] in insitu:
            ispairs.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
            if '_flav' in key[0] or '_flav' in key[1] or '_tnfa' in key[0] or '_tnfa' in key[1]: continue
            if ('IMR90' in key[0] and 'IMR90' in key[1]) or ('hESC' in key[0] and 'hESC' in key[1]) or ('HFF-hTERT' in key[0] and 'HFF-hTERT' in key[1]):
                samecelltypesameexp[0].extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
        else:
            isdil.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
            if '_flav' in key[0] or '_flav' in key[1] or '_tnfa' in key[0] or '_tnfa' in key[1]: continue
            if ('IMR90' in key[0] and 'IMR90' in key[1]) or ('hESC' in key[0] and 'hESC' in key[1]) or ('HFF-hTERT' in key[0] and 'HFF-hTERT' in key[1]):
                samecelltypediffexp.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])

    diffcelltype = []
    for key,jinums in nonrepji.iteritems():
        jinums = np.array(jinums)
        ct1split = key[0].split('_')
        ct2split = key[1].split('_')
        if ct1split[0] != ct2split[0]:
            diffcelltype.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])

    tad.plotISDilPairs(ispairs,dilpairs,isdil,measure,outputloc+'ji_protocolpairboxplot_withnew_violin.pdf')
    tad.plotSameCellTypeBoxplots(samecelltypediffexp, samecelltypesameexp, diffcelltype, measure, outputloc+'ji_samecelltype_protocol_comparison_withnew_violin.pdf')

    dilreps = []
    isreps = []
    for key,jinums in repjidict.iteritems():
        jinums = np.array(jinums)
        if key[0] in dil:
            dilreps.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
        elif key[0] in insitu:
            isreps.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
    
    tad.protocolRepBoxplots(isreps,dilreps,measure,outputloc+'ji_protocolrepboxplot_withnew_violin.pdf')

def tissueAnalysis(nonrepji,donorjivals,repjivals,measure,outputloc):

    tissuesamps = ['SkelMuscle','TransColon','Adrenal','Bladder','DPC','Hippocampus','Lung','Ovary','Pancreas','Psoas','RightVentricle','SmallBowel','Spleen']
    acrosstissue = []
    background = []
    for key,jinums in nonrepji.iteritems():
        jinums = np.array(jinums)
        if key[0] in tissuesamps and key[1] in tissuesamps:
            acrosstissue.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
        else:
            background.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
    print 'average JI across tissue types =', np.mean(acrosstissue)
    tad.plotTissueBoxplots(donorjivals,acrosstissue,repjivals,background,measure,outputloc+'ji_tissueboxplots_withnew_violin.pdf')

def trioAnalysis(nonrepji,repjidict,measure,outputloc):

    withintrio = []
    acrosstrio = []
    background = []
    bloodlymph = []
    parentchild = []
    parentparent = []
    children = ['HG00733', 'HG00514', 'GM19240']
    bloodlymphtypes = ['GM12878','GM20431','GM19240','GM19239','GM19238','HG00733','HG00732','HG00731','HG00514','HG00513','HG00512']
    for key,jinums in nonrepji.iteritems():
        jinums = np.array(jinums)
        if key[0][:5] == key[1][:5] and (key[0][:4] == 'HG00' or key[0][:4] == 'GM19'):
            #print 'within trio,',key
            withintrio.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
            if key[0] in children or key[1] in children:
                parentchild.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
                #print 'parentchild',key,np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))
            else:
                parentparent.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
                #print 'parentparent',key,np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))
        elif key[0][:5] != key[1][:5] and (key[0][:4] == 'HG00' or key[0][:4] == 'GM19') and (key[1][:4] == 'HG00' or key[1][:4] == 'GM19'):
            acrosstrio.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
            bloodlymph.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
            #print 'across trio,', key
        elif key[0] in bloodlymphtypes and key[1] in bloodlymphtypes:
            #print key
            bloodlymph.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
        else:
            background.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
    trioreps = []
    for key,jinums in repjidict.iteritems():
        jinums = np.array(jinums)
        if key[0][:4] == 'GM19' or key[0][:4] == 'HG00':
            trioreps.extend([np.sum(jinums[:,0])/float(np.sum(jinums[:,1]))])
    tad.plotTrioBoxplots(parentchild,parentparent,trioreps,bloodlymph,background,measure,outputloc+'ji_trioboxplots_byfamily_withnew_violin.pdf')

def main(files_rep, files_nonrep, files_tissuedonor, rawcountsfile, res, outputloc):

    measure = 'Jaccard Index'

    repfilelist = readFilenames(files_rep)
    nonrepfilelist = readFilenames(files_nonrep)

    repjidict,repjivals = computeAllRepJI(repfilelist,res)
    nonrepji,nonrepjivals = computeAllNonRepJI(nonrepfilelist,res)

    #reppairs = []
    #with open('go/testing/allhicrep_replicates.txt','rb') as f:
    #    freader = csv.reader(f,delimiter=' ')
    #    for line in freader:
    #        reppairs.extend([(line[0][:-5], int(line[0][-1])-1, int(line[1][-1])-1)])
    #for pair in reppairs:
    #    if pair not in repjidict and (pair[0],pair[2],pair[1]) not in repjidict:
    #        print 'in hicrep, not ji:', pair
    #for key in repjidict.keys():
    #    if key not in reppairs and (key[0],key[2],key[1]) not in reppairs:
    #        print 'in ji, not hicrep:',key
    #sys.exit()

    tad.plotReplicateVNonReplicate(repjivals, nonrepjivals, [], measure, outputloc+'ji_replicateVnonreplicate.pdf')
    if len(rawcountsfile) > 0:
        rawcontactcounts = tad.readRawCountsFile(rawcountsfile)
        tad.plotSimVContactCounts(rawcontactcounts, repjidict, measure, outputloc+'ji_rawcountsVrep.pdf')

    resFragAnalysis(nonrepji, measure, outputloc)
    acrossLabAnalysis(nonrepji,repjivals,measure,outputloc)
    simmat,labels = tad.generateHeatMap(nonrepji, measure, outputloc+'ji_fullheatmap.png')    
    #plot2Dembed(simmat,labels,outputloc+'ji_2Dembedding.pdf')
    insituDilutionAnalysis(nonrepji, repjidict, measure, outputloc)
    trioAnalysis(nonrepji,repjidict,measure,outputloc)

    donorfilelist = readFilenames(files_tissuedonor)
    donorjivals = computeAllDonorJI(donorfilelist,res)
    tissueAnalysis(nonrepji,donorjivals,repjivals,measure,outputloc)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-rf', type=str, help='File containing filenames of TAD sets from replicates')
    parser.add_argument('-nrf', type=str, help='File containing filenames of TAD sets from non-replicates')
    parser.add_argument('-td',type=str, help='File containing filenames of TAD sets from different donors of the same tissue')
    parser.add_argument('-rcf', type=str, default='', help='Raw count file')
    parser.add_argument('-res', type=int, help='Resolution of Hi-C data')
    parser.add_argument('-o', type=str, help='Location for output files')

    args = parser.parse_args()
    main(args.rf, args.nrf, args.td, args.rcf, args.res, args.o)
