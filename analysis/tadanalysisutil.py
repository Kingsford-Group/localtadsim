import csv
import numpy as np
import scipy.stats
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
import sys
from sklearn import manifold, cluster, metrics
import umap
from scipy.spatial.distance import pdist, squareform
from sklearn import datasets
from scipy.cluster.hierarchy import linkage

newpalette = ['#8dd3c7','#58508d','#bc5090','#003f5c','#ff6361']
#newpalette = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3']
sns.set(style='whitegrid',palette=sns.color_palette(newpalette),color_codes=False,font='Ubuntu')

def plotReplicateVNonReplicate(repvals, nonrepvals, intralabvals, measure, filename):

    newpalette = ['#8dd3c7','#58508d','#bc5090','#003f5c','#ff6361']
    sns.set(style='whitegrid',palette=sns.color_palette(newpalette),color_codes=False,font='Ubuntu')

    if len(intralabvals) == 0:
        alldata = [repvals,nonrepvals]
    else:
        alldata = [repvals,nonrepvals,intralabvals]
    mwuval,pval = scipy.stats.mannwhitneyu(repvals, nonrepvals,alternative='greater')
    print 'pvalue of replicates v nonreplicates =', pval
    print 'average replicate similarity =', np.mean(repvals)
    print 'average non-replicate similarity =', np.mean(nonrepvals)
    #print alldata                                                           

    fig = plt.figure()
    ax = plt.gca()
    #plt.boxplot(alldata,widths=0.5,medianprops=dict(color='red'))             

    alldatavals = nonrepvals+repvals
    alllabels =  ['Non-Replicate']*len(nonrepvals) + ['Replicate']*len(repvals)
    datadict = {'labels':alllabels, 'values':alldatavals}
    alldatadf = pd.DataFrame(datadict)
    sns.violinplot(x=alldatadf['labels'],y=alldatadf['values'],cut=1)
    plt.tick_params(top=False,right=False,left=False)
    plt.xticks([0,1],['Non-Replicates','Replicates'],size=12)
    ax.set_aspect(3)
    ax.set_ylim(0,1)
    plt.ylabel(measure)
    plt.xlabel('')
    plt.savefig(filename,bbox_inches='tight')
    plt.close(fig)
    print 'Saved replicate comparison figure to',filename


def plotResFragResults(pointvals, backgroundvals, measure, filename):

    newpalette = ['#8dd3c7','#58508d','#bc5090','#003f5c','#ff6361']
    sns.set(style='whitegrid',palette=sns.color_palette(newpalette),color_codes=False,font='Ubuntu')

    fig = plt.figure()
    ax = plt.gca()
    datadict = {'xloc':[0]*len(backgroundvals), 'values':backgroundvals}
    datadf = pd.DataFrame(datadict)
    #plt.boxplot([backgroundvals,[]], widths=0.5,medianprops=dict(color='red'))#,showfliers=False)                          

    #order points by yvals
    labels = ['hESC']*len(pointvals[0])+['HFF-hTERT']*len(pointvals[1])
    allvals =  pointvals[0]+pointvals[1]
    labelledvals = [list(x) for x in zip(labels,allvals)]
    labelledvals.sort(key=lambda x: x[1])
    labels = [l[0] for l in labelledvals]
    allvals = [l[1] for l in labelledvals]

    if labels[0] != 'hESC':
        ptpalette = newpalette[3:5]
    else:
        ptpalette = [newpalette[4],newpalette[3]]

    xvals = np.linspace(1-.03*len(allvals),1+.03*len(allvals),len(allvals))
    resfragdata = {'labels': labels, 'xvals':xvals, 'yvals':allvals}
    sns.scatterplot(x='xvals',y='yvals',hue='labels',data=resfragdata,palette=ptpalette,legend=False,s=35)
    sns.violinplot(x=datadf['xloc'],y=datadf['values'],cut=1)
    plt.xticks([0,1],['Background','Same cell \n type, different \n restriction enzyme'],size=12)
    ax.set_aspect(3)
    ax.set_ylim(0,1)
    ax.set_xlim(-0.5,1.5)
    plt.ylabel(measure)
    plt.xlabel('')
    plt.tick_params(top=False,right=False,left=False)
    #ax.legend(handles,['hESC','HFF-hTERT'],loc='lower right',frameon=False)
    #plt.title(plottitle)                                                      

    plt.savefig(filename,bbox_inches='tight')
    plt.close(fig)
    print 'Saved restriction fragment figure to', filename

def plotAcrossLabResults(pointvals, backgroundvals, repvals, measure, filename):

    newpalette = ['#8dd3c7','#58508d','#bc5090','#003f5c','#ff6361']
    sns.set(style='whitegrid',palette=sns.color_palette(newpalette),color_codes=False,font='Ubuntu')

    #order points by yvals
    labels = ['hESC']*len(pointvals[0])+['IMR90']*len(pointvals[1])
    allvals =  pointvals[0]+pointvals[1]
    labelledvals = [list(x) for x in zip(labels,allvals)]
    labelledvals.sort(key=lambda x: x[1])
    labels = [l[0] for l in labelledvals]
    allvals = [l[1] for l in labelledvals]

    fig = plt.figure()
    ax = plt.gca()
    xvals = np.linspace(2-.03*len(allvals),2+.03*len(allvals),len(allvals))
    crosslabdata = {'labels': labels,'xvals': xvals,'yvals': allvals}
    if labels[0] == 'hESC':
        ptpalette = newpalette[3:5]
    else:
        ptpalette = [newpalette[4],newpalette[3]]
    sns.scatterplot(x='xvals',y='yvals',hue='labels',data=crosslabdata,palette=ptpalette,legend=False,s=35)
    #alldatavals = nonrepvals+repvals
    datadict = {'labels':['Non-Replicates']*len(backgroundvals) + ['Replicate']*len(repvals), 'values':backgroundvals+repvals}
    alldatadf = pd.DataFrame(datadict)
    sns.violinplot(x=alldatadf['labels'],y=alldatadf['values'],cut=1)
    ax.set_ylim(0,1)
    ax.set_xlim(-0.5,2.5)
    plt.ylabel(measure)
    plt.xticks([0,1,2],['All \n Non-replicates','All \n Replicates','Same cell type, \n different study'],size=12)
    ax.set_aspect(4)
    plt.xlabel('')
    plt.tick_params(top=False,right=False,left=False)
    #plt.title(plottitle)                                                     
    #ax.legend(handles,['hESC','IMR90'],loc='lower right',frameon=False)        
    plt.savefig(filename,bbox_inches='tight')
    plt.close(fig)
    print 'Saved across lab figure to', filename

def generateHeatMap(nonrepdict, measure, filename):

    labelledvals = []
    for key,data in nonrepdict.iteritems():
        if isinstance(data,list):
            data = np.array(data)
            jival = np.sum(data[:,0])/float(np.sum(data[:,1]))
            labelledvals.append([key,jival])
        else:
            labelledvals.append([key,data])

    # need to order matrix rows by dataset/lab of origin                     
    raodata = ['GM12878', 'IMR90rao', 'NHEK', 'HMEC', 'HUVEC','K562rao','KBM7']
    encodedata = ['A549','Caki2','G401','LNCaP-FGC','NCI-H460','Panc1','RPMI-7951','SKMEL5','SKNDZ','SKNMC','T47D','GM20431', 'SkelMuscle', 'TransColon', 'BrainMicroEndo','AstrocyteCerebellum','AstrocyteSpinalCord','DLD1','BrainPericyte','EndomMicroEndoth','HepSin', 'ACHN']
    dixon2012data = ['IMR90dixon','hESC']
    dixon2015data = ['hESC_Dixon2015']
    fourdndata = ['hESC_NcoI','hESC_DpnII_rep','hESC_HindIII','hESC_DpnII','HFF-hTERT_HindIII_techreps','HFF-hTERT_DpnII','HFF-hTERT_NcoI','HFF-hTERT_HindIII_beads','HFF-hTERT_HindIII_plate','HFF-hTERT_MboI','HFF-c6','HG00733','HG00732', 'HG00731', 'HG00514', 'HG00513', 'HG00512', 'GM19238', 'GM19239', 'GM19240']
    jindata = ['hESC_Jin', 'IMR90_flav', 'IMR90_normal', 'IMR90_tnfa']
    ridata = ['IMR90_RI']
    fpdata = ['hESC_FP']
    schmittdata = ['Adrenal','Bladder','DPC','Hippocampus','Lung','Ovary','Pancreas','Psoas','RightVentricle','SmallBowel','Spleen']

    linelocs = [len(raodata), len(encodedata), len(dixon2012data), len(dixon2015data), len(fourdndata), len(jindata), len(ridata), len(fpdata),len(schmittdata)]

    celltypelist = raodata+encodedata+dixon2012data+dixon2015data+fourdndata+jindata+ridata+fpdata+schmittdata

    #put percsim data into matrix                                        
    simmat = np.ones((len(celltypelist),len(celltypelist)))
    totalsim = sorted(labelledvals, key=lambda x:x[1], reverse=True)

    for idx,ctpair in enumerate(totalsim):
        listidx1 = celltypelist.index(ctpair[0][0])
        listidx2 = celltypelist.index(ctpair[0][1])
 #       print ctpair[0]                                 
        simmat[listidx1, listidx2] = ctpair[1]
        simmat[listidx2, listidx1] = ctpair[1]

    plotHeatMap(simmat,celltypelist,linelocs,measure,filename)

    # order by clustering                                    
    clusorder = compute_matrix_clus(1-simmat,method="single")
    #print clusorder                          
    #sys.exit()                                
    celltypelist2 = []
    simmat2 = np.ones((len(celltypelist),len(celltypelist)))
    a,b = np.triu_indices(len(celltypelist),k=1)
    simmat2[a,b] = simmat[ [ clusorder[i] for i in a], [clusorder[j] for j in b] ]
    simmat2[b,a] = simmat2[a,b]
    for oldidx in clusorder:
        celltypelist2.extend([celltypelist[oldidx]])
    #plotHeatMap(simmat2, celltypelist2,[],filename)

    return simmat,celltypelist


def plotHeatMap(data,labels,linelocs,measure,filename):

    newpalette = ['#8dd3c7','#58508d','#bc5090','#003f5c','#ff6361']
    sns.set(style='whitegrid',palette=sns.color_palette(newpalette),color_codes=False,font='Ubuntu')

    hescnum = 1
    imr90num = 1
    hffhtertnum = 1
    for i,l in enumerate(labels):
        if 'hESC' in l:
            labels[i] = 'hESC ('+str(hescnum)+ ')'
            hescnum += 1
        elif 'IMR90' in l:
            labels[i] = 'IMR90 ('+str(imr90num)+ ')'
            imr90num += 1
        elif 'HFF-hTERT' in l:
            labels[i] = 'HFF-hTERT ('+str(hffhtertnum)+ ')'
            hffhtertnum += 1
        elif l == 'K562rao':
            labels[i] = 'K562'

    fig,ax1 = plt.subplots(figsize=(15,15))
    #sns.set(font_scale=0.9)
    sns.heatmap(data,cmap='YlGnBu',cbar_kws = {'label':measure,'shrink':0.5},xticklabels=labels,yticklabels=labels,square=True)

    ax1.grid(False)
    plt.tick_params(top=False,right=False,left=False,bottom=False)
    linestart = 0
    for num in linelocs:
        #plt.axhline(y=num+linestart,color='black')    
        #plt.axvline(x=num+linestart,color='black') 
        ax1.plot([linestart,linestart+num],[linestart,linestart],color='black',linestyle='dashed',linewidth=3)
        ax1.plot([linestart+num,linestart+num],[linestart,linestart+num],color='black',linestyle='dashed',linewidth=3)
        linestart+=num
    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    plt.close(fig)
    print 'Heatmap saved as',filename

def seriation(Z,N,cur_index):
    '''                                                   
        input:     
            - Z is a hierarchical tree (dendrogram)
            - N is the number of points given to the clustering process
            - cur_index is the position in the tree for the recursive traversal
        output:
            - order implied by the hierarchical tree Z

        seriation computes the order implied by a hierarchical tree (dendrogram)
    '''
    if cur_index < N:
        return [cur_index]
    else:
        left = int(Z[cur_index-N,0])
        right = int(Z[cur_index-N,1])
        return (seriation(Z,N,left) + seriation(Z,N,right))

def compute_matrix_clus(dist_mat,method="ward"):
    '''                                            
        input:                                                       
            - dist_mat is a distance matrix                      
            - method = ["ward","single","average","complete"]         
        output:                                                  
            - seriated_dist is the input dist_mat,                     
              but with re-ordered rows and columns         
              according to the seriation, i.e. the           
              order implied by the hierarchical tree         
            - res_order is the order implied by             
              the hierarhical tree                          
            - res_linkage is the hierarhical tree (dendrogram)       
                                                                     
        compute_serial_matrix transforms a distance matrix into          
        a sorted distance matrix according to the order implied   
        by the hierarchical tree (dendrogram)                         
    '''
    N = len(dist_mat)
    flat_dist_mat = squareform(dist_mat)
    res_linkage = linkage(flat_dist_mat, method=method,optimal_ordering=True)
    res_order = seriation(res_linkage, N, N + N-2)
    #seriated_dist = np.zeros((N,N))  
    #a,b = np.triu_indices(N,k=1)          
    #seriated_dist[a,b] = dist_mat[ [res_order[i] for i in a], [res_order[j] for j in b]]         
    #seriated_dist[b,a] = seriated_dist[a,b]

    #return seriated_dist, res_order, res_linkage       
    return res_order

def readRawCountsFile(filename):
    with open(filename,'rb') as f:
        freader = csv.reader(f, delimiter='\t')
        rawcontactcounts = {}
        freader.next()
        for line in freader:
            rawcontactcounts[line[0]] = [float(x) for x in line[1:]]
    return rawcontactcounts

def plotSimVContactCounts(rawcontactcounts, repdata, measure, figname):

    newpalette = ['#8dd3c7','#58508d','#bc5090','#003f5c','#ff6361']
    sns.set(style='whitegrid',palette=sns.color_palette(newpalette),color_codes=False,font='Ubuntu')

    rawcountsx = []
    simvalsy = []
    labels = []

    for celltype,rawcounts in rawcontactcounts.iteritems():
        if celltype in ['HepG2','SJCRH30','HeLaS3']:
            continue
        numreps = len(rawcounts)
        for i in xrange(numreps-1):
            for j in xrange(i+1,numreps):
                if measure == 'Jaccard Index':
                    if (celltype,i,j) in repdata:
                        rawcountsx.extend([min([rawcounts[i],rawcounts[j]])])
                        jivals = np.array(repdata[(celltype,i,j)])
                        simvalsy.extend([np.sum(jivals[:,0])/float(np.sum(jivals[:,1]))])
                        labels.extend(celltype)
                    else:
                        print 'cant find', (celltype,i,j)
                elif measure == 'HiCRep':
                    key1 = (celltype+'_rep'+str(i+1), celltype+'_rep'+str(j+1))
                    key2 = (celltype+'_rep'+str(j+1), celltype+'_rep'+str(i+1))
                    if key1 in repdata:
                        rawcountsx.extend([min([rawcounts[i],rawcounts[j]])])
                        simvalsy.extend([repdata[key1]])
                        labels.extend(celltype)
                    elif key2 in repdata:
                        rawcountsx.extend([min([rawcounts[i],rawcounts[j]])])
                        simvalsy.extend([repdata[key2]])
                        labels.extend(celltype)
                    else:
                        print 'cant find', celltype,i,j
    sns.scatterplot(rawcountsx,simvalsy,hue=['a']*len(rawcountsx),palette=sns.color_palette(['#58508d']),legend=False,s=35)
    plt.tick_params(top=False,right=False,left=False,bottom=False)
    #plt.scatter(rawcountsx,simvalsy, s=5,c='b')                  
    #plt.xlim((0,250000000))                      
    plt.xlim(0,2e9)
    plt.ylim(0.2,1)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.xlabel('Hi-C coverage',size=14)
    plt.ylabel('Replicate '+measure,size=14)
    plt.savefig(figname, bbox_inches='tight')
    print 'Saved raw counts scatter plot figure to',figname
    #slope,intcpt,rval,pval,stderr = scipy.stats.linregress(rawcountsx,simvalsy)
    #print 'R^2 =',rval**2                                   
    r,pval = scipy.stats.spearmanr(rawcountsx,simvalsy)
    print 'rho =', r
    print 'pval =',pval


def protocolRepBoxplots(isreps,dilreps,measure,filename):

    newpalette = ['#8dd3c7','#58508d','#bc5090','#003f5c','#ff6361']
    sns.set(style='whitegrid',palette=sns.color_palette(newpalette),color_codes=False,font='Ubuntu')

    mwuval,pval = scipy.stats.mannwhitneyu(isreps,dilreps,alternative='greater')
    print 'pvalue for in situ/dilution  replicate comparison =',pval
    fig = plt.figure()
    ax = plt.gca()
    protdata = {'labels':['insitu']*len(isreps)+['dilution']*len(dilreps), 'vals':isreps+dilreps}
    sns.violinplot(x=protdata['labels'],y=protdata['vals'],cut=1)
    #plt.boxplot([isreps,dilreps], widths = 0.5,medianprops=dict(color='red'))   
    ax.set_aspect(4)
    ax.set_ylim(0,1)
    plt.xticks([0,1],['In situ \n replicates', 'Dilution \n replicates'],size=11)
    plt.ylabel(measure)
    plt.xlabel('')
    plt.tick_params(top=False,right=False,left=False)
    plt.savefig(filename, bbox_inches = 'tight')
    plt.close(fig)
    print 'saved protocol replicate comparison to',filename

def plotISDilPairs(ispairs,dilpairs,isdil,measure,filename):

    mwuval,pval = scipy.stats.mannwhitneyu(ispairs, dilpairs)
    print 'pvalue for in situ pairs vs dilution pairs =', pval
    mwuval,pval = scipy.stats.mannwhitneyu(ispairs,isdil)
    print 'pvalue for in situ pairs vs is-dilution pairs =', pval
    mwuval,pval = scipy.stats.mannwhitneyu(dilpairs, isdil)
    print 'pvalue for dilution pairs vs is-dilution pairs =',pval

    fig = plt.figure()
    ax = plt.gca()
    protdata = {'labels': ['is']*len(ispairs)+['dil']*len(dilpairs)+['isdil']*len(isdil), 'vals':ispairs+dilpairs+isdil}
    #plt.boxplot(,widths=0.5,medianprops=dict(color='red'))
    sns.violinplot(x=protdata['labels'],y=protdata['vals'],cut=1)
    ax.set_aspect(4)
    ax.set_ylim(0,1)
    plt.xticks([0,1,2],['In situ \n pairs','Dilution \n pairs', 'In situ-dilution \n pairs'],size=11)
    plt.ylabel(measure)
    plt.xlabel('')
    plt.tick_params(top=False,right=False,left=False)
    plt.savefig(filename,bbox_inches='tight')
    plt.close(fig)
    print 'Saved protocol pair figure to',filename


def plotSameCellTypeBoxplots(samecelltypediffexp, samecelltypesameexp, diffcelltype, measure, filename):

    newpalette = ['#8dd3c7','#58508d','#bc5090','#003f5c','#ff6361']
    sns.set(style='whitegrid',palette=sns.color_palette(newpalette),color_codes=False,font='Ubuntu')

    samedistall = [x for lst in samecelltypesameexp for x in lst]
    mwuval,pval = scipy.stats.mannwhitneyu(samedistall, samecelltypediffexp, alternative='greater')
    print 'pvalue for same cell type protocol comparison =', pval
    mwuval,pval = scipy.stats.mannwhitneyu(samedistall, diffcelltype,alternative='greater')
    print 'pvalue for same protocol v different cell types =', pval
    mwuval,pval = scipy.stats.mannwhitneyu(diffcelltype, samecelltypediffexp,alternative='less')
    print 'pvalue for diff protocol v diff cell type =', pval
    #print alldata                       
    fig = plt.figure()
    ax = plt.gca()
    samecelltypesameexp = [x for subl in samecelltypesameexp for x in subl]
    protdata = {'labels':['diffct']*len(diffcelltype)+['samectdiffprot']*len(samecelltypediffexp)+['samectsameprot']*len(samecelltypesameexp), 'vals': diffcelltype+samecelltypediffexp+samecelltypesameexp}
    sns.violinplot(x=protdata['labels'],y=protdata['vals'],cut=1)
    #plt.boxplot([diffcelltype,samedistall,samecelltypediffexp],widths=0.5,medianprops=dict(color='red'))                                                                                          
    ax.set_aspect(4)
    ax.set_ylim(0,1)
    #for idx,pts in enumerate(samedist):                     
    #    x = np.random.normal(1,0.04,size=len(pts))   
    #    plt.scatter(x,pts)   

    plt.xticks([0,1,2],['Different \n cell types','Same cell \n type, different \n protocol','Same cell \n type, same \n protocol'],size=11)
    plt.ylabel(measure)
    plt.xlabel('')
    plt.tick_params(top=False,right=False,left=False)
    plt.savefig(filename,bbox_inches='tight')
    plt.close(fig)
    print 'Saved same cell type protocol figure to',filename


def plotTissueBoxplots(withintissue,acrosstissue,replicates,background, measure, filename):

    newpalette = ['#8dd3c7','#58508d','#bc5090','#003f5c','#ff6361']
    sns.set(style='whitegrid',palette=sns.color_palette(newpalette),color_codes=False,font='Ubuntu')


    mwuval,pvalwithinvsacross = scipy.stats.mannwhitneyu(withintissue, acrosstissue,alternative='greater')
    mwuval,pvalrepvswithin = scipy.stats.mannwhitneyu(withintissue,replicates)
    mwuval,pvalrepvsacross = scipy.stats.mannwhitneyu(acrosstissue,replicates)
    mwuval,pvalbgvswithin = scipy.stats.mannwhitneyu(background, withintissue, alternative='less')
    mwuval,pvalbgvsacross = scipy.stats.mannwhitneyu(background, acrosstissue, alternative='less')

    print 'pvalue within vs across =', pvalwithinvsacross
    print 'pvalue replicates vs within =', pvalrepvswithin
    print 'pvalue replicates vs across =', pvalrepvsacross
    print 'pvalue background vs within =', pvalbgvswithin
    print 'pvalue background vs across =', pvalbgvsacross
    print 'average value among different tissue types:', np.mean(acrosstissue)

    fig = plt.figure()
    ax = plt.gca()
    tissuedata = {'labels':['bg']*len(background)+['reps']*len(replicates)+['intissue']*len(withintissue)+['difftissue']*len(acrosstissue), 'vals':background+replicates+withintissue+acrosstissue}
    sns.violinplot(x=tissuedata['labels'],y=tissuedata['vals'],cut=1)
    #plt.boxplot([background,replicates,withintissue,acrosstissue],widths=0.5,medianprops=dict(color='red'))            
    #for idx,pts in enumerate(sameindiv):              
    #    x = np.random.normal(4, 0.04, size=1)
    #    plt.scatter(x,pts) 
    ax.set_aspect(3)
    ax.set_ylim(0,1)
    plt.xticks([0,1,2,3],['Background','Replicates','Same tissue, \n different donor','Different tissue'],size=13.5)
    plt.ylabel(measure,size=14)
    plt.ylim((0,1))
    plt.xlabel('')
    plt.tick_params(top=False,right=False,left=False)
    plt.savefig(filename,bbox_inches='tight')
    plt.close(fig)
    print 'Saved tissue boxplot figure to',filename

def plotTrioBoxplots(parentchild,parentparent,trioreps,bloodlymph,background, measure, filename):

    newpalette = ['#8dd3c7','#58508d','#bc5090','#003f5c','#ff6361']
    sns.set(style='whitegrid',palette=sns.color_palette(newpalette),color_codes=False,font='Ubuntu')

    fig = plt.figure()
    ax = plt.gca()

    #order points by yvals
    labels = ['pp']*len(parentparent)+['pc']*len(parentchild)
    allvals =  parentparent+parentchild
    labelledvals = [list(x) for x in zip(labels,allvals)]
    labelledvals.sort(key=lambda x: x[1])
    labels = [l[0] for l in labelledvals]
    allvals = [l[1] for l in labelledvals]

    if labels[0] == 'pp':
        ptpalette = newpalette[3:5]
    else:
        ptpalette = [newpalette[4],newpalette[3]]

    #plt.boxplot([background,bloodlymph,trioreps,[]],widths=0.5,medianprops=dict(color='red'))    
    comppts = [parentparent,parentchild]
    xvals = np.linspace(3-.04*len(allvals),3+.04*len(allvals),len(allvals))
    triodata = {'labels': labels,'xvals': xvals,'yvals': allvals}
    sns.scatterplot(x='xvals',y='yvals',hue='labels',data=triodata,legend=False,palette=ptpalette,s=35)
    triodata = {'labels':['bg']*len(background)+['blp']*len(bloodlymph)+['treps']*len(trioreps), 'vals':background+bloodlymph+trioreps}
    sns.violinplot(x=triodata['labels'],y=triodata['vals'],cut=1)
    ax.set_aspect(3)
    ax.set_ylim(0,1)
    plt.xlim(-0.5,3.5)
    plt.xticks([0,1,2,3],['Background','Blood \n lymphocyte \n pairs','Trio \n replicates','Within trio \n pairs'],size=14)
    plt.ylabel(measure,size=14)
    plt.tick_params(axis='y', which='major', labelsize=14)
    plt.ylim((0,1))
    plt.xlabel('')
    plt.tick_params(top=False,right=False,left=False)
    #ax.legend(handles,['parent-parent','parent-child'],loc='lower right',prop={'size': 14},frameon=False)                                                                                         
    plt.savefig(filename,bbox_inches='tight')
    plt.close(fig)
    print 'Saved trio boxplot figure to',filename

