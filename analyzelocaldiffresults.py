import argparse
import numpy as np
import math
import csv
import glob
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
import scipy.stats


def readLocalDiffFiles(filelist,celltypelist):
    compdata = {}
    for filename in filelist:
        splitname = filename.split('_')
        celltype1 = ''
        celltype2 = ''
        for i,split in enumerate(splitname):
            if len(celltype1) > 0 and len(celltype2) > 0: break
            if len(celltype1) == 0 and split in celltypelist:
                celltype1 = split
            elif len(celltype1) == 0 and split+'_D' == splitname[i]+'_'+splitname[i+1]:
                celltype1 = split+'_D'
            elif len(celltype1) == 0 and split+'_LA' == splitname[i]+'_'+splitname[i+1]:
                celltype1 = split+'_LA'
            elif len(celltype1) == 0 and split+'_R' == splitname[i]+'_'+splitname[i+1]:
                celltype1 = split+'_R'
            elif split in celltypelist:
                celltype2 = split
            elif split+'_D' in celltypelist and split+'_D' == splitname[i]+'_'+splitname[i+1]:
                celltype2 = split+'_D'
            elif split+'_LA' in celltypelist and split+'_LA' == splitname[i]+'_'+splitname[i+1]:
                celltype2 = split+'_LA'
            elif split+'_R' in celltypelist and split+'_R' == splitname[i]+'_'+splitname[i+1]:
                celltype2 = split+'_R'
        chrloc = filename.find('chr')
        chrnum = int(filename[chrloc+3:-4])
        dictkey = (celltype1, celltype2, chrnum)
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

def readGeneLocFile(filename='hg19_genelocs.txt'):
    genelocs = []
    with open(filename,'rb') as f:
        freader = csv.reader(f, delimiter='\t')
        freader.next()
        for line in freader:
            try:
                chrnum = int(line[1])
            except:
                continue
            genestart = int(line[2])
            geneend = int(line[3])
            genelocs.append([chrnum, genestart, geneend])
    return genelocs

def countNumGenes(genelocs,cancercolsums,numcancer, numnormal,res):

    numlowcancercons_chr = np.zeros(22)
    totalgenes_chr = np.zeros(22)
    deltadists = [[] for i in xrange(22)]
    for gene in genelocs:
        chrnum = gene[0]
        genestart = gene[1]
        geneend = gene[2]
        genebin = int(round((genestart + geneend)/float(2)/res))
        if genebin > len(cancercolsums['nn'][chrnum-1])-1 or genebin > len(cancercolsums['cc'][chrnum-1])-1: continue
        percsimnormal = cancercolsums['nn'][chrnum-1][genebin]/numnormal
        percsimcancer = cancercolsums['cc'][chrnum-1][genebin]/numcancer
        if percsimnormal > percsimcancer:
            numlowcancercons_chr[chrnum-1] += 1
        totalgenes_chr[chrnum-1] += 1
        deltadists[chrnum-1].extend([percsimnormal - percsimcancer])
    return numlowcancercons_chr,totalgenes_chr,deltadists

def calcPvals(probbychr, totalnumgenes, totallowcancercons, cancergenechrs, cancergenelocs, numlowcancercons):

    # for overall genome p-value:
    pval1 = scipy.stats.hypergeom.sf(numlowcancercons-1,totalnumgenes,totallowcancercons,len(cancergenechrs))

    # for p-value using same chromosomes as top ten cancer genes
    # prob of any 9/10 being higher:
    pval2 = 0
    for i in xrange(len(cancergenechrs)):
        currprob = 1
        for genenum,chrnum in enumerate(cancergenechrs):
            if i != genenum:
                currprob = currprob*probbychr[chrnum-1]
            else:
                currprob = currprob*(1-probbychr[chrnum-1])
        pval2 += currprob
    # prob of all 10/10 being higher:
    proball10 = 1.0
    for i,chrnum in enumerate(cancergenechrs):
        proball10 = proball10*probbychr[chrnum-1]
    pval2 += proball10
    return pval1, pval2

def computeBasicStats(compdata):
    
    avglength = 0
    intcount = 0
    minlength = [(), 1e10]
    maxlength = [(), 0]
    for key,data in compdata.iteritems():
        if len(data) > 0:
            for intdata in data:
                intlen = intdata[1]-intdata[0]+1
                intcount += 1
                avglength += intlen
                if intlen < minlength[1] and intlen > 1:
                    minlength = [key, intlen]
                if intlen > maxlength[1]:
                    maxlength = [key, intlen]
    avglength = avglength/intcount
    print 'Total number of intervals = ', intcount
    print 'Average interval length = ', avglength
    print 'Max interval length = ', maxlength
    print 'Min interval length = ', minlength

def calcBinCons(compdata, chrlengths, res = 100000):
    # format data for Figure 4, calculate average conservation per bin

    chrdict = {} # keys will be chromosome numbers
    for key,data in compdata.iteritems(): # compdata[(celltype1,celltype2,chrnum)] = [[start,end,VI,qval]...]
        chrnum = key[2]
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
    avgperbin = 0
    for chrnum,chrmat in chrdict.iteritems():
        bincons[chrnum-1] = np.sum(chrmat, axis=0)
        avgperbin += np.sum(bincons[chrnum-1])
    print 'Average conservation per bin = ', avgperbin/np.sum(chrlengths)
    return bincons

def calcBinConsbyCondition(compdata, cancertypes, normaltypes, chrlengths, res = 100000):
    # format data for Figure 6, also calculate average conservation for cancer/normal conditions

    conditiondict = {} # keys will be ['cc'](cancer-cancer), ['cn'] (cancer-normal), and ['nn'] (normal-normal)
    conditiondict['cc'] = {}
    conditiondict['cn'] = {}
    conditiondict['nn'] = {}
    for key,data in compdata.iteritems():
        chrnum = key[2]
        if key[0] in cancertypes and key[1] in cancertypes:
            condition = 'cc'
        elif key[0] in normaltypes and key[1] in normaltypes:
            condition = 'nn'
        else:
            condition = 'cn'
        if chrnum not in conditiondict[condition]:
            chrmat = np.zeros((len(compdata),chrlengths[chrnum-1]))
        else:
            chrmat = conditiondict[condition][chrnum]
        rowidx = np.where(~chrmat.any(axis=1))[0][0]
        for interval in data:
            chrmat[rowidx, int(interval[0]):int(interval[1])] = 1
        conditiondict[condition][chrnum] = chrmat

    binconsbycondition = {}
    for condition in conditiondict:
        avgperbin = 0
        bincons = [[] for i in range(22)]
        for chrnum,chrmat in conditiondict[condition].iteritems():
            bincons[chrnum-1] = np.sum(chrmat, axis=0)
            avgperbin += np.sum(bincons[chrnum-1])
        binconsbycondition[condition] = bincons
        if condition == 'cc':
            numpairs = 0.5*len(cancertypes)*(len(cancertypes)-1)
        elif condition == 'cn':
            numpairs = len(cancertypes)*len(normaltypes)
        else:
            numpairs = 0.5*len(normaltypes)*(len(normaltypes)-1)
        print 'Average conservation per bin for '+condition+' cell type pairs =', avgperbin/np.sum(chrlengths)/numpairs
    return binconsbycondition

def plotCancerGeneData(binconsbycondition,filename,topten, chrnums, binnums,numcancer,numnormal,probbychr,totalprob,deltadists):
    # Figure 7 in manuscript

    nnpairs = [binconsbycondition['nn'][chrnums[i]-1][binnum]/numnormal for i,binnum in enumerate(binnums)]
    ccpairs = [binconsbycondition['cc'][chrnums[i]-1][binnum]/numcancer for i,binnum in enumerate(binnums)]

    cancergenedelta = np.array(nnpairs) - np.array(ccpairs)
    fig = plt.figure()
    ax = plt.gca()
    plt.boxplot(deltadists)
    plt.scatter(chrnums, cancergenedelta, s=15, color='red', marker='*')
    #ax.set_ylim(0,1)
    plt.xlabel('Chromosome')
    plt.ylabel('% normal-normal similarity - % cancer-cancer similarity')
    distfilename = filename[:-21] + 'deltadistplot_allchr.eps'
    plt.savefig(distfilename,bbox_inches='tight')
    plt.close(fig)
    print 'Saved deltadistribution fig for all chr to',distfilename

    idx = np.arange(len(nnpairs))
    fig,ax = plt.subplots()
    width=0.3
    bars1 = ax.bar(idx,nnpairs,width,color='purple')
    bars3 = ax.bar(idx+width,ccpairs,width,color='blue')
    ax.set_ylabel('Fraction of structurally similar pairs')
    ax.set_xticks(idx+.5*width)
    ax.set_xticklabels(topten)
    ax.legend((bars1[0],bars3[0]),('normal-normal','cancer-cancer'))
    plt.savefig(filename,bbox_inches='tight')
    plt.close(fig)
    print 'Saved cancer gene plot to',filename
    
def plotBinConservation(bincons, res, centromerelocs, filename):
    #Figure 4 in manuscript
    #bincons is a list with 22 rows where each row represents the conservation vector of a chromosome, and each value in the vector represents a genomic bin, and the number of pairwise comparisons in which that bin was within a significant, dominating interval
    for i,row in enumerate(centromerelocs):
        centromerelocs[i] = [np.floor(row[0]/res), np.ceil(row[1]/res)]
    for i,chrdata in enumerate(bincons):
        fig = plt.figure(figsize=(10,3.3))
        ax = plt.gca()
        x_pos = np.arange(len(chrdata))
        plt.bar(x_pos, chrdata, 1)
        plt.axvline(x=centromerelocs[i][0],color='r')
        plt.axvline(x=centromerelocs[i][1],color='r')
        ax.set_xlim(0,len(chrdata))
        ax.set_ylim(0,253)
        ax.tick_params(axis='both',labelsize=8)
        plt.ylabel('Number of similar pairs of cell types',fontsize=10)
        plt.xlabel('Genomic position (100kb)',fontsize=10)
        figname = filename+'_chr'+str(i+1)+'.eps'
        plt.savefig(figname,bbox_inches='tight')
        plt.close(fig)
        print 'Wrote num pairs histogram for chr'+str(i+1)+' to',figname

def plotChrLevelBoxplot(chrdata, totalpairs, filename):
    # Figure 6 in manuscript
    # input: chrdata is a list of 22 lists, one for each chromosome, where each list is the length of the number of bins on that chromosome, and the values represent the number of pairs in which that particular bin was in a significant, dominating interval
    # totalpairs: integer representing the total number of pairs for the given condition of chrdata (for us, 28 for normal-normal pairs, 105 for cancer-cancer pairs)

    chrdata = np.array(chrdata)/float(totalpairs)
    fig = plt.figure()
    ax = plt.gca()
    plt.boxplot(chrdata)
    ax.set_ylim(0,1)
    plt.xlabel('Chromosome')
    plt.ylabel('Fraction of structurally similar pairs')
    plt.savefig(filename,bbox_inches='tight')
    plt.close(fig)
    print 'Saved chromosome-level similarity fig to',filename

def calcPercSimilarity(compdata, chrlengths):
    # for data in Table 2, other values cited in mansucript (such as comparison to Dixon/Rao)
    # compdata: dictionary with ((celltype1,celltype2),chrnum) where each entry is the list of significant, dominating intervals [[start,end,VI,pval],...]

    totallen = np.sum(chrlengths)
    totalsim = [] #will be list of [(celltype1, celltype2), sum] values
    simbychr = {} #dict keys will be chrnum, entry will be list of same format as totalsim
    for key,data in compdata.iteritems():
        chrnum = key[2]
        if len(data) > 0:
            sigints = np.zeros(int(max(data[:,1])))
            for row in data:
                sigints[int(row[0]):int(row[1])] = 1
            chrsum = np.sum(sigints)
        else:
            chrsum = 0
        celltypes = (key[0], key[1])
        if chrnum in simbychr:
            simbychr[chrnum].append([celltypes, chrsum/chrlengths[chrnum-1]])
        else:
            simbychr[chrnum] = [[celltypes, chrsum/chrlengths[chrnum-1]]]
        idx = -1
        for i,sublist in enumerate(totalsim):
            if sublist[0] == celltypes:
                idx = i
        if idx > -1:
            totalsim[idx][1] += chrsum/totallen
        else:
            totalsim.append([celltypes, chrsum/totallen])
    #totalsimsorted = sorted(totalsim,key=lambda x: x[1], reverse=True)
    #for idx,x in enumerate(totalsimsorted):
    #    if idx < 10:
    #        print x
    #    if x[0] == ('IMR90_D', 'IMR90_R'):
    #        print x, idx
    #    if x[0] == ('K562_LA', 'K562_R'):
    #        print x, idx
    #sys.exit()
    return totalsim,simbychr

def analyzePercSimByChr(simbychr,numpairs):
    # for values mentioned in section 3.2, chromosome-level similarity

    #sort data within each chromosome
    for chrnum,data in simbychr.iteritems():
        data = sorted(data, key=lambda x:x[1], reverse=True)
        simbychr[chrnum] = data
    # calculate average sim per chromosome, highest max, lowest min
    avgperchr = np.zeros(22)
    maxsim = [(),0,0]
    minsim = [(),1000,1000]
    zeropairs = []
    for chrnum,data in simbychr.iteritems():
        for pair in data:
            avgperchr[chrnum-1] += pair[1]/numpairs
            if pair[1] > maxsim[2]:
                maxsim[0] = pair[0]
                maxsim[1] = chrnum
                maxsim[2] = pair[1]
            if pair[1] < minsim[2]:
                minsim[0] = pair[0]
                minsim[1] = chrnum
                minsim[2] = pair[1]
            if pair[1] == 0:
                zeropairs.append([pair[0], chrnum, pair[1]])
    print 'Maximally similar pair:', maxsim
    print 'Minimally similar pair:', minsim
    print 'Number of zero pairs:',len(zeropairs)
    # find which chromosomes have highest and lowest average similarity
    minavg = [0,1]
    maxavg = [0,0]
    for i,avgval in enumerate(avgperchr):
        if avgval > maxavg[1]:
            maxavg = [i+1, avgval]
        if avgval < minavg[1]:
            minavg = [i+1, avgval]
    print 'Max avg on chr'+str(maxavg[0])+' = ',maxavg[1]
    print 'Min avg on chr'+str(minavg[0])+' = ',minavg[1]

def printDixonRaoPairs(percsimlist):

    for ctpair in percsimlist:
        if ctpair[0] == ('IMR90dixon','hESC') or ctpair[0] == ('IMR90_D', 'hESC'):
            print ctpair
        if ctpair[0] == ('IMR90rao','GM12878') or ctpair[0] == ('IMR90_R', 'GM12878'):
            print ctpair
        if ctpair[0][0] == 'GM12878':
            if ctpair[0][1] == 'HMEC':
                print ctpair
            if ctpair[0][1] == 'HUVEC':
                print ctpair
            if ctpair[0][1] == 'K562rao' or ctpair[0][1] == 'K562_R':
                print ctpair
            if ctpair[0][1] == 'KBM7':
                print ctpair
            if ctpair[0][1] == 'NHEK':
                print ctpair

def plotPercSimHeatMap(percsimlist, filename):
    # Figure 5 in manuscript

    celltypelist = []
    for line in percsimlist:
        celltypelist.extend([line[0][0], line[0][1]])
    celltypelist = list(set(celltypelist))
    print celltypelist
    #put percsim data into matrix
    simmat = np.ones((len(celltypelist),len(celltypelist)))
    totalsim = sorted(percsimlist, key=lambda x:x[1], reverse=True)
    # print ranks of IMR90 pairs and K562 pairs
    for idx,ctpair in enumerate(totalsim):
 #       print ctpair[0]
        if ctpair[0] == ('IMR90rao','IMR90dixon') or ctpair[0] == ('IMR90_D', 'IMR90_R'):
            print ctpair, idx
        if ctpair[0] == ('K562rao', 'K562la') or ctpair[0] == ('K562_LA', 'K562_R'):
            print ctpair, idx
        if ctpair[0] == ('K562rao', 'KBM7') or ctpair[0] == ('K562_R', 'KBM7'):
            print ctpair, idx
        simmat[celltypelist.index(ctpair[0][0]), celltypelist.index(ctpair[0][1])] = ctpair[1]
        simmat[celltypelist.index(ctpair[0][1]), celltypelist.index(ctpair[0][0])] = ctpair[1]
    # sort matrix by cell types with highest average similarity
    simsum = np.sum(simmat, axis=0)
    idx = (-simsum).argsort()
    simmat = simmat[:,idx]
    simmat = simmat[idx,:]
    # fix some of the cell type names
    for i,celltype in enumerate(celltypelist):
        if celltype == 'IMR90_R' or celltype == 'IMR90rao':
            celltypelist[i] = 'IMR90(R)'
        if celltype == 'K562_R' or celltype == 'K562rao':
            celltypelist[i] = 'K562(R)'
        if celltype == 'IMR90_D' or celltype == 'IMR90dixon':
            celltypelist[i] = 'IMR90(D)'
        if celltype == 'K562_LA' or celltype == 'K562la':
            celltypelist[i] = 'K562(LA)'
    celllist = [celltypelist[i] for i in idx]
    plotHeatMap(simmat,celllist,filename)

def writeToFile(data, filename):
    
    with open(filename,'wb') as f:
        fwriter = csv.writer(f, delimiter='\t')
        for row in data:
            fwriter.writerow(row)
    print 'Wrote data to', filename

def plotHeatMap(data,labels,filename):
    
    fig,ax1 = plt.subplots(1,1)
    plt.imshow(data, interpolation='none',cmap='YlGnBu')
    cbar = plt.colorbar()
    cbar.set_label('% similarity', rotation=270,labelpad=25)
    ax1.set_xticks(list(range(len(labels))))
    ax1.set_yticks(list(range(len(labels))))
    ax1.set_xticklabels(labels,size='small',rotation=90)
    ax1.set_yticklabels(labels,size='small')
    ax1.grid(False)
    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    plt.close(fig)
    print 'Heatmap saved as',filename

def plotIntervalResults(compdata,tadlists,chrlength,newdata,figname):
    # used to make parts of Figures 1 and 3

    matsize = int(chrlength) # should already be divided by resolution
    # plot TAD sets
    for interval in tadlists[0]:
        plt.plot((interval[0]-1,interval[1]), (matsize-interval[0] + 1, matsize-interval[0]+1),'g')
        plt.plot((interval[1],interval[1]),(matsize-interval[0]+1,matsize-interval[1]),'g')
        dom1Artist = plt.Line2D((0,1),(0,0), color='green', linestyle='solid')
    for interval in tadlists[1]:
        plt.plot((interval[0]-1,interval[1]),(matsize-interval[1],matsize-interval[1]),'b')
        plt.plot((interval[0]-1,interval[0]-1),(matsize-interval[0]+1,matsize-interval[1]),'b')
        dom2Artist = plt.Line2D((0,1),(0,0), color='blue', linestyle='solid')
    ax = plt.gca()
    # plot stars over newdata points (for Figure 1, used all boundary pts, then just significant points, then just dominating ones)
    if len(newdata)>0:
        for interval in newdata:
            plt.plot((interval[0]-1,interval[1]),(matsize-interval[1],matsize-interval[1]),'r')
            plt.plot((interval[0]-1,interval[0]-1),(matsize-interval[0]+1,matsize-interval[1]),'r')
            dom3Artist = plt.Line2D((0,1),(0,0), color='red', linestyle='solid')
#        ax.scatter(newdata[:,1], matsize-newdata[:,0], s=15, color='red', marker='*')
    else: # to just plot TAD sets alone
        ax.set_facecolor('white')
        plt.tick_params(axis='x',labelbottom='off')
    ax.set_aspect('equal', adjustable='box')
    plt.tick_params(axis='y',labelleft='off')
    plt.savefig(figname,bbox_inches='tight')
    plt.close()
    print 'Saved figure as',figname


def readTADfile(filename, res):
    # reads Armatus files
    tadbdys = []
    with open(filename, 'rb') as tf:
        freader = csv.reader(tf, delimiter = '\t')
        for line in freader:
            tadbdys.append([int(line[1]), int(line[2])])
    if tadbdys[0][0] > tadbdys[-1][1]:
        tadbdys = np.flipud(np.array(tadbdys))
    tadbdys = np.array(tadbdys)
    tadbdys[:,0] = np.floor(np.divide(tadbdys[:,0],res))
    tadbdys[:,1] = np.floor(np.divide((tadbdys[:,1]+1),res)-1)
    return tadbdys

def readIntermedFiles(filename):
    # used to read the files with the data for Figure 1 - intermediate results from the method

    newdata = []
    with open(filename,'rb') as f:
        freader = csv.reader(f,delimiter='\t')
        freader.next()
        for line in freader:
            newdata.append([int(line[0]), int(line[1]), float(line[2])])
    return np.array(newdata)


def main(fileseed, res, cancertypes, normaltypes, armatusfilename, chrlengthfile, centromerefile, genelocfile, scatterplotfile, outputloc):

    if len(cancertypes) == 0:
        cancer = ['K562_LA','K562_R','KBM7','A549','Caki2','G401','LNCaP-FGC','NCI-H460','Panc1','RPMI-7951','SJCRH30','SKMEL5','SKNDZ','SKNMC','T47D']
    else:
        cancer = cancertypes
    if len(normaltypes) == 0:
        noncancer = ['IMR90_R','GM12878','HMEC','HUVEC','NHEK','IMR90_D','hESC','GM06990']
    else:
        noncancer = normaltypes

    # if no file for chrlengths and centromere locations are input, use these (hg19)
    if len(chrlengthfile) == 0:
        chrlengths = [248956422, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566] # hg19 chr lengths from https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
        for i,chrlen in enumerate(chrlengths):
            chrlengths[i] = int(np.ceil(float(chrlen)/res))
    else:
        pass
        #read file

    if len(centromerefile) == 0:
        centromerelocs = [[121535434, 124535434],[92326171, 95326171],[90504854,93504854],[49660117, 52660117],[46405641, 49405641],[58830166, 61830166],[58054331, 61054331],[43838887, 46838887],[47367679, 50367679],[39254935,42254935],[51644205,54644205],[34856694,37856694],[16000000,19000000],[16000000,19000000],[17000000,20000000],[35335801,38335801],[22263006,25263006],[15460898,18460898],[24681782,27681782],[26369569,29369569],[11288129,14288129],[13000000,16000000]] #hg19 centromere locations from http://genome.ucsc.edu/cgi-bin/hgTables
    else:
        pass
        #read file

    if len(genelocfile) == 0:
        genelocfile = 'hg19_genelocs.txt'

    filelist = glob.glob(fileseed)
    compdata = readLocalDiffFiles(filelist,cancer+noncancer)
    numpairs = len(compdata)/22 # should be 253 in our case

    # to get data for sec 3.1 of manuscript
    computeBasicStats(compdata)
    bincons = calcBinCons(compdata, chrlengths, res)
    binconsconditional = calcBinConsbyCondition(compdata, cancer, noncancer, chrlengths, res)
    # plot all chromosomes for Figure 4
    plotBinConservation(bincons, res, centromerelocs, outputloc+'binconservationplots')
    # calculate percent similarities
    totalsim, simbychr = calcPercSimilarity(compdata, chrlengths)
    totalsimsorted = sorted(totalsim,key=lambda x: x[1], reverse=True)

    # compare percsim to Dixon and Rao results
    printDixonRaoPairs(totalsim)
    # analyze percent similarity data
    analyzePercSimByChr(simbychr, numpairs)
    # make heat map (Figure 5)
    plotPercSimHeatMap(totalsim, outputloc+'percsimheatmap.eps')
    # make boxplots for Figure 6
    numcancerpairs = 0.5*(len(cancer))*(len(cancer)-1)
    numnormalpairs = 0.5*(len(noncancer))*(len(noncancer)-1)
    for condition,condbincons in binconsconditional.iteritems():
        if condition == 'cc':
            totalpairs = numcancerpairs
            filename = outputloc+'chrlevel_boxplot_cancercancer.eps'
        elif condition == 'nn':
            totalpairs = numnormalpairs
            filename = outputloc+'chrlevel_boxplot_normalnormal.eps'
        else:
            continue
        plotChrLevelBoxplot(condbincons, totalpairs, filename)
    

    # analyze cancer v noncancer data
    toptencancer = ('TP53','PIK3CA','PTEN','APC','VHL','KRAS','KMT2C','KMT2D','ARID1A','PBRM1')
    cancerchrnums = [17, 3, 10, 5, 3, 12, 7, 12, 1, 3]
    # gene locations from ghr.nlm.nih.gov
    cancergenelocs = [[7668402, 7687550], [179148114, 179240093], [87863438, 87971930], [112707505, 112846239], [10141635, 10153670], [25204789, 25252093], [152134925, 152436642], [49018975, 49060884], [26696031, 26782110], [52545352, 52685933]]
    cancerbinnums = [(loc[0]+loc[1])/2/res for loc in cancergenelocs]

    genelocs = readGeneLocFile(genelocfile)
    probbychr,genesperchr,deltadists = countNumGenes(genelocs, binconsconditional, numcancerpairs, numnormalpairs, res)
    totalprob = np.sum(probbychr)/np.sum(genesperchr)
    print 'Percent of all human genes with higher structural conservation in normal-normal pairs = ',totalprob
    overallpval, chrpval = calcPvals(np.divide(probbychr,genesperchr), np.sum(genesperchr), np.sum(probbychr), cancerchrnums, cancerbinnums, 9)
    print 'Probability of at least 9/10 random genes with higher similarity among normal pairs =',overallpval
    print 'Probability of at least 9/10 genes from the same chromosomes as top 10 with higher similarity among normal pairs =',chrpval
    # plot Figure 7
    filename = outputloc+'cancergenebarplot.eps'
    plotCancerGeneData(binconsconditional, filename, toptencancer, cancerchrnums, cancerbinnums, numcancerpairs, numnormalpairs, probbychr, totalprob, deltadists)


    # make figures used for Figs 1 and 3
    if len(armatusfilename) > 0:
        tadlists = []
        for filename in armatusfilename:
            tadlists.append(readTADfile(filename,res))
        # plot TAD sets alone
        matsize = np.max([np.max(tadlists[0]), np.max(tadlists[1])])
        plotIntervalResults(compdata,tadlists,matsize,[],outputloc+'TADsetplot.eps')
        if len(scatterplotfile) > 0:
            newdata = readIntermedFiles(scatterplotfile)
            plotIntervalResults(compdata,tadlists,matsize,newdata,outputloc+'TADset_plusscatterpts.eps')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='fileseed of files to analyze, with wildcard characters')
    # fileseed, ie '/home/go/outputs/localdiff/localdiff_*optgamma_chr*.txt'
    parser.add_argument('-r', type=int, help='resolution of Hi-C data')
    parser.add_argument('-a', default=[], nargs='+', help='Armatus files for plotting TADs')
    parser.add_argument('-cl', default='', type=str, help='File containing chromosome lengths')
    parser.add_argument('-cm', default='', type=str, help='File containing centromere locations')
    parser.add_argument('-gl', default='', type=str, help='File listing all human gene locations')
    parser.add_argument('-p', default='', type=str, help='File containing points to plot in matrix (like in Fig 1)')
    parser.add_argument('-c', default=[], nargs='+', help='Cancer cell types')
    parser.add_argument('-n', default=[], nargs='+', help='Normal cell types')
    parser.add_argument('-o', default='', type=str, help='Path to location for output files/figures to be written')

    args = parser.parse_args()

    main(args.i, args.r, args.c, args.n, args.a, args.cl, args.cm, args.gl, args.p, args.o)

