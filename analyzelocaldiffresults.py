import argparse
import numpy as np
import csv
import glob
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
import scipy.stats

# input should be fileseed, ie /mnt/disk34/user/nsauerwa/go/outputs/localdiff/localdiff_*_100kb_pvals_optgamma_chr*.txt

def readLocalDiffFiles(filelist):
    compdata = {}
    for filename in filelist:
        splitname = filename.split('_')
        celltype1 = splitname[1]
        celltype2 = splitname[2]
        chrloc = filename.find('chr')
        chrnum = int(filename[chrloc+3:-4])
        dictkey = (celltype1, celltype2, chrnum)
        dictdata = []
        with open(filename, 'rb') as f:
            freader = csv.reader(f, delimiter='\t')
            freader.next() #skip header line
            for line in freader:
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

def countNumGenes(genelocs,cancercolsums,res=100000):

    numlowcancersim = 0
    totalgenes = 0
    for gene in genelocs:
        chrnum = gene[0]
        genestart = gene[1]
        geneend = gene[2]
        genebin = int(round((genestart + geneend)/float(2)/res))
        if genebin > len(cancercolsums[chrnum]['nvn'])-1: continue
        #print genebin
        #print cancercolsums[chrnum]['nvn'][genebin]/28
        #print cancercolsums[chrnum]['cvc'][genebin]/105
        if cancercolsums[chrnum]['nvn'][genebin]/28 > cancercolsums[chrnum]['cvc'][genebin]/105:
            numlowcancersim += 1
        totalgenes += 1
    perclowcancersim = numlowcancersim/float(totalgenes)
    return perclowcancersim

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


def intervalHistogram(compdata, chrlengths, res = 100000):

    # colorcode pairs by cancer-cancer, cancer-noncancer, noncancer-noncancer
    cancer = ['K562la','K562rao','KBM7','A549','Caki2','G401','LNCaP-FGC','NCI-H460','Panc1','RPMI-7951','SJCRH30','SKMEL5','SKNDZ','SKNMC','T47D']
    noncancer = ['IMR90rao','GM12878','HMEC','HUVEC','NHEK','IMR90dixon','hESC','GM06990']
    cancerdict = {}
    chrdict = {} # keys are chrnum
    pairdict = {}
    for key,data in compdata.iteritems(): # compdata[(celltype1,celltype2,chrnum)] = [[start,end,VI,qval]...]
        chrnum = key[2]
        # for each chromosome, make a matrix with a row for each key, and columns for each bin, indicator if bin is included in overlaps
        if chrnum not in chrdict:
            chrmat = np.zeros((len(compdata),chrlengths[chrnum-1]))
            keylist = [[] for i in range(len(compdata))]
        else:
            chrmat = chrdict[chrnum][0]
            keylist = chrdict[chrnum][1]
        if chrnum not in pairdict:
            pairlist = [set() for i in range(chrlengths[chrnum-1])]
        else:
            pairlist = pairdict[chrnum]
        rowidx = np.where(~chrmat.any(axis=1))[0][0]
        keylist[rowidx] = key
        #keep list of which key each row belongs to? per chromosome? then classify by cancer/noncancer?
        for interval in data:
            chrmat[rowidx,int(interval[0]):int(interval[1])] = 1
            for i in xrange(int(interval[0]),int(interval[1])):
#                print interval[0:2], chrlengths[chrnum-1], chrnum
                pairlist[i].add((key[0],key[1]))
            # add (celltype1, celltype2) to set list?
        chrdict[chrnum] = [chrmat,keylist]
        pairdict[chrnum] = pairlist
    colsums = [[] for i in range(22)] # initialize empty list of length 22 for each chromosome #
    setdiffs = [[] for i in range(22)]
    cancercolsums = {}
    avgperbin = 0
    avgpctsim_cc = 0
    avgpctsim_cn = 0
    avgpctsim_nn = 0
    for chrnum,data in chrdict.iteritems():
        chrmat = data[0]
        keylist = data[1]
        cvc_idx, cvn_idx, nvn_idx = [],[],[]
        for i,key in enumerate(keylist):
            if len(key)>0:
                if key[0] in cancer and key[1] in cancer:
                    cvc_idx.extend([i])
                elif (key[0] in cancer and key[1] in noncancer) or (key[1] in cancer and key[0] in noncancer):
                    cvn_idx.extend([i])
                elif key[0] in noncancer and key[1] in noncancer:
                    nvn_idx.extend([i])
                else:
                    print 'that wasnt supposed to happen'
                    print key
                    sys.exit()
        cancercolsums[chrnum] = {}
#        print cvc_idx
        cancercolsums[chrnum]['cvc'] = np.sum(chrmat[cvc_idx,:], axis=0)
        avgpctsim_cc += np.sum(cancercolsums[chrnum]['cvc'])/105
        cancercolsums[chrnum]['cvn'] = np.sum(chrmat[cvn_idx,:], axis=0)
        avgpctsim_cn += np.sum(cancercolsums[chrnum]['cvn'])/120
        cancercolsums[chrnum]['nvn'] = np.sum(chrmat[nvn_idx,:], axis=0)
        avgpctsim_nn += np.sum(cancercolsums[chrnum]['nvn'])/28
        colsums[chrnum-1] = np.sum(chrmat, axis=0)
        avgperbin += np.sum(colsums[chrnum-1])
    if False:
        print 'at TP53 gene location (chr 17, bin 76-77)'
        print 'prop of normal pairs:',cancercolsums[17]['nvn'][76:78]/28
        print 'prop of n-c pairs:',cancercolsums[17]['cvn'][76:78]/120
        print 'prop of cancer pairs:',cancercolsums[17]['cvc'][76:78]/105
        print 'at RB1 gene location (chr 13, bin 483-485)'
        print 'prop of normal pairs:',cancercolsums[13]['nvn'][483:486]/28
        print 'prop of n-c pairs:',cancercolsums[13]['cvn'][483:486]/120
        print 'prop of cancer pairs:',cancercolsums[13]['cvc'][483:486]/105
        print 'at PTEN gene location (chr 10, bin 878-879)'
        print 'prop of normal pairs:',cancercolsums[10]['nvn'][878:880]/28
        print 'prop of n-c pairs:',cancercolsums[10]['cvn'][878:880]/120
        print 'prop of cancer pairs:',cancercolsums[10]['cvc'][878:880]/105
        print 'at PIK3CA gene location (chr 3, bin 1791-1792)'
        print 'prop of normal pairs:',cancercolsums[3]['nvn'][1791:1793]/28
        print 'prop of n-c pairs:',cancercolsums[3]['cvn'][1791:1793]/120
        print 'prop of cancer pairs:',cancercolsums[3]['cvc'][1791:1793]/105
        print 'at MYC gene location (chr 8, bin 1277)'
        print 'prop of normal pairs:',cancercolsums[8]['nvn'][1277]/28
        print 'prop of n-c pairs:',cancercolsums[8]['cvn'][1277]/120
        print 'prop of cancer pairs:',cancercolsums[8]['cvc'][1277]/105
        print('')
        print 'Average percent similarity for cancer-cancer pairs =',avgpctsim_cc/np.sum(chrlengths)
        print 'Average percent similarity for cancer-normal pairs =',avgpctsim_cn/np.sum(chrlengths)
        print 'Average percent similarity for normal-normal pairs =',avgpctsim_nn/np.sum(chrlengths)
#    sys.exit()
#    print 'Average conservation per bin = ', avgperbin/np.sum(chrlengths)
    for chrnum,pairlist in pairdict.iteritems():
        setdiffs[chrnum-1] = [len(pairlist[i]^pairlist[i+1]) for i in xrange(len(pairlist)-1)]
    return colsums, setdiffs, cancercolsums

def plotCancerGeneData(cancercolsums):

    topten = ('TP53','PIK3CA','PTEN','APC','VHL','KRAS','KMT2C','KMT2D','ARID1A','PBRM1')
    chrnums = [17, 3, 10, 5, 3, 12, 7, 12, 1, 3]
    binnums = [76, 1791, 879, 1127, 101, 252, 1522, 490, 267, 526]

    nnpairs = [cancercolsums[chrnums[i]]['nvn'][binnum]/28 for i,binnum in enumerate(binnums)]
    ccpairs = [cancercolsums[chrnums[i]]['cvc'][binnum]/105 for i,binnum in enumerate(binnums)]

#    nnpairs = [cancercolsums[17]['nvn'][76]/28, cancercolsums[13]['nvn'][483]/28, cancercolsums[10]['nvn'][879]/28, cancercolsums[3]['nvn'][1791]/28, cancercolsums[8]['nvn'][1277]/28]
#    cnpairs = [cancercolsums[17]['cvn'][76]/120, cancercolsums[13]['cvn'][483]/120, cancercolsums[10]['cvn'][879]/120, cancercolsums[3]['cvn'][1791]/120, cancercolsums[8]['cvn'][1277]/120]
#    ccpairs = [cancercolsums[17]['cvc'][76]/105, cancercolsums[13]['cvc'][483]/105, cancercolsums[10]['cvc'][879]/105, cancercolsums[3]['cvc'][1791]/105, cancercolsums[8]['cvc'][1277]/105]

 #   xlabels=('TP53','RB1','PTEN','PIK3CA','Myc')
    idx = np.arange(len(nnpairs))
    fig,ax = plt.subplots()
    width=0.3
    bars1 = ax.bar(idx,nnpairs,width,color='purple')
 #   bars2 = ax.bar(idx+width,cnpairs,width,color='green')
    bars3 = ax.bar(idx+width,ccpairs,width,color='blue')
    ax.set_ylabel('Fraction of structurally similar pairs')
    ax.set_xticks(idx+.5*width)
    ax.set_xticklabels(topten)
    #ax.legend((bars1[0],bars2[0],bars3[0]),('normal-normal','cancer-normal','cancer-cancer'))
    ax.legend((bars1[0],bars3[0]),('normal-normal','cancer-cancer'))
    filename = '/mnt/disk34/user/nsauerwa/go/outputs/localdiff/figures/cancergenebarplot.eps'
    plt.savefig(filename,bbox_inches='tight')
    print 'Saved cancer gene plot to',filename


def histogramByCellType(compdata, chrlengths):
    # go through all 23 cell types, for each bin on each chromosome create list of other cell types that were similar in this bin
    # dict[celltype][chr][bins[other celltypes]]
    
    celltypehist = {}
    for key,data in compdata.iteritems():
        chrnum = key[2]
        celltype1 = key[0]
        celltype2 = key[1]
        if celltype1 not in celltypehist:
            celltypehist[celltype1] = {}
        if chrnum not in celltypehist[celltype1]:
            matchlist1 = [[] for i in range(chrlengths[chrnum-1])]
        else:
            matchlist1 = celltypehist[celltype1][chrnum]
        
        if celltype2 not in celltypehist:
            celltypehist[celltype2] = {}
        if chrnum not in celltypehist[celltype2]:
            matchlist2 = [[] for i in range(chrlengths[chrnum-1])]
        else:
            matchlist2 = celltypehist[celltype2][chrnum]
        
        # loop through and add celltype1 to all matchlist2 bins in intervals, and celltype2 to all matchlist1 bins in intervals
        for interval in data:
            for i in range(int(interval[0]), int(interval[1])):
                matchlist1[i].extend([celltype2])
                matchlist2[i].extend([celltype1])

        celltypehist[celltype1][chrnum] = matchlist1
        celltypehist[celltype2][chrnum] = matchlist2
    return celltypehist
    
    
def plotHistograms(colsums, setdiffs, cancercolsums, celltypehists):

    # colsums is a list of lists for each chromosome
    if True:
        res = 100000
        centromerelocs = [[121535434, 124535434],[92326171, 95326171],[90504854,93504854],[49660117, 52660117],[46405641, 49405641],[58830166, 61830166],[58054331, 61054331],[43838887, 46838887],[47367679, 50367679],[39254935,42254935],[51644205,54644205],[34856694,37856694],[16000000,19000000],[16000000,19000000],[17000000,20000000],[35335801,38335801],[22263006,25263006],[15460898,18460898],[24681782,27681782],[26369569,29369569],[11288129,14288129],[13000000,16000000]]
        for i,row in enumerate(centromerelocs):
            centromerelocs[i] = [np.floor(row[0]/res), np.ceil(row[1]/res)]
        for i,chrdata in enumerate(colsums):
            fig = plt.figure(figsize=(10,3.3))
            ax = plt.gca()
            x_pos = np.arange(len(chrdata))
            plt.bar(x_pos, chrdata, 1)
            plt.axvline(x=centromerelocs[i][0],color='r')
            plt.axvline(x=centromerelocs[i][1],color='r')
            ax.set_xlim(0,len(chrdata))
            ax.set_ylim(0,253)
            ax.tick_params(axis='both',labelsize=4)
            plt.ylabel('Number of similar pairs of cell types',fontsize=6)
            plt.xlabel('Genomic position (100kb)',fontsize=6)
            #plt.title('chromosome'+str(i+1))
            #figname = 'testplt_chr'+str(i+1)+'.png'
            figname = '/mnt/disk34/user/nsauerwa/go/outputs/localdiff/figures/numpairshistogram_chr'+str(i+1)+'.eps'
            plt.savefig(figname,bbox_inches='tight')
            plt.close(fig)
            print 'Wrote num pairs histogram for chr'+str(i+1)+' to',figname
            #if i==2:
            #    sys.exit()
    sys.exit()

    if False:
        for i,chrdata in enumerate(setdiffs):
            fig = plt.figure()
            x_pos = np.arange(len(chrdata))
            plt.bar(x_pos,chrdata)
            plt.ylabel('size of set difference')
            plt.xlabel('between genomic bins')
            plt.title('set differences, chromosome'+str(i+1))
            figname = '/mnt/disk34/user/nsauerwa/go/outputs/localdiff/figures/setdiffhistogram_chr'+str(i+1)+'.png'
            plt.savefig(figname)
            plt.close(fig)
            print 'Wrote set diff histogram for chr'+str(i+1)+' to',figname
   

    if False:
        #for each celltype and chromosome, go through and create a binary list for each other celltype, plot them on top of each other
        #need 22 different colors
        colorlist = ['aqua','black','blue','brown','coral','darkblue','darkgreen','fuchsia','gold','green','grey','indigo','lavender','lime','maroon','olive','orange','pink','purple','red','yellow','yellowgreen','teal']
        celltypelist = []
        for celltype in celltypehists:
            celltypelist.extend([celltype])
        for maincelltype,chrmap in celltypehists.iteritems():
            for chrnum,data in chrmap.iteritems():
                # go through and create a list of zeros and ones for presence of each other cell type
                # should be 22 other cell types
                compdict = {}
                for celltype in celltypelist:
                    compdict[celltype] = np.zeros(len(data))
                    #for all other cell types
                for binnum,sublist in enumerate(data):
                    for matchtype in sublist:
                        compdict[matchtype][binnum] = 1
                    fig = plt.figure()
                    x_pos = np.arange(len(data))
                    prev = np.zeros(len(data))
                for celltype in celltypelist:
                    #print colorlist[celltypelist.index(celltype)]
                    plt.bar(x_pos, compdict[celltype], 1, color=colorlist[celltypelist.index(celltype)], bottom=list(prev))
                    plt.title(maincelltype+' vs all, chr'+str(chrnum+1))
                    prev += compdict[celltype]
                    plt.legend(celltypelist)
                filename = '/mnt/disk34/user/nsauerwa/go/outputs/localdiff/figures/one_v_all_plot_'+maincelltype+'_chr'+str(chrnum)+'.png'
                print 'Wrote barplot to',filename
                plt.savefig(filename,bbox_inches='tight')
                plt.close(fig)
                #            sys.exit()

    #print cancer v non cancer histograms
    if False:
        for chrnum,chrdata in cancercolsums.iteritems():
            fig = plt.figure()
            ax = plt.gca()
            x_pos = np.arange(len(chrdata['cvc']))
            #plt.bar(x_pos, chrdata['cvc'], 1, color='blue')
            #plt.bar(x_pos, chrdata['cvn'], 1, bottom=chrdata['cvc'], color='green')
            #plt.bar(x_pos, chrdata['nvn'], 1, bottom=chrdata['cvc']+chrdata['cvn'], color='purple')
            labels=['both cancer', 'one cancer, one normal', 'both normal']
            totalnnpairs = 28 
            totalcnpairs = 120
            totalccpairs = 105
            plt.plot(x_pos, (chrdata['cvc']/totalccpairs), color='blue')
            plt.plot(x_pos, (chrdata['cvn']/totalcnpairs), color='green')
            plt.plot(x_pos, (chrdata['nvn']/totalnnpairs), color='purple')
            #labels=['both cancer','both normal']
            plt.legend(labels)
            plt.ylabel('Fraction of structurally similar pairs')
            #plt.ylabel('Number of cell type pairs')
            plt.xlabel('Genomic location (100kb)')
            ax.set_xlim(0,len(chrdata['cvc']))
            ax.set_ylim(0,1)
            #        plt.title('Cancer v noncancer, chromosome'+str(chrnum))
            #figname = 'testplt_chr'+str(i+1)+'.png'
            figname = '/mnt/disk34/user/nsauerwa/go/outputs/localdiff/figures/cancervnormalconservation_chr'+str(chrnum)+'.eps'
            #figname = '/mnt/disk34/user/nsauerwa/go/outputs/localdiff/figures/cancerhistogram_chr'+str(chrnum)+'.eps'
            plt.savefig(figname,bbox_inches='tight')
            plt.close(fig)
            print 'Wrote cancer v noncancer histogram for chr'+str(chrnum)+' to',figname
            #sys.exit()

    # plot chromosome-level box and whisker plots for cancer-cancer and normal-normal pairs
    totalnnpairs = 28 
    totalcnpairs = 120
    totalccpairs = 105
    #first put all chromosomes into one list for each condition
    allnndata = []
    allccdata = []
    for chrnum,chrdata in cancercolsums.iteritems():
        allnndata.append(chrdata['nvn']/totalnnpairs)
        allccdata.append(chrdata['cvc']/totalccpairs)
    fig = plt.figure()
    ax = plt.gca()
    plt.boxplot(allnndata)
    ax.set_ylim(0,1)
    plt.xlabel('Chromosome')
    plt.ylabel('Fraction of structurally similar pairs')
    figname = '/mnt/disk34/user/nsauerwa/go/outputs/localdiff/figures/chromosomelevel_normalpair_boxplot.eps'
    plt.savefig(figname,bbox_inches='tight')
    plt.close(fig)
    print 'Saved chromosome-level similarity fig to',figname
    fig = plt.figure()
    ax = plt.gca()
    plt.boxplot(allccdata)
    ax.set_ylim(0,1)
    plt.xlabel('Chromosome')
    plt.ylabel('Fraction of structurally similar pairs')
    figname = '/mnt/disk34/user/nsauerwa/go/outputs/localdiff/figures/chromosomelevel_cancerpair_boxplot.eps'
    plt.savefig(figname,bbox_inches='tight')
    plt.close(fig)    


def computeSignificantChromosomes(cancercolsums):

    totalnnpairs = 28
    totalccpairs = 105
    chrlvl_normal = []
    chrlvl_cancer = []
    for chrnum,chrdata in cancercolsums.iteritems():
        chrlvl_normal.append(chrdata['nvn']/totalnnpairs)
        chrlvl_cancer.append(chrdata['cvc']/totalccpairs)
    background_normal = [x for sublist in chrlvl_normal for x in sublist]
    background_cancer = [x for sublist in chrlvl_cancer for x in sublist]
    for chrnum in range(len(chrlvl_normal)):
        print 'Chromosome',chrnum+1
        kstest_normal = scipy.stats.ks_2samp(chrlvl_normal[chrnum], background_normal)
        andtest_normal = scipy.stats.anderson_ksamp([chrlvl_normal[chrnum], background_normal])
        #if kstest_normal[1] < 0.05:
        #    print 'normal-normal pairs KS test:',kstest_normal
        print 'normal:',kstest_normal,andtest_normal
        kstest_cancer = scipy.stats.ks_2samp(chrlvl_cancer[chrnum], background_cancer)
        andtest_cancer = scipy.stats.anderson_ksamp([chrlvl_cancer[chrnum], background_cancer])
        print 'cancer:',kstest_cancer,andtest_normal
        #if kstest_cancer[1] < 0.05:
        #    print 'cancer-cancer pairs KS test:',kstest_cancer

def oneValldata(celltypehists):

    onevalldata = {}
    for celltype,chrdict in celltypehists.iteritems():
        onevalldata[celltype] = [[] for i in range(22)]
        for chrnum,data in chrdict.iteritems():
            onevalldata[celltype][chrnum-1] = [len(sublist) for sublist in data]
    return onevalldata

def calcPercSimilarity(compdata, chrlengths):

    totallen = np.sum(chrlengths)
    totalsim = [] #will be list of [(celltype1, celltype2), sum] values
    celltypelist = []
    totalsimbychr = {} #dict keys will be chrnum, value will be list of same format as totalsim
    for key,data in compdata.iteritems():
        chrnum = key[2]
        celltypelist.extend(key[:2])
        if len(data) > 0:
            #chrsum = np.sum(data[:,1]-data[:,0]+1) # doesn't take into account overlapping intervals - need to be smarter than this
            overlaps = np.zeros(int(max(data[:,1])))
            for row in data:
                overlaps[int(row[0]):int(row[1])] = 1
            chrsum = np.sum(overlaps)
        else:
            chrsum = 0
        celltypes = (key[0], key[1])
        if chrnum in totalsimbychr:
            totalsimbychr[chrnum].append([celltypes, chrsum/chrlengths[chrnum-1]])
        else:
            totalsimbychr[chrnum] = [[celltypes, chrsum/chrlengths[chrnum-1]]]
        idx = -1
        for i,sublist in enumerate(totalsim):
            if sublist[0] == celltypes:
                idx = i
        if idx > -1:
            totalsim[idx][1] += chrsum/totallen
        else:
            totalsim.append([celltypes, chrsum/totallen])
    #loop through totalsim to find pairs to compare to Dixon and Rao
    if False:
        for ctpair in totalsim:
            if ctpair[0] == ('IMR90dixon','hESC'):
                print ctpair
            if ctpair[0] == ('IMR90rao','GM12878'):
                print ctpair
            if ctpair[0][0] == 'GM12878':
                #if ctpair[0][1] == 'IMR90rao':
                #    print ctpair
                if ctpair[0][1] == 'HMEC':
                    print ctpair
                if ctpair[0][1] == 'HUVEC':
                    print ctpair
                if ctpair[0][1] == 'K562rao':
                    print ctpair
                if ctpair[0][1] == 'KBM7':
                    print ctpair
                if ctpair[0][1] == 'NHEK':
                    print ctpair
        sys.exit()

        
        for chrnum,data in totalsimbychr.iteritems():
            data = sorted(data, key=lambda x:x[1], reverse=True)
            totalsimbychr[chrnum] = data
        # calculate average sim per chromosome, highest max, lowest min
        avgperchr = np.zeros(22)
        maxsim = [(),0,0]
        zeropairs = []
        for chrnum,data in totalsimbychr.iteritems():
            for pair in data:
                avgperchr[chrnum-1] += pair[1]/253
                if pair[1] > maxsim[2]:
                    maxsim[0] = pair[0]
                    maxsim[1] = chrnum
                    maxsim[2] = pair[1]
                    #if pair[1] < minsim[2]:
                    #minsim[0] = pair[0]
                    #minsim[1] = chrnum
                    #minsim[2] = pair[1]
                if pair[1] == 0:
                    zeropairs.append([pair[0], chrnum, pair[1]])
        print 'Maximally similar pair:', maxsim
        print 'Number of zero pairs:',len(zeropairs)
        #    print 'Minimally similar pair:', minsim
        #print 'Averages by chromosome:'
        minavg = [0,1]
        maxavg = [0,0]
        for i,avgval in enumerate(avgperchr):
            if avgval > maxavg[1]:
                maxavg = [i+1, avgval]
            if avgval < minavg[1]:
                minavg = [i+1, avgval]
        print 'Max avg on chr'+str(maxavg[0])+' = ',maxavg[1]
        print 'Min avg on chr'+str(minavg[0])+' = ',minavg[1]
        sys.exit()

    # sort by overall similarity
    celltypelist = list(set(celltypelist))
    #print celltypelist
    simmat = np.ones((len(celltypelist),len(celltypelist)))
    totalsim = sorted(totalsim, key=lambda x:x[1], reverse=True)
    for ctpair in totalsim:
        simmat[celltypelist.index(ctpair[0][0]), celltypelist.index(ctpair[0][1])] = ctpair[1]
        simmat[celltypelist.index(ctpair[0][1]), celltypelist.index(ctpair[0][0])] = ctpair[1]
    # sort matrix by cell types with highest average similarity
    simsum = np.sum(simmat, axis=0)
    idx = (-simsum).argsort()
    #print idx
    #sys.exit()
    simmat = simmat[:,idx]
    simmat = simmat[idx,:]


    for i,celltype in enumerate(celltypelist):
        if celltype == 'IMR90rao':
            celltypelist[i] = 'IMR90(R)'
        if celltype == 'K562rao':
            celltypelist[i] = 'K562(R)'
        if celltype == 'IMR90dixon':
            celltypelist[i] = 'IMR90(D)'
        if celltype == 'K562la':
            celltypelist[i] = 'K562(LA)'
    celllist = [celltypelist[i] for i in idx]
    #celllist = celltypelist
    return totalsim,simmat,celllist

def writeToFile(colsums, filename):
    
    with open(filename,'wb') as f:
        fwriter = csv.writer(f, delimiter='\t')
        for row in colsums:
            fwriter.writerow(row)
    print 'Wrote colsums to', filename

def writeSetDifftoNarrowPeak(setdiffs):
    
    filename = 'structsimchanges.narrowPeak'
    name = '.'
    score = 1000
    strand = '.'
    pvalue = -1
    qvalue = -1
    peak = -1
    with open(filename, 'wb') as f:
        fwriter = csv.writer(f, delimiter='\t',quotechar="'")
        fwriter.writerow(['track type=narrowPeak', 'db=hg19', 'name=\"structsimPks\"', 'description=\"Changes in sets of pairwise structural similarity\"'])
        for i,data in enumerate(setdiffs):
            chrlabel = 'chr'+str(i+1)
            for genloc,sigval in enumerate(data):
                if sigval > 0:
                    chromstart = (genloc*100000) + 100000/2
                    chromend = ((genloc+1)*100000) + 100000/2
                    newpk = [chrlabel, chromstart, chromend, name, score, strand, sigval, pvalue, qvalue, peak]
                    fwriter.writerow(newpk)
    print 'Wrote narrowPeak file to', filename


def writeSetDifftoBedGraph(setdiffs,filename):

    with open(filename, 'wb') as f:
        fwriter = csv.writer(f, delimiter='\t')
        fwriter.writerow(['track type=bedGraph'])
        for i,data in enumerate(setdiffs):
            chrlabel = 'chr'+str(i+1)
            for genloc,dataval in enumerate(data):
                if dataval > 0:
                    chromstart = genloc*100000 + 100000/2
                    chromend = (genloc+1)*100000 + 100000/2
                    newline = [chrlabel, chromstart, chromend, dataval]
                    fwriter.writerow(newline)
    print 'Wrote bedGraph file to', filename


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
    print 'Heatmap saved as',filename

def plotIntervalResults(compdata,tadlists,chrlengths,key,newdata):

    chrnum = key[2]
    matsize = int(np.ceil(chrlengths[chrnum]/100000)) #812 #90354753/res 81195210/res
    for interval in tadlists[0]:
        plt.plot((interval[0]-1,interval[1]), (matsize-interval[0] + 1, matsize-interval[0]+1),'g')
        plt.plot((interval[1],interval[1]),(matsize-interval[0]+1,matsize-interval[1]),'g')
        dom1Artist = plt.Line2D((0,1),(0,0), color='green', linestyle='solid')
    for interval in tadlists[1]:
        plt.plot((interval[0]-1,interval[1]),(matsize-interval[1],matsize-interval[1]),'b')
        plt.plot((interval[0]-1,interval[0]-1),(matsize-interval[0]+1,matsize-interval[1]),'b')
        dom2Artist = plt.Line2D((0,1),(0,0), color='blue', linestyle='solid')
    #for interval in compdata[key]: #compdata[('T47D','HUVEC',16)]:
    #    plt.plot((interval[0]-1,interval[1]),(matsize-interval[1],matsize-interval[1]),'r')
    #    plt.plot((interval[0]-1,interval[0]-1),(matsize-interval[0]+1,matsize-interval[1]),'r')
    #    dom3Artist = plt.Line2D((0,1),(0,0), color='red', linestyle='solid')
    ax = plt.gca()
    if len(newdata)>0:
        ax.scatter(newdata[:,1], matsize-newdata[:,0], s=15, color='red', marker='*')
    figname = '/mnt/disk34/user/nsauerwa/go/outputs/localdiff/figures/methodoverview_allbdypts.eps'
#    ax.set_facecolor('white')
    ax.set_aspect('equal', adjustable='box')
#    plt.tick_params(axis='x',labelbottom='off')
    plt.tick_params(axis='y',labelleft='off')
#    figname = '/mnt/disk34/user/nsauerwa/go/outputs/localdiff/figures/methodoverview_KBM7tads.eps'
#   figname = '/mnt/disk34/user/nsauerwa/go/outputs/localdiff/figures/methodoutput_example_'+key[0]+'_'+key[1]+'_chr'+str(chrnum)+'_justTADs.eps'
    plt.savefig(figname,bbox_inches='tight')
    print 'saved figure as',figname
#    ax.scatter(localmins[:,1], matsize-localmins[:,0], s=25, marker='*', color='red')

def readTADfile(filename, res):
    # written for Armatus files for now
    tadbdys = []
    #print filename
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

    newdata = []
    with open(filename,'rb') as f:
        freader = csv.reader(f,delimiter='\t')
        freader.next()
        for line in freader:
            newdata.append([int(line[0]), int(line[1]), float(line[2])])
    return np.array(newdata)


def main(fileseed):
    
    filelist = glob.glob(fileseed)

    res = 100000
    chrlengths = [248956422, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566] # hg19 chr lengths from https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
    for i,chrlen in enumerate(chrlengths):
        chrlengths[i] = int(np.ceil(float(chrlen)/res))

    compdata = readLocalDiffFiles(filelist)

    #    tadfiles = ['/mnt/disk17/armatusresults/A549/100kb/EncodeChr18_A549_combo_IC_100kb_gmax1.0.gamma.1.0.txt','/mnt/disk17/armatusresults/KBM7/100kb/RaoChr18_KBM7_IC_100kb_gmax1.0.gamma.1.0.txt']
    #['/mnt/disk17/armatusresults/A549/100kb/EncodeChr21_A549_combo_IC_100kb_gmax1.0.gamma.0.9.0.txt','/mnt/disk17/armatusresults/hESC/hicpro/100kb/DivxonChr21_hESC_combined_100kb_gmax1.0.gamma.0.9.0.txt']
    #['/mnt/disk17/armatusresults/IMR90/100kb/RaoChr17_IMR90_IC_100kb_gmax1.0.gamma.0.5.0.txt', '/mnt/disk17/armatusresults/HMEC/100kb/RaoChr17_HMEC_IC_100kb_gmax1.0.gamma.0.5.0.txt']
    tadfiles = ['mini_example_a549.txt','mini_example_kbm7.txt']
    tadlists = []
    for filename in tadfiles:
        tadlists.append(readTADfile(filename,100000))
    allvis = readIntermedFiles('mini_example_bdyvis.txt')
    plotIntervalResults(compdata,tadlists,[165],('A549','KBM7',0),allvis)
    sys.exit()

#    computeBasicStats(compdata)
   # percsim,simmat,celltypelist = calcPercSimilarity(compdata, chrlengths)
    #percsim heatmaps
  #  plotHeatMap(simmat,celltypelist,'outputs/localdiff/figures/global_percsim_heatmap.eps')
 #   sys.exit()


#    for i,row in enumerate(percsim):
#        if 'K562la' == row[0][1] or 'GM06990' == row[0][1] or 'GM06990' == row[0][0]:
#            print row, 'rank', i+1
        #if i < 10:
        #    print row, 'rank', i+1, 'out of', len(percsim)
        #if i == 11: print('')
        #if len(percsim)-i < 11:
        #    print row, 'rank', i+1, 'out of', len(percsim)
    #print('')
    #for i,row in enumerate(percsim):   
    #    if row[0] == ('K562rao','K562la') or row[0] == ('IMR90rao','IMR90dixon') or row[0] == ('K562rao','KBM7') or row[0] == ('GM12878','GM06990'):
    #        print row, 'rank',i+1,'out of',len(percsim)
 #   colsums,setdiffs, cancercolsums = intervalHistogram(compdata,chrlengths,res)
#    computeSignificantChromosomes(cancercolsums)
#    genelocs = readGeneLocFile()
#    print 'Percent of genes with higher similarity in normal-normal than cancer-cancer pairs:',countNumGenes(genelocs,cancercolsums)
#    plotCancerGeneData(cancercolsums)
   #writeSetDifftoBedGraph(colsums,'similaritybyloc.bedGraph')
    #writeSetDifftoNarrowPeak(setdiffs)
    #writeSetDifftoBedGraph(setdiffs, 'structsimchanges.bedGraph')
    #celltypehists = histogramByCellType(compdata,chrlengths)
    #onevalldata = oneValldata(celltypehists)
    #for celltype in ['GM12878','HUVEC','K562rao','NHEK']:
    #    filename = 'onevall_'+celltype+'.bedGraph'
    #    writeSetDifftoBedGraph(onevalldata[celltype],filename)
    # want: bar plots of colsums and setdiffs for each chromosome
    #       color-coded bar plots for each celltype hist
#    plotHistograms(colsums,setdiffs,cancercolsums,[])#celltypehists)

#    writeToFile(colsums,'lochistograms_allchr.txt')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='fileseed of files to analyze, with wildcard characters')
    
    args = parser.parse_args()
    main(args.i)
