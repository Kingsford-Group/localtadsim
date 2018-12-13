# turn 3 column sparse Hi-C matrix format into n x (n+3) matrix for input into HiCRep
# usage: R run_hicrep.r filename1 celltype1 filename2 celltype2 resolution
# filename1 and filename2 should be partial paths to locations of Hi-C data, ie "GM12878/matrix/combined/iced/100000/combined_100000_iced_chr"

#options(scipen=999)

args = commandArgs(trailingOnly=TRUE)
#print(args)
filename1 = args[1]
label1 = args[2]
filename2 = args[3]
label2 = args[4]	 
res = as.numeric(args[5])

# create data frame to store chr values in, then take mean at the end
chrrepvals <- vector(mode="numeric",length=22)
for (chr in 1:22) {
    chrnum <- paste("chr",as.character(chr),sep="")
    #print(chrnum)

    fullfilename1 <- paste(paste(filename1,chrnum,sep=""),".matrix",sep="")
    fullfilename2 <- paste(paste(filename2,chrnum,sep=""),".matrix",sep="")

    rawdata1 = read.table(file = fullfilename1, sep = '\t', header = FALSE)
    rawdata1[,1:2] <- rawdata1[,1:2]/res

    rawdata2 = read.table(file = fullfilename2, sep = '\t', header = FALSE)
    rawdata2[,1:2] <- rawdata2[,1:2]/res

    nbins1 <- max(rawdata1$V2)
    nbins2 <- max(rawdata2$V2)
    nbins <- max(c(nbins1,nbins2))

    fullmatrix <- matrix(data=0, nrow=nbins, ncol=nbins+2)
    for (row in 1:nbins) {
    	fullmatrix[row,1] = (row-1)*res
    	fullmatrix[row,2] = row*res
    }   

    for (row in 1:nrow(rawdata1)) {
    	fullmatrix[rawdata1[row,1],rawdata1[row,2]+2] = rawdata1[row,3]
    }
    fullmatrix = data.frame(fullmatrix)

    chrlist <- rep(chrnum,nbins)
    hicmat1 <- data.frame(V1=chrlist,fullmatrix)

    fullmatrix <- matrix(data=0, nrow=nbins, ncol=nbins+2)
    #fullmatrix[,1] = chrnum
    for (row in 1:nbins) {
    	fullmatrix[row,1] = (row-1)*res
    	fullmatrix[row,2] = row*res
    }   

    for (row in 1:nrow(rawdata2)) {
    	fullmatrix[rawdata2[row,1],rawdata2[row,2]+2] = rawdata2[row,3]
    }
    fullmatrix = data.frame(fullmatrix)
    hicmat2 <- data.frame(V1=chrlist,fullmatrix)

    library(hicrep)

    h_hat <- htrain(hicmat1,hicmat2,100000,5000000,0:3)
    #print(h_hat)
    pre_hic <- prep(hicmat1,hicmat2,100000,h_hat,5000000)

    SCC.out = get.scc(pre_hic,100000,5000000)

    chrrepvals[chr] = SCC.out$scc
}
totalrepval <- mean(chrrepvals)

cat(sprintf("%s %s %f \n",label1,label2,totalrepval))
