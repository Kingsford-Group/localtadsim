import sys
import csv
import numpy as np

hicfilename = sys.argv[1]

#read HiCfile (blah_raw.matrix) into n x 3 matrix
def readHiCProFile(filename):
    with open(filename,'rb') as f:
        freader = csv.reader(f,delimiter='\t')
        rawdata = []
        for row in freader:
            rawdata.append([int(row[0]), int(row[1]), int(row[2])])
    rawdata = np.array(rawdata)
    return rawdata

def sumContactsByBin(rawdata):
    veclength = np.max(rawdata[:,:2])
    bincontacts = np.zeros(veclength)
    for row in rawdata:
        bincontacts[row[0]-1] += row[2]
        bincontacts[row[1]-1] += row[2]
    return bincontacts

def checkPercAbove1000(bincontacts):
    numAbove1000 = (bincontacts > 1000).sum()
    percAbove1000 = float(numAbove1000)/len(bincontacts)
    return percAbove1000

rawdata = readHiCProFile(hicfilename)
print 'read data...'
bincontacts = sumContactsByBin(rawdata)
print 'summed contacts...'
percAbove1000 = checkPercAbove1000(bincontacts)
print 'calculated percentage over 1000...'
if percAbove1000 >= 0.8:
    print 'Passed:', hicfilename
else:
    print 'Failed:', hicfilename
