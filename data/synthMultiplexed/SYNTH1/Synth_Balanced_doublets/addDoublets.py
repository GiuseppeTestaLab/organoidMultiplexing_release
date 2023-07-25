#!/usr/bin/python3

import os
import pysam
import pandas as pd
from multiprocessing import Pool
import subprocess
import random
import shutil
import numpy as np



barcodeList=pd.read_csv("/SYNTH1/Synth_Balanced_noDoublets/Barcodes-IDs.tsv", header= None , sep = "\t", names = ["Barcode","ID"])
targetDblRatio=0.1
inputBam="/SYNTH1/Synth_Balanced_noDoublets/SilicoMultiplexed.bam"
outfile="/SYNTH1/Synth_Balanced_doublets/SilicoMultiplexed.dblAdded.bam"
nThreads=10
WD=os.getcwd()+'/'
TEMPDIR=WD+'tempdir/'
contigs=list(pysam.AlignmentFile(inputBam, "rb").header.references)



try:
    os.mkdir(TEMPDIR)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass

chunkList=["Chunk_"+ contig + ".bam" for contig in contigs]
with open(TEMPDIR+"ChunkList.txt", "w") as chunks:
    for chunk in chunkList:
        chunks.write(TEMPDIR+chunk+"\n")

sampleSize =int((len(barcodeList) * targetDblRatio * 2/(1 + targetDblRatio))-((len(barcodeList) * targetDblRatio * 2/(1 + targetDblRatio))%2))
print("After rounding actual dbl ratio will be "+str(round((sampleSize/2)/(len(barcodeList)-(sampleSize/2)), 5)))
print("Beads are being merged .  .  .")

def replacingLister(barcodeList, targetDblRatio, sampleSize):
    sample=random.sample(list(barcodeList["Barcode"]), int(sampleSize))
    dblAcceptor=sample[0::2]
    dblDonor=sample[1::2]
    doubletsMap=pd.DataFrame({'dblAcceptor': dblAcceptor, 'dblDonor': dblDonor})
    doubletsMap.to_csv(WD+"doubletsMap.tsv", sep="\t", header=True, line_terminator="\n", index = False)
    conversionDict = pd.Series(doubletsMap.dblDonor.values, index=doubletsMap.dblAcceptor).to_dict()
    doubletsMap["DBLstatus"] = "Doublet"
    #Produce df containing new mixed genotypes per barcode
    barcodeListMod=barcodeList.rename(columns={'Barcode':'dblAcceptor'}, inplace=False)
    doubletsMap=doubletsMap.merge(barcodeListMod, on = "dblAcceptor")
    barcodeListMod=barcodeListMod.rename(columns={'dblAcceptor':'dblDonor', 'ID':'ID2'}, inplace=False)
    doubletsMap=doubletsMap.merge(barcodeListMod, on = "dblDonor")
    doubletsMap["ID"] = doubletsMap["ID"]+ "," +doubletsMap["ID2"]
    doubletsMap=doubletsMap.drop(columns = ["dblAcceptor" ,    "ID2" ])
    doubletsMap=doubletsMap.rename(columns={'dblDonor':'Barcode'}, inplace=False)
    return conversionDict, dblAcceptor, doubletsMap


def DBLfactory(INBAM, chr):
    bamFile = bamFile=pysam.AlignmentFile(INBAM, "rb")
    out= TEMPDIR+"Chunk_" + str(chr)+".bam"
    outBam = pysam.AlignmentFile(out,"wb", template=bamFile)
    for read in bamFile.fetch(str(chr)):
        try:
            CBtag = read.get_tag('CB')
        except:
            outBam.write(read)
        if CBtag in list(conversionDict.keys()):
            read.set_tag(tag='CB', value=conversionDict[CBtag] , value_type="Z")
            read.set_tag(tag='CR', value=str(conversionDict[CBtag]).split("-")[0], value_type="Z")
            outBam.write(read)
        else:
            outBam.write(read)
    bamFile.close()
    outBam.close()

conversionDict, dblAcceptor, doubletsMap=replacingLister(barcodeList, targetDblRatio, sampleSize)


temp=barcodeList[(~ barcodeList.Barcode.isin(conversionDict.keys())) & (~ barcodeList.Barcode.isin(conversionDict.values()))]
temp.insert(1, "DBLstatus", "Singlet", True)
doubletsMap=pd.concat([doubletsMap, temp], axis = 0)






pool=Pool(nThreads)
for contig in contigs:
    pool.apply_async(DBLfactory, (inputBam, contig))

pool.close()
pool.join()



print("Merging chunks .  .  .")
subprocess.run(["samtools","cat","-o",TEMPDIR+"Merged.temp.bam","-b",TEMPDIR+"ChunkList.txt"])
print("Sorting outbam .  .  .")
subprocess.run(["samtools","sort","-o",outfile,"-@"+str(nThreads),TEMPDIR+"Merged.temp.bam"])
print("Indexing outbam .  .  .")
subprocess.run(["samtools","index",outfile])
print("Deleting temps .  .  .")
shutil.rmtree(TEMPDIR)
print("Writing new barcodes file in "+WD+"Barcodes.new.tsv .  .  .")
with open(WD+"Barcodes.new.tsv", "w") as barcodes:
    for newbarcode in list(np.setdiff1d(barcodeList.Barcode, dblAcceptor)):
        barcodes.write(newbarcode+"\n")

print("Writing new barcodes file in "+WD+"doubletsMap.tsv .  .  .")
doubletsMap.to_csv(WD+"doubletsMap.tsv", sep = "\t", header = True, index = False)
