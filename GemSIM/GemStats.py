#!/usr/bin/python

# Copyright (C) 2011, Kerensa McElroy.
# kerensa@unsw.edu.au

# This file is part of the sequence simulator GemSIM. 
# It is used to calculate some statistics from error
# model files created with GemErr.py

# GemSIM is free software; it may be redistributed and 
# modified under the terms of the GNU General Public 
# License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option)
# any later version.

# GemSIM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, without even the implied 
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more 
# details.

# You should have recieved a copy of the GNU General Public
# License along with GemSIM. If not, see 
# http://www.gnu.org/licenses/. 


import sys
import getopt
import numpy as np
import cPickle
import gzip

def processGZIP(gzipFile,paired):
    f=gzip.open(gzipFile)
    if paired:
       readLen=cPickle.load(f)
       mx1=cPickle.load(f)
       mx2=cPickle.load(f)
       insD1=cPickle.load(f)
       insD2=cPickle.load(f)
       delD1=cPickle.load(f)
       delD2=cPickle.load(f)
       intD=cPickle.load(f)
       gQualL=cPickle.load(f)
       bQualL=cPickle.load(f)
       iQualL=cPickle.load(f)
       mates=cPickle.load(f)
       rds=cPickle.load(f)
       rdLenD=cPickle.load(f)
       f.close()
       return mx1,mx2,insD1,insD2,delD1,delD2,intD,gQualL,bQualL,iQualL,mates,rds,rdLenD 
    else:
       readLen=cPickle.load(f)
       mx=cPickle.load(f)
       insD=cPickle.load(f)
       delD=cPickle.load(f)
       gQualL=cPickle.load(f)
       bQualL=cPickle.load(f)
       iQualL=cPickle.load(f)
       readCount=cPickle.load(f)
       rdLenD=cPickle.load(f)
       f.close()
       return mx,insD,delD,gQualL,bQualL,iQualL,readCount,rdLenD

def parseMX(mx,insD, delD):
    tot=np.apply_over_axes(np.sum,mx,[0,1,2,3,4,5])[0][0][0][0][0][0][5]
    muts=np.sum(mx)-tot    
    bases=np.apply_over_axes(np.sum,mx,[0,2,3,4,5])
    aTot=bases[0][0][0][0][0][0][5]
    tTot=bases[0][1][0][0][0][0][5]
    gTot=bases[0][2][0][0][0][0][5]
    cTot=bases[0][3][0][0][0][0][5]
    nTot=bases[0][4][0][0][0][0][5] 
    aMut=np.sum(bases[0][0][0][0][0][0][:5])
    tMut=np.sum(bases[0][1][0][0][0][0][:5])
    gMut=np.sum(bases[0][2][0][0][0][0][:5])
    cMut=np.sum(bases[0][3][0][0][0][0][:5])
    posM=np.add.reduce(np.apply_over_axes(np.sum, mx, [1,2,3,4,5])[:,0,0,0,0,0,:5],axis=1)
    posT=np.apply_over_axes(np.sum, mx, [1,2,3,4,5])[:,0,0,0,0,0,5]
    rPos=posM.astype(float)/posT.astype(float)
    mx=np.add.reduce(mx,axis=0)
    mxTot=np.repeat(mx[:,:,:,:,:,5],5)
    mxTot=mxTot.reshape([5,5,5,5,5,5])
    mx=mx[:,:,:,:,:,:5]
    mx=mx.astype(float)/mxTot.astype(float) 
    mutants=np.ndenumerate(mx)
    list=[]
    for pos,val in mutants:
        if 4 not in pos:
            if str(val)!='NaN' and str(val)!='nan':
                if str(val)=='inf':
                    val=1.0
                list.append(val)
    avg=np.mean(list)
    stDv=np.std(list)
    highErrs=[]
    mutants=np.ndenumerate(mx)
    P=float(muts)/float(tot)
    for pos, val in mutants:
        if 4 not in pos:
            if val>(avg+stDv*2) and val>P:
                if str(val)!='NaN' and str(val)!='nan':
                    if val=='inf':
                        val=1.0
                    highErrs.append((val,pos))
    A=float(aMut)/float(aTot)
    T=float(tMut)/float(tTot) 
    G=float(gMut)/float(gTot)
    C=float(cMut)/float(cTot)
    ikeys=insD.keys()
    iDict={}
    for k in ikeys:
        seq=k.partition('.')[2]
        if seq in iDict:
            for s in insD[k].keys():
                if s in iDict[seq]:
                    iDict[seq][s]+=insD[k][s] 
                else:
                    iDict[seq][s]=insD[k][s]
        else:
            iDict[seq]=insD[k]
    insD=iDict
    ikeys=insD.keys()
    kList=[]
    kDict={} 
    for key in ikeys:
        ksum=0
        kL=[]
        ke = insD[key].keys()
        for k in ke:
            kList.append(insD[key][k])
            ksum+=insD[key][k]
            kL.append(k)
        if ksum in kDict:
            kDict[ksum].append((key, kL))
        else:
            kDict[ksum]=[(key,kL)] 
    iSums=kDict.keys()
    iSums.sort()
    iSums=iSums[::-1]
    i=0
    totI=sum(kList)
    iP=float(totI)/float(tot)
    dkeys=delD.keys()
    dList=[]
    dDict={}
    for key in dkeys:
        dsum=0
        dL=[]
        for k in range(len(delD[key])):
            dList.append(delD[key][k])
            dsum+=delD[key][k]
            dL.append(k)
        if dsum in dDict:
            dDict[dsum].append((key,dL))
        else:
            dDict[dsum]=[(key,dL)]
    dSums=dDict.keys()
    dSums.sort()
    dSums=dSums[::-1]
    totD=sum(dList)
    dP=float(totD)/float(tot)
    return rPos,tot, P, A, T, G, C, highErrs, avg, stDv, iP,iSums,kDict,dP,dSums,dDict

def prtPrint(name,read, bP,Tot, P, A, T, G, C, highErrs,avg,stDv,totI,iSums,iDict,totD,dSums,dDict):
    inds={0:'A',1:'T',2:'G',3:'C',4:'N'}
    g=open(name+'_'+read+'_stats.txt','w')
    g.write('\nError analysis for '+read+'\n\n')
    g.write('Overall error rate:\t'+str(P)+'\n\n')
    g.write('Bases covered by '+read+':\n')
    g.write('\n'+str(Tot)+'\n\n') 
    g.write('Error rate along length of read:\n')
    x=1
    for i in bP.flatten():
        g.write(str(x)+'\t'+str(i)+'\n')
        x+=1
    g.write('\n\nError rate for nucleotides:\n')
    g.write('A:\t'+str(A))
    g.write('\nT:\t'+str(T))
    g.write('\nG:\t'+str(G))
    g.write('\nC:\t'+str(C))
    g.write('\n\nSequences with unusually high error rates:\n')
    highErrs.sort()
    for i in highErrs:
        g.write(inds[i[1][3]]+inds[i[1][2]]+inds[i[1][1]]+' '+inds[i[1][0]]+' '+inds[i[1][4]]+' -> '+inds[i[1][3]]+inds[i[1][2]]+inds[i[1][1]]+' '+inds[i[1][5]]+' '+inds[i[1][4]]+':\t'+str(i[0])+'\n')
    g.write("\n\nProbability of an insertion: "+str(totI))
    g.write("\n\nProbability of a deletion: "+str(totD)+'\n')

def usage():
    print '\n\n########################################################################'
    print '# GemSIM - Generic Error Model based SIMulator of N.G. sequencing data #'
    print '########################################################################\n'
    print '\nGemStats.py:\n'
    print 'Takes error model files produce by GemErr.py, and generates statistics'
    print 'for a particular error model. Output saved as .txt file.'
    print '\nOptions:'
    print '      -h prints these instructions.'
    print '      -m error model file *_single.gzip or *_paired.gzip.'
    print '      -p use if model is for paired end reads.'
    print '      -n prefix for output files.\n\n'

def main(argv):
    models=''
    paired=False
    name=''
    try:
        opts, args = getopt.getopt(argv, "hm:pn:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt =='-h':
            usage()
            sys.exit()
        elif opt =='-m':
            models=arg
        elif opt =='-p':
            paired=True
        elif opt =='-n':
            name=arg
    if models=='' or name=='':
        usage()
        sys.exit()
    if paired:
        mx1,mx2,insD1,insD2,delD1,delD2,inD,gQualL,bQualL,iQualL,mates,rds,rdLenD=processGZIP(models,True)
        pos,tot, P, A, T, G, C, highErrs, avg, std, totI,iSums,iDict,totD,dSums,dDict=parseMX(mx1,insD1,delD1)
        prtPrint(name, 'read_1', pos, tot, P, A, T, G, C, highErrs,avg,std,totI,iSums,iDict,totD,dSums,dDict)
        pos,tot, P, A, T, G, C, highErrs, avg, std,totI,iSums,iDict,totD,dSums,dDict=parseMX(mx2,insD2,delD2)
        prtPrint(name, 'read_2', pos, tot, P, A, T, G, C, highErrs,avg,std,totI,iSums,iDict,totD,dSums,dDict)
    else:
        mx, insD, delD, gQualL, bQualL, iQualL,rds,rdLenD=processGZIP(models,False)
        pos,tot,P, A, T, G, C, highErrs, avg, std,totI,iSums,iDict,totD, dSums,dDict=parseMX(mx,insD,delD)
        prtPrint(name,'single_read', pos, tot,P, A, T, G, C, highErrs,avg,std,totI,iSums,iDict,totD,dSums,dDict) 

if __name__=="__main__":
    main(sys.argv[1:])


