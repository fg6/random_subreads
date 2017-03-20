#!/usr/bin/env python
from __future__ import print_function  
import sys, getopt
from Bio import SeqIO
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os.path

plot=0

totall=0
readall=[]
subreads=[]
meantot=0
mean=0
std=1000
warning=0.6 # about min ratio of coverage needed/initial coverage (avg prob for a read to get chosen if inside the gaussian peak)

def read(inputfile,faqtype,cov,refsize):
  records = (r for r in SeqIO.parse(inputfile, faqtype)) 
  global totall
  global readall
  global meantot
  global inin
  global bins
  global maxN
  global totN
  global weigth
  weigth=1.
  totN=0
  loctot=0
  for r in records:
     totall+=len(r.seq) 
     readall.append(len(r.seq))
     if mean > 0 and len(r.seq)< mean+2*std and len(r.seq)> mean-2*std: 
        meantot+= len(r.seq)
        loctot+= 1
     totN+=1

  if mean == 0:
    pcov=totall/refsize
    pmean=totall/totN
    thisTot=totN
  else:
    pcov=meantot/refsize
    pmean=meantot/loctot
    thisTot=loctot

  # delete this?? 
  [inin,bins]=np.histogram(readall, 20)
  maxN=max(inin)

  #neededN=thisTot*(cov/pcov)*(mean/pmean)/warning
  #neededN=thisTot*(cov/mean)*(pmean/pcov)
  #neededN=(warning+0.05)*thisTot

  weigth=1  #(neededN/totN)+warning # slightly reduce weigth to make sure you almost reach the end of the reads, to account for non-randomness distri in the read order
  #print(" need ",neededN," reads, percentage of total is",neededN/totN,"weigth=",weigth," tot-needed:",totN-neededN)
  #weigth=0.8 

  if(cov > pcov): 
    print ("\n Sorry your original file does not have enough coverage! \n  Total coverage in the selected area=%0.1fX, desired coverage=%0.1fX" % (pcov, cov))
    return 1
  elif cov/pcov >= warning:
    print ("\n Probably not enough coverage in the selected area: initial coverage for selected area=%0.1fX, can extract about %0.1fX-%0.1fX" % (pcov,(warning)*pcov,(warning+0.05)*pcov))
    return 0
  else: 
    print ("\n Initial coverage for selected area=%0.1fX" % pcov)
    return 0


def plotnow():
    plt.rc('xtick', labelsize=8)  # x axis labels
    plt.rc('ytick', labelsize=8)   # y axis labels
    font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 10}
    plt.rc('font', **font)
     
    fig = plt.figure()
    plt.hist(readall,100,histtype='step',color='b',label="All")
    plt.hist(subreads,100,histtype='step',color='r',label="Sub")


     
    plt.tight_layout()
    plt.show()

    pdfname='subsample.pdf'
    myfig = PdfPages(pdfname)
    myfig.savefig(fig)
    myfig.close()

def gaussian(x, mu, sig):
    return 1.*np.exp(-np.power((x - mu)/sig, 2.)/2)

def findN(length):
   global inin
   global bins


   mybin=-1
   for b in range(1,len(bins)):
     if length > bins[b-1] and length < bins[b]: mybin=b-1
     if mybin != -1: break
   return inin[mybin]


def readnrandom(inputfile,faqtype,ofaqtype,cov,refsize): 

  global subreads


  records = list(SeqIO.parse(inputfile, faqtype))
  readlist=[]
 
  
  if(mean==0):
    mytotal=totall
  else:
    mytotal=meantot

  mysum=0.
  allmean=0
  i=-1
  ii=0

  # choose a random, non-repetitive order:
  randord=[]
  listord=range(0,len(records))
#  listord=range(0,10)
  nn=len(listord)
  
  for ir in range(0,nn):
    thisr=random.choice(listord)      #random.randint(0,len(records))
    randord.append(thisr)
    listord.remove(thisr)


  for ir in randord:
   
   thisl=len(records[ir].seq)
   allmean+=thisl   #len(r.seq)
   i+=1
   if i<10: print(i,ir,thisl)

   if mysum/refsize > cov:
    if(1<2): print (cov, "% hit!",i,totN)
    break
   else: 
    test=0

    #thisN=findN(thisl)   #len(r.seq))


    if(mean==0):
      value=1-cov/(mytotal/refsize)
      ran=random.randint(0,100)
      if ran > value: test=1
    else:
       model=gaussian(thisl,mean,std) #+gaussian(thisl,mean-3000,std) 
       sigma_up=model*weigth  
       ran=random.uniform(0,1)

       if ran <= sigma_up: test=1

    if test:  
     readlist.append(records[ir])
     mysum+=thisl   
     subreads.append(thisl) 
     ii+=1


  print("Min subread:",min(subreads), " tot subreads: ", ii)

  if(float(mysum)/refsize < cov): 
    print (" Finished: could not extract enough coverage, extracted coverage=%0.2fX" % (float(mysum)/refsize))
    return 0 
  else:
    print ("  Final number of bases=%0.0f,  Coverage=%0.1fX" % (mysum,float(mysum)/refsize))
    print ("  Initial read length mean:%0.1f  Final read length mean=%0.1f" %(allmean/i,mysum/ii))

  ofile="subreads_%0.0fX.%s" % (cov,ofaqtype)
  if mean > 0: ofile="subreads_%0.0fX_mean%d.%s" % (cov,mysum/ii,ofaqtype)
  output_handle = open(ofile,"w")  
  SeqIO.write(readlist, output_handle,ofaqtype )
  output_handle.close()
  print (" Subsample written in",ofile,"\n")


def usage():
  print (' Missing required variables ! Usage: ')
  print (' subrand.py -i <inputfile> -t <faqtype> -c <coverage> -r <refsize> -m <mean> -s <std>')
  print ('  input file: fasta or fastq file of reads from which extract a subsample')
  print ('  faqtype: format of output file, fasta or fastq')
  print ('  coverage: desired coverage for subsample')
  print ('  refsize: reference size or expected genome size needed to calculate coverage, in Mb')
  print ('  mean,std: mean and standard deviation of read lengths desired for subsample. If not defined, a uniform random selection of reads is performed')
  

def main(argv):
   inputfile = ''
   ofaqtype = 'fasta'
   cov = 0
   refsize = 0
   global mean
   global std
   try:
      (opts, args) = getopt.getopt(argv,"i:t:c:r:m:s:",["ifile=","faqtype=","cov=","refsize=","mean=","std="])   
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('subrand.py -i <inputfile> -t <faqtype>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-t", "--faqtype"):
         ofaqtype = arg
      elif opt in ("-c", "--cov"):
         cov = float(arg)
      elif opt in ("-r", "--refsize"):
         refsize = float(arg)*1000000
      elif opt in ("-m", "--mean"):
         mean = float(arg)
      elif opt in ("-s", "--std"):
         std = float(arg)

   if len(inputfile) == 0 or cov ==0 or refsize ==0:
      usage()
      sys.exit(2)
   if not os.path.exists(inputfile): 
      print("Sorry, file ", inputfile, "does not exists")
      sys.exit(2)


   type=inputfile.split(".")[-1]
   if type == "fq":
    faqtype="fastq"
   elif type == "fa":
    faqtype="fasta"
   else:
    faqtype=type

      
   print ("\n Getting a",cov,"X subsample of reads from", inputfile)
   print ("   (Coverage calculated for a reference size of %0.1f Mb)" % (refsize/1000000))
  
   
   ok=read(inputfile,faqtype,cov,refsize)
   if(ok==0): 
     readnrandom(inputfile,faqtype,ofaqtype,cov,refsize)    
     if(plot):
      plotnow()
   

if __name__ == "__main__":
   main(sys.argv[1:])
