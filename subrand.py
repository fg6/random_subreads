#!/usr/bin/env python

import sys, getopt
from Bio import SeqIO
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


plot=1

totall=0
readall=[]
subreads=[]
meantot=0
mean=0
std=1000
warning=0.6

def read(inputfile,faqtype,cov,refsize):
  records = (r for r in SeqIO.parse(inputfile, faqtype))
  global totall
  global readall
  global meantot


  for r in records:
     totall+=len(r.seq)
     readall.append(len(r.seq))
     if mean > 0 and len(r.seq)< mean+2*std and len(r.seq)> mean-2*std: meantot+= len(r.seq)

  if mean == 0:
    pcov=totall/refsize
  else:
    pcov=meantot/refsize

  if(cov > pcov):
    print "\n Sorry your original file does not have enough coverage! \n  Total coverage in the selected area=%0.1fX, desired coverage=%0.1fX" % (pcov, cov)
    return 1
  elif cov/pcov >= warning:
    print "\n Probably not enough coverage in the selected area: initial coverage for selected area=%0.1fX, can extract about %0.1fX-%0.1fX" % (pcov,(warning)*pcov,(warning+0.05)*pcov)
    return 0
  else:
    print "\n Initial coverage for selected area=%0.1fX" % pcov
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

def readnrandom(inputfile,faqtype,ofaqtype,cov,refsize):

  global subreads
  records = (r for r in SeqIO.parse(inputfile, faqtype))
  readlist=[]

  allhist, edges =  np.histogram(readall, bins=100)

  if(mean==0):
    mytotal=totall
  else:
    mytotal=meantot

  mysum=0.
  allmean=0
  i=-1
  ii=0
  for r in records:
   allmean+=len(r.seq)
   i+=1
   if mysum/refsize > cov:
    if(1<0): print cov, "% hit!"
   else:
    test=0
    if(mean==0):
      value=1-cov/(mytotal/refsize)
      ran=random.randint(0,100)
      if ran > value: test=1
    else:
       #value=cov/(mytotal/refsize)
       model=gaussian(len(r.seq),mean,std)
       sigma_up=model
       ran=random.uniform(0,1)
       #if i<20 : print len(r.seq), peak,"thisvalue:", value, "ran:",ran, "sigma:",sigma_up," model:",model
       if ran <= sigma_up: test=1

    if test:
     readlist.append(r)
     mysum+=len(r.seq)
     subreads.append(len(r.seq))
     ii+=1

  #print "Total=",i, "Selected=",ii,"Ratio=",float(ii)/float(i)
  if(float(mysum)/refsize < cov):
    print " Finished: could not extract enough coverage, extracted coverage=%0.1fX" % (float(mysum)/refsize)
    return 0
  else:
    print "  Final number of bases=%0.0f,  Coverage=%0.1fX" % (mysum,float(mysum)/refsize)
    print "  Initial read length mean:%0.1f  Final read length mean=%0.1f" %(allmean/i,mysum/ii)

  ofile="subreads_%sX.%s" % (cov,ofaqtype)
  output_handle = open(ofile,"w")
  SeqIO.write(readlist, output_handle,ofaqtype )
  output_handle.close()
  print " Subsample written in",ofile,"\n"


def main(argv):
   inputfile = ''
   ofaqtype = 'fasta'
   cov = 0
   refsize = 0
   global mean
   global std
   try:
      opts, args = getopt.getopt(argv,"i:t:c:r:m:s:",["ifile=","faqtype=","cov=","refsize=","mean=","std="])
   except getopt.GetoptError:
      print 'subrand.py -i <inputfile> -t <faqtype> -c <coverage> -r <refsize> -m <mean> -s <std>'
      print '  input file: fasta or fastq file of reads from which extract a subsample'
      print '  faqtype: format of output file, fasta or fastq'
      print '  coverage: desired coverage for subsample'
      print '  refsize: reference size or expected genome size needed to calculate coverage, in Mb'
      print '  mean,std: mean and standard deviation of read lengths desired for subsample. If not defined, a uniform random selection of reads is performed'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'test.py -i <inputfile> -t <faqtype>'
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

   type=inputfile.split(".")[-1]
   if type == "fq":
    faqtype="fastq"
   elif type == "fa":
    faqtype="fasta"
   else:
    faqtype=type

   print "\n Getting a",cov,"X subsample of reads from", inputfile
   print "   (Coverage calculated for a reference size of %0.1f Mb)" % (refsize/1000000)



   ok=read(inputfile,faqtype,cov,refsize)
   if(ok==0):
     readnrandom(inputfile,faqtype,ofaqtype,cov,refsize)
     if(plot):
      plotnow()


if __name__ == "__main__":
   main(sys.argv[1:])

 
