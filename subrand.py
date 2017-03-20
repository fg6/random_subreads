#!/usr/bin/env python
from __future__ import print_function  
import sys, getopt
from Bio import SeqIO
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os.path
from scipy.stats import norm

plot=1

totall=0
readall=[]

subreads=[]
meantot=0
mean=0
std=1000
warning=0.6 # about min ratio of coverage needed/initial coverage (avg prob for a read to get chosen if inside the gaussian peak)

compall=[]

def read(inputfile,faqtype,cov,refsize):
  records = (r for r in SeqIO.parse(inputfile, faqtype)) 
  global totall
  global readall
  global meantot
  global inin

  global totN
  global weigth


  global warning
  if mean is 0:
    warning=0.94
  else:
    warning=0.4
  
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


  if(cov >= pcov): 
    print ("\n Sorry your original file does not have enough coverage! \n  Total coverage in the selected area=", "{:.0f}".format(pcov),"X, desired coverage=", "{:.0f}".format(cov),"X")
    return 1
  elif cov/pcov >= warning:
    print ("\n !!! Warning !!!")
    if mean == 0:
      print(" Probably not enough coverage in the originl file: initial coverage is",
            "{:.0f}".format(pcov),"X, can extract about",
            "{:.1f}".format((warning)*pcov),"X-","{:.1f}".format((warning+0.05)*pcov),"X")
    else:
      print(" Probably not enough coverage in the selected area ("
            ,mean,"+/-",std,") initial coverage for selected area=",
            "{:.0f}".format(pcov),"X, for a Gaussian shape select a coverage of about",
            "{:.0f}".format((warning)*pcov),"X  (Higher coverage => less Gaussian shape)")
      

    return 0
  else: 
    print ("\n Initial coverage for selected area=", "{:.0f}".format(pcov),"X")
    return 0

def compare(cfile,faqtype):
  print("reading comp file")
  comparec = (r for r in SeqIO.parse(cfile, faqtype)) 
  for r in comparec:
     compall.append(len(r.seq))
  return 0

def plotnow():
    plt.rc('xtick', labelsize=8)  # x axis labels
    plt.rc('ytick', labelsize=8)   # y axis labels
    font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 10}
    plt.rc('font', **font)
     
    fig = plt.figure()
    plt.hist(readall,100,histtype='step',color='b',label="Initial reads")
    values, bins, _ = plt.hist(subreads,100,histtype='step',color='r',label="Sub-selected reads")
    if len(compall) > 0:
      plt.hist(compall,100,histtype='step',color='g',label="Comparison reads")


    plt.legend(bbox_to_anchor=(0, 0, 1, 1), loc="upper right", borderaxespad=0.,fontsize=10, title='Read length distribution')


    # plot gaussian
    if(mean is not 0):
      sub_area = sum(np.diff(bins)*values)
      x = np.linspace(0, max(readall), 100)
      plt.plot(x,norm.pdf(x, mean, std)*sub_area,color='r')
      
    plt.ylabel('N')
    plt.xlabel('Read lengths')
    plt.tight_layout()
    plt.show()

    pdfname='lengths.pdf'
    myfig = PdfPages(pdfname)
    myfig.savefig(fig)
    myfig.close()

def gaussian(x, mu, sig):
    return 1.*np.exp(-np.power((x - mu)/sig, 2.)/2)


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

  # choose a random, non-repetitive order to go through reads:
  randord=[]
  listord=range(0,len(records))
  nn=len(listord)
  for ir in range(0,nn):
    thisr=random.choice(listord)  
    randord.append(thisr)
    listord.remove(thisr)


  # go through reads and subselect them
  for ir in randord:
   thislen=len(records[ir].seq)
   allmean+=thislen
   i+=1

   if mysum/refsize > cov:
    if(1<2): print (" ... {:.0f}".format(cov),"X reached!")
    break
   else: 
    chosen=0

    if(mean==0):  # random selection
      value=1-cov/(mytotal/refsize)
      ran=random.randint(0,100)
      if ran > value: chosen=1
    else:
      weigth=1
      if(thislen<3000): weigth=0.8
      model=gaussian(thislen,mean,std)
      sigma_up=model*weigth  
      ran=random.uniform(0,1)

      if ran <= sigma_up: chosen=1

    if chosen: ## add to selected reads  
     readlist.append(records[ir])
     mysum+=thislen   
     subreads.append(thislen) 
     ii+=1


  print("\n Min length selected:",min(subreads), " Total number of reads selected: ", ii)

  if(float(mysum)/refsize < cov): 
    print (" Finished: could not extract enough coverage, extracted coverage=%0.2fX" % (float(mysum)/refsize))
    return 0 
  else:
    print (" Final number of bases=%0.0f,  Coverage=%0.0fX" % (mysum,float(mysum)/refsize))
    print (" Initial read mean-length: %0.1f  Final read mean-length: %0.1f" %(allmean/i,mysum/ii))

  ofile="subreads_%0.0fX.%s" % (cov,ofaqtype)
  if mean > 0: ofile="subreads_%0.0fX_mean%d.%s" % (cov,mysum/ii,ofaqtype)
  output_handle = open(ofile,"w")  
  SeqIO.write(readlist, output_handle,ofaqtype )
  output_handle.close()
  print (" Subsample written in",ofile,", read length distribution plotted in lengths.pdf\n")


def usage():
  print (' Missing required variables ! Usage: ')
  print (' subrand.py -i <inputfile> -t <faqtype> -c <coverage> -r <refsize> -m <mean> -s <std>  -x <compfile>')
  print ('\n  input file: fasta or fastq file of reads from which extract a subsample')
  print ('  faqtype: format of output file, fasta or fastq [fasta]')
  print ('  coverage: desired coverage for subsample')
  print ('  refsize: reference size or expected genome size needed to calculate coverage, in Mb ')
  print ('  mean,std: mean and standard deviation of read lengths desired for subsample. If not defined, a uniform random selection of reads is performed')
  print ('  compfile: comparison file:  fasta or fastq file of reads to compare with the subsample')


def main(argv):
   inputfile = ''
   compfile = ''
   ofaqtype = 'fasta'
   cov = 0
   refsize = 0
   global mean
   global std

   try:
      (opts, args) = getopt.getopt(argv,"i:t:c:r:m:s:p:",["ifile=","faqtype=","cov=","refsize=","mean=","std=","compfile="])   
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
     if opt == '-h':
       usage()
       sys.exit()
     elif opt in ("-i", "--ifile"):
       inputfile = arg
     elif opt in ("-p", "--cfile"):
       compfile = arg
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
         
   print(compfile)
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

      
   print ("\n Getting a", "{:.0f}".format(cov),"X subsample of reads from", inputfile)
   print ("   (Coverage calculated for a reference size of %0.1f Mb)" % (refsize/1000000))
  
   
   # read file and check if enough data for desired subsample
   ok=read(inputfile,faqtype,cov,refsize)



   if len(compfile) is not 0 and os.path.exists(compfile):
     ctype=compfile.split(".")[-1]
     if ctype == "fq":
       faqctype="fastq"
     elif ctype == "fa":
       faqctype="fasta"
     else:
       faqctype=ctype

     compare(compfile,faqctype)
   if(ok==0): 
     readnrandom(inputfile,faqtype,ofaqtype,cov,refsize)    
     if(plot):
      plotnow()
   

if __name__ == "__main__":
   main(sys.argv[1:])
