#!/usr/bin/env python
from __future__ import print_function  
import sys, getopt
from Bio import SeqIO
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os.path


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

  global totall
  global readall
  global meantot
  global inin
  global randord
  global lrecords

  global totN
  global weigth


  lrecords = list(SeqIO.parse(inputfile, faqtype))

  listord=range(0,len(lrecords))
  nn=len(listord)
  for ir in range(0,nn):
    thisr=random.choice(listord)  
    randord.append(thisr)
    listord.remove(thisr)

 

  global warning
  if mean is 0:
    warning=0.94
  else:
    warning=0.4
  
  weigth=1.
  totN=0
  loctot=0

  for ir in range(0,len(lrecords)):
     totall+=len(lrecords[ir].seq)
     readall.append(len(lrecords[ir].seq))
     if mean > 0 and len(lrecords[ir].seq)< mean+2*std and len(lrecords[ir].seq)> mean-2*std: 
        meantot+= len(lrecords[ir].seq)
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
     
    binwidth=400
    fig = plt.figure()
    plt.hist(readall, bins=np.arange(0, max(readall) + binwidth, binwidth),histtype='step',color='b',label="Initial PacBio reads")
    values, bins, _ = plt.hist(subreads,bins=np.arange(0, max(subreads) + binwidth, binwidth),histtype='step',color='r',label="Sub-selected PacBio reads")
    if len(compall) > 0:
      plt.hist(compall,bins=np.arange(0, max(compall) + binwidth, binwidth),histtype='step',color='g',label="ONT-to-emulate reads")


    plt.legend(bbox_to_anchor=(0, 0, 1, 1), loc="upper right", borderaxespad=0.,fontsize=10, title='Read length distribution')


      
    plt.ylabel('N')
    plt.xlabel('Read lengths')
    plt.tight_layout()
    plt.show()

    pdfname='lengths.pdf'
    myfig = PdfPages(pdfname)
    myfig.savefig(fig)
    myfig.close()

def gaussian(x, mu, sig):
   n=1./np.sqrt(2*np.pi*sig*sig)
   exp=(x - mu)*1./sig
   return n*np.exp(-exp*exp/2)


def readnrandom(inputfile,faqtype,ofaqtype,cov,refsize): 

  global subreads
  readlist=[]


  if mean>0:
    xx = np.linspace(0, max(readall), 100)
    vv=[gaussian(x,mean,std) for x in xx]
    gmax=max(vv)   
    vv2=[gaussian(x,mean2,std2) for x in xx]
    gmax2=max(vv2)   
    vv3=[gaussian(x,omean,ostd) for x in xx]
    gmax3=max(vv3)   

 
  
  if(mean==0):
    mytotal=totall
  else:
    mytotal=meantot

  mysum=0.
  allmean=0
  i=-1
  ii=0

  # go through reads and subselect them
  for ir in randord:
   thislen=len(lrecords[ir].seq)
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
     
      model1=1*gaussian(thislen,mean,std)/gmax
      if model1>1: model1=1.
      model2=0.28*gaussian(thislen,mean2,std2)/gmax2 
      model3=corrfact*gaussian(thislen,omean,ostd)/gmax3


      model=model1+model2-model3
      ran=random.uniform(0,1)

      if ran <= model: chosen=1

    if chosen: ## add to selected reads  
     readlist.append(lrecords[ir])
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
  print (' subrand.py -i <inputfile>  -c <coverage> -x <compfile>')
  print ('\n  input file: fasta or fastq file of reads from which extract a subsample')
  print ('  coverage: desired coverage for subsample [31X]')
  print ('  compfile: comparison file:  fasta or fastq file of reads to compare with the subsample')

  print (' This branch is to run only the PacBio ont-emu case described in paper [add_ref_when_available]')
  print (' For a more general case, check the master branch in https://github.com/fg6/random_subreads.git')
  
def main(argv):
   global cov
   global mean
   global std
   global omean
   global ostd
   global mean2
   global std2
   global randord
   global corrfact
   randord=[]

   inputfile = ''
   ofaqtype = 'fastq'
   cov = 31
   refsize = 12*1000000
 
   mean=25000
   std=12000
   omean = 2500
   ostd = 5000
   corrfact = 0.3
   mean2=1000
   std2=4000
           
   
   try:
     (opts, args) = getopt.getopt(argv,"i:c:p:",["ifile=","cov=","compfile="])
 
   except getopt.GetoptError:
     usage()
     sys.exit(2)
   for opt, arg in opts:
     if opt == '-h':
       usage()
       sys.exit()
     elif opt in ("-i", "--ifile"):
       inputfile = arg
     elif opt in ("-c", "--cov"):
       cov = float(arg)
     elif opt in ("-p", "--cfile"):
       compfile = arg


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
