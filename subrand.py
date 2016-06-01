#!/usr/bin/env python

import sys, getopt
from Bio import SeqIO
import random


def read(inputfile,faqtype,cov,refsize):
  records = (r for r in SeqIO.parse(inputfile, faqtype))
  totall=0.0
  totreads=0.
  for r in records:
     totall+=len(r.seq)
     totreads+=1
  pcov=totall/refsize
  if(cov > pcov): 
    print "\n Sorry your original file does not have the needed coverage! \n  Total coverage=%0.1fX, desired coverage=%0.1fX" % (pcov, cov)
    return 0
  else: 
    print "  Initial coverage=%0.1fX" % pcov
    return totall


def readnrandom(inputfile,faqtype,ofaqtype,cov,refsize,totall,mean,std): 

  records = (r for r in SeqIO.parse(inputfile, faqtype))

  readlist=[]
 
  value=1-cov/(totall/refsize)  # for uniform selection 
  # for gaussian selection  std=5000  


  mysum=0
  allmean=0
  i=-1
  ii=0
  for r in records:
   allmean+=len(r.seq)
   i+=1
   if mysum*1./refsize > cov:
    if(1<0): print cov, "% hit!"
   else: 
    test=0
    if(mean==0):
      ran=random.randint(0,100)
      if ran > value: test=1
    else:
      ran=random.gauss(mean,std)
      if len(r.seq) > ran: test=1

    if test:  
     readlist.append(r)
     mysum+=len(r.seq)
     ii+=1
 
#  print "Total=",i, "Selected=",ii,"Ratio=",float(ii)/float(i)
  if(float(mysum)/refsize < cov): 
    print "  Not enough coverage selected:",float(mysum)/refsize 
    return 0 
  else:
    print "  Final number of bases=%0.0f,  Coverage=%0.1fX" % (mysum,float(mysum)/refsize)
    print "  Initial read length mean:",allmean/i,",  Final read length mean=",mysum/ii

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
   mean= 0
   std=0
   try:
      opts, args = getopt.getopt(argv,"i:t:c:r:m:s",["ifile=","faqtype=","cov=","refsize=","mean=","std="])
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

   print "\n Getting a",cov,"X subsample of reads from", inputfile, "wrt a reference size of ", refsize 

   totall=read(inputfile,faqtype,cov,refsize)
   if(totall): readnrandom(inputfile,faqtype,ofaqtype,cov,refsize,totall,mean,std)

if __name__ == "__main__":
   main(sys.argv[1:])


