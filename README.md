# random_subreads

Python script to randomly select a subsample of reads with desired coverage from a fasta or fastq file.
The subsample can be selected having the read lengths followin a gaussian distribution a chosen mean and standard deviation, for instance to match a different data sample.


### Usage:
#### subrand.py -i \<inputfile\> -t \<faqtype\> -c \<coverage\> -r \<refsize\> -m \<mean\> -s \<std\>
 
where:

  input file: fasta or fastq file of reads from which extract a subsample
 
   faqtype: format of output file, fasta or fastq
   
   coverage: desired coverage for subsample
   
   refsize: reference size or expected genome size needed to calculate coverage, in Mb
   
   mean,std: mean and standard deviation of read lengths distribution desired for subsample. If not defined, reads are randomly selected uniformly in read length
