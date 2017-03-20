# random_subreads

Python script to randomly select a subsample of reads with desired coverage from a fasta or fastq file.
The subsample can be randomly selected or selected having the read lengths following a gaussian distribution with a chosen mean and standard deviation, for instance to match a different data sample.


### Usage:
#### subrand.py -i \<inputfile\> -t \<faqtype\> -c \<coverage\> -r \<refsize\> -m \<mean\> -s \<std\>  -p \<compfile\> 
 
where:

  input file: fasta or fastq file of reads from which extract a subsample
 
   faqtype: format of output file, fasta or fastq  [fasta]
   
   coverage: desired coverage for subsample
   
   refsize: reference size or expected genome size needed to calculate coverage, in Mb
   
   mean,std: mean and standard deviation of read lengths distribution desired for subsample. If not defined, reads are randomly selected 

   compfile: fasta or fastq file of reads to compare to the subsample read length distribution

### Output:

   subreads...fasta:  fasta or fastq file with the selected reads
  
   lengths.pdf: plot comparing (a) the read length distribution of the original file, 
			       (b) the read length distribution of the selected subsample,
                               [c] if mean and std are selected: a gaussian with mean and std, having the same area as (b)
			       [d] if a compfile is defined: the read length distribution of the cmopfile 

### Warnings:
   When mean and std are defined, the shape of the subsample will be roughly gaussian, but the actual shape is affected by the
   original distribution. Generally, the more the desired coverage is close to the original one, the more the shape will differ from a gaussian.
   Will try generalize the script to become more independent from the original shape sometime in the future.

### Requirements:
Python 2 or 3, modules needed: biopython, numpy, matplotlib
