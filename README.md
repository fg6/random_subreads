# random_subreads

Python script to randomly select a subsample of reads with desired coverage from a fasta or fastq file.
The subsample can be randomly selected or selected having the read lengths following a gaussian distribution with a chosen mean and standard deviation, for instance to match a different data sample.


### Usage:
#### subrand.py -i \<inputfile\> -t \<faqtype\> -c \<coverage\> -r \<refsize\> -m \<mean\> -s \<std\> -o \<omean\> -u \<ostd\>  -f \<corrfact\>    -p \<compfile\>  
 
where:

  input file: fasta or fastq file of reads from which extract a subsample
 
   faqtype: format of output file, fasta or fastq  [fasta]
   
   coverage: desired coverage for subsample
   
   refsize: reference size or expected genome size needed to calculate coverage, in Mb
   
   mean,std: mean and standard deviation of read lengths distribution desired for subsample. If not defined, reads are randomly selected 

   omean,ostd,corrfact: to get a gaussian shape you likely need to suppress the area around the original read-length peak: set here original mean, original std, and correction factor (0:1), 0 being not corrected, maximally corrected. If not defined, oringinal shape not corrected for, subset gaussian shape will likely be deformed 

   compfile: fasta or fastq file of reads to compare to the subsample read length distribution

### Output:

   subreads...fasta:  fasta or fastq file with the selected reads
  
   lengths.pdf: plot comparing (a) the read length distribution of the original file, 
			       (b) the read length distribution of the selected subsample,
                               [c] if mean and std are selected: a gaussian with mean and std, having the same area as (b)
			       [d] if a compfile is defined: the read length distribution of the cmopfile 

### Warnings: 
   With the random selection, each read has the same probability to be picked, but as the original read-length distribution will follow a peaked shape,
   read-lengths around the peak will have higher chances to get picked, as there are more reads in that area. Thus the random selection will reproduce
   the original shape but at a lower depth (chosen by the parameter -c). 
   When mean and std are defined, the probability to get picked has been modified to follow the shape of a gaussian, but the actual shape is still affected by the
   original distribution. Generally the more the desired coverage is close to the original one, or the further from a flat distribution is the original shape
   the more the final shape will differ from a gaussian.
   The gaussian option is not optimized, you can try to correct for the original shape by using the parameters -o, -u and -f, but a manual correction 
   of the original shape might be necessary.

### Requirements:
Python 2 or 3, modules needed: biopython, numpy, matplotlib
