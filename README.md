# random_subreads

Python script for subsampling PacBio reads for the study ERP021971 (will add reference when available). 
This script need as input the location of the PacBio file and the coverage to extract. The final shape of the
subsample is selected according to the shape of the ONT dataset (2D-Pass 31X S288c) (will add reference when available).

### Usage:
#### subrand.py -i \<inputfile\>  -c \<coverage\> -p \<compfile\>  
 
where:

  input file: PacBio fasta or fastq file of reads from which to extract a subsample
    
  coverage: desired coverage for subsample [31]  

  compfile: (optional) fastq file of ONT reads 

### Output:

   subreads...fasta:  fasta or fastq file with the selected reads
  
   lengths.pdf: plot comparing (a) the read length distribution of the original file, 
			       (b) the read length distribution of the selected subsample,
                               [c] if mean and std are selected: a gaussian with mean and std, having the same area as (b)
			       [d] if a compfile is defined: the read length distribution of the cmopfile 

### Data
   The PacBio data used in the paper are from the S288C strain of Saccharomyces cerevisiae, accession numbers: ERR1655125 ERR1655118 ERR1655119
   To use this this script download the PacBio folders for the mentioned accession numbers and extract fastq reads in a pacbio.fastq file, then run
   the script as:  python subrand.py -i pacbio.fastq 
   To compare in a figure the extracted subsample distribution and the ONT distribution, download the data from accession numbers
   ERR1883389 ERR1883398 ERR1883399 ERR1883400 ERR1883401 ERR1883402 and generate a fastq file from the 2D Pass reads in a ont.fastq file, then run
   the script as: python subrand.py -i pacbio.fastq -p ont.fastq

### Warnings
   This branch works only for the special case described above, as the final subsample shape is fixed. For a more general case use the master branch.

### Requirements:
Python 2 or 3, modules needed: biopython, numpy, matplotlib
