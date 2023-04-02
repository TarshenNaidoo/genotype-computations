
Data
----

There are three sets of sample data given

a. small

   I am also including two text files small.tped and small.tfam which is a more human friendly way of
   representing genomic data. Each line of small.tfam contains the genotypes for a particular, with
   each individual being represented by the two variants inherited from father and mother. A 0 represents
   missingness (there is always two 0s or no 0s, can't have partially missing data). The reason for including
   this is so that you can see what the data is and design tests for your code.
b. mid (1000 SNPs, 100 individuals)
c. large (~2m SNPs; ~5k individuals)

I am also including the output files .imiss and .freq for the larger (computing the results for the small one is a good exercise to make sure you undersand)


The following is a set of clarifications for assignmment 1

1. The --aaf and --missing flags only need to work for SNP-order data
2. For simplicity you may assume that the bed file maximum size will be 4GB. However, you should think about
   what changes would be required to change to terabyte size data (ideally changing two letters in your code and
   recompiling)
3. You may find the xxd program useful This allows viewing/editing of binary files. Read the manual page but here are
   some hints
      xxd see.bed  (show the bed file in hexadecimal)
      xxd -s 3 see.bed (show the bed file in hex, skipping the magic number)
      xxd -s 3 -b see.bed (show the bed file in binary, skipping the magic number)
      xxd -s 3 -c 5 -b see.bed (as above but show 5 bytes per line)
