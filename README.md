# genotype-computations

To compile:

rustc wink

To run:

wink <command> <dataset name>

commands:

missing: calcuates the per individual missingness. Expects SNP-major format .bed file
aaf: calculates the alternate allele frequency. Expects SNP-major format
trans10: transposes .bed file to individual-major format. Expects SNP-major format
trans01: transposes .bed file to SNP-major format. Expects individual-major format

Make sure <dataset name>.bed <dataset name>.bim, and <dataset name>.fam are located in the same directory as the wink.rs file

.check file should be the same as .bed file

Lots of curses we uttered in the making of this

Detailed explanation of each command

missing:

reads in the bed, bim and fam files

Uses number of lines of bim file to determine number of SNPS, and the same for fam and individuals

The vector missing_nums will store the sum of all missing genotypes for each individual

The first for loop (i) to extract the SNPs is simple but the second (for the individuals) have to be grouped into groups of 4 with a secondary loop to loop through the byte itself. 4 iterations with 2 bits each. This is so that the loop can be terminated and transitioned to the next byte for the next SNP.

The specific genotype retrieved through the k loop is done by right shifting the buffer 0 bits followed by 2, 4, and 6 for each subsequent iteration. Afterwards, bitwise AND operation with 0b11 is performed to pull just the 2 rightmost bits. If the value is 0b01, then the index in the missing_nums vector increments by 1.

The vector is divided by num_snps and then written to file

aaf:

This function is very similar to 'missing', however, the aaf vector tracks the sum of alternate allele in the bed file, 2 for ob00, or 1 for 0b10. allele_count tracks occurences of 0b00 (missing genotypes).

aaf is then divided by allele_counts and then written to file


trans10:

Also similar to 'missing', the difference being: the contents of genotype is stored into a 2d vector bed_vec, which is then transposed and stored in tbed_vec. bytes are then reassembled from the elements of tbed_vec using the same logic as the disassembly from earlier and written to file.

trans01:

Exactly the same as 'trans10', but from individual-major to SNP-major, i.e. SNP and individual variables such as num_snps and num_individuals are essentially swapped
