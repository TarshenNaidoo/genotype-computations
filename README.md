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