use std::env;
use std::fs::File;
use std::io::{BufReader, Read};

fn main() {
    // Get the command line arguments
    let args: Vec<String> = env::args().collect();

    // Check that the input file names were provided
    if args.len() != 4 {
        eprintln!("Usage: {} <bed_file> <bim_file> <fam_file>", args[0]);
        std::process::exit(1);
    }

    // Open the input files
    let bed_file = File::open(&args[1]).expect("Failed to open bed file");
    let bim_file = File::open(&args[2]).expect("Failed to open bim file");
    let fam_file = File::open(&args[3]).expect("Failed to open fam file");

    // Get the number of SNPs from the bim file
    let num_snps = BufReader::new(bim_file).lines().count();

    // Get the number of individuals from the fam file
    let num_inds = BufReader::new(fam_file).lines().count();

    // Wrap the bed file in a buffered reader
    let mut input_reader = BufReader::new(bed_file);

    // Read the header bytes
    let mut header = [0u8; 3];
    input_reader.read_exact(&mut header).expect("Failed to read header bytes");

    // Extract the genotype organization flag
    let snp_major = (header[2] & 0b0000_0001) != 0;

    let mut genotype_data = vec![vec![0u8; num_individuals as usize]; num_snps as usize];

    if (!snp_major) {
        genotype_data = vec![vec![0u8; num_individuals as usize]; num_snps as usize];
    }

    // Read the rest of the file
    let mut buf = [0u8; 1];
    while let Ok(_) = input_reader.read_exact(&mut buf) {
        // Extract the genotype value from the buffer
        let mut bit_idx = 6;
        for _ in 0..4 {

            // Will shift the buffer by bit_idx to extract each genotype
            let genotype = (buf[0] >> bit_idx) & 0b11;
            println!("genotype: {}", genotype);
            if snp_major {
                // Store the genotype value by SNP

                genotype_data
            } else {
                // Store the genotype value by individual
            }
            bit_idx -= 2;
        }
    }
}
