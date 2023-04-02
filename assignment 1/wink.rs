use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

const BUF_SIZE: usize = 1000;

fn missing_arthritis<P: AsRef<Path>>(&bed_file_reference, &bim_file_reference, &fam_file_reference) -> std::io::Result<Vec<Vec<u8>>> {

    
    let bed_file = File::open(&bed_file_reference).expect("Failed to open bed file");
    let bim_file = File::open(&bim_file_reference).expect("Failed to open bim file");
    let fam_file = File::open(&fam_file_reference).expect("Failed to open fam file");

    // read the number of individuals and SNPs
    let mut num_individuals = BufReader::new(fam_file).lines().count();
    let mut num_snps = BufReader::new(bim_file).lines().count();

    // Open the bed file and read in the genotypes
    let mut buf = [0u8; BUF_SIZE];
    let mut genotypes = Vec::new();
    loop {
        let bytes_read = bed_file.read(&mut buf).unwrap();
        if bytes_read == 0 {
            break;
        }
        for i in 0..bytes_read {
            let byte = buf[i];
            for bit_idx in (0..8).step_by(2) {
                let genotype = (byte >> bit_idx) & 0b11;
                genotypes.push(genotype);
            }
        }
    }

    // Compute the per-individual missingness
    let mut missingness = Vec::new();
    for i in 0..num_individuals {
        let start_idx = i * num_snps;
        let end_idx = start_idx + num_snps;
        let num_missing = genotypes[start_idx..end_idx]
            .iter()
            .filter(|&&genotype| genotype == 0b11)
            .count();
        let missingness_rate = num_missing as f64 / num_snps as f64;
        missingness.push(missingness_rate);
    }

    // Write the missingness results to a file
    let mut imiss_file = BufWriter::new(File::create("arthritis.imiss").unwrap());
    for m in missingness {
        let _ = write!(imiss_file, "{:.4}\n", m);
    }
}

fn main() {
    
    // Get the command line arguments
    let args: Vec<String> = env::args().collect();

    // Check that the input file names were provided
    if args.len() != 4 {
        eprintln!("Usage: {} <bed_file> <bim_file> <fam_file>", args[0]);
        std::process::exit(1);
    }

    // Open the input files
    // let bed_file = File::open(&args[1]).expect("Failed to open bed file");
    // let bim_file = File::open(&args[2]).expect("Failed to open bim file");
    // let fam_file = File::open(&args[3]).expect("Failed to open fam file");

    let bed_file_reference = &args[1];
    let bim_file_reference = &args[2];
    let fam_file_reference = &args[3];
    
    read_plink_bed_file(bed_file_reference, bim_file_reference, fam_file_reference)?;
}
