use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::env;
use std::io::Read;

const BUF_SIZE: usize = 1000;

fn missing(filename: &str) {

    let file_base_name = filename; // replace with the actual file name without the extension
    let bed_file_name = format!("{}.bed", file_base_name);
    let fam_file_name = format!("{}.fam", file_base_name);
    let bim_file_name = format!("{}.bim", file_base_name);

    let mut bed_file = File::open(&bed_file_name).expect("Failed to open bed file");
    let bim_file = File::open(&fam_file_name).expect("Failed to open bim file");
    let fam_file = File::open(&bim_file_name).expect("Failed to open fam file");

    // read the number of individuals and SNPs
    let num_individuals = BufReader::new(fam_file).lines().count();
    let num_snps = BufReader::new(bim_file).lines().count();

    
    // read the magic number and file format version
    let mut reader = BufReader::new(bed_file);
    let mut magic_number = [0u8; 2];
    reader.read_exact(&mut magic_number);
    let third_byte = reader.read_u8();

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

    println!(third_byte)

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


fn aaf(filename: &str) {
    println!("Called aaf with argument: {}", filename);
}

fn trans10(filename: &str) {
    println!("Called trans10 with argument: {}", filename);
}

fn trans01(filename: &str) {
    println!("Called trans01 with argument: {}", filename);
}

fn main() {
    
    // Get the command line arguments
    let args: Vec<String> = env::args().collect();

    // Check that the input file names were provided
    if args.len() != 3 {
        eprintln!("Usage: {} <command> <dataset>", args[0]);
        std::process::exit(1);
    }

    let option = &args[1];
    let string = &args[2];

    match option.as_str() {
        "missing" => missing(string),
        "aaf" => aaf(string),
        "trans10" => trans10(string),
        "trans01" => trans01(string),
        _ => println!("Invalid option: {}", option),
    }
}
