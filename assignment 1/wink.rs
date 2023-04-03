use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::env;
use std::io::prelude::*;

fn missing(prefix: &str){
    // Open bed, bim, and fam files
    let bed_file = File::open(format!("{}.bed", prefix)).unwrap();
    let mut bed_reader = BufReader::new(bed_file);

    let bim_file = File::open(format!("{}.bim", prefix)).unwrap();
    let bim_reader = BufReader::new(bim_file);

    let fam_file = File::open(format!("{}.fam", prefix)).unwrap();
    let fam_reader = BufReader::new(fam_file);

    // Ignore first two bytes of bed file
    let mut buf = [0; 2];
    bed_reader.read_exact(&mut buf).unwrap();

    // Read third byte to determine snp major
    let mut snp_major_buf = [0; 1];
    bed_reader.read_exact(&mut snp_major_buf).unwrap();
    let snp_major = snp_major_buf[0] == 1;

    // Calculate number of individuals and number of SNPs
    let num_individuals = fam_reader.lines().count();
    let num_snps = bim_reader.lines().count();
    println!("num_individuals: {}", num_individuals);
    println!("num_snps: {}", num_snps);

    let mut missing_nums = Vec::with_capacity(num_individuals);
    missing_nums.extend(vec![0.0; num_individuals]);

    let mut j_iter = (num_individuals/4) as usize;

    if num_individuals % 4 > 0 {
        j_iter += 1;
    }

    if snp_major {
        for _i in 0..num_snps {
            for j in 0..(j_iter) {
                // let byte_idx = if is_snp_major == 1 {
                //     j / 4
                // } else {
                //     i / 4
                // } as usize;
                let mut buf = [0u8; 1];
                // println!("{}", (_i * j));
                bed_reader.read_exact(&mut buf).unwrap();
                // if j == 0 {
                //     println!("binary: {}", format!("{:b}", buf[0]));
                // }
                // println!("{:?}", buf);
                for k in 0..4 {
                    if j * 4 + k > num_individuals {
                        break;
                    }
                    // let bit_idx = j % 4;

                    // if j == 0 {
                    //     println!("shift: {}", format!("{:#010b}", buf[0] >> k*2));
                    // }
                    let genotype = (buf[0] >> k*2) & 0b11;
                    // println!("{:?}", genotype & 0b11);
                    if genotype == 0b01 {
                        missing_nums[j * 4 + k] += 1.0;
                    }

                    // print!("[");
                    // for element in missing_nums.iter_mut() {
                    //      print!(" {} ", element);
                    // }
                    // println!("]");
    
                }
            }
        }
    } else {

    }

    for element in missing_nums.iter_mut() {
        *element /= num_snps as f64;
    }

    let mut file = File::create(format!("{}.imiss", prefix)).unwrap();
    for element in &missing_nums {
        let element_str = element.to_string();
        let element_bytes = element_str.as_bytes();
        file.write_all(element_bytes).unwrap();
        file.write_all(b"\n").unwrap();
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
