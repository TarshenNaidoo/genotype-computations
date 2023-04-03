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
                    if j * 4 + k > num_individuals - 1{
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

fn aaf(prefix: &str) {
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

    let mut aaf = Vec::with_capacity(num_snps);
    aaf.extend(vec![0.0; num_snps]);
    let mut allelle_count = Vec::with_capacity(num_snps);
    allelle_count.extend(vec![0.0; num_snps]);

    let mut j_iter = (num_individuals/4) as usize;

    if num_individuals % 4 > 0 {
        j_iter += 1;
    }

    if snp_major {
        for i in 0..num_snps {
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
                    if j * 4 + k > num_individuals - 1 {
                        break;
                    }
                    // let bit_idx = j % 4;

                    // if j == 0 {
                    //     println!("shift: {}", format!("{:#010b}", buf[0] >> k*2));
                    // }
                    // println!("{:?}", genotype & 0b11);

                    //Homo aa 1, hetero aa +0.5 to alternate allele count. Missing allele +1 to missing allele
                    if (buf[0] >> k*2) & 0b11 == 0b00 {
                        aaf[i] += 2.0;
                    } else if (buf[0] >> k*2) & 0b11 == 0b10 {
                        aaf[i] += 1.0;
                    } else if (buf[0] >> k*2) & 0b01 == 0b00 {
                        allelle_count[i] += 1.0;
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

    for ind in 0..aaf.len() {
        aaf[ind] *= 10000.0;
        aaf[ind] /= 2.0 * (num_individuals as f64 - allelle_count[0]);
        aaf[ind] = aaf[ind].round() / 10000.0;
    }

    let mut file = File::create(format!("{}.freq", prefix)).unwrap();
    for element in &aaf {
        let element_str = element.to_string();
        let element_bytes = element_str.as_bytes();
        file.write_all(element_bytes).unwrap();
        file.write_all(b"\n").unwrap();
    }
}

fn trans10(prefix: &str) {
    // Open bed, bim, and fam files
    let bed_file = File::open(format!("{}.bed", prefix)).unwrap();
    let mut bed_reader = BufReader::new(bed_file);

    let bim_file = File::open(format!("{}.bim", prefix)).unwrap();
    let bim_reader = BufReader::new(bim_file);

    let fam_file = File::open(format!("{}.fam", prefix)).unwrap();
    let fam_reader = BufReader::new(fam_file);

    // Ignore first two bytes of bed file
    let mut magic_buf = [0; 2];
    bed_reader.read_exact(&mut magic_buf).unwrap();

    // Read third byte to determine snp major
    let mut snp_major_buf = [0; 1];
    bed_reader.read_exact(&mut snp_major_buf).unwrap();
    let snp_major = snp_major_buf[0] == 1;

    if !snp_major {
        panic!("bed file is already in individual-major format!");
    }

    // Calculate number of individuals and number of SNPs
    let num_individuals = fam_reader.lines().count();
    let num_snps = bim_reader.lines().count();
    println!("num_individuals: {}", num_individuals);
    println!("num_snps: {}", num_snps);

    let mut bed_vec = vec![vec![0 ; num_individuals] ; num_snps];

    let mut j_iter = (num_individuals/4) as usize;

    if num_individuals % 4 > 0 {
        j_iter += 1;
    }

    if snp_major {
        for i in 0..num_snps {
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
                    if j * 4 + k > num_individuals - 1 {
                        break;
                    }
                    // let bit_idx = j % 4;
                    // if j == 0 {
                    //     println!("shift: {}", format!("{:#010b}", buf[0] >> k*2));
                    // }
                    let genotype = (buf[0] >> k*2) & 0b11;
                    // println!(" {:?} ", j * 4 + k);
                    bed_vec[i][j * 4 + k] = genotype;
                    // print!(" {:?} ", bed_vec[i][j * 4 + k]);

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

    let mut tbed_vec = vec![vec![0; num_snps]; num_individuals];
    // println!("{:?}", bed_vec);

    for i in 0..num_snps {
        for j in 0..num_individuals {
            tbed_vec[j][i] = bed_vec[i][j];
        }
    }

    // let mut ind_major = Vec::new();

    let mut j_iter = (num_snps/4) as usize;

    if num_snps % 4 > 0 {
        j_iter += 1;
    }
    // println!("j_iter {:?}", j_iter );

    let mut file = File::create(format!("{}.tbed", prefix)).unwrap();


    file.write(&magic_buf).unwrap();
    file.write(&[0u8; 1]).unwrap();

    for i in 0..num_individuals {
        for j in 0..(j_iter) {
            let mut buf = Vec::new();
            let mut byte: u8 = 0;
            for k in 0..4 {
                if j * 4 + k > num_snps - 1 {
                    break;
                }
                // println!("{:?}", j * 4 + k );
                byte += tbed_vec[i][j * 4 + k] << 2 * (3 - k);
            }

            print!(" {} ", j);
            println!("shift: {} ", format!("{:#010b}", byte));
            buf.push(byte);

            file.write(&buf).unwrap();

        }
    }
    
    // for i in 0..num_snps {
    //     for j in 0..num_individuals {
    //         tbed_vec[j][i] = bed_vec[i][j];
    //     }
    // }
    
    // let mut file = File::create(format!("{}.tbed", prefix)).unwrap();
    // println!("{:?}", bed_vec[0][0]);
    // for i in 0..num_individuals {
    //     for j in 0..num_snps {
    //         file.write(&[tbed_vec[i][j]]).unwrap();
    //         // transposed_vec2d[j][i] = vec2d[i][j];
    //     }
    // }

    // let mut file = File::create(format!("{}.tbed", prefix)).unwrap();
    // for element in &tbed_vec {
    //     let element_str = element.to_string();
    //     let element_bytes = element_str.as_bytes();
    //     file.write_all(element_bytes).unwrap();
    //     file.write_all(b"\n").unwrap();
    // }
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
