use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

fn read_plink_bed_file<P: AsRef<Path>>(bed_file_reference, bim_file_reference, fam_file_reference) -> std::io::Result<Vec<Vec<u8>>> {
    // open the .bed file
    let bed_file = File::open(&bed_file_reference).expect("Failed to open bed file");
    let bim_file = File::open(&bim_file_reference).expect("Failed to open bim file");
    let fam_file = File::open(&fam_file_reference).expect("Failed to open fam file");
    let mut reader = BufReader::new(bed_file);

    // read the number of individuals and SNPs
    let num_inds = BufReader::new(fam_file).lines().count();
    let num_snps = BufReader::new(bim_file).lines().count();

    // read the magic number and file format version
    let mut magic_number = [0u8; 2];
    reader.read_exact(&mut magic_number)?;
    let third_byte = reader.read_u8()?;

    // create a 2D array to store the genotype data
    let mut genotype_data = vec![vec![0u8; num_snps as usize]; num_inds as usize];

    // read in the genotype data
    for i in 0..num_inds {
        for j in 0..num_snps {

            //if byte and bit index is determined whether we are in SNP major or individual major format
            let byte_idx = if third_byte == 1 {
                j / 4
            } else {
                i / 4
            } as usize;
            let bit_idx = if third_byte == 1 {
                (j % 4) * 2
            } else {
                (i % 4) * 2
            };
            
            // reads in exactly 1 byte
            let mut buf = [0u8; 1];
            reader.read_exact(&mut buf)?;

            //this line reads in the shifts the bits right by the index value and then takes the two right-most bits
            let genotype = (buf[0] >> bit_idx) & 0b11;

            //stores genotype
            genotype_data[i as usize][j as usize] = genotype;
        }
    }

    println!("{:?}", genotype_data);

    let mut file = File::create("output.txt").expect("Unable to create file");
    for row in &genotype_data {
        let row_str = row.iter()
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(file, "{}", row_str).expect("Unable to write to file");
    }

    Ok(genotype_data)
}

fn main() -> std::io::Result<()> {

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
    
    let genotype_data = read_plink_bed_file(bed_file_reference, bim_file_reference, fam_file_reference)?;
    // let genotype_data = read_plink_bed_file("example.bed")?;
    Ok(())
}
