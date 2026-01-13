use csv::ReaderBuilder;
use std::fs::File;

/*
Gaurav Sablok
codeprog@icloud.com
*/

pub fn load_variants(path: &str) -> Vec<Variant> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .unwrap();
    let mut variants = Vec::new();
    for result in rdr.records() {
        let record = result.unwrap();
        let chrom = record[0].to_string();
        let pos: u32 = record[1].parse().unwrap();
        let ref_allele = record[3].to_string();
        let alt_allele = record[4].to_string();
        let info: HashMap<_, _> = record[7]
            .split(';')
            .filter_map(|kv| kv.split_once('='))
            .collect();

        let clinvar = info.get("CLNSIG").and_then(|s| {
            if s.contains("Pathogenic") {
                Some(1)
            } else if s.contains("Benign") {
                Some(0)
            } else {
                None
            }
        });

        let cadd = info.get("CADD").and_then(|s| s.parse().ok()).unwrap_or(0.0);
        let af = info
            .get("gnomAD_AF")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);

        variants.push(Variant {
            chrom,
            pos,
            ref_allele,
            alt_allele,
            clinvar_label: clinvar,
            cadd_phred: cadd,
            gnomad_af: af,
        });
    }
    variants
}
