/*
Gaurav Sablok
codeprog@icloud.com
*/

use petgraph::Graph;
use petgraph::graph::NodeIndex;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum NodeType {
    Variant,
    Gene,
    Regulator,
}

#[derive(Debug)]
struct Variant {
    chrom: String,
    pos: u32,
    ref_allele: String,
    alt_allele: String,
    clinvar_label: Option<i32>, // 1 = pathogenic, 0 = benign, None = unknown
    cadd_phred: f32,
    gnomad_af: f32,
}

#[derive(Debug)]
struct Gene {
    symbol: String,
    entrez_id: String,
}
