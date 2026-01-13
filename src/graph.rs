/*
Gaurav Sablok
codeprog@icloud.com
*/

pub fn build_genomic_graph(
    variants: Vec<Variant>,
    genes: Vec<Gene>,
) -> (GenomicGraph, HashMap<String, NodeIndex>) {
    let mut graph = GenomicGraph::new();
    let mut id_to_idx: HashMap<String, NodeIndex> = HashMap::new();

    // Add variants
    for v in &variants {
        let id = format!("{}_{}_{}_{}", v.chrom, v.pos, v.ref_allele, v.alt_allele);
        let idx = graph.add_node(Node::Variant(v.clone()));
        id_to_idx.insert(id, idx);
    }

    // Add genes
    for g in &genes {
        let idx = graph.add_node(Node::Gene(g.clone()));
        id_to_idx.insert(g.entrez_id.clone(), idx);
    }

    // Add edges (example: in-gene, co-occurrence in 1kb)
    let mut var_list: Vec<_> = variants.iter().enumerate().collect();
    var_list.sort_by_key(|(_, v)| (v.chrom.clone(), v.pos));

    for i in 0..var_list.len() {
        let (i_idx, v1) = var_list[i];
        let v1_node = id_to_idx[&format!(
            "{}_{}_{}_{}",
            v1.chrom, v1.pos, v1.ref_allele, v1.alt_allele
        )];

        // Co-occurrence in 1kb
        for j in (i + 1)..var_list.len() {
            let (_, v2) = var_list[j];
            if v1.chrom != v2.chrom {
                break;
            }
            if v2.pos - v1.pos > 1000 {
                break;
            }

            let v2_node = id_to_idx[&format!(
                "{}_{}_{}_{}",
                v2.chrom, v2.pos, v2.ref_allele, v2.alt_allele
            )];
            graph.add_edge(v1_node, v2_node, EdgeType::CoOccurs);
            graph.add_edge(v2_node, v1_node, EdgeType::CoOccurs);
        }
    }

    (graph, id_to_idx)
}

pub fn node_features(node: &Node) -> Vec<f32> {
    match node {
        Node::Variant(v) => {
            vec![
                v.cadd_phred,
                v.gnomad_af.ln_1p(),                             // log(1 + AF)
                if v.ref_allele.len() == 1 { 1.0 } else { 0.0 }, // SNP?
                if v.alt_allele.len() == 1 { 1.0 } else { 0.0 },
            ]
        }
        Node::Gene(_) => vec![0.0, 0.0, 0.0, 0.0], // placeholder
        Node::Regulator { .. } => vec![0.0, 0.0, 0.0, 0.0],
    }
}
