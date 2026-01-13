/*
Gaurav Sablok
codeprog@icloud.com
*/

type GenomicGraph = Graph<Node, EdgeType, petgraph::Directed>;

#[derive(Debug, Clone, Copy)]
enum EdgeType {
    CoOccurs,
    Regulates,
    InGene,
    HiCContact,
}

#[derive(Debug)]
enum Node {
    Variant(Variant),
    Gene(Gene),
    Regulator { name: String },
}
