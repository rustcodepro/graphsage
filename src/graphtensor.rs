use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

/*
Gaurav Sablok
codeprog@icloud.com
*/

use tch::{Device, Kind, Tensor};

fn graph_to_tensors(
    graph: &GenomicGraph,
    train_mask: &[bool],
) -> (Tensor, Tensor, Tensor, Tensor, Tensor) {
    let n = graph.node_count();

    // Node features
    let mut x_vec = Vec::with_capacity(n * 4);
    let mut y_vec = Vec::with_capacity(n);
    for (i, node) in graph.node_indices().map(|i| (i.index(), &graph[i])) {
        x_vec.extend(node_features(node));
        let label = match node {
            Node::Variant(v) => v.clinvar_label.unwrap_or(-1),
            _ => -1,
        };
        y_vec.push(label as f32);
    }

    let x = Tensor::of_slice(&x_vec).f_reshape(&[n as i64, 4]).unwrap();
    let y = Tensor::of_slice(&y_vec);

    // Edge index
    let mut src = Vec::new();
    let mut dst = Vec::new();
    for edge in graph.edge_indices() {
        let (s, t) = graph.edge_endpoints(edge).unwrap();
        src.push(s.index() as i64);
        dst.push(t.index() as i64);
    }
    let edge_index = Tensor::cat(
        &[
            Tensor::of_slice(&src).unsqueeze(0),
            Tensor::of_slice(&dst).unsqueeze(0),
        ],
        0,
    );

    // Masks
    let train_mask: Vec<i8> = train_mask.iter().map(|&b| b as i8).collect();
    let train_mask = Tensor::of_slice(&train_mask);

    (
        x,
        edge_index,
        y,
        train_mask,
        Tensor::zeros(&[n as i64], (Kind::Int64, Device::Cpu)),
    )
}
