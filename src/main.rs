mod genomic;
mod graph;
mod graphsage;
mod graphtensor;
mod loadvariants;
mod sage;

/*
Gaurav Sablok
codeprog@icloud.com
*/

fn main() -> tch::Result<()> {
    let device = Device::Cpu;
    let vs = nn::VarStore::new(device);

    // Load data
    let variants = load_variants("clinvar_snvs.tsv");
    let genes = vec![Gene {
        symbol: "BRCA1".to_string(),
        entrez_id: "672".to_string(),
    }]; // extend
    let (graph, _) = build_genomic_graph(variants, genes);

    // Split: 70% train, 30% test
    let n = graph.node_count();
    let mut train_mask = vec![false; n];
    for (i, node) in graph.node_indices().enumerate() {
        if let Node::Variant(v) = &graph[node] {
            if v.clinvar_label.is_some() {
                train_mask[i] = i % 10 < 7; // 70%
            }
        }
    }

    let (x, edge_index, y, train_mask_tensor, _) = graph_to_tensors(&graph, &train_mask);

    let x = x.to_device(device);
    let edge_index = edge_index.to_device(device);
    let y = y.to_device(device);
    let train_mask = train_mask_tensor.to_device(device);

    // Model
    let model = GraphSAGE::new(vs.root(), 4, 64, 2);
    let mut opt = nn::Optimizer::adam(&vs, 1e-3, Default::default());

    // Train
    for epoch in 1..=300 {
        opt.zero_grad();
        let logits = model.forward_t(&x, &edge_index, true);
        let loss = logits
            .masked_select(&train_mask)
            .cross_entropy_for_logits(&y.masked_select(&train_mask));

        loss.backward();
        opt.step();

        if epoch % 50 == 0 {
            let pred = logits.argmax(-1, false);
            let acc = pred
                .masked_select(&train_mask)
                .eq(&y.masked_select(&train_mask).argmax(-1, false))
                .f_mean(Kind::Float)?
                .double_value(&[]);
            println!(
                "Epoch {} | Loss: {:.4} | Train Acc: {:.1}%",
                epoch,
                loss.double_value(&[]),
                acc * 100.0
            );
        }
    }

    vs.save("graphsage_genomic.pt")?;
    Ok(())
}
