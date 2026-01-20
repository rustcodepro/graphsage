use tch::nn::{self, ModuleT, Optimizer};

/*
Gaurav Sablok
codeprog@icloud.com
*/

pub struct GraphSAGE {
    conv1: nn::GraphSageConv,
    conv2: nn::GraphSageConv,
    lin: nn::Linear,
}

impl GraphSAGE {
    pub fn new(vs: nn::Path, in_dim: i64, hidden: i64, out_dim: i64) -> Self {
        let c1 = nn::graph_sage_conv(&vs / "conv1", in_dim, hidden, "mean", true);
        let c2 = nn::graph_sage_conv(&vs / "conv2", hidden, hidden, "mean", true);
        let lin = nn::linear(&vs / "lin", hidden, out_dim, Default::default());
        GraphSAGE {
            conv1: c1,
            conv2: c2,
            lin,
        }
    }
}

impl ModuleT for GraphSAGE {
    fn forward_t(&self, x: &Tensor, edge_index: &Tensor, train: bool) -> Tensor {
        let h = self.conv1.forward_t(x, edge_index, train).relu();
        let h = self.conv2.forward_t(&h, edge_index, train).relu();
        self.lin.forward(&h)
    }
}
