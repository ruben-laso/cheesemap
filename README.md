# `cheesemap` (`🧀map`)

This is the `cheesemap` project.
This project provides a set of data structures for indexing points in point clouds,
inspired by the amazing properties of a collection of slices of cheese (**I'm not kidding**).

Our lemma is "🧀 is the best way to know your neighbors".

# Flavours of `🧀map`
Our 🧀 is available in the following flavours:
- `chs::Dense`: A dense grid of voxels storing the points in each cell.

- `chs::Sparse`: A sparse grid of voxels storing the points in each cell.
Empty voxels are not stored.

- `chs::Mixed`: A combination of dense and sparse grids.
The space is divided into slices, which can be either dense or sparse.
Note that different slices can have different sparsity patterns, resembling the holes in different slices of cheese.

**Yes, the name of the data structure comes from this analogy!**.

Of course, which flavour is best for you depends on the specific use case.

# Citation

If you use `🧀map` in your research, please cite the latest paper(s):

```
@misc{laso2025cheesemaphighperformancepointindexingdata,
      title={Cheesemap: A High-Performance Point-Indexing Data Structure for Neighbor Search in LiDAR Data},
      author={Ruben Laso and Miguel Yermo},
      year={2025},
      eprint={2502.11602},
      archivePrefix={arXiv},
      primaryClass={cs.DS},
      url={https://arxiv.org/abs/2502.11602},
}
```

# Building and installing

See the [BUILDING](BUILDING.md) document.

# Contributing

See the [CONTRIBUTING](CONTRIBUTING.md) document.

# Licensing
The code in this project is licensed under GNU Affero General Public License v3.0. See the [LICENSE](LICENSE) document for more information.
