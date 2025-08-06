# `cheesemap` (`🧀map`)

This is the `cheesemap` project.
This project provides a set of data structures for indexing points in point clouds,
inspired by the amazing properties of a collection of slices of cheese (**Not kidding**).

Our lemma is "🧀 is the best way to know your neighbors".

# Features

- **High performance**: `🧀map` is designed for high-performance point indexing, making it suitable for large-scale point cloud applications.

- **$\mathcal{O}(1)$ voxel access**: All the different flavours provide (amortised) constant-time access to voxel data, ensuring fast query performance.

- **Efficient memory usage**: The sparse and mixed flavours are designed to use memory efficiently, reducing the overall memory footprint.

- **Use of custom kernel queries**: `🧀map` supports custom kernel queries, allowing users to define their own search criteria and operations on the point data as long as they comply with the `chs::concepts::Kernel` concept (file [include/cheesemap/concepts/Kernel.hpp](https://github.com/ruben-laso/cheesemap/blob/main/include/cheesemap/concepts/Kernel.hpp)).

We encourage you to check the papers down below for more information.

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

```bibtex
@article{fgcs2025cheesemap,
      title = {Cheesemap: A high-performance point-indexing data structure for neighbor search in {LiDAR} data},
      journal = {Future Generation Computer Systems},
      pages = {108060},
      year = {2025},
      issn = {0167-739X},
      doi = {https://doi.org/10.1016/j.future.2025.108060},
      url = {https://www.sciencedirect.com/science/article/pii/S0167739X25003553},
      author = {Ruben Laso and Miguel Yermo},
      keywords = {Point cloud, Data structure, Nearest neighbors, LiDAR},
}

@misc{arxiv2025cheesemap,
      title={Cheesemap: A High-Performance Point-Indexing Data Structure for Neighbor Search in {LiDAR} Data},
      author={Ruben Laso and Miguel Yermo},
      year={2025},
      eprint={2502.11602},
      archivePrefix={arXiv},
      primaryClass={cs.DS},
      url={https://arxiv.org/abs/2502.11602},
}
```

# Building and installing

TLDR: This is a header-only library, clone the source code, and include the headers in your project.

For more details, see the [BUILDING](BUILDING.md) document.

# Contributing

See the [CONTRIBUTING](CONTRIBUTING.md) document.

# Licensing
The code in this project is licensed under GNU Affero General Public License v3.0. See the [LICENSE](LICENSE) document for more information.
