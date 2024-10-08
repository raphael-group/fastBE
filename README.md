## fastBE: An ultrafast method for phylogenetic reconstruction from bulk DNA sequencing data

*fastBE* is a method for inferring the evolutionary history
of tumors from multi-sample bulk DNA sequencing data.
Our method uses ideas from distance based phylogenetics and 
a handcrafted solver of the variant allele frequency 
$\ell_1$-regression problem.

![overview](docs/overview.png)

If you find this tool useful in your research, please cite us at:
```
@article{schmidt2024regression,
  title={A regression based approach to phylogenetic reconstruction from multi-sample bulk DNA sequencing of tumors},
  author={Schmidt, Henri and Raphael, Benjamin J},
  journal={bioRxiv},
  pages={2024--04},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Installation

#### Option 1: conda install

To install using `conda` execute the following command:
```
$ conda install schmidt73::fastbe
```
Currently, `fastBE` is only available through `conda` for Linux. Additional
platforms will be supported upon request.

#### Option 2: Pre-compiled binaries

Pre-compiled binaries for `fastBE` are available for both macOS and Linux under the releases section. 
To install, simply copy the binary into a directory which is recognized by your shell.

#### Option 3: Build from source

`fastBE` is implemented in C++ and is packaged with the dependencies
needed to execute the program. The only dependencies are
a recent version of CMAKE and a modern C++17 compliant compiler.

To build `fastBE` from source, first clone the repository and its submodules:
```
$ git clone --recurse-submodules https://github.com/schmidt73/fastbe.git
```

Then from the root of the project directory, execute the following sequence o
commands:
```
$ mkdir build; cd build
$ cmake ..
$ make
```
The output binary will be located at `build/src/fastbe`.

## Usage

To run *fastbe*, simply execute the binary. 
```
$ fastbe
Usage: fastbe [--help] [--version] {cluster,regress,search}

Optional arguments:
  -h, --help     shows help message and exits
  -v, --version  prints version information and exits

Subcommands:
  cluster       clusters the mutations (columns) of a frequency matrix.
  regress       regresses a clone tree onto a frequency matrix.
  search        searches for a clone tree that best fits a frequency matrix.
```

The three modes of fastbe are `search`, `regress`, and `cluster`. 
The `search` mode infers the phylogenetic tree best fitting the input 
frequency matrix, the `regress` mode infers the fit of the input 
tree to the frequency matrix, and the `cluster` mode clusters mutations
into clones using phylogeny constrained clustering. If one is interested in 
inferring a phylogenetic tree from a frequency matrix, the `search` mode 
should be used. 

Formally, the `search` mode
solves the variant allele frequency $\ell_1$-factorization problem, 
while the `regress` mode solves the variant allele frequency 
$\ell_1$-regression problem, both of which are defined in 
our manuscript. 
The `search` mode takes as input an $m \times n$ frequency matrix $F$ 
and outputs an $n$-clonal tree that best fits the frequency matrix. 
The `regress` mode
takes as input an $m \times n$ frequency matrix $F$ and an $n$-clonal
tree $\mathcal{T}$ and outputs the minimum value of 
$\lVert F - UB_{\mathcal{T}} \rVert_1$ over all usage matrices $U$.
The `cluster` mode takes as input an $n$-mutation tree $\mathcal{T}$,
an $m \times n$ frequency matrix $F$, and the number of clusters $k$. 
The `cluster` mode then outputs a clustering of the $n$ mutations into
$k$ clusters.

> [!IMPORTANT] 
> The search command requires a root vertex specified with the 
> `-f/--assigned-root` flag, corresponding to the appropriate column of the 
> frequency matrix. By default this root vertex is set to be $0$, or the
> first column of the frequency matrix. When the root vertex is unknown,
> it suffices to append an extra column to the beginning of the frequency matrix 
> and specify the root as $0$. 

> [!TIP]
> Use the first column of the frequency matrix as the mutation off the root.

### Input format

The input format for the `search` mode of *fastbe* consists of a frequency 
matrix $F$ in `.txt` format. Rows are separated by newlines
and columns are separated by spaces. Rows correspond
to distinct samples and columns correspond to distinct mutations 
(or mutation clusters).
More formally, $F_{ij}$ is the frequency of the $j^{\text{th}}$ mutation
in the $i^{\text{th}}$ sample. As an example, a frequency matrix $F$ 
describing $20$ samples and $10$ mutations is:
```
1.0000 0.3217 0.3425 0.0000 0.0000 0.0000 0.4444 0.0000 0.6857 0.0000
1.0000 0.5141 0.5128 0.0165 0.0000 0.0000 0.5417 0.3842 0.0891 0.0000
1.0000 0.0563 0.3483 0.0000 0.5000 0.0833 0.3566 0.4167 0.1527 0.4953
1.0000 0.1644 0.3261 0.0000 0.4679 0.0000 0.1826 0.4906 0.0000 0.4500
1.0000 0.0000 0.0000 0.0000 1.0000 0.0000 0.0000 1.0000 0.0000 1.0000
1.0000 0.0909 0.2800 0.0000 0.3047 0.0000 0.3220 0.2986 0.2510 0.2857
1.0000 1.0000 1.0000 0.0000 0.0000 0.0000 1.0000 0.0000 0.0000 0.0000
1.0000 0.1953 0.2206 0.1354 0.1917 0.0896 0.1339 0.5255 0.0996 0.3945
1.0000 0.0815 0.3684 0.0505 0.1500 0.0238 0.2056 0.5045 0.1148 0.2400
1.0000 0.0000 0.0476 0.0000 0.0000 0.0000 0.0000 0.8865 0.0000 0.9300
1.0000 0.1437 0.8081 0.0000 0.0152 0.0000 0.1610 0.0046 0.0000 0.0175
1.0000 0.4094 0.5625 0.0000 0.0000 0.0000 0.4054 0.0000 0.0000 0.0000
1.0000 0.0000 0.2018 0.2866 0.5050 0.0000 0.0000 0.4650 0.0000 0.4623
1.0000 0.0000 0.5275 0.4424 0.0000 0.0000 0.0000 0.0000 0.1310 0.0000
1.0000 0.0000 0.0000 0.0000 0.7236 0.0000 0.0000 0.8315 0.0545 0.7200
1.0000 0.4534 0.4459 0.0000 0.0143 0.0000 0.4766 0.1818 0.2152 0.1282
1.0000 1.0000 1.0000 0.0000 0.0000 0.0000 1.0000 0.0000 0.0000 0.0000
1.0000 0.0704 0.1757 0.0000 0.0432 0.1094 0.2000 0.2020 0.5352 0.0000
1.0000 0.2138 0.2877 0.6715 0.0175 0.0000 0.2736 0.0711 0.0000 0.0804
1.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.8703 0.0000 0.0000
```

The above frequency matrix $F$ is provided as an input file at `examples/sim_obs_frequency_matrix.txt`.

The input format for the `regress` mode of *fastbe* is the aforementioned
frequency matrix $F$ and an $n$-clonal tree $\mathcal{T}$. The tree is specified
as an adjacency list in `.txt` format. An example of a clonal tree
which consists of $10$ clones rooted at the $0$ vertex is:
```
0 2 3 7 8
1
2 6
3
4
5
6 1
7 9
8 5
9 4
```

The above clonal tree $\mathcal{T}$ is provided as an input file at `examples/sim_tree.txt`. Each
vertex in the adjacency list corresponds to a clone. The name of the vertex corresponds
to the index of the mutation in the frequency matrix on the edge leading to the vertex.

The above frequency matrix and clonal tree were generated using the command,
`python scripts/simulation.py --clones 10 --samples 20 --coverage 40 --seed 0 --mutations 40 --output examples/sim`
which simulates the evolution of a tumor with $10$ mutation clusters (equivalently, clones) and $20$ samples at a 
read depth of $40\times$. The $40$ mutations are distributed across the $10$ mutation clusters.
Several other files such as the mutation to clone mapping, the ground truth usage matrix $U$, clonal
matrix $B$, and read count matrices are also provided in the `examples/` directory.

### Example

As an example, we infer a phylogenetic tree from the simulated
data with $20$ samples and $10$ clones and $40$ mutations. 
To run `fastbe` on this data, assuming the clones are known we 
can pass in the *clonal frequency matrix* by executing:
```
fastbe search examples/sim_obs_frequency_matrix.txt -o examples/fastbe
```
This command will output an adjacency list describing the clonal tree 
at `examples/fastbe_tree.txt` and a `.json` file containing metadata
at `examples/fastbe_results.json`. Note that we did not specify the root
of the tree, so the root is assumed to be the first column of the 
frequency matrix.

When the clones are unknown, we run `fastbe` with the *mutation
frequency matrix* and specify the root as the second mutation in the
frequency matrix with the command:
```
fastbe search examples/sim_obs_full_frequency_matrix.txt -o examples/fastbe -f 1
```
As one often wants the clones, one can use the inferred tree to 
infer the clones with phylogenetically constrained clustering. To do this, 
we run the *cluster* command with the inferred *mutation tree* and 
the frequency matrix:
```
fastbe cluster -k 10 examples/fastbe_tree.txt examples/sim_obs_full_frequency_matrix.txt -o examples/fastbe -l L2
```
This command will output the inferred clones at 
`examples/fastbe_clustering.csv` and a `examples/fastbe_results.json` file 
containing important metadata. 

> [!TIP] For selecting the number of mutation
> clusters automatically, we recommend the usage of the Python
> package [kneed](https://github.com/raphael-group/fastBE/tree/v1.0.1).

### Benchmarks

The runtime of `fastbe` scales as $\mathcal{O}(m)$ with the number 
$m$ of samples and as $\mathcal{O}(n^4\log^2(n))$ with the number $n$ of mutations. As a 
reference, we provide the runtime of `fastbe` on the simulated data
with $m = 100$ samples and a varying number of mutations $n$ with
default parameters. The system used for benchmarking is an 
Intel Core i7-8665U CPU @ 1.90GHz × 8 with 8GB of RAM. 

| Mutations ($n$) | Samples ($m$) | Runtime (s) |
|-----------------|---------------|-------------|
| 10              | 100           | 4           |
| 25              | 100           | 41          |
| 50              | 100           | 148         |
| 100             | 100           | 486         |
| 200             | 100           | 923         |
| 500             | 100           | 3161        |
| 1000            | 100           | 4145        |

The memory footprint of `fastbe` is $O(mn + nb)$ where $b$ is the
beam width for almost all applications. The memory footprint of 
`fastbe` is negligible for almost all applications.

When running the `search` mode of `fastbe`, the beam width of the
search algorithm can be specified with the `-b/--beam_width` flag.
The beam width controls the number of candidate trees that are
considered at each iteration of the search algorithm. The runtime
of the search algorithm scales approximately linearly with the beam 
width, and the quality of the inferred tree improves with increased 
beam width. 

> [!TIP]
> Use the default beam width chosen by `fastbe` for most applications.
