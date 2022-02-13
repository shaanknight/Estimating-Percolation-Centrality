# Estimating-Percolation-Centrality
This work presents a new randomized approximation algorithm, called Algorithm Dual-Brandes, for the source-destination variant of percolation centrality for undirected graphs. Our algorithm is inspired by the source-sampling technique developed by Brandes and Pich. For baselines, we generalize existing ideas used for approximation of the similar betweenness centrality measure to estimate Percolation centrality. We experimentally evaluate our algorithm against these baselines on a collection of 12 real-world graphs. Our algorithm outperforms baselines in accuracy and other metrics, achieving a mean squared error 61\% lower than the nearest baseline.

## Dependencies

- g++ version 7.3.0 compiler
- C++

## Dataset

The graph instances used in the experiment can be found in the Dataset folder.

## Specifics

- The baseline algorithms and Dual-Brandes are referred by codes. Use the following codes : Algorithm-A, Algorithm-B, Algorithm-C, Dual-Brandes
- The code accepts a sampling set percentage as a parameter, provide a suitable integer value in [0,100].

## Run

To compile the code use :
```
g++ -O3 -static-libstdc++ approximatePercolationCentrality.cpp
```

To run the CPU codes use :
```
./approximatePC <path_to_graph_file> <sampling_percentage> <algorithm_code>
```
