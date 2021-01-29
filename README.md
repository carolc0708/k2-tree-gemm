# k2-tree-gemm
## CPU ver.
#### Tree Representation 
* To compile the program: `g++ -o k2tree k2tree.cpp`
* Execute: `./k2tree`
* Note: to change the k value of the tree, the `BLOCK_SIZE` in `k2tree.h` should be adjust accordingly
    * **fixed-k mode**: enter a `k` value, the `BLOCK_SIZE` should be set as k^2
    * **hybrid-k mode**: enter 0 for the `k` value, the `BLOCK_SIZE` should be set as the square of the smallest prime factor of matrix length 

#### GEMM
* spmv: sparse matrix (k2-tree) * dense vector (1d vector) = dense vector (1d vector)
* spgemm: sparse matrix (k2-tree) * sparse matrix (k2-tree) = sparse matrix (k2-tree)
    * `spgemm()` only support multiplication of trees with same k value and same matrix length


#### Graph Algorithm
* bfs: `g++ -o bfs bfs.cpp`, then `./bfs`

## GPU ver.
* To compile the program in `k2tree-gpu`: `nvcc -std=c++11 -O3 -w -arch=sm_60 -maxrregcount=64 -rdc=true -o k2tree k2tree.cu`
* Note: it is verified on only CUDA 10.0