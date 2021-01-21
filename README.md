# k2-tree-gemm
#### Tree Representation 
* To compile the program: `g++ -o k2tree k2tree.cpp`
* Execute: `./k2tree`
* Note: to change the k value of the tree, the `BLOCK_SIZE` in `k2tree.h` should be adjust accordingly
    * **fixed-k mode**: enter a `k` value, the `BLOCK_SIZE` should be set as k^2
    * **hybrid-k mode**: enter 0 for the `k` value, the `BLOCK_SIZE` should be set as the square of the smallest prime factor of matrix length 

#### GEMM
* spmv: sparse matrix (k2-tree) * dense vector (1d vector) = dense vector (1d vector)
* spgemm: sparse matrix (k2-tree) * sparse matrix (k2-tree) = sparse matrix (k2-tree)
* Note: `spgemm()` result hasn't been verified

#### Graph Algorithm
* bfs: `g++ -o bfs bfs.cpp`, then `./bfs`