# k2-tree-gemm
#### Tree Representation 
* To compile the program: `g++ -o k2tree k2tree.cpp`
* Execute: `./k2tree`
* Note: to change the k value of the tree, the `BLOCK_SIZE` in `k2tree.h` should be adjust to k^2 accordingly

#### GEMM
* spmm: sparse matrix (k2-tree) * dense matrix (2d vector) = dense matrix (2d vector)
* spmv: sparse matrix (k2-tree) * dense vector (1d vector) = dense vector (1d vector)
* spgemm: sparse matrix (k2-tree) * sparse matrix (k2-tree) = sparse matrix (k2-tree)
* Note: the current gemm model is still under revision.