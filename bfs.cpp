#include "k2tree.h"
/**
* graph algorithm template - bfs
*
*/
int main() {

    // generate matrix representation
    char* filename = "matrix/ash85.mtx";
    k2tree *A = new k2tree(filename, 2);

    // initialize frontier vector
    int mat_len = A->mat_height; // a newly generated len with padding
    std::vector<int> f(mat_len, 0);

    // initialize source point e.g. 1
    f[1] = 1;

    // traverse
    int count = std::count(f.begin(), f.end(), 0);
    int iter = 0;
    int newcount = count-1;
    while (newcount < count) {
        std::cout << "=== iteration: " << iter << " === not visited: " << count << std::endl;
        count = newcount;
        f = A->spmv(f);
        newcount = std::count(f.begin(), f.end(), 0);
        iter ++;
    }

    return 0;
}