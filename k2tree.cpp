#include "k2tree.h"

int main() {

    // k2-tree representation
    //char* filename = "matrix/hollywood-2009.mtx";
    char* filename = "matrix/ash85.mtx"; //85x85
    //char *filename = "matrix/test4.mtx"; //4x4
    k2tree *tree = new k2tree(filename, 2);

    // printout leaf group
    for (auto it : tree->leafgroup) {
        std::cout << "[" << it.first << "]" << " ---- " << std::endl;
        for (auto iv : (it.second)) {
            std::cout << iv.first << "," << iv.second << std::endl;
        }
    }

    // multiply with dense vector ---
    // an example dense vector (assume vector only contains 0 and 1)
    int rows = tree->mat_height, cols = tree->mat_width;
    std::vector<int> dv;
	for(int i = 0; i < rows; i++){
	    if (i < 85) dv.push_back(1);
	    else dv.push_back(0);
	}

    // multiply with sparse matrix ---
    k2tree *tree2 = new k2tree(filename, 2);

    tree->spgemm(tree2);

    // release
    delete tree;
    delete tree2;

    return 0;
}