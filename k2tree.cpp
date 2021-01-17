#include "k2tree.h"

int main() {

    // k2-tree representation
    //char* filename = "matrix/road_usa.mtx";
    char* filename = "matrix/ash85.mtx"; //84x84
    k2tree *tree = new k2tree(filename, 2);

    // printout leaf group
    for (auto it : tree->leafgroup) {
        std::cout << "[" << it.first << "]" << " ---- " << std::endl;
        for (auto iv : (it.second)) {
            std::cout << iv.first << "," << iv.second << std::endl;
        }
    }

    // multiply with dense matrix ---
    // an example dense matrix (assume matrix only contains 0 and 1)
    std::vector<std::vector<int>> dm;
    int rows = 84, cols = 84;
	for(int i = 0; i < rows; i++){
		std::vector<int> temp;
		for(int j = 0; j < cols; j++){
			temp.push_back(1);//(i * cols + j + 1);
		}
		dm.push_back(temp);
	}

    tree->spmm(dm);

    // multiply with dense vector ---
    // an example dense vector (assume vector only contains 0 and 1)
    std::vector<int> dv;
	for(int i = 0; i < rows; i++){
		dv.push_back(1);
	}

    tree->spmv(dv);

    // multiply with sparse matrix ---
    k2tree *tree2 = new k2tree(filename, 2);

    tree->spgemm(tree2);

    // release
    delete tree;
    delete tree2;

    return 0;
}