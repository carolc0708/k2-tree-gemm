#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <utility>
#include <bits/stdc++.h> // stringstream, bitset

#include "mmio.h"
using namespace std;

/* TODOs:
* 1. change the std::string representation of binary code to bitset
* 2. use class to wrap up k2-tree representation
* 3. multiply sparse matrix, multiply sparse vector (k2-tree, csr)
* 4. column-group-wise multiplication
* 5. rectangular block
* 6. verify if block size can be generalize to 32,
*    add some test case for special number
* 7. k2-tree to .mtx output
* 8. spMM, spMV on GPU ver.
*/


// given code (block location) and indices range, split the block and output indices
int *code2ind(std::string code, int k, int rmin_ind, int rmax_ind, int cmin_ind, int cmax_ind) {

    static int result[4];
    if (code == "0_0") { result[0] = rmin_ind; result[1] = (rmin_ind+rmax_ind)/k; result[2] = cmin_ind; result[3] = (cmin_ind+cmax_ind)/k; }
    if (code == "0_1") { result[0] = rmin_ind; result[1] = (rmin_ind+rmax_ind)/k; result[2] = ceil(float(cmin_ind+cmax_ind)/k); result[3] = cmax_ind; }
    if (code == "1_0") { result[0] = ceil(float(rmin_ind+rmax_ind)/k); result[1] = rmax_ind; result[2] = cmin_ind; result[3] = (cmin_ind+cmax_ind)/k; }
    if (code == "1_1") { result[0] = ceil(float(rmin_ind+rmax_ind)/k); result[1] = rmax_ind; result[2] = ceil(float(cmin_ind+cmax_ind)/k); result[3] = cmax_ind; }

    return result;
}

// gather leafs with the same rows
// [rmin-rmax] <cmin_cmax, 0010>
unordered_map<std::string, vector<pair<std::string, std::string>>> leafgroup;


// input block range and k
// output result of k^2 bit
std::string getBlockRep(int *csrRowIdx_tmp, int *csrColIdx_tmp, int nnz_mtx_report, int rowmin, int rowmax,
                        int colmin, int colmax, int sub_block_len, int k, int rmin_ind, int rmax_ind, int cmin_ind, int cmax_ind) {

    //std::cout << "sub_block_len: " << sub_block_len << std::endl; // sub_block_len of the same length means at same level
    //std::cout << "rowmin, rowmax, colmin, colmax: " << rowmin << " " << rowmax << " " << colmin << " " << colmax << std::endl;

    bool returnflag = false;
    if (sub_block_len == 1) returnflag = true; // supposed to be leaf value

    // init block bucket // this cannot be generalized to block size > 2
    unordered_map<string, int> block_bucket;
    for (int l=0; l<k; l++){
        for (int r=0; r<k; r++) {
            std::string code = "";
            code += std::to_string(l);
            code += "_";
            code += std::to_string(r);
            block_bucket[code] = 0;
        }
    }

    // start recording
    for (int i=0; i<nnz_mtx_report; i++) {

        //std::cout << "(" << csrRowIdx_tmp[i] << "," << csrColIdx_tmp[i] << ")" << std::endl;
        int vrowid = csrRowIdx_tmp[i]/sub_block_len;
        int vcolid = csrColIdx_tmp[i]/sub_block_len;

        if (rowmin == rowmax && colmin == colmax) { // first level
            std::string code = "";
            code += std::to_string(vrowid);
            code += "_";
            code += std::to_string(vcolid);
            block_bucket[code] += 1;

        } else { // levels other than first
            if ((rowmin*k <= vrowid && vrowid < rowmax*k) && (colmin*k <= vcolid && vcolid < colmax*k)) {
                std::string code = "";
                code += std::to_string(vrowid-rowmin*k);
                code += "_";
                code += std::to_string(vcolid-colmin*k);
                //std::cout << code << " " << vrowid << " " << vcolid << std::endl;
                block_bucket[code] += 1;
            }
        }
    }

    // output level result
    std::string result = "";
    std::string result_tmp = "";

    for (int l=0; l<k; l++){
        for (int r=0; r<k; r++) {
            std::string code = "";
            code += std::to_string(l);
            code += "_";
            code += std::to_string(r);

            if (block_bucket[code] == 0) result += "0";
            else {
                result += "1";

                //process next level ind: code 00, 01, 10, 11 => ids
                int *ids = code2ind(code, k, rmin_ind, rmax_ind, cmin_ind, cmax_ind);
                //std::cout << "next level ids to pass down: " << ids[0] << "," << ids[1] << "," << ids[2] << "," << ids[3] << std::endl;
                if (!returnflag) result_tmp += getBlockRep(csrRowIdx_tmp, csrColIdx_tmp, nnz_mtx_report, code[0]-'0', (code[0]-'0')+1, code[2]-'0', (code[2]-'0')+1,
                                                sub_block_len/k, k, ids[0], ids[1], ids[2], ids[3]);
            }

            //std::cout << "code, result: " << code << " " << result << std::endl;
        }
    }

    // leaf return
    if (returnflag) {
        //std::cout << "L: " << result << std::endl;
        std::string code(1, '0' + rowmin);
        code += "_";
        code += ('0' + colmin);
        //std::cout << code << std::endl;
        //std::cout << "id pass to this level: "<< rmin_ind << "," << rmax_ind << "," << cmin_ind << "," << cmax_ind << std::endl;
        int *lids = code2ind(code, k, rmin_ind, rmax_ind, cmin_ind, cmax_ind);
        //std::cout << "Lid: " << lids[0] << "," << lids[1] << "," << lids[2] << "," << lids[3] << std::endl;

        //record the result
        if (result != "0000") { // we care only non-empty leaf block

            std::string cind = "";
            cind += std::to_string(lids[2]);
            cind += "-";
            cind += std::to_string(lids[3]);

            std::string rind = "";
            rind += std::to_string(lids[0]);
            rind += "-";
            rind += std::to_string(lids[1]);

            pair<std::string, std::string> content;
            content.first = cind;
            content.second = result;
            //std::cout << cind << "," << rind << std::endl;

            leafgroup[rind].push_back(content);
        }


        return ""; }

    //std::cout << "result: " << result << std::endl;

    if (result == "0000") return "";
    else return result + " " + result_tmp;
}


int buildK2Tree(const char *filename, int k) {
    int m_tmp, n_tmp, nnz_tmp;

    int ret_code;
    MM_typecode matcode;
    FILE *f;

    int nnz_mtx_report;
    int isInteger = 0, isReal = 0, isPattern = 0, isSymmetric_tmp = 0, isComplex = 0;

        // load matrix
    if ((f = fopen(filename, "r")) == NULL)
    {
        printf("Error loading matrix file.\n");
        return -1;
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        return -2;
    }

    if ( mm_is_pattern( matcode ) )  { isPattern = 1; /*printf("type = Pattern\n");*/ }
    if ( mm_is_real ( matcode) )     { isReal = 1; /*printf("type = real\n");*/ }
    if ( mm_is_complex( matcode ) )  { isComplex = 1; /*printf("type = real\n");*/ }
    if ( mm_is_integer ( matcode ) ) { isInteger = 1; /*printf("type = integer\n");*/ }

    /* find out size of sparse matrix .... */
    ret_code = mm_read_mtx_crd_size(f, &m_tmp, &n_tmp, &nnz_mtx_report);
    if (ret_code != 0)
        return -4;


    if ( mm_is_symmetric( matcode ) || mm_is_hermitian( matcode ) )
    {
        isSymmetric_tmp = 1;
        //printf("input matrix is symmetric = true\n");
    }
    else
    {
        //printf("input matrix is symmetric = false\n");
    }


    // read from matrix -----------------------
    int *csrRowIdx_tmp = (int *)malloc(nnz_mtx_report * sizeof(int));
    int *csrColIdx_tmp = (int *)malloc(nnz_mtx_report * sizeof(int));

    vector<vector<int>> temp(n_tmp, std::vector<int>(n_tmp, 0)); // for verifying

    for (int i=0; i<nnz_mtx_report; i++) {
        int idxi, idxj;
        double fval, fval_im;
        int ival;
        int returnvalue;

        if (isReal)
        {
            returnvalue = fscanf(f, "%d %d %lg\n", &idxi, &idxj, &fval);
        }
        else if (isComplex)
        {
            returnvalue = fscanf(f, "%d %d %lg %lg\n", &idxi, &idxj, &fval, &fval_im);
        }
        else if (isInteger)
        {
            returnvalue = fscanf(f, "%d %d %d\n", &idxi, &idxj, &ival);
            fval = ival;
        }
        else if (isPattern)
        {
            returnvalue = fscanf(f, "%d %d\n", &idxi, &idxj);
            fval = 1.0;
        }

        // adjust from 1-based to 0-based
        idxi--;
        idxj--;

        // store to temporary location
        csrRowIdx_tmp[i] = idxi;
        csrColIdx_tmp[i] = idxj;

        temp[idxi][idxj] = 1;
    }

    // process tree representation -----------------------
    std::string result = "";
    result += getBlockRep(csrRowIdx_tmp, csrColIdx_tmp, nnz_mtx_report, 0, 0, 0, 0, ceil(m_tmp/k), k, 0, m_tmp, 0, m_tmp);
    std::cout << "T result: " << result << std::endl;

//    // print out the matrix (for verifying)
//    for(int i=0; i<n_tmp; i++) {
//        for(int j=0; j<n_tmp; j++) {
//            if (temp[i][j] == 1) printf("1 ");
//            else printf("0 ");
//        }
//        printf("\n");
//    }


    // free tmp space
    free(csrColIdx_tmp);
    free(csrRowIdx_tmp);

    return 0;
}

// tokenize ids that represent as "XX-XX"
std::vector<int> tokenizeIDs(std::string id) {

    // Vector of string to save tokens
    std::vector<int> tokens;

    // stringstream class check1
    stringstream check(id);

    std::string intermediate;

    // Tokenizing w.r.t. '-'
    while(getline(check, intermediate, '-'))
    {
        int val;
        std::istringstream ss(intermediate);
        ss >> val;
        tokens.push_back(val);
    }

    return tokens;
}

int main() {

    // k2-tree representation
    //char* filename = "matrix/road_usa.mtx";
    char* filename = "matrix/ash85.mtx"; //84x84
    buildK2Tree(filename, 2);

    // printout leaf group
    for (auto it : leafgroup) {
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

	// 1: block-wise multiplication
	// init output vector
	std::vector<std::vector<int>> output;
    for(int i = 0; i < rows; i++){
        std::vector<int> temp;
        for(int j = 0; j < cols; j++){
            temp.push_back(0);
        }
        output.push_back(temp);
    }

    for (auto it : leafgroup) {
        // split rids
        std::vector<int> rids  = tokenizeIDs(it.first);
        //std::cout << rids[0] << " " << rids[1] << std::endl;

        // for each rids
        if (rids[0] >= rows || rids[1] >= rows) continue;
        for (int i=rids[0]; i<=rids[1]; i++) { // for each rids
            for (auto iv : (it.second)) { // for each cids
               std::vector<int> cids = tokenizeIDs(iv.first);
               //std::cout << cids[0] << " " << cids[1] << std::endl;
               std::string content = iv.second;
               int cnt = 0;

               if (cids[0] >= cols || cids[1] >= cols) continue;
               for (int j=cids[0]; j<=cids[1]; j++) {
                    int temp = (content[cnt]-'0');
                    temp *= dm[j][i];
                    output[i][j] += temp;
                    cnt ++;
               }
            }
        }
    }

    // print the output for debug
    std::cout << "--- output result: ---" << std::endl;
    for(int i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            std::cout << output[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // multiply with dense vector ---
    // an example dense vector (assume vector only contains 0 and 1)
    std::vector<int> dv;
	for(int i = 0; i < rows; i++){
		dv.push_back(1);
	}

    // init output vector
	std::vector<int> outputv;
    for(int i = 0; i < rows; i++){
        outputv.push_back(0);
    }

    for (auto it : leafgroup) {
        // split rids
        std::vector<int> rids  = tokenizeIDs(it.first);
        //std::cout << rids[0] << " " << rids[1] << std::endl;

        // for each rids
        if (rids[0] >= rows || rids[1] >= rows) continue;
        for (int i=rids[0]; i<=rids[1]; i++) { // for each rids
            for (auto iv : (it.second)) { // for each cids
               std::vector<int> cids = tokenizeIDs(iv.first);
               //std::cout << cids[0] << " " << cids[1] << std::endl;
               std::string content = iv.second;
               int cnt = 0;

               if (cids[0] >= cols || cids[1] >= cols) continue;
               for (int j=cids[0]; j<=cids[1]; j++) {
                    int temp = (content[cnt]-'0');
                    temp *= dv[j];
                    outputv[i] += temp;
                    cnt ++;
               }
            }
        }
    }

    std::cout << "--- vector output result: ---" << std::endl;
    for(int i=0; i<rows; i++) {
        std::cout << outputv[i] << " ";
    }
    std::cout << std::endl;

}