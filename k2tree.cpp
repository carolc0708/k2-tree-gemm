#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include <vector>

#include "mmio.h"
using namespace std;

// input block range and k
// output result of k^2 bit
std::string getBlockRep(int *csrRowIdx_tmp, int *csrColIdx_tmp, int nnz_mtx_report, int rowmin, int rowmax,
                        int colmin, int colmax, int sub_block_len, int k) {

    std::cout << "sub_block_len: " << sub_block_len << std::endl; // sub_block_len of the same length means at same level
    //std::cout << "rowmin, rowmax, colmin, colmax: " << rowmin << " " << rowmax << " " << colmin << " " << colmax << std::endl;

    if (sub_block_len <= k) return ""; // supposed to be leave value

    // init block bucket
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
                //std::cout << code << std::endl;
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
                result_tmp += getBlockRep(csrRowIdx_tmp, csrColIdx_tmp, nnz_mtx_report, code[0]-'0', (code[0]-'0')+1, code[2]-'0', (code[2]-'0')+1, sub_block_len/k, k);
            }

            //std::cout << "code, result: " << code << " " << result << std::endl;
        }
    }
    //std::cout << "result: " << result << std::endl;
    result = result + result_tmp;


    return result;
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
    result += getBlockRep(csrRowIdx_tmp, csrColIdx_tmp, nnz_mtx_report, 0, 0, 0, 0, ceil(m_tmp/k), k);
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

int main() {

    // k2-tree representation
    //char* filename = "matrix/road_usa.mtx";
    char* filename = "matrix/ash85.mtx"; //84x84
    buildK2Tree(filename, 2);

}