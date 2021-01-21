#ifndef K2TREE_H
#define K2TREE_H

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <utility>
#include <bits/stdc++.h> // stringstream, bitset
#include "mmio.h"

#define BLOCK_SIZE 4 // k*k, need to be predefined
using namespace std;

// binary matrix multiplication
bitset<BLOCK_SIZE> bmm(bitset<BLOCK_SIZE> a, bitset<BLOCK_SIZE> b) {

    bitset<BLOCK_SIZE> result;

    int n = sqrt(BLOCK_SIZE);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {

            bitset<BLOCK_SIZE> value;
            for (int m=0; m<n; m++) {
                value |= (a[BLOCK_SIZE-1-(i*n+m)] & b[BLOCK_SIZE-1-(m*n+j)]);
                if (value[0]) break;
            }
            result[BLOCK_SIZE-1-(i*n+j)] = value[0];
        }

    }

    return result;
}

class k2tree
{
public:
    // sparse matrix metadata
    int mat_height = 0;
    int mat_width = 0; // assume the same as mat_height
    int mat_nnz = 0;
    int k = 0;

    int* csrRowIdx_tmp;
    int* csrColIdx_tmp;

    vector<int> prime; // in hybrid-k mode, record the prime factors of mat_height

    // representation result
    std::string T_string;
    unordered_map<std::string, unordered_map<std::string, bitset<BLOCK_SIZE>>> leafgroup; // [rowidrange] : ([colidrange]: leaf_bset)

public:
    // ************* constructor and destructor *************
    k2tree(const char* filename, int k) {
        if (filename == "") {
                T_string = "";
                std::cout << "Construct an empty tree. " << std::endl;
            } else {
            buildK2Tree(filename, k);
        }

    }
    ~k2tree() {} // leave empty for now

        // ************* tree build up functions *************
 /**
 * Read the matrix from .mtx file and build the k2-tree
 *
 * @param filename  The filename (.mtx extension) of the ultra-sparse matrix
 * @param k  For fixed-k mode, the specified k value; For hybrid-k mode, remains 0
 * @return   1 or -1 to indicate successfully read in the matrix or not
 */
    int buildK2Tree(const char *filename, int k) {

        // load matrix -----------------------
        int m_tmp, n_tmp, nnz_tmp;

        int ret_code;
        MM_typecode matcode;
        FILE *f;

        int nnz_mtx_report;
        int isInteger = 0, isReal = 0, isPattern = 0, isSymmetric_tmp = 0, isComplex = 0;

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


        this->csrRowIdx_tmp = (int *)malloc(nnz_mtx_report * sizeof(int));
        this->csrColIdx_tmp = (int *)malloc(nnz_mtx_report * sizeof(int));

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
            this->csrRowIdx_tmp[i] = idxi;
            this->csrColIdx_tmp[i] = idxj;
        }

        // k2-tree metadata -----------------------
        this->mat_height = m_tmp;
        this->mat_width = n_tmp;
        this->mat_nnz = nnz_mtx_report;


        // process tree representation -----------------------
        // choice of k-split
        if (k > 0) { // fixed-k, get padding
            this->k = k;
            this->getPadding();
            this->T_string = getBlockRep(0, this->mat_height, 0, this->mat_height, ceil(float(this->mat_height)/this->k), -1);
        } else { // hybrid-k, get prime factors for a sequence of k
            this->getPrimeFactor(this->mat_height);
            this->k = this->prime[(this->prime).size()-1];
            this->T_string = getBlockRep(0, this->mat_height, 0, this->mat_height, ceil(float(this->mat_height)/this->k), (this->prime).size()-1);
        }

        // free tmp space
        free(this->csrRowIdx_tmp);
        free(this->csrColIdx_tmp);

        return 0;
    }

 /**
 * Recursive function to record the binary string of block at each level
 *
 * @param rmin_ind, rmax_ind  The row index range of the block to be processed
 * @param cmin_ind, cmax_ind  The col index range of the block to be processed
 * @param sub_block_len  The sub-block size to split the block
 * @param k_ind  For hybrid-k mode, k2tree::prime[k_ind] indicates k value at this level;
 *               For fixed-k mode, the value remains -1
 * @return      part of T_string for non-leaf blocks; "" for leaf blocks
 */
    std::string getBlockRep(int rmin_ind, int rmax_ind, int cmin_ind, int cmax_ind, int sub_block_len, int k_ind) {

        bool returnflag = false;
        int parent_k;
        if (k_ind == -1) { // fixed-k mode
            if (sub_block_len < this->k) returnflag = true; // record leaf block result when next-level len smaller than k
        } else { // hybrid-k mode
            if (k_ind == 0) returnflag = true;
            parent_k = this->k; // parent k
            this->k = this->prime[k_ind];
            //std::cout << "sub block len: " << sub_block_len << ", this level k: " << this->k << std::endl;
        }


        // init block bucket
        unordered_map<string, int> block_bucket;
        for (int l=0; l<this->k; l++){
            for (int r=0; r<this->k; r++) {

                std::string code = genRangeCode(l, r);
                block_bucket[code] = 0;
            }
        }

        // start recording
        for (int i=0; i<this->mat_nnz; i++) {
            if ((rmin_ind <= this->csrRowIdx_tmp[i] && this->csrRowIdx_tmp[i] < rmax_ind)
            && (cmin_ind <= this->csrColIdx_tmp[i] && this->csrColIdx_tmp[i] < cmax_ind)) {

                std::string code = genRangeCode((this->csrRowIdx_tmp[i]-rmin_ind)/sub_block_len,
                                                (this->csrColIdx_tmp[i]-cmin_ind)/sub_block_len);
                block_bucket[code] += 1;
            }
        }

        // output level result
        std::string result = "";
        std::string result_tmp = "";

        for (int l=0; l<this->k; l++){
            for (int r=0; r<this->k; r++) {

                std::string code = genRangeCode(l, r);
                if (block_bucket[code] == 0) result += "0";
                else {
                    result += "1";
                    if (!returnflag) {
                        vector<int> ids = code2ind(code, rmin_ind, rmax_ind, cmin_ind, cmax_ind);
                        if (k_ind == -1) // fixed-k mode
                            result_tmp += getBlockRep(ids[0], ids[1], ids[2], ids[3], ceil(float(sub_block_len)/this->k), -1);
                        else // hybrid-k mode
                            result_tmp += getBlockRep(ids[0], ids[1], ids[2], ids[3], ceil(float(sub_block_len)/this->prime[k_ind-1]), k_ind-1);
                    }
                }
            }
        }

        // convert result bitstring to bitset
        bitset<BLOCK_SIZE> result_bset(result);

        // return the leaf block
        if (returnflag) {
            if (result_bset.count() != 0) { // we care only non-empty leaf block

                std::string rind = genRangeCode(rmin_ind, rmax_ind);
                std::string cind = genRangeCode(cmin_ind, cmax_ind);
                this->leafgroup[rind][cind] = result_bset;
            }

            if (k_ind != -1) this->k = parent_k; // for hybrid-k mode
            return "";
        }

        // return result
        if (k_ind != -1) this->k = parent_k; // for hybrid-k mode
        if (result_bset.count() == 0) return "";
        else return result + " " + result_tmp;
    }

    // ************* multiplication functions *************
    std::vector<int> spmv(std::vector<int> dv) {
        // init output vector
        std::vector<int> outputv;
        int outrows = this->mat_height;
        for(int i = 0; i < outrows; i++){
            outputv.push_back(0);
        }

        // multiply dense vector
        for (auto it : this->leafgroup) {

            std::vector<int> rids  = tokenizeIDs(it.first);
            for (int r=rids[0]; r<rids[1]; r++) { // for each rids
                for (auto iv : (it.second)) { // for each cids
                   std::vector<int> cids = tokenizeIDs(iv.first);
                   bitset<BLOCK_SIZE> block_bset = iv.second;
                   int cnt = (r-rids[0]) * sqrt(BLOCK_SIZE);
                   for (int c=cids[0]; c<cids[1]; c++) {
                        int temp = block_bset[BLOCK_SIZE-1-cnt];
                        temp *= dv[c];
                        outputv[r] += temp;
                        cnt ++;
                   }
                }
            }
        }

        // print out vector result
        for(int i=0; i<outrows; i++) {
            std::cout << outputv[i] << " ";
        }
        std::cout << std::endl;

        return outputv;
    }

    k2tree* spgemm(k2tree* B) { // Cik = Sum(A_ij * B_jk), considering only non-empty leaf blocks
        // warning message
        if (this->mat_height != B->mat_height) std::cout << "[k2tree::spgemm()] Invalid matrix size for matrix multiplication!" << std::endl;
        if (this->k != B->k) std::cout << "[k2tree::spgemm()] k value mismatch for matrix multiplication!" << std::endl;

        // initialize a new tree
        int block_len = sqrt(BLOCK_SIZE);
        k2tree *output = new k2tree("", block_len); // construct a empty structure
        output->mat_height = this->mat_height;

        // multiplication
        int n = this->mat_height;

        for (int i=0; i<n; i+=block_len) {
            std::string i_ind = genRangeCode(i, i+block_len);
            for (int j=0; j<n; j+=block_len) {
                std::string j_ind = genRangeCode(j, j+block_len);
                bitset<BLOCK_SIZE> value;
                for (int m=0; m<n; m+=block_len) {
                    std::string m_ind = genRangeCode(m, m+block_len);

                    // add to global sum
                    if (this->leafBlockExist(i_ind, m_ind) && B->leafBlockExist(m_ind, j_ind))
                        value |= bmm(this->leafgroup[i_ind][m_ind], B->leafgroup[m_ind][j_ind]);
                }
                if (value.count() != 0) output->leafgroup[i_ind][j_ind] = value;
            }
        }

        return output;
    }


    // ************* utility functions *************
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

    // given code (block location) and indices range, split the block and output indices
    vector<int> code2ind(std::string code, int rmin_ind, int rmax_ind, int cmin_ind, int cmax_ind) {
        vector<int> result{0, 0, 0, 0};

        vector<int> codeid = tokenizeIDs(code);
        int r = codeid[0], c = codeid[1];
        int r_interval = ceil(float(rmax_ind - rmin_ind)/this->k), c_interval = ceil(float(cmax_ind - cmin_ind)/this->k);

        result[0] = (rmin_ind + r_interval * r);
        result[1] = (rmin_ind + r_interval * (r+1));
        result[2] = (cmin_ind + c_interval * c);
        result[3] = (cmin_ind + c_interval * (c+1));
        return result;
    }

    // check if leaf block at certain index exist
    bool leafBlockExist(std::string rids, std::string cids) {
        if (this->leafgroup.find(rids) == this->leafgroup.end()) return false;
        else if (this->leafgroup[rids].find(cids) == this->leafgroup[rids].end()) return false;
        return true;
    }

    // calculate the padding and replace metadata
    void getPadding() {
        int mat_len = pow(this->k, ceil(log(this->mat_height)/log(this->k)));
        this->mat_height = mat_len;
        this->mat_width = mat_len;

        std::cout << "Add paddings to the matrix, now mat_len = " << mat_len << std::endl;
    }

    // generate range code based on given ids
    std::string genRangeCode(int min, int max) {
        std::string code = "";
        code += std::to_string(min);
        code += "-";
        code += std::to_string(max);
        return code;
    }

    // prime factorization for hybrid-k tree
    void getPrimeFactor(int n) {
        vector<int> prime;
        // Print the number of 2s that divide n
        while (n%2 == 0)
        {
            prime.push_back(2);
            n = n/2;
        }

        // n must be odd at this point.  So we can skip
        // one element (Note i = i +2)
        for (int i = 3; i <= sqrt(n); i = i+2)
        {
            // While i divides n, print i and divide n
            while (n%i == 0)
            {
                prime.push_back(i);
                n = n/i;
            }
        }

        // This condition is to handle the case when n
        // is a prime number greater than 2
        if (n > 2) prime.push_back(n);

        // print out result
        std::cout << "hybrid-k values by level: ";
        for (int i=0; i<prime.size(); i++) std::cout << prime[prime.size()-1-i] << " ";
        std::cout << std::endl;

        this->prime = prime;
    }

    // ************* print functions *************
    void printTree() {
        std::cout << "*********************************" << std::endl;
        std::cout << "T_string: " << T_string << std::endl;

        // printout leaf group
        for (auto it : this->leafgroup) {
            std::cout << "[" << it.first << "]" << " ---- " << std::endl;
            for (auto iv : (it.second)) {
                std::cout << iv.first << "," << iv.second << std::endl;
            }
        }

        std::cout << "*********************************" << std::endl;

    }
};

#endif