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

#define BLOCK_LEN 4
using namespace std;

class k2tree
{
    public:
        // sparse matrix metadata
        int mat_height = 0;
        int mat_width = 0;
        int mat_nnz = 0;
        int k = k;

        // representation result
        std::string T_string;
        unordered_map<std::string, vector<pair<std::string, bitset<BLOCK_LEN>>>> leafgroup; // [rowidrange] : (colidrange, leaf_bset)

    public:
        // ************* constructor and destructor *************
        k2tree(const char* filename, int k) {
            buildK2Tree(filename, k);
            std::cout << "T_string: " << T_string << std::endl;
        }
        ~k2tree() {} // leave empty for now

        // ************* tree build up functions *************
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
            }

            // k2-tree metadata
            this->mat_height = m_tmp;
            this->mat_width = n_tmp;
            this->mat_nnz = nnz_mtx_report;

            // process tree representation -----------------------
            this->T_string = getBlockRep(csrRowIdx_tmp, csrColIdx_tmp, 0, 0, 0, 0, ceil(m_tmp/k), k, 0, m_tmp, 0, m_tmp);

            // free tmp space
            free(csrColIdx_tmp);
            free(csrRowIdx_tmp);

            return 0;
        }
        std::string getBlockRep(int *csrRowIdx_tmp, int *csrColIdx_tmp, int rowmin, int rowmax,
                        int colmin, int colmax, int sub_block_len, int k, int rmin_ind, int rmax_ind, int cmin_ind, int cmax_ind) {

            bool returnflag = false;
            if (sub_block_len == 1) returnflag = true; // record leaf block result when next-level len is 1

            // init block bucket // TODO: fix! this cannot be generalized to block size > 2
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
            for (int i=0; i<this->mat_nnz; i++) {
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
                        if (!returnflag) result_tmp += getBlockRep(csrRowIdx_tmp, csrColIdx_tmp, code[0]-'0', (code[0]-'0')+1, code[2]-'0', (code[2]-'0')+1,
                                                        sub_block_len/k, k, ids[0], ids[1], ids[2], ids[3]); //TODO: fix! should fix this as when blen > 2, code is not single char
                    }
                }
            }

            // convert result bitstring to bitset
            bitset<BLOCK_LEN> result_bset(result);

            // leaf return
            if (returnflag) {
                std::string code(1, '0' + rowmin);
                code += "_";
                code += ('0' + colmin);

                int *lids = code2ind(code, k, rmin_ind, rmax_ind, cmin_ind, cmax_ind);

                //record the result
                if (result_bset.count() != 0) { // we care only non-empty leaf block

                    std::string cind = "";
                    cind += std::to_string(lids[2]);
                    cind += "-";
                    cind += std::to_string(lids[3]);

                    std::string rind = "";
                    rind += std::to_string(lids[0]);
                    rind += "-";
                    rind += std::to_string(lids[1]);

                    pair<std::string, bitset<BLOCK_LEN>> content;
                    content.first = cind;
                    content.second = result_bset;

                    this->leafgroup[rind].push_back(content);
                }

                return ""; }

            // return result
            if (result_bset.count() == 0) return "";
            else return result + " " + result_tmp;
        }

        // ************* multiplication functions *************
        void spmm(std::vector<std::vector<int>> dm) {
            // 1: block-wise multiplication
            // init output vector
            std::vector<std::vector<int>> output;
            int outrows = this->mat_height, outcols = dm[0].size();
            for(int i = 0; i < outrows; i++){
                std::vector<int> temp;
                for(int j = 0; j < outcols; j++){
                    temp.push_back(0);
                }
                output.push_back(temp);
            }

            // multiply dense matrix
            for (auto it : this->leafgroup) {
                // split rids
                std::vector<int> rids = tokenizeIDs(it.first);
                //std::cout << rids[0] << " " << rids[1] << std::endl;

                // for each rids
                if (rids[0] >= this->mat_height || rids[1] >= this->mat_height) continue;
                for (int i=rids[0]; i<=rids[1]; i++) { // for each rids
                    for (auto iv : (it.second)) { // for each cids
                       std::vector<int> cids = tokenizeIDs(iv.first);
                       //std::cout << cids[0] << " " << cids[1] << std::endl;
                       bitset<BLOCK_LEN> block_bset = iv.second;
                       int cnt = 0;

                       if (cids[0] >= this->mat_width || cids[1] >= this->mat_width) continue;
                       for (int j=cids[0]; j<=cids[1]; j++) {
                            int temp = block_bset[cnt];
                            temp *= dm[j][i];
                            output[i][j] += temp;
                            cnt ++;
                       }
                    }
                }
            }

            // print the output for debug
            std::cout << "--- spmm output result: ---" << std::endl;
            for(int i=0; i<outrows; i++) {
                for(int j=0; j<outcols; j++) {
                    std::cout << output[i][j] << " ";
                }
                std::cout << std::endl;
            }
        }

        void spmv(std::vector<int> dv) {
            // init output vector
            std::vector<int> outputv;
            int outrows = this->mat_height;
            for(int i = 0; i < outrows; i++){
                outputv.push_back(0);
            }

            // multiply dense vector
            for (auto it : this->leafgroup) {
                // split rids
                std::vector<int> rids  = tokenizeIDs(it.first);
                //std::cout << rids[0] << " " << rids[1] << std::endl;

                // for each rids
                if (rids[0] >= this->mat_height || rids[1] >= this->mat_height) continue;
                for (int i=rids[0]; i<=rids[1]; i++) { // for each rids
                    for (auto iv : (it.second)) { // for each cids
                       std::vector<int> cids = tokenizeIDs(iv.first);
                       //std::cout << cids[0] << " " << cids[1] << std::endl;
                       bitset<BLOCK_LEN> block_bset = iv.second;
                       int cnt = 0;

                       if (cids[0] >= this->mat_width || cids[1] >= this->mat_width) continue;
                       for (int j=cids[0]; j<=cids[1]; j++) {
                            int temp = block_bset[cnt];
                            temp *= dv[j];
                            outputv[i] += temp;
                            cnt ++;
                       }
                    }
                }
            }

            std::cout << "--- spmv output result: ---" << std::endl;
            for(int i=0; i<outrows; i++) {
                std::cout << outputv[i] << " ";
            }
            std::cout << std::endl;


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

        // given code (block location) and indices range, split the block and output indices //TODO: Fix! <- code range should be generalized
        int *code2ind(std::string code, int k, int rmin_ind, int rmax_ind, int cmin_ind, int cmax_ind) {

            static int result[4];
            if (code == "0_0") { result[0] = rmin_ind; result[1] = (rmin_ind+rmax_ind)/k; result[2] = cmin_ind; result[3] = (cmin_ind+cmax_ind)/k; }
            if (code == "0_1") { result[0] = rmin_ind; result[1] = (rmin_ind+rmax_ind)/k; result[2] = ceil(float(cmin_ind+cmax_ind)/k); result[3] = cmax_ind; }
            if (code == "1_0") { result[0] = ceil(float(rmin_ind+rmax_ind)/k); result[1] = rmax_ind; result[2] = cmin_ind; result[3] = (cmin_ind+cmax_ind)/k; }
            if (code == "1_1") { result[0] = ceil(float(rmin_ind+rmax_ind)/k); result[1] = rmax_ind; result[2] = ceil(float(cmin_ind+cmax_ind)/k); result[3] = cmax_ind; }

            return result;
        }

};

#endif