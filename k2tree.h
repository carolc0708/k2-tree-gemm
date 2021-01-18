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

class k2tree
{
    public:
        // sparse matrix metadata
        int mat_height = 0;
        int mat_width = 0;
        int mat_nnz = 0;
        int k = 0;

        // representation result
        std::string T_string;
        unordered_map<std::string, unordered_map<std::string, bitset<BLOCK_SIZE>>> leafgroup; // [rowidrange] : ([colidrange]: leaf_bset)

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
            this->k = k;

            // process tree representation -----------------------
            this->T_string = getBlockRep(csrRowIdx_tmp, csrColIdx_tmp, 0, 0, 0, 0, ceil(float(this->mat_height)/this->k), 0, this->mat_height, 0, this->mat_height);

            // free tmp space
            free(csrColIdx_tmp);
            free(csrRowIdx_tmp);

            return 0;
        }

        std::string getBlockRep(int *csrRowIdx_tmp, int *csrColIdx_tmp, int rowmin, int rowmax,
                        int colmin, int colmax, int sub_block_len, int rmin_ind, int rmax_ind, int cmin_ind, int cmax_ind) {

            //std::cout << sub_block_len << std::endl;
            bool returnflag = false;
            //if (this->k <= sub_block_len && sub_block_len < this->k*this->k) returnflag = true; // for cases that the next-level block become too small
            if (sub_block_len < this->k) returnflag = true; // record leaf block result when next-level len smaller than k

            // init block bucket
            unordered_map<string, int> block_bucket;
            for (int l=0; l<this->k; l++){
                for (int r=0; r<this->k; r++) {
                    std::string code = "";
                    code += std::to_string(l);
                    code += "-";
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
                    code += "-";
                    code += std::to_string(vcolid);
                    block_bucket[code] += 1;

                } else { // levels other than first
                    if ((rowmin*this->k <= vrowid && vrowid < rowmax*this->k) && (colmin*this->k <= vcolid && vcolid < colmax*this->k)) {
                        std::string code = "";
                        code += std::to_string(vrowid-rowmin*this->k);
                        code += "-";
                        code += std::to_string(vcolid-colmin*this->k);
                        block_bucket[code] += 1;
                    }
                }
            }

            // output level result
            std::string result = "";
            std::string result_tmp = "";

            for (int l=0; l<this->k; l++){
                for (int r=0; r<this->k; r++) {
                    std::string code = "";
                    code += std::to_string(l);
                    code += "-";
                    code += std::to_string(r);

                    if (block_bucket[code] == 0) result += "0";
                    else {
                        result += "1";

                        //process next level ind: code 00, 01, 10, 11 => ids
                        vector<int> ids = code2ind(code, rmin_ind, rmax_ind, cmin_ind, cmax_ind);
                        vector<int> codeinds = tokenizeIDs(code);
                        if (!returnflag) result_tmp += getBlockRep(csrRowIdx_tmp, csrColIdx_tmp, codeinds[0], codeinds[0]+1, codeinds[1], codeinds[1]+1,
                                                        ceil(float(sub_block_len)/this->k), ids[0], ids[1], ids[2], ids[3]);
                    }
                }
            }



            // convert result bitstring to bitset
            bitset<BLOCK_SIZE> result_bset(result);

            // leaf return
            if (returnflag) {

                //std::cout << rowmin << ", " << rowmax << ", " << colmin << ", " << colmax << ", " << sub_block_len
                //<< ", " << rmin_ind << ", " << rmax_ind << ", " << cmin_ind << ", " << cmax_ind << ", " << result << std::endl;

                //record the result
                if (result_bset.count() != 0) { // we care only non-empty leaf block
//                    // the last level split
//                    std::string code = "";
//                    code += std::to_string(rowmin);
//                    code += "-";
//                    code += std::to_string(colmin);
//
//                    vector<int> lids = code2ind(code, rmin_ind, rmax_ind, cmin_ind, cmax_ind);
//
//                    std::string cind = "";
//                    cind += std::to_string(lids[2]);
//                    cind += "-";
//                    cind += std::to_string(lids[3]);
//
//                    std::string rind = "";
//                    rind += std::to_string(lids[0]);
//                    rind += "-";
//                    rind += std::to_string(lids[1]);

                    std::string rind = "";
                    rind += std::to_string(rmin_ind);
                    rind += "-";
                    rind += std::to_string(rmax_ind);

                    std::string cind = "";
                    cind += std::to_string(cmin_ind);
                    cind += "-";
                    cind += std::to_string(cmax_ind);


                    this->leafgroup[rind][cind] = result_bset;
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
//            for (auto it : this->leafgroup) {
//                std::vector<int> rids = tokenizeIDs(it.first);
//                for (int r=rids[0]; r<rids[1]; r++) { // for each row
//                    for(auto iv : (it.second)) { // for each blocks
//                        std::vector<int> cids = tokenizeIDs(iv.first);
//                        bitset<BLOCK_SIZE> block_bset = iv.second;
//                        int cnt = (r-rids[0]) * this->k; // first row starts at 0, second row starts at 2
//                        for (int c=cids[0]; c<cids[1]; c++) { // for each cids
//
//                            int temp = block_bset[BLOCK_SIZE-1-cnt]; std::cout << r << ", " << c << ", " << temp << std::endl;
//                            temp *= dm[c][c];
//                            output[r][c] += temp;
//                            cnt ++;
//                        }
//                    }
//                }
//            }

            // without particular assumption,
            // virtually all value locations should be considered

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

                std::vector<int> rids  = tokenizeIDs(it.first);
                for (int r=rids[0]; r<rids[1]; r++) { // for each rids
                    for (auto iv : (it.second)) { // for each cids
                       std::vector<int> cids = tokenizeIDs(iv.first);
                       bitset<BLOCK_SIZE> block_bset = iv.second;
                       int cnt = (r-rids[0]) * this->k;
                       for (int c=cids[0]; c<cids[1]; c++) {
                            int temp = block_bset[BLOCK_SIZE-1-cnt];
                            temp *= dv[c];
                            outputv[r] += temp;
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

        void spgemm(k2tree* tree) { // with k2-tree
            // result storage: store back to original tree

            // multiply tree
            for (auto it : this->leafgroup) {
                // split rids
                std::vector<int> rids = tokenizeIDs(it.first);

                // for each rids
                if (rids[0] >= this->mat_height || rids[1] >= this->mat_height) continue;
                for (int i=rids[0]; i<=rids[1]; i++) { // for each rids
                    for (auto iv : (it.second)) { // for each cids
                       std::vector<int> cids = tokenizeIDs(iv.first);
                       bitset<BLOCK_SIZE> block_bset = iv.second;
                       int cnt = 0;

                       if (cids[0] >= this->mat_width || cids[1] >= this->mat_width) continue;
                       for (int j=cids[0]; j<=cids[1]; j++) {

                            /* bmm here */
                            if (tree->leafblockexist(it.first, iv.first)) {
                                // TODO: suppose to be bmm operation. Just to simulate it now.
                                /* (iv.second) = bmm();*/
                                (iv.second) ^= tree->leafgroup[it.first][iv.first];
                            }
                            else {(iv.second).reset();}
                       }
                    }
                }
            }

            // print out result for debug
            std::cout << "(temp) finish tree * tree !" << std::endl;
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

        bool leafblockexist(std::string rids, std::string cids) {
            if (this->leafgroup.find(rids) == this->leafgroup.end()) return false;
            else if (this->leafgroup[rids].find(cids) == this->leafgroup[rids].end()) return false;
            return true;
        }

};

#endif