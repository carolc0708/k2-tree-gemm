
#include <stdio.h>
#include <assert.h>
#include "k2tree.h"
#include "k2tree.cuh"

using namespace std;

__global__ void bfs(MulAParam* p)
{
    //grid_group grid = this_grid();
    MatMul(p);
    //grid.sync();

}

int main()
{
    // config input
    int dev = 0;
    cudaSetDevice(dev);
    const unsigned mat_length = 128;

    printf("feil: %u \n", FEIL(mat_length));
    printf("ceil: %u \n", CEIL(mat_length));

    // set input vector/matrix
    float* f = NULL;
    float* f_gpu = NULL;

    // f
    srand(3333);
    cudaMallocHost((void**) &f, sizeof(unsigned)*mat_length*mat_length);

    for (int i = 0; i < mat_length; ++i) {
        for (int j = 0; j < mat_length; ++j) {
            f[i * mat_length + j] = rand() % 2;
        }
    }

    // print f
    printf(" \n input: \n");
    for (int i = 0; i < mat_length; ++i) {
        for (int j = 0; j < mat_length; ++j) {
            printf("%d", (f[i * mat_length + j])>0?1:0);
        }
        printf("\n");
    }

//    // print f transpose
//    printf("\n f transpose: \n");
//    for (int i = 0; i < mat_length; ++i) {
//        for (int j = 0; j < mat_length; ++j) {
//            printf("%d", (f[j * mat_length + i])>0?1:0);
//        }
//        printf("\n");
//    }

    // set iterations
    MulAParam* p = new MulAParam(mat_length);
    FILE* config_file = fopen("matrix.csv", "r"); // read in A (128 x 128) // better keep the content as 1D array
    MulAParam* p_gpu = p->initialize(config_file, f);

    // setup kernel
    int numThreads = 1024;
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties (&deviceProp, dev);
    int numBlocksPerSm;
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, bfs, numThreads, 0);
    void* args[] = {&p_gpu};


    START_TIMER;
    cudaLaunchCooperativeKernel((void*)bfs, numBlocksPerSm*deviceProp.multiProcessorCount, numThreads, args);
    STOP_TIMER;

    // output
    float* output = p->download_full_output();
    validate(output, mat_length); // only printout for now


    // release
    delete p;
    cudaFreeHost(f);

    return 0;
}