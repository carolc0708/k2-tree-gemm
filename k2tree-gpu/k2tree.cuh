#ifndef K2TREE_CUH
#define K2TREE_CUH

// the row-by-row model
__device__ __inline__ void MatMul(MulAParam* p)
{
    GET_LANEID;
    const int gdx = (CEIL(p->mat_length));
    const int gdy = (CEIL(p->mat_length));
    for (int bid=blockIdx.x*32+warpid; bid<gdx*gdy; bid+=gridDim.x*32)
    {
        unsigned bx = bid / gdy;
        unsigned by = bid % gdy;
        const unsigned* input_sub = &(p->input_gpu[bx*32]);
        const unsigned* weight_sub = &(p->A_gpu[by*32]); // each col of weight (A)
        unsigned* output_sub = &(p->output_gpu[by*gdx*32+bx*32]);
        //
        register int Cm[32] = {0};
        for (int i=0; (i*32)<(p->mat_length); i++)
        {
            unsigned r0 = input_sub[i*32*gdx+laneid];
            unsigned r1 = weight_sub[i*32*gdy+laneid];

            #pragma unroll
            for (int j=0; j<32; j++)
            {
                unsigned r2 = __shfl(r1, j); //from lane-j, r1 of weight matrix
                Cm[j] += __popc(r0 & r2);
            }
        }
        //
        unsigned C = 0;
        if ((bx*32+laneid)<(p->mat_length))
        {
            for (int i=0; i<32; i++)
            {
                if (by*32+i<(p->mat_length))
                {
                    C = (C << 1) | (Cm[i] > 0);
                }
            }
        }
        output_sub[laneid] = C;
    }
}

#endif