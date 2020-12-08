__global__ void add_vec( double * v1, const double * v2 ) 
{
    int idx = threadIdx.x;
    v1[idx] += v2[idx];
}