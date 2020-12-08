__global__ void smoothing_vec( double * x1, const double * x2, const double * x3 , const int numx)
{
    int index = blockIdx.x + blockIdx.y * gridDim.x 
                + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    while (index < numx)
    {
        x1[index] = x1[index] * x2[index] + x3[index];
        index = index + num_threads;
    }
}