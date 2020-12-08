__global__ void calc_average2(double * Xi1, const double * Xi2, const double * V, const int Nx1, const int Nx2, const double dx1, const double dx2, const double minx1, const double minx2, const int * corners_index)
{
    int index = blockIdx.x + blockIdx.y * gridDim.x 
                + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    
    while (index < Nx1)
    {
        Xi1[index] = (Xi1[index] - minx1) / dx1;
        double xi2 = (Xi2[index] - minx2) / dx2;
        
        int loc1 = min(Nx1 - 2, int(floor(Xi1[index])));
        int loc2 = min(Nx2 - 2, int(floor(xi2)));
        
        Xi1[index] = Xi1[index] - loc1; // Weights from the 0 corner
        xi2        = xi2 - loc2;
        
        int scalar_loc = loc1 + Nx1 * loc2;
        
        Xi1[index] = (1 - Xi1[index]) * ((1 - xi2) * V[scalar_loc + corners_index[0]]
                                             + xi2 * V[scalar_loc + corners_index[1]])
                         + Xi1[index] * ((1 - xi2) * V[scalar_loc + corners_index[2]]
                                             + xi2 * V[scalar_loc + corners_index[3]]);
        index = index + num_threads;
    }
        
}