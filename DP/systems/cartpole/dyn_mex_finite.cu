/*
 * Example of how to use the mxGPUArray API in a MEX file.  This example shows
 * how to write a MEX function that takes a gpuArray input and returns a
 * gpuArray output, e.g. B=mexFunction(A).
 *
 * Copyright 2012 The MathWorks, Inc.
 */

#include "mex.h"
#include "matrix.h"
#include "gpu/mxGPUArray.h"
#include <iostream>

 /*
  * Device code
  */

__global__ void dyn3_mex_continuous(double* const d_dx, // outputs
    double* const x, double* const dx,
    double* const th, double* const dth, // input states
    double const* const in1, double const* const in2, // input actions
    const double mc, const double mp,
    const double l, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2]
                       * grid_size[3] * grid_size[4];
    double x3, x4, u1, u2;

    int index3 = (grid_size[2] > 1) ? index : 0;
    int index4 = (grid_size[3] > 1) ? index : 0;

    int uindex1 = (active_actions[0] == 1) ? index : 0;
    int uindex2 = (active_actions[1] == 1) ? index : 0;

    while (index < num_elements)
    {
        x3 = th[index3];
        x4 = dth[index4];

        u1 = in1[uindex1];
        u2 = in2[uindex2];

        double t2 = cos(x3);
        double t3 = sin(x3);
        double t4 = x4 * x4;
        double t6 = 1.0 / l;
        double t5 = t3 * t3;
        double t7 = mp * t5;
        double t8 = mc + t7;
        double t9 = 1.0 / t8;

        d_dx[index] = t6*t9*(l*u1 - t2*u2 + mp*t3*t4*1.0/(t6*t6) + g*l*mp*t2*t3);

        index = index + num_threads;

        index3 = (grid_size[2] > 1) ? index : 0;
        index4 = (grid_size[3] > 1) ? index : 0;

        uindex1 = (active_actions[0] == 1) ? index : 0;
        uindex2 = (active_actions[1] == 1) ? index : 0;
    }
}

__global__ void dyn4_mex_continuous(double* const d_dth, // outputs
    double* const x, double* const dx,
    double* const th, double* const dth, // input states
    double const* const in1, double const* const in2, // input actions
    const double mc, const double mp,
    const double l, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2]
                       * grid_size[3] * grid_size[4];
    double x3, x4, u1, u2;

    int index3 = (grid_size[2] > 1) ? index : 0;
    int index4 = (grid_size[3] > 1) ? index : 0;

    int uindex1 = (active_actions[0] == 1) ? index : 0;
    int uindex2 = (active_actions[1] == 1) ? index : 0;

    while (index < num_elements)
    {
        x3 = th[index3];
        x4 = dth[index4];

        u1 = in1[uindex1];
        u2 = in2[uindex2];

        double t2 = cos(x3);
        double t3 = sin(x3);
        double t4 = x4 * x4;
        double t6 = 1.0 / l;
        double t5 = t3 * t3;
        double t7 = mp * t5;
        double t8 = mc + t7;
        double t9 = 1.0 / t8;

        d_dth[index] = -t6*t9*(t2*u1 - t6*u2*(mc/mp + 1.0)) - t6*t9*(g*t3*(mc + mp) + l*mp*t2*t3*t4);

        index = index + num_threads;

        index3 = (grid_size[2] > 1) ? index : 0;
        index4 = (grid_size[3] > 1) ? index : 0;

        uindex1 = (active_actions[0] == 1) ? index : 0;
        uindex2 = (active_actions[1] == 1) ? index : 0;
    }
}

__global__ void step(double* const x1n, // outputs
                     double* const x1, // current states
                     double* const d_x1, // current derivatives
                     const double delta, int32_t const* const grid_size)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2]
                       * grid_size[3] * grid_size[4];
    
    while (index < num_elements)
    {
        x1n[index] = x1[index] + delta * d_x1[index];
        index = index + num_threads;
    }
}

__global__ void rk4_wbounds(double* const x1n, // outputs
                    double* const x1, // current states
                    double* const k1_x1, // K1
                    double* const k2_x1, // K2
                    double* const k3_x1, // K3
                    double* const k4_x1, // K4
                    double dt, double lower_limit, double upper_limit, int32_t const* const grid_size)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2]
                       * grid_size[3] * grid_size[4];
    while (index < num_elements)
    {
        x1n[index] = x1[index] + dt / 6.0 * (k1_x1[index] + 2.0 * k2_x1[index] + 2.0 * k3_x1[index] + k4_x1[index]);        
        x1n[index] = max(lower_limit, min(upper_limit, x1n[index]));
        index = index + num_threads;
    }
}
/*
 * Host code
 */
void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, mxArray const* prhs[])
{
    /* Declare all variables*/
    // Inputs
    // GPU Arrays
    mxGPUArray const* X1;
    mxGPUArray const* X2;
    mxGPUArray const* dX1;
    mxGPUArray const* dX2;
    mxGPUArray const* in1;
    mxGPUArray const* in2;
    mxGPUArray const* grid_size;
    mxGPUArray const* active_actions;
    // Pointers for GPU Arrays
    double* p_X1; 
    double* p_X2; 
    double* p_dX1;
    double* p_dX2;
    double const* p_in1; 
    double const* p_in2; 
    int32_t const* p_grid_size;
    int32_t const* p_active_actions;
    // Pointers for normal inputs
    double* p_m;
    double* p_l;
    double* p_g;
    double* p_dt;
    double* p_limits;
    double* p_x_dims_free;

    // Intermediate variables
    mxGPUArray* k1_dX1;
    mxGPUArray* k1_dX2;
    mxGPUArray* k2_X1;
    mxGPUArray* k2_X2;
    mxGPUArray* k2_dX1;
    mxGPUArray* k2_dX2;
    mxGPUArray* k3_X1;
    mxGPUArray* k3_X2;
    mxGPUArray* k3_dX1;
    mxGPUArray* k3_dX2;
    mxGPUArray* k4_X1;
    mxGPUArray* k4_X2;
    mxGPUArray* k4_dX1;
    mxGPUArray* k4_dX2;
    // Pointers for intermediate variables
    double* p_k1_X1;
    double* p_k1_X2;
    double* p_k1_dX1;
    double* p_k1_dX2;
    double* p_k2_X1;
    double* p_k2_X2;
    double* p_k2_dX1;
    double* p_k2_dX2;
    double* p_k3_X1;
    double* p_k3_X2;
    double* p_k3_dX1;
    double* p_k3_dX2;
    double* p_k4_X1;
    double* p_k4_X2;
    double* p_k4_dX1;
    double* p_k4_dX2;

    // Outputs
    mxGPUArray* X1n;
    mxGPUArray* X2n;
    mxGPUArray* dX1n;
    mxGPUArray* dX2n;
    // Pointers for outputs
    double* p_X1n;
    double* p_X2n;
    double* p_dX1n;
    double* p_dX2n;

    char const* const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const* const errMsg = "Invalid input to MEX file.";

    /* Choose a reasonably sized number of threads for the block. */
    int const threadsPerBlock = 256;
    int const blocksPerGrid = 1024;

    /* Initialize the MathWorks GPU API. */
    mxInitGPU();

    /* Throw an error if the input is not a GPU array. */
    if ((nrhs != 16) || !(mxIsGPUArray(prhs[0])) || !(mxIsGPUArray(prhs[1])) || !(mxIsGPUArray(prhs[2])) 
                     || !(mxIsGPUArray(prhs[3])) || !(mxIsGPUArray(prhs[4])) || !(mxIsGPUArray(prhs[5]))
                     || !(mxIsGPUArray(prhs[6]))) {
        mexErrMsgIdAndTxt(errId, errMsg);
    }

    X1 = mxGPUCreateFromMxArray(prhs[0]);
    X2 = mxGPUCreateFromMxArray(prhs[2]);
    dX1 = mxGPUCreateFromMxArray(prhs[1]);
    dX2 = mxGPUCreateFromMxArray(prhs[3]);
    in1 = mxGPUCreateFromMxArray(prhs[5]);
    in2 = mxGPUCreateFromMxArray(prhs[6]);
    grid_size = mxGPUCreateFromMxArray(prhs[12]);
    active_actions = mxGPUCreateFromMxArray(prhs[13]);

    p_m = mxGetDoubles(prhs[7]); 
    double const mc = p_m[0];
    p_m = mxGetDoubles(prhs[8]); 
    double const mp = p_m[0];

    p_l = mxGetDoubles(prhs[9]); 
    double const l = p_l[0];

    p_g = mxGetDoubles(prhs[10]);
    double const g = p_g[0];
    
    p_dt = mxGetDoubles(prhs[11]);
    double const dt = p_dt[0];
    
    p_limits = mxGetDoubles(prhs[14]);
    
    p_x_dims_free = mxGetDoubles(prhs[15]);
    mwSize const* num_x_dims = mxGetDimensions(prhs[15]);
    if ((mxGetNumberOfDimensions(prhs[15]) != 2) || (num_x_dims[1] > 1)) {
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    
    /*
     * Verify that inputs are of appropriate type before extracting the pointer.
     */
    if ((mxGPUGetClassID(X1) != mxDOUBLE_CLASS) || (mxGPUGetClassID(X2) != mxDOUBLE_CLASS)
        || (mxGPUGetClassID(dX1) != mxDOUBLE_CLASS) || (mxGPUGetClassID(dX2) != mxDOUBLE_CLASS)     
        || (mxGPUGetClassID(in1) != mxDOUBLE_CLASS) || (mxGPUGetClassID(in2) != mxDOUBLE_CLASS)
        || (mxGPUGetClassID(grid_size) != mxINT32_CLASS) || (mxGPUGetClassID(active_actions) != mxINT32_CLASS)) {
        mexErrMsgIdAndTxt(errId, errMsg);
    }

    /*
     * Now that we have verified the data type, extract a pointer to the input
     * data on the device.
     */
    p_X1 = (double*)(mxGPUGetDataReadOnly(X1));
    p_X2 = (double*)(mxGPUGetDataReadOnly(X2));
    p_dX1 = (double*)(mxGPUGetDataReadOnly(dX1));
    p_dX2 = (double*)(mxGPUGetDataReadOnly(dX2));
    p_in1 = (double const*)(mxGPUGetDataReadOnly(in1));
    p_in2 = (double const*)(mxGPUGetDataReadOnly(in2));
    p_grid_size = (int32_t const*)(mxGPUGetDataReadOnly(grid_size));
    p_active_actions = (int32_t const*)(mxGPUGetDataReadOnly(active_actions));
    
    /* Create output arrays*/
    X1n = mxGPUCopyGPUArray(X1);
    p_X1n = (double*)(mxGPUGetData(X1n));
    
    X2n = mxGPUCopyGPUArray(X2);
    p_X2n = (double*)(mxGPUGetData(X2n));
    
    dX1n = mxGPUCopyGPUArray(dX1);
    p_dX1n = (double*)(mxGPUGetData(dX1n));
    
    dX2n = mxGPUCopyGPUArray(dX2);
    p_dX2n = (double*)(mxGPUGetData(dX2n));
    
    // RK4 - Step 1
    int32_t curr_free_dim;
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            p_k1_X1 = p_dX1;
            step << <blocksPerGrid, threadsPerBlock >> > (p_X1n, p_X1, p_k1_X1, 0.5 * dt, p_grid_size);
    
        } else if (curr_free_dim == 3) {
            p_k1_X2 = p_dX2;
            step << <blocksPerGrid, threadsPerBlock >> > (p_X2n, p_X2, p_k1_X2, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 2) {
            k1_dX1 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dX1), mxGPUGetDimensions(dX1), mxGPUGetClassID(dX1), mxGPUGetComplexity(dX1), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_dX1 = (double*)(mxGPUGetData(k1_dX1));
            dyn3_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_dX1, p_X1, p_dX1, p_X2, p_dX2, 
                                                                         p_in1, p_in2, mc, mp, l, 
                                                                         g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_dX1n, p_dX1, p_k1_dX1, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 4) {
            k1_dX2 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dX2), mxGPUGetDimensions(dX2), mxGPUGetClassID(dX2), mxGPUGetComplexity(dX2), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_dX2 = (double*)(mxGPUGetData(k1_dX2));
            dyn4_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_dX2, p_X1, p_dX1, p_X2, p_dX2, 
                                                                         p_in1, p_in2, mc, mp, l, 
                                                                         g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_dX2n, p_dX2, p_k1_dX2, 0.5 * dt, p_grid_size);   
        }
    }
    
    // RK4 - Step 2
    // Continuous Dynamics
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            k2_X1 = mxGPUCopyGPUArray(dX1n);
            p_k2_X1 = (double*)(mxGPUGetData(k2_X1));
    
        } else if (curr_free_dim == 3) {
            k2_X2 = mxGPUCopyGPUArray(dX2n);
            p_k2_X2 = (double*)(mxGPUGetData(k2_X2));
            
        } else if (curr_free_dim == 2) {
            k2_dX1 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dX1), mxGPUGetDimensions(dX1), mxGPUGetClassID(dX1), mxGPUGetComplexity(dX1), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_dX1 = (double*)(mxGPUGetData(k2_dX1));
            dyn3_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_dX1, p_X1n, p_dX1n, p_X2n, p_dX2n, 
                                                                         p_in1, p_in2, mc, mp, l, 
                                                                         g, p_grid_size, p_active_actions);
        } else if (curr_free_dim == 4) {
            k2_dX2 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dX2), mxGPUGetDimensions(dX2), mxGPUGetClassID(dX2), mxGPUGetComplexity(dX2), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_dX2 = (double*)(mxGPUGetData(k2_dX2));
            dyn4_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_dX2, p_X1n, p_dX1n, p_X2n, p_dX2n, 
                                                                         p_in1, p_in2, mc, mp, l, 
                                                                         g, p_grid_size, p_active_actions);            
        }
    }
    // Discrete step
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_X1n, p_X1, p_k2_X1, 0.5 * dt, p_grid_size);
    
        } else if (curr_free_dim == 3) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_X2n, p_X2, p_k2_X2, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 2) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dX1n, p_dX1, p_k2_dX1, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 4) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dX2n, p_dX2, p_k2_dX2, 0.5 * dt, p_grid_size);
            
        } 
    }
    
    // RK4 - Step 3
    // Continuous Dynamics
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            k3_X1 = mxGPUCopyGPUArray(dX1n);
            p_k3_X1 = (double*)(mxGPUGetData(k3_X1));
    
        } else if (curr_free_dim == 3) {
            k3_X2 = mxGPUCopyGPUArray(dX2n);
            p_k3_X2 = (double*)(mxGPUGetData(k3_X2));
            
        } else if (curr_free_dim == 2) {
            k3_dX1 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dX1), mxGPUGetDimensions(dX1), mxGPUGetClassID(dX1), mxGPUGetComplexity(dX1), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_dX1 = (double*)(mxGPUGetData(k3_dX1));
            dyn3_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_dX1, p_X1n, p_dX1n, p_X2n, p_dX2n, 
                                                                         p_in1, p_in2, mc, mp, l,
                                                                         g, p_grid_size, p_active_actions);
        } else if (curr_free_dim == 4) {
            k3_dX2 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dX2), mxGPUGetDimensions(dX2), mxGPUGetClassID(dX2), mxGPUGetComplexity(dX2), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_dX2 = (double*)(mxGPUGetData(k3_dX2));
            dyn4_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_dX2, p_X1n, p_dX1n, p_X2n, p_dX2n,
                                                                         p_in1, p_in2, mc, mp, l,
                                                                         g, p_grid_size, p_active_actions);            
        }
    }
    // Discrete step
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_X1n, p_X1, p_k3_X1, dt, p_grid_size);
    
        } else if (curr_free_dim == 3) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_X2n, p_X2, p_k3_X2, dt, p_grid_size);
            
        } else if (curr_free_dim == 2) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dX1n, p_dX1, p_k3_dX1, dt, p_grid_size);
            
        } else if (curr_free_dim == 4) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dX2n, p_dX2, p_k3_dX2, dt, p_grid_size);
            
        }
    }
    
    // RK4 - Step 4
    // Continuous Dynamics
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            k4_X1 = mxGPUCopyGPUArray(dX1n);
            p_k4_X1 = (double*)(mxGPUGetData(k4_X1));
    
        } else if (curr_free_dim == 3) {
            k4_X2 = mxGPUCopyGPUArray(dX2n);
            p_k4_X2 = (double*)(mxGPUGetData(k4_X2));
            
        } else if (curr_free_dim == 2) {
            k4_dX1 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dX1), mxGPUGetDimensions(dX1), mxGPUGetClassID(dX1), mxGPUGetComplexity(dX1), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_dX1 = (double*)(mxGPUGetData(k4_dX1));
            dyn3_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_dX1, p_X1n, p_dX1n, p_X2n, p_dX2n, 
                                                                         p_in1, p_in2, mc, mp, l, 
                                                                         g, p_grid_size, p_active_actions);
        } else if (curr_free_dim == 4) {
            k4_dX2 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dX2), mxGPUGetDimensions(dX2), mxGPUGetClassID(dX2), mxGPUGetComplexity(dX2), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_dX2 = (double*)(mxGPUGetData(k4_dX2));
            dyn4_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_dX2, p_X1n, p_dX1n, p_X2n, p_dX2n, 
                                                                         p_in1, p_in2, mc, mp, l, 
                                                                         g, p_grid_size, p_active_actions);            
        }
    }

    // RK4 Ouput
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_X1n, p_X1, p_k1_X1, p_k2_X1, p_k3_X1, p_k4_X1, dt, p_limits[0], p_limits[4], p_grid_size);
            mxGPUDestroyGPUArray(k2_X1);
            mxGPUDestroyGPUArray(k3_X1);
            mxGPUDestroyGPUArray(k4_X1);
        } else if (curr_free_dim == 3) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_X2n, p_X2, p_k1_X2, p_k2_X2, p_k3_X2, p_k4_X2, dt, p_limits[2], p_limits[6], p_grid_size);
            mxGPUDestroyGPUArray(k2_X2);
            mxGPUDestroyGPUArray(k3_X2);
            mxGPUDestroyGPUArray(k4_X2);
        } else if (curr_free_dim == 2) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_dX1n, p_dX1, p_k1_dX1, p_k2_dX1, p_k3_dX1, p_k4_dX1, dt, p_limits[1], p_limits[5], p_grid_size);
            mxGPUDestroyGPUArray(k1_dX1);
            mxGPUDestroyGPUArray(k2_dX1);
            mxGPUDestroyGPUArray(k3_dX1);
            mxGPUDestroyGPUArray(k4_dX1);
        } else if (curr_free_dim == 4) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_dX2n, p_dX2, p_k1_dX2, p_k2_dX2, p_k3_dX2, p_k4_dX2, dt, p_limits[3], p_limits[7], p_grid_size);
            mxGPUDestroyGPUArray(k1_dX2);
            mxGPUDestroyGPUArray(k2_dX2);
            mxGPUDestroyGPUArray(k3_dX2);
            mxGPUDestroyGPUArray(k4_dX2);
        }
    }
 
    /* Wrap the result up as a MATLAB gpuArray for return. */
    plhs[0] = mxGPUCreateMxArrayOnGPU(X1n);
    plhs[2] = mxGPUCreateMxArrayOnGPU(X2n);
    plhs[1] = mxGPUCreateMxArrayOnGPU(dX1n);
    plhs[3] = mxGPUCreateMxArrayOnGPU(dX2n);

    /*
     * The mxGPUArray pointers are host-side structures that refer to device
     * data. These must be destroyed before leaving the MEX function.
     */
    mxGPUDestroyGPUArray(X1);
    mxGPUDestroyGPUArray(X2);
    mxGPUDestroyGPUArray(dX1);
    mxGPUDestroyGPUArray(dX2);
    mxGPUDestroyGPUArray(in1);
    mxGPUDestroyGPUArray(in2);
    mxGPUDestroyGPUArray(grid_size);
    mxGPUDestroyGPUArray(active_actions);
    mxGPUDestroyGPUArray(X1n);
    mxGPUDestroyGPUArray(X2n);
    mxGPUDestroyGPUArray(dX1n);
    mxGPUDestroyGPUArray(dX2n);
}
