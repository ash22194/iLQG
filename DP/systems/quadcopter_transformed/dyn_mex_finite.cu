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

__global__ void dyn1_mex_continuous(double* const d_z, // outputs
    double* const z, double* const roll, double* const pitch, double* const yaw,
    double* const vx, double* const vy, double* const vz, double* const vroll, double* const vpitch, double* const vyaw, // input states
    double const* const in1, double const* const in2, double const* const in3, double const* const in4, // input actions
    const double m, const double I1, const double I2, const double I3, 
    const double l, const double bk, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3] * grid_size[4] 
                       * grid_size[5] * grid_size[6] * grid_size[7] * grid_size[8] * grid_size[9];
    double dz;
    int index7 = (grid_size[6] > 1) ? index : 0;

    while (index < num_elements)
    {
        dz = vz[index7];
        d_z[index] = dz;
        
        index = index + num_threads;
        index7 = (grid_size[6] > 1) ? index : 0;
    }
}

__global__ void dyn2_mex_continuous(double* const d_roll, // outputs
    double* const z, double* const roll, double* const pitch, double* const yaw,
    double* const vx, double* const vy, double* const vz, double* const vroll, double* const vpitch, double* const vyaw, // input states
    double const* const in1, double const* const in2, double const* const in3, double const* const in4, // input actions
    const double m, const double I1, const double I2, const double I3, 
    const double l, const double bk, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3] * grid_size[4] 
                       * grid_size[5] * grid_size[6] * grid_size[7] * grid_size[8] * grid_size[9];
    double droll;
    int index8 = (grid_size[7] > 1) ? index : 0;

    while (index < num_elements)
    {
        droll = vroll[index8];
        d_roll[index] = droll;
        
        index = index + num_threads;
        index8 = (grid_size[7] > 1) ? index : 0;
    }
}

__global__ void dyn3_mex_continuous(double* const d_pitch, // outputs
    double* const z, double* const roll, double* const pitch, double* const yaw,
    double* const vx, double* const vy, double* const vz, double* const vroll, double* const vpitch, double* const vyaw, // input states
    double const* const in1, double const* const in2, double const* const in3, double const* const in4, // input actions
    const double m, const double I1, const double I2, const double I3, 
    const double l, const double bk, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3] * grid_size[4] 
                       * grid_size[5] * grid_size[6] * grid_size[7] * grid_size[8] * grid_size[9];
    double dpitch;
    int index9 = (grid_size[8] > 1) ? index : 0;

    while (index < num_elements)
    {
        dpitch = vpitch[index9];
        d_pitch[index] = dpitch;
        
        index = index + num_threads;
        index9 = (grid_size[8] > 1) ? index : 0;
    }
}

__global__ void dyn4_mex_continuous(double* const d_yaw, // outputs
    double* const z, double* const roll, double* const pitch, double* const yaw,
    double* const vx, double* const vy, double* const vz, double* const vroll, double* const vpitch, double* const vyaw, // input states
    double const* const in1, double const* const in2, double const* const in3, double const* const in4, // input actions
    const double m, const double I1, const double I2, const double I3, 
    const double l, const double bk, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3] * grid_size[4] 
                       * grid_size[5] * grid_size[6] * grid_size[7] * grid_size[8] * grid_size[9];
    double dyaw;
    int index10 = (grid_size[9] > 1) ? index : 0;

    while (index < num_elements)
    {
        dyaw = vyaw[index10];
        d_yaw[index] = dyaw;
        
        index = index + num_threads;
        index10 = (grid_size[9] > 1) ? index : 0;
    }
}

__global__ void dyn5_mex_continuous(double* const d_dx, // outputs
    double* const z, double* const roll, double* const pitch, double* const yaw,
    double* const vx, double* const vy, double* const vz, double* const vroll, double* const vpitch, double* const vyaw, // input states
    double const* const in1, double const* const in2, double const* const in3, double const* const in4, // input actions
    const double m, const double I1, const double I2, const double I3, 
    const double l, const double bk, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3] * grid_size[4] 
                       * grid_size[5] * grid_size[6] * grid_size[7] * grid_size[8] * grid_size[9];
    
    double x2, x3, x4, u1;

    int index2  = (grid_size[1] > 1) ? index : 0;
    int index3  = (grid_size[2] > 1) ? index : 0;
    int index4  = (grid_size[3] > 1) ? index : 0;

    int uindex1 = (active_actions[0] == 1) ? index : 0;

    while (index < num_elements)
    {
        x2  = roll[index2];
        x3  = pitch[index3];
        x4  = yaw[index4];
        
        u1 = in1[uindex1];
        
        double t2 = cos(x2);
        double t4 = cos(x4);
        double t5 = sin(x2);
        double t6 = sin(x3);
        double t7 = sin(x4);
        double t23 = 1.0 / m;
        
        d_dx[index] = t23 * u1 * (t5 * t7 + t2 * t4 * t6);
        
        index = index + num_threads;
        
        index2  = (grid_size[1] > 1) ? index : 0;
        index3  = (grid_size[2] > 1) ? index : 0;
        index4  = (grid_size[3] > 1) ? index : 0;

        uindex1 = (active_actions[0] == 1) ? index : 0;
    }
}

__global__ void dyn6_mex_continuous(double* const d_dy, // outputs
    double* const z, double* const roll, double* const pitch, double* const yaw,
    double* const vx, double* const vy, double* const vz, double* const vroll, double* const vpitch, double* const vyaw, // input states
    double const* const in1, double const* const in2, double const* const in3, double const* const in4, // input actions
    const double m, const double I1, const double I2, const double I3, 
    const double l, const double bk, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3] * grid_size[4] 
                       * grid_size[5] * grid_size[6] * grid_size[7] * grid_size[8] * grid_size[9];
    
    double x2, x3, x4, u1;

    int index2  = (grid_size[1] > 1) ? index : 0;
    int index3  = (grid_size[2] > 1) ? index : 0;
    int index4  = (grid_size[3] > 1) ? index : 0;

    int uindex1 = (active_actions[0] == 1) ? index : 0;

    while (index < num_elements)
    {
        x2  = roll[index2];
        x3  = pitch[index3];
        x4  = yaw[index4];
        
        u1 = in1[uindex1];
        
        double t2 = cos(x2);
        double t4 = cos(x4);
        double t5 = sin(x2);
        double t6 = sin(x3);
        double t7 = sin(x4);
        double t23 = 1.0 / m;
        
        d_dy[index] = -t23 * u1 * (t4 * t5 - t2 * t6 * t7);
        
        index = index + num_threads;
        
        index2  = (grid_size[1] > 1) ? index : 0;
        index3  = (grid_size[2] > 1) ? index : 0;
        index4  = (grid_size[3] > 1) ? index : 0;

        uindex1 = (active_actions[0] == 1) ? index : 0;
    }
}

__global__ void dyn7_mex_continuous(double* const d_dz, // outputs
    double* const z, double* const roll, double* const pitch, double* const yaw,
    double* const vx, double* const vy, double* const vz, double* const vroll, double* const vpitch, double* const vyaw, // input states
    double const* const in1, double const* const in2, double const* const in3, double const* const in4, // input actions
    const double m, const double I1, const double I2, const double I3, 
    const double l, const double bk, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3] * grid_size[4] 
                       * grid_size[5] * grid_size[6] * grid_size[7] * grid_size[8] * grid_size[9];
    
    double x2, x3, u1;

    int index2  = (grid_size[1] > 1) ? index : 0;
    int index3  = (grid_size[2] > 1) ? index : 0;

    int uindex1 = (active_actions[0] == 1) ? index : 0;

    while (index < num_elements)
    {
        x2  = roll[index2];
        x3  = pitch[index3];
        
        u1 = in1[uindex1];
        
        double t2 = cos(x2);
        double t3 = cos(x3);
        double t23 = 1.0 / m;
        
        d_dz[index] = -g + t2 * t3 * t23 * u1;
        
        index = index + num_threads;
        
        index2  = (grid_size[1] > 1) ? index : 0;
        index3  = (grid_size[2] > 1) ? index : 0;

        uindex1 = (active_actions[0] == 1) ? index : 0;
    }
}

__global__ void dyn8_mex_continuous(double* const d_droll, // outputs
    double* const z, double* const roll, double* const pitch, double* const yaw,
    double* const vx, double* const vy, double* const vz, double* const vroll, double* const vpitch, double* const vyaw, // input states
    double const* const in1, double const* const in2, double const* const in3, double const* const in4, // input actions
    const double m, const double I1, const double I2, const double I3, 
    const double l, const double bk, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3] * grid_size[4] 
                       * grid_size[5] * grid_size[6] * grid_size[7] * grid_size[8] * grid_size[9];
    
    double x2, x3, x4, x8, x9, x10, u2, u3, u4;

    int index2  = (grid_size[1] > 1) ? index : 0;
    int index3  = (grid_size[2] > 1) ? index : 0;
    int index4  = (grid_size[3] > 1) ? index : 0;
    int index8  = (grid_size[7] > 1) ? index : 0;
    int index9  = (grid_size[8] > 1) ? index : 0;
    int index10 = (grid_size[9] > 1) ? index : 0;

    int uindex2 = (active_actions[1] == 1) ? index : 0;
    int uindex3 = (active_actions[2] == 1) ? index : 0;
    int uindex4 = (active_actions[3] == 1) ? index : 0;

    while (index < num_elements)
    {
        x2  = roll[index2];
        x3  = pitch[index3];
        x4  = yaw[index4];
        x8  = vroll[index8];
        x9  = vpitch[index9];
        x10 = vyaw[index10];
        
        u2 = in2[uindex2];
        u3 = in3[uindex3];
        u4 = in4[uindex4];
        
        double t2 = cos(x2);
        double t3 = cos(x3);
        double t4 = cos(x4);
        double t5 = sin(x2);
        double t6 = sin(x3);
        double t7 = sin(x4);
        double t8 = I1 * I1;
        double t9 = I2 * I2;
        double t10 = I3 * I3;
        double t11 = x2 * 2.0;
        double t12 = x3 * 2.0;
        double t13 = x9 * x9;
        double t14 = x10 * x10;
        double t21 = 1.0 / I2;
        double t22 = 1.0 / I3;
        double t23 = 1.0 / m;
        double t15 = t2 * t2;
        double t16 = t3 * t3;
        double t17 = t3 * t16;
        double t18 = sin(t11);
        double t19 = sin(t12);
        double t20 = t5 * t5;
        double t24 = 1.0 / t3;
        double t25 = 1.0 / t16;
        
        double et1 = I1 * t10 * x9 * x10-I3 * t8 * x9 * x10+I1 * I2 * I3 * x9 * x10+I2 * I3 * l * t3 * u2-I1 * t6 * t10 * x8 * x9+I3 * t6 * t8 * x8 * x9+I1 * t9 * t15 * x9 * x10-I2 * t8 * t15 * x9 * x10-I1 * t10 * t15 * x9 * x10+I3 * t8 * t15 * x9 * x10-I1 * t10 * t16 * x9 * x10+I3 * t8 * t16 * x9 * x10+I2 * t10 * t16 * x9 * x10-I3 * t9 * t16 * x9 * x10+I1 * I2 * I3 * t16 * x8 * x9+I1 * I2 * bk * t2 * t6 * u4+I1 * I3 * l * t5 * t6 * u3+I1 * t2 * t3 * t5 * t9 * t14-I2 * t2 * t3 * t5 * t8 * t14-I1 * t2 * t3 * t5 * t10 * t14+I2 * t2 * t3 * t5 * t10 * t13+I3 * t2 * t3 * t5 * t8 * t14-I3 * t2 * t3 * t5 * t9 * t13-I1 * t2 * t5 * t9 * t14 * t17+I2 * t2 * t5 * t8 * t14 * t17+I1 * t2 * t5 * t10 * t14 * t17-I3 * t2 * t5 * t8 * t14 * t17;
        
        double et2 = -I2 * t2 * t5 * t10 * t14 * t17+I3 * t2 * t5 * t9 * t14 * t17-I1 * t6 * t9 * t15 * x8 * x9+I2 * t6 * t8 * t15 * x8 * x9+I1 * t6 * t10 * t15 * x8 * x9-I3 * t6 * t8 * t15 * x8 * x9-I1 * t9 * t15 * t16 * x9 * x10+I2 * t8 * t15 * t16 * x9 * x10+I1 * t10 * t15 * t16 * x9 * x10-I3 * t8 * t15 * t16 * x9 * x10-I2 * t10 * t15 * t16 * x9 * x10 * 2.0+I3 * t9 * t15 * t16 * x9 * x10 * 2.0+I1 * I2 * I3 * t6 * t15 * x8 * x9-I1 * I2 * I3 * t15 * t16 * x8 * x9-I1 * I2 * I3 * t2 * t5 * t17 * x8 * x10-I1 * t2 * t3 * t5 * t6 * t9 * x8 * x10+I2 * t2 * t3 * t5 * t6 * t8 * x8 * x10+I1 * t2 * t3 * t5 * t6 * t10 * x8 * x10-I3 * t2 * t3 * t5 * t6 * t8 * x8 * x10+I1 * I2 * I3 * t2 * t3 * t5 * t6 * x8 * x10;
        
        d_droll[index] = (t21 * t22 * t24 * (et1 + et2)) / I1;
        
        index = index + num_threads;
        
        index2  = (grid_size[1] > 1) ? index : 0;
        index3  = (grid_size[2] > 1) ? index : 0;
        index4  = (grid_size[3] > 1) ? index : 0;
        index8  = (grid_size[7] > 1) ? index : 0;
        index9  = (grid_size[8] > 1) ? index : 0;
        index10 = (grid_size[9] > 1) ? index : 0;

        uindex2 = (active_actions[1] == 1) ? index : 0;
        uindex3 = (active_actions[2] == 1) ? index : 0;
        uindex4 = (active_actions[3] == 1) ? index : 0;
    }
}

__global__ void dyn9_mex_continuous(double* const d_dpitch, // outputs
    double* const z, double* const roll, double* const pitch, double* const yaw,
    double* const vx, double* const vy, double* const vz, double* const vroll, double* const vpitch, double* const vyaw, // input states
    double const* const in1, double const* const in2, double const* const in3, double const* const in4, // input actions
    const double m, const double I1, const double I2, const double I3, 
    const double l, const double bk, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3] * grid_size[4] 
                       * grid_size[5] * grid_size[6] * grid_size[7] * grid_size[8] * grid_size[9];
    
    double x2, x3, x4, x8, x9, x10, u3, u4;

    int index2  = (grid_size[1] > 1) ? index : 0;
    int index3  = (grid_size[2] > 1) ? index : 0;
    int index4  = (grid_size[3] > 1) ? index : 0;
    int index8  = (grid_size[7] > 1) ? index : 0;
    int index9  = (grid_size[8] > 1) ? index : 0;
    int index10 = (grid_size[9] > 1) ? index : 0;

    int uindex3 = (active_actions[2] == 1) ? index : 0;
    int uindex4 = (active_actions[3] == 1) ? index : 0;

    while (index < num_elements)
    {
        x2  = roll[index2];
        x3  = pitch[index3];
        x4  = yaw[index4];
        x8  = vroll[index8];
        x9  = vpitch[index9];
        x10 = vyaw[index10];
        
        u3 = in3[uindex3];
        u4 = in4[uindex4];
        
        double t2 = cos(x2);
        double t3 = cos(x3);
        double t4 = cos(x4);
        double t5 = sin(x2);
        double t6 = sin(x3);
        double t7 = sin(x4);
        double t8 = I1 * I1;
        double t9 = I2 * I2;
        double t10 = I3 * I3;
        double t11 = x2 * 2.0;
        double t12 = x3 * 2.0;
        double t13 = x9 * x9;
        double t14 = x10 * x10;
        double t21 = 1.0 / I2;
        double t22 = 1.0 / I3;
        double t23 = 1.0 / m;
        double t15 = t2 * t2;
        double t16 = t3 * t3;
        double t17 = t3 * t16;
        double t18 = sin(t11);
        double t19 = sin(t12);
        double t20 = t5 * t5;
        double t24 = 1.0 / t3;
        double t25 = 1.0 / t16;
        
        double mt2 = t21 * t22 * ((t9 * t14 * t19) - (t3 * t9 * x8 * x10 * 2.0) - (t9 * t18 * x8 * x9) + (t10 * t18 * x8 * x9) - (I1 * I2 * t14 * t19) + (I2 * bk * t5 * u4 * 2.0) - (I3 * l * t2 * u3 * 2.0) + (I1 * I2 * t3 * x8 * x10 * 2.0) + (I2 * I3 * t3 * x8 * x10 * 2.0) + (I1 * I2 * t18 * x8 * x9) - (I1 * I3 * t18 * x8 * x9) - (t3 * t6 * t9 * t14 * t15 * 2.0) + (t3 * t6 * t10 * t14 * t15 * 2.0) + (t3 * t9 * t15 * x8 * x10 * 2.0) - (t3 * t10 * t15 * x8 * x10 * 2.0) + (t2 * t5 * t6 * t9 * x9 * x10 * 2.0) - (t2 * t5 * t6 * t10 * x9 * x10 * 2.0) + (I1 * I2 * t3 * t6 * t14 * t15 * 2.0) - (I1 * I3 * t3 * t6 * t14 * t15 * 2.0) - (I1 * I2 * t3 * t15 * x8 * x10 * 2.0) + (I1 * I3 * t3 * t15 * x8 * x10 * 2.0) - (I1 * I2 * t2 * t5 * t6 * x9 * x10 * 2.0) + (I1 * I3 * t2 * t5 * t6 * x9 * x10 * 2.0)) * (-0.5);
        
        d_dpitch[index] = mt2;
        
        index = index + num_threads;
        
        index2  = (grid_size[1] > 1) ? index : 0;
        index3  = (grid_size[2] > 1) ? index : 0;
        index4  = (grid_size[3] > 1) ? index : 0;
        index8  = (grid_size[7] > 1) ? index : 0;
        index9  = (grid_size[8] > 1) ? index : 0;
        index10 = (grid_size[9] > 1) ? index : 0;

        uindex3 = (active_actions[2] == 1) ? index : 0;
        uindex4 = (active_actions[3] == 1) ? index : 0;
    }
}

__global__ void dyn10_mex_continuous(double* const d_dyaw, // outputs
    double* const z, double* const roll, double* const pitch, double* const yaw,
    double* const vx, double* const vy, double* const vz, double* const vroll, double* const vpitch, double* const vyaw, // input states
    double const* const in1, double const* const in2, double const* const in3, double const* const in4, // input actions
    const double m, const double I1, const double I2, const double I3, 
    const double l, const double bk, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3] * grid_size[4] 
                       * grid_size[5] * grid_size[6] * grid_size[7] * grid_size[8] * grid_size[9];
    
    double x2, x3, x4, x8, x9, x10, u3, u4;

    int index2  = (grid_size[1] > 1) ? index : 0;
    int index3  = (grid_size[2] > 1) ? index : 0;
    int index4  = (grid_size[3] > 1) ? index : 0;
    int index8  = (grid_size[7] > 1) ? index : 0;
    int index9  = (grid_size[8] > 1) ? index : 0;
    int index10 = (grid_size[9] > 1) ? index : 0;

    int uindex3 = (active_actions[2] == 1) ? index : 0;
    int uindex4 = (active_actions[3] == 1) ? index : 0;

    while (index < num_elements)
    {
        x2  = roll[index2];
        x3  = pitch[index3];
        x4  = yaw[index4];
        x8  = vroll[index8];
        x9  = vpitch[index9];
        x10 = vyaw[index10];
        
        u3 = in3[uindex3];
        u4 = in4[uindex4];
        
        double t2 = cos(x2);
        double t3 = cos(x3);
        double t4 = cos(x4);
        double t5 = sin(x2);
        double t6 = sin(x3);
        double t7 = sin(x4);
        double t8 = I1 * I1;
        double t9 = I2 * I2;
        double t10 = I3 * I3;
        double t11 = x2 * 2.0;
        double t12 = x3 * 2.0;
        double t13 = x9 * x9;
        double t14 = x10 * x10;
        double t21 = 1.0 / I2;
        double t22 = 1.0 / I3;
        double t23 = 1.0 / m;
        double t15 = t2 * t2;
        double t16 = t3 * t3;
        double t17 = t3 * t16;
        double t18 = sin(t11);
        double t19 = sin(t12);
        double t20 = t5 * t5;
        double t24 = 1.0 / t3;
        double t25 = 1.0 / t16;
        
        double mt3 = t25 * (t5 * x9-t2 * t3 * x10) * (t3 * t5 * x8-t2 * t6 * x9)+t21 * t22 * t24 * (-t9 * t15 * x8 * x9-t10 * t20 * x8 * x9+I2 * bk * t2 * u4+I3 * l * t5 * u3+I1 * I2 * t15 * x8 * x9+I1 * I3 * t20 * x8 * x9+t6 * t9 * t15 * x9 * x10+t6 * t10 * t20 * x9 * x10+t2 * t3 * t5 * t6 * t9 * t14-t2 * t3 * t5 * t6 * t10 * t14-t2 * t3 * t5 * t9 * x8 * x10+t2 * t3 * t5 * t10 * x8 * x10-I1 * I2 * t6 * t15 * x9 * x10-I1 * I3 * t6 * t20 * x9 * x10-I1 * I2 * t2 * t3 * t5 * t6 * t14+I1 * I3 * t2 * t3 * t5 * t6 * t14+I1 * I2 * t2 * t3 * t5 * x8 * x10-I1 * I3 * t2 * t3 * t5 * x8 * x10)+t25 * x8 * cos(x2-x3) * (t2 * x9+t3 * t5 * x10);
        
        d_dyaw[index] = mt3;
        
        index = index + num_threads;
        
        index2  = (grid_size[1] > 1) ? index : 0;
        index3  = (grid_size[2] > 1) ? index : 0;
        index4  = (grid_size[3] > 1) ? index : 0;
        index8  = (grid_size[7] > 1) ? index : 0;
        index9  = (grid_size[8] > 1) ? index : 0;
        index10 = (grid_size[9] > 1) ? index : 0;

        uindex3 = (active_actions[2] == 1) ? index : 0;
        uindex4 = (active_actions[3] == 1) ? index : 0;
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
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3] * grid_size[4] 
                       * grid_size[5] * grid_size[6] * grid_size[7] * grid_size[8] * grid_size[9];
    
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
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3] * grid_size[4] 
                       * grid_size[5] * grid_size[6] * grid_size[7] * grid_size[8] * grid_size[9];
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
    /* Declare all variables. */
    // Inputs
    // GPU Arrays
    mxGPUArray const* Z;
    mxGPUArray const* Roll;
    mxGPUArray const* Pitch;
    mxGPUArray const* Yaw;
    mxGPUArray const* VX;
    mxGPUArray const* VY;
    mxGPUArray const* VZ;
    mxGPUArray const* VRoll;
    mxGPUArray const* VPitch;
    mxGPUArray const* VYaw;
    mxGPUArray const* in1;
    mxGPUArray const* in2;
    mxGPUArray const* in3;
    mxGPUArray const* in4;
    mxGPUArray const* grid_size;
    mxGPUArray const* active_actions;
    // Pointers for GPU Arrays
    double* p_Z; 
    double* p_Roll;
    double* p_Pitch;
    double* p_Yaw; 
    double* p_VX;
    double* p_VY;
    double* p_VZ;
    double* p_VRoll;
    double* p_VPitch;
    double* p_VYaw;
    double const* p_in1; 
    double const* p_in2;
    double const* p_in3;
    double const* p_in4; 
    int32_t const* p_grid_size;
    int32_t const* p_active_actions;
    // Pointers for normal inputs
    double* p_m;
    double* p_I;
    double* p_l;
    double* p_bk;
    double* p_g;
    double* p_dt;
    double* p_limits;
    double* p_x_dims_free;
    
    // Intermediate variables
//     mxGPUArray* k1_X1;
//     mxGPUArray* k1_X2;
//     mxGPUArray* k1_X3;
//     mxGPUArray* k1_X4;
    mxGPUArray* k1_VX;
    mxGPUArray* k1_VY;
    mxGPUArray* k1_VZ;
    mxGPUArray* k1_VRoll;
    mxGPUArray* k1_VPitch;
    mxGPUArray* k1_VYaw;
    mxGPUArray* k2_Z;
    mxGPUArray* k2_Roll;
    mxGPUArray* k2_Pitch;
    mxGPUArray* k2_Yaw;
    mxGPUArray* k2_VX;
    mxGPUArray* k2_VY;
    mxGPUArray* k2_VZ;
    mxGPUArray* k2_VRoll;
    mxGPUArray* k2_VPitch;
    mxGPUArray* k2_VYaw;
    mxGPUArray* k3_Z;
    mxGPUArray* k3_Roll;
    mxGPUArray* k3_Pitch;
    mxGPUArray* k3_Yaw;
    mxGPUArray* k3_VX;
    mxGPUArray* k3_VY;
    mxGPUArray* k3_VZ;
    mxGPUArray* k3_VRoll;
    mxGPUArray* k3_VPitch;
    mxGPUArray* k3_VYaw;
    mxGPUArray* k4_Z;
    mxGPUArray* k4_Roll;
    mxGPUArray* k4_Pitch;
    mxGPUArray* k4_Yaw;
    mxGPUArray* k4_VX;
    mxGPUArray* k4_VY;
    mxGPUArray* k4_VZ;
    mxGPUArray* k4_VRoll;
    mxGPUArray* k4_VPitch;
    mxGPUArray* k4_VYaw;
    // Pointers for intermediate variables
    double* p_k1_Z;
    double* p_k1_Roll;
    double* p_k1_Pitch;
    double* p_k1_Yaw;
    double* p_k1_VX;
    double* p_k1_VY;
    double* p_k1_VZ;
    double* p_k1_VRoll;
    double* p_k1_VPitch;
    double* p_k1_VYaw;
    double* p_k2_Z;
    double* p_k2_Roll;
    double* p_k2_Pitch;
    double* p_k2_Yaw;
    double* p_k2_VX;
    double* p_k2_VY;
    double* p_k2_VZ;
    double* p_k2_VRoll;
    double* p_k2_VPitch;
    double* p_k2_VYaw;
    double* p_k3_Z;
    double* p_k3_Roll;
    double* p_k3_Pitch;
    double* p_k3_Yaw;
    double* p_k3_VX;
    double* p_k3_VY;
    double* p_k3_VZ;
    double* p_k3_VRoll;
    double* p_k3_VPitch;
    double* p_k3_VYaw;
    double* p_k4_Z;
    double* p_k4_Roll;
    double* p_k4_Pitch;
    double* p_k4_Yaw;
    double* p_k4_VX;
    double* p_k4_VY;
    double* p_k4_VZ;
    double* p_k4_VRoll;
    double* p_k4_VPitch;
    double* p_k4_VYaw;

    // Outputs
    mxGPUArray* Zn;
    mxGPUArray* Rolln;
    mxGPUArray* Pitchn;
    mxGPUArray* Yawn;
    mxGPUArray* VXn;
    mxGPUArray* VYn;
    mxGPUArray* VZn;
    mxGPUArray* VRolln;
    mxGPUArray* VPitchn;
    mxGPUArray* VYawn;
    // Pointers for outputs
    double* p_Zn;
    double* p_Rolln;
    double* p_Pitchn;
    double* p_Yawn;
    double* p_VXn;
    double* p_VYn;
    double* p_VZn;
    double* p_VRolln;
    double* p_VPitchn;
    double* p_VYawn;

    char const* const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const* const errMsg = "Invalid input to MEX file.";

    /* Choose a reasonably sized number of threads for the block. */
    int const threadsPerBlock = 256;
    int const blocksPerGrid = 1024;

    /* Initialize the MathWorks GPU API. */
    mxInitGPU();

    /* Throw an error if the input is not a GPU array. */
    if ((nrhs != 26) || !(mxIsGPUArray(prhs[0])) || !(mxIsGPUArray(prhs[1])) || !(mxIsGPUArray(prhs[2])) || !(mxIsGPUArray(prhs[3])) 
                     || !(mxIsGPUArray(prhs[4])) || !(mxIsGPUArray(prhs[5])) || !(mxIsGPUArray(prhs[6])) || !(mxIsGPUArray(prhs[7]))
                     || !(mxIsGPUArray(prhs[8])) || !(mxIsGPUArray(prhs[9])) || !(mxIsGPUArray(prhs[10])) || !(mxIsGPUArray(prhs[11]))
                     || !(mxIsGPUArray(prhs[12])) || !(mxIsGPUArray(prhs[13]))) {
        mexErrMsgIdAndTxt(errId, errMsg);
    }

    Z      = mxGPUCreateFromMxArray(prhs[0]);
    Roll   = mxGPUCreateFromMxArray(prhs[1]);
    Pitch  = mxGPUCreateFromMxArray(prhs[2]);
    Yaw    = mxGPUCreateFromMxArray(prhs[3]);
    VX     = mxGPUCreateFromMxArray(prhs[4]);
    VY     = mxGPUCreateFromMxArray(prhs[5]);
    VZ     = mxGPUCreateFromMxArray(prhs[6]);
    VRoll  = mxGPUCreateFromMxArray(prhs[7]);
    VPitch = mxGPUCreateFromMxArray(prhs[8]);
    VYaw   = mxGPUCreateFromMxArray(prhs[9]);
    in1 = mxGPUCreateFromMxArray(prhs[10]);
    in2 = mxGPUCreateFromMxArray(prhs[11]);
    in3 = mxGPUCreateFromMxArray(prhs[12]);
    in4 = mxGPUCreateFromMxArray(prhs[13]);
    grid_size = mxGPUCreateFromMxArray(prhs[22]);
    active_actions = mxGPUCreateFromMxArray(prhs[23]);
    
    p_m = mxGetDoubles(prhs[14]); 
    double const m = p_m[0];
    
    p_I = mxGetDoubles(prhs[15]); 
    double const I1 = p_I[0];
    p_I = mxGetDoubles(prhs[16]); 
    double const I2 = p_I[0];
    p_I = mxGetDoubles(prhs[17]); 
    double const I3 = p_I[0];

    p_l = mxGetDoubles(prhs[18]); 
    double const l = p_l[0];

    p_bk = mxGetDoubles(prhs[19]); 
    double const bk = p_bk[0];

    p_g = mxGetDoubles(prhs[20]);
    double const g = p_g[0];
    
    p_dt = mxGetDoubles(prhs[21]);
    double const dt = p_dt[0];
    
    p_limits = mxGetDoubles(prhs[24]);
    
    p_x_dims_free = mxGetDoubles(prhs[25]);
    mwSize const* num_x_dims = mxGetDimensions(prhs[25]);
    if ((mxGetNumberOfDimensions(prhs[25]) != 2) || (num_x_dims[1] > 1)) {
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    
    /*
     * Verify that inputs are of appropriate type before extracting the pointer.
     */
    if ((mxGPUGetClassID(Z) != mxDOUBLE_CLASS) || (mxGPUGetClassID(Roll) != mxDOUBLE_CLASS)
        || (mxGPUGetClassID(Pitch) != mxDOUBLE_CLASS) || (mxGPUGetClassID(Yaw) != mxDOUBLE_CLASS)
        || (mxGPUGetClassID(VX) != mxDOUBLE_CLASS) || (mxGPUGetClassID(VY) != mxDOUBLE_CLASS)
        || (mxGPUGetClassID(VZ) != mxDOUBLE_CLASS) || (mxGPUGetClassID(VRoll) != mxDOUBLE_CLASS) 
        || (mxGPUGetClassID(VPitch) != mxDOUBLE_CLASS) || (mxGPUGetClassID(VYaw) != mxDOUBLE_CLASS)        
        || (mxGPUGetClassID(in1) != mxDOUBLE_CLASS) || (mxGPUGetClassID(in2) != mxDOUBLE_CLASS)
        || (mxGPUGetClassID(in3) != mxDOUBLE_CLASS) || (mxGPUGetClassID(in4) != mxDOUBLE_CLASS)
        || (mxGPUGetClassID(grid_size) != mxINT32_CLASS) || (mxGPUGetClassID(active_actions) != mxINT32_CLASS)) {
        mexErrMsgIdAndTxt(errId, errMsg);
    }

    /*
     * Now that we have verified the data type, extract a pointer to the input
     * data on the device.
     */
    p_Z      = (double*)(mxGPUGetDataReadOnly(Z));
    p_Roll   = (double*)(mxGPUGetDataReadOnly(Roll));
    p_Pitch  = (double*)(mxGPUGetDataReadOnly(Pitch));
    p_Yaw    = (double*)(mxGPUGetDataReadOnly(Yaw));
    p_VX     = (double*)(mxGPUGetDataReadOnly(VX));
    p_VY     = (double*)(mxGPUGetDataReadOnly(VY));
    p_VZ     = (double*)(mxGPUGetDataReadOnly(VZ));
    p_VRoll  = (double*)(mxGPUGetDataReadOnly(VRoll));
    p_VPitch = (double*)(mxGPUGetDataReadOnly(VPitch));
    p_VYaw   = (double*)(mxGPUGetDataReadOnly(VYaw));
    p_in1 = (double const*)(mxGPUGetDataReadOnly(in1));
    p_in2 = (double const*)(mxGPUGetDataReadOnly(in2));
    p_in3 = (double const*)(mxGPUGetDataReadOnly(in3));
    p_in4 = (double const*)(mxGPUGetDataReadOnly(in4));
    p_grid_size = (int32_t const*)(mxGPUGetDataReadOnly(grid_size));
    p_active_actions = (int32_t const*)(mxGPUGetDataReadOnly(active_actions));
    
    /* Create output arrays*/
    Zn = mxGPUCopyGPUArray(Z);
    p_Zn = (double*)(mxGPUGetData(Zn));
    
    Rolln = mxGPUCopyGPUArray(Roll);
    p_Rolln = (double*)(mxGPUGetData(Rolln));
    
    Pitchn = mxGPUCopyGPUArray(Pitch);
    p_Pitchn = (double*)(mxGPUGetData(Pitchn));
    
    Yawn = mxGPUCopyGPUArray(Yaw);
    p_Yawn = (double*)(mxGPUGetData(Yawn));
    
    VXn = mxGPUCopyGPUArray(VX);
    p_VXn = (double*)(mxGPUGetData(VXn));
    
    VYn = mxGPUCopyGPUArray(VY);
    p_VYn = (double*)(mxGPUGetData(VYn));
    
    VZn = mxGPUCopyGPUArray(VZ);
    p_VZn = (double*)(mxGPUGetData(VZn));
    
    VRolln = mxGPUCopyGPUArray(VRoll);
    p_VRolln = (double*)(mxGPUGetData(VRolln));

    VPitchn = mxGPUCopyGPUArray(VPitch);
    p_VPitchn = (double*)(mxGPUGetData(VPitchn));

    VYawn = mxGPUCopyGPUArray(VYaw);
    p_VYawn = (double*)(mxGPUGetData(VYawn));
    
    // RK4 - Step 1
    int32_t curr_free_dim;
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            p_k1_Z = p_VZ;
            step << <blocksPerGrid, threadsPerBlock >> > (p_Zn, p_Z, p_k1_Z, 0.5 * dt, p_grid_size);
    
        } else if (curr_free_dim == 2) {
            p_k1_Roll = p_VRoll;
            step << <blocksPerGrid, threadsPerBlock >> > (p_Rolln, p_Roll, p_k1_Roll, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 3) {
            p_k1_Pitch = p_VPitch;
            step << <blocksPerGrid, threadsPerBlock >> > (p_Pitchn, p_Pitch, p_k1_Pitch, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 4) {
            p_k1_Yaw = p_VYaw;
            step << <blocksPerGrid, threadsPerBlock >> > (p_Yawn, p_Yaw, p_k1_Yaw, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 5) {
            k1_VX = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VX), mxGPUGetDimensions(VX), mxGPUGetClassID(VX), mxGPUGetComplexity(VX), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_VX = (double*)(mxGPUGetData(k1_VX));
            dyn5_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_VX, p_Z, p_Roll, p_Pitch, p_Yaw, p_VX, p_VY, p_VZ, p_VRoll, p_VPitch, p_VYaw, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_VXn, p_VX, p_k1_VX, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 6) {
            k1_VY = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VY), mxGPUGetDimensions(VY), mxGPUGetClassID(VY), mxGPUGetComplexity(VY), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_VY = (double*)(mxGPUGetData(k1_VY));
            dyn6_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_VY, p_Z, p_Roll, p_Pitch, p_Yaw, p_VX, p_VY, p_VZ, p_VRoll, p_VPitch, p_VYaw, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_VYn, p_VY, p_k1_VY, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 7) {
            k1_VZ = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VZ), mxGPUGetDimensions(VZ), mxGPUGetClassID(VZ), mxGPUGetComplexity(VZ), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_VZ = (double*)(mxGPUGetData(k1_VZ));
            dyn7_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_VZ, p_Z, p_Roll, p_Pitch, p_Yaw, p_VX, p_VY, p_VZ, p_VRoll, p_VPitch, p_VYaw, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_VZn, p_VZ, p_k1_VZ, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 8) {
            k1_VRoll = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VRoll), mxGPUGetDimensions(VRoll), mxGPUGetClassID(VRoll), mxGPUGetComplexity(VRoll), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_VRoll = (double*)(mxGPUGetData(k1_VRoll));
            dyn8_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_VRoll, p_Z, p_Roll, p_Pitch, p_Yaw, p_VX, p_VY, p_VZ, p_VRoll, p_VPitch, p_VYaw, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_VRolln, p_VRoll, p_k1_VRoll, 0.5 * dt, p_grid_size);
        } else if (curr_free_dim == 9) {
            k1_VPitch = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VPitch), mxGPUGetDimensions(VPitch), mxGPUGetClassID(VPitch), mxGPUGetComplexity(VPitch), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_VPitch = (double*)(mxGPUGetData(k1_VPitch));
            dyn9_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_VPitch, p_Z, p_Roll, p_Pitch, p_Yaw, p_VX, p_VY, p_VZ, p_VRoll, p_VPitch, p_VYaw, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_VPitchn, p_VPitch, p_k1_VPitch, 0.5 * dt, p_grid_size);
        } else if (curr_free_dim == 10) {
            k1_VYaw = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VYaw), mxGPUGetDimensions(VYaw), mxGPUGetClassID(VYaw), mxGPUGetComplexity(VYaw), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_VYaw = (double*)(mxGPUGetData(k1_VYaw));
            dyn10_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_VYaw, p_Z, p_Roll, p_Pitch, p_Yaw, p_VX, p_VY, p_VZ, p_VRoll, p_VPitch, p_VYaw, 
                                                                          p_in1, p_in2, p_in3, p_in4,
                                                                          m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_VYawn, p_VYaw, p_k1_VYaw, 0.5 * dt, p_grid_size);
        }
    }
    
    // RK4 - Step 2
    // Continuous Dynamics
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            k2_Z = mxGPUCopyGPUArray(VZn);
            p_k2_Z = (double*)(mxGPUGetData(k2_Z));
    
        } else if (curr_free_dim == 2) {
            k2_Roll = mxGPUCopyGPUArray(VRolln);
            p_k2_Roll = (double*)(mxGPUGetData(k2_Roll));
            
        } else if (curr_free_dim == 3) {
            k2_Pitch = mxGPUCopyGPUArray(VPitchn);
            p_k2_Pitch = (double*)(mxGPUGetData(k2_Pitch));
            
        } else if (curr_free_dim == 4) {
            k2_Yaw = mxGPUCopyGPUArray(VYawn);
            p_k2_Yaw = (double*)(mxGPUGetData(k2_Yaw));
            
        } else if (curr_free_dim == 5) {
            k2_VX = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VX), mxGPUGetDimensions(VX), mxGPUGetClassID(VX), mxGPUGetComplexity(VX), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_VX = (double*)(mxGPUGetData(k2_VX));
            dyn5_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_VX, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
        } else if (curr_free_dim == 6) {
            k2_VY = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VY), mxGPUGetDimensions(VY), mxGPUGetClassID(VY), mxGPUGetComplexity(VY), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_VY = (double*)(mxGPUGetData(k2_VY));
            dyn6_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_VY, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);            
        } else if (curr_free_dim == 7) {
            k2_VZ = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VZ), mxGPUGetDimensions(VZ), mxGPUGetClassID(VZ), mxGPUGetComplexity(VZ), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_VZ = (double*)(mxGPUGetData(k2_VZ));
            dyn7_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_VZ, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions); 
        } else if (curr_free_dim == 8) {
            k2_VRoll = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VRoll), mxGPUGetDimensions(VRoll), mxGPUGetClassID(VRoll), mxGPUGetComplexity(VRoll), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_VRoll = (double*)(mxGPUGetData(k2_VRoll));
            dyn8_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_VRoll, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
        } else if (curr_free_dim == 9) {
            k2_VPitch = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VPitch), mxGPUGetDimensions(VPitch), mxGPUGetClassID(VPitch), mxGPUGetComplexity(VPitch), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_VPitch = (double*)(mxGPUGetData(k2_VPitch));
            dyn9_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_VPitch, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
        } else if (curr_free_dim == 10) {
            k2_VYaw = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VYaw), mxGPUGetDimensions(VYaw), mxGPUGetClassID(VYaw), mxGPUGetComplexity(VYaw), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_VYaw = (double*)(mxGPUGetData(k2_VYaw));
            dyn10_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_VYaw, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                          p_in1, p_in2, p_in3, p_in4,
                                                                          m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
        }

    }
    // Discrete step
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_Zn, p_Z, p_k2_Z, 0.5 * dt, p_grid_size);
    
        } else if (curr_free_dim == 2) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_Rolln, p_Roll, p_k2_Roll, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 3) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_Pitchn, p_Pitch, p_k2_Pitch, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 4) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_Yawn, p_Yaw, p_k2_Yaw, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 5) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_VXn, p_VX, p_k2_VX, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 6) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_VYn, p_VY, p_k2_VY, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 7) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_VZn, p_VZ, p_k2_VZ, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 8) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_VRolln, p_VRoll, p_k2_VRoll, 0.5 * dt, p_grid_size);

        } else if (curr_free_dim == 9) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_VPitchn, p_VPitch, p_k2_VPitch, 0.5 * dt, p_grid_size);

        } else if (curr_free_dim == 10) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_VYawn, p_VYaw, p_k2_VYaw, 0.5 * dt, p_grid_size);

        }
    }
    
    // RK4 - Step 3
    // Continuous Dynamics
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            k3_Z = mxGPUCopyGPUArray(VZn);
            p_k3_Z = (double*)(mxGPUGetData(k3_Z));
    
        } else if (curr_free_dim == 2) {
            k3_Roll = mxGPUCopyGPUArray(VRolln);
            p_k3_Roll = (double*)(mxGPUGetData(k3_Roll));
            
        } else if (curr_free_dim == 3) {
            k3_Pitch = mxGPUCopyGPUArray(VPitchn);
            p_k3_Pitch = (double*)(mxGPUGetData(k3_Pitch));
            
        } else if (curr_free_dim == 4) {
            k3_Yaw = mxGPUCopyGPUArray(VYawn);
            p_k3_Yaw = (double*)(mxGPUGetData(k3_Yaw));
            
        } else if (curr_free_dim == 5) {
            k3_VX = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VX), mxGPUGetDimensions(VX), mxGPUGetClassID(VX), mxGPUGetComplexity(VX), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_VX = (double*)(mxGPUGetData(k3_VX));
            dyn5_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_VX, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
        } else if (curr_free_dim == 6) {
            k3_VY = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VY), mxGPUGetDimensions(VY), mxGPUGetClassID(VY), mxGPUGetComplexity(VY), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_VY = (double*)(mxGPUGetData(k3_VY));
            dyn6_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_VY, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);            
        } else if (curr_free_dim == 7) {
            k3_VZ = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VZ), mxGPUGetDimensions(VZ), mxGPUGetClassID(VZ), mxGPUGetComplexity(VZ), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_VZ = (double*)(mxGPUGetData(k3_VZ));
            dyn7_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_VZ, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions); 
        } else if (curr_free_dim == 8) {
            k3_VRoll = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VRoll), mxGPUGetDimensions(VRoll), mxGPUGetClassID(VRoll), mxGPUGetComplexity(VRoll), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_VRoll = (double*)(mxGPUGetData(k3_VRoll));
            dyn8_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_VRoll, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
        } else if (curr_free_dim == 9) {
            k3_VPitch = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VPitch), mxGPUGetDimensions(VPitch), mxGPUGetClassID(VPitch), mxGPUGetComplexity(VPitch), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_VPitch = (double*)(mxGPUGetData(k3_VPitch));
            dyn9_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_VPitch, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
        } else if (curr_free_dim == 10) {
            k3_VYaw = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VYaw), mxGPUGetDimensions(VYaw), mxGPUGetClassID(VYaw), mxGPUGetComplexity(VYaw), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_VYaw = (double*)(mxGPUGetData(k3_VYaw));
            dyn10_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_VYaw, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                          p_in1, p_in2, p_in3, p_in4,
                                                                          m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
        }
    }
    // Discrete step
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_Zn, p_Z, p_k3_Z, dt, p_grid_size);
    
        } else if (curr_free_dim == 2) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_Rolln, p_Roll, p_k3_Roll, dt, p_grid_size);
            
        } else if (curr_free_dim == 3) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_Pitchn, p_Pitch, p_k3_Pitch, dt, p_grid_size);
            
        } else if (curr_free_dim == 4) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_Yawn, p_Yaw, p_k3_Yaw, dt, p_grid_size);
            
        } else if (curr_free_dim == 5) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_VXn, p_VX, p_k3_VX, dt, p_grid_size);
            
        } else if (curr_free_dim == 6) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_VYn, p_VY, p_k3_VY, dt, p_grid_size);
            
        } else if (curr_free_dim == 7) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_VZn, p_VZ, p_k3_VZ, dt, p_grid_size);
            
        } else if (curr_free_dim == 8) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_VRolln, p_VRoll, p_k3_VRoll, dt, p_grid_size);

        } else if (curr_free_dim == 9) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_VPitchn, p_VPitch, p_k3_VPitch, dt, p_grid_size);

        } else if (curr_free_dim == 10) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_VYawn, p_VYaw, p_k3_VYaw, dt, p_grid_size);

        }
    }
    
    // RK4 - Step 4
    // Continuous Dynamics
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            k4_Z = mxGPUCopyGPUArray(VZn);
            p_k4_Z = (double*)(mxGPUGetData(k4_Z));
    
        } else if (curr_free_dim == 2) {
            k4_Roll = mxGPUCopyGPUArray(VRolln);
            p_k4_Roll = (double*)(mxGPUGetData(k4_Roll));
            
        } else if (curr_free_dim == 3) {
            k4_Pitch = mxGPUCopyGPUArray(VPitchn);
            p_k4_Pitch = (double*)(mxGPUGetData(k4_Pitch));
            
        } else if (curr_free_dim == 4) {
            k4_Yaw = mxGPUCopyGPUArray(VYawn);
            p_k4_Yaw = (double*)(mxGPUGetData(k4_Yaw));
            
        } else if (curr_free_dim == 5) {
            k4_VX = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VX), mxGPUGetDimensions(VX), mxGPUGetClassID(VX), mxGPUGetComplexity(VX), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_VX = (double*)(mxGPUGetData(k4_VX));
            dyn5_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_VX, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
        } else if (curr_free_dim == 6) {
            k4_VY = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VY), mxGPUGetDimensions(VY), mxGPUGetClassID(VY), mxGPUGetComplexity(VY), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_VY = (double*)(mxGPUGetData(k4_VY));
            dyn6_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_VY, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);            
        } else if (curr_free_dim == 7) {
            k4_VZ = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VZ), mxGPUGetDimensions(VZ), mxGPUGetClassID(VZ), mxGPUGetComplexity(VZ), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_VZ = (double*)(mxGPUGetData(k4_VZ));
            dyn7_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_VZ, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions); 
        } else if (curr_free_dim == 8) {
            k4_VRoll = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VRoll), mxGPUGetDimensions(VRoll), mxGPUGetClassID(VRoll), mxGPUGetComplexity(VRoll), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_VRoll = (double*)(mxGPUGetData(k4_VRoll));
            dyn8_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_VRoll, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
        } else if (curr_free_dim == 9) {
            k4_VPitch = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VPitch), mxGPUGetDimensions(VPitch), mxGPUGetClassID(VPitch), mxGPUGetComplexity(VPitch), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_VPitch = (double*)(mxGPUGetData(k4_VPitch));
            dyn9_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_VPitch, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                         p_in1, p_in2, p_in3, p_in4,
                                                                         m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
        } else if (curr_free_dim == 10) {
            k4_VYaw = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(VYaw), mxGPUGetDimensions(VYaw), mxGPUGetClassID(VYaw), mxGPUGetComplexity(VYaw), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_VYaw = (double*)(mxGPUGetData(k4_VYaw));
            dyn10_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_VYaw, p_Zn, p_Rolln, p_Pitchn, p_Yawn, p_VXn, p_VYn, p_VZn, p_VRolln, p_VPitchn, p_VYawn, 
                                                                          p_in1, p_in2, p_in3, p_in4,
                                                                          m, I1, I2, I3, l, bk, g, p_grid_size, p_active_actions);
        }
    }

    // RK4 Ouput
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_Zn, p_Z, p_k1_Z, p_k2_Z, p_k3_Z, p_k4_Z, dt, p_limits[0], p_limits[10], p_grid_size);
            mxGPUDestroyGPUArray(k2_Z);
            mxGPUDestroyGPUArray(k3_Z);
            mxGPUDestroyGPUArray(k4_Z);
        } else if (curr_free_dim == 2) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_Rolln, p_Roll, p_k1_Roll, p_k2_Roll, p_k3_Roll, p_k4_Roll, dt, p_limits[1], p_limits[11], p_grid_size);
            mxGPUDestroyGPUArray(k2_Roll);
            mxGPUDestroyGPUArray(k3_Roll);
            mxGPUDestroyGPUArray(k4_Roll);
        } else if (curr_free_dim == 3) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_Pitchn, p_Pitch, p_k1_Pitch, p_k2_Pitch, p_k3_Pitch, p_k4_Pitch, dt, p_limits[2], p_limits[12], p_grid_size);
            mxGPUDestroyGPUArray(k2_Pitch);
            mxGPUDestroyGPUArray(k3_Pitch);
            mxGPUDestroyGPUArray(k4_Pitch);
        } else if (curr_free_dim == 4) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_Yawn, p_Yaw, p_k1_Yaw, p_k2_Yaw, p_k3_Yaw, p_k4_Yaw, dt, p_limits[3], p_limits[13], p_grid_size);
            mxGPUDestroyGPUArray(k2_Yaw);
            mxGPUDestroyGPUArray(k3_Yaw);
            mxGPUDestroyGPUArray(k4_Yaw);
        } else if (curr_free_dim == 5) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_VXn, p_VX, p_k1_VX, p_k2_VX, p_k3_VX, p_k4_VX, dt, p_limits[4], p_limits[14], p_grid_size);
            mxGPUDestroyGPUArray(k1_VX);
            mxGPUDestroyGPUArray(k2_VX);
            mxGPUDestroyGPUArray(k3_VX);
            mxGPUDestroyGPUArray(k4_VX);
        } else if (curr_free_dim == 6) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_VYn, p_VY, p_k1_VY, p_k2_VY, p_k3_VY, p_k4_VY, dt, p_limits[5], p_limits[15], p_grid_size);
            mxGPUDestroyGPUArray(k1_VY);
            mxGPUDestroyGPUArray(k2_VY);
            mxGPUDestroyGPUArray(k3_VY);
            mxGPUDestroyGPUArray(k4_VY);
        } else if (curr_free_dim == 7) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_VZn, p_VZ, p_k1_VZ, p_k2_VZ, p_k3_VZ, p_k4_VZ, dt, p_limits[6], p_limits[16], p_grid_size);
            mxGPUDestroyGPUArray(k1_VZ);
            mxGPUDestroyGPUArray(k2_VZ);
            mxGPUDestroyGPUArray(k3_VZ);
            mxGPUDestroyGPUArray(k4_VZ);
        } else if (curr_free_dim == 8) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_VRolln, p_VRoll, p_k1_VRoll, p_k2_VRoll, p_k3_VRoll, p_k4_VRoll, dt, p_limits[7], p_limits[17], p_grid_size);
            mxGPUDestroyGPUArray(k1_VRoll);
            mxGPUDestroyGPUArray(k2_VRoll);
            mxGPUDestroyGPUArray(k3_VRoll);
            mxGPUDestroyGPUArray(k4_VRoll);
        } else if (curr_free_dim == 9) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_VPitchn, p_VPitch, p_k1_VPitch, p_k2_VPitch, p_k3_VPitch, p_k4_VPitch, dt, p_limits[8], p_limits[18], p_grid_size);
            mxGPUDestroyGPUArray(k1_VPitch);
            mxGPUDestroyGPUArray(k2_VPitch);
            mxGPUDestroyGPUArray(k3_VPitch);
            mxGPUDestroyGPUArray(k4_VPitch);
        } else if (curr_free_dim == 10) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_VYawn, p_VYaw, p_k1_VYaw, p_k2_VYaw, p_k3_VYaw, p_k4_VYaw, dt, p_limits[9], p_limits[19], p_grid_size);
            mxGPUDestroyGPUArray(k1_VYaw);
            mxGPUDestroyGPUArray(k2_VYaw);
            mxGPUDestroyGPUArray(k3_VYaw);
            mxGPUDestroyGPUArray(k4_VYaw);
        } 
    }
 
    /* Wrap the result up as a MATLAB gpuArray for return. */
    plhs[0] = mxGPUCreateMxArrayOnGPU(Zn);
    plhs[1] = mxGPUCreateMxArrayOnGPU(Rolln);
    plhs[2] = mxGPUCreateMxArrayOnGPU(Pitchn);
    plhs[3] = mxGPUCreateMxArrayOnGPU(Yawn);
    plhs[4] = mxGPUCreateMxArrayOnGPU(VXn);
    plhs[5] = mxGPUCreateMxArrayOnGPU(VYn);
    plhs[6] = mxGPUCreateMxArrayOnGPU(VZn);
    plhs[7] = mxGPUCreateMxArrayOnGPU(VRolln);
    plhs[8] = mxGPUCreateMxArrayOnGPU(VPitchn);
    plhs[9] = mxGPUCreateMxArrayOnGPU(VYawn);

    /*
     * The mxGPUArray pointers are host-side structures that refer to device
     * data. These must be destroyed before leaving the MEX function.
     */
    mxGPUDestroyGPUArray(Z);
    mxGPUDestroyGPUArray(Roll);
    mxGPUDestroyGPUArray(Pitch);
    mxGPUDestroyGPUArray(Yaw);
    mxGPUDestroyGPUArray(VX);
    mxGPUDestroyGPUArray(VY);
    mxGPUDestroyGPUArray(VZ);
    mxGPUDestroyGPUArray(VRoll);
    mxGPUDestroyGPUArray(VPitch);
    mxGPUDestroyGPUArray(VYaw);
    mxGPUDestroyGPUArray(in1);
    mxGPUDestroyGPUArray(in2);
    mxGPUDestroyGPUArray(in3);
    mxGPUDestroyGPUArray(in4);
    mxGPUDestroyGPUArray(grid_size);
    mxGPUDestroyGPUArray(active_actions);
    mxGPUDestroyGPUArray(Zn);
    mxGPUDestroyGPUArray(Rolln);
    mxGPUDestroyGPUArray(Pitchn);
    mxGPUDestroyGPUArray(Yawn);
    mxGPUDestroyGPUArray(VXn);
    mxGPUDestroyGPUArray(VYn);
    mxGPUDestroyGPUArray(VZn);
    mxGPUDestroyGPUArray(VRolln);
    mxGPUDestroyGPUArray(VPitchn);
    mxGPUDestroyGPUArray(VYawn);
}

