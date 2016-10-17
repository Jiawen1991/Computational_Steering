//
// Created by Yonghong Yan on 4/26/15.
//
#ifndef OFFOMP_STENCIL2D_H
#define OFFOMP_STENCIL2D_H
#ifdef __cplusplus
extern "C" {
#endif
// flexible between REAL and double
#define REAL float
#include "homp.h"

/* this should go to a compiler generated header file and to be included in the souce files (launcher, off pragma call, etc)*/
struct stencil2d_off_args {
    long n; long m; REAL *u; int radius; REAL *coeff; int num_its;

    long u_dimX;
    long u_dimY;
    int coeff_dimX;
    REAL * coeff_center;
    REAL * uold;
};

extern int dist_dim;
extern int dist_policy;
extern double stencil2d_omp_mdev_iterate(int ndevs, int *targets, long n, long m, REAL *u, int radius, REAL *coeff, int num_its);
extern double stencil2d_omp_mdev(int ndevs, int *targets, long n, long m, REAL *u, int radius, REAL *coeff, int num_its);

extern void stencil2d_cpu_omp_wrapper(omp_offloading_t *off, int start_n, int len_n, long n, long m, int u_dimX, int u_dimY, REAL *u, REAL *uold, int radius, int coeff_dimX, REAL *coeff);
extern void stencil2d_nvgpu_cuda_wrapper(omp_offloading_t *off, int start_n, int len_n, long n, long m, int u_dimX, int u_dimY, REAL *u, REAL *uold, int radius, int coeff_dimX, REAL *coeff);
extern void stencil2d_itlmic_wrapper(omp_offloading_t *off, int start, int len, long n, long m, int u_dimX, int u_dimY, REAL *u, REAL *uold, int radius, int coeff_dimX, REAL *coeff);

#ifdef __cplusplus
 }
#endif
#endif //OFFOMP_STENCIL2D_H

