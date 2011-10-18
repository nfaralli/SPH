#ifndef __MATHE_H__
#define __MATHE_H__

#include <stdlib.h>
#include <math.h>

#define PI (3.141592653589793)

float  myrand();
//float* invert3x3(float *a, float *dest, float *determinant);
//float* mult3x3_3x1(float *a, float *b, float *dest);

float  kernelPoly6    (float *r, float h);
float* kernelPoly6Grad(float *r, float h);
float  kernelPoly6Lap (float *r, float h);

float  kernelSpiky    (float *r, float h);
float* kernelSpikyGrad(float *r, float h);
float  kernelSpikyLap (float *r, float h);

float  kernelVisco    (float *r, float h);
float* kernelViscoGrad(float *r, float h);
float  kernelViscoLap (float *r, float h);

#endif
