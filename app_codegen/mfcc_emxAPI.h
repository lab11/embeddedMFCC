/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mfcc_emxAPI.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Dec-2016 14:53:47
 */

#ifndef MFCC_EMXAPI_H
#define MFCC_EMXAPI_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "mfcc_types.h"

/* Function Declarations */
extern emxArray_creal_T *emxCreateND_creal_T(int numDimensions, int *size);
extern emxArray_creal_T *emxCreateWrapperND_creal_T(creal_T *data, int
  numDimensions, int *size);
extern emxArray_creal_T *emxCreateWrapper_creal_T(creal_T *data, int rows, int
  cols);
extern emxArray_creal_T *emxCreate_creal_T(int rows, int cols);
extern void emxDestroyArray_creal_T(emxArray_creal_T *emxArray);
extern void emxInitArray_creal_T(emxArray_creal_T **pEmxArray, int numDimensions);

#endif

/*
 * File trailer for mfcc_emxAPI.h
 *
 * [EOF]
 */
