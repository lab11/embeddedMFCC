/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_mfcc_api.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Dec-2016 14:53:47
 */

#ifndef _CODER_MFCC_API_H
#define _CODER_MFCC_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_mfcc_api.h"

/* Type Definitions */
#ifndef struct_emxArray_creal_T
#define struct_emxArray_creal_T

struct emxArray_creal_T
{
  creal_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_creal_T*/

#ifndef typedef_emxArray_creal_T
#define typedef_emxArray_creal_T

typedef struct emxArray_creal_T emxArray_creal_T;

#endif                                 /*typedef_emxArray_creal_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void mfcc(real_T currFrame, emxArray_creal_T *mfcc_coeff);
extern void mfcc_api(const mxArray * const prhs[1], const mxArray *plhs[1]);
extern void mfcc_atexit(void);
extern void mfcc_initialize(void);
extern void mfcc_terminate(void);
extern void mfcc_xil_terminate(void);

#endif

/*
 * File trailer for _coder_mfcc_api.h
 *
 * [EOF]
 */
