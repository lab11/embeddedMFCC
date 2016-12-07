/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: bluestein_setup.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Dec-2016 14:53:47
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "mfcc.h"
#include "bluestein_setup.h"

/* Function Definitions */

/*
 * Arguments    : creal_T wwc[2645]
 * Return Type  : void
 */
void bluestein_setup(creal_T wwc[2645])
{
  int idx;
  int rt;
  int k;
  int y;
  double nt_im;
  double nt_re;
  idx = 1321;
  rt = 0;
  wwc[1322].re = 1.0;
  wwc[1322].im = 0.0;
  for (k = 0; k < 1322; k++) {
    y = ((k + 1) << 1) - 1;
    if (2646 - rt <= y) {
      rt = (y + rt) - 2646;
    } else {
      rt += y;
    }

    nt_im = -3.1415926535897931 * (double)rt / 1323.0;
    if (nt_im == 0.0) {
      nt_re = 1.0;
      nt_im = 0.0;
    } else {
      nt_re = cos(nt_im);
      nt_im = sin(nt_im);
    }

    wwc[idx].re = nt_re;
    wwc[idx].im = -nt_im;
    idx--;
  }

  idx = 0;
  for (k = 1321; k >= 0; k += -1) {
    wwc[k + 1323] = wwc[idx];
    idx++;
  }
}

/*
 * File trailer for bluestein_setup.c
 *
 * [EOF]
 */
