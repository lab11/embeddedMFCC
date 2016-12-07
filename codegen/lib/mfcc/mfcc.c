/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mfcc.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Dec-2016 14:53:47
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "mfcc.h"
#include "mfcc_emxutil.h"
#include "bluestein_setup.h"

/* Function Declarations */
static void b_r2br_r2dit_trig(const creal_T x[4096], const double costab[2049],
  const double sintab[2049], creal_T y[4096]);
static void b_r2br_r2dit_trig_impl(const emxArray_creal_T *x, int unsigned_nRows,
  const emxArray_real_T *costab, const emxArray_real_T *sintab, emxArray_creal_T
  *y);
static void dobluesteinfft(const emxArray_real_T *x, int N2, int n1, const
  emxArray_real_T *costab, const emxArray_real_T *sintab, const emxArray_real_T *
  sintabinv, emxArray_creal_T *y);
static void eml_fft(const emxArray_real_T *x, int n, emxArray_creal_T *y);
static void r2br_r2dit_trig(const creal_T x[2645], const double costab[2049],
  const double sintab[2049], creal_T y[4096]);
static void r2br_r2dit_trig_impl(const creal_T x[1323], const double costab[2049],
  const double sintab[2049], creal_T y[4096]);
static double rt_hypotd_snf(double u0, double u1);

/* Function Definitions */

/*
 * Arguments    : const creal_T x[4096]
 *                const double costab[2049]
 *                const double sintab[2049]
 *                creal_T y[4096]
 * Return Type  : void
 */
static void b_r2br_r2dit_trig(const creal_T x[4096], const double costab[2049],
  const double sintab[2049], creal_T y[4096])
{
  int ix;
  int ju;
  int iy;
  int i;
  boolean_T tst;
  double temp_re;
  double temp_im;
  int iheight;
  int istart;
  int j;
  double twid_re;
  double twid_im;
  int ihi;
  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 0; i < 4095; i++) {
    y[iy] = x[ix];
    iy = 4096;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }

    iy = ju;
    ix++;
  }

  y[iy] = x[ix];
  for (i = 0; i <= 4095; i += 2) {
    temp_re = y[i + 1].re;
    temp_im = y[i + 1].im;
    y[i + 1].re = y[i].re - y[i + 1].re;
    y[i + 1].im = y[i].im - y[i + 1].im;
    y[i].re += temp_re;
    y[i].im += temp_im;
  }

  iy = 2;
  ix = 4;
  ju = 1024;
  iheight = 4093;
  while (ju > 0) {
    for (i = 0; i < iheight; i += ix) {
      temp_re = y[i + iy].re;
      temp_im = y[i + iy].im;
      y[i + iy].re = y[i].re - temp_re;
      y[i + iy].im = y[i].im - temp_im;
      y[i].re += temp_re;
      y[i].im += temp_im;
    }

    istart = 1;
    for (j = ju; j < 2048; j += ju) {
      twid_re = costab[j];
      twid_im = sintab[j];
      i = istart;
      ihi = istart + iheight;
      while (i < ihi) {
        temp_re = twid_re * y[i + iy].re - twid_im * y[i + iy].im;
        temp_im = twid_re * y[i + iy].im + twid_im * y[i + iy].re;
        y[i + iy].re = y[i].re - temp_re;
        y[i + iy].im = y[i].im - temp_im;
        y[i].re += temp_re;
        y[i].im += temp_im;
        i += ix;
      }

      istart++;
    }

    ju /= 2;
    iy = ix;
    ix <<= 1;
    iheight -= iy;
  }

  for (iy = 0; iy < 4096; iy++) {
    y[iy].re *= 0.000244140625;
    y[iy].im *= 0.000244140625;
  }
}

/*
 * Arguments    : const emxArray_creal_T *x
 *                int unsigned_nRows
 *                const emxArray_real_T *costab
 *                const emxArray_real_T *sintab
 *                emxArray_creal_T *y
 * Return Type  : void
 */
static void b_r2br_r2dit_trig_impl(const emxArray_creal_T *x, int unsigned_nRows,
  const emxArray_real_T *costab, const emxArray_real_T *sintab, emxArray_creal_T
  *y)
{
  int j;
  int nRowsD2;
  int nRowsD4;
  int iy;
  int iDelta;
  int ix;
  int ju;
  int i;
  boolean_T tst;
  double temp_re;
  double temp_im;
  double twid_re;
  double twid_im;
  int ihi;
  if (x->size[0] <= unsigned_nRows) {
    j = x->size[0];
  } else {
    j = unsigned_nRows;
  }

  nRowsD2 = unsigned_nRows / 2;
  nRowsD4 = nRowsD2 / 2;
  iy = y->size[0];
  y->size[0] = unsigned_nRows;
  emxEnsureCapacity((emxArray__common *)y, iy, (int)sizeof(creal_T));
  iy = x->size[0];
  if (unsigned_nRows > iy) {
    iDelta = y->size[0];
    iy = y->size[0];
    y->size[0] = iDelta;
    emxEnsureCapacity((emxArray__common *)y, iy, (int)sizeof(creal_T));
    for (iy = 0; iy < iDelta; iy++) {
      y->data[iy].re = 0.0;
      y->data[iy].im = 0.0;
    }
  }

  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 1; i < j; i++) {
    y->data[iy] = x->data[ix];
    iDelta = unsigned_nRows;
    tst = true;
    while (tst) {
      iDelta >>= 1;
      ju ^= iDelta;
      tst = ((ju & iDelta) == 0);
    }

    iy = ju;
    ix++;
  }

  y->data[iy] = x->data[ix];
  if (unsigned_nRows > 1) {
    for (i = 0; i <= unsigned_nRows - 2; i += 2) {
      temp_re = y->data[i + 1].re;
      temp_im = y->data[i + 1].im;
      y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
      y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }
  }

  iDelta = 2;
  iy = 4;
  ix = 1 + ((nRowsD4 - 1) << 2);
  while (nRowsD4 > 0) {
    for (i = 0; i < ix; i += iy) {
      temp_re = y->data[i + iDelta].re;
      temp_im = y->data[i + iDelta].im;
      y->data[i + iDelta].re = y->data[i].re - temp_re;
      y->data[i + iDelta].im = y->data[i].im - temp_im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }

    ju = 1;
    for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
      twid_re = costab->data[j];
      twid_im = sintab->data[j];
      i = ju;
      ihi = ju + ix;
      while (i < ihi) {
        temp_re = twid_re * y->data[i + iDelta].re - twid_im * y->data[i +
          iDelta].im;
        temp_im = twid_re * y->data[i + iDelta].im + twid_im * y->data[i +
          iDelta].re;
        y->data[i + iDelta].re = y->data[i].re - temp_re;
        y->data[i + iDelta].im = y->data[i].im - temp_im;
        y->data[i].re += temp_re;
        y->data[i].im += temp_im;
        i += iy;
      }

      ju++;
    }

    nRowsD4 /= 2;
    iDelta = iy;
    iy <<= 1;
    ix -= iDelta;
  }
}

/*
 * Arguments    : const emxArray_real_T *x
 *                int N2
 *                int n1
 *                const emxArray_real_T *costab
 *                const emxArray_real_T *sintab
 *                const emxArray_real_T *sintabinv
 *                emxArray_creal_T *y
 * Return Type  : void
 */
static void dobluesteinfft(const emxArray_real_T *x, int N2, int n1, const
  emxArray_real_T *costab, const emxArray_real_T *sintab, const emxArray_real_T *
  sintabinv, emxArray_creal_T *y)
{
  emxArray_creal_T *wwc;
  int nInt2m1;
  int nInt2;
  int idx;
  int rt;
  int j;
  int ihi;
  double twid_im;
  double r;
  emxArray_creal_T *fy;
  emxArray_creal_T *fv;
  int nRowsD2;
  int nRowsD4;
  int i;
  boolean_T tst;
  double temp_re;
  double temp_im;
  double fv_re;
  double fv_im;
  double wwc_im;
  double b_fv_re;
  emxInit_creal_T(&wwc, 1);
  nInt2m1 = (n1 + n1) - 1;
  nInt2 = wwc->size[0];
  wwc->size[0] = nInt2m1;
  emxEnsureCapacity((emxArray__common *)wwc, nInt2, (int)sizeof(creal_T));
  idx = n1;
  rt = 0;
  wwc->data[n1 - 1].re = 1.0;
  wwc->data[n1 - 1].im = 0.0;
  nInt2 = n1 << 1;
  for (j = 1; j < n1; j++) {
    ihi = (j << 1) - 1;
    if (nInt2 - rt <= ihi) {
      rt += ihi - nInt2;
    } else {
      rt += ihi;
    }

    twid_im = -3.1415926535897931 * (double)rt / (double)n1;
    if (twid_im == 0.0) {
      r = 1.0;
      twid_im = 0.0;
    } else {
      r = cos(twid_im);
      twid_im = sin(twid_im);
    }

    wwc->data[idx - 2].re = r;
    wwc->data[idx - 2].im = -twid_im;
    idx--;
  }

  idx = 0;
  for (j = nInt2m1 - 1; j >= n1; j--) {
    wwc->data[j] = wwc->data[idx];
    idx++;
  }

  if (n1 <= x->size[0]) {
    rt = n1;
  } else {
    rt = x->size[0];
  }

  nInt2 = y->size[0];
  y->size[0] = n1;
  emxEnsureCapacity((emxArray__common *)y, nInt2, (int)sizeof(creal_T));
  if (n1 > x->size[0]) {
    ihi = y->size[0];
    nInt2 = y->size[0];
    y->size[0] = ihi;
    emxEnsureCapacity((emxArray__common *)y, nInt2, (int)sizeof(creal_T));
    for (nInt2 = 0; nInt2 < ihi; nInt2++) {
      y->data[nInt2].re = 0.0;
      y->data[nInt2].im = 0.0;
    }
  }

  idx = 0;
  for (j = 0; j + 1 <= rt; j++) {
    r = wwc->data[(n1 + j) - 1].re;
    twid_im = wwc->data[(n1 + j) - 1].im;
    y->data[j].re = r * x->data[idx];
    y->data[j].im = twid_im * -x->data[idx];
    idx++;
  }

  while (rt + 1 <= n1) {
    y->data[rt].re = 0.0;
    y->data[rt].im = 0.0;
    rt++;
  }

  emxInit_creal_T(&fy, 1);
  b_r2br_r2dit_trig_impl(y, N2, costab, sintab, fy);
  if (wwc->size[0] <= N2) {
    nInt2m1 = wwc->size[0];
  } else {
    nInt2m1 = N2;
  }

  emxInit_creal_T(&fv, 1);
  nRowsD2 = N2 / 2;
  nRowsD4 = nRowsD2 / 2;
  nInt2 = fv->size[0];
  fv->size[0] = N2;
  emxEnsureCapacity((emxArray__common *)fv, nInt2, (int)sizeof(creal_T));
  if (N2 > wwc->size[0]) {
    idx = fv->size[0];
    nInt2 = fv->size[0];
    fv->size[0] = idx;
    emxEnsureCapacity((emxArray__common *)fv, nInt2, (int)sizeof(creal_T));
    for (nInt2 = 0; nInt2 < idx; nInt2++) {
      fv->data[nInt2].re = 0.0;
      fv->data[nInt2].im = 0.0;
    }
  }

  rt = 0;
  nInt2 = 0;
  idx = 0;
  for (i = 1; i < nInt2m1; i++) {
    fv->data[idx] = wwc->data[rt];
    idx = N2;
    tst = true;
    while (tst) {
      idx >>= 1;
      nInt2 ^= idx;
      tst = ((nInt2 & idx) == 0);
    }

    idx = nInt2;
    rt++;
  }

  fv->data[idx] = wwc->data[rt];
  if (N2 > 1) {
    for (i = 0; i <= N2 - 2; i += 2) {
      temp_re = fv->data[i + 1].re;
      temp_im = fv->data[i + 1].im;
      fv->data[i + 1].re = fv->data[i].re - fv->data[i + 1].re;
      fv->data[i + 1].im = fv->data[i].im - fv->data[i + 1].im;
      fv->data[i].re += temp_re;
      fv->data[i].im += temp_im;
    }
  }

  idx = 2;
  rt = 4;
  nInt2 = 1 + ((nRowsD4 - 1) << 2);
  while (nRowsD4 > 0) {
    for (i = 0; i < nInt2; i += rt) {
      temp_re = fv->data[i + idx].re;
      temp_im = fv->data[i + idx].im;
      fv->data[i + idx].re = fv->data[i].re - temp_re;
      fv->data[i + idx].im = fv->data[i].im - temp_im;
      fv->data[i].re += temp_re;
      fv->data[i].im += temp_im;
    }

    nInt2m1 = 1;
    for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
      r = costab->data[j];
      twid_im = sintab->data[j];
      i = nInt2m1;
      ihi = nInt2m1 + nInt2;
      while (i < ihi) {
        temp_re = r * fv->data[i + idx].re - twid_im * fv->data[i + idx].im;
        temp_im = r * fv->data[i + idx].im + twid_im * fv->data[i + idx].re;
        fv->data[i + idx].re = fv->data[i].re - temp_re;
        fv->data[i + idx].im = fv->data[i].im - temp_im;
        fv->data[i].re += temp_re;
        fv->data[i].im += temp_im;
        i += rt;
      }

      nInt2m1++;
    }

    nRowsD4 /= 2;
    idx = rt;
    rt <<= 1;
    nInt2 -= idx;
  }

  nInt2 = fy->size[0];
  emxEnsureCapacity((emxArray__common *)fy, nInt2, (int)sizeof(creal_T));
  idx = fy->size[0];
  for (nInt2 = 0; nInt2 < idx; nInt2++) {
    r = fy->data[nInt2].re;
    twid_im = fy->data[nInt2].im;
    fv_re = fv->data[nInt2].re;
    fv_im = fv->data[nInt2].im;
    fy->data[nInt2].re = r * fv_re - twid_im * fv_im;
    fy->data[nInt2].im = r * fv_im + twid_im * fv_re;
  }

  if (fy->size[0] <= N2) {
    nInt2m1 = fy->size[0];
  } else {
    nInt2m1 = N2;
  }

  nRowsD2 = N2 / 2;
  nRowsD4 = nRowsD2 / 2;
  nInt2 = fv->size[0];
  fv->size[0] = N2;
  emxEnsureCapacity((emxArray__common *)fv, nInt2, (int)sizeof(creal_T));
  if (N2 > fy->size[0]) {
    idx = fv->size[0];
    nInt2 = fv->size[0];
    fv->size[0] = idx;
    emxEnsureCapacity((emxArray__common *)fv, nInt2, (int)sizeof(creal_T));
    for (nInt2 = 0; nInt2 < idx; nInt2++) {
      fv->data[nInt2].re = 0.0;
      fv->data[nInt2].im = 0.0;
    }
  }

  rt = 0;
  nInt2 = 0;
  idx = 0;
  for (i = 1; i < nInt2m1; i++) {
    fv->data[idx] = fy->data[rt];
    idx = N2;
    tst = true;
    while (tst) {
      idx >>= 1;
      nInt2 ^= idx;
      tst = ((nInt2 & idx) == 0);
    }

    idx = nInt2;
    rt++;
  }

  fv->data[idx] = fy->data[rt];
  emxFree_creal_T(&fy);
  if (N2 > 1) {
    for (i = 0; i <= N2 - 2; i += 2) {
      temp_re = fv->data[i + 1].re;
      temp_im = fv->data[i + 1].im;
      fv->data[i + 1].re = fv->data[i].re - fv->data[i + 1].re;
      fv->data[i + 1].im = fv->data[i].im - fv->data[i + 1].im;
      fv->data[i].re += temp_re;
      fv->data[i].im += temp_im;
    }
  }

  idx = 2;
  rt = 4;
  nInt2 = 1 + ((nRowsD4 - 1) << 2);
  while (nRowsD4 > 0) {
    for (i = 0; i < nInt2; i += rt) {
      temp_re = fv->data[i + idx].re;
      temp_im = fv->data[i + idx].im;
      fv->data[i + idx].re = fv->data[i].re - temp_re;
      fv->data[i + idx].im = fv->data[i].im - temp_im;
      fv->data[i].re += temp_re;
      fv->data[i].im += temp_im;
    }

    nInt2m1 = 1;
    for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
      r = costab->data[j];
      twid_im = sintabinv->data[j];
      i = nInt2m1;
      ihi = nInt2m1 + nInt2;
      while (i < ihi) {
        temp_re = r * fv->data[i + idx].re - twid_im * fv->data[i + idx].im;
        temp_im = r * fv->data[i + idx].im + twid_im * fv->data[i + idx].re;
        fv->data[i + idx].re = fv->data[i].re - temp_re;
        fv->data[i + idx].im = fv->data[i].im - temp_im;
        fv->data[i].re += temp_re;
        fv->data[i].im += temp_im;
        i += rt;
      }

      nInt2m1++;
    }

    nRowsD4 /= 2;
    idx = rt;
    rt <<= 1;
    nInt2 -= idx;
  }

  if (fv->size[0] > 1) {
    r = 1.0 / (double)fv->size[0];
    nInt2 = fv->size[0];
    emxEnsureCapacity((emxArray__common *)fv, nInt2, (int)sizeof(creal_T));
    idx = fv->size[0];
    for (nInt2 = 0; nInt2 < idx; nInt2++) {
      fv->data[nInt2].re *= r;
      fv->data[nInt2].im *= r;
    }
  }

  idx = 0;
  for (j = n1 - 1; j + 1 <= wwc->size[0]; j++) {
    r = wwc->data[j].re;
    fv_re = fv->data[j].re;
    twid_im = wwc->data[j].im;
    fv_im = fv->data[j].im;
    temp_re = wwc->data[j].re;
    temp_im = fv->data[j].im;
    wwc_im = wwc->data[j].im;
    b_fv_re = fv->data[j].re;
    y->data[idx].re = r * fv_re + twid_im * fv_im;
    y->data[idx].im = temp_re * temp_im - wwc_im * b_fv_re;
    idx++;
  }

  emxFree_creal_T(&fv);
  emxFree_creal_T(&wwc);
}

/*
 * Arguments    : const emxArray_real_T *x
 *                int n
 *                emxArray_creal_T *y
 * Return Type  : void
 */
static void eml_fft(const emxArray_real_T *x, int n, emxArray_creal_T *y)
{
  boolean_T useRadix2;
  int pmin;
  int nn1m1;
  int pmax;
  emxArray_real_T *costab1q;
  double e;
  boolean_T exitg1;
  int nRowsD4;
  int istart;
  int pow2p;
  emxArray_real_T *costab;
  emxArray_real_T *sintab;
  emxArray_real_T *sintabinv;
  int b_x[1];
  emxArray_real_T c_x;
  int nRowsD2;
  int i;
  double temp_re;
  double temp_im;
  double twid_im;
  int ihi;
  useRadix2 = ((n & (n - 1)) == 0);
  pmin = 1;
  if (useRadix2) {
    nn1m1 = n;
  } else {
    nn1m1 = (n + n) - 1;
    pmax = 31;
    if (nn1m1 <= 1) {
      pmax = 0;
    } else {
      pmin = 0;
      exitg1 = false;
      while ((!exitg1) && (pmax - pmin > 1)) {
        istart = (pmin + pmax) >> 1;
        pow2p = 1 << istart;
        if (pow2p == nn1m1) {
          pmax = istart;
          exitg1 = true;
        } else if (pow2p > nn1m1) {
          pmax = istart;
        } else {
          pmin = istart;
        }
      }
    }

    pmin = 1 << pmax;
    nn1m1 = pmin;
  }

  emxInit_real_T(&costab1q, 2);
  e = 6.2831853071795862 / (double)nn1m1;
  nRowsD4 = nn1m1 / 2 / 2;
  pmax = costab1q->size[0] * costab1q->size[1];
  costab1q->size[0] = 1;
  costab1q->size[1] = nRowsD4 + 1;
  emxEnsureCapacity((emxArray__common *)costab1q, pmax, (int)sizeof(double));
  costab1q->data[0] = 1.0;
  nn1m1 = nRowsD4 / 2;
  for (pmax = 1; pmax <= nn1m1; pmax++) {
    costab1q->data[pmax] = cos(e * (double)pmax);
  }

  for (pmax = nn1m1 + 1; pmax < nRowsD4; pmax++) {
    costab1q->data[pmax] = sin(e * (double)(nRowsD4 - pmax));
  }

  costab1q->data[nRowsD4] = 0.0;
  emxInit_real_T(&costab, 2);
  emxInit_real_T(&sintab, 2);
  emxInit_real_T(&sintabinv, 2);
  if (!useRadix2) {
    pow2p = costab1q->size[1] - 1;
    nn1m1 = (costab1q->size[1] - 1) << 1;
    pmax = costab->size[0] * costab->size[1];
    costab->size[0] = 1;
    costab->size[1] = nn1m1 + 1;
    emxEnsureCapacity((emxArray__common *)costab, pmax, (int)sizeof(double));
    pmax = sintab->size[0] * sintab->size[1];
    sintab->size[0] = 1;
    sintab->size[1] = nn1m1 + 1;
    emxEnsureCapacity((emxArray__common *)sintab, pmax, (int)sizeof(double));
    costab->data[0] = 1.0;
    sintab->data[0] = 0.0;
    pmax = sintabinv->size[0] * sintabinv->size[1];
    sintabinv->size[0] = 1;
    sintabinv->size[1] = nn1m1 + 1;
    emxEnsureCapacity((emxArray__common *)sintabinv, pmax, (int)sizeof(double));
    for (pmax = 1; pmax <= pow2p; pmax++) {
      sintabinv->data[pmax] = costab1q->data[pow2p - pmax];
    }

    for (pmax = costab1q->size[1]; pmax <= nn1m1; pmax++) {
      sintabinv->data[pmax] = costab1q->data[pmax - pow2p];
    }

    for (pmax = 1; pmax <= pow2p; pmax++) {
      costab->data[pmax] = costab1q->data[pmax];
      sintab->data[pmax] = -costab1q->data[pow2p - pmax];
    }

    for (pmax = costab1q->size[1]; pmax <= nn1m1; pmax++) {
      costab->data[pmax] = -costab1q->data[nn1m1 - pmax];
      sintab->data[pmax] = -costab1q->data[pmax - pow2p];
    }
  } else {
    pow2p = costab1q->size[1] - 1;
    nn1m1 = (costab1q->size[1] - 1) << 1;
    pmax = costab->size[0] * costab->size[1];
    costab->size[0] = 1;
    costab->size[1] = nn1m1 + 1;
    emxEnsureCapacity((emxArray__common *)costab, pmax, (int)sizeof(double));
    pmax = sintab->size[0] * sintab->size[1];
    sintab->size[0] = 1;
    sintab->size[1] = nn1m1 + 1;
    emxEnsureCapacity((emxArray__common *)sintab, pmax, (int)sizeof(double));
    costab->data[0] = 1.0;
    sintab->data[0] = 0.0;
    for (pmax = 1; pmax <= pow2p; pmax++) {
      costab->data[pmax] = costab1q->data[pmax];
      sintab->data[pmax] = -costab1q->data[pow2p - pmax];
    }

    for (pmax = costab1q->size[1]; pmax <= nn1m1; pmax++) {
      costab->data[pmax] = -costab1q->data[nn1m1 - pmax];
      sintab->data[pmax] = -costab1q->data[pmax - pow2p];
    }

    pmax = sintabinv->size[0] * sintabinv->size[1];
    sintabinv->size[0] = 1;
    sintabinv->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)sintabinv, pmax, (int)sizeof(double));
  }

  emxFree_real_T(&costab1q);
  if (useRadix2) {
    istart = x->size[0];
    if (istart <= n) {
      istart = x->size[0];
    } else {
      istart = n;
    }

    nRowsD2 = n / 2;
    nRowsD4 = nRowsD2 / 2;
    pmax = y->size[0];
    y->size[0] = n;
    emxEnsureCapacity((emxArray__common *)y, pmax, (int)sizeof(creal_T));
    nn1m1 = x->size[0];
    if (n > nn1m1) {
      nn1m1 = y->size[0];
      pmax = y->size[0];
      y->size[0] = nn1m1;
      emxEnsureCapacity((emxArray__common *)y, pmax, (int)sizeof(creal_T));
      for (pmax = 0; pmax < nn1m1; pmax++) {
        y->data[pmax].re = 0.0;
        y->data[pmax].im = 0.0;
      }
    }

    pmax = 0;
    pmin = 0;
    nn1m1 = 0;
    for (i = 1; i < istart; i++) {
      y->data[nn1m1].re = x->data[pmax];
      y->data[nn1m1].im = 0.0;
      pow2p = n;
      useRadix2 = true;
      while (useRadix2) {
        pow2p >>= 1;
        pmin ^= pow2p;
        useRadix2 = ((pmin & pow2p) == 0);
      }

      nn1m1 = pmin;
      pmax++;
    }

    y->data[nn1m1].re = x->data[pmax];
    y->data[nn1m1].im = 0.0;
    if (n > 1) {
      for (i = 0; i <= n - 2; i += 2) {
        temp_re = y->data[i + 1].re;
        temp_im = y->data[i + 1].im;
        y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
        y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
        y->data[i].re += temp_re;
        y->data[i].im += temp_im;
      }
    }

    nn1m1 = 2;
    pmax = 4;
    pmin = 1 + ((nRowsD4 - 1) << 2);
    while (nRowsD4 > 0) {
      for (i = 0; i < pmin; i += pmax) {
        temp_re = y->data[i + nn1m1].re;
        temp_im = y->data[i + nn1m1].im;
        y->data[i + nn1m1].re = y->data[i].re - temp_re;
        y->data[i + nn1m1].im = y->data[i].im - temp_im;
        y->data[i].re += temp_re;
        y->data[i].im += temp_im;
      }

      istart = 1;
      for (pow2p = nRowsD4; pow2p < nRowsD2; pow2p += nRowsD4) {
        e = costab->data[pow2p];
        twid_im = sintab->data[pow2p];
        i = istart;
        ihi = istart + pmin;
        while (i < ihi) {
          temp_re = e * y->data[i + nn1m1].re - twid_im * y->data[i + nn1m1].im;
          temp_im = e * y->data[i + nn1m1].im + twid_im * y->data[i + nn1m1].re;
          y->data[i + nn1m1].re = y->data[i].re - temp_re;
          y->data[i + nn1m1].im = y->data[i].im - temp_im;
          y->data[i].re += temp_re;
          y->data[i].im += temp_im;
          i += pmax;
        }

        istart++;
      }

      nRowsD4 /= 2;
      nn1m1 = pmax;
      pmax <<= 1;
      pmin -= nn1m1;
    }
  } else {
    b_x[0] = x->size[0];
    c_x = *x;
    c_x.size = (int *)&b_x;
    c_x.numDimensions = 1;
    dobluesteinfft(&c_x, pmin, n, costab, sintab, sintabinv, y);
  }

  emxFree_real_T(&sintabinv);
  emxFree_real_T(&sintab);
  emxFree_real_T(&costab);
}

/*
 * Arguments    : const creal_T x[2645]
 *                const double costab[2049]
 *                const double sintab[2049]
 *                creal_T y[4096]
 * Return Type  : void
 */
static void r2br_r2dit_trig(const creal_T x[2645], const double costab[2049],
  const double sintab[2049], creal_T y[4096])
{
  int i;
  int ix;
  int ju;
  int iy;
  boolean_T tst;
  double temp_re;
  double temp_im;
  int iheight;
  int istart;
  int j;
  double twid_re;
  double twid_im;
  int ihi;
  for (i = 0; i < 4096; i++) {
    y[i].re = 0.0;
    y[i].im = 0.0;
  }

  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 0; i < 2644; i++) {
    y[iy] = x[ix];
    iy = 4096;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }

    iy = ju;
    ix++;
  }

  y[iy] = x[ix];
  for (i = 0; i <= 4095; i += 2) {
    temp_re = y[i + 1].re;
    temp_im = y[i + 1].im;
    y[i + 1].re = y[i].re - y[i + 1].re;
    y[i + 1].im = y[i].im - y[i + 1].im;
    y[i].re += temp_re;
    y[i].im += temp_im;
  }

  iy = 2;
  ix = 4;
  ju = 1024;
  iheight = 4093;
  while (ju > 0) {
    for (i = 0; i < iheight; i += ix) {
      temp_re = y[i + iy].re;
      temp_im = y[i + iy].im;
      y[i + iy].re = y[i].re - temp_re;
      y[i + iy].im = y[i].im - temp_im;
      y[i].re += temp_re;
      y[i].im += temp_im;
    }

    istart = 1;
    for (j = ju; j < 2048; j += ju) {
      twid_re = costab[j];
      twid_im = sintab[j];
      i = istart;
      ihi = istart + iheight;
      while (i < ihi) {
        temp_re = twid_re * y[i + iy].re - twid_im * y[i + iy].im;
        temp_im = twid_re * y[i + iy].im + twid_im * y[i + iy].re;
        y[i + iy].re = y[i].re - temp_re;
        y[i + iy].im = y[i].im - temp_im;
        y[i].re += temp_re;
        y[i].im += temp_im;
        i += ix;
      }

      istart++;
    }

    ju /= 2;
    iy = ix;
    ix <<= 1;
    iheight -= iy;
  }
}

/*
 * Arguments    : const creal_T x[1323]
 *                const double costab[2049]
 *                const double sintab[2049]
 *                creal_T y[4096]
 * Return Type  : void
 */
static void r2br_r2dit_trig_impl(const creal_T x[1323], const double costab[2049],
  const double sintab[2049], creal_T y[4096])
{
  int i;
  int ix;
  int ju;
  int iy;
  boolean_T tst;
  double temp_re;
  double temp_im;
  int iheight;
  int istart;
  int j;
  double twid_re;
  double twid_im;
  int ihi;
  for (i = 0; i < 4096; i++) {
    y[i].re = 0.0;
    y[i].im = 0.0;
  }

  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 0; i < 1322; i++) {
    y[iy] = x[ix];
    iy = 4096;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }

    iy = ju;
    ix++;
  }

  y[iy] = x[ix];
  for (i = 0; i <= 4095; i += 2) {
    temp_re = y[i + 1].re;
    temp_im = y[i + 1].im;
    y[i + 1].re = y[i].re - y[i + 1].re;
    y[i + 1].im = y[i].im - y[i + 1].im;
    y[i].re += temp_re;
    y[i].im += temp_im;
  }

  iy = 2;
  ix = 4;
  ju = 1024;
  iheight = 4093;
  while (ju > 0) {
    for (i = 0; i < iheight; i += ix) {
      temp_re = y[i + iy].re;
      temp_im = y[i + iy].im;
      y[i + iy].re = y[i].re - temp_re;
      y[i + iy].im = y[i].im - temp_im;
      y[i].re += temp_re;
      y[i].im += temp_im;
    }

    istart = 1;
    for (j = ju; j < 2048; j += ju) {
      twid_re = costab[j];
      twid_im = sintab[j];
      i = istart;
      ihi = istart + iheight;
      while (i < ihi) {
        temp_re = twid_re * y[i + iy].re - twid_im * y[i + iy].im;
        temp_im = twid_re * y[i + iy].im + twid_im * y[i + iy].re;
        y[i + iy].re = y[i].re - temp_re;
        y[i + iy].im = y[i].im - temp_im;
        y[i].re += temp_re;
        y[i].im += temp_im;
        i += ix;
      }

      istart++;
    }

    ju /= 2;
    iy = ix;
    ix <<= 1;
    iheight -= iy;
  }
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = b;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

/*
 * Arguments    : double currFrame
 *                emxArray_creal_T *mfcc_coeff
 * Return Type  : void
 */
void mfcc(double currFrame, emxArray_creal_T *mfcc_coeff)
{
  double mag_dft_currFrame[1323];
  int i0;
  creal_T wwc[2645];
  static const double dv0[1323] = { 0.080000000000000016, 0.0800051954580786,
    0.080020781714954414, 0.080046758418550357, 0.080083124982079867,
    0.080129880584060731, 0.080187024168332954, 0.080254554444083082,
    0.080332469885873226, 0.080420768733675541, 0.080519448992911913,
    0.080628508434499035, 0.080747944594898924, 0.08087775477617426,
    0.081017936046049621, 0.0811684852379776, 0.081329398951210274,
    0.08150067355087609, 0.081682305168061964, 0.081874289699900793,
    0.082076622809663935, 0.082289299926859238, 0.0825123162473343,
    0.082745666733385015, 0.082989346113869433, 0.083243348884326551,
    0.083507669307100962, 0.083782301411472415, 0.084067238993790538,
    0.084362475617614952, 0.084668004613860925, 0.084983819080949818,
    0.085309911884964951, 0.085646275659812976, 0.085992902807389859,
    0.086349785497753018, 0.0867169156692979, 0.087094285028940066,
    0.087481885052302588, 0.087879706983908623, 0.088287741837379141,
    0.0887059803956361, 0.089134413211110319, 0.089573030605955117,
    0.090021822672264884, 0.090480779272298917, 0.090949890038710457,
    0.091429144374780558, 0.091918531454658, 0.092418040223603337,
    0.092927659398238882, 0.093447377466803427, 0.093977182689412686,
    0.09451706309832375, 0.095067006498206041, 0.095627000466416756,
    0.096197032353281031, 0.09677708928237827, 0.097367158150832356,
    0.097967225629608357, 0.098577278163813065, 0.099197301973001417,
    0.099827283051487692, 0.10046720716866198, 0.10111705986931158,
    0.10177682647394742, 0.10244649207913603, 0.10312604155783567,
    0.10381545955973837, 0.10451473051161669, 0.10522383861767531,
    0.10594276785990792, 0.10667150199845921, 0.1074100245719915,
    0.1081583188980566, 0.10891636807347271, 0.10968415497470646,
    0.11046166225825926, 0.11124887236105935, 0.11204576750085854,
    0.11285232967663378, 0.11366854066899396, 0.11449438204059126,
    0.11532983513653777, 0.11617488108482676, 0.11702950079675906,
    0.11789367496737441, 0.11876738407588716, 0.11965060838612762,
    0.12054332794698758, 0.12144552259287117, 0.1223571719441503,
    0.12327825540762499, 0.12420875217698862, 0.125148641233298,
    0.12609790134544785, 0.12705651107065075, 0.12802444875492125,
    0.12900169253356525, 0.12998822033167362, 0.13098400986462111,
    0.13198903863856953, 0.133003283950976, 0.13402672289110579,
    0.13505933234054968, 0.13610108897374629, 0.13715196925850914,
    0.1382119494565579, 0.13928100562405477, 0.14035911361214548,
    0.14144624906750441, 0.14254238743288539, 0.14364750394767539,
    0.14476157364845493, 0.14588457136956123, 0.14701647174365695,
    0.14815724920230328, 0.14930687797653724, 0.15046533209745405,
    0.15163258539679375, 0.15280861150753194, 0.15399338386447575,
    0.15518687570486389, 0.15638906006897096, 0.15759990980071675,
    0.1588193975482794, 0.16004749576471339, 0.16128417670857176,
    0.16252941244453278, 0.16378317484403088, 0.1650454355858923,
    0.16631616615697453, 0.16759533785281056, 0.16888292177825726,
    0.17017888884814814, 0.17148320978795006, 0.17279585513442514,
    0.17411679523629553, 0.17544600025491375, 0.17678344016493663,
    0.17812908475500316, 0.17948290362841768, 0.18084486620383572,
    0.18221494171595531, 0.18359309921621159, 0.1849793075734763,
    0.18637353547476082, 0.18777575142592312, 0.18918592375237991,
    0.19060402059982162, 0.19203000993493202, 0.19346385954611195,
    0.19490553704420677, 0.19635500986323834, 0.19781224526114022,
    0.19927721032049756, 0.20074987194929045, 0.20223019688164179,
    0.20371815167856849, 0.20521370272873662, 0.20671681624922111,
    0.20822745828626849, 0.20974559471606419, 0.21127119124550298,
    0.21280421341296402, 0.21434462658908893, 0.21589239597756438,
    0.2174474866159079, 0.21900986337625755, 0.2205794909661658,
    0.22215633392939638, 0.2237403566467252, 0.22533152333674528,
    0.22692979805667479, 0.22853514470316888, 0.23014752701313529,
    0.23176690856455368, 0.23339325277729805, 0.23502652291396336,
    0.23666668208069513, 0.23831369322802287, 0.23996751915169723,
    0.24162812249353016, 0.2432954657422387, 0.24496951123429284,
    0.24665022115476554, 0.24833755753818759, 0.25003148226940491,
    0.25173195708443952, 0.25343894357135405, 0.25515240317111904,
    0.25687229717848453, 0.25859858674285391, 0.26033123286916177,
    0.2620701964187544, 0.26381543811027441, 0.2655669185205477,
    0.26732459808547404, 0.26908843710092079, 0.27085839572361997,
    0.272634433972068, 0.2744165117274287, 0.27620458873444037,
    0.277998624602324, 0.27979857880569631, 0.28160441068548508,
    0.28341607944984776, 0.28523354417509239, 0.28705676380660261,
    0.28888569715976448, 0.2907203029208974, 0.292560539648187,
    0.29440636577262108, 0.29625773959892909, 0.2981146193065235,
    0.29997696295044496, 0.30184472846230936, 0.30371787365125846,
    0.30559635620491277, 0.3074801336903269, 0.30936916355494914,
    0.31126340312758172, 0.31316280961934517, 0.31506734012464455,
    0.31697695162213912, 0.31889160097571367, 0.32081124493545343,
    0.32273584013862011, 0.32466534311063266, 0.32659971026604839,
    0.32853889790954788, 0.330482862236922, 0.33243155933606139,
    0.33438494518794815, 0.33634297566765048, 0.33830560654531955,
    0.34027279348718809, 0.34224449205657226, 0.3442206577148752,
    0.34620124582259332, 0.34818621164032454, 0.350175510329779,
    0.35216909695479143, 0.35416692648233705, 0.35616895378354829,
    0.35817513363473408, 0.36018542071840187, 0.36219976962428091,
    0.364218134850348, 0.3662404708038558, 0.36826673180236225,
    0.3702968720747627, 0.37233084576232389, 0.3743686069197194,
    0.37641010951606813, 0.37845530743597378, 0.38050415448056663,
    0.38255660436854677, 0.38461261073723019, 0.38667212714359561,
    0.38873510706533365, 0.39080150390189788, 0.39287127097555707,
    0.39494436153245016, 0.39702072874364214, 0.39910032570618148,
    0.40118310544416036, 0.40326902090977523, 0.40535802498438961,
    0.4074500704795988, 0.40954511013829559, 0.41164309663573784,
    0.41374398258061695, 0.4158477205161295, 0.4179542629210482,
    0.420063562210796, 0.4221755707385203, 0.42429024079617006,
    0.426407524615573, 0.42852737436951471, 0.43064974217281882,
    0.43277458008342917, 0.43490184010349237, 0.43703147418044191,
    0.43916343420808407, 0.44129767202768433, 0.4434341394290553,
    0.44557278815164536, 0.44771356988562955, 0.44985643627300015,
    0.45200133890865951, 0.454148229341513, 0.456297059075564, 0.458447779571009,
    0.4606003422453343, 0.4627546984744132, 0.46491079959360437,
    0.46706859689885161, 0.4692280416477832, 0.47138908506081373,
    0.47355167832224559, 0.47571577258137177, 0.47788131895357894,
    0.4800482685214526, 0.48221657233588117, 0.48438618141716228,
    0.48655704675610867, 0.48872911931515584, 0.4909023500294693,
    0.49307668980805308, 0.49525208953485872, 0.49742850006989442,
    0.49960587225033559, 0.50178415689163458, 0.50396330478863283,
    0.50614326671667131, 0.50832399343270318, 0.51050543567640549,
    0.51268754417129248, 0.5148702696258286, 0.51705356273454162,
    0.51923737417913662, 0.5214216546296101, 0.5236063547453641,
    0.52579142517632116, 0.52797681656403839, 0.53016247954282314,
    0.53234836474084757, 0.53453442278126428, 0.53672060428332147,
    0.53890685986347853, 0.54109314013652154, 0.5432793957166786,
    0.54546557721873579, 0.5476516352591525, 0.54983752045717682,
    0.55202318343596146, 0.5542085748236788, 0.55639364525463586,
    0.55857834537039, 0.56076262582086356, 0.56294643726545845,
    0.56512973037417136, 0.56731245582870737, 0.56949456432359447,
    0.57167600656729678, 0.57385673328332865, 0.57603669521136713,
    0.57821584310836549, 0.58039412774966459, 0.58257149993010549,
    0.58474791046514119, 0.58692331019194688, 0.58909764997053071,
    0.59127088068484412, 0.5934429532438914, 0.5956138185828378,
    0.59778342766411874, 0.59995173147854741, 0.60211868104642108,
    0.60428422741862831, 0.60644832167775442, 0.60861091493918629,
    0.61077195835221687, 0.61293140310114846, 0.61508920040639548,
    0.61724530152558676, 0.61939965775466566, 0.62155222042899094,
    0.623702940924436, 0.62585177065848707, 0.62799866109134062,
    0.63014356372699976, 0.63228643011437047, 0.6344272118483546,
    0.63656586057094466, 0.63870232797231563, 0.640836565791916,
    0.64296852581955821, 0.64509815989650776, 0.64722541991657079,
    0.649350257827181, 0.65147262563048536, 0.653592475384427,
    0.6557097592038299, 0.65782442926147988, 0.65993643778920419,
    0.66204573707895176, 0.66415227948387057, 0.666256017419383,
    0.66835690336426223, 0.67045488986170443, 0.67254992952040116,
    0.67464197501561052, 0.676730979090225, 0.6788168945558396,
    0.68089967429381848, 0.68297927125635793, 0.6850556384675498,
    0.687128729024443, 0.68919849609810224, 0.69126489293466631,
    0.69332787285640429, 0.69538738926276977, 0.69744339563145319,
    0.69949584551943333, 0.7015446925640263, 0.703589890483932,
    0.70563139308028067, 0.70766915423767618, 0.70970312792523715,
    0.7117332681976376, 0.71375952919614427, 0.715781865149652,
    0.71780023037571916, 0.7198145792815982, 0.72182486636526588,
    0.72383104621645167, 0.72583307351766291, 0.72783090304520859,
    0.72982448967022107, 0.73181378835967548, 0.73379875417740681,
    0.73577934228512487, 0.73775550794342792, 0.73972720651281187,
    0.74169439345468047, 0.74365702433234948, 0.74561505481205193,
    0.74756844066393868, 0.74951713776307793, 0.75146110209045214,
    0.75340028973395157, 0.75533465688936741, 0.75726415986137985,
    0.75918875506454664, 0.76110839902428629, 0.76302304837786106,
    0.76493265987535541, 0.766837190380655, 0.76873659687241824,
    0.77063083644505093, 0.772519866309673, 0.7744036437950873,
    0.7762821263487415, 0.77815527153769048, 0.780023037049555,
    0.78188538069347646, 0.783742260401071, 0.78559363422737882,
    0.787439460351813, 0.78927969707910262, 0.79111430284023565,
    0.79294323619339746, 0.79476645582490768, 0.79658392055015215,
    0.79839558931451493, 0.80020142119430371, 0.802001375397676,
    0.80379541126555965, 0.8055834882725712, 0.80736556602793208,
    0.80914160427637993, 0.81091156289907917, 0.81267540191452614,
    0.81443308147945226, 0.81618456188972566, 0.81792980358124567,
    0.81966876713083825, 0.82140141325714611, 0.82312770282151537,
    0.824847596828881, 0.826561056428646, 0.82826804291556044,
    0.8299685177305951, 0.83166244246181242, 0.83334977884523453,
    0.83503048876570718, 0.83670453425776126, 0.83837187750646991,
    0.84003248084830273, 0.84168630677197709, 0.843333317919305,
    0.84497347708603665, 0.846606747222702, 0.84823309143544634,
    0.84985247298686473, 0.85146485529683114, 0.85307020194332517,
    0.85466847666325485, 0.85625964335327476, 0.85784366607060369,
    0.859420509033834, 0.86099013662374246, 0.86255251338409211,
    0.86410760402243558, 0.86565537341091114, 0.86719578658703611,
    0.86872880875449709, 0.87025440528393583, 0.87177254171373142,
    0.873283183750779, 0.87478629727126334, 0.87628184832143141,
    0.87776980311835828, 0.87925012805070946, 0.88072278967950246,
    0.88218775473885991, 0.88364499013676168, 0.8850944629557933,
    0.88653614045388807, 0.88796999006506794, 0.88939597940017845,
    0.8908140762476201, 0.892224248574077, 0.89362646452523919,
    0.89502069242652371, 0.89640690078378837, 0.89778505828404476,
    0.8991551337961643, 0.90051709637158228, 0.90187091524499685,
    0.90321655983506344, 0.90455399974508621, 0.90588320476370454,
    0.90720414486557488, 0.9085167902120499, 0.90982111115185194,
    0.9111170782217427, 0.91240466214718952, 0.91368383384302554,
    0.91495456441410772, 0.916216825155969, 0.91747058755546729,
    0.91871582329142831, 0.91995250423528674, 0.92118060245172062,
    0.92240009019928326, 0.923610939931029, 0.92481312429513618,
    0.92600661613552426, 0.92719138849246807, 0.92836741460320638,
    0.929534667902546, 0.93069312202346288, 0.93184275079769674,
    0.932983528256343, 0.93411542863043873, 0.93523842635154519,
    0.93635249605232462, 0.9374576125671148, 0.9385537509324956,
    0.9396408863878547, 0.94071899437594531, 0.94178805054344217,
    0.94284803074149082, 0.94389891102625367, 0.9449406676594504,
    0.94597327710889423, 0.946996716049024, 0.94801096136143048,
    0.949015990135379, 0.95001177966832639, 0.95099830746643477,
    0.95197555124507882, 0.95294348892934932, 0.95390209865455222,
    0.954851358766702, 0.95579124782301139, 0.956721744592375,
    0.95764282805584977, 0.95855447740712885, 0.95945667205301244,
    0.96034939161387245, 0.96123261592411291, 0.96210632503262561,
    0.962970499203241, 0.96382511891517331, 0.96467016486346224,
    0.96550561795940881, 0.96633145933100617, 0.96714767032336624,
    0.96795423249914148, 0.96875112763894067, 0.96953833774174081,
    0.97031584502529356, 0.97108363192652725, 0.97184168110194347,
    0.97258997542800851, 0.97332849800154086, 0.9740572321400921,
    0.9747761613823247, 0.97548526948838332, 0.97618454044026159,
    0.9768739584421644, 0.977553507920864, 0.97822317352605259,
    0.97888294013068844, 0.979532792831338, 0.98017271694851238,
    0.9808026980269986, 0.981422721836187, 0.98203277437039171,
    0.98263284184916766, 0.9832229107176218, 0.983802967646719,
    0.98437299953358326, 0.984932993501794, 0.98548293690167632,
    0.98602281731058739, 0.98655262253319664, 0.98707234060176119,
    0.98758195977639673, 0.98808146854534207, 0.98857085562521951,
    0.98905010996128961, 0.989519220727701, 0.98997817732773519,
    0.990426969394045, 0.99086558678888981, 0.991294019604364,
    0.99171225816262087, 0.9921202930160915, 0.99251811494769748,
    0.99290571497106006, 0.99328308433070212, 0.993650214502247,
    0.99400709719261027, 0.9943537243401871, 0.99469008811503512,
    0.99501618091905031, 0.99533199538613915, 0.99563752438238518,
    0.99593276100620953, 0.9962176985885276, 0.996492330692899,
    0.99675665111567358, 0.99701065388613064, 0.99725433326661506,
    0.99748768375266583, 0.99771070007314089, 0.99792337719033608,
    0.99812571030009933, 0.99831769483193811, 0.998499326449124,
    0.9986706010487898, 0.99883151476202248, 0.99898206395395039,
    0.99912224522382576, 0.9992520554051012, 0.999371491565501,
    0.9994805510070881, 0.99957923126632453, 0.99966753011412679,
    0.99974544555591693, 0.99981297583166717, 0.9998701194159394,
    0.99991687501792015, 0.99995324158144971, 0.99997921828504566,
    0.99999480454192147, 1.0, 0.99999480454192147, 0.99997921828504566,
    0.99995324158144971, 0.99991687501792015, 0.9998701194159394,
    0.99981297583166717, 0.99974544555591693, 0.99966753011412679,
    0.99957923126632453, 0.9994805510070881, 0.999371491565501,
    0.9992520554051012, 0.99912224522382576, 0.99898206395395039,
    0.99883151476202248, 0.9986706010487898, 0.998499326449124,
    0.99831769483193811, 0.99812571030009933, 0.99792337719033608,
    0.99771070007314089, 0.99748768375266583, 0.99725433326661506,
    0.99701065388613064, 0.99675665111567358, 0.996492330692899,
    0.9962176985885276, 0.99593276100620953, 0.99563752438238518,
    0.99533199538613915, 0.99501618091905031, 0.99469008811503512,
    0.9943537243401871, 0.99400709719261027, 0.993650214502247,
    0.99328308433070212, 0.99290571497106006, 0.99251811494769748,
    0.9921202930160915, 0.99171225816262087, 0.991294019604364,
    0.99086558678888981, 0.990426969394045, 0.98997817732773519,
    0.989519220727701, 0.98905010996128961, 0.98857085562521951,
    0.98808146854534207, 0.98758195977639673, 0.98707234060176119,
    0.98655262253319664, 0.98602281731058739, 0.98548293690167632,
    0.984932993501794, 0.98437299953358326, 0.983802967646719,
    0.9832229107176218, 0.98263284184916766, 0.98203277437039171,
    0.981422721836187, 0.9808026980269986, 0.98017271694851238,
    0.979532792831338, 0.97888294013068844, 0.97822317352605259,
    0.977553507920864, 0.9768739584421644, 0.97618454044026159,
    0.97548526948838332, 0.9747761613823247, 0.9740572321400921,
    0.97332849800154086, 0.97258997542800851, 0.97184168110194347,
    0.97108363192652725, 0.97031584502529356, 0.96953833774174081,
    0.96875112763894067, 0.96795423249914148, 0.96714767032336624,
    0.96633145933100617, 0.96550561795940881, 0.96467016486346224,
    0.96382511891517331, 0.962970499203241, 0.96210632503262561,
    0.96123261592411291, 0.96034939161387245, 0.95945667205301244,
    0.95855447740712885, 0.95764282805584977, 0.956721744592375,
    0.95579124782301139, 0.954851358766702, 0.95390209865455222,
    0.95294348892934932, 0.95197555124507882, 0.95099830746643477,
    0.95001177966832639, 0.949015990135379, 0.94801096136143048,
    0.946996716049024, 0.94597327710889423, 0.9449406676594504,
    0.94389891102625367, 0.94284803074149082, 0.94178805054344217,
    0.94071899437594531, 0.9396408863878547, 0.9385537509324956,
    0.9374576125671148, 0.93635249605232462, 0.93523842635154519,
    0.93411542863043873, 0.932983528256343, 0.93184275079769674,
    0.93069312202346288, 0.929534667902546, 0.92836741460320638,
    0.92719138849246807, 0.92600661613552426, 0.92481312429513618,
    0.923610939931029, 0.92240009019928326, 0.92118060245172062,
    0.91995250423528674, 0.91871582329142831, 0.91747058755546729,
    0.916216825155969, 0.91495456441410772, 0.91368383384302554,
    0.91240466214718952, 0.9111170782217427, 0.90982111115185194,
    0.9085167902120499, 0.90720414486557488, 0.90588320476370454,
    0.90455399974508621, 0.90321655983506344, 0.90187091524499685,
    0.90051709637158228, 0.8991551337961643, 0.89778505828404476,
    0.89640690078378837, 0.89502069242652371, 0.89362646452523919,
    0.892224248574077, 0.8908140762476201, 0.88939597940017845,
    0.88796999006506794, 0.88653614045388807, 0.8850944629557933,
    0.88364499013676168, 0.88218775473885991, 0.88072278967950246,
    0.87925012805070946, 0.87776980311835828, 0.87628184832143141,
    0.87478629727126334, 0.873283183750779, 0.87177254171373142,
    0.87025440528393583, 0.86872880875449709, 0.86719578658703611,
    0.86565537341091114, 0.86410760402243558, 0.86255251338409211,
    0.86099013662374246, 0.859420509033834, 0.85784366607060369,
    0.85625964335327476, 0.85466847666325485, 0.85307020194332517,
    0.85146485529683114, 0.84985247298686473, 0.84823309143544634,
    0.846606747222702, 0.84497347708603665, 0.843333317919305,
    0.84168630677197709, 0.84003248084830273, 0.83837187750646991,
    0.83670453425776126, 0.83503048876570718, 0.83334977884523453,
    0.83166244246181242, 0.8299685177305951, 0.82826804291556044,
    0.826561056428646, 0.824847596828881, 0.82312770282151537,
    0.82140141325714611, 0.81966876713083825, 0.81792980358124567,
    0.81618456188972566, 0.81443308147945226, 0.81267540191452614,
    0.81091156289907917, 0.80914160427637993, 0.80736556602793208,
    0.8055834882725712, 0.80379541126555965, 0.802001375397676,
    0.80020142119430371, 0.79839558931451493, 0.79658392055015215,
    0.79476645582490768, 0.79294323619339746, 0.79111430284023565,
    0.78927969707910262, 0.787439460351813, 0.78559363422737882,
    0.783742260401071, 0.78188538069347646, 0.780023037049555,
    0.77815527153769048, 0.7762821263487415, 0.7744036437950873,
    0.772519866309673, 0.77063083644505093, 0.76873659687241824,
    0.766837190380655, 0.76493265987535541, 0.76302304837786106,
    0.76110839902428629, 0.75918875506454664, 0.75726415986137985,
    0.75533465688936741, 0.75340028973395157, 0.75146110209045214,
    0.74951713776307793, 0.74756844066393868, 0.74561505481205193,
    0.74365702433234948, 0.74169439345468047, 0.73972720651281187,
    0.73775550794342792, 0.73577934228512487, 0.73379875417740681,
    0.73181378835967548, 0.72982448967022107, 0.72783090304520859,
    0.72583307351766291, 0.72383104621645167, 0.72182486636526588,
    0.7198145792815982, 0.71780023037571916, 0.715781865149652,
    0.71375952919614427, 0.7117332681976376, 0.70970312792523715,
    0.70766915423767618, 0.70563139308028067, 0.703589890483932,
    0.7015446925640263, 0.69949584551943333, 0.69744339563145319,
    0.69538738926276977, 0.69332787285640429, 0.69126489293466631,
    0.68919849609810224, 0.687128729024443, 0.6850556384675498,
    0.68297927125635793, 0.68089967429381848, 0.6788168945558396,
    0.676730979090225, 0.67464197501561052, 0.67254992952040116,
    0.67045488986170443, 0.66835690336426223, 0.666256017419383,
    0.66415227948387057, 0.66204573707895176, 0.65993643778920419,
    0.65782442926147988, 0.6557097592038299, 0.653592475384427,
    0.65147262563048536, 0.649350257827181, 0.64722541991657079,
    0.64509815989650776, 0.64296852581955821, 0.640836565791916,
    0.63870232797231563, 0.63656586057094466, 0.6344272118483546,
    0.63228643011437047, 0.63014356372699976, 0.62799866109134062,
    0.62585177065848707, 0.623702940924436, 0.62155222042899094,
    0.61939965775466566, 0.61724530152558676, 0.61508920040639548,
    0.61293140310114846, 0.61077195835221687, 0.60861091493918629,
    0.60644832167775442, 0.60428422741862831, 0.60211868104642108,
    0.59995173147854741, 0.59778342766411874, 0.5956138185828378,
    0.5934429532438914, 0.59127088068484412, 0.58909764997053071,
    0.58692331019194688, 0.58474791046514119, 0.58257149993010549,
    0.58039412774966459, 0.57821584310836549, 0.57603669521136713,
    0.57385673328332865, 0.57167600656729678, 0.56949456432359447,
    0.56731245582870737, 0.56512973037417136, 0.56294643726545845,
    0.56076262582086356, 0.55857834537039, 0.55639364525463586,
    0.5542085748236788, 0.55202318343596146, 0.54983752045717682,
    0.5476516352591525, 0.54546557721873579, 0.5432793957166786,
    0.54109314013652154, 0.53890685986347853, 0.53672060428332147,
    0.53453442278126428, 0.53234836474084757, 0.53016247954282314,
    0.52797681656403839, 0.52579142517632116, 0.5236063547453641,
    0.5214216546296101, 0.51923737417913662, 0.51705356273454162,
    0.5148702696258286, 0.51268754417129248, 0.51050543567640549,
    0.50832399343270318, 0.50614326671667131, 0.50396330478863283,
    0.50178415689163458, 0.49960587225033559, 0.49742850006989442,
    0.49525208953485872, 0.49307668980805308, 0.4909023500294693,
    0.48872911931515584, 0.48655704675610867, 0.48438618141716228,
    0.48221657233588117, 0.4800482685214526, 0.47788131895357894,
    0.47571577258137177, 0.47355167832224559, 0.47138908506081373,
    0.4692280416477832, 0.46706859689885161, 0.46491079959360437,
    0.4627546984744132, 0.4606003422453343, 0.458447779571009, 0.456297059075564,
    0.454148229341513, 0.45200133890865951, 0.44985643627300015,
    0.44771356988562955, 0.44557278815164536, 0.4434341394290553,
    0.44129767202768433, 0.43916343420808407, 0.43703147418044191,
    0.43490184010349237, 0.43277458008342917, 0.43064974217281882,
    0.42852737436951471, 0.426407524615573, 0.42429024079617006,
    0.4221755707385203, 0.420063562210796, 0.4179542629210482,
    0.4158477205161295, 0.41374398258061695, 0.41164309663573784,
    0.40954511013829559, 0.4074500704795988, 0.40535802498438961,
    0.40326902090977523, 0.40118310544416036, 0.39910032570618148,
    0.39702072874364214, 0.39494436153245016, 0.39287127097555707,
    0.39080150390189788, 0.38873510706533365, 0.38667212714359561,
    0.38461261073723019, 0.38255660436854677, 0.38050415448056663,
    0.37845530743597378, 0.37641010951606813, 0.3743686069197194,
    0.37233084576232389, 0.3702968720747627, 0.36826673180236225,
    0.3662404708038558, 0.364218134850348, 0.36219976962428091,
    0.36018542071840187, 0.35817513363473408, 0.35616895378354829,
    0.35416692648233705, 0.35216909695479143, 0.350175510329779,
    0.34818621164032454, 0.34620124582259332, 0.3442206577148752,
    0.34224449205657226, 0.34027279348718809, 0.33830560654531955,
    0.33634297566765048, 0.33438494518794815, 0.33243155933606139,
    0.330482862236922, 0.32853889790954788, 0.32659971026604839,
    0.32466534311063266, 0.32273584013862011, 0.32081124493545343,
    0.31889160097571367, 0.31697695162213912, 0.31506734012464455,
    0.31316280961934517, 0.31126340312758172, 0.30936916355494914,
    0.3074801336903269, 0.30559635620491277, 0.30371787365125846,
    0.30184472846230936, 0.29997696295044496, 0.2981146193065235,
    0.29625773959892909, 0.29440636577262108, 0.292560539648187,
    0.2907203029208974, 0.28888569715976448, 0.28705676380660261,
    0.28523354417509239, 0.28341607944984776, 0.28160441068548508,
    0.27979857880569631, 0.277998624602324, 0.27620458873444037,
    0.2744165117274287, 0.272634433972068, 0.27085839572361997,
    0.26908843710092079, 0.26732459808547404, 0.2655669185205477,
    0.26381543811027441, 0.2620701964187544, 0.26033123286916177,
    0.25859858674285391, 0.25687229717848453, 0.25515240317111904,
    0.25343894357135405, 0.25173195708443952, 0.25003148226940491,
    0.24833755753818759, 0.24665022115476554, 0.24496951123429284,
    0.2432954657422387, 0.24162812249353016, 0.23996751915169723,
    0.23831369322802287, 0.23666668208069513, 0.23502652291396336,
    0.23339325277729805, 0.23176690856455368, 0.23014752701313529,
    0.22853514470316888, 0.22692979805667479, 0.22533152333674528,
    0.2237403566467252, 0.22215633392939638, 0.2205794909661658,
    0.21900986337625755, 0.2174474866159079, 0.21589239597756438,
    0.21434462658908893, 0.21280421341296402, 0.21127119124550298,
    0.20974559471606419, 0.20822745828626849, 0.20671681624922111,
    0.20521370272873662, 0.20371815167856849, 0.20223019688164179,
    0.20074987194929045, 0.19927721032049756, 0.19781224526114022,
    0.19635500986323834, 0.19490553704420677, 0.19346385954611195,
    0.19203000993493202, 0.19060402059982162, 0.18918592375237991,
    0.18777575142592312, 0.18637353547476082, 0.1849793075734763,
    0.18359309921621159, 0.18221494171595531, 0.18084486620383572,
    0.17948290362841768, 0.17812908475500316, 0.17678344016493663,
    0.17544600025491375, 0.17411679523629553, 0.17279585513442514,
    0.17148320978795006, 0.17017888884814814, 0.16888292177825726,
    0.16759533785281056, 0.16631616615697453, 0.1650454355858923,
    0.16378317484403088, 0.16252941244453278, 0.16128417670857176,
    0.16004749576471339, 0.1588193975482794, 0.15759990980071675,
    0.15638906006897096, 0.15518687570486389, 0.15399338386447575,
    0.15280861150753194, 0.15163258539679375, 0.15046533209745405,
    0.14930687797653724, 0.14815724920230328, 0.14701647174365695,
    0.14588457136956123, 0.14476157364845493, 0.14364750394767539,
    0.14254238743288539, 0.14144624906750441, 0.14035911361214548,
    0.13928100562405477, 0.1382119494565579, 0.13715196925850914,
    0.13610108897374629, 0.13505933234054968, 0.13402672289110579,
    0.133003283950976, 0.13198903863856953, 0.13098400986462111,
    0.12998822033167362, 0.12900169253356525, 0.12802444875492125,
    0.12705651107065075, 0.12609790134544785, 0.125148641233298,
    0.12420875217698862, 0.12327825540762499, 0.1223571719441503,
    0.12144552259287117, 0.12054332794698758, 0.11965060838612762,
    0.11876738407588716, 0.11789367496737441, 0.11702950079675906,
    0.11617488108482676, 0.11532983513653777, 0.11449438204059126,
    0.11366854066899396, 0.11285232967663378, 0.11204576750085854,
    0.11124887236105935, 0.11046166225825926, 0.10968415497470646,
    0.10891636807347271, 0.1081583188980566, 0.1074100245719915,
    0.10667150199845921, 0.10594276785990792, 0.10522383861767531,
    0.10451473051161669, 0.10381545955973837, 0.10312604155783567,
    0.10244649207913603, 0.10177682647394742, 0.10111705986931158,
    0.10046720716866198, 0.099827283051487692, 0.099197301973001417,
    0.098577278163813065, 0.097967225629608357, 0.097367158150832356,
    0.09677708928237827, 0.096197032353281031, 0.095627000466416756,
    0.095067006498206041, 0.09451706309832375, 0.093977182689412686,
    0.093447377466803427, 0.092927659398238882, 0.092418040223603337,
    0.091918531454658, 0.091429144374780558, 0.090949890038710457,
    0.090480779272298917, 0.090021822672264884, 0.089573030605955117,
    0.089134413211110319, 0.0887059803956361, 0.088287741837379141,
    0.087879706983908623, 0.087481885052302588, 0.087094285028940066,
    0.0867169156692979, 0.086349785497753018, 0.085992902807389859,
    0.085646275659812976, 0.085309911884964951, 0.084983819080949818,
    0.084668004613860925, 0.084362475617614952, 0.084067238993790538,
    0.083782301411472415, 0.083507669307100962, 0.083243348884326551,
    0.082989346113869433, 0.082745666733385015, 0.0825123162473343,
    0.082289299926859238, 0.082076622809663935, 0.081874289699900793,
    0.081682305168061964, 0.08150067355087609, 0.081329398951210274,
    0.0811684852379776, 0.081017936046049621, 0.08087775477617426,
    0.080747944594898924, 0.080628508434499035, 0.080519448992911913,
    0.080420768733675541, 0.080332469885873226, 0.080254554444083082,
    0.080187024168332954, 0.080129880584060731, 0.080083124982079867,
    0.080046758418550357, 0.080020781714954414, 0.0800051954580786,
    0.080000000000000016 };

  int xidx;
  creal_T b_y1[1323];
  int center;
  static const double costab[2049] = { 1.0, 0.99999882345170188,
    0.99999529380957619, 0.9999894110819284, 0.99998117528260111,
    0.99997058643097414, 0.9999576445519639, 0.99994234967602391,
    0.9999247018391445, 0.9999047010828529, 0.99988234745421256,
    0.99985764100582386, 0.9998305817958234, 0.99980116988788426,
    0.99976940535121528, 0.99973528826056168, 0.99969881869620425,
    0.99965999674395922, 0.99961882249517864, 0.99957529604674922,
    0.99952941750109314, 0.999481186966167, 0.99943060455546173,
    0.99937767038800285, 0.99932238458834954, 0.99926474728659442,
    0.99920475861836389, 0.99914241872481691, 0.99907772775264536,
    0.99901068585407338, 0.99894129318685687, 0.99886954991428356,
    0.99879545620517241, 0.99871901223387294, 0.99864021818026527,
    0.99855907422975931, 0.99847558057329477, 0.99838973740734016,
    0.99830154493389289, 0.99821100336047819, 0.99811811290014918,
    0.99802287377148624, 0.997925286198596, 0.99782535041111164,
    0.99772306664419164, 0.99761843513851955, 0.99751145614030345,
    0.9974021299012753, 0.99729045667869021, 0.99717643673532619,
    0.997060070339483, 0.99694135776498216, 0.99682029929116567,
    0.99669689520289606, 0.99657114579055484, 0.99644305135004263,
    0.996312612182778, 0.996179828595697, 0.996044700901252, 0.99590722941741172,
    0.99576741446765982, 0.99562525638099431, 0.99548075549192694,
    0.99533391214048228, 0.99518472667219693, 0.99503319943811863,
    0.99487933079480562, 0.9947231211043257, 0.99456457073425542,
    0.9944036800576791, 0.9942404494531879, 0.99407487930487937,
    0.99390697000235606, 0.9937367219407246, 0.9935641355205953,
    0.99338921114808065, 0.9932119492347945, 0.99303235019785141,
    0.9928504144598651, 0.992666142448948, 0.99247953459871, 0.99229059134825737,
    0.9920993131421918, 0.99190570043060933, 0.99170975366909953,
    0.9915114733187439, 0.99131085984611544, 0.99110791372327689,
    0.99090263542778, 0.99069502544266463, 0.99048508425645709,
    0.99027281236316911, 0.99005821026229712, 0.98984127845882053,
    0.98962201746320089, 0.98940042779138038, 0.989176509964781,
    0.988950264510303, 0.98872169196032378, 0.98849079285269659,
    0.98825756773074946, 0.98802201714328353, 0.98778414164457218,
    0.98754394179435923, 0.98730141815785843, 0.987056571305751,
    0.98680940181418553, 0.98655991026477541, 0.98630809724459867,
    0.98605396334619544, 0.98579750916756748, 0.98553873531217606,
    0.98527764238894122, 0.98501423101223984, 0.98474850180190421,
    0.98448045538322093, 0.984210092386929, 0.98393741344921892,
    0.98366241921173025, 0.98338511032155118, 0.98310548743121629,
    0.98282355119870524, 0.98253930228744124, 0.98225274136628937,
    0.98196386910955524, 0.98167268619698311, 0.98137919331375456,
    0.98108339115048671, 0.98078528040323043, 0.98048486177346938,
    0.98018213596811743, 0.97987710369951764, 0.97956976568544052,
    0.979260122649082, 0.9789481753190622, 0.97863392442942321,
    0.97831737071962765, 0.97799851493455714, 0.97767735782450993,
    0.97735390014520007, 0.97702814265775439, 0.97670008612871184,
    0.97636973133002114, 0.976037079039039, 0.97570213003852857,
    0.975364885116657, 0.97502534506699412, 0.97468351068851067,
    0.97433938278557586, 0.97399296216795583, 0.973644249650812,
    0.97329324605469825, 0.97293995220556018, 0.97258436893473221,
    0.97222649707893627, 0.9718663374802794, 0.97150389098625178,
    0.97113915844972509, 0.97077214072895035, 0.9704028386875555,
    0.970031253194544, 0.96965738512429245, 0.96928123535654853,
    0.96890280477642887, 0.96852209427441727, 0.96813910474636244,
    0.96775383709347551, 0.96736629222232851, 0.96697647104485207,
    0.96658437447833312, 0.9661900034454125, 0.96579335887408368,
    0.9653944416976894, 0.96499325285492032, 0.96458979328981276,
    0.96418406395174583, 0.96377606579543984, 0.963365799780954,
    0.96295326687368388, 0.96253846804435916, 0.96212140426904158,
    0.96170207652912254, 0.96128048581132064, 0.96085663310767966,
    0.96043051941556579, 0.960002145737666, 0.95957151308198452,
    0.95913862246184189, 0.9587034748958716, 0.95826607140801767,
    0.95782641302753291, 0.95738450078897586, 0.95694033573220882,
    0.9564939189023951, 0.95604525134999641, 0.95559433413077111,
    0.95514116830577078, 0.95468575494133834, 0.95422809510910567,
    0.95376818988599033, 0.95330604035419386, 0.95284164760119872,
    0.95237501271976588, 0.95190613680793235, 0.95143502096900834,
    0.95096166631157508, 0.9504860739494817, 0.950008245001843,
    0.94952818059303667, 0.94904588185270056, 0.94856134991573027,
    0.94807458592227623, 0.94758559101774109, 0.94709436635277722,
    0.94660091308328353, 0.94610523237040345, 0.94560732538052128,
    0.94510719328526061, 0.94460483726148026, 0.94410025849127266,
    0.94359345816196039, 0.94308443746609349, 0.94257319760144687,
    0.94205973977101731, 0.94154406518302081, 0.94102617505088926,
    0.9405060705932683, 0.93998375303401394, 0.93945922360218992,
    0.9389324835320646, 0.93840353406310806, 0.93787237643998989,
    0.937339011912575, 0.93680344173592156, 0.93626566717027826,
    0.93572568948108037, 0.93518350993894761, 0.93463912981968078,
    0.93409255040425887, 0.93354377297883617, 0.932992798834739,
    0.93243962926846236, 0.93188426558166815, 0.93132670908118043,
    0.93076696107898371, 0.93020502289221907, 0.92964089584318121,
    0.92907458125931586, 0.92850608047321559, 0.92793539482261789,
    0.92736252565040111, 0.92678747430458175, 0.92621024213831138,
    0.92563083050987272, 0.92504924078267758, 0.9244654743252626,
    0.92387953251128674, 0.92329141671952764, 0.92270112833387863,
    0.92210866874334518, 0.9215140393420419, 0.92091724152918952,
    0.92031827670911059, 0.91971714629122736, 0.91911385169005777,
    0.91850839432521225, 0.9179007756213905, 0.91729099700837791,
    0.9166790599210427, 0.91606496579933172, 0.91544871608826783,
    0.9148303122379462, 0.91420975570353069, 0.91358704794525081,
    0.91296219042839821, 0.91233518462332275, 0.91170603200542988,
    0.91107473405517636, 0.91044129225806725, 0.90980570810465222,
    0.90916798309052238, 0.90852811871630612, 0.90788611648766626,
    0.90724197791529582, 0.90659570451491533, 0.90594729780726846,
    0.90529675931811882, 0.90464409057824624, 0.90398929312344334,
    0.90333236849451182, 0.90267331823725883, 0.90201214390249318,
    0.901348847046022, 0.900683429228647, 0.90001589201616017,
    0.89934623697934157, 0.89867446569395382, 0.89800057974073988,
    0.89732458070541832, 0.89664647017868015, 0.89596624975618522,
    0.89528392103855758, 0.8945994856313827, 0.89391294514520325,
    0.89322430119551532, 0.89253355540276458, 0.89184070939234272,
    0.89114576479458318, 0.89044872324475788, 0.88974958638307278,
    0.88904835585466457, 0.88834503330959635, 0.88763962040285393,
    0.88693211879434219, 0.88622253014888064, 0.8855108561362,
    0.88479709843093779, 0.884081258712635, 0.88336333866573158,
    0.88264333997956279, 0.881921264348355, 0.88119711347122209,
    0.88047088905216075, 0.87974259280004741, 0.87901222642863353,
    0.87827979165654158, 0.87754529020726135, 0.87680872380914565,
    0.8760700941954066, 0.87532940310411089, 0.87458665227817611,
    0.87384184346536686, 0.87309497841829009, 0.87234605889439154,
    0.87159508665595109, 0.870842063470079, 0.87008699110871146,
    0.86932987134860684, 0.8685707059713409, 0.86780949676330332,
    0.86704624551569265, 0.866280954024513, 0.86551362409056909,
    0.86474425751946238, 0.8639728561215867, 0.86319942171212416,
    0.8624239561110405, 0.8616464611430813, 0.86086693863776731,
    0.86008539042939014, 0.85930181835700847, 0.85851622426444274,
    0.85772861000027212, 0.85693897741782876, 0.85614732837519447,
    0.855353664735196, 0.85455798836540053, 0.85376030113811141,
    0.85296060493036363, 0.85215890162391983, 0.8513551931052652,
    0.85054948126560348, 0.84974176800085255, 0.84893205521163961,
    0.84812034480329723, 0.84730663868585832, 0.84649093877405213,
    0.84567324698729907, 0.84485356524970712, 0.84403189549006641,
    0.84320823964184544, 0.84238259964318585, 0.84155497743689844,
    0.84072537497045807, 0.83989379419599952, 0.83906023707031274,
    0.83822470555483808, 0.83738720161566194, 0.836547727223512,
    0.8357062843537526, 0.83486287498638, 0.83401750110601813,
    0.83317016470191319, 0.83232086776792968, 0.83146961230254524,
    0.83061640030884631, 0.829761233794523, 0.82890411477186487,
    0.8280450452577558, 0.82718402727366913, 0.82632106284566353,
    0.82545615400437755, 0.82458930278502529, 0.82372051122739143,
    0.82284978137582643, 0.82197711527924155, 0.82110251499110465,
    0.82022598256943469, 0.819347520076797, 0.81846712958029866,
    0.81758481315158371, 0.81670057286682785, 0.81581441080673378,
    0.81492632905652662, 0.81403632970594841, 0.81314441484925359,
    0.81225058658520388, 0.81135484701706373, 0.81045719825259477,
    0.80955764240405126, 0.808656181588175, 0.80775281792619036,
    0.80684755354379933, 0.80594039057117628, 0.80503133114296366,
    0.8041203773982657, 0.80320753148064494, 0.80229279553811572,
    0.80137617172314024, 0.80045766219262282, 0.799537269107905,
    0.79861499463476093, 0.79769084094339116, 0.79676481020841883,
    0.79583690460888357, 0.794907126328237, 0.79397547755433717,
    0.79304196047944364, 0.79210657730021239, 0.7911693302176902,
    0.79023022143731, 0.78928925316888565, 0.78834642762660634,
    0.78740174702903143, 0.78645521359908577, 0.78550682956405393,
    0.78455659715557524, 0.78360451860963831, 0.78265059616657573,
    0.78169483207105939, 0.78073722857209449, 0.77977778792301455,
    0.778816512381476, 0.77785340420945315, 0.77688846567323244,
    0.77592169904340769, 0.77495310659487393, 0.7739826906068229,
    0.773010453362737, 0.77203639715038452, 0.77106052426181382,
    0.7700828369933479, 0.7691033376455797, 0.76812202852336542,
    0.7671389119358204, 0.76615399019631292, 0.765167265622459,
    0.76417874053611679, 0.76318841726338127, 0.7621962981345789,
    0.76120238548426178, 0.76020668165120242, 0.759209188978388,
    0.75820990981301528, 0.75720884650648457, 0.75620600141439454,
    0.75520137689653655, 0.75419497531688917, 0.75318679904361252,
    0.7521768504490427, 0.75116513190968637, 0.75015164580621507,
    0.74913639452345937, 0.7481193804504036, 0.74710060598018013,
    0.74608007351006378, 0.745057785441466, 0.74403374417992929,
    0.74300795213512172, 0.74198041172083107, 0.74095112535495922,
    0.7399200954595162, 0.73888732446061511, 0.737852814788466,
    0.73681656887736979, 0.73577858916571359, 0.7347388780959635,
    0.73369743811466037, 0.73265427167241282, 0.73160938122389263,
    0.73056276922782759, 0.729514438146997, 0.7284643904482252,
    0.72741262860237577, 0.726359155084346, 0.72530397237306077,
    0.724247082951467, 0.72318848930652746, 0.72212819392921535,
    0.72106619931450811, 0.72000250796138165, 0.71893712237280449,
    0.71787004505573171, 0.71680127852109954, 0.71573082528381859,
    0.71465868786276909, 0.71358486878079352, 0.71250937056469243,
    0.71143219574521643, 0.71035334685706242, 0.70927282643886569,
    0.7081906370331954, 0.70710678118654757, 0.70602126144933974,
    0.70493408037590488, 0.70384524052448494, 0.7027547444572253,
    0.70166259474016845, 0.70056879394324834, 0.69947334464028377,
    0.69837624940897292, 0.69727751083088652, 0.696177131491463,
    0.69507511398000088, 0.69397146088965389, 0.69286617481742463,
    0.69175925836415775, 0.6906507141345346, 0.68954054473706683,
    0.68842875278409044, 0.687315340891759, 0.68620031168003859,
    0.68508366777270036, 0.6839654117973154, 0.68284554638524808,
    0.68172407417164971, 0.680600997795453, 0.679476319899365,
    0.67835004312986147, 0.67722217013718033, 0.67609270357531592,
    0.674961646102012, 0.673829000378756, 0.67269476907077286,
    0.67155895484701833, 0.67042156038017309, 0.669282588346636,
    0.66814204142651845, 0.66699992230363747, 0.66585623366550972,
    0.66471097820334479, 0.66356415861203977, 0.66241577759017178,
    0.66126583783999227, 0.66011434206742048, 0.65896129298203732,
    0.65780669329707864, 0.65665054572942894, 0.65549285299961535,
    0.65433361783180044, 0.65317284295377676, 0.6520105310969595,
    0.650846684996381, 0.64968130739068319, 0.64851440102211244,
    0.64734596863651206, 0.64617601298331628, 0.64500453681554393,
    0.64383154288979139, 0.64265703396622686, 0.641481012808583,
    0.64030348218415167, 0.63912444486377573, 0.637943903621844,
    0.6367618612362842, 0.63557832048855611, 0.63439328416364549,
    0.63320675505005719, 0.63201873593980906, 0.63082922962842447,
    0.629638238914927, 0.6284457666018326, 0.62725181549514408,
    0.62605638840434352, 0.62485948814238634, 0.62366111752569453,
    0.62246127937415, 0.62125997651108755, 0.6200572117632891,
    0.61885298796097632, 0.61764730793780387, 0.61644017453085365,
    0.61523159058062682, 0.61402155893103849, 0.61281008242940971,
    0.61159716392646191, 0.61038280627630948, 0.60916701233645321,
    0.60794978496777363, 0.60673112703452448, 0.60551104140432555,
    0.604289530948156, 0.60306659854034816, 0.60184224705858, 0.600616479383869,
    0.59938929840056454, 0.59816070699634238, 0.59693070806219639,
    0.59569930449243336, 0.59446649918466443, 0.5932322950397998,
    0.591996694962041, 0.59075970185887416, 0.58952131864106394,
    0.58828154822264522, 0.587040393520918, 0.58579785745643886,
    0.58455394295301533, 0.58330865293769829, 0.58206199034077544,
    0.58081395809576453, 0.57956455913940563, 0.57831379641165559,
    0.57706167285567944, 0.57580819141784534, 0.57455335504771576,
    0.5732971666980422, 0.572039629324757, 0.57078074588696726,
    0.56952051934694714, 0.56825895267013149, 0.56699604882510868,
    0.56573181078361312, 0.5644662415205195, 0.56319934401383409,
    0.56193112124468936, 0.560661576197336, 0.55939071185913614,
    0.5581185312205561, 0.5568450372751601, 0.55557023301960218,
    0.55429412145362, 0.55301670558002747, 0.55173798840470734,
    0.55045797293660481, 0.54917666218771966, 0.54789405917310019,
    0.54661016691083486, 0.54532498842204646, 0.54403852673088382,
    0.54275078486451589, 0.54146176585312344, 0.54017147272989285,
    0.53887990853100842, 0.53758707629564539, 0.53629297906596318,
    0.53499761988709715, 0.533701001807153, 0.5324031278771979,
    0.531104001151255, 0.52980362468629461, 0.52850200154222848,
    0.52719913478190128, 0.52589502747108463, 0.524589682678469,
    0.52328310347565643, 0.52197529293715439, 0.52066625414036716,
    0.51935599016558964, 0.51804450409599934, 0.51673179901764987,
    0.51541787801946293, 0.51410274419322166, 0.512786400633563,
    0.5114688504379703, 0.51015009670676681, 0.508830142543107,
    0.50750899105297087, 0.50618664534515523, 0.50486310853126759,
    0.50353838372571758, 0.50221247404571079, 0.50088538261124071,
    0.49955711254508184, 0.49822766697278181, 0.49689704902265447,
    0.49556526182577254, 0.49423230851595967, 0.49289819222978404,
    0.4915629161065499, 0.49022648328829116, 0.48888889691976317,
    0.487550160148436, 0.48621027612448642, 0.48486924800079106,
    0.48352707893291874, 0.48218377207912272, 0.48083933060033396,
    0.47949375766015295, 0.478147056424843, 0.47679923006332209,
    0.47545028174715587, 0.47410021465054997, 0.47274903195034279,
    0.47139673682599764, 0.47004333245959562, 0.46868882203582796,
    0.46733320874198842, 0.46597649576796618, 0.46461868630623782,
    0.46325978355186015, 0.46189979070246273, 0.46053871095824,
    0.45917654752194409, 0.45781330359887717, 0.45644898239688392,
    0.45508358712634384, 0.45371712100016387, 0.45234958723377089,
    0.45098098904510386, 0.44961132965460654, 0.44824061228521989,
    0.44686884016237416, 0.44549601651398174, 0.4441221445704292,
    0.44274722756457, 0.44137126873171667, 0.43999427130963326,
    0.43861623853852766, 0.43723717366104409, 0.43585707992225547,
    0.43447596056965565, 0.43309381885315196, 0.43171065802505726,
    0.43032648134008261, 0.42894129205532949, 0.42755509343028208,
    0.42616788872679962, 0.42477968120910881, 0.42339047414379605,
    0.42200027079979968, 0.42060907444840251, 0.41921688836322391,
    0.41782371582021227, 0.41642956009763715, 0.41503442447608163,
    0.4136383122384345, 0.41224122666988289, 0.41084317105790391,
    0.40944414869225759, 0.40804416286497869, 0.40664321687036903,
    0.40524131400498986, 0.40383845756765407, 0.40243465085941843,
    0.40102989718357562, 0.39962419984564679, 0.39821756215337356,
    0.39680998741671031, 0.39540147894781635, 0.3939920400610481,
    0.39258167407295147, 0.39117038430225387, 0.38975817406985641,
    0.38834504669882625, 0.38693100551438858, 0.38551605384391885,
    0.38410019501693504, 0.38268343236508978, 0.38126576922216238,
    0.37984720892405116, 0.37842775480876556, 0.37700741021641826,
    0.37558617848921722, 0.37416406297145793, 0.37274106700951576,
    0.37131719395183749, 0.3698924471489341, 0.36846682995337232,
    0.36704034571976718, 0.36561299780477385, 0.36418478956707989,
    0.36275572436739723, 0.36132580556845428, 0.35989503653498811,
    0.35846342063373654, 0.35703096123343, 0.35559766170478385,
    0.35416352542049034, 0.35272855575521073, 0.35129275608556709,
    0.34985612979013492, 0.34841868024943456, 0.34698041084592368,
    0.34554132496398909, 0.34410142598993881, 0.34266071731199438,
    0.34121920232028236, 0.33977688440682685, 0.33833376696554113,
    0.33688985339222005, 0.3354451470845316, 0.33399965144200938,
    0.33255336986604422, 0.33110630575987643, 0.32965846252858749,
    0.3282098435790925, 0.32676045232013173, 0.32531029216226293,
    0.32385936651785285, 0.32240767880106985, 0.32095523242787521,
    0.31950203081601569, 0.31804807738501495, 0.31659337555616585,
    0.31513792875252244, 0.31368174039889152, 0.31222481392182488,
    0.31076715274961147, 0.30930876031226873, 0.30784964004153487,
    0.30638979537086092, 0.30492922973540237, 0.30346794657201132,
    0.30200594931922808, 0.30054324141727345, 0.29907982630804048,
    0.2976157074350862, 0.29615088824362379, 0.29468537218051433,
    0.29321916269425863, 0.29175226323498926, 0.29028467725446233,
    0.28881640820604948, 0.28734745954472951, 0.28587783472708062,
    0.28440753721127188, 0.28293657045705539, 0.28146493792575794,
    0.27999264308027322, 0.27851968938505306, 0.2770460803060999,
    0.27557181931095814, 0.27409690986870638, 0.272621355449949,
    0.271145159526808, 0.26966832557291509, 0.26819085706340318,
    0.26671275747489837, 0.26523403028551179, 0.26375467897483135,
    0.26227470702391359, 0.26079411791527551, 0.25931291513288623,
    0.257831102162159, 0.25634868248994291, 0.25486565960451457,
    0.25338203699557016, 0.25189781815421697, 0.25041300657296522,
    0.24892760574572015, 0.24744161916777327, 0.24595505033579459,
    0.24446790274782415, 0.24298017990326387, 0.24149188530286933,
    0.2400030224487415, 0.23851359484431842, 0.2370236059943672,
    0.23553305940497549, 0.23404195858354343, 0.23255030703877524,
    0.23105810828067111, 0.22956536582051887, 0.22807208317088573,
    0.22657826384561, 0.22508391135979283, 0.22358902922979, 0.22209362097320351,
    0.22059769010887351, 0.2191012401568698, 0.21760427463848364,
    0.21610679707621952, 0.21460881099378676, 0.21311031991609136,
    0.21161132736922755, 0.21011183688046961, 0.20861185197826349,
    0.20711137619221856, 0.20561041305309924, 0.20410896609281687,
    0.20260703884442113, 0.2011046348420919, 0.19960175762113097,
    0.19809841071795356, 0.19659459767008022, 0.19509032201612825,
    0.19358558729580361, 0.19208039704989244, 0.19057475482025274,
    0.18906866414980619, 0.1875621285825296, 0.18605515166344663,
    0.18454773693861962, 0.18303988795514095, 0.18153160826112497,
    0.18002290140569951, 0.17851377093899751, 0.17700422041214875,
    0.17549425337727143, 0.17398387338746382, 0.17247308399679595,
    0.17096188876030122, 0.16945029123396796, 0.16793829497473117,
    0.1664259035404641, 0.16491312048996992, 0.16339994938297323,
    0.16188639378011183, 0.16037245724292828, 0.15885814333386145,
    0.15734345561623825, 0.15582839765426523, 0.1543129730130201,
    0.15279718525844344, 0.15128103795733022, 0.14976453467732151,
    0.14824767898689603, 0.14673047445536175, 0.14521292465284746,
    0.14369503315029447, 0.14217680351944803, 0.14065823933284921,
    0.1391393441638262, 0.13762012158648604, 0.1361005751757062,
    0.13458070850712617, 0.13306052515713906, 0.13154002870288312,
    0.13001922272223335, 0.12849811079379317, 0.12697669649688587,
    0.12545498341154623, 0.12393297511851216, 0.1224106751992162,
    0.12088808723577708, 0.11936521481099135, 0.11784206150832498,
    0.11631863091190475, 0.11479492660651008, 0.11327095217756435,
    0.11174671121112659, 0.11022220729388306, 0.10869744401313872,
    0.10717242495680884, 0.10564715371341062, 0.10412163387205459,
    0.10259586902243628, 0.10106986275482782, 0.099543618660069319,
    0.0980171403295606, 0.096490431355252593, 0.094963495329638992,
    0.093436335845747787, 0.091908956497132724, 0.090381360877864983,
    0.0888535525825246, 0.087325535206192059, 0.0857973123444399,
    0.084268887593324071, 0.082740264549375692, 0.081211446809592441,
    0.079682437971430126, 0.078153241632794232, 0.076623861392031492,
    0.0750943008479213, 0.073564563599667426, 0.072034653246889332,
    0.070504573389613856, 0.068974327628266746, 0.067443919563664051,
    0.0659133527970038, 0.064382630929857465, 0.0628517575641614,
    0.061320736302208578, 0.059789570746639868, 0.058258264500435752,
    0.056726821166907748, 0.055195244349689941, 0.05366353765273052,
    0.052131704680283324, 0.050599749036899282, 0.049067674327418015,
    0.0475354841569593, 0.046003182130914623, 0.044470771854938668,
    0.04293825693494082, 0.041405640977076739, 0.039872927587739811,
    0.038340120373552694, 0.036807222941358832, 0.035274238898213947,
    0.03374117185137758, 0.032208025408304586, 0.030674803176636626,
    0.029141508764193722, 0.02760814577896574, 0.0260747178291039,
    0.024541228522912288, 0.023007681468839369, 0.021474080275469508,
    0.019940428551514441, 0.01840672990580482, 0.01687298794728171,
    0.0153392062849881, 0.013805388528060391, 0.012271538285719925,
    0.010737659167264491, 0.00920375478205982, 0.007669828739531097,
    0.0061358846491544753, 0.0046019261204485705, 0.0030679567629659761,
    0.0015339801862847655, 0.0, -0.0015339801862847655, -0.0030679567629659761,
    -0.0046019261204485705, -0.0061358846491544753, -0.007669828739531097,
    -0.00920375478205982, -0.010737659167264491, -0.012271538285719925,
    -0.013805388528060391, -0.0153392062849881, -0.01687298794728171,
    -0.01840672990580482, -0.019940428551514441, -0.021474080275469508,
    -0.023007681468839369, -0.024541228522912288, -0.0260747178291039,
    -0.02760814577896574, -0.029141508764193722, -0.030674803176636626,
    -0.032208025408304586, -0.03374117185137758, -0.035274238898213947,
    -0.036807222941358832, -0.038340120373552694, -0.039872927587739811,
    -0.041405640977076739, -0.04293825693494082, -0.044470771854938668,
    -0.046003182130914623, -0.0475354841569593, -0.049067674327418015,
    -0.050599749036899282, -0.052131704680283324, -0.05366353765273052,
    -0.055195244349689941, -0.056726821166907748, -0.058258264500435752,
    -0.059789570746639868, -0.061320736302208578, -0.0628517575641614,
    -0.064382630929857465, -0.0659133527970038, -0.067443919563664051,
    -0.068974327628266746, -0.070504573389613856, -0.072034653246889332,
    -0.073564563599667426, -0.0750943008479213, -0.076623861392031492,
    -0.078153241632794232, -0.079682437971430126, -0.081211446809592441,
    -0.082740264549375692, -0.084268887593324071, -0.0857973123444399,
    -0.087325535206192059, -0.0888535525825246, -0.090381360877864983,
    -0.091908956497132724, -0.093436335845747787, -0.094963495329638992,
    -0.096490431355252593, -0.0980171403295606, -0.099543618660069319,
    -0.10106986275482782, -0.10259586902243628, -0.10412163387205459,
    -0.10564715371341062, -0.10717242495680884, -0.10869744401313872,
    -0.11022220729388306, -0.11174671121112659, -0.11327095217756435,
    -0.11479492660651008, -0.11631863091190475, -0.11784206150832498,
    -0.11936521481099135, -0.12088808723577708, -0.1224106751992162,
    -0.12393297511851216, -0.12545498341154623, -0.12697669649688587,
    -0.12849811079379317, -0.13001922272223335, -0.13154002870288312,
    -0.13306052515713906, -0.13458070850712617, -0.1361005751757062,
    -0.13762012158648604, -0.1391393441638262, -0.14065823933284921,
    -0.14217680351944803, -0.14369503315029447, -0.14521292465284746,
    -0.14673047445536175, -0.14824767898689603, -0.14976453467732151,
    -0.15128103795733022, -0.15279718525844344, -0.1543129730130201,
    -0.15582839765426523, -0.15734345561623825, -0.15885814333386145,
    -0.16037245724292828, -0.16188639378011183, -0.16339994938297323,
    -0.16491312048996992, -0.1664259035404641, -0.16793829497473117,
    -0.16945029123396796, -0.17096188876030122, -0.17247308399679595,
    -0.17398387338746382, -0.17549425337727143, -0.17700422041214875,
    -0.17851377093899751, -0.18002290140569951, -0.18153160826112497,
    -0.18303988795514095, -0.18454773693861962, -0.18605515166344663,
    -0.1875621285825296, -0.18906866414980619, -0.19057475482025274,
    -0.19208039704989244, -0.19358558729580361, -0.19509032201612825,
    -0.19659459767008022, -0.19809841071795356, -0.19960175762113097,
    -0.2011046348420919, -0.20260703884442113, -0.20410896609281687,
    -0.20561041305309924, -0.20711137619221856, -0.20861185197826349,
    -0.21011183688046961, -0.21161132736922755, -0.21311031991609136,
    -0.21460881099378676, -0.21610679707621952, -0.21760427463848364,
    -0.2191012401568698, -0.22059769010887351, -0.22209362097320351,
    -0.22358902922979, -0.22508391135979283, -0.22657826384561,
    -0.22807208317088573, -0.22956536582051887, -0.23105810828067111,
    -0.23255030703877524, -0.23404195858354343, -0.23553305940497549,
    -0.2370236059943672, -0.23851359484431842, -0.2400030224487415,
    -0.24149188530286933, -0.24298017990326387, -0.24446790274782415,
    -0.24595505033579459, -0.24744161916777327, -0.24892760574572015,
    -0.25041300657296522, -0.25189781815421697, -0.25338203699557016,
    -0.25486565960451457, -0.25634868248994291, -0.257831102162159,
    -0.25931291513288623, -0.26079411791527551, -0.26227470702391359,
    -0.26375467897483135, -0.26523403028551179, -0.26671275747489837,
    -0.26819085706340318, -0.26966832557291509, -0.271145159526808,
    -0.272621355449949, -0.27409690986870638, -0.27557181931095814,
    -0.2770460803060999, -0.27851968938505306, -0.27999264308027322,
    -0.28146493792575794, -0.28293657045705539, -0.28440753721127188,
    -0.28587783472708062, -0.28734745954472951, -0.28881640820604948,
    -0.29028467725446233, -0.29175226323498926, -0.29321916269425863,
    -0.29468537218051433, -0.29615088824362379, -0.2976157074350862,
    -0.29907982630804048, -0.30054324141727345, -0.30200594931922808,
    -0.30346794657201132, -0.30492922973540237, -0.30638979537086092,
    -0.30784964004153487, -0.30930876031226873, -0.31076715274961147,
    -0.31222481392182488, -0.31368174039889152, -0.31513792875252244,
    -0.31659337555616585, -0.31804807738501495, -0.31950203081601569,
    -0.32095523242787521, -0.32240767880106985, -0.32385936651785285,
    -0.32531029216226293, -0.32676045232013173, -0.3282098435790925,
    -0.32965846252858749, -0.33110630575987643, -0.33255336986604422,
    -0.33399965144200938, -0.3354451470845316, -0.33688985339222005,
    -0.33833376696554113, -0.33977688440682685, -0.34121920232028236,
    -0.34266071731199438, -0.34410142598993881, -0.34554132496398909,
    -0.34698041084592368, -0.34841868024943456, -0.34985612979013492,
    -0.35129275608556709, -0.35272855575521073, -0.35416352542049034,
    -0.35559766170478385, -0.35703096123343, -0.35846342063373654,
    -0.35989503653498811, -0.36132580556845428, -0.36275572436739723,
    -0.36418478956707989, -0.36561299780477385, -0.36704034571976718,
    -0.36846682995337232, -0.3698924471489341, -0.37131719395183749,
    -0.37274106700951576, -0.37416406297145793, -0.37558617848921722,
    -0.37700741021641826, -0.37842775480876556, -0.37984720892405116,
    -0.38126576922216238, -0.38268343236508978, -0.38410019501693504,
    -0.38551605384391885, -0.38693100551438858, -0.38834504669882625,
    -0.38975817406985641, -0.39117038430225387, -0.39258167407295147,
    -0.3939920400610481, -0.39540147894781635, -0.39680998741671031,
    -0.39821756215337356, -0.39962419984564679, -0.40102989718357562,
    -0.40243465085941843, -0.40383845756765407, -0.40524131400498986,
    -0.40664321687036903, -0.40804416286497869, -0.40944414869225759,
    -0.41084317105790391, -0.41224122666988289, -0.4136383122384345,
    -0.41503442447608163, -0.41642956009763715, -0.41782371582021227,
    -0.41921688836322391, -0.42060907444840251, -0.42200027079979968,
    -0.42339047414379605, -0.42477968120910881, -0.42616788872679962,
    -0.42755509343028208, -0.42894129205532949, -0.43032648134008261,
    -0.43171065802505726, -0.43309381885315196, -0.43447596056965565,
    -0.43585707992225547, -0.43723717366104409, -0.43861623853852766,
    -0.43999427130963326, -0.44137126873171667, -0.44274722756457,
    -0.4441221445704292, -0.44549601651398174, -0.44686884016237416,
    -0.44824061228521989, -0.44961132965460654, -0.45098098904510386,
    -0.45234958723377089, -0.45371712100016387, -0.45508358712634384,
    -0.45644898239688392, -0.45781330359887717, -0.45917654752194409,
    -0.46053871095824, -0.46189979070246273, -0.46325978355186015,
    -0.46461868630623782, -0.46597649576796618, -0.46733320874198842,
    -0.46868882203582796, -0.47004333245959562, -0.47139673682599764,
    -0.47274903195034279, -0.47410021465054997, -0.47545028174715587,
    -0.47679923006332209, -0.478147056424843, -0.47949375766015295,
    -0.48083933060033396, -0.48218377207912272, -0.48352707893291874,
    -0.48486924800079106, -0.48621027612448642, -0.487550160148436,
    -0.48888889691976317, -0.49022648328829116, -0.4915629161065499,
    -0.49289819222978404, -0.49423230851595967, -0.49556526182577254,
    -0.49689704902265447, -0.49822766697278181, -0.49955711254508184,
    -0.50088538261124071, -0.50221247404571079, -0.50353838372571758,
    -0.50486310853126759, -0.50618664534515523, -0.50750899105297087,
    -0.508830142543107, -0.51015009670676681, -0.5114688504379703,
    -0.512786400633563, -0.51410274419322166, -0.51541787801946293,
    -0.51673179901764987, -0.51804450409599934, -0.51935599016558964,
    -0.52066625414036716, -0.52197529293715439, -0.52328310347565643,
    -0.524589682678469, -0.52589502747108463, -0.52719913478190128,
    -0.52850200154222848, -0.52980362468629461, -0.531104001151255,
    -0.5324031278771979, -0.533701001807153, -0.53499761988709715,
    -0.53629297906596318, -0.53758707629564539, -0.53887990853100842,
    -0.54017147272989285, -0.54146176585312344, -0.54275078486451589,
    -0.54403852673088382, -0.54532498842204646, -0.54661016691083486,
    -0.54789405917310019, -0.54917666218771966, -0.55045797293660481,
    -0.55173798840470734, -0.55301670558002747, -0.55429412145362,
    -0.55557023301960218, -0.5568450372751601, -0.5581185312205561,
    -0.55939071185913614, -0.560661576197336, -0.56193112124468936,
    -0.56319934401383409, -0.5644662415205195, -0.56573181078361312,
    -0.56699604882510868, -0.56825895267013149, -0.56952051934694714,
    -0.57078074588696726, -0.572039629324757, -0.5732971666980422,
    -0.57455335504771576, -0.57580819141784534, -0.57706167285567944,
    -0.57831379641165559, -0.57956455913940563, -0.58081395809576453,
    -0.58206199034077544, -0.58330865293769829, -0.58455394295301533,
    -0.58579785745643886, -0.587040393520918, -0.58828154822264522,
    -0.58952131864106394, -0.59075970185887416, -0.591996694962041,
    -0.5932322950397998, -0.59446649918466443, -0.59569930449243336,
    -0.59693070806219639, -0.59816070699634238, -0.59938929840056454,
    -0.600616479383869, -0.60184224705858, -0.60306659854034816,
    -0.604289530948156, -0.60551104140432555, -0.60673112703452448,
    -0.60794978496777363, -0.60916701233645321, -0.61038280627630948,
    -0.61159716392646191, -0.61281008242940971, -0.61402155893103849,
    -0.61523159058062682, -0.61644017453085365, -0.61764730793780387,
    -0.61885298796097632, -0.6200572117632891, -0.62125997651108755,
    -0.62246127937415, -0.62366111752569453, -0.62485948814238634,
    -0.62605638840434352, -0.62725181549514408, -0.6284457666018326,
    -0.629638238914927, -0.63082922962842447, -0.63201873593980906,
    -0.63320675505005719, -0.63439328416364549, -0.63557832048855611,
    -0.6367618612362842, -0.637943903621844, -0.63912444486377573,
    -0.64030348218415167, -0.641481012808583, -0.64265703396622686,
    -0.64383154288979139, -0.64500453681554393, -0.64617601298331628,
    -0.64734596863651206, -0.64851440102211244, -0.64968130739068319,
    -0.650846684996381, -0.6520105310969595, -0.65317284295377676,
    -0.65433361783180044, -0.65549285299961535, -0.65665054572942894,
    -0.65780669329707864, -0.65896129298203732, -0.66011434206742048,
    -0.66126583783999227, -0.66241577759017178, -0.66356415861203977,
    -0.66471097820334479, -0.66585623366550972, -0.66699992230363747,
    -0.66814204142651845, -0.669282588346636, -0.67042156038017309,
    -0.67155895484701833, -0.67269476907077286, -0.673829000378756,
    -0.674961646102012, -0.67609270357531592, -0.67722217013718033,
    -0.67835004312986147, -0.679476319899365, -0.680600997795453,
    -0.68172407417164971, -0.68284554638524808, -0.6839654117973154,
    -0.68508366777270036, -0.68620031168003859, -0.687315340891759,
    -0.68842875278409044, -0.68954054473706683, -0.6906507141345346,
    -0.69175925836415775, -0.69286617481742463, -0.69397146088965389,
    -0.69507511398000088, -0.696177131491463, -0.69727751083088652,
    -0.69837624940897292, -0.69947334464028377, -0.70056879394324834,
    -0.70166259474016845, -0.7027547444572253, -0.70384524052448494,
    -0.70493408037590488, -0.70602126144933974, -0.70710678118654757,
    -0.7081906370331954, -0.70927282643886569, -0.71035334685706242,
    -0.71143219574521643, -0.71250937056469243, -0.71358486878079352,
    -0.71465868786276909, -0.71573082528381859, -0.71680127852109954,
    -0.71787004505573171, -0.71893712237280449, -0.72000250796138165,
    -0.72106619931450811, -0.72212819392921535, -0.72318848930652746,
    -0.724247082951467, -0.72530397237306077, -0.726359155084346,
    -0.72741262860237577, -0.7284643904482252, -0.729514438146997,
    -0.73056276922782759, -0.73160938122389263, -0.73265427167241282,
    -0.73369743811466037, -0.7347388780959635, -0.73577858916571359,
    -0.73681656887736979, -0.737852814788466, -0.73888732446061511,
    -0.7399200954595162, -0.74095112535495922, -0.74198041172083107,
    -0.74300795213512172, -0.74403374417992929, -0.745057785441466,
    -0.74608007351006378, -0.74710060598018013, -0.7481193804504036,
    -0.74913639452345937, -0.75015164580621507, -0.75116513190968637,
    -0.7521768504490427, -0.75318679904361252, -0.75419497531688917,
    -0.75520137689653655, -0.75620600141439454, -0.75720884650648457,
    -0.75820990981301528, -0.759209188978388, -0.76020668165120242,
    -0.76120238548426178, -0.7621962981345789, -0.76318841726338127,
    -0.76417874053611679, -0.765167265622459, -0.76615399019631292,
    -0.7671389119358204, -0.76812202852336542, -0.7691033376455797,
    -0.7700828369933479, -0.77106052426181382, -0.77203639715038452,
    -0.773010453362737, -0.7739826906068229, -0.77495310659487393,
    -0.77592169904340769, -0.77688846567323244, -0.77785340420945315,
    -0.778816512381476, -0.77977778792301455, -0.78073722857209449,
    -0.78169483207105939, -0.78265059616657573, -0.78360451860963831,
    -0.78455659715557524, -0.78550682956405393, -0.78645521359908577,
    -0.78740174702903143, -0.78834642762660634, -0.78928925316888565,
    -0.79023022143731, -0.7911693302176902, -0.79210657730021239,
    -0.79304196047944364, -0.79397547755433717, -0.794907126328237,
    -0.79583690460888357, -0.79676481020841883, -0.79769084094339116,
    -0.79861499463476093, -0.799537269107905, -0.80045766219262282,
    -0.80137617172314024, -0.80229279553811572, -0.80320753148064494,
    -0.8041203773982657, -0.80503133114296366, -0.80594039057117628,
    -0.80684755354379933, -0.80775281792619036, -0.808656181588175,
    -0.80955764240405126, -0.81045719825259477, -0.81135484701706373,
    -0.81225058658520388, -0.81314441484925359, -0.81403632970594841,
    -0.81492632905652662, -0.81581441080673378, -0.81670057286682785,
    -0.81758481315158371, -0.81846712958029866, -0.819347520076797,
    -0.82022598256943469, -0.82110251499110465, -0.82197711527924155,
    -0.82284978137582643, -0.82372051122739143, -0.82458930278502529,
    -0.82545615400437755, -0.82632106284566353, -0.82718402727366913,
    -0.8280450452577558, -0.82890411477186487, -0.829761233794523,
    -0.83061640030884631, -0.83146961230254524, -0.83232086776792968,
    -0.83317016470191319, -0.83401750110601813, -0.83486287498638,
    -0.8357062843537526, -0.836547727223512, -0.83738720161566194,
    -0.83822470555483808, -0.83906023707031274, -0.83989379419599952,
    -0.84072537497045807, -0.84155497743689844, -0.84238259964318585,
    -0.84320823964184544, -0.84403189549006641, -0.84485356524970712,
    -0.84567324698729907, -0.84649093877405213, -0.84730663868585832,
    -0.84812034480329723, -0.84893205521163961, -0.84974176800085255,
    -0.85054948126560348, -0.8513551931052652, -0.85215890162391983,
    -0.85296060493036363, -0.85376030113811141, -0.85455798836540053,
    -0.855353664735196, -0.85614732837519447, -0.85693897741782876,
    -0.85772861000027212, -0.85851622426444274, -0.85930181835700847,
    -0.86008539042939014, -0.86086693863776731, -0.8616464611430813,
    -0.8624239561110405, -0.86319942171212416, -0.8639728561215867,
    -0.86474425751946238, -0.86551362409056909, -0.866280954024513,
    -0.86704624551569265, -0.86780949676330332, -0.8685707059713409,
    -0.86932987134860684, -0.87008699110871146, -0.870842063470079,
    -0.87159508665595109, -0.87234605889439154, -0.87309497841829009,
    -0.87384184346536686, -0.87458665227817611, -0.87532940310411089,
    -0.8760700941954066, -0.87680872380914565, -0.87754529020726135,
    -0.87827979165654158, -0.87901222642863353, -0.87974259280004741,
    -0.88047088905216075, -0.88119711347122209, -0.881921264348355,
    -0.88264333997956279, -0.88336333866573158, -0.884081258712635,
    -0.88479709843093779, -0.8855108561362, -0.88622253014888064,
    -0.88693211879434219, -0.88763962040285393, -0.88834503330959635,
    -0.88904835585466457, -0.88974958638307278, -0.89044872324475788,
    -0.89114576479458318, -0.89184070939234272, -0.89253355540276458,
    -0.89322430119551532, -0.89391294514520325, -0.8945994856313827,
    -0.89528392103855758, -0.89596624975618522, -0.89664647017868015,
    -0.89732458070541832, -0.89800057974073988, -0.89867446569395382,
    -0.89934623697934157, -0.90001589201616017, -0.900683429228647,
    -0.901348847046022, -0.90201214390249318, -0.90267331823725883,
    -0.90333236849451182, -0.90398929312344334, -0.90464409057824624,
    -0.90529675931811882, -0.90594729780726846, -0.90659570451491533,
    -0.90724197791529582, -0.90788611648766626, -0.90852811871630612,
    -0.90916798309052238, -0.90980570810465222, -0.91044129225806725,
    -0.91107473405517636, -0.91170603200542988, -0.91233518462332275,
    -0.91296219042839821, -0.91358704794525081, -0.91420975570353069,
    -0.9148303122379462, -0.91544871608826783, -0.91606496579933172,
    -0.9166790599210427, -0.91729099700837791, -0.9179007756213905,
    -0.91850839432521225, -0.91911385169005777, -0.91971714629122736,
    -0.92031827670911059, -0.92091724152918952, -0.9215140393420419,
    -0.92210866874334518, -0.92270112833387863, -0.92329141671952764,
    -0.92387953251128674, -0.9244654743252626, -0.92504924078267758,
    -0.92563083050987272, -0.92621024213831138, -0.92678747430458175,
    -0.92736252565040111, -0.92793539482261789, -0.92850608047321559,
    -0.92907458125931586, -0.92964089584318121, -0.93020502289221907,
    -0.93076696107898371, -0.93132670908118043, -0.93188426558166815,
    -0.93243962926846236, -0.932992798834739, -0.93354377297883617,
    -0.93409255040425887, -0.93463912981968078, -0.93518350993894761,
    -0.93572568948108037, -0.93626566717027826, -0.93680344173592156,
    -0.937339011912575, -0.93787237643998989, -0.93840353406310806,
    -0.9389324835320646, -0.93945922360218992, -0.93998375303401394,
    -0.9405060705932683, -0.94102617505088926, -0.94154406518302081,
    -0.94205973977101731, -0.94257319760144687, -0.94308443746609349,
    -0.94359345816196039, -0.94410025849127266, -0.94460483726148026,
    -0.94510719328526061, -0.94560732538052128, -0.94610523237040345,
    -0.94660091308328353, -0.94709436635277722, -0.94758559101774109,
    -0.94807458592227623, -0.94856134991573027, -0.94904588185270056,
    -0.94952818059303667, -0.950008245001843, -0.9504860739494817,
    -0.95096166631157508, -0.95143502096900834, -0.95190613680793235,
    -0.95237501271976588, -0.95284164760119872, -0.95330604035419386,
    -0.95376818988599033, -0.95422809510910567, -0.95468575494133834,
    -0.95514116830577078, -0.95559433413077111, -0.95604525134999641,
    -0.9564939189023951, -0.95694033573220882, -0.95738450078897586,
    -0.95782641302753291, -0.95826607140801767, -0.9587034748958716,
    -0.95913862246184189, -0.95957151308198452, -0.960002145737666,
    -0.96043051941556579, -0.96085663310767966, -0.96128048581132064,
    -0.96170207652912254, -0.96212140426904158, -0.96253846804435916,
    -0.96295326687368388, -0.963365799780954, -0.96377606579543984,
    -0.96418406395174583, -0.96458979328981276, -0.96499325285492032,
    -0.9653944416976894, -0.96579335887408368, -0.9661900034454125,
    -0.96658437447833312, -0.96697647104485207, -0.96736629222232851,
    -0.96775383709347551, -0.96813910474636244, -0.96852209427441727,
    -0.96890280477642887, -0.96928123535654853, -0.96965738512429245,
    -0.970031253194544, -0.9704028386875555, -0.97077214072895035,
    -0.97113915844972509, -0.97150389098625178, -0.9718663374802794,
    -0.97222649707893627, -0.97258436893473221, -0.97293995220556018,
    -0.97329324605469825, -0.973644249650812, -0.97399296216795583,
    -0.97433938278557586, -0.97468351068851067, -0.97502534506699412,
    -0.975364885116657, -0.97570213003852857, -0.976037079039039,
    -0.97636973133002114, -0.97670008612871184, -0.97702814265775439,
    -0.97735390014520007, -0.97767735782450993, -0.97799851493455714,
    -0.97831737071962765, -0.97863392442942321, -0.9789481753190622,
    -0.979260122649082, -0.97956976568544052, -0.97987710369951764,
    -0.98018213596811743, -0.98048486177346938, -0.98078528040323043,
    -0.98108339115048671, -0.98137919331375456, -0.98167268619698311,
    -0.98196386910955524, -0.98225274136628937, -0.98253930228744124,
    -0.98282355119870524, -0.98310548743121629, -0.98338511032155118,
    -0.98366241921173025, -0.98393741344921892, -0.984210092386929,
    -0.98448045538322093, -0.98474850180190421, -0.98501423101223984,
    -0.98527764238894122, -0.98553873531217606, -0.98579750916756748,
    -0.98605396334619544, -0.98630809724459867, -0.98655991026477541,
    -0.98680940181418553, -0.987056571305751, -0.98730141815785843,
    -0.98754394179435923, -0.98778414164457218, -0.98802201714328353,
    -0.98825756773074946, -0.98849079285269659, -0.98872169196032378,
    -0.988950264510303, -0.989176509964781, -0.98940042779138038,
    -0.98962201746320089, -0.98984127845882053, -0.99005821026229712,
    -0.99027281236316911, -0.99048508425645709, -0.99069502544266463,
    -0.99090263542778, -0.99110791372327689, -0.99131085984611544,
    -0.9915114733187439, -0.99170975366909953, -0.99190570043060933,
    -0.9920993131421918, -0.99229059134825737, -0.99247953459871,
    -0.992666142448948, -0.9928504144598651, -0.99303235019785141,
    -0.9932119492347945, -0.99338921114808065, -0.9935641355205953,
    -0.9937367219407246, -0.99390697000235606, -0.99407487930487937,
    -0.9942404494531879, -0.9944036800576791, -0.99456457073425542,
    -0.9947231211043257, -0.99487933079480562, -0.99503319943811863,
    -0.99518472667219693, -0.99533391214048228, -0.99548075549192694,
    -0.99562525638099431, -0.99576741446765982, -0.99590722941741172,
    -0.996044700901252, -0.996179828595697, -0.996312612182778,
    -0.99644305135004263, -0.99657114579055484, -0.99669689520289606,
    -0.99682029929116567, -0.99694135776498216, -0.997060070339483,
    -0.99717643673532619, -0.99729045667869021, -0.9974021299012753,
    -0.99751145614030345, -0.99761843513851955, -0.99772306664419164,
    -0.99782535041111164, -0.997925286198596, -0.99802287377148624,
    -0.99811811290014918, -0.99821100336047819, -0.99830154493389289,
    -0.99838973740734016, -0.99847558057329477, -0.99855907422975931,
    -0.99864021818026527, -0.99871901223387294, -0.99879545620517241,
    -0.99886954991428356, -0.99894129318685687, -0.99901068585407338,
    -0.99907772775264536, -0.99914241872481691, -0.99920475861836389,
    -0.99926474728659442, -0.99932238458834954, -0.99937767038800285,
    -0.99943060455546173, -0.999481186966167, -0.99952941750109314,
    -0.99957529604674922, -0.99961882249517864, -0.99965999674395922,
    -0.99969881869620425, -0.99973528826056168, -0.99976940535121528,
    -0.99980116988788426, -0.9998305817958234, -0.99985764100582386,
    -0.99988234745421256, -0.9999047010828529, -0.9999247018391445,
    -0.99994234967602391, -0.9999576445519639, -0.99997058643097414,
    -0.99998117528260111, -0.9999894110819284, -0.99999529380957619,
    -0.99999882345170188, -1.0 };

  static const double sintab[2049] = { 0.0, -0.0015339801862847655,
    -0.0030679567629659761, -0.0046019261204485705, -0.0061358846491544753,
    -0.007669828739531097, -0.00920375478205982, -0.010737659167264491,
    -0.012271538285719925, -0.013805388528060391, -0.0153392062849881,
    -0.01687298794728171, -0.01840672990580482, -0.019940428551514441,
    -0.021474080275469508, -0.023007681468839369, -0.024541228522912288,
    -0.0260747178291039, -0.02760814577896574, -0.029141508764193722,
    -0.030674803176636626, -0.032208025408304586, -0.03374117185137758,
    -0.035274238898213947, -0.036807222941358832, -0.038340120373552694,
    -0.039872927587739811, -0.041405640977076739, -0.04293825693494082,
    -0.044470771854938668, -0.046003182130914623, -0.0475354841569593,
    -0.049067674327418015, -0.050599749036899282, -0.052131704680283324,
    -0.05366353765273052, -0.055195244349689941, -0.056726821166907748,
    -0.058258264500435752, -0.059789570746639868, -0.061320736302208578,
    -0.0628517575641614, -0.064382630929857465, -0.0659133527970038,
    -0.067443919563664051, -0.068974327628266746, -0.070504573389613856,
    -0.072034653246889332, -0.073564563599667426, -0.0750943008479213,
    -0.076623861392031492, -0.078153241632794232, -0.079682437971430126,
    -0.081211446809592441, -0.082740264549375692, -0.084268887593324071,
    -0.0857973123444399, -0.087325535206192059, -0.0888535525825246,
    -0.090381360877864983, -0.091908956497132724, -0.093436335845747787,
    -0.094963495329638992, -0.096490431355252593, -0.0980171403295606,
    -0.099543618660069319, -0.10106986275482782, -0.10259586902243628,
    -0.10412163387205459, -0.10564715371341062, -0.10717242495680884,
    -0.10869744401313872, -0.11022220729388306, -0.11174671121112659,
    -0.11327095217756435, -0.11479492660651008, -0.11631863091190475,
    -0.11784206150832498, -0.11936521481099135, -0.12088808723577708,
    -0.1224106751992162, -0.12393297511851216, -0.12545498341154623,
    -0.12697669649688587, -0.12849811079379317, -0.13001922272223335,
    -0.13154002870288312, -0.13306052515713906, -0.13458070850712617,
    -0.1361005751757062, -0.13762012158648604, -0.1391393441638262,
    -0.14065823933284921, -0.14217680351944803, -0.14369503315029447,
    -0.14521292465284746, -0.14673047445536175, -0.14824767898689603,
    -0.14976453467732151, -0.15128103795733022, -0.15279718525844344,
    -0.1543129730130201, -0.15582839765426523, -0.15734345561623825,
    -0.15885814333386145, -0.16037245724292828, -0.16188639378011183,
    -0.16339994938297323, -0.16491312048996992, -0.1664259035404641,
    -0.16793829497473117, -0.16945029123396796, -0.17096188876030122,
    -0.17247308399679595, -0.17398387338746382, -0.17549425337727143,
    -0.17700422041214875, -0.17851377093899751, -0.18002290140569951,
    -0.18153160826112497, -0.18303988795514095, -0.18454773693861962,
    -0.18605515166344663, -0.1875621285825296, -0.18906866414980619,
    -0.19057475482025274, -0.19208039704989244, -0.19358558729580361,
    -0.19509032201612825, -0.19659459767008022, -0.19809841071795356,
    -0.19960175762113097, -0.2011046348420919, -0.20260703884442113,
    -0.20410896609281687, -0.20561041305309924, -0.20711137619221856,
    -0.20861185197826349, -0.21011183688046961, -0.21161132736922755,
    -0.21311031991609136, -0.21460881099378676, -0.21610679707621952,
    -0.21760427463848364, -0.2191012401568698, -0.22059769010887351,
    -0.22209362097320351, -0.22358902922979, -0.22508391135979283,
    -0.22657826384561, -0.22807208317088573, -0.22956536582051887,
    -0.23105810828067111, -0.23255030703877524, -0.23404195858354343,
    -0.23553305940497549, -0.2370236059943672, -0.23851359484431842,
    -0.2400030224487415, -0.24149188530286933, -0.24298017990326387,
    -0.24446790274782415, -0.24595505033579459, -0.24744161916777327,
    -0.24892760574572015, -0.25041300657296522, -0.25189781815421697,
    -0.25338203699557016, -0.25486565960451457, -0.25634868248994291,
    -0.257831102162159, -0.25931291513288623, -0.26079411791527551,
    -0.26227470702391359, -0.26375467897483135, -0.26523403028551179,
    -0.26671275747489837, -0.26819085706340318, -0.26966832557291509,
    -0.271145159526808, -0.272621355449949, -0.27409690986870638,
    -0.27557181931095814, -0.2770460803060999, -0.27851968938505306,
    -0.27999264308027322, -0.28146493792575794, -0.28293657045705539,
    -0.28440753721127188, -0.28587783472708062, -0.28734745954472951,
    -0.28881640820604948, -0.29028467725446233, -0.29175226323498926,
    -0.29321916269425863, -0.29468537218051433, -0.29615088824362379,
    -0.2976157074350862, -0.29907982630804048, -0.30054324141727345,
    -0.30200594931922808, -0.30346794657201132, -0.30492922973540237,
    -0.30638979537086092, -0.30784964004153487, -0.30930876031226873,
    -0.31076715274961147, -0.31222481392182488, -0.31368174039889152,
    -0.31513792875252244, -0.31659337555616585, -0.31804807738501495,
    -0.31950203081601569, -0.32095523242787521, -0.32240767880106985,
    -0.32385936651785285, -0.32531029216226293, -0.32676045232013173,
    -0.3282098435790925, -0.32965846252858749, -0.33110630575987643,
    -0.33255336986604422, -0.33399965144200938, -0.3354451470845316,
    -0.33688985339222005, -0.33833376696554113, -0.33977688440682685,
    -0.34121920232028236, -0.34266071731199438, -0.34410142598993881,
    -0.34554132496398909, -0.34698041084592368, -0.34841868024943456,
    -0.34985612979013492, -0.35129275608556709, -0.35272855575521073,
    -0.35416352542049034, -0.35559766170478385, -0.35703096123343,
    -0.35846342063373654, -0.35989503653498811, -0.36132580556845428,
    -0.36275572436739723, -0.36418478956707989, -0.36561299780477385,
    -0.36704034571976718, -0.36846682995337232, -0.3698924471489341,
    -0.37131719395183749, -0.37274106700951576, -0.37416406297145793,
    -0.37558617848921722, -0.37700741021641826, -0.37842775480876556,
    -0.37984720892405116, -0.38126576922216238, -0.38268343236508978,
    -0.38410019501693504, -0.38551605384391885, -0.38693100551438858,
    -0.38834504669882625, -0.38975817406985641, -0.39117038430225387,
    -0.39258167407295147, -0.3939920400610481, -0.39540147894781635,
    -0.39680998741671031, -0.39821756215337356, -0.39962419984564679,
    -0.40102989718357562, -0.40243465085941843, -0.40383845756765407,
    -0.40524131400498986, -0.40664321687036903, -0.40804416286497869,
    -0.40944414869225759, -0.41084317105790391, -0.41224122666988289,
    -0.4136383122384345, -0.41503442447608163, -0.41642956009763715,
    -0.41782371582021227, -0.41921688836322391, -0.42060907444840251,
    -0.42200027079979968, -0.42339047414379605, -0.42477968120910881,
    -0.42616788872679962, -0.42755509343028208, -0.42894129205532949,
    -0.43032648134008261, -0.43171065802505726, -0.43309381885315196,
    -0.43447596056965565, -0.43585707992225547, -0.43723717366104409,
    -0.43861623853852766, -0.43999427130963326, -0.44137126873171667,
    -0.44274722756457, -0.4441221445704292, -0.44549601651398174,
    -0.44686884016237416, -0.44824061228521989, -0.44961132965460654,
    -0.45098098904510386, -0.45234958723377089, -0.45371712100016387,
    -0.45508358712634384, -0.45644898239688392, -0.45781330359887717,
    -0.45917654752194409, -0.46053871095824, -0.46189979070246273,
    -0.46325978355186015, -0.46461868630623782, -0.46597649576796618,
    -0.46733320874198842, -0.46868882203582796, -0.47004333245959562,
    -0.47139673682599764, -0.47274903195034279, -0.47410021465054997,
    -0.47545028174715587, -0.47679923006332209, -0.478147056424843,
    -0.47949375766015295, -0.48083933060033396, -0.48218377207912272,
    -0.48352707893291874, -0.48486924800079106, -0.48621027612448642,
    -0.487550160148436, -0.48888889691976317, -0.49022648328829116,
    -0.4915629161065499, -0.49289819222978404, -0.49423230851595967,
    -0.49556526182577254, -0.49689704902265447, -0.49822766697278181,
    -0.49955711254508184, -0.50088538261124071, -0.50221247404571079,
    -0.50353838372571758, -0.50486310853126759, -0.50618664534515523,
    -0.50750899105297087, -0.508830142543107, -0.51015009670676681,
    -0.5114688504379703, -0.512786400633563, -0.51410274419322166,
    -0.51541787801946293, -0.51673179901764987, -0.51804450409599934,
    -0.51935599016558964, -0.52066625414036716, -0.52197529293715439,
    -0.52328310347565643, -0.524589682678469, -0.52589502747108463,
    -0.52719913478190128, -0.52850200154222848, -0.52980362468629461,
    -0.531104001151255, -0.5324031278771979, -0.533701001807153,
    -0.53499761988709715, -0.53629297906596318, -0.53758707629564539,
    -0.53887990853100842, -0.54017147272989285, -0.54146176585312344,
    -0.54275078486451589, -0.54403852673088382, -0.54532498842204646,
    -0.54661016691083486, -0.54789405917310019, -0.54917666218771966,
    -0.55045797293660481, -0.55173798840470734, -0.55301670558002747,
    -0.55429412145362, -0.55557023301960218, -0.5568450372751601,
    -0.5581185312205561, -0.55939071185913614, -0.560661576197336,
    -0.56193112124468936, -0.56319934401383409, -0.5644662415205195,
    -0.56573181078361312, -0.56699604882510868, -0.56825895267013149,
    -0.56952051934694714, -0.57078074588696726, -0.572039629324757,
    -0.5732971666980422, -0.57455335504771576, -0.57580819141784534,
    -0.57706167285567944, -0.57831379641165559, -0.57956455913940563,
    -0.58081395809576453, -0.58206199034077544, -0.58330865293769829,
    -0.58455394295301533, -0.58579785745643886, -0.587040393520918,
    -0.58828154822264522, -0.58952131864106394, -0.59075970185887416,
    -0.591996694962041, -0.5932322950397998, -0.59446649918466443,
    -0.59569930449243336, -0.59693070806219639, -0.59816070699634238,
    -0.59938929840056454, -0.600616479383869, -0.60184224705858,
    -0.60306659854034816, -0.604289530948156, -0.60551104140432555,
    -0.60673112703452448, -0.60794978496777363, -0.60916701233645321,
    -0.61038280627630948, -0.61159716392646191, -0.61281008242940971,
    -0.61402155893103849, -0.61523159058062682, -0.61644017453085365,
    -0.61764730793780387, -0.61885298796097632, -0.6200572117632891,
    -0.62125997651108755, -0.62246127937415, -0.62366111752569453,
    -0.62485948814238634, -0.62605638840434352, -0.62725181549514408,
    -0.6284457666018326, -0.629638238914927, -0.63082922962842447,
    -0.63201873593980906, -0.63320675505005719, -0.63439328416364549,
    -0.63557832048855611, -0.6367618612362842, -0.637943903621844,
    -0.63912444486377573, -0.64030348218415167, -0.641481012808583,
    -0.64265703396622686, -0.64383154288979139, -0.64500453681554393,
    -0.64617601298331628, -0.64734596863651206, -0.64851440102211244,
    -0.64968130739068319, -0.650846684996381, -0.6520105310969595,
    -0.65317284295377676, -0.65433361783180044, -0.65549285299961535,
    -0.65665054572942894, -0.65780669329707864, -0.65896129298203732,
    -0.66011434206742048, -0.66126583783999227, -0.66241577759017178,
    -0.66356415861203977, -0.66471097820334479, -0.66585623366550972,
    -0.66699992230363747, -0.66814204142651845, -0.669282588346636,
    -0.67042156038017309, -0.67155895484701833, -0.67269476907077286,
    -0.673829000378756, -0.674961646102012, -0.67609270357531592,
    -0.67722217013718033, -0.67835004312986147, -0.679476319899365,
    -0.680600997795453, -0.68172407417164971, -0.68284554638524808,
    -0.6839654117973154, -0.68508366777270036, -0.68620031168003859,
    -0.687315340891759, -0.68842875278409044, -0.68954054473706683,
    -0.6906507141345346, -0.69175925836415775, -0.69286617481742463,
    -0.69397146088965389, -0.69507511398000088, -0.696177131491463,
    -0.69727751083088652, -0.69837624940897292, -0.69947334464028377,
    -0.70056879394324834, -0.70166259474016845, -0.7027547444572253,
    -0.70384524052448494, -0.70493408037590488, -0.70602126144933974,
    -0.70710678118654757, -0.7081906370331954, -0.70927282643886569,
    -0.71035334685706242, -0.71143219574521643, -0.71250937056469243,
    -0.71358486878079352, -0.71465868786276909, -0.71573082528381859,
    -0.71680127852109954, -0.71787004505573171, -0.71893712237280449,
    -0.72000250796138165, -0.72106619931450811, -0.72212819392921535,
    -0.72318848930652746, -0.724247082951467, -0.72530397237306077,
    -0.726359155084346, -0.72741262860237577, -0.7284643904482252,
    -0.729514438146997, -0.73056276922782759, -0.73160938122389263,
    -0.73265427167241282, -0.73369743811466037, -0.7347388780959635,
    -0.73577858916571359, -0.73681656887736979, -0.737852814788466,
    -0.73888732446061511, -0.7399200954595162, -0.74095112535495922,
    -0.74198041172083107, -0.74300795213512172, -0.74403374417992929,
    -0.745057785441466, -0.74608007351006378, -0.74710060598018013,
    -0.7481193804504036, -0.74913639452345937, -0.75015164580621507,
    -0.75116513190968637, -0.7521768504490427, -0.75318679904361252,
    -0.75419497531688917, -0.75520137689653655, -0.75620600141439454,
    -0.75720884650648457, -0.75820990981301528, -0.759209188978388,
    -0.76020668165120242, -0.76120238548426178, -0.7621962981345789,
    -0.76318841726338127, -0.76417874053611679, -0.765167265622459,
    -0.76615399019631292, -0.7671389119358204, -0.76812202852336542,
    -0.7691033376455797, -0.7700828369933479, -0.77106052426181382,
    -0.77203639715038452, -0.773010453362737, -0.7739826906068229,
    -0.77495310659487393, -0.77592169904340769, -0.77688846567323244,
    -0.77785340420945315, -0.778816512381476, -0.77977778792301455,
    -0.78073722857209449, -0.78169483207105939, -0.78265059616657573,
    -0.78360451860963831, -0.78455659715557524, -0.78550682956405393,
    -0.78645521359908577, -0.78740174702903143, -0.78834642762660634,
    -0.78928925316888565, -0.79023022143731, -0.7911693302176902,
    -0.79210657730021239, -0.79304196047944364, -0.79397547755433717,
    -0.794907126328237, -0.79583690460888357, -0.79676481020841883,
    -0.79769084094339116, -0.79861499463476093, -0.799537269107905,
    -0.80045766219262282, -0.80137617172314024, -0.80229279553811572,
    -0.80320753148064494, -0.8041203773982657, -0.80503133114296366,
    -0.80594039057117628, -0.80684755354379933, -0.80775281792619036,
    -0.808656181588175, -0.80955764240405126, -0.81045719825259477,
    -0.81135484701706373, -0.81225058658520388, -0.81314441484925359,
    -0.81403632970594841, -0.81492632905652662, -0.81581441080673378,
    -0.81670057286682785, -0.81758481315158371, -0.81846712958029866,
    -0.819347520076797, -0.82022598256943469, -0.82110251499110465,
    -0.82197711527924155, -0.82284978137582643, -0.82372051122739143,
    -0.82458930278502529, -0.82545615400437755, -0.82632106284566353,
    -0.82718402727366913, -0.8280450452577558, -0.82890411477186487,
    -0.829761233794523, -0.83061640030884631, -0.83146961230254524,
    -0.83232086776792968, -0.83317016470191319, -0.83401750110601813,
    -0.83486287498638, -0.8357062843537526, -0.836547727223512,
    -0.83738720161566194, -0.83822470555483808, -0.83906023707031274,
    -0.83989379419599952, -0.84072537497045807, -0.84155497743689844,
    -0.84238259964318585, -0.84320823964184544, -0.84403189549006641,
    -0.84485356524970712, -0.84567324698729907, -0.84649093877405213,
    -0.84730663868585832, -0.84812034480329723, -0.84893205521163961,
    -0.84974176800085255, -0.85054948126560348, -0.8513551931052652,
    -0.85215890162391983, -0.85296060493036363, -0.85376030113811141,
    -0.85455798836540053, -0.855353664735196, -0.85614732837519447,
    -0.85693897741782876, -0.85772861000027212, -0.85851622426444274,
    -0.85930181835700847, -0.86008539042939014, -0.86086693863776731,
    -0.8616464611430813, -0.8624239561110405, -0.86319942171212416,
    -0.8639728561215867, -0.86474425751946238, -0.86551362409056909,
    -0.866280954024513, -0.86704624551569265, -0.86780949676330332,
    -0.8685707059713409, -0.86932987134860684, -0.87008699110871146,
    -0.870842063470079, -0.87159508665595109, -0.87234605889439154,
    -0.87309497841829009, -0.87384184346536686, -0.87458665227817611,
    -0.87532940310411089, -0.8760700941954066, -0.87680872380914565,
    -0.87754529020726135, -0.87827979165654158, -0.87901222642863353,
    -0.87974259280004741, -0.88047088905216075, -0.88119711347122209,
    -0.881921264348355, -0.88264333997956279, -0.88336333866573158,
    -0.884081258712635, -0.88479709843093779, -0.8855108561362,
    -0.88622253014888064, -0.88693211879434219, -0.88763962040285393,
    -0.88834503330959635, -0.88904835585466457, -0.88974958638307278,
    -0.89044872324475788, -0.89114576479458318, -0.89184070939234272,
    -0.89253355540276458, -0.89322430119551532, -0.89391294514520325,
    -0.8945994856313827, -0.89528392103855758, -0.89596624975618522,
    -0.89664647017868015, -0.89732458070541832, -0.89800057974073988,
    -0.89867446569395382, -0.89934623697934157, -0.90001589201616017,
    -0.900683429228647, -0.901348847046022, -0.90201214390249318,
    -0.90267331823725883, -0.90333236849451182, -0.90398929312344334,
    -0.90464409057824624, -0.90529675931811882, -0.90594729780726846,
    -0.90659570451491533, -0.90724197791529582, -0.90788611648766626,
    -0.90852811871630612, -0.90916798309052238, -0.90980570810465222,
    -0.91044129225806725, -0.91107473405517636, -0.91170603200542988,
    -0.91233518462332275, -0.91296219042839821, -0.91358704794525081,
    -0.91420975570353069, -0.9148303122379462, -0.91544871608826783,
    -0.91606496579933172, -0.9166790599210427, -0.91729099700837791,
    -0.9179007756213905, -0.91850839432521225, -0.91911385169005777,
    -0.91971714629122736, -0.92031827670911059, -0.92091724152918952,
    -0.9215140393420419, -0.92210866874334518, -0.92270112833387863,
    -0.92329141671952764, -0.92387953251128674, -0.9244654743252626,
    -0.92504924078267758, -0.92563083050987272, -0.92621024213831138,
    -0.92678747430458175, -0.92736252565040111, -0.92793539482261789,
    -0.92850608047321559, -0.92907458125931586, -0.92964089584318121,
    -0.93020502289221907, -0.93076696107898371, -0.93132670908118043,
    -0.93188426558166815, -0.93243962926846236, -0.932992798834739,
    -0.93354377297883617, -0.93409255040425887, -0.93463912981968078,
    -0.93518350993894761, -0.93572568948108037, -0.93626566717027826,
    -0.93680344173592156, -0.937339011912575, -0.93787237643998989,
    -0.93840353406310806, -0.9389324835320646, -0.93945922360218992,
    -0.93998375303401394, -0.9405060705932683, -0.94102617505088926,
    -0.94154406518302081, -0.94205973977101731, -0.94257319760144687,
    -0.94308443746609349, -0.94359345816196039, -0.94410025849127266,
    -0.94460483726148026, -0.94510719328526061, -0.94560732538052128,
    -0.94610523237040345, -0.94660091308328353, -0.94709436635277722,
    -0.94758559101774109, -0.94807458592227623, -0.94856134991573027,
    -0.94904588185270056, -0.94952818059303667, -0.950008245001843,
    -0.9504860739494817, -0.95096166631157508, -0.95143502096900834,
    -0.95190613680793235, -0.95237501271976588, -0.95284164760119872,
    -0.95330604035419386, -0.95376818988599033, -0.95422809510910567,
    -0.95468575494133834, -0.95514116830577078, -0.95559433413077111,
    -0.95604525134999641, -0.9564939189023951, -0.95694033573220882,
    -0.95738450078897586, -0.95782641302753291, -0.95826607140801767,
    -0.9587034748958716, -0.95913862246184189, -0.95957151308198452,
    -0.960002145737666, -0.96043051941556579, -0.96085663310767966,
    -0.96128048581132064, -0.96170207652912254, -0.96212140426904158,
    -0.96253846804435916, -0.96295326687368388, -0.963365799780954,
    -0.96377606579543984, -0.96418406395174583, -0.96458979328981276,
    -0.96499325285492032, -0.9653944416976894, -0.96579335887408368,
    -0.9661900034454125, -0.96658437447833312, -0.96697647104485207,
    -0.96736629222232851, -0.96775383709347551, -0.96813910474636244,
    -0.96852209427441727, -0.96890280477642887, -0.96928123535654853,
    -0.96965738512429245, -0.970031253194544, -0.9704028386875555,
    -0.97077214072895035, -0.97113915844972509, -0.97150389098625178,
    -0.9718663374802794, -0.97222649707893627, -0.97258436893473221,
    -0.97293995220556018, -0.97329324605469825, -0.973644249650812,
    -0.97399296216795583, -0.97433938278557586, -0.97468351068851067,
    -0.97502534506699412, -0.975364885116657, -0.97570213003852857,
    -0.976037079039039, -0.97636973133002114, -0.97670008612871184,
    -0.97702814265775439, -0.97735390014520007, -0.97767735782450993,
    -0.97799851493455714, -0.97831737071962765, -0.97863392442942321,
    -0.9789481753190622, -0.979260122649082, -0.97956976568544052,
    -0.97987710369951764, -0.98018213596811743, -0.98048486177346938,
    -0.98078528040323043, -0.98108339115048671, -0.98137919331375456,
    -0.98167268619698311, -0.98196386910955524, -0.98225274136628937,
    -0.98253930228744124, -0.98282355119870524, -0.98310548743121629,
    -0.98338511032155118, -0.98366241921173025, -0.98393741344921892,
    -0.984210092386929, -0.98448045538322093, -0.98474850180190421,
    -0.98501423101223984, -0.98527764238894122, -0.98553873531217606,
    -0.98579750916756748, -0.98605396334619544, -0.98630809724459867,
    -0.98655991026477541, -0.98680940181418553, -0.987056571305751,
    -0.98730141815785843, -0.98754394179435923, -0.98778414164457218,
    -0.98802201714328353, -0.98825756773074946, -0.98849079285269659,
    -0.98872169196032378, -0.988950264510303, -0.989176509964781,
    -0.98940042779138038, -0.98962201746320089, -0.98984127845882053,
    -0.99005821026229712, -0.99027281236316911, -0.99048508425645709,
    -0.99069502544266463, -0.99090263542778, -0.99110791372327689,
    -0.99131085984611544, -0.9915114733187439, -0.99170975366909953,
    -0.99190570043060933, -0.9920993131421918, -0.99229059134825737,
    -0.99247953459871, -0.992666142448948, -0.9928504144598651,
    -0.99303235019785141, -0.9932119492347945, -0.99338921114808065,
    -0.9935641355205953, -0.9937367219407246, -0.99390697000235606,
    -0.99407487930487937, -0.9942404494531879, -0.9944036800576791,
    -0.99456457073425542, -0.9947231211043257, -0.99487933079480562,
    -0.99503319943811863, -0.99518472667219693, -0.99533391214048228,
    -0.99548075549192694, -0.99562525638099431, -0.99576741446765982,
    -0.99590722941741172, -0.996044700901252, -0.996179828595697,
    -0.996312612182778, -0.99644305135004263, -0.99657114579055484,
    -0.99669689520289606, -0.99682029929116567, -0.99694135776498216,
    -0.997060070339483, -0.99717643673532619, -0.99729045667869021,
    -0.9974021299012753, -0.99751145614030345, -0.99761843513851955,
    -0.99772306664419164, -0.99782535041111164, -0.997925286198596,
    -0.99802287377148624, -0.99811811290014918, -0.99821100336047819,
    -0.99830154493389289, -0.99838973740734016, -0.99847558057329477,
    -0.99855907422975931, -0.99864021818026527, -0.99871901223387294,
    -0.99879545620517241, -0.99886954991428356, -0.99894129318685687,
    -0.99901068585407338, -0.99907772775264536, -0.99914241872481691,
    -0.99920475861836389, -0.99926474728659442, -0.99932238458834954,
    -0.99937767038800285, -0.99943060455546173, -0.999481186966167,
    -0.99952941750109314, -0.99957529604674922, -0.99961882249517864,
    -0.99965999674395922, -0.99969881869620425, -0.99973528826056168,
    -0.99976940535121528, -0.99980116988788426, -0.9998305817958234,
    -0.99985764100582386, -0.99988234745421256, -0.9999047010828529,
    -0.9999247018391445, -0.99994234967602391, -0.9999576445519639,
    -0.99997058643097414, -0.99998117528260111, -0.9999894110819284,
    -0.99999529380957619, -0.99999882345170188, -1.0, -0.99999882345170188,
    -0.99999529380957619, -0.9999894110819284, -0.99998117528260111,
    -0.99997058643097414, -0.9999576445519639, -0.99994234967602391,
    -0.9999247018391445, -0.9999047010828529, -0.99988234745421256,
    -0.99985764100582386, -0.9998305817958234, -0.99980116988788426,
    -0.99976940535121528, -0.99973528826056168, -0.99969881869620425,
    -0.99965999674395922, -0.99961882249517864, -0.99957529604674922,
    -0.99952941750109314, -0.999481186966167, -0.99943060455546173,
    -0.99937767038800285, -0.99932238458834954, -0.99926474728659442,
    -0.99920475861836389, -0.99914241872481691, -0.99907772775264536,
    -0.99901068585407338, -0.99894129318685687, -0.99886954991428356,
    -0.99879545620517241, -0.99871901223387294, -0.99864021818026527,
    -0.99855907422975931, -0.99847558057329477, -0.99838973740734016,
    -0.99830154493389289, -0.99821100336047819, -0.99811811290014918,
    -0.99802287377148624, -0.997925286198596, -0.99782535041111164,
    -0.99772306664419164, -0.99761843513851955, -0.99751145614030345,
    -0.9974021299012753, -0.99729045667869021, -0.99717643673532619,
    -0.997060070339483, -0.99694135776498216, -0.99682029929116567,
    -0.99669689520289606, -0.99657114579055484, -0.99644305135004263,
    -0.996312612182778, -0.996179828595697, -0.996044700901252,
    -0.99590722941741172, -0.99576741446765982, -0.99562525638099431,
    -0.99548075549192694, -0.99533391214048228, -0.99518472667219693,
    -0.99503319943811863, -0.99487933079480562, -0.9947231211043257,
    -0.99456457073425542, -0.9944036800576791, -0.9942404494531879,
    -0.99407487930487937, -0.99390697000235606, -0.9937367219407246,
    -0.9935641355205953, -0.99338921114808065, -0.9932119492347945,
    -0.99303235019785141, -0.9928504144598651, -0.992666142448948,
    -0.99247953459871, -0.99229059134825737, -0.9920993131421918,
    -0.99190570043060933, -0.99170975366909953, -0.9915114733187439,
    -0.99131085984611544, -0.99110791372327689, -0.99090263542778,
    -0.99069502544266463, -0.99048508425645709, -0.99027281236316911,
    -0.99005821026229712, -0.98984127845882053, -0.98962201746320089,
    -0.98940042779138038, -0.989176509964781, -0.988950264510303,
    -0.98872169196032378, -0.98849079285269659, -0.98825756773074946,
    -0.98802201714328353, -0.98778414164457218, -0.98754394179435923,
    -0.98730141815785843, -0.987056571305751, -0.98680940181418553,
    -0.98655991026477541, -0.98630809724459867, -0.98605396334619544,
    -0.98579750916756748, -0.98553873531217606, -0.98527764238894122,
    -0.98501423101223984, -0.98474850180190421, -0.98448045538322093,
    -0.984210092386929, -0.98393741344921892, -0.98366241921173025,
    -0.98338511032155118, -0.98310548743121629, -0.98282355119870524,
    -0.98253930228744124, -0.98225274136628937, -0.98196386910955524,
    -0.98167268619698311, -0.98137919331375456, -0.98108339115048671,
    -0.98078528040323043, -0.98048486177346938, -0.98018213596811743,
    -0.97987710369951764, -0.97956976568544052, -0.979260122649082,
    -0.9789481753190622, -0.97863392442942321, -0.97831737071962765,
    -0.97799851493455714, -0.97767735782450993, -0.97735390014520007,
    -0.97702814265775439, -0.97670008612871184, -0.97636973133002114,
    -0.976037079039039, -0.97570213003852857, -0.975364885116657,
    -0.97502534506699412, -0.97468351068851067, -0.97433938278557586,
    -0.97399296216795583, -0.973644249650812, -0.97329324605469825,
    -0.97293995220556018, -0.97258436893473221, -0.97222649707893627,
    -0.9718663374802794, -0.97150389098625178, -0.97113915844972509,
    -0.97077214072895035, -0.9704028386875555, -0.970031253194544,
    -0.96965738512429245, -0.96928123535654853, -0.96890280477642887,
    -0.96852209427441727, -0.96813910474636244, -0.96775383709347551,
    -0.96736629222232851, -0.96697647104485207, -0.96658437447833312,
    -0.9661900034454125, -0.96579335887408368, -0.9653944416976894,
    -0.96499325285492032, -0.96458979328981276, -0.96418406395174583,
    -0.96377606579543984, -0.963365799780954, -0.96295326687368388,
    -0.96253846804435916, -0.96212140426904158, -0.96170207652912254,
    -0.96128048581132064, -0.96085663310767966, -0.96043051941556579,
    -0.960002145737666, -0.95957151308198452, -0.95913862246184189,
    -0.9587034748958716, -0.95826607140801767, -0.95782641302753291,
    -0.95738450078897586, -0.95694033573220882, -0.9564939189023951,
    -0.95604525134999641, -0.95559433413077111, -0.95514116830577078,
    -0.95468575494133834, -0.95422809510910567, -0.95376818988599033,
    -0.95330604035419386, -0.95284164760119872, -0.95237501271976588,
    -0.95190613680793235, -0.95143502096900834, -0.95096166631157508,
    -0.9504860739494817, -0.950008245001843, -0.94952818059303667,
    -0.94904588185270056, -0.94856134991573027, -0.94807458592227623,
    -0.94758559101774109, -0.94709436635277722, -0.94660091308328353,
    -0.94610523237040345, -0.94560732538052128, -0.94510719328526061,
    -0.94460483726148026, -0.94410025849127266, -0.94359345816196039,
    -0.94308443746609349, -0.94257319760144687, -0.94205973977101731,
    -0.94154406518302081, -0.94102617505088926, -0.9405060705932683,
    -0.93998375303401394, -0.93945922360218992, -0.9389324835320646,
    -0.93840353406310806, -0.93787237643998989, -0.937339011912575,
    -0.93680344173592156, -0.93626566717027826, -0.93572568948108037,
    -0.93518350993894761, -0.93463912981968078, -0.93409255040425887,
    -0.93354377297883617, -0.932992798834739, -0.93243962926846236,
    -0.93188426558166815, -0.93132670908118043, -0.93076696107898371,
    -0.93020502289221907, -0.92964089584318121, -0.92907458125931586,
    -0.92850608047321559, -0.92793539482261789, -0.92736252565040111,
    -0.92678747430458175, -0.92621024213831138, -0.92563083050987272,
    -0.92504924078267758, -0.9244654743252626, -0.92387953251128674,
    -0.92329141671952764, -0.92270112833387863, -0.92210866874334518,
    -0.9215140393420419, -0.92091724152918952, -0.92031827670911059,
    -0.91971714629122736, -0.91911385169005777, -0.91850839432521225,
    -0.9179007756213905, -0.91729099700837791, -0.9166790599210427,
    -0.91606496579933172, -0.91544871608826783, -0.9148303122379462,
    -0.91420975570353069, -0.91358704794525081, -0.91296219042839821,
    -0.91233518462332275, -0.91170603200542988, -0.91107473405517636,
    -0.91044129225806725, -0.90980570810465222, -0.90916798309052238,
    -0.90852811871630612, -0.90788611648766626, -0.90724197791529582,
    -0.90659570451491533, -0.90594729780726846, -0.90529675931811882,
    -0.90464409057824624, -0.90398929312344334, -0.90333236849451182,
    -0.90267331823725883, -0.90201214390249318, -0.901348847046022,
    -0.900683429228647, -0.90001589201616017, -0.89934623697934157,
    -0.89867446569395382, -0.89800057974073988, -0.89732458070541832,
    -0.89664647017868015, -0.89596624975618522, -0.89528392103855758,
    -0.8945994856313827, -0.89391294514520325, -0.89322430119551532,
    -0.89253355540276458, -0.89184070939234272, -0.89114576479458318,
    -0.89044872324475788, -0.88974958638307278, -0.88904835585466457,
    -0.88834503330959635, -0.88763962040285393, -0.88693211879434219,
    -0.88622253014888064, -0.8855108561362, -0.88479709843093779,
    -0.884081258712635, -0.88336333866573158, -0.88264333997956279,
    -0.881921264348355, -0.88119711347122209, -0.88047088905216075,
    -0.87974259280004741, -0.87901222642863353, -0.87827979165654158,
    -0.87754529020726135, -0.87680872380914565, -0.8760700941954066,
    -0.87532940310411089, -0.87458665227817611, -0.87384184346536686,
    -0.87309497841829009, -0.87234605889439154, -0.87159508665595109,
    -0.870842063470079, -0.87008699110871146, -0.86932987134860684,
    -0.8685707059713409, -0.86780949676330332, -0.86704624551569265,
    -0.866280954024513, -0.86551362409056909, -0.86474425751946238,
    -0.8639728561215867, -0.86319942171212416, -0.8624239561110405,
    -0.8616464611430813, -0.86086693863776731, -0.86008539042939014,
    -0.85930181835700847, -0.85851622426444274, -0.85772861000027212,
    -0.85693897741782876, -0.85614732837519447, -0.855353664735196,
    -0.85455798836540053, -0.85376030113811141, -0.85296060493036363,
    -0.85215890162391983, -0.8513551931052652, -0.85054948126560348,
    -0.84974176800085255, -0.84893205521163961, -0.84812034480329723,
    -0.84730663868585832, -0.84649093877405213, -0.84567324698729907,
    -0.84485356524970712, -0.84403189549006641, -0.84320823964184544,
    -0.84238259964318585, -0.84155497743689844, -0.84072537497045807,
    -0.83989379419599952, -0.83906023707031274, -0.83822470555483808,
    -0.83738720161566194, -0.836547727223512, -0.8357062843537526,
    -0.83486287498638, -0.83401750110601813, -0.83317016470191319,
    -0.83232086776792968, -0.83146961230254524, -0.83061640030884631,
    -0.829761233794523, -0.82890411477186487, -0.8280450452577558,
    -0.82718402727366913, -0.82632106284566353, -0.82545615400437755,
    -0.82458930278502529, -0.82372051122739143, -0.82284978137582643,
    -0.82197711527924155, -0.82110251499110465, -0.82022598256943469,
    -0.819347520076797, -0.81846712958029866, -0.81758481315158371,
    -0.81670057286682785, -0.81581441080673378, -0.81492632905652662,
    -0.81403632970594841, -0.81314441484925359, -0.81225058658520388,
    -0.81135484701706373, -0.81045719825259477, -0.80955764240405126,
    -0.808656181588175, -0.80775281792619036, -0.80684755354379933,
    -0.80594039057117628, -0.80503133114296366, -0.8041203773982657,
    -0.80320753148064494, -0.80229279553811572, -0.80137617172314024,
    -0.80045766219262282, -0.799537269107905, -0.79861499463476093,
    -0.79769084094339116, -0.79676481020841883, -0.79583690460888357,
    -0.794907126328237, -0.79397547755433717, -0.79304196047944364,
    -0.79210657730021239, -0.7911693302176902, -0.79023022143731,
    -0.78928925316888565, -0.78834642762660634, -0.78740174702903143,
    -0.78645521359908577, -0.78550682956405393, -0.78455659715557524,
    -0.78360451860963831, -0.78265059616657573, -0.78169483207105939,
    -0.78073722857209449, -0.77977778792301455, -0.778816512381476,
    -0.77785340420945315, -0.77688846567323244, -0.77592169904340769,
    -0.77495310659487393, -0.7739826906068229, -0.773010453362737,
    -0.77203639715038452, -0.77106052426181382, -0.7700828369933479,
    -0.7691033376455797, -0.76812202852336542, -0.7671389119358204,
    -0.76615399019631292, -0.765167265622459, -0.76417874053611679,
    -0.76318841726338127, -0.7621962981345789, -0.76120238548426178,
    -0.76020668165120242, -0.759209188978388, -0.75820990981301528,
    -0.75720884650648457, -0.75620600141439454, -0.75520137689653655,
    -0.75419497531688917, -0.75318679904361252, -0.7521768504490427,
    -0.75116513190968637, -0.75015164580621507, -0.74913639452345937,
    -0.7481193804504036, -0.74710060598018013, -0.74608007351006378,
    -0.745057785441466, -0.74403374417992929, -0.74300795213512172,
    -0.74198041172083107, -0.74095112535495922, -0.7399200954595162,
    -0.73888732446061511, -0.737852814788466, -0.73681656887736979,
    -0.73577858916571359, -0.7347388780959635, -0.73369743811466037,
    -0.73265427167241282, -0.73160938122389263, -0.73056276922782759,
    -0.729514438146997, -0.7284643904482252, -0.72741262860237577,
    -0.726359155084346, -0.72530397237306077, -0.724247082951467,
    -0.72318848930652746, -0.72212819392921535, -0.72106619931450811,
    -0.72000250796138165, -0.71893712237280449, -0.71787004505573171,
    -0.71680127852109954, -0.71573082528381859, -0.71465868786276909,
    -0.71358486878079352, -0.71250937056469243, -0.71143219574521643,
    -0.71035334685706242, -0.70927282643886569, -0.7081906370331954,
    -0.70710678118654757, -0.70602126144933974, -0.70493408037590488,
    -0.70384524052448494, -0.7027547444572253, -0.70166259474016845,
    -0.70056879394324834, -0.69947334464028377, -0.69837624940897292,
    -0.69727751083088652, -0.696177131491463, -0.69507511398000088,
    -0.69397146088965389, -0.69286617481742463, -0.69175925836415775,
    -0.6906507141345346, -0.68954054473706683, -0.68842875278409044,
    -0.687315340891759, -0.68620031168003859, -0.68508366777270036,
    -0.6839654117973154, -0.68284554638524808, -0.68172407417164971,
    -0.680600997795453, -0.679476319899365, -0.67835004312986147,
    -0.67722217013718033, -0.67609270357531592, -0.674961646102012,
    -0.673829000378756, -0.67269476907077286, -0.67155895484701833,
    -0.67042156038017309, -0.669282588346636, -0.66814204142651845,
    -0.66699992230363747, -0.66585623366550972, -0.66471097820334479,
    -0.66356415861203977, -0.66241577759017178, -0.66126583783999227,
    -0.66011434206742048, -0.65896129298203732, -0.65780669329707864,
    -0.65665054572942894, -0.65549285299961535, -0.65433361783180044,
    -0.65317284295377676, -0.6520105310969595, -0.650846684996381,
    -0.64968130739068319, -0.64851440102211244, -0.64734596863651206,
    -0.64617601298331628, -0.64500453681554393, -0.64383154288979139,
    -0.64265703396622686, -0.641481012808583, -0.64030348218415167,
    -0.63912444486377573, -0.637943903621844, -0.6367618612362842,
    -0.63557832048855611, -0.63439328416364549, -0.63320675505005719,
    -0.63201873593980906, -0.63082922962842447, -0.629638238914927,
    -0.6284457666018326, -0.62725181549514408, -0.62605638840434352,
    -0.62485948814238634, -0.62366111752569453, -0.62246127937415,
    -0.62125997651108755, -0.6200572117632891, -0.61885298796097632,
    -0.61764730793780387, -0.61644017453085365, -0.61523159058062682,
    -0.61402155893103849, -0.61281008242940971, -0.61159716392646191,
    -0.61038280627630948, -0.60916701233645321, -0.60794978496777363,
    -0.60673112703452448, -0.60551104140432555, -0.604289530948156,
    -0.60306659854034816, -0.60184224705858, -0.600616479383869,
    -0.59938929840056454, -0.59816070699634238, -0.59693070806219639,
    -0.59569930449243336, -0.59446649918466443, -0.5932322950397998,
    -0.591996694962041, -0.59075970185887416, -0.58952131864106394,
    -0.58828154822264522, -0.587040393520918, -0.58579785745643886,
    -0.58455394295301533, -0.58330865293769829, -0.58206199034077544,
    -0.58081395809576453, -0.57956455913940563, -0.57831379641165559,
    -0.57706167285567944, -0.57580819141784534, -0.57455335504771576,
    -0.5732971666980422, -0.572039629324757, -0.57078074588696726,
    -0.56952051934694714, -0.56825895267013149, -0.56699604882510868,
    -0.56573181078361312, -0.5644662415205195, -0.56319934401383409,
    -0.56193112124468936, -0.560661576197336, -0.55939071185913614,
    -0.5581185312205561, -0.5568450372751601, -0.55557023301960218,
    -0.55429412145362, -0.55301670558002747, -0.55173798840470734,
    -0.55045797293660481, -0.54917666218771966, -0.54789405917310019,
    -0.54661016691083486, -0.54532498842204646, -0.54403852673088382,
    -0.54275078486451589, -0.54146176585312344, -0.54017147272989285,
    -0.53887990853100842, -0.53758707629564539, -0.53629297906596318,
    -0.53499761988709715, -0.533701001807153, -0.5324031278771979,
    -0.531104001151255, -0.52980362468629461, -0.52850200154222848,
    -0.52719913478190128, -0.52589502747108463, -0.524589682678469,
    -0.52328310347565643, -0.52197529293715439, -0.52066625414036716,
    -0.51935599016558964, -0.51804450409599934, -0.51673179901764987,
    -0.51541787801946293, -0.51410274419322166, -0.512786400633563,
    -0.5114688504379703, -0.51015009670676681, -0.508830142543107,
    -0.50750899105297087, -0.50618664534515523, -0.50486310853126759,
    -0.50353838372571758, -0.50221247404571079, -0.50088538261124071,
    -0.49955711254508184, -0.49822766697278181, -0.49689704902265447,
    -0.49556526182577254, -0.49423230851595967, -0.49289819222978404,
    -0.4915629161065499, -0.49022648328829116, -0.48888889691976317,
    -0.487550160148436, -0.48621027612448642, -0.48486924800079106,
    -0.48352707893291874, -0.48218377207912272, -0.48083933060033396,
    -0.47949375766015295, -0.478147056424843, -0.47679923006332209,
    -0.47545028174715587, -0.47410021465054997, -0.47274903195034279,
    -0.47139673682599764, -0.47004333245959562, -0.46868882203582796,
    -0.46733320874198842, -0.46597649576796618, -0.46461868630623782,
    -0.46325978355186015, -0.46189979070246273, -0.46053871095824,
    -0.45917654752194409, -0.45781330359887717, -0.45644898239688392,
    -0.45508358712634384, -0.45371712100016387, -0.45234958723377089,
    -0.45098098904510386, -0.44961132965460654, -0.44824061228521989,
    -0.44686884016237416, -0.44549601651398174, -0.4441221445704292,
    -0.44274722756457, -0.44137126873171667, -0.43999427130963326,
    -0.43861623853852766, -0.43723717366104409, -0.43585707992225547,
    -0.43447596056965565, -0.43309381885315196, -0.43171065802505726,
    -0.43032648134008261, -0.42894129205532949, -0.42755509343028208,
    -0.42616788872679962, -0.42477968120910881, -0.42339047414379605,
    -0.42200027079979968, -0.42060907444840251, -0.41921688836322391,
    -0.41782371582021227, -0.41642956009763715, -0.41503442447608163,
    -0.4136383122384345, -0.41224122666988289, -0.41084317105790391,
    -0.40944414869225759, -0.40804416286497869, -0.40664321687036903,
    -0.40524131400498986, -0.40383845756765407, -0.40243465085941843,
    -0.40102989718357562, -0.39962419984564679, -0.39821756215337356,
    -0.39680998741671031, -0.39540147894781635, -0.3939920400610481,
    -0.39258167407295147, -0.39117038430225387, -0.38975817406985641,
    -0.38834504669882625, -0.38693100551438858, -0.38551605384391885,
    -0.38410019501693504, -0.38268343236508978, -0.38126576922216238,
    -0.37984720892405116, -0.37842775480876556, -0.37700741021641826,
    -0.37558617848921722, -0.37416406297145793, -0.37274106700951576,
    -0.37131719395183749, -0.3698924471489341, -0.36846682995337232,
    -0.36704034571976718, -0.36561299780477385, -0.36418478956707989,
    -0.36275572436739723, -0.36132580556845428, -0.35989503653498811,
    -0.35846342063373654, -0.35703096123343, -0.35559766170478385,
    -0.35416352542049034, -0.35272855575521073, -0.35129275608556709,
    -0.34985612979013492, -0.34841868024943456, -0.34698041084592368,
    -0.34554132496398909, -0.34410142598993881, -0.34266071731199438,
    -0.34121920232028236, -0.33977688440682685, -0.33833376696554113,
    -0.33688985339222005, -0.3354451470845316, -0.33399965144200938,
    -0.33255336986604422, -0.33110630575987643, -0.32965846252858749,
    -0.3282098435790925, -0.32676045232013173, -0.32531029216226293,
    -0.32385936651785285, -0.32240767880106985, -0.32095523242787521,
    -0.31950203081601569, -0.31804807738501495, -0.31659337555616585,
    -0.31513792875252244, -0.31368174039889152, -0.31222481392182488,
    -0.31076715274961147, -0.30930876031226873, -0.30784964004153487,
    -0.30638979537086092, -0.30492922973540237, -0.30346794657201132,
    -0.30200594931922808, -0.30054324141727345, -0.29907982630804048,
    -0.2976157074350862, -0.29615088824362379, -0.29468537218051433,
    -0.29321916269425863, -0.29175226323498926, -0.29028467725446233,
    -0.28881640820604948, -0.28734745954472951, -0.28587783472708062,
    -0.28440753721127188, -0.28293657045705539, -0.28146493792575794,
    -0.27999264308027322, -0.27851968938505306, -0.2770460803060999,
    -0.27557181931095814, -0.27409690986870638, -0.272621355449949,
    -0.271145159526808, -0.26966832557291509, -0.26819085706340318,
    -0.26671275747489837, -0.26523403028551179, -0.26375467897483135,
    -0.26227470702391359, -0.26079411791527551, -0.25931291513288623,
    -0.257831102162159, -0.25634868248994291, -0.25486565960451457,
    -0.25338203699557016, -0.25189781815421697, -0.25041300657296522,
    -0.24892760574572015, -0.24744161916777327, -0.24595505033579459,
    -0.24446790274782415, -0.24298017990326387, -0.24149188530286933,
    -0.2400030224487415, -0.23851359484431842, -0.2370236059943672,
    -0.23553305940497549, -0.23404195858354343, -0.23255030703877524,
    -0.23105810828067111, -0.22956536582051887, -0.22807208317088573,
    -0.22657826384561, -0.22508391135979283, -0.22358902922979,
    -0.22209362097320351, -0.22059769010887351, -0.2191012401568698,
    -0.21760427463848364, -0.21610679707621952, -0.21460881099378676,
    -0.21311031991609136, -0.21161132736922755, -0.21011183688046961,
    -0.20861185197826349, -0.20711137619221856, -0.20561041305309924,
    -0.20410896609281687, -0.20260703884442113, -0.2011046348420919,
    -0.19960175762113097, -0.19809841071795356, -0.19659459767008022,
    -0.19509032201612825, -0.19358558729580361, -0.19208039704989244,
    -0.19057475482025274, -0.18906866414980619, -0.1875621285825296,
    -0.18605515166344663, -0.18454773693861962, -0.18303988795514095,
    -0.18153160826112497, -0.18002290140569951, -0.17851377093899751,
    -0.17700422041214875, -0.17549425337727143, -0.17398387338746382,
    -0.17247308399679595, -0.17096188876030122, -0.16945029123396796,
    -0.16793829497473117, -0.1664259035404641, -0.16491312048996992,
    -0.16339994938297323, -0.16188639378011183, -0.16037245724292828,
    -0.15885814333386145, -0.15734345561623825, -0.15582839765426523,
    -0.1543129730130201, -0.15279718525844344, -0.15128103795733022,
    -0.14976453467732151, -0.14824767898689603, -0.14673047445536175,
    -0.14521292465284746, -0.14369503315029447, -0.14217680351944803,
    -0.14065823933284921, -0.1391393441638262, -0.13762012158648604,
    -0.1361005751757062, -0.13458070850712617, -0.13306052515713906,
    -0.13154002870288312, -0.13001922272223335, -0.12849811079379317,
    -0.12697669649688587, -0.12545498341154623, -0.12393297511851216,
    -0.1224106751992162, -0.12088808723577708, -0.11936521481099135,
    -0.11784206150832498, -0.11631863091190475, -0.11479492660651008,
    -0.11327095217756435, -0.11174671121112659, -0.11022220729388306,
    -0.10869744401313872, -0.10717242495680884, -0.10564715371341062,
    -0.10412163387205459, -0.10259586902243628, -0.10106986275482782,
    -0.099543618660069319, -0.0980171403295606, -0.096490431355252593,
    -0.094963495329638992, -0.093436335845747787, -0.091908956497132724,
    -0.090381360877864983, -0.0888535525825246, -0.087325535206192059,
    -0.0857973123444399, -0.084268887593324071, -0.082740264549375692,
    -0.081211446809592441, -0.079682437971430126, -0.078153241632794232,
    -0.076623861392031492, -0.0750943008479213, -0.073564563599667426,
    -0.072034653246889332, -0.070504573389613856, -0.068974327628266746,
    -0.067443919563664051, -0.0659133527970038, -0.064382630929857465,
    -0.0628517575641614, -0.061320736302208578, -0.059789570746639868,
    -0.058258264500435752, -0.056726821166907748, -0.055195244349689941,
    -0.05366353765273052, -0.052131704680283324, -0.050599749036899282,
    -0.049067674327418015, -0.0475354841569593, -0.046003182130914623,
    -0.044470771854938668, -0.04293825693494082, -0.041405640977076739,
    -0.039872927587739811, -0.038340120373552694, -0.036807222941358832,
    -0.035274238898213947, -0.03374117185137758, -0.032208025408304586,
    -0.030674803176636626, -0.029141508764193722, -0.02760814577896574,
    -0.0260747178291039, -0.024541228522912288, -0.023007681468839369,
    -0.021474080275469508, -0.019940428551514441, -0.01840672990580482,
    -0.01687298794728171, -0.0153392062849881, -0.013805388528060391,
    -0.012271538285719925, -0.010737659167264491, -0.00920375478205982,
    -0.007669828739531097, -0.0061358846491544753, -0.0046019261204485705,
    -0.0030679567629659761, -0.0015339801862847655, -0.0 };

  static creal_T fy[4096];
  creal_T fv[4096];
  static const double sintabinv[2049] = { 0.0, 0.0015339801862847655,
    0.0030679567629659761, 0.0046019261204485705, 0.0061358846491544753,
    0.007669828739531097, 0.00920375478205982, 0.010737659167264491,
    0.012271538285719925, 0.013805388528060391, 0.0153392062849881,
    0.01687298794728171, 0.01840672990580482, 0.019940428551514441,
    0.021474080275469508, 0.023007681468839369, 0.024541228522912288,
    0.0260747178291039, 0.02760814577896574, 0.029141508764193722,
    0.030674803176636626, 0.032208025408304586, 0.03374117185137758,
    0.035274238898213947, 0.036807222941358832, 0.038340120373552694,
    0.039872927587739811, 0.041405640977076739, 0.04293825693494082,
    0.044470771854938668, 0.046003182130914623, 0.0475354841569593,
    0.049067674327418015, 0.050599749036899282, 0.052131704680283324,
    0.05366353765273052, 0.055195244349689941, 0.056726821166907748,
    0.058258264500435752, 0.059789570746639868, 0.061320736302208578,
    0.0628517575641614, 0.064382630929857465, 0.0659133527970038,
    0.067443919563664051, 0.068974327628266746, 0.070504573389613856,
    0.072034653246889332, 0.073564563599667426, 0.0750943008479213,
    0.076623861392031492, 0.078153241632794232, 0.079682437971430126,
    0.081211446809592441, 0.082740264549375692, 0.084268887593324071,
    0.0857973123444399, 0.087325535206192059, 0.0888535525825246,
    0.090381360877864983, 0.091908956497132724, 0.093436335845747787,
    0.094963495329638992, 0.096490431355252593, 0.0980171403295606,
    0.099543618660069319, 0.10106986275482782, 0.10259586902243628,
    0.10412163387205459, 0.10564715371341062, 0.10717242495680884,
    0.10869744401313872, 0.11022220729388306, 0.11174671121112659,
    0.11327095217756435, 0.11479492660651008, 0.11631863091190475,
    0.11784206150832498, 0.11936521481099135, 0.12088808723577708,
    0.1224106751992162, 0.12393297511851216, 0.12545498341154623,
    0.12697669649688587, 0.12849811079379317, 0.13001922272223335,
    0.13154002870288312, 0.13306052515713906, 0.13458070850712617,
    0.1361005751757062, 0.13762012158648604, 0.1391393441638262,
    0.14065823933284921, 0.14217680351944803, 0.14369503315029447,
    0.14521292465284746, 0.14673047445536175, 0.14824767898689603,
    0.14976453467732151, 0.15128103795733022, 0.15279718525844344,
    0.1543129730130201, 0.15582839765426523, 0.15734345561623825,
    0.15885814333386145, 0.16037245724292828, 0.16188639378011183,
    0.16339994938297323, 0.16491312048996992, 0.1664259035404641,
    0.16793829497473117, 0.16945029123396796, 0.17096188876030122,
    0.17247308399679595, 0.17398387338746382, 0.17549425337727143,
    0.17700422041214875, 0.17851377093899751, 0.18002290140569951,
    0.18153160826112497, 0.18303988795514095, 0.18454773693861962,
    0.18605515166344663, 0.1875621285825296, 0.18906866414980619,
    0.19057475482025274, 0.19208039704989244, 0.19358558729580361,
    0.19509032201612825, 0.19659459767008022, 0.19809841071795356,
    0.19960175762113097, 0.2011046348420919, 0.20260703884442113,
    0.20410896609281687, 0.20561041305309924, 0.20711137619221856,
    0.20861185197826349, 0.21011183688046961, 0.21161132736922755,
    0.21311031991609136, 0.21460881099378676, 0.21610679707621952,
    0.21760427463848364, 0.2191012401568698, 0.22059769010887351,
    0.22209362097320351, 0.22358902922979, 0.22508391135979283, 0.22657826384561,
    0.22807208317088573, 0.22956536582051887, 0.23105810828067111,
    0.23255030703877524, 0.23404195858354343, 0.23553305940497549,
    0.2370236059943672, 0.23851359484431842, 0.2400030224487415,
    0.24149188530286933, 0.24298017990326387, 0.24446790274782415,
    0.24595505033579459, 0.24744161916777327, 0.24892760574572015,
    0.25041300657296522, 0.25189781815421697, 0.25338203699557016,
    0.25486565960451457, 0.25634868248994291, 0.257831102162159,
    0.25931291513288623, 0.26079411791527551, 0.26227470702391359,
    0.26375467897483135, 0.26523403028551179, 0.26671275747489837,
    0.26819085706340318, 0.26966832557291509, 0.271145159526808,
    0.272621355449949, 0.27409690986870638, 0.27557181931095814,
    0.2770460803060999, 0.27851968938505306, 0.27999264308027322,
    0.28146493792575794, 0.28293657045705539, 0.28440753721127188,
    0.28587783472708062, 0.28734745954472951, 0.28881640820604948,
    0.29028467725446233, 0.29175226323498926, 0.29321916269425863,
    0.29468537218051433, 0.29615088824362379, 0.2976157074350862,
    0.29907982630804048, 0.30054324141727345, 0.30200594931922808,
    0.30346794657201132, 0.30492922973540237, 0.30638979537086092,
    0.30784964004153487, 0.30930876031226873, 0.31076715274961147,
    0.31222481392182488, 0.31368174039889152, 0.31513792875252244,
    0.31659337555616585, 0.31804807738501495, 0.31950203081601569,
    0.32095523242787521, 0.32240767880106985, 0.32385936651785285,
    0.32531029216226293, 0.32676045232013173, 0.3282098435790925,
    0.32965846252858749, 0.33110630575987643, 0.33255336986604422,
    0.33399965144200938, 0.3354451470845316, 0.33688985339222005,
    0.33833376696554113, 0.33977688440682685, 0.34121920232028236,
    0.34266071731199438, 0.34410142598993881, 0.34554132496398909,
    0.34698041084592368, 0.34841868024943456, 0.34985612979013492,
    0.35129275608556709, 0.35272855575521073, 0.35416352542049034,
    0.35559766170478385, 0.35703096123343, 0.35846342063373654,
    0.35989503653498811, 0.36132580556845428, 0.36275572436739723,
    0.36418478956707989, 0.36561299780477385, 0.36704034571976718,
    0.36846682995337232, 0.3698924471489341, 0.37131719395183749,
    0.37274106700951576, 0.37416406297145793, 0.37558617848921722,
    0.37700741021641826, 0.37842775480876556, 0.37984720892405116,
    0.38126576922216238, 0.38268343236508978, 0.38410019501693504,
    0.38551605384391885, 0.38693100551438858, 0.38834504669882625,
    0.38975817406985641, 0.39117038430225387, 0.39258167407295147,
    0.3939920400610481, 0.39540147894781635, 0.39680998741671031,
    0.39821756215337356, 0.39962419984564679, 0.40102989718357562,
    0.40243465085941843, 0.40383845756765407, 0.40524131400498986,
    0.40664321687036903, 0.40804416286497869, 0.40944414869225759,
    0.41084317105790391, 0.41224122666988289, 0.4136383122384345,
    0.41503442447608163, 0.41642956009763715, 0.41782371582021227,
    0.41921688836322391, 0.42060907444840251, 0.42200027079979968,
    0.42339047414379605, 0.42477968120910881, 0.42616788872679962,
    0.42755509343028208, 0.42894129205532949, 0.43032648134008261,
    0.43171065802505726, 0.43309381885315196, 0.43447596056965565,
    0.43585707992225547, 0.43723717366104409, 0.43861623853852766,
    0.43999427130963326, 0.44137126873171667, 0.44274722756457,
    0.4441221445704292, 0.44549601651398174, 0.44686884016237416,
    0.44824061228521989, 0.44961132965460654, 0.45098098904510386,
    0.45234958723377089, 0.45371712100016387, 0.45508358712634384,
    0.45644898239688392, 0.45781330359887717, 0.45917654752194409,
    0.46053871095824, 0.46189979070246273, 0.46325978355186015,
    0.46461868630623782, 0.46597649576796618, 0.46733320874198842,
    0.46868882203582796, 0.47004333245959562, 0.47139673682599764,
    0.47274903195034279, 0.47410021465054997, 0.47545028174715587,
    0.47679923006332209, 0.478147056424843, 0.47949375766015295,
    0.48083933060033396, 0.48218377207912272, 0.48352707893291874,
    0.48486924800079106, 0.48621027612448642, 0.487550160148436,
    0.48888889691976317, 0.49022648328829116, 0.4915629161065499,
    0.49289819222978404, 0.49423230851595967, 0.49556526182577254,
    0.49689704902265447, 0.49822766697278181, 0.49955711254508184,
    0.50088538261124071, 0.50221247404571079, 0.50353838372571758,
    0.50486310853126759, 0.50618664534515523, 0.50750899105297087,
    0.508830142543107, 0.51015009670676681, 0.5114688504379703,
    0.512786400633563, 0.51410274419322166, 0.51541787801946293,
    0.51673179901764987, 0.51804450409599934, 0.51935599016558964,
    0.52066625414036716, 0.52197529293715439, 0.52328310347565643,
    0.524589682678469, 0.52589502747108463, 0.52719913478190128,
    0.52850200154222848, 0.52980362468629461, 0.531104001151255,
    0.5324031278771979, 0.533701001807153, 0.53499761988709715,
    0.53629297906596318, 0.53758707629564539, 0.53887990853100842,
    0.54017147272989285, 0.54146176585312344, 0.54275078486451589,
    0.54403852673088382, 0.54532498842204646, 0.54661016691083486,
    0.54789405917310019, 0.54917666218771966, 0.55045797293660481,
    0.55173798840470734, 0.55301670558002747, 0.55429412145362,
    0.55557023301960218, 0.5568450372751601, 0.5581185312205561,
    0.55939071185913614, 0.560661576197336, 0.56193112124468936,
    0.56319934401383409, 0.5644662415205195, 0.56573181078361312,
    0.56699604882510868, 0.56825895267013149, 0.56952051934694714,
    0.57078074588696726, 0.572039629324757, 0.5732971666980422,
    0.57455335504771576, 0.57580819141784534, 0.57706167285567944,
    0.57831379641165559, 0.57956455913940563, 0.58081395809576453,
    0.58206199034077544, 0.58330865293769829, 0.58455394295301533,
    0.58579785745643886, 0.587040393520918, 0.58828154822264522,
    0.58952131864106394, 0.59075970185887416, 0.591996694962041,
    0.5932322950397998, 0.59446649918466443, 0.59569930449243336,
    0.59693070806219639, 0.59816070699634238, 0.59938929840056454,
    0.600616479383869, 0.60184224705858, 0.60306659854034816, 0.604289530948156,
    0.60551104140432555, 0.60673112703452448, 0.60794978496777363,
    0.60916701233645321, 0.61038280627630948, 0.61159716392646191,
    0.61281008242940971, 0.61402155893103849, 0.61523159058062682,
    0.61644017453085365, 0.61764730793780387, 0.61885298796097632,
    0.6200572117632891, 0.62125997651108755, 0.62246127937415,
    0.62366111752569453, 0.62485948814238634, 0.62605638840434352,
    0.62725181549514408, 0.6284457666018326, 0.629638238914927,
    0.63082922962842447, 0.63201873593980906, 0.63320675505005719,
    0.63439328416364549, 0.63557832048855611, 0.6367618612362842,
    0.637943903621844, 0.63912444486377573, 0.64030348218415167,
    0.641481012808583, 0.64265703396622686, 0.64383154288979139,
    0.64500453681554393, 0.64617601298331628, 0.64734596863651206,
    0.64851440102211244, 0.64968130739068319, 0.650846684996381,
    0.6520105310969595, 0.65317284295377676, 0.65433361783180044,
    0.65549285299961535, 0.65665054572942894, 0.65780669329707864,
    0.65896129298203732, 0.66011434206742048, 0.66126583783999227,
    0.66241577759017178, 0.66356415861203977, 0.66471097820334479,
    0.66585623366550972, 0.66699992230363747, 0.66814204142651845,
    0.669282588346636, 0.67042156038017309, 0.67155895484701833,
    0.67269476907077286, 0.673829000378756, 0.674961646102012,
    0.67609270357531592, 0.67722217013718033, 0.67835004312986147,
    0.679476319899365, 0.680600997795453, 0.68172407417164971,
    0.68284554638524808, 0.6839654117973154, 0.68508366777270036,
    0.68620031168003859, 0.687315340891759, 0.68842875278409044,
    0.68954054473706683, 0.6906507141345346, 0.69175925836415775,
    0.69286617481742463, 0.69397146088965389, 0.69507511398000088,
    0.696177131491463, 0.69727751083088652, 0.69837624940897292,
    0.69947334464028377, 0.70056879394324834, 0.70166259474016845,
    0.7027547444572253, 0.70384524052448494, 0.70493408037590488,
    0.70602126144933974, 0.70710678118654757, 0.7081906370331954,
    0.70927282643886569, 0.71035334685706242, 0.71143219574521643,
    0.71250937056469243, 0.71358486878079352, 0.71465868786276909,
    0.71573082528381859, 0.71680127852109954, 0.71787004505573171,
    0.71893712237280449, 0.72000250796138165, 0.72106619931450811,
    0.72212819392921535, 0.72318848930652746, 0.724247082951467,
    0.72530397237306077, 0.726359155084346, 0.72741262860237577,
    0.7284643904482252, 0.729514438146997, 0.73056276922782759,
    0.73160938122389263, 0.73265427167241282, 0.73369743811466037,
    0.7347388780959635, 0.73577858916571359, 0.73681656887736979,
    0.737852814788466, 0.73888732446061511, 0.7399200954595162,
    0.74095112535495922, 0.74198041172083107, 0.74300795213512172,
    0.74403374417992929, 0.745057785441466, 0.74608007351006378,
    0.74710060598018013, 0.7481193804504036, 0.74913639452345937,
    0.75015164580621507, 0.75116513190968637, 0.7521768504490427,
    0.75318679904361252, 0.75419497531688917, 0.75520137689653655,
    0.75620600141439454, 0.75720884650648457, 0.75820990981301528,
    0.759209188978388, 0.76020668165120242, 0.76120238548426178,
    0.7621962981345789, 0.76318841726338127, 0.76417874053611679,
    0.765167265622459, 0.76615399019631292, 0.7671389119358204,
    0.76812202852336542, 0.7691033376455797, 0.7700828369933479,
    0.77106052426181382, 0.77203639715038452, 0.773010453362737,
    0.7739826906068229, 0.77495310659487393, 0.77592169904340769,
    0.77688846567323244, 0.77785340420945315, 0.778816512381476,
    0.77977778792301455, 0.78073722857209449, 0.78169483207105939,
    0.78265059616657573, 0.78360451860963831, 0.78455659715557524,
    0.78550682956405393, 0.78645521359908577, 0.78740174702903143,
    0.78834642762660634, 0.78928925316888565, 0.79023022143731,
    0.7911693302176902, 0.79210657730021239, 0.79304196047944364,
    0.79397547755433717, 0.794907126328237, 0.79583690460888357,
    0.79676481020841883, 0.79769084094339116, 0.79861499463476093,
    0.799537269107905, 0.80045766219262282, 0.80137617172314024,
    0.80229279553811572, 0.80320753148064494, 0.8041203773982657,
    0.80503133114296366, 0.80594039057117628, 0.80684755354379933,
    0.80775281792619036, 0.808656181588175, 0.80955764240405126,
    0.81045719825259477, 0.81135484701706373, 0.81225058658520388,
    0.81314441484925359, 0.81403632970594841, 0.81492632905652662,
    0.81581441080673378, 0.81670057286682785, 0.81758481315158371,
    0.81846712958029866, 0.819347520076797, 0.82022598256943469,
    0.82110251499110465, 0.82197711527924155, 0.82284978137582643,
    0.82372051122739143, 0.82458930278502529, 0.82545615400437755,
    0.82632106284566353, 0.82718402727366913, 0.8280450452577558,
    0.82890411477186487, 0.829761233794523, 0.83061640030884631,
    0.83146961230254524, 0.83232086776792968, 0.83317016470191319,
    0.83401750110601813, 0.83486287498638, 0.8357062843537526, 0.836547727223512,
    0.83738720161566194, 0.83822470555483808, 0.83906023707031274,
    0.83989379419599952, 0.84072537497045807, 0.84155497743689844,
    0.84238259964318585, 0.84320823964184544, 0.84403189549006641,
    0.84485356524970712, 0.84567324698729907, 0.84649093877405213,
    0.84730663868585832, 0.84812034480329723, 0.84893205521163961,
    0.84974176800085255, 0.85054948126560348, 0.8513551931052652,
    0.85215890162391983, 0.85296060493036363, 0.85376030113811141,
    0.85455798836540053, 0.855353664735196, 0.85614732837519447,
    0.85693897741782876, 0.85772861000027212, 0.85851622426444274,
    0.85930181835700847, 0.86008539042939014, 0.86086693863776731,
    0.8616464611430813, 0.8624239561110405, 0.86319942171212416,
    0.8639728561215867, 0.86474425751946238, 0.86551362409056909,
    0.866280954024513, 0.86704624551569265, 0.86780949676330332,
    0.8685707059713409, 0.86932987134860684, 0.87008699110871146,
    0.870842063470079, 0.87159508665595109, 0.87234605889439154,
    0.87309497841829009, 0.87384184346536686, 0.87458665227817611,
    0.87532940310411089, 0.8760700941954066, 0.87680872380914565,
    0.87754529020726135, 0.87827979165654158, 0.87901222642863353,
    0.87974259280004741, 0.88047088905216075, 0.88119711347122209,
    0.881921264348355, 0.88264333997956279, 0.88336333866573158,
    0.884081258712635, 0.88479709843093779, 0.8855108561362, 0.88622253014888064,
    0.88693211879434219, 0.88763962040285393, 0.88834503330959635,
    0.88904835585466457, 0.88974958638307278, 0.89044872324475788,
    0.89114576479458318, 0.89184070939234272, 0.89253355540276458,
    0.89322430119551532, 0.89391294514520325, 0.8945994856313827,
    0.89528392103855758, 0.89596624975618522, 0.89664647017868015,
    0.89732458070541832, 0.89800057974073988, 0.89867446569395382,
    0.89934623697934157, 0.90001589201616017, 0.900683429228647,
    0.901348847046022, 0.90201214390249318, 0.90267331823725883,
    0.90333236849451182, 0.90398929312344334, 0.90464409057824624,
    0.90529675931811882, 0.90594729780726846, 0.90659570451491533,
    0.90724197791529582, 0.90788611648766626, 0.90852811871630612,
    0.90916798309052238, 0.90980570810465222, 0.91044129225806725,
    0.91107473405517636, 0.91170603200542988, 0.91233518462332275,
    0.91296219042839821, 0.91358704794525081, 0.91420975570353069,
    0.9148303122379462, 0.91544871608826783, 0.91606496579933172,
    0.9166790599210427, 0.91729099700837791, 0.9179007756213905,
    0.91850839432521225, 0.91911385169005777, 0.91971714629122736,
    0.92031827670911059, 0.92091724152918952, 0.9215140393420419,
    0.92210866874334518, 0.92270112833387863, 0.92329141671952764,
    0.92387953251128674, 0.9244654743252626, 0.92504924078267758,
    0.92563083050987272, 0.92621024213831138, 0.92678747430458175,
    0.92736252565040111, 0.92793539482261789, 0.92850608047321559,
    0.92907458125931586, 0.92964089584318121, 0.93020502289221907,
    0.93076696107898371, 0.93132670908118043, 0.93188426558166815,
    0.93243962926846236, 0.932992798834739, 0.93354377297883617,
    0.93409255040425887, 0.93463912981968078, 0.93518350993894761,
    0.93572568948108037, 0.93626566717027826, 0.93680344173592156,
    0.937339011912575, 0.93787237643998989, 0.93840353406310806,
    0.9389324835320646, 0.93945922360218992, 0.93998375303401394,
    0.9405060705932683, 0.94102617505088926, 0.94154406518302081,
    0.94205973977101731, 0.94257319760144687, 0.94308443746609349,
    0.94359345816196039, 0.94410025849127266, 0.94460483726148026,
    0.94510719328526061, 0.94560732538052128, 0.94610523237040345,
    0.94660091308328353, 0.94709436635277722, 0.94758559101774109,
    0.94807458592227623, 0.94856134991573027, 0.94904588185270056,
    0.94952818059303667, 0.950008245001843, 0.9504860739494817,
    0.95096166631157508, 0.95143502096900834, 0.95190613680793235,
    0.95237501271976588, 0.95284164760119872, 0.95330604035419386,
    0.95376818988599033, 0.95422809510910567, 0.95468575494133834,
    0.95514116830577078, 0.95559433413077111, 0.95604525134999641,
    0.9564939189023951, 0.95694033573220882, 0.95738450078897586,
    0.95782641302753291, 0.95826607140801767, 0.9587034748958716,
    0.95913862246184189, 0.95957151308198452, 0.960002145737666,
    0.96043051941556579, 0.96085663310767966, 0.96128048581132064,
    0.96170207652912254, 0.96212140426904158, 0.96253846804435916,
    0.96295326687368388, 0.963365799780954, 0.96377606579543984,
    0.96418406395174583, 0.96458979328981276, 0.96499325285492032,
    0.9653944416976894, 0.96579335887408368, 0.9661900034454125,
    0.96658437447833312, 0.96697647104485207, 0.96736629222232851,
    0.96775383709347551, 0.96813910474636244, 0.96852209427441727,
    0.96890280477642887, 0.96928123535654853, 0.96965738512429245,
    0.970031253194544, 0.9704028386875555, 0.97077214072895035,
    0.97113915844972509, 0.97150389098625178, 0.9718663374802794,
    0.97222649707893627, 0.97258436893473221, 0.97293995220556018,
    0.97329324605469825, 0.973644249650812, 0.97399296216795583,
    0.97433938278557586, 0.97468351068851067, 0.97502534506699412,
    0.975364885116657, 0.97570213003852857, 0.976037079039039,
    0.97636973133002114, 0.97670008612871184, 0.97702814265775439,
    0.97735390014520007, 0.97767735782450993, 0.97799851493455714,
    0.97831737071962765, 0.97863392442942321, 0.9789481753190622,
    0.979260122649082, 0.97956976568544052, 0.97987710369951764,
    0.98018213596811743, 0.98048486177346938, 0.98078528040323043,
    0.98108339115048671, 0.98137919331375456, 0.98167268619698311,
    0.98196386910955524, 0.98225274136628937, 0.98253930228744124,
    0.98282355119870524, 0.98310548743121629, 0.98338511032155118,
    0.98366241921173025, 0.98393741344921892, 0.984210092386929,
    0.98448045538322093, 0.98474850180190421, 0.98501423101223984,
    0.98527764238894122, 0.98553873531217606, 0.98579750916756748,
    0.98605396334619544, 0.98630809724459867, 0.98655991026477541,
    0.98680940181418553, 0.987056571305751, 0.98730141815785843,
    0.98754394179435923, 0.98778414164457218, 0.98802201714328353,
    0.98825756773074946, 0.98849079285269659, 0.98872169196032378,
    0.988950264510303, 0.989176509964781, 0.98940042779138038,
    0.98962201746320089, 0.98984127845882053, 0.99005821026229712,
    0.99027281236316911, 0.99048508425645709, 0.99069502544266463,
    0.99090263542778, 0.99110791372327689, 0.99131085984611544,
    0.9915114733187439, 0.99170975366909953, 0.99190570043060933,
    0.9920993131421918, 0.99229059134825737, 0.99247953459871, 0.992666142448948,
    0.9928504144598651, 0.99303235019785141, 0.9932119492347945,
    0.99338921114808065, 0.9935641355205953, 0.9937367219407246,
    0.99390697000235606, 0.99407487930487937, 0.9942404494531879,
    0.9944036800576791, 0.99456457073425542, 0.9947231211043257,
    0.99487933079480562, 0.99503319943811863, 0.99518472667219693,
    0.99533391214048228, 0.99548075549192694, 0.99562525638099431,
    0.99576741446765982, 0.99590722941741172, 0.996044700901252,
    0.996179828595697, 0.996312612182778, 0.99644305135004263,
    0.99657114579055484, 0.99669689520289606, 0.99682029929116567,
    0.99694135776498216, 0.997060070339483, 0.99717643673532619,
    0.99729045667869021, 0.9974021299012753, 0.99751145614030345,
    0.99761843513851955, 0.99772306664419164, 0.99782535041111164,
    0.997925286198596, 0.99802287377148624, 0.99811811290014918,
    0.99821100336047819, 0.99830154493389289, 0.99838973740734016,
    0.99847558057329477, 0.99855907422975931, 0.99864021818026527,
    0.99871901223387294, 0.99879545620517241, 0.99886954991428356,
    0.99894129318685687, 0.99901068585407338, 0.99907772775264536,
    0.99914241872481691, 0.99920475861836389, 0.99926474728659442,
    0.99932238458834954, 0.99937767038800285, 0.99943060455546173,
    0.999481186966167, 0.99952941750109314, 0.99957529604674922,
    0.99961882249517864, 0.99965999674395922, 0.99969881869620425,
    0.99973528826056168, 0.99976940535121528, 0.99980116988788426,
    0.9998305817958234, 0.99985764100582386, 0.99988234745421256,
    0.9999047010828529, 0.9999247018391445, 0.99994234967602391,
    0.9999576445519639, 0.99997058643097414, 0.99998117528260111,
    0.9999894110819284, 0.99999529380957619, 0.99999882345170188, 1.0,
    0.99999882345170188, 0.99999529380957619, 0.9999894110819284,
    0.99998117528260111, 0.99997058643097414, 0.9999576445519639,
    0.99994234967602391, 0.9999247018391445, 0.9999047010828529,
    0.99988234745421256, 0.99985764100582386, 0.9998305817958234,
    0.99980116988788426, 0.99976940535121528, 0.99973528826056168,
    0.99969881869620425, 0.99965999674395922, 0.99961882249517864,
    0.99957529604674922, 0.99952941750109314, 0.999481186966167,
    0.99943060455546173, 0.99937767038800285, 0.99932238458834954,
    0.99926474728659442, 0.99920475861836389, 0.99914241872481691,
    0.99907772775264536, 0.99901068585407338, 0.99894129318685687,
    0.99886954991428356, 0.99879545620517241, 0.99871901223387294,
    0.99864021818026527, 0.99855907422975931, 0.99847558057329477,
    0.99838973740734016, 0.99830154493389289, 0.99821100336047819,
    0.99811811290014918, 0.99802287377148624, 0.997925286198596,
    0.99782535041111164, 0.99772306664419164, 0.99761843513851955,
    0.99751145614030345, 0.9974021299012753, 0.99729045667869021,
    0.99717643673532619, 0.997060070339483, 0.99694135776498216,
    0.99682029929116567, 0.99669689520289606, 0.99657114579055484,
    0.99644305135004263, 0.996312612182778, 0.996179828595697, 0.996044700901252,
    0.99590722941741172, 0.99576741446765982, 0.99562525638099431,
    0.99548075549192694, 0.99533391214048228, 0.99518472667219693,
    0.99503319943811863, 0.99487933079480562, 0.9947231211043257,
    0.99456457073425542, 0.9944036800576791, 0.9942404494531879,
    0.99407487930487937, 0.99390697000235606, 0.9937367219407246,
    0.9935641355205953, 0.99338921114808065, 0.9932119492347945,
    0.99303235019785141, 0.9928504144598651, 0.992666142448948, 0.99247953459871,
    0.99229059134825737, 0.9920993131421918, 0.99190570043060933,
    0.99170975366909953, 0.9915114733187439, 0.99131085984611544,
    0.99110791372327689, 0.99090263542778, 0.99069502544266463,
    0.99048508425645709, 0.99027281236316911, 0.99005821026229712,
    0.98984127845882053, 0.98962201746320089, 0.98940042779138038,
    0.989176509964781, 0.988950264510303, 0.98872169196032378,
    0.98849079285269659, 0.98825756773074946, 0.98802201714328353,
    0.98778414164457218, 0.98754394179435923, 0.98730141815785843,
    0.987056571305751, 0.98680940181418553, 0.98655991026477541,
    0.98630809724459867, 0.98605396334619544, 0.98579750916756748,
    0.98553873531217606, 0.98527764238894122, 0.98501423101223984,
    0.98474850180190421, 0.98448045538322093, 0.984210092386929,
    0.98393741344921892, 0.98366241921173025, 0.98338511032155118,
    0.98310548743121629, 0.98282355119870524, 0.98253930228744124,
    0.98225274136628937, 0.98196386910955524, 0.98167268619698311,
    0.98137919331375456, 0.98108339115048671, 0.98078528040323043,
    0.98048486177346938, 0.98018213596811743, 0.97987710369951764,
    0.97956976568544052, 0.979260122649082, 0.9789481753190622,
    0.97863392442942321, 0.97831737071962765, 0.97799851493455714,
    0.97767735782450993, 0.97735390014520007, 0.97702814265775439,
    0.97670008612871184, 0.97636973133002114, 0.976037079039039,
    0.97570213003852857, 0.975364885116657, 0.97502534506699412,
    0.97468351068851067, 0.97433938278557586, 0.97399296216795583,
    0.973644249650812, 0.97329324605469825, 0.97293995220556018,
    0.97258436893473221, 0.97222649707893627, 0.9718663374802794,
    0.97150389098625178, 0.97113915844972509, 0.97077214072895035,
    0.9704028386875555, 0.970031253194544, 0.96965738512429245,
    0.96928123535654853, 0.96890280477642887, 0.96852209427441727,
    0.96813910474636244, 0.96775383709347551, 0.96736629222232851,
    0.96697647104485207, 0.96658437447833312, 0.9661900034454125,
    0.96579335887408368, 0.9653944416976894, 0.96499325285492032,
    0.96458979328981276, 0.96418406395174583, 0.96377606579543984,
    0.963365799780954, 0.96295326687368388, 0.96253846804435916,
    0.96212140426904158, 0.96170207652912254, 0.96128048581132064,
    0.96085663310767966, 0.96043051941556579, 0.960002145737666,
    0.95957151308198452, 0.95913862246184189, 0.9587034748958716,
    0.95826607140801767, 0.95782641302753291, 0.95738450078897586,
    0.95694033573220882, 0.9564939189023951, 0.95604525134999641,
    0.95559433413077111, 0.95514116830577078, 0.95468575494133834,
    0.95422809510910567, 0.95376818988599033, 0.95330604035419386,
    0.95284164760119872, 0.95237501271976588, 0.95190613680793235,
    0.95143502096900834, 0.95096166631157508, 0.9504860739494817,
    0.950008245001843, 0.94952818059303667, 0.94904588185270056,
    0.94856134991573027, 0.94807458592227623, 0.94758559101774109,
    0.94709436635277722, 0.94660091308328353, 0.94610523237040345,
    0.94560732538052128, 0.94510719328526061, 0.94460483726148026,
    0.94410025849127266, 0.94359345816196039, 0.94308443746609349,
    0.94257319760144687, 0.94205973977101731, 0.94154406518302081,
    0.94102617505088926, 0.9405060705932683, 0.93998375303401394,
    0.93945922360218992, 0.9389324835320646, 0.93840353406310806,
    0.93787237643998989, 0.937339011912575, 0.93680344173592156,
    0.93626566717027826, 0.93572568948108037, 0.93518350993894761,
    0.93463912981968078, 0.93409255040425887, 0.93354377297883617,
    0.932992798834739, 0.93243962926846236, 0.93188426558166815,
    0.93132670908118043, 0.93076696107898371, 0.93020502289221907,
    0.92964089584318121, 0.92907458125931586, 0.92850608047321559,
    0.92793539482261789, 0.92736252565040111, 0.92678747430458175,
    0.92621024213831138, 0.92563083050987272, 0.92504924078267758,
    0.9244654743252626, 0.92387953251128674, 0.92329141671952764,
    0.92270112833387863, 0.92210866874334518, 0.9215140393420419,
    0.92091724152918952, 0.92031827670911059, 0.91971714629122736,
    0.91911385169005777, 0.91850839432521225, 0.9179007756213905,
    0.91729099700837791, 0.9166790599210427, 0.91606496579933172,
    0.91544871608826783, 0.9148303122379462, 0.91420975570353069,
    0.91358704794525081, 0.91296219042839821, 0.91233518462332275,
    0.91170603200542988, 0.91107473405517636, 0.91044129225806725,
    0.90980570810465222, 0.90916798309052238, 0.90852811871630612,
    0.90788611648766626, 0.90724197791529582, 0.90659570451491533,
    0.90594729780726846, 0.90529675931811882, 0.90464409057824624,
    0.90398929312344334, 0.90333236849451182, 0.90267331823725883,
    0.90201214390249318, 0.901348847046022, 0.900683429228647,
    0.90001589201616017, 0.89934623697934157, 0.89867446569395382,
    0.89800057974073988, 0.89732458070541832, 0.89664647017868015,
    0.89596624975618522, 0.89528392103855758, 0.8945994856313827,
    0.89391294514520325, 0.89322430119551532, 0.89253355540276458,
    0.89184070939234272, 0.89114576479458318, 0.89044872324475788,
    0.88974958638307278, 0.88904835585466457, 0.88834503330959635,
    0.88763962040285393, 0.88693211879434219, 0.88622253014888064,
    0.8855108561362, 0.88479709843093779, 0.884081258712635, 0.88336333866573158,
    0.88264333997956279, 0.881921264348355, 0.88119711347122209,
    0.88047088905216075, 0.87974259280004741, 0.87901222642863353,
    0.87827979165654158, 0.87754529020726135, 0.87680872380914565,
    0.8760700941954066, 0.87532940310411089, 0.87458665227817611,
    0.87384184346536686, 0.87309497841829009, 0.87234605889439154,
    0.87159508665595109, 0.870842063470079, 0.87008699110871146,
    0.86932987134860684, 0.8685707059713409, 0.86780949676330332,
    0.86704624551569265, 0.866280954024513, 0.86551362409056909,
    0.86474425751946238, 0.8639728561215867, 0.86319942171212416,
    0.8624239561110405, 0.8616464611430813, 0.86086693863776731,
    0.86008539042939014, 0.85930181835700847, 0.85851622426444274,
    0.85772861000027212, 0.85693897741782876, 0.85614732837519447,
    0.855353664735196, 0.85455798836540053, 0.85376030113811141,
    0.85296060493036363, 0.85215890162391983, 0.8513551931052652,
    0.85054948126560348, 0.84974176800085255, 0.84893205521163961,
    0.84812034480329723, 0.84730663868585832, 0.84649093877405213,
    0.84567324698729907, 0.84485356524970712, 0.84403189549006641,
    0.84320823964184544, 0.84238259964318585, 0.84155497743689844,
    0.84072537497045807, 0.83989379419599952, 0.83906023707031274,
    0.83822470555483808, 0.83738720161566194, 0.836547727223512,
    0.8357062843537526, 0.83486287498638, 0.83401750110601813,
    0.83317016470191319, 0.83232086776792968, 0.83146961230254524,
    0.83061640030884631, 0.829761233794523, 0.82890411477186487,
    0.8280450452577558, 0.82718402727366913, 0.82632106284566353,
    0.82545615400437755, 0.82458930278502529, 0.82372051122739143,
    0.82284978137582643, 0.82197711527924155, 0.82110251499110465,
    0.82022598256943469, 0.819347520076797, 0.81846712958029866,
    0.81758481315158371, 0.81670057286682785, 0.81581441080673378,
    0.81492632905652662, 0.81403632970594841, 0.81314441484925359,
    0.81225058658520388, 0.81135484701706373, 0.81045719825259477,
    0.80955764240405126, 0.808656181588175, 0.80775281792619036,
    0.80684755354379933, 0.80594039057117628, 0.80503133114296366,
    0.8041203773982657, 0.80320753148064494, 0.80229279553811572,
    0.80137617172314024, 0.80045766219262282, 0.799537269107905,
    0.79861499463476093, 0.79769084094339116, 0.79676481020841883,
    0.79583690460888357, 0.794907126328237, 0.79397547755433717,
    0.79304196047944364, 0.79210657730021239, 0.7911693302176902,
    0.79023022143731, 0.78928925316888565, 0.78834642762660634,
    0.78740174702903143, 0.78645521359908577, 0.78550682956405393,
    0.78455659715557524, 0.78360451860963831, 0.78265059616657573,
    0.78169483207105939, 0.78073722857209449, 0.77977778792301455,
    0.778816512381476, 0.77785340420945315, 0.77688846567323244,
    0.77592169904340769, 0.77495310659487393, 0.7739826906068229,
    0.773010453362737, 0.77203639715038452, 0.77106052426181382,
    0.7700828369933479, 0.7691033376455797, 0.76812202852336542,
    0.7671389119358204, 0.76615399019631292, 0.765167265622459,
    0.76417874053611679, 0.76318841726338127, 0.7621962981345789,
    0.76120238548426178, 0.76020668165120242, 0.759209188978388,
    0.75820990981301528, 0.75720884650648457, 0.75620600141439454,
    0.75520137689653655, 0.75419497531688917, 0.75318679904361252,
    0.7521768504490427, 0.75116513190968637, 0.75015164580621507,
    0.74913639452345937, 0.7481193804504036, 0.74710060598018013,
    0.74608007351006378, 0.745057785441466, 0.74403374417992929,
    0.74300795213512172, 0.74198041172083107, 0.74095112535495922,
    0.7399200954595162, 0.73888732446061511, 0.737852814788466,
    0.73681656887736979, 0.73577858916571359, 0.7347388780959635,
    0.73369743811466037, 0.73265427167241282, 0.73160938122389263,
    0.73056276922782759, 0.729514438146997, 0.7284643904482252,
    0.72741262860237577, 0.726359155084346, 0.72530397237306077,
    0.724247082951467, 0.72318848930652746, 0.72212819392921535,
    0.72106619931450811, 0.72000250796138165, 0.71893712237280449,
    0.71787004505573171, 0.71680127852109954, 0.71573082528381859,
    0.71465868786276909, 0.71358486878079352, 0.71250937056469243,
    0.71143219574521643, 0.71035334685706242, 0.70927282643886569,
    0.7081906370331954, 0.70710678118654757, 0.70602126144933974,
    0.70493408037590488, 0.70384524052448494, 0.7027547444572253,
    0.70166259474016845, 0.70056879394324834, 0.69947334464028377,
    0.69837624940897292, 0.69727751083088652, 0.696177131491463,
    0.69507511398000088, 0.69397146088965389, 0.69286617481742463,
    0.69175925836415775, 0.6906507141345346, 0.68954054473706683,
    0.68842875278409044, 0.687315340891759, 0.68620031168003859,
    0.68508366777270036, 0.6839654117973154, 0.68284554638524808,
    0.68172407417164971, 0.680600997795453, 0.679476319899365,
    0.67835004312986147, 0.67722217013718033, 0.67609270357531592,
    0.674961646102012, 0.673829000378756, 0.67269476907077286,
    0.67155895484701833, 0.67042156038017309, 0.669282588346636,
    0.66814204142651845, 0.66699992230363747, 0.66585623366550972,
    0.66471097820334479, 0.66356415861203977, 0.66241577759017178,
    0.66126583783999227, 0.66011434206742048, 0.65896129298203732,
    0.65780669329707864, 0.65665054572942894, 0.65549285299961535,
    0.65433361783180044, 0.65317284295377676, 0.6520105310969595,
    0.650846684996381, 0.64968130739068319, 0.64851440102211244,
    0.64734596863651206, 0.64617601298331628, 0.64500453681554393,
    0.64383154288979139, 0.64265703396622686, 0.641481012808583,
    0.64030348218415167, 0.63912444486377573, 0.637943903621844,
    0.6367618612362842, 0.63557832048855611, 0.63439328416364549,
    0.63320675505005719, 0.63201873593980906, 0.63082922962842447,
    0.629638238914927, 0.6284457666018326, 0.62725181549514408,
    0.62605638840434352, 0.62485948814238634, 0.62366111752569453,
    0.62246127937415, 0.62125997651108755, 0.6200572117632891,
    0.61885298796097632, 0.61764730793780387, 0.61644017453085365,
    0.61523159058062682, 0.61402155893103849, 0.61281008242940971,
    0.61159716392646191, 0.61038280627630948, 0.60916701233645321,
    0.60794978496777363, 0.60673112703452448, 0.60551104140432555,
    0.604289530948156, 0.60306659854034816, 0.60184224705858, 0.600616479383869,
    0.59938929840056454, 0.59816070699634238, 0.59693070806219639,
    0.59569930449243336, 0.59446649918466443, 0.5932322950397998,
    0.591996694962041, 0.59075970185887416, 0.58952131864106394,
    0.58828154822264522, 0.587040393520918, 0.58579785745643886,
    0.58455394295301533, 0.58330865293769829, 0.58206199034077544,
    0.58081395809576453, 0.57956455913940563, 0.57831379641165559,
    0.57706167285567944, 0.57580819141784534, 0.57455335504771576,
    0.5732971666980422, 0.572039629324757, 0.57078074588696726,
    0.56952051934694714, 0.56825895267013149, 0.56699604882510868,
    0.56573181078361312, 0.5644662415205195, 0.56319934401383409,
    0.56193112124468936, 0.560661576197336, 0.55939071185913614,
    0.5581185312205561, 0.5568450372751601, 0.55557023301960218,
    0.55429412145362, 0.55301670558002747, 0.55173798840470734,
    0.55045797293660481, 0.54917666218771966, 0.54789405917310019,
    0.54661016691083486, 0.54532498842204646, 0.54403852673088382,
    0.54275078486451589, 0.54146176585312344, 0.54017147272989285,
    0.53887990853100842, 0.53758707629564539, 0.53629297906596318,
    0.53499761988709715, 0.533701001807153, 0.5324031278771979,
    0.531104001151255, 0.52980362468629461, 0.52850200154222848,
    0.52719913478190128, 0.52589502747108463, 0.524589682678469,
    0.52328310347565643, 0.52197529293715439, 0.52066625414036716,
    0.51935599016558964, 0.51804450409599934, 0.51673179901764987,
    0.51541787801946293, 0.51410274419322166, 0.512786400633563,
    0.5114688504379703, 0.51015009670676681, 0.508830142543107,
    0.50750899105297087, 0.50618664534515523, 0.50486310853126759,
    0.50353838372571758, 0.50221247404571079, 0.50088538261124071,
    0.49955711254508184, 0.49822766697278181, 0.49689704902265447,
    0.49556526182577254, 0.49423230851595967, 0.49289819222978404,
    0.4915629161065499, 0.49022648328829116, 0.48888889691976317,
    0.487550160148436, 0.48621027612448642, 0.48486924800079106,
    0.48352707893291874, 0.48218377207912272, 0.48083933060033396,
    0.47949375766015295, 0.478147056424843, 0.47679923006332209,
    0.47545028174715587, 0.47410021465054997, 0.47274903195034279,
    0.47139673682599764, 0.47004333245959562, 0.46868882203582796,
    0.46733320874198842, 0.46597649576796618, 0.46461868630623782,
    0.46325978355186015, 0.46189979070246273, 0.46053871095824,
    0.45917654752194409, 0.45781330359887717, 0.45644898239688392,
    0.45508358712634384, 0.45371712100016387, 0.45234958723377089,
    0.45098098904510386, 0.44961132965460654, 0.44824061228521989,
    0.44686884016237416, 0.44549601651398174, 0.4441221445704292,
    0.44274722756457, 0.44137126873171667, 0.43999427130963326,
    0.43861623853852766, 0.43723717366104409, 0.43585707992225547,
    0.43447596056965565, 0.43309381885315196, 0.43171065802505726,
    0.43032648134008261, 0.42894129205532949, 0.42755509343028208,
    0.42616788872679962, 0.42477968120910881, 0.42339047414379605,
    0.42200027079979968, 0.42060907444840251, 0.41921688836322391,
    0.41782371582021227, 0.41642956009763715, 0.41503442447608163,
    0.4136383122384345, 0.41224122666988289, 0.41084317105790391,
    0.40944414869225759, 0.40804416286497869, 0.40664321687036903,
    0.40524131400498986, 0.40383845756765407, 0.40243465085941843,
    0.40102989718357562, 0.39962419984564679, 0.39821756215337356,
    0.39680998741671031, 0.39540147894781635, 0.3939920400610481,
    0.39258167407295147, 0.39117038430225387, 0.38975817406985641,
    0.38834504669882625, 0.38693100551438858, 0.38551605384391885,
    0.38410019501693504, 0.38268343236508978, 0.38126576922216238,
    0.37984720892405116, 0.37842775480876556, 0.37700741021641826,
    0.37558617848921722, 0.37416406297145793, 0.37274106700951576,
    0.37131719395183749, 0.3698924471489341, 0.36846682995337232,
    0.36704034571976718, 0.36561299780477385, 0.36418478956707989,
    0.36275572436739723, 0.36132580556845428, 0.35989503653498811,
    0.35846342063373654, 0.35703096123343, 0.35559766170478385,
    0.35416352542049034, 0.35272855575521073, 0.35129275608556709,
    0.34985612979013492, 0.34841868024943456, 0.34698041084592368,
    0.34554132496398909, 0.34410142598993881, 0.34266071731199438,
    0.34121920232028236, 0.33977688440682685, 0.33833376696554113,
    0.33688985339222005, 0.3354451470845316, 0.33399965144200938,
    0.33255336986604422, 0.33110630575987643, 0.32965846252858749,
    0.3282098435790925, 0.32676045232013173, 0.32531029216226293,
    0.32385936651785285, 0.32240767880106985, 0.32095523242787521,
    0.31950203081601569, 0.31804807738501495, 0.31659337555616585,
    0.31513792875252244, 0.31368174039889152, 0.31222481392182488,
    0.31076715274961147, 0.30930876031226873, 0.30784964004153487,
    0.30638979537086092, 0.30492922973540237, 0.30346794657201132,
    0.30200594931922808, 0.30054324141727345, 0.29907982630804048,
    0.2976157074350862, 0.29615088824362379, 0.29468537218051433,
    0.29321916269425863, 0.29175226323498926, 0.29028467725446233,
    0.28881640820604948, 0.28734745954472951, 0.28587783472708062,
    0.28440753721127188, 0.28293657045705539, 0.28146493792575794,
    0.27999264308027322, 0.27851968938505306, 0.2770460803060999,
    0.27557181931095814, 0.27409690986870638, 0.272621355449949,
    0.271145159526808, 0.26966832557291509, 0.26819085706340318,
    0.26671275747489837, 0.26523403028551179, 0.26375467897483135,
    0.26227470702391359, 0.26079411791527551, 0.25931291513288623,
    0.257831102162159, 0.25634868248994291, 0.25486565960451457,
    0.25338203699557016, 0.25189781815421697, 0.25041300657296522,
    0.24892760574572015, 0.24744161916777327, 0.24595505033579459,
    0.24446790274782415, 0.24298017990326387, 0.24149188530286933,
    0.2400030224487415, 0.23851359484431842, 0.2370236059943672,
    0.23553305940497549, 0.23404195858354343, 0.23255030703877524,
    0.23105810828067111, 0.22956536582051887, 0.22807208317088573,
    0.22657826384561, 0.22508391135979283, 0.22358902922979, 0.22209362097320351,
    0.22059769010887351, 0.2191012401568698, 0.21760427463848364,
    0.21610679707621952, 0.21460881099378676, 0.21311031991609136,
    0.21161132736922755, 0.21011183688046961, 0.20861185197826349,
    0.20711137619221856, 0.20561041305309924, 0.20410896609281687,
    0.20260703884442113, 0.2011046348420919, 0.19960175762113097,
    0.19809841071795356, 0.19659459767008022, 0.19509032201612825,
    0.19358558729580361, 0.19208039704989244, 0.19057475482025274,
    0.18906866414980619, 0.1875621285825296, 0.18605515166344663,
    0.18454773693861962, 0.18303988795514095, 0.18153160826112497,
    0.18002290140569951, 0.17851377093899751, 0.17700422041214875,
    0.17549425337727143, 0.17398387338746382, 0.17247308399679595,
    0.17096188876030122, 0.16945029123396796, 0.16793829497473117,
    0.1664259035404641, 0.16491312048996992, 0.16339994938297323,
    0.16188639378011183, 0.16037245724292828, 0.15885814333386145,
    0.15734345561623825, 0.15582839765426523, 0.1543129730130201,
    0.15279718525844344, 0.15128103795733022, 0.14976453467732151,
    0.14824767898689603, 0.14673047445536175, 0.14521292465284746,
    0.14369503315029447, 0.14217680351944803, 0.14065823933284921,
    0.1391393441638262, 0.13762012158648604, 0.1361005751757062,
    0.13458070850712617, 0.13306052515713906, 0.13154002870288312,
    0.13001922272223335, 0.12849811079379317, 0.12697669649688587,
    0.12545498341154623, 0.12393297511851216, 0.1224106751992162,
    0.12088808723577708, 0.11936521481099135, 0.11784206150832498,
    0.11631863091190475, 0.11479492660651008, 0.11327095217756435,
    0.11174671121112659, 0.11022220729388306, 0.10869744401313872,
    0.10717242495680884, 0.10564715371341062, 0.10412163387205459,
    0.10259586902243628, 0.10106986275482782, 0.099543618660069319,
    0.0980171403295606, 0.096490431355252593, 0.094963495329638992,
    0.093436335845747787, 0.091908956497132724, 0.090381360877864983,
    0.0888535525825246, 0.087325535206192059, 0.0857973123444399,
    0.084268887593324071, 0.082740264549375692, 0.081211446809592441,
    0.079682437971430126, 0.078153241632794232, 0.076623861392031492,
    0.0750943008479213, 0.073564563599667426, 0.072034653246889332,
    0.070504573389613856, 0.068974327628266746, 0.067443919563664051,
    0.0659133527970038, 0.064382630929857465, 0.0628517575641614,
    0.061320736302208578, 0.059789570746639868, 0.058258264500435752,
    0.056726821166907748, 0.055195244349689941, 0.05366353765273052,
    0.052131704680283324, 0.050599749036899282, 0.049067674327418015,
    0.0475354841569593, 0.046003182130914623, 0.044470771854938668,
    0.04293825693494082, 0.041405640977076739, 0.039872927587739811,
    0.038340120373552694, 0.036807222941358832, 0.035274238898213947,
    0.03374117185137758, 0.032208025408304586, 0.030674803176636626,
    0.029141508764193722, 0.02760814577896574, 0.0260747178291039,
    0.024541228522912288, 0.023007681468839369, 0.021474080275469508,
    0.019940428551514441, 0.01840672990580482, 0.01687298794728171,
    0.0153392062849881, 0.013805388528060391, 0.012271538285719925,
    0.010737659167264491, 0.00920375478205982, 0.007669828739531097,
    0.0061358846491544753, 0.0046019261204485705, 0.0030679567629659761,
    0.0015339801862847655, 0.0 };

  double sum;
  emxArray_real_T *mfc_spectrum_bins;
  int i;
  emxArray_real_T *log_spectral;
  static const signed char iv0[18] = { 2, 5, 7, 10, 13, 16, 20, 24, 28, 32, 38,
    43, 49, 56, 63, 71, 80, 90 };

  int j;
  emxArray_creal_T *c_y1;
  int b_log_spectral[2];
  emxArray_real_T c_log_spectral;

  /*  EECS351 Final Project */
  /*  rohrer, malinas */
  /*  November 6, 2016 */
  /*  mfcc: compute the MFCC coefficients for a 30 ms window and store all of  */
  /*        these results in matrix where each row is a vector of MFCC coefficients */
  /*  framing the signal to prevent spectral leakage */
  /*  calculate the hamming window */
  /*  apply the hamming window */
  /*  take the DFT of the current frame */
  for (i0 = 0; i0 < 1323; i0++) {
    mag_dft_currFrame[i0] = currFrame * dv0[i0];
  }

  bluestein_setup(wwc);
  xidx = 0;
  for (center = 0; center < 1323; center++) {
    b_y1[center].re = wwc[center + 1322].re * mag_dft_currFrame[xidx];
    b_y1[center].im = wwc[center + 1322].im * -mag_dft_currFrame[xidx];
    xidx++;
  }

  r2br_r2dit_trig_impl(b_y1, costab, sintab, fy);
  r2br_r2dit_trig(wwc, costab, sintab, fv);
  for (i0 = 0; i0 < 4096; i0++) {
    sum = fy[i0].re;
    fy[i0].re = fy[i0].re * fv[i0].re - fy[i0].im * fv[i0].im;
    fy[i0].im = sum * fv[i0].im + fy[i0].im * fv[i0].re;
  }

  b_r2br_r2dit_trig(fy, costab, sintabinv, fv);
  xidx = 0;
  for (center = 0; center < 1323; center++) {
    b_y1[xidx].re = wwc[center + 1322].re * fv[center + 1322].re + wwc[center +
      1322].im * fv[center + 1322].im;
    b_y1[xidx].im = wwc[center + 1322].re * fv[center + 1322].im - wwc[center +
      1322].im * fv[center + 1322].re;
    xidx++;
  }

  /*  put into mel freq banks */
  /*  calculate magnitude of DFT coefficients */
  for (center = 0; center < 1323; center++) {
    sum = rt_hypotd_snf(b_y1[center].re, b_y1[center].im);
    mag_dft_currFrame[center] = sum * sum;
  }

  emxInit_real_T(&mfc_spectrum_bins, 2);

  /*  num freq bins in DFT */
  /*  remove redundant information from DFT */
  /*  Compute the Mel filters */
  /*  in hertz */
  /*  in hertz */
  /*  somewhat arbitrarily chosen */
  /*  compute max mel freq */
  /*  compute min mel freq */
  /*  linearly space the mel bins */
  /*  find the mel filters */
  /*  these now correspond to the index in the DFT of each bin */
  /*  now the triangular binning */
  /*  iterate over every bin */
  i0 = mfc_spectrum_bins->size[0] * mfc_spectrum_bins->size[1];
  mfc_spectrum_bins->size[0] = 1;
  mfc_spectrum_bins->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)mfc_spectrum_bins, i0, (int)sizeof
                    (double));
  mfc_spectrum_bins->data[0] = 0.0;
  for (i = 0; i < 16; i++) {
    xidx = iv0[i];
    center = iv0[i + 1];
    sum = 0.0;

    /*  iterate over left side individual bin to the center */
    i0 = iv0[i + 1] - iv0[i];
    for (j = 0; j < i0; j++) {
      sum += ((1.0 + (double)j) - 1.0) / (double)(center - xidx) *
        mag_dft_currFrame[(xidx + j) - 1];
    }

    /*  add the middle term */
    sum += mag_dft_currFrame[iv0[i + 1] - 1];

    /*  iterate over left side individual bin to the center */
    xidx = iv0[2 + i] - iv0[i + 1];
    for (j = 0; j < xidx; j++) {
      sum += ((double)(xidx - j) - 1.0) / (double)xidx *
        mag_dft_currFrame[center + j];
    }

    xidx = mfc_spectrum_bins->size[1];
    i0 = mfc_spectrum_bins->size[0] * mfc_spectrum_bins->size[1];
    mfc_spectrum_bins->size[1] = xidx + 1;
    emxEnsureCapacity((emxArray__common *)mfc_spectrum_bins, i0, (int)sizeof
                      (double));
    mfc_spectrum_bins->data[xidx] = sum;
  }

  emxInit_real_T(&log_spectral, 2);

  /*  compute the log-spectral vector */
  i0 = log_spectral->size[0] * log_spectral->size[1];
  log_spectral->size[0] = 1;
  log_spectral->size[1] = mfc_spectrum_bins->size[1];
  emxEnsureCapacity((emxArray__common *)log_spectral, i0, (int)sizeof(double));
  xidx = mfc_spectrum_bins->size[0] * mfc_spectrum_bins->size[1];
  for (i0 = 0; i0 < xidx; i0++) {
    log_spectral->data[i0] = mfc_spectrum_bins->data[i0];
  }

  for (center = 0; center + 1 <= mfc_spectrum_bins->size[1]; center++) {
    log_spectral->data[center] = log(log_spectral->data[center]);
  }

  emxFree_real_T(&mfc_spectrum_bins);
  emxInit_creal_T(&c_y1, 1);

  /*  get the mel coefficients by take the IFFT */
  /*    in this case reduces to DCT as signal real */
  b_log_spectral[0] = log_spectral->size[1];
  b_log_spectral[1] = 1;
  c_log_spectral = *log_spectral;
  c_log_spectral.size = (int *)&b_log_spectral;
  c_log_spectral.numDimensions = 1;
  eml_fft(&c_log_spectral, log_spectral->size[1], c_y1);
  i0 = mfcc_coeff->size[0] * mfcc_coeff->size[1];
  mfcc_coeff->size[0] = 1;
  mfcc_coeff->size[1] = log_spectral->size[1];
  emxEnsureCapacity((emxArray__common *)mfcc_coeff, i0, (int)sizeof(creal_T));
  xidx = log_spectral->size[1];
  emxFree_real_T(&log_spectral);
  for (i0 = 0; i0 < xidx; i0++) {
    mfcc_coeff->data[i0] = c_y1->data[i0];
  }

  emxFree_creal_T(&c_y1);

  /*  mfcc */
}

/*
 * File trailer for mfcc.c
 *
 * [EOF]
 */
