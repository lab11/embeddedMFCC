/* vim: set sw=2 expandtab tw=80: */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// add fft libraries
#include "fft_lib/fft.h"
// add mfcc libraries
#include "libmfcc.h"

#define numel 128

void nop() {}

fft_complex_t data[numel];
double spectrum[numel];

int main() {

  // test out the FFT
  int bits = log2(numel);

  // fill the array with data
  for (int i = 0; i < numel; i++) {
    int32_t real = numel - i;
    if (i % 2 == 0) {
      data[i].r = real;
    } else {
      data[i].r = -real;
    }
  }

  printf("data created\n");
  fft_permutate(data, bits);
  fft_forward(data, bits);

  printf("after FFT\n");
  // make FFT results purely real
  for (int i = 0; i < numel; i++) {
    spectrum[i] = data[i].r;
  }
  printf("after making spectrum real\n");

  // compute the first 13 coefficients
  int numbins = 10;
  for (int coeff = 0; coeff < numbins; coeff++) {
    double mfcc_result = GetCoefficient(spectrum, 22050, numbins, numel, coeff);
    printf("%i %i\n", coeff, (int)(1000 * mfcc_result));
  }
  printf("after mfcc\n");

  return 0;
}
