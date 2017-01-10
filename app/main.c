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

#define numel 256

void nop() {}

fft_complex_t data[numel];
float spectrum[numel];

int main() {

  // for timing
  timer_oneshot(10000);
  unsigned int start = timer_read();

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

  // printf("data created\n");
  fft_permutate(data, bits);
  fft_forward(data, bits);

  // make FFT results purely real
  for (int i = 0; i < numel; i++) {
    spectrum[i] = data[i].r;
  }

  // compute the first 13 coefficients
  int numbins = 10;
  for (int coeff = 0; coeff < numbins; coeff++) {
    float mfcc_result = GetCoefficient(spectrum, 44100, numbins, numel, coeff);
  }
  // print time
  unsigned int end = timer_read();
  unsigned int freq = 16000;
  printf("end: %i, start: %i, time elapsed is (divide by 16 kHz for seconds): "
         "%i \n", end, start, (end - start));

  return 0;
}
