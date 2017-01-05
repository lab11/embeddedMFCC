/* vim: set sw=2 expandtab tw=80: */

#include <console.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// add fft libraries
#include "fft_lib/fft.h"
#include "fft_lib/fpmath.h"

#define numel 8

char hello[] = "Hello World!\r\n";
char bye[] = "Farewell world! I have computed a FFT!!\r\n";

void nop() {}

fft_complex_t data[numel];

void printArr() {
  for (int i = 0; i < numel; i++) {
    printf("arr[%d], real=%d, imag=%d\n", i, data[i].r, data[i].i);
  }
  return;
}

int main() {

  // test out the FFT
  int bits = log2(numel);

  // fill the array with data
  /*for (int i = 0; i < numel; i++) {
    int32_t real = 1;
    data[i].r = real;
  }*/
  data[0].r = 98;
  data[1].r = 8;
  data[2].r = 798;
  data[3].r = 90;
  data[4].r = 45;
  data[5].r = 103;
  data[6].r = 421;
  data[7].r = 59;

  printf("data created\n");
  printArr();
  fft_permutate(data, bits);
  printf("data permuted \n");
  printArr();

  fft_forward(data, bits);

  printf("after FFT\n");
  printArr();
  putnstr(bye, sizeof(bye));
  return 0;
}
