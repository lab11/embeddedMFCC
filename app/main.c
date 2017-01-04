/* vim: set sw=2 expandtab tw=80: */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <console.h>

// add fft libraries
#include "fft_lib/fft.h"
#include "fft_lib/fpmath.h"


char hello[] = "Hello World!\r\n";
char bye[] = "Farewell world! I have computed a FFT!!\r\n";

void nop() {}

int main() {

  putnstr(hello, sizeof(hello));
  // test out the FFT
  unsigned bits = 256;
  
  int numel = 214;

  fft_complex_t data[numel];
  // fill the array with data
  for (int i = 0; i < numel; i++){
    int32_t real = i;
    int32_t imag = numel - i;
    data[i].r = real;
    data[i].i = imag;
 }
  
  fft_forward(data, bits); 

    
  putnstr(bye, sizeof(bye));
  return 0;
}
