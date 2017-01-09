MFCC Benchmarks
================
Using the following code: [SYLT-FFT](https://github.com/stg/SYLT-FFT) and 
[libmfcc](https://github.com/jsawruk/libmfcc)

## Algorithm Overview
1. Call `fft_permutate` from SYLT-FFT 
2. Call `fft_forward` from SYLT-FFT
3. Keep only real elements of FFT
4. Call `GetCoefficient` from libmfcc for each desired coefficient

## Results
This serves a collection of MFCC computation benchmarks on a cortex-m4 based
system, the hail board.

- **N=256, Fs=44100 Hz, 20 coefficients**: for 0.005 second sample takes ~45 seconds.
- **N=256, Fs=22050 Hz, 20 coefficients**: for 0.010 second sample takes ~45 seconds.
- **N=256, Fs=44100 Hz, 13 coefficients**: for 0.005 second sample takes ~4.5 seconds.
- **N=256, Fs=22050 Hz, 13 coefficients**: for 0.010 second sample takes ~4.5 seconds.
- **N=256, Fs=44100 Hz, 10 coefficients**: for 0.005 second sample takes ~3 seconds.
- **N=256, Fs=22050 Hz, 10 coefficients**: for 0.010 second sample takes ~3 seconds.
- **N=128, Fs=44100 Hz, 10 coefficients**: for 0.0025 second sample takes ~1.8 seconds.
- **N=128, Fs=22050 Hz, 10 coefficients**: for 0.005 second sample takes ~1.8 seconds.

### After switching `double` to `float`
- **N=256, Fs=44100 Hz, 10 coefficients**: for 0.005 second sample takes ~2.5 seconds.
- **N=256, Fs=22050 Hz, 10 coefficients**: for 0.010 second sample takes ~2.5 seconds.
