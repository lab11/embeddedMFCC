# embeddedMFCC
- Can successfully run MFCCs on cortex-m4. 
- Would take about 3 minutes to compute FFT and MFCC on embedded device for a 
1 second audio sample.
- See complete results in the `app/` README

- Benchmark of FFT algorithm is 3 ms for 512 samples (32 ms at 16 kHz) run on
  a hail board running tock.
