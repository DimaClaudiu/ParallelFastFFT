# Parallel Fast Fourier Transform
A multithreaded implementation of the fast Fourier transform, written in C, using pthreads.

#### - Project Status: [Completed]

### Built with
  * C and pthreads 
  
## Project Description
  The program aims to calculate the fast fourier transform for a given input as fast as possible. The implementation starts from the
iterative version of the Fast Fourier Transform, where is loop is split on threads at different times in the execution. This "sweet-spot" approach stops bottle-necks from many processor cores or low input values.
  
  It scales quite well on any number of threads, not just powers of 2, and can handle really big inputs since recursion isn't used at all.

#### Time and CPU Usage
Average of 3 runs per test. Input size of 2^24:

| Threads       | Time (s)      | CPU Usage |
| ------------- |:-------------:| -----:|
| 1 Thread      | 14.83s | 99%
| 2 Threads      | 9.51s      |   161% |
| 3 Threads | 8.11s  |    200% |
| 4 Threads | 7.33s  |    244% |
| 8 Threads | 6.93s  |    391% |
    
    
In comparison, the Rossetta code version of FFT can't handle this input due to it's use of recursivity. On the largest input I could find (2^14), this implementation is about 346% faster using just 4 threads.

## Usage
  Compile with `gcc` and `-lm -lpthread` flags.
  Run with command line arguments: `inputfile outputfile nrOfThreads`
  
  `./a.out largeInput.in fft.out 8`

## Completion Date
2019, October 26
