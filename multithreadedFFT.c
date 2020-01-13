#include <complex.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

int N;
int P;
int size;
double *x;
double complex *bitrev_data;

int rev_bits(unsigned int index) {
  int rev = 0;
  int mySize = size;
  for (; mySize > 1; mySize >>= 1) {
    rev = (rev << 1) | (index & 1);
    index >>= 1;
  }
  return rev;
}

void *thRev(void *var) {
  int thread_id = *(int *)var;
  unsigned long start, end;
  start = thread_id * ceil((double)size / P);
  end = fmin(size, (thread_id + 1) * ceil((double)size / P));

  for (int i = start; i < end; i++) {
    int rb = rev_bits(i);
    bitrev_data[rb] = x[i];
  }
}

unsigned long imax, istep;
double *wr, *wi;

void initConsts() {
  double wtemp, twr, wpr, wpi, twi, theta;
  theta = -2.0 * M_PI / istep;
  wtemp = sin(0.5 * theta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(theta);
  twr = 1.0;
  twi = 0.0;

  wr[0] = twr;
  wi[0] = twi;

  for (int i = 1; i <= size / 2; i++) {
    wtemp = twr;
    twr = wtemp * wpr - twi * wpi + twr;
    twi = twi * wpr + wtemp * wpi + twi;
    wr[i] = twr;
    wi[i] = twi;
  }
}

void *thOutward(void *var) {
  int thread_id = *(int *)var;

  unsigned long start, end;
  start = thread_id * ceil((double)imax / P);
  end = fmin(size, (thread_id + 1) * ceil((double)imax / P));
  end = fmin(end, imax);

  double complex tc;

  for (unsigned long m = start; m < end; m++) {
    for (unsigned long i = m; i < size; i += istep) {
      unsigned long j = i + imax;
      tc = wr[m] * creal(bitrev_data[j]) - wi[m] * cimag(bitrev_data[j]) +
           (wr[m] * cimag(bitrev_data[j]) + wi[m] * creal(bitrev_data[j])) * I;
      bitrev_data[j] = bitrev_data[i] - tc;
      bitrev_data[i] += tc;
    }
  }
}

int inner_m;

void *thInward(void *var) {
  int thread_id = *(int *)var;
  unsigned long start, end;
  start = thread_id * ceil((double)size / P);
  end = fmin(size, (thread_id + 1) * ceil((double)size / P));
  end = fmin(end, size);

  // printf("FROM TID%d -> start: %ld, end:%ld\n", thread_id, start, end);

  for (int i = start; i < end; i += istep) {
    unsigned long j = i + imax;
    unsigned complex tc =
        wr[0] * creal(bitrev_data[j]) - wi[0] * cimag(bitrev_data[j]) +
        (wr[0] * cimag(bitrev_data[j]) + wi[0] * creal(bitrev_data[j])) * I;
    bitrev_data[j] = bitrev_data[i] - tc;
    bitrev_data[i] += tc;
  }
}

void fft(pthread_t tid[], int thread_id[]) {
  double complex tc;
  unsigned long i, j, m;
  double wtemp, wr, wpr, wpi, wi, theta;

  imax = 1;
  istep = 2;

  while (imax < size) {
    istep = 2 * imax;

    theta = -2.0 * M_PI / istep;
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;

    if (P > 1 && imax >= P) {
      for (int i = 0; i < P; i++) {
        pthread_create(&(tid[i]), NULL, thOutward, &(thread_id[i]));
      }

      for (int i = 0; i < P; i++) pthread_join(tid[i], NULL);

    } else {
      for (m = 0; m < imax; m++) {
        if (P > 1 && imax == 1) {
          inner_m = m;
          for (int i = 0; i < P; i++)
            pthread_create(&(tid[i]), NULL, thInward, &(thread_id[i]));

          for (int i = 0; i < P; i++) pthread_join(tid[i], NULL);
        } else {
          for (i = m; i < size; i += istep) {
            j = i + imax;
            tc = wr * creal(bitrev_data[j]) - wi * cimag(bitrev_data[j]) +
                 (wr * cimag(bitrev_data[j]) + wi * creal(bitrev_data[j])) * I;
            bitrev_data[j] = bitrev_data[i] - tc;
            bitrev_data[i] += tc;
          }

          wtemp = wr;
          wr = wtemp * wpr - wi * wpi + wr;
          wi = wi * wpr + wtemp * wpi + wi;
        }
      }
    }
    imax = istep;
  }
}

int main(int argc, char *argv[]) {
  int N = 0, r, i;

  P = atoi(argv[3]);

  pthread_t tid[P];
  int thread_id[P];
  for (i = 0; i < P; i++) thread_id[i] = i;

  FILE *f = fopen(argv[1], "r");

  r = fscanf(f, "%d", &N);
  size = N;

  x = (double *)malloc(N * sizeof(double));
  bitrev_data = (complex double *)malloc(N * sizeof(complex double));
  if (P > 1) {
    wr = (double *)malloc(N * sizeof(double));
    wi = (double *)malloc(N * sizeof(double));
    initConsts();
  }

  double myX = 0;
  for (i = 0; i < N; i++) {
    r = fscanf(f, "%lf", &myX);
    x[i] = myX;
  }
  fclose(f);

  FILE *fo = fopen(argv[2], "w+");

  fprintf(fo, "%d\n", N);

  for (i = 0; i < P; i++)
    pthread_create(&(tid[i]), NULL, thRev, &(thread_id[i]));

  for (i = 0; i < P; i++) pthread_join(tid[i], NULL);

  fft(tid, thread_id);

  for (i = 0; i < N; i++)
    fprintf(fo, "%lf %lf\n", creal(bitrev_data[i]), cimag(bitrev_data[i]));

  fclose(fo);
  return 0;
}