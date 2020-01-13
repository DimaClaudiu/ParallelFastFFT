#include <complex.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>

int P;
int N;
double complex x[5000];
double complex out[5000];

void FT(const double complex input[], double complex output[], int N) {
  int k, n;
  for (k = 0; k < N; k++) {
    output[k] = 0;
    for (n = 0; n < N; n++) {
      double expn = 2 * M_PI * n * k / N;
      output[k] += input[n] * cexp(-expn * I);
    }
  }
}

void *threadFunction(void *var) {
  int thread_id = *(int *)var;

  int start = thread_id * ceil((double)N / P);
  int end = fmin(N, (thread_id + 1) * ceil((double)N / P));
  end = fmin(end, N);
  int k, n;
  for (k = start; k < end; k++)
    for (n = 0; n < N; n++) {
      double expn = 2 * M_PI * n * k / N;
      out[k] += x[n] * cexp(-expn * I);
    }
}

int main(int argc, char *argv[]) {
  int r;

  P = atoi(argv[3]);
  pthread_t tid[P];
  int thread_id[P];
  int i;
  for (i = 0; i < P; i++) thread_id[i] = i;

  FILE *f = fopen(argv[1], "r");
  r = fscanf(f, "%d", &N);
  double myX = 0;
  for (i = 0; i < N; i++) {
    r = fscanf(f, "%lf", &myX);
    x[i] = myX + 0 * I;
  }
  fclose(f);

  // FT(x, out, N);
  for (i = 0; i < P; i++)
    pthread_create(&(tid[i]), NULL, threadFunction, &(thread_id[i]));

  for (i = 0; i < P; i++) pthread_join(tid[i], NULL);

  FILE *fo = fopen(argv[2], "w+");
  fprintf(fo, "%d\n", N);
  for (i = 0; i < N; i++)
    fprintf(fo, "%lf %lf\n", creal(out[i]), cimag(out[i]));
  fclose(fo);

  return 0;
}
