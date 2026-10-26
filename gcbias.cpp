#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gcbias.h"

gc_bias_t *gc_bias_read(const char *filename) {
  if (filename == NULL) return NULL;

  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "WARNING: failed to open GC bias file '%s', continuing without GC bias\n", filename);
    return NULL;
  }

  int cap = 256, n = 0;
  double *f = (double *)malloc(sizeof(double) * cap);
  double *w = (double *)malloc(sizeof(double) * cap);
  char line[4096];

  while (fgets(line, sizeof(line), fp)) {
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
    double a = 0.0, b = 0.0;
    if (sscanf(line, "%lf %lf", &a, &b) != 2) continue;
    if (a < 0.0) a = 0.0;
    if (a > 1.0) a = 1.0;
    if (b < 0.0) b = 0.0;
    if (n == cap) {
      cap *= 2;
      f = (double *)realloc(f, sizeof(double) * cap);
      w = (double *)realloc(w, sizeof(double) * cap);
    }
    f[n] = a;
    w[n] = b;
    n++;
  }
  fclose(fp);

  if (n < 2) {
    fprintf(stderr, "WARNING: GC bias file '%s' has fewer than 2 usable rows, ignoring\n", filename);
    free(f);
    free(w);
    return NULL;
  }

  // sort by GC fraction (insertion sort, tables are tiny)
  for (int i = 1; i < n; i++) {
    double fa = f[i], wa = w[i];
    int j = i - 1;
    while (j >= 0 && f[j] > fa) {
      f[j + 1] = f[j];
      w[j + 1] = w[j];
      j--;
    }
    f[j + 1] = fa;
    w[j + 1] = wa;
  }

  gc_bias_t *gb = (gc_bias_t *)malloc(sizeof(gc_bias_t));
  gb->frac = f;
  gb->weight = w;
  gb->n = n;

  double wmin = w[0], wmax = w[0];
  for (int i = 1; i < n; i++) {
    if (w[i] < wmin) wmin = w[i];
    if (w[i] > wmax) wmax = w[i];
  }
  fprintf(stderr,
          "\t-> GC bias table '%s': %d bins, GC fraction %.3f-%.3f, weight %.3f-%.3f (mean acceptance is reported at exit)\n",
          filename, n, f[0], f[n - 1], wmin, wmax);
  return gb;
}

double gc_bias_weight(const gc_bias_t *gb, double gc) {
  if (gb == NULL || gb->n < 2) return 1.0;
  if (gc <= gb->frac[0]) return gb->weight[0];
  if (gc >= gb->frac[gb->n - 1]) return gb->weight[gb->n - 1];

  int lo = 0, hi = gb->n - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (gb->frac[mid] <= gc) lo = mid;
    else hi = mid;
  }
  double x0 = gb->frac[lo], x1 = gb->frac[hi];
  if (x1 <= x0) return gb->weight[lo];
  double t = (gc - x0) / (x1 - x0);
  return gb->weight[lo] + t * (gb->weight[hi] - gb->weight[lo]);
}

void gc_bias_destroy(gc_bias_t *gb) {
  if (gb == NULL) return;
  free(gb->frac);
  free(gb->weight);
  free(gb);
}
