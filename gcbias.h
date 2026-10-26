#ifndef GCBIAS_H
#define GCBIAS_H

#ifdef __cplusplus
extern "C" {
#endif

/*
  gc_bias_t - lookup table mapping a fragment's GC fraction to an acceptance weight.

  File format (whitespace separated, '#' starts a comment line):
      <GC_fraction 0..1>   <weight 0..1>

  Rows must be sortable by GC fraction; they need not be sorted in the file,
  gc_bias_read() sorts them. Linear interpolation is used between rows and the
  first/last weight is clamped outside the tabulated range.
*/
typedef struct {
  double *frac;
  double *weight;
  int n;
} gc_bias_t;

gc_bias_t *gc_bias_read(const char *filename);
double gc_bias_weight(const gc_bias_t *gb, double gc_frac);
void gc_bias_destroy(gc_bias_t *gb);

#ifdef __cplusplus
}
#endif

#endif
