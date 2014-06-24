#include <stddef.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

#include "likelihood.h"

int main(int argc, char *argv[])
{

  unsigned int rng_seed = 334671461;

  // Declare and configure GSL RNG
  gsl_rng * rng;
  const gsl_rng_type * T;

  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  gsl_rng_set(rng, rng_seed);

  size_t NumParam = 17, InitSamples = 100;

  gsl_matrix * Xmat = gsl_matrix_alloc(InitSamples, NumParam);

  sample_prior(rng, InitSamples, Xmat);

  printf("rel.priminf <- c(");
  for(size_t i = 0; i < InitSamples; i++){
    printf("%f", gsl_matrix_get(Xmat, i, NumParam - 1));
    if(i != InitSamples - 1) printf(", ");
  }

  printf(")\n\nprior <- c(");
  for(size_t i = 0; i < InitSamples; i++){
    gsl_vector_view v = gsl_matrix_row(Xmat, i);
    printf("%f", prior(&v.vector));
    if(i != InitSamples - 1) printf(", ");
  }

  printf(")\n\nll <- c(");
  for(size_t i = 0; i < InitSamples; i++){
    gsl_vector_view v = gsl_matrix_row(Xmat, i);
    printf("%f", log(likelihood(&v.vector)));
    if(i != InitSamples - 1) printf(", ");
  }

  printf(")\n\n");

  return 0;
}
