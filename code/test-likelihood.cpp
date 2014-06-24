#include <stddef.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

#include "likelihood.h"

int main(int argc, char *argv[])
{

  const size_t NumParam = 17;

  double theta[NumParam] = {1986.989968056878979, 0.248537739624382, 0.484384635623043, 0.359295751121344, 0.361065867767938, 2.313501394133991, 67.547262910640569, 46.497647430714409, 0.706899939540192, 0.194808796557237, 0.418877554277150, 0.585004044217372, 0.048652671124887, 0.313264144122407, 2000.445334431683023, 2009.465323351676943, 26.928118355200720};

  gsl_vector_view v = gsl_vector_view_array(theta, NumParam);

  printf("prior = %f\n", prior(&v.vector));
  printf("ll = %f\n", log(likelihood(&v.vector)));

  return 0;
}
