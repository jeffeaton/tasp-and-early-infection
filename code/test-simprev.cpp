#include <stddef.h>
#include <stdio.h>
#include <math.h>

#include "simprev.h"

// $ g++ -lgsl -lgslcblas test-simprev.cpp simprev.cpp model.cpp states.cpp 

int main(int argc, char *argv[])
{

  double theta[numParam] = {1986.989968056878979, 0.248537739624382, 0.484384635623043, 0.359295751121344, 0.361065867767938, 2.313501394133991, 67.547262910640569, 46.497647430714409, 0.706899939540192, 0.194808796557237, 0.418877554277150, 0.585004044217372, 0.048652671124887, 0.313264144122407, 2000.445334431683023, 2009.465323351676943, 26.928118355200720};    

  double evec[30];
  printf("R0 = %.4f\n\n", R0(theta, -INFINITY, evec));
  printf("r = %.4f\n\n", epidemicGrowthRate(theta, -INFINITY, evec));
  
  for(size_t i = 0; i < 30; i++)
    printf("%8.5f ", evec[i]);
  printf("\n");

  struct modprev out;
  SimPrev(theta, &out);
  printf("\n$prev\n");
  for(size_t i = 0; i < 21; i++)
    printf("%.9f %.9f\n", out.male[i], out.female[i]);
  
  printf("\n$on.art\n");
  
  for(size_t i = 0; i < 21; i++){
    if(i == 15) printf("\n");
    printf("%.10f ", out.artcov[i]);
  }

  printf("\n");


  return 0;
}
