#include <stddef.h>
#include "parameters.h"
#include "states.h"

states::states() {}

states states::operator+(const states& a)
{
  states out;

  for(int g = 0; g < NG; g++)
    for(int r = 0; r <= RG; r++){
      for(int m  = 0; m < DS; m++){
	for(int u  = 0; u < art_index[m]; u++){
	  out.X[g][r][m][u] = X[g][r][m][u] + a.X[g][r][m][u];
	  #if fullmod
	  if(r < RG)
	    out.trans[g][r][m][u] = trans[g][r][m][u] + a.trans[g][r][m][u];
	  #endif
	}
	#if fullmod
	out.art_init[g][r][m] = art_init[g][r][m] + a.art_init[g][r][m];
	#endif
      }
      #if fullmod
      if(r < RG)
	out.inc[g][r] = inc[g][r] + a.inc[g][r];
      #endif
    }

  for(int m  = 0; m < DS; m++)
    out.art_index[m] = art_index[m];

  return out;
}


states states::operator*(const double& c)
{
  states out;

  for(int g = 0; g < NG; g++)
    for(int r = 0; r <= RG; r++){
      for(int m  = 0; m < DS; m++){
	for(int u  = 0; u < art_index[m]; u++){
	  out.X[g][r][m][u] = X[g][r][m][u] * c;
	  #if fullmod
	  if(r < RG)
	    out.trans[g][r][m][u] = trans[g][r][m][u] * c;
	  #endif
	}
	#if fullmod
	out.art_init[g][r][m] = art_init[g][r][m] * c;
	#endif
      }
      #if fullmod
      if(r < RG)
	out.inc[g][r] = inc[g][r] * c;
      #endif
    }

  for(int m  = 0; m < DS; m++)
    out.art_index[m] = art_index[m];

  return out;
}

states & states::operator+=(const states& a)
{
  *this = *this + a;
  return *this;
}


