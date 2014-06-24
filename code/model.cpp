#include <stddef.h>
#include <math.h>
#include "parameters.h"
#include "model.h"
#include "states.h"

///////////////////////////////////////////////////////////////

/* numerical solvers */
states rk4(states y, double dt, const struct parameters * param)
{
  states k[4], out;
  k[0] = grad(y, param);
  k[1] = grad(y + k[0]*(.5*dt), param);
  k[2] = grad(y + k[1]*(.5*dt), param);
  k[3] = grad(y + k[2]*dt, param);

  out = y + (k[0] + k[1]*2 + k[2]*2 + k[3])*(dt/6);

  return out;
}

states euler(states y, const double dt, const struct parameters * param)
{
  states out = y + grad(y, param)*dt;
  return out;
}

///////////////////////////////////////////////////////////////

/* the model gradient (differential equations) */
states grad(const states y, const struct parameters * param)
{
  states out; // initialise states object to store output
  for(int m  = 0; m < DS; m++)
    out.art_index[m] = y.art_index[m];

  // sum population
  double Xr[NG][RG], Xg[NG], Xtot = 0;
  for(size_t g = 0; g < NG; g++){
    Xg[g] = 0;
    for(size_t r = 0; r < RG; r++){
      Xr[g][r] = y.X[g][r][0][0];
      for(size_t m = 1; m < DS; m++)
        for(size_t u = 0; u < y.art_index[m]; u++)
          Xr[g][r] += y.X[g][r][m][u];
      Xg[g] += Xr[g][r];
    }
    Xtot += Xg[g];
  }

  // calculate the force of infection
  double W[NG][RG];
  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r < RG; r++)
      W[g][r] = param->c[g][r]*Xr[g][r];

  double rho[NG][RG][RG];
  for(size_t g = 0; g < NG; g++){
    size_t gprime = (g==0)?1:0;
    for(size_t r1 = 0; r1 < RG; r1++){
      double sumW = 0;
      for(size_t rsum = 0; rsum < RG; rsum++)
        sumW += W[gprime][rsum];
      for(size_t r2 = 0; r2 < RG; r2++)
        rho[g][r1][r2] = param->epsilon*((r1==r2)?1:0) + (1-param->epsilon)*W[gprime][r2]/sumW;
    }
  }
  
  double D[RG][RG];
  for(size_t rm = 0; rm < RG; rm++)
    for(size_t rf = 0; rf < RG; rf++)
      D[rm][rf] = W[0][rm]*rho[0][rm][rf]/(W[1][rf]*rho[1][rf][rm]);

  for(size_t rm = 0; rm < RG; rm++)
    for(size_t rf = 0; rf < RG; rf++){
      rho[0][rm][rf] = rho[0][rm][rf]*pow(D[rm][rf], theta_g - 1.0);
      rho[1][rf][rm] = rho[1][rf][rm]*pow(D[rm][rf], theta_g);
    }

  // calculate incident infections
  #if fullmod
  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r <= RG; r++){
      if(r < RG){
        out.inc[g][r] = 0;
        out.trans[g][r][0][0] = 0;
      }
      for(size_t m = 0; m < DS; m++){
        out.art_init[g][r][m] = 0;
        if(m > 0 && r < RG)
          for(size_t u = 0; u < y.art_index[m]; u++)
            out.trans[g][r][m][u] = 0;
      }
    }
  #else
  double inc[NG][RG];
  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r < RG; r++)
      inc[g][r] = 0;
  #endif

  for(size_t g = 0; g < NG; g++){
    size_t gprime = (g==0)?1:0;
    for(size_t r = 0; r < RG; r++)
      for(size_t rprime = 0; rprime < RG; rprime++)
        for(size_t mprime = 1; mprime < DS; mprime++){
          for(size_t uprime = 0; uprime < y.art_index[mprime]; uprime++){
            #if fullmod
            double tmp_inf = y.X[g][r][0][0]*param->c[g][r]*rho[g][r][rprime]*y.X[gprime][rprime][mprime][uprime]/Xr[gprime][rprime]*(1-exp(-param->beta[gprime][mprime][uprime]*param->nL[r][rprime]));
            out.inc[g][r] += tmp_inf;
            out.trans[gprime][rprime][mprime][uprime] += tmp_inf;
            #else
            inc[g][r] += y.X[g][r][0][0]*param->c[g][r]*rho[g][r][rprime]*y.X[gprime][rprime][mprime][uprime]/Xr[gprime][rprime]*(1-exp(-param->beta[gprime][mprime][uprime]*param->nL[r][rprime]));
            #endif
          }
        }
  }

  // do model
  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r < RG; r++){
  
      // susceptibles: add mortality and disease mortality, remove infections and natural mortality
      #if fullmod
      out.X[g][r][0][0] = (alpha+nu)*param->prop_new[g][r]*Xtot/2 - out.inc[g][r] - nu*y.X[g][r][0][0];
      #else
      out.X[g][r][0][0] = (alpha+nu)*param->prop_new[g][r]*Xtot/2 - inc[g][r] - nu*y.X[g][r][0][0];
      #endif

      // disease stage 1 (early infection): add new infections, remove natural mortality
      #if fullmod
      out.X[g][r][1][0] = (1.0-param->frac_scr[g])*out.inc[g][r] - nu*y.X[g][r][1][0];
      out.X[g][r][1][1] = param->frac_scr[g]*out.inc[g][r] - nu*y.X[g][r][1][1];
      #else
      out.X[g][r][1][0] = (1.0-param->frac_scr[g])*inc[g][r] - nu*y.X[g][r][1][0];
      out.X[g][r][1][1] = param->frac_scr[g]*inc[g][r] - nu*y.X[g][r][1][1];
      #endif

      // disease progression for stages 2:DS: natural mortality
      for(size_t m = 2; m < DS; m++)
        for(size_t u = 0; u < 2; u++)
          out.X[g][r][m][u] =  -nu*y.X[g][r][m][u];

      // ART stage progression for 1:ART_ST, natural mortality
      if(ART_ST > 1)
        for(size_t m = 1; m < DS; m++)
          for(size_t u = 2; u < y.art_index[m]; u++){
            out.X[g][r][m][u] = -nu*y.X[g][r][m][u];
          }
    }

  // move risk groups
  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r < RG; r++){
      for(size_t rprime = 0; rprime < RG; rprime++)
        if(param->psi[g][r][rprime] > 0){
          for(size_t m = 0; m < DS; m++)
            for(size_t u = 0; u < y.art_index[m]; u++){
              out.X[g][r][m][u] -= param->psi[g][r][rprime]*y.X[g][r][m][u];
              out.X[g][rprime][m][u] += param->psi[g][r][rprime]*y.X[g][r][m][u];
            }
        }
    }

  // model the removed group
  for(size_t g = 0; g < NG; g++){
    for(size_t m = 0; m < DS; m++)
      for(size_t u = 0; u < y.art_index[m]; u++){
        out.X[g][RG][m][u] = - mu*y.X[g][RG][m][u];
        for(size_t r = 0; r < RG; r++)
          out.X[g][RG][m][u] += nu*y.X[g][r][m][u];
      }
  }

  // HIV stage progression and ART initiation
  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r <= RG; r++){
      for(size_t m = 1; m < DS; m++){
        out.X[g][r][m][0] += sigma[g][m-1][0]*y.X[g][r][m-1][0]  - (sigma[g][m][0] + param->lambda[g][m])*y.X[g][r][m][0];
        #if fullmod
        out.art_init[g][r][m] += param->lambda[g][m]*y.X[g][r][m][0] + param->lambda_scr[g][m]*y.X[g][r][m][1];
        #endif
        if(y.art_index[m] > 1)
          out.X[g][r][m][1] += sigma[g][m-1][0]*((y.art_index[m-1] > 1)?y.X[g][r][m-1][1]:0) - (sigma[g][m][0] + param->lambda_scr[g][m])*y.X[g][r][m][1];

        if(y.art_index[m] > 2){

          // ART initiation and disease stage progression
          out.X[g][r][m][2] += param->lambda[g][m]*y.X[g][r][m][0] + param->lambda_scr[g][m]*y.X[g][r][m][1] - (sigma[g][m][1] + eta[g][m][1])*y.X[g][r][m][2];
          out.X[g][r][m][3] += (1-xi[g][m])*sigma[g][m][1]*y.X[g][r][m][2] - (sigma[g][m][2] + eta[g][m][2])*y.X[g][r][m][3];
          out.X[g][r][m][4] += sigma[g][m][2]*y.X[g][r][m][3] - (sigma[g][m][3] + eta[g][m][3])*y.X[g][r][m][4];
          out.X[g][r][m][5] += sigma[g][m][3]*y.X[g][r][m][4] - (sigma[g][m][4] + eta[g][m][4])*y.X[g][r][m][5];
          out.X[g][r][m][6] += xi[g][m]*sigma[g][m][1]*y.X[g][r][m][2] + sigma[g][m][4]*y.X[g][r][m][5] - (sigma[g][m][5] + eta[g][m][5])*y.X[g][r][m][6];

          // disease progression for dropped out individuals
          out.X[g][r][m][7] += sigma[g][m-1][6]*((y.art_index[m-1] > 7)?y.X[g][r][m-1][7]:0) - (sigma[g][m][6] + lambda_reinit[g][m] + 0.5)*y.X[g][r][m][7];
          out.X[g][r][m][8] += sigma[g][m-1][6]*((y.art_index[m-1] > 7)?y.X[g][r][m-1][8]:0) + 0.5*y.X[g][r][m][7] - (sigma[g][m][6] + lambda_reinit_late[g][m])*y.X[g][r][m][8];

          // ART reinitiation and ART stage progression for re-initiated persons
          out.X[g][r][m][9] += lambda_reinit[g][m]*y.X[g][r][m][7] + lambda_reinit_late[g][m]*y.X[g][r][m][8]  - (sigma[g][m][1] + eta_reinit[g][m][1])*y.X[g][r][m][9];
          out.X[g][r][m][10] += (1-xi[g][m])*sigma[g][m][1]*y.X[g][r][m][9] - (sigma[g][m][2] + eta_reinit[g][m][2])*y.X[g][r][m][10];
          out.X[g][r][m][11] += sigma[g][m][2]*y.X[g][r][m][10] - (sigma[g][m][3] + eta_reinit[g][m][3])*y.X[g][r][m][11];
          out.X[g][r][m][12] += sigma[g][m][3]*y.X[g][r][m][11] - (sigma[g][m][4] + eta_reinit[g][m][4])*y.X[g][r][m][12];
          out.X[g][r][m][13] += xi[g][m]*sigma[g][m][1]*y.X[g][r][m][9] + sigma[g][m][4]*y.X[g][r][m][12] - (sigma[g][m][5] + eta_reinit[g][m][5])*y.X[g][r][m][13];

          //disease progression for second time (permanent) dropouts
          out.X[g][r][m][14] += sigma[g][m-1][6]*((y.art_index[m-1] > 2)?y.X[g][r][m-1][14]:0) - sigma[g][m][6]*y.X[g][r][m][14];
        }
      }

      // drop out from tx during suppressing stage
      for(size_t m = 1; m < DS; m++)
        if(y.art_index[m] > 2){
          out.X[g][r][m][7] += eta[g][m][1]*y.X[g][r][m][2]; // in suppressing stage, go back to same CD4 category
          out.X[g][r][m][14] += eta_reinit[g][m][1]*y.X[g][r][m][9]; // in suppressing stage, go back to same CD4 category
        }

      // drop out from successful treatment < 2 years
      out.X[g][r][3][7] += 0.5*eta[g][5][2]*y.X[g][r][5][3];
      out.X[g][r][4][7] += 0.5*eta[g][5][2]*y.X[g][r][5][3];
      out.X[g][r][2][7] += 0.5*eta[g][4][2]*y.X[g][r][4][3];
      out.X[g][r][3][7] += 0.5*eta[g][4][2]*y.X[g][r][4][3];
      out.X[g][r][2][7] += eta[g][3][2]*y.X[g][r][3][3];
      out.X[g][r][2][7] += eta[g][2][2]*y.X[g][r][2][3];

      out.X[g][r][3][14] += 0.5*eta_reinit[g][5][2]*y.X[g][r][5][10];
      out.X[g][r][4][14] += 0.5*eta_reinit[g][5][2]*y.X[g][r][5][10];
      out.X[g][r][2][14] += 0.5*eta_reinit[g][4][2]*y.X[g][r][4][10];
      out.X[g][r][3][14] += 0.5*eta_reinit[g][4][2]*y.X[g][r][4][10];
      out.X[g][r][2][14] += eta_reinit[g][3][2]*y.X[g][r][3][10];
      out.X[g][r][2][14] += eta_reinit[g][2][2]*y.X[g][r][2][10];

      // drop out from successful treatment > 2 years
      out.X[g][r][2][7] += eta[g][2][3]*y.X[g][r][2][4] + eta[g][3][3]*y.X[g][r][3][4] + eta[g][4][3]*y.X[g][r][4][4];
      out.X[g][r][3][7] += eta[g][5][3]*y.X[g][r][5][4];

      out.X[g][r][2][14] += eta[g][2][3]*y.X[g][r][2][11] + eta[g][3][3]*y.X[g][r][3][11] + eta[g][4][3]*y.X[g][r][4][11];
      out.X[g][r][3][14] += eta[g][5][3]*y.X[g][r][5][11];
    }

  return out;
}
