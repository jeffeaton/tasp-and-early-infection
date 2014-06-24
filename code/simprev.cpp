#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include "parameters.h"
#include "states.h"
#include "model.h"
#include "simprev.h"


/////////////////////////////////////////
////  Define fixed model parameters  ////
/////////////////////////////////////////

const unsigned int breaks = 10;
const double dt = 1.0/breaks;

/////////////////////////////
////  Declare functions  ////
/////////////////////////////

void SetCbarLogistic(const double theta[numParam], const double currTime, struct parameters * param);
void SetDelta(const double theta[numParam], struct parameters * param);
void SetARTRates(const double ARTRate, const size_t ARTEligStage, struct parameters * param);
void SetFracScr(const double fracScr, struct parameters * param);
void InitializePopulation(states * y, const struct parameters * param);
void PrepareART(states * current);
void InitiateART(double art_frac, states * y);

void RecordPrevOutput(const states * mod, struct modprev * out, const size_t outIdx);
void RecordFullOutput(states * y, const size_t outIdx, const size_t numOutDates, struct modout * out);
void PrintParameters(const struct parameters * param);


////////////////////////////
////  Define functions  ////
////////////////////////////

void SimPrev(const double theta[numParam], struct modprev * out)
{

  // fixed model parameters for calibration simulation
  const size_t nARTdates = 60;
  const double art_fraction[nARTdates] = {0.000183501693954, 0.000367003387909, 0.000550505081863, 0.000734006775818, 0.000917508469772, 0.001101010163726, 0.001284511857681, 0.001468013551635, 0.001651515245589, 0.001835016939544,
                                          0.002135802936128, 0.002436588932711, 0.002737374929295, 0.003038160925879, 0.003338946922462, 0.003639732919046, 0.003940518915630, 0.004241304912214, 0.004542090908797, 0.004842876905381,
                                          0.005238861672595, 0.005634846439809, 0.006030831207024, 0.006426815974238, 0.006822800741452, 0.007218785508666, 0.007614770275881, 0.008010755043095, 0.008406739810309, 0.008802724577523,
                                          0.009401255559511, 0.009999786541500, 0.010598317523488, 0.011196848505476, 0.011795379487464, 0.012393910469453, 0.012992441451441, 0.013590972433429, 0.014189503415417, 0.014788034397405,
                                          0.015650731438506, 0.016513428479608, 0.017376125520709, 0.018238822561810, 0.019101519602911, 0.019964216644012, 0.020826913685113, 0.021689610726214, 0.022552307767315, 0.023415004808416,
                                          0.023978065526726, 0.024541126245035, 0.025104186963344, 0.025667247681653, 0.026230308399962, 0.026793369118272, 0.027356429836581, 0.027919490554890, 0.028482551273199, 0.029045611991509};
  const double art_dates[nARTdates] = {2004.6, 2004.7, 2004.8, 2004.9, 2005.0, 2005.1, 2005.2, 2005.3, 2005.4, 2005.5, 2005.6, 2005.7, 2005.8, 2005.9, 2006.0, 2006.1, 2006.2, 2006.3, 2006.4, 2006.5,
                                       2006.6, 2006.7, 2006.8, 2006.9, 2007.0, 2007.1, 2007.2, 2007.3, 2007.4, 2007.5, 2007.6, 2007.7, 2007.8, 2007.9, 2008.0, 2008.1, 2008.2, 2008.3, 2008.4, 2008.5,
                                       2008.6, 2008.7, 2008.8, 2008.9, 2009.0, 2009.1, 2009.2, 2009.3, 2009.4, 2009.5, 2009.6, 2009.7, 2009.8, 2009.9, 2010.0, 2010.1, 2010.2, 2010.3, 2010.4, 2010.5};
  const double calibOutDates[nOutDates] = {1990.5, 1991.5, 1992.5, 1993.5, 1994.5, 1995.5, 1996.5, 1997.5, 1998.5, 1999.5, 2000.5, 2001.5, 2002.5, 2003.5, 2004.5, 2005.5, 2006.5, 2007.5, 2008.5, 2009.5, 2010.5};
  const double stopYear = 2010.6;

  // set parameter values
  struct parameters param;
  SetParameters(theta, &param);
  SetDelta(theta, &param);

  // set simulation time steps
  size_t startTS = param.startYear * breaks;
  size_t stopTS =  stopYear * breaks;

  // determine the time steps to output
  unsigned int currOutDate = 0, *out_ts = new unsigned int [nOutDates];
  for(size_t i = 0; i < nOutDates; i++)
    out_ts[i] = calibOutDates[i]*breaks; // add breaks/2 to get mid-year results

  // determine ART scale up dates
  size_t currARTdate = 0, *ART_ts = new size_t [nARTdates];
  for(size_t i = 0; i < nARTdates; i++)
    ART_ts[i] = art_dates[i] * breaks;

  // initialise the population
  states current;
  InitializePopulation(&current, &param);


  // numerically solve the model
  for(size_t ts = startTS; ts < stopTS; ts++){

    // set contact rate
    SetCbarLogistic(theta, param.startYear + dt*(ts - startTS), &param);

    current = euler(current, dt, &param);

    // do ART initiation
    if(currARTdate < nARTdates && ART_ts[currARTdate] == ts){
      if(currARTdate == 0)
        PrepareART(&current);
      InitiateART(art_fraction[currARTdate], &current);
      currARTdate++;
    }

    // record the outputs
    if(currOutDate < nOutDates && ts == out_ts[currOutDate]){
      RecordPrevOutput(&current, out, currOutDate);
      currOutDate++;
    }
  } // end running the model

  delete[] out_ts;
  delete[] ART_ts;

  return;
}


void SimMod(const double theta[numParam], const double * outDates, const size_t numOutDates, const double rateARTinit, const size_t stageARTelig, const double dateARTelig, const double fracScr, struct modout * out)
{

  // set the parameters
  struct parameters param;
  SetParameters(theta, &param);
  SetDelta(theta, &param);

  // determine the time steps to output
  size_t startTS = param.startYear * breaks, stopTS = (outDates[numOutDates - 1] + dt) * breaks;
  unsigned int currOutDate = 0, * out_ts = (unsigned int *) malloc(numOutDates * sizeof(unsigned int));
  for(size_t i = 0; i < numOutDates; i++)
    out_ts[i] = outDates[i] * breaks;

  // determine ART introduction time step and set fraction who will be 
  size_t ARTts = dateARTelig * breaks;
  SetFracScr(fracScr, &param);

  // initialise the population 
  states current;
  InitializePopulation(&current, &param);


  // numerically solve the model
  for(size_t ts = startTS; ts < stopTS; ts++){

    // set contact rate
    SetCbarLogistic(theta, param.startYear + dt*(ts - startTS), &param);

    // set ART eligibility
    if(ts == ARTts){
      PrepareART(&current);
      SetARTRates(rateARTinit, stageARTelig, &param);
    }

    current = euler(current, dt, &param);

    // record the outputs 
    if(currOutDate < numOutDates && ts == out_ts[currOutDate])
      RecordFullOutput(&current, currOutDate++, numOutDates, out);

  } // end running the model

  free(out_ts);

  return;
}



size_t get_NGM_idx(size_t g, size_t r, size_t m)
{
  // gender g
  // risk group r
  // disease stage m

  return g*(DS-1)*RG + (m-1)*RG + r;
}



void createNGM(const double theta[numParam], const double r0time, gsl_matrix * Fmat, gsl_matrix * Vmat)
{
  struct parameters param;
  SetParameters(theta, &param);
  SetCbarLogistic(theta, r0time, &param);

  // transmission matrix (T)
  for(size_t g = 0; g < NG; g++){
    size_t gprime = (g == 0) ? 1 : 0;
    for(size_t r = 0; r < RG; r++){
      for(size_t rprime = 0; rprime < RG; rprime++){
        for(size_t mprime = 1; mprime < DS; mprime++){
          gsl_matrix_set(Fmat, get_NGM_idx(g, r, 1), get_NGM_idx(gprime, rprime, mprime),
                         param.pi[g][r] / param.pi[gprime][rprime] * param.c[g][r] *
                         (param.epsilon * ((r==rprime)?1:0) + (1.0 - param.epsilon) * (param.c[gprime][rprime]/param.cbar[gprime]) * param.pi[gprime][rprime]) *
                         (1.0- exp(-param.beta[gprime][mprime][0]*param.nL[rprime][r])));
        }
      }
    }
  }
  
  // state transition matrix (Sigma)
  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r < RG; r++)
      for(size_t m = 1; m < DS; m++){
        gsl_matrix_set(Vmat, get_NGM_idx(g, r, m), get_NGM_idx(g, r, m), nu);
        for(size_t rprime = 0; rprime < RG; rprime++){
          gsl_matrix_set(Vmat, get_NGM_idx(g, r, m), get_NGM_idx(g, rprime, m),
                         gsl_matrix_get(Vmat, get_NGM_idx(g, r, m), get_NGM_idx(g, rprime, m)) - param.psi[g][rprime][r]);
          gsl_matrix_set(Vmat, get_NGM_idx(g, r, m), get_NGM_idx(g, r, m),
                         gsl_matrix_get(Vmat, get_NGM_idx(g, r, m), get_NGM_idx(g, r, m)) + param.psi[g][r][rprime]);
        }
        gsl_matrix_set(Vmat, get_NGM_idx(g, r, m), get_NGM_idx(g, r, m),
                       gsl_matrix_get(Vmat, get_NGM_idx(g, r, m), get_NGM_idx(g, r, m)) + sigma[g][m][0]);
        if(m < (DS-1))
          gsl_matrix_set(Vmat, get_NGM_idx(g, r, m+1), get_NGM_idx(g, r, m),
                         gsl_matrix_get(Vmat, get_NGM_idx(g, r, m+1), get_NGM_idx(g, r, m)) - sigma[g][m][0]);
      }

  return;
}

double R0(const double theta[numParam], const double r0time, double * eigenvec)
{

  gsl_matrix * Fmat = gsl_matrix_calloc(NG*(DS-1)*RG, NG*(DS-1)*RG);
  gsl_matrix * Vmat = gsl_matrix_calloc(NG*(DS-1)*RG, NG*(DS-1)*RG);
  gsl_matrix * VmatInv = gsl_matrix_calloc(NG*(DS-1)*RG, NG*(DS-1)*RG);
  gsl_matrix * ngm = gsl_matrix_calloc(NG*(DS-1)*RG, NG*(DS-1)*RG);

  createNGM(theta, r0time, Fmat, Vmat);

  gsl_permutation * p = gsl_permutation_alloc(NG*(DS-1)*RG);
  int s;
  gsl_linalg_LU_decomp(Vmat, p, &s);
  gsl_linalg_LU_invert(Vmat, p, VmatInv);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Fmat, VmatInv, 0.0, ngm);

  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(NG*(DS-1)*RG);
  gsl_vector_complex * eval = gsl_vector_complex_alloc(NG*(DS-1)*RG);
  gsl_matrix_complex * evec = gsl_matrix_complex_alloc(NG*(DS-1)*RG, NG*(DS-1)*RG);

  gsl_set_error_handler_off();
  gsl_eigen_nonsymmv(ngm, eval, evec, w);

  size_t r0_idx = 0;
  double r0 = 0.0;
  for(size_t i = 0; i < NG*(DS-1)*RG; i++){
    if(GSL_REAL(gsl_vector_complex_get(eval, i)) > r0){
      r0_idx = i;
      r0 = GSL_REAL(gsl_vector_complex_get(eval, i));
    }
  }

  if(eigenvec != NULL){
    for(size_t i = 0; i < NG*(DS-1)*RG; i++)
      eigenvec[i] = GSL_REAL(gsl_matrix_complex_get(evec, i, r0_idx));
  }

  gsl_matrix_free(Fmat);
  gsl_matrix_free(Vmat);
  gsl_matrix_free(VmatInv);
  gsl_matrix_free(ngm);
  gsl_permutation_free(p);
  gsl_eigen_nonsymmv_free(w);
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);

  return r0;
}

double epidemicGrowthRate(const double theta[numParam], const double r0time, double * eigenvec)
{

  gsl_matrix * Fmat = gsl_matrix_calloc(NG*(DS-1)*RG, NG*(DS-1)*RG);
  gsl_matrix * Vmat = gsl_matrix_calloc(NG*(DS-1)*RG, NG*(DS-1)*RG);

  createNGM(theta, r0time, Fmat, Vmat);
  gsl_matrix_sub(Fmat, Vmat);

  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(NG*(DS-1)*RG);
  gsl_vector_complex * eval = gsl_vector_complex_alloc(NG*(DS-1)*RG);
  gsl_matrix_complex * evec = gsl_matrix_complex_alloc(NG*(DS-1)*RG, NG*(DS-1)*RG);

  gsl_set_error_handler_off();
  gsl_eigen_nonsymmv(Fmat, eval, evec, w);

  size_t growth_rate_idx = 0;
  double growth_rate = -INFINITY;
  for(size_t i = 0; i < NG*(DS-1)*RG; i++){
    if(GSL_REAL(gsl_vector_complex_get(eval, i)) > growth_rate){
      growth_rate_idx = i;
      growth_rate = GSL_REAL(gsl_vector_complex_get(eval, i));
    }
  }

  if(eigenvec != NULL){
    for(size_t i = 0; i < NG*(DS-1)*RG; i++)
      eigenvec[i] = GSL_REAL(gsl_matrix_complex_get(evec, i, growth_rate_idx));
  }

  gsl_matrix_free(Fmat);
  gsl_matrix_free(Vmat);
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
  gsl_eigen_nonsymmv_free(w);

  return growth_rate;
}


void SetParameters(const double theta[numParam], struct parameters * param)
{

  // start time
  param->startYear = theta[0];

  // prop_new
  param->prop_new[0][0] = theta[1]*theta[3];
  param->prop_new[1][0] = theta[2]*theta[4];
  param->prop_new[0][1] = theta[1] - theta[1]*theta[3];
  param->prop_new[1][1] = theta[2] - theta[2]*theta[4];
  param->prop_new[0][2] = 1.0 - theta[1];
  param->prop_new[1][2] = 1.0 - theta[2];

  // psi
  for(size_t i = 0; i < NG; i++)
    for(size_t j = 0; j < RG; j++)
      for(size_t k = 0; k < RG; k++)
        if(k > j)
          param->psi[i][j][k] = theta[12];
        else
          param->psi[i][j][k] = 0.0;

  // pi
  double tmpSum;
  param->pi[0][0] = (alpha + nu + theta[12]) / ((param->prop_new[0][1]/param->prop_new[0][0])*(alpha + nu + 2*theta[12]) + theta[12]);
  param->pi[0][1] = 1;
  param->pi[0][2] = ((param->prop_new[0][2] + theta[12]/(alpha+nu)) * (alpha + nu + 2*theta[12])) / (param->prop_new[0][1] * (alpha + nu + theta[12]) + (1.0 - param->prop_new[0][2])*theta[12]);
  tmpSum = param->pi[0][0] + param->pi[0][1] + param->pi[0][2];
  param->pi[0][0] /= tmpSum;
  param->pi[0][1] /= tmpSum;
  param->pi[0][2] /= tmpSum;

  param->pi[1][0] = (alpha + nu + theta[12]) / ((param->prop_new[1][1]/param->prop_new[1][0])*(alpha + nu + 2*theta[12]) + theta[12]);
  param->pi[1][1] = 1;
  param->pi[1][2] = ((param->prop_new[1][2] + theta[12]/(alpha+nu)) * (alpha + nu + 2*theta[12])) / (param->prop_new[1][1] * (alpha + nu + theta[12]) + (1.0 - param->prop_new[1][2])*theta[12]);
  tmpSum = param->pi[1][0] + param->pi[1][1] + param->pi[1][2];
  param->pi[1][0] /= tmpSum;
  param->pi[1][1] /= tmpSum;
  param->pi[1][2] /= tmpSum;

  // omega
  param->omega[1][0] = theta[6] + theta[7];
  param->omega[1][1] = theta[6];
  param->omega[1][2] = 1.0;
  param->omega[0][0] = (param->omega[1][0]*param->pi[1][0]/param->pi[0][0]) / (param->omega[1][2]*param->pi[1][2] / param->pi[0][2]);
  param->omega[0][1] = (param->omega[1][1]*param->pi[1][1]/param->pi[0][1]) / (param->omega[1][2]*param->pi[1][2] / param->pi[0][2]);
  param->omega[0][2] = 1.0;

  // nL (same as kappa)
  for(size_t r1 = 0; r1 < RG; r1++)
    for(size_t r2 = 0; r2 < RG; r2++)
      if(r1 == 0 | r2 == 0)
        param->nL[r1][r2] = theta[9];
      else if(r1 == 1 | r2 == 1)
        param->nL[r1][r2] = theta[10];
      else
        param->nL[r1][r2] = theta[11];

  // epsilon
  param->epsilon = theta[8];

  // lambda
  for(size_t i = 0; i < NG; i++)
    for(size_t j = 0; j < DS; j++){
      param->lambda[i][j] = 0.0;
      param->lambda_scr[i][j] = 0.0;
    }

  // frac_scr
  for(size_t g = 0; g < NG; g++)
    param->frac_scr[g] = 0.0;

  // beta
  for(size_t g = 0; g < NG; g++)
    for(size_t m = 0; m < DS; m++)
      for(size_t u = 0; u < ART_ST; u++)
	param->beta[g][m][u] = beta_init[g][m][u];

  for(size_t g = 0; g < NG; g++){
    param->beta[g][1][0] = param->beta[g][2][0]*theta[16];
    param->beta[g][1][1] = param->beta[g][2][1]*theta[16];
    param->beta[g][1][2] = param->beta[g][2][2]*theta[16];
    param->beta[g][1][7] = param->beta[g][2][7]*theta[16];
    param->beta[g][1][8] = param->beta[g][2][8]*theta[16];
    param->beta[g][1][9] = param->beta[g][2][9]*theta[16];
    param->beta[g][1][14] = param->beta[g][2][14]*theta[16];
  }    

  return;
}

void SetCbarLogistic(const double theta[numParam], const double currTime, struct parameters * param)
{

  // set cbar
  double cx = theta[13];
  double cbar = theta[5];
  double cx_start = theta[14];
  double cx_dur = (theta[15] - theta[14]);
  double cbar_t = (1 - cx)*theta[5] + theta[5]*cx/(1+exp((currTime - (cx_start + cx_dur/2))*10/cx_dur));
  param->cbar[0] = param->cbar[1] = cbar_t;

  // calculate within risk group contact rate
  double pi_omega[NG] = {0, 0};
  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r < RG; r++)
      pi_omega[g] += param->pi[g][r]*param->omega[g][r];

  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r < RG; r++)
      param->c[g][r] = param->omega[g][r]*cbar_t/pi_omega[g];

  return;
}

void SetARTRates(const double ARTRate, const size_t ARTEligStage, struct parameters * param)
{
  for(size_t g = 0; g < NG; g++)
    for(size_t m = ARTEligStage; m < DS; m++){
      param->lambda[g][m] = ARTRate;
      param->lambda_scr[g][m] = 1.0;
    }

  return;
}

void SetFracScr(const double fracScr, struct parameters * param)
{
  for(size_t g = 0; g < NG; g++)
    param->frac_scr[g] = fracScr;

  return;
}

void PrepareART(states * current)
{

  for(size_t mup = 1; mup < DS; mup++){  // hard coded to change art_index for stages 1+
    if(current->art_index[mup] < ART_ST){
      for(size_t g = 0; g < NG; g++)
        for(size_t r = 0; r <= RG; r++)
          for(size_t u = current->art_index[mup]; u < ART_ST; u++){
            current->X[g][r][mup][u] = 0;
           #if fullmod
            if(r < RG)
              current->trans[g][r][mup][u] = 0;
           #endif
          }
      current->art_index[mup] = ART_ST;
    }
  }

  return;
}

void InitiateART(double art_frac, states * y)
{
#define ART_ELIG 4

  const double relative_ART[NG][DS] = {{0.0, 0.0, 0.0, 0.0, 1.0/8, 1.0},
                                       {0.0, 0.0, 0.0, 0.0, 1.4/8, 1.4}};

  // Current ART coverage
  double Xart = 0.0, XnoTx = 0.0, Xtot = 0.0, weighted_elig = 0.0;
  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r <= RG; r++)
      for(size_t m = 0; m < DS; m++){
        if(m >= ART_ELIG){
          XnoTx += y->X[g][r][m][0];
          weighted_elig += relative_ART[g][m]*y->X[g][r][m][0];
        }
        for(size_t u = 0; u < y->art_index[m]; u++){
          Xtot += y->X[g][r][m][u];
          if((u >= 2 && u <= 6) || (u >= 9 && u <= 13))
            Xart += y->X[g][r][m][u];
        }
      }

  double curr_cov = Xart/Xtot;
  double art_init = art_frac*Xtot - Xart;
  if(art_init > XnoTx)
    art_init = XnoTx;

  double numinit = 0.0;
  if(art_init > 0){

    for(size_t g = 0; g < NG; g++)
      for(size_t r = 0; r <= RG; r++)
        for(size_t m = ART_ELIG; m < DS; m++){
          double tmp_init = art_init * relative_ART[g][m] * y->X[g][r][m][0]/weighted_elig;
          if(tmp_init > y->X[g][r][m][0])
            tmp_init = y->X[g][r][m][0];
          y->X[g][r][m][0] -= tmp_init;
          y->X[g][r][m][2] += tmp_init;
          numinit += tmp_init;
        }
  }

  return;
}

void SetDelta(const double theta[numParam], struct parameters * param)
{

  double eigenvec[NG*(DS-1)*RG];
  epidemicGrowthRate(theta, -INFINITY, eigenvec);
  double sumEvec = 0.0;

  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r < RG; r++)
      for(size_t m = 0; m < DS-1; m++)
	sumEvec += eigenvec[g*(DS-1)*RG + m*RG + r];

  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r < RG; r++){
      param->delta[g][r][0] = 0.0;
      for(size_t m = 1; m < DS; m++){
	param->delta[g][r][m] = 0.0005*eigenvec[g*(DS-1)*RG + (m-1)*RG + r] / sumEvec;
      }
    }

  return;
}

void InitializePopulation(states * y, const struct parameters * param)
{

  // initialise population
  states current;
  for(size_t g = 0; g < NG; g++)
    for(size_t m = 0; m < DS; m++)
      for(size_t u = 0; u < ART_ST; u++){
        for(size_t r = 0; r < RG; r++)
          y->X[g][r][m][u] = 0;
        y->X[g][RG][m][u] = 0;
      }

  for(size_t g = 0; g < NG; g++){
    for(size_t r = 0; r < RG; r++){
      double rg_inf = 0.0;
      for(size_t m = 1; m < DS; m++){
	y->X[g][r][m][0] = n_init[g] * param->delta[g][r][m];
	rg_inf += y->X[g][r][m][0];
      }
      y->X[g][r][0][0] = n_init[g]*param->pi[g][r] - rg_inf;
    }
    y->X[g][RG][0][0] = n_init[g]*nu/(mu+alpha);
  }

  for(size_t m = 0; m < DS; m++)
    if(param->frac_scr[0] > 0.0 || param->frac_scr[1] > 0.0)
      y->art_index[m] = 2;
    else
      y->art_index[m] = 1;
  
  #if fullmod
  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r < RG; r++){
      y->inc[g][r] = 0;
      for(size_t m = 0; m < DS; m++)
	for(size_t u = 0; u < ART_ST; u++)
	  y->trans[g][r][m][u] = 0;
    }
  #endif

  return;
}

void RecordPrevOutput(const states * mod, struct modprev * out, const size_t outIdx)
{

  double Xg, Xneg, XonART = 0, Xtot = 0;
  for(size_t g = 0; g < NG; g++){
    Xg = Xneg = 0;
    for(size_t r = 0; r < RG; r++){
      Xneg += mod->X[g][r][0][0];
      Xg += mod->X[g][r][0][0];
      for(size_t m = 1; m < DS; m++){
        Xg += mod->X[g][r][m][0];
        for(size_t u = 1; u < mod->art_index[m]; u++){
          if((u >= 2 && u <=6) || (u >= 9 && u <= 13)) // !!!! THIS IS HARDCODED ART STAGE PATTERN !!!!
            XonART += mod->X[g][r][m][u];
          Xg += mod->X[g][r][m][u];
        }
      }
    }
    if(g == 0)
      out->male[outIdx] = (Xg-Xneg)/Xg;
    if(g == 1)
      out->female[outIdx] = (Xg-Xneg)/Xg;
    Xtot += Xg;
  }
  out->artcov[outIdx] = XonART/Xtot;

  return;
}


void RecordFullOutput(states * y, const size_t outIdx, const size_t numOutDates, struct modout * out)
{
  for(size_t g = 0; g < NG; g++)
    for(size_t r = 0; r <= RG; r++){
      for(size_t m = 0; m < DS; m++){
	for(size_t u = 0; u < y->art_index[m]; u++){
	  out->X_out[outIdx + (g + NG*(r + (RG+1)*(m + DS*u)))*numOutDates] = y->X[g][r][m][u];
	  #if fullmod
	  if(r < RG){
	    out->trans_out[outIdx + (g + NG*(r + RG*(m + DS*u)))*numOutDates] = y->trans[g][r][m][u];
	    y->trans[g][r][m][u] = 0;
	  }
	  #endif
	}
	for(size_t u = y->art_index[m]; u < ART_ST; u++){
	  out->X_out[outIdx + (g + NG*(r + (RG+1)*(m + DS*u)))*numOutDates] = 0;
	  #if fullmod
	  if(r < RG)
	    out->trans_out[outIdx + (g + NG*(r + RG*(m + DS*u)))*numOutDates] = 0;
	  #endif
	}
	#if fullmod
	out->art_init_out[outIdx + (g + NG*(r + (RG+1)*m))*numOutDates] = y->art_init[g][r][m];
	y->art_init[g][r][m] = 0;
	#endif
      }
      #if fullmod
      if(r < RG){
	out->inc_out[outIdx + (g + NG*r)*numOutDates] = y->inc[g][r];
	y->inc[g][r] = 0;
      }
      #endif
    }

  return;
}


void PrintParameters(const struct parameters * param)
{

  printf("\npi:\n");
  for(size_t i = 0; i < NG; i++){
    for(size_t j = 0; j < RG; j++)
      printf("%f ", param->pi[i][j]);
    printf("\n");
  }

  printf("\nprop_new:\n");
  for(size_t i = 0; i < NG; i++){
    for(size_t j = 0; j < RG; j++)
      printf("%f ", param->prop_new[i][j]);
    printf("\n");
  }

  printf("\nomega:\n");
  for(size_t i = 0; i < NG; i++){
    for(size_t j = 0; j < RG; j++)
      printf("%f ", param->omega[i][j]);
    printf("\n");
  }

  printf("\nc:\n");
  for(size_t i = 0; i < NG; i++){
    for(size_t j = 0; j < RG; j++)
      printf("%f ", param->c[i][j]);
    printf("\n");
  }

  printf("\nepsilon = %f\n", param->epsilon);

  printf("\npsi:\n");
  for(size_t i = 0; i < NG; i++){
    for(size_t j = 0; j < RG; j++){
      for(size_t k = 0; k < RG; k++)
        printf("%f ", param->psi[i][j][k]);
      printf("\n");
    }
    printf("\n");
  }

  printf("\nlambda:\n");
  for(size_t i = 0; i < NG; i++){
    for(size_t j = 0; j < DS; j++){
      printf("%f ", param->lambda[i][j]);
    }
    printf("\n");
  }

  return;
}
