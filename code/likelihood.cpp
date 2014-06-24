#include <stddef.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>

#include "likelihood.h"
#include "simprev.h"

/////////////////////////////
////  Declare functions  ////
/////////////////////////////

double ll(struct modprev * prev);
double logit(const double x);

////  model parameters  ////

const size_t modSimYears = 21; 

////////////////////////////
////  Prior parameters  ////
////////////////////////////

const double ancBiasPrior[2] = {0, 1};

const double theta_unif_lower[numParam] = {1983.0, 0.05, 0.05, 0.2, 0.2, 0.5,  1.0,  0.0, 0.2, 0.0, 0.0, 0.0, 0.00, 0.0, 1990.0, 2002.0};
const double theta_unif_upper[numParam] = {1988.0, 0.70, 0.70, 0.8, 0.8, 4.0, 70.0, 50.0, 0.8, 1.0, 1.0, 1.0, 0.15, 0.7, 2002.0, 2010.0};

const double priminf_rel_prior_mu = 3.2;
const double priminf_rel_prior_sigma = 0.34;


///////////////////////////////////////
////  Declare data for likelihood  ////
///////////////////////////////////////

const double ancLogitVar[] = {0.0174009254,  0.0084540456,  0.0053964641,  0.0032544019,  0.0016378424,  0.0010016838,  0.0009123278,  0.0016000000,  0.0020186588,  0.0011394100,
			      0.0009205834,  0.0011694585,  0.0007565123,  0.0007092568,  0.0006018407,  0.0006458887,  0.0003913861,  0.0003867049,  0.0003882484};
const double ancLogitPrev[] = {-4.8719780, -4.2914736, -3.6969050, -3.1148216, -2.5022585, -2.1492642, -1.8012415, -1.5856273, -1.2196389, -1.2425065,
			       -1.1254595, -1.1093076, -1.0201407, -0.9494274, -0.8712224, -0.8377921, -0.8905323, -0.8760355, -0.8808581};

const double maleLogitVar[] = {0.010017457, 0.007895279, 0.009358600};
const double maleLogitPrev[] = {-1.918759, -2.021151, -2.030867};
const double femaleLogitVar[] = {0.008607043, 0.003814583, 0.003172518};
const double femaleLogitPrev[] = {-1.536806, -1.373841, -1.306936};
const size_t hsrc_idx[] = {12, 15, 18};

const size_t ancDataYears = 19;
const size_t numHSRCYears = 3;


/////////////////////
////  Functions  ////
/////////////////////

double logit(const double x) 
{
  return log(x / (1.0 - x));
};

double ll(struct modprev * prev)
{

  // ll for ANC data
  double ll_S2 = 0.0, ll_dbar = 0.0, ll_d2bar = 0.0, ll_anc = 0.0;
  for(size_t i = 0; i < ancDataYears; i++){
    double logitPi = logit(prev->female[i]);
    ll_S2 += 1.0/ancLogitVar[i];
    ll_dbar += (ancLogitPrev[i] - logitPi)/ancLogitVar[i];
    ll_d2bar += pow(ancLogitPrev[i] - logitPi, 2.0)/ancLogitVar[i];
    ll_anc += log(ancLogitVar[i]);
  }
  ll_S2 = 1.0/ll_S2;
  ll_dbar *= ll_S2;
  ll_d2bar *= ll_S2;

  ll_anc = log(ll_S2)/2.0 - ll_anc/2.0 + (pow(ll_dbar, 2.0) - ll_d2bar)/ll_S2 + log(gsl_cdf_gaussian_P(ancBiasPrior[1] - ll_dbar, sqrt(ll_S2)) - gsl_cdf_gaussian_P(ancBiasPrior[0] - ll_dbar, sqrt(ll_S2)));

  // ll for HSRC survey data
  double ll_male = 0.0, ll_female = 0.0;
  for(size_t i = 0; i < numHSRCYears; i++){
    ll_male += -log(maleLogitVar[i]) - pow(maleLogitPrev[i] - logit(prev->male[hsrc_idx[i]]), 2.0)/maleLogitVar[i];
    ll_female += -log(femaleLogitVar[i]) - pow(femaleLogitPrev[i] - logit(prev->female[hsrc_idx[i]]), 2.0)/femaleLogitVar[i];
  }
  ll_male /= 2.0;
  ll_female /= 2.0;

  return ll_anc + ll_male + ll_female;
 }; 

double likelihood(const gsl_vector * theta)
{
  if(prior(theta) == 0.0)
    return 0.0;

  double th[numParam];
  for(size_t i = 0; i < numParam; i++)
    th[i] = gsl_vector_get(theta, i);

  if(R0(th, -INFINITY, NULL) < 1.0)
    return 0.0;

  // simulate the model
  struct modprev out;
  SimPrev(th, &out);
  for(size_t i = 0; i < modSimYears; i++)
    if(isnan(out.male[i]) || isnan(out.female[i]) || isnan(out.artcov[i])) 
      return 0;

  // calculate the log-likelihood
  return exp(ll(&out));
}

double prior(const gsl_vector * theta)
{

  for(size_t i = 0; i < numParam - 1; i++)
    if(gsl_vector_get(theta, i) < theta_unif_lower[i] || gsl_vector_get(theta, i) > theta_unif_upper[i])
      return 0.0;

  if(gsl_vector_get(theta, 9) > gsl_vector_get(theta, 10) || gsl_vector_get(theta, 10) > gsl_vector_get(theta, 11))
    return 0.0;

  return gsl_ran_lognormal_pdf(gsl_vector_get(theta, numParam - 1), priminf_rel_prior_mu, priminf_rel_prior_sigma);
}

void sample_prior(gsl_rng * r, size_t numSamples, gsl_matrix * storeSamples)
{

  for(size_t i = 0; i < numSamples; i++){
    for(size_t j = 0; j < (numParam - 1); j++)
      gsl_matrix_set(storeSamples, i, j, gsl_ran_flat(r, theta_unif_lower[j], theta_unif_upper[j]));

    gsl_matrix_set(storeSamples, i, numParam - 1, gsl_ran_lognormal(r, priminf_rel_prior_mu, priminf_rel_prior_sigma));
    
    gsl_matrix_set(storeSamples, i, 11, 
		   pow(gsl_matrix_get(storeSamples, i, 11), 1.0/3));
    gsl_matrix_set(storeSamples, i, 10, 
		   pow(gsl_matrix_get(storeSamples, i, 10), 1.0/2) * 
		   gsl_matrix_get(storeSamples, i, 11));
    gsl_matrix_set(storeSamples, i, 9, 
		   gsl_matrix_get(storeSamples, i, 9) * 
		   gsl_matrix_get(storeSamples, i, 10));
  }

  return;
}
