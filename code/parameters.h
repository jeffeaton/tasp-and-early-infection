// parameters
const int NG = 2; // Number of genders
const int RG = 3; // Number of sexual risk groups
const int DS = 6; // Number of disease stages (including uninfected)
const int ART_ST = 15; // Number of stages of ART (including not on ART and screened)

#ifndef fullmod
  #define fullmod 1 // set to 0 if don't need to output count of incident or transmitted infections (model calibration)
#endif

// ART_ST = 0: no treatment
// ART_ST = 1: screened, no treatment
// ART_ST = 2: initiating treatment
// ART_ST = 3: successful treatment, first 2 years
// ART_ST = 4: successful treatment, 2+ years
// ART_ST = 5: treatment failing
// ART_ST = 6: very sick, death
// ART_ST = 7: dropped out from treatment, early dropout (first 2 years)
// ART_ST = 8: dropped out from treatment, long-term dropout (2+ years)
// ART_ST = 9: re-initiating treatment
// ART_ST = 10: successful treatment, first 2 years
// ART_ST = 11: successful treatment, 2+ years
// ART_ST = 12: treatment failing
// ART_ST = 13: very sick, death
// ART_ST = 14: dropped out, permanently off treatment

struct parameters {

  //  Declare varying parameters
  double c[NG][RG];
  double cbar[NG];
  double epsilon;
  double nL[RG][RG];
  double prop_new[NG][RG];
  double psi[NG][RG][RG];
  double lambda[NG][DS];
  double lambda_scr[NG][DS];
  double frac_scr[NG];
  double pi[NG][RG];
  double omega[NG][RG];
  double startYear;
  double delta[NG][RG][DS];
  double beta[NG][DS][ART_ST];

};

// fixed model parameters (global)
const double beta_init[NG][DS][ART_ST] = {{{0.00000000, 0.00000000, 0.00000000, 0.0000000, 0.0000000, 0.0000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.0000000, 0.0000000, 0.0000000, 0.00000000, 0.00000000},
					   {-1.0      , -1.0      , -1.0      , 0.0057714, 0.0057714, 0.2706632, 0.05107432, -1.0      , -1.0      , -1.0      , 0.0057714, 0.0057714, 0.2706632, 0.05107432, -1.0      },
					   {0.04395779, 0.04395779, 0.02197890, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.04395779, 0.04395779, 0.02197890, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.04395779},
					   {0.07214250, 0.07214250, 0.03607125, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.07214250, 0.07214250, 0.03607125, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.07214250},
					   {0.27066315, 0.27066315, 0.13533158, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.27066315, 0.27066315, 0.13533158, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.27066315},
					   {0.05107432, 0.05107432, 0.02553716, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.05107432, 0.05107432, 0.02553716, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.05107432}},
					  {{0.00000000, 0.00000000, 0.00000000, 0.0000000, 0.0000000, 0.0000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.0000000, 0.0000000, 0.0000000, 0.00000000, 0.00000000},
					   {-1.0      , -1.0      , -1.0      , 0.0057714, 0.0057714, 0.2706632, 0.05107432, -1.0      , -1.0      , -1.0      , 0.0057714, 0.0057714, 0.2706632, 0.05107432, -1.0      },
					   {0.04395779, 0.04395779, 0.02197890, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.04395779, 0.04395779, 0.02197890, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.04395779},
					   {0.07214250, 0.07214250, 0.03607125, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.07214250, 0.07214250, 0.03607125, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.07214250},
					   {0.27066315, 0.27066315, 0.13533158, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.27066315, 0.27066315, 0.13533158, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.27066315},
					   {0.05107432, 0.05107432, 0.02553716, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.05107432, 0.05107432, 0.02553716, 0.0057714, 0.0057714, 0.2706632, 0.05107432, 0.05107432}}};

const double alpha = 0.02261967;
const double nu = 1.0/35;
const double mu = 0.087318;
const double theta_g = 0.5;

const double n_init[NG] = {1000, 1000};

const double sigma[NG][DS][ART_ST-8] = {{{0.0000000, 0.0, 0.0000000, 0.00000000, 0.0000000, 0.000000, 0.0000000},
                                         {4.1666667, 4.0, 0.5714286, 0.02191630, 0.4792898, 1.947118, 8.3333333},
                                         {0.2192982, 4.0, 0.5714286, 0.02191630, 0.4792898, 1.947118, 0.4385965},
                                         {0.2173913, 4.0, 0.5714286, 0.02581109, 0.4792898, 1.947118, 0.4347826},
                                         {0.2396449, 4.0, 0.5714286, 0.03030303, 0.4792898, 1.947118, 0.4792898},
                                         {0.9735592, 4.0, 0.5714286, 0.03030303, 0.4792898, 1.947118, 1.9471184}},
                                        {{0.0000000, 0.0, 0.0000000, 0.00000000, 0.0000000, 0.000000, 0.0000000},
                                         {4.1666667, 4.0, 0.5714286, 0.02191630, 0.4792898, 1.947118, 8.3333333},
                                         {0.2192982, 4.0, 0.5714286, 0.02191630, 0.4792898, 1.947118, 0.4385965},
                                         {0.2173913, 4.0, 0.5714286, 0.02581109, 0.4792898, 1.947118, 0.4347826},
                                         {0.2396449, 4.0, 0.5714286, 0.03030303, 0.4792898, 1.947118, 0.4792898},
                                         {0.9735592, 4.0, 0.5714286, 0.03030303, 0.4792898, 1.947118, 1.9471184}}};

const double xi[NG][DS] = {{0.0, 0.0, 0, 0.025, 0.06717725, 0.1890904},
                           {0.0, 0.0, 0, 0.025, 0.06717725, 0.1890904}};

const double eta[NG][DS][ART_ST-9] = {{{0.0, 0.000, 0.000, 0.000, 0.0, 0.0},
                                       {0.0, 0.168, 0.168, 0.088, 0.0, 0.0},
                                       {0.0, 0.168, 0.168, 0.088, 0.0, 0.0},
                                       {0.0, 0.168, 0.168, 0.088, 0.0, 0.0},
                                       {0.0, 0.156, 0.156, 0.088, 0.0, 0.0},
                                       {0.0, 0.120, 0.120, 0.088, 0.0, 0.0}},
                                      {{0.0, 0.000, 0.000, 0.000, 0.0, 0.0},
                                       {0.0, 0.168, 0.168, 0.088, 0.0, 0.0},
                                       {0.0, 0.168, 0.168, 0.088, 0.0, 0.0},
                                       {0.0, 0.168, 0.168, 0.088, 0.0, 0.0},
                                       {0.0, 0.156, 0.156, 0.088, 0.0, 0.0},
                                       {0.0, 0.120, 0.120, 0.088, 0.0, 0.0}}};

const double eta_reinit[NG][DS][ART_ST-9] = {{{0.0, 0.00, 0.00, 0.00, 0.0, 0.0},
                                              {0.0, 0.06, 0.06, 0.06, 0.0, 0.0},
                                              {0.0, 0.06, 0.06, 0.06, 0.0, 0.0},
                                              {0.0, 0.06, 0.06, 0.06, 0.0, 0.0},
                                              {0.0, 0.06, 0.06, 0.06, 0.0, 0.0},
                                              {0.0, 0.06, 0.06, 0.06, 0.0, 0.0}},
                                             {{0.0, 0.00, 0.00, 0.00, 0.0, 0.0},
                                              {0.0, 0.06, 0.06, 0.06, 0.0, 0.0},
                                              {0.0, 0.06, 0.06, 0.06, 0.0, 0.0},
                                              {0.0, 0.06, 0.06, 0.06, 0.0, 0.0},
                                              {0.0, 0.06, 0.06, 0.06, 0.0, 0.0},
                                              {0.0, 0.06, 0.06, 0.06, 0.0, 0.0}}};

const double lambda_reinit[NG][DS] = {{0.0, 0.0, 0.0, 0.04830918, 0.1597633, 2.920678},
                                      {0.0, 0.0, 0.0, 0.04830918, 0.1597633, 2.920678}};
const double lambda_reinit_late[NG][DS] = {{0.0, 0.0, 0.0, 0.04830918, 0.1597633, 2.920678},
                                           {0.0, 0.0, 0.0, 0.04830918, 0.1597633, 2.920678}};