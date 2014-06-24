#include <gsl/gsl_matrix.h>

#define nOutDates 21
#define numParam 17

struct modprev {
  double male[nOutDates];
  double female[nOutDates];
  double artcov[nOutDates];
};

struct modout {
  double * X_out;
  double * trans_out;
  double * inc_out;
  double * art_init_out;
};

void SimPrev(const double theta[numParam], struct modprev * out);
void SimMod(const double theta[numParam], const double * outDates, const size_t numOutDates, const double rateARTinit, const size_t stageARTelig, const double dateARTelig, const double fracScr, struct modout * out);
double R0(const double theta[numParam], const double r0time, double * eigenvec);
double epidemicGrowthRate(const double theta[numParam], const double r0time, double * eigenvec);
void createNGM(const double theta[numParam], const double r0time, gsl_matrix * Fmat, gsl_matrix * Vmat);
void SetParameters(const double theta[numParam], struct parameters * param);
