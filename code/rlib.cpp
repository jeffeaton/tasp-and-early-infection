#include <R.h>
#include <Rinternals.h>

#include <gsl/gsl_matrix.h>

#include "simprev.h"
#include "parameters.h"

extern "C" {

  void callCreateNGM(double * theta, int * afterCbarChange, double * Fmat, double * Vmat){
    // Note: this does NOT check that theta is valid parameter set, length, etc.


    gsl_matrix_view vwFmat = gsl_matrix_view_array(Fmat, NG*RG*(DS-1), NG*RG*(DS-1));
    gsl_matrix_view vwVmat = gsl_matrix_view_array(Vmat, NG*RG*(DS-1), NG*RG*(DS-1));

    if(*afterCbarChange)
      createNGM(theta, INFINITY, &vwFmat.matrix, &vwVmat.matrix);
    else
      createNGM(theta, -INFINITY, &vwFmat.matrix, &vwVmat.matrix);

    return;
  }

  void callR0(double * theta, int * afterCbarChange, double * retR0, double * retVec){
    // Note: this does NOT check that theta is valid parameter set, length, etc.

    if(*afterCbarChange)
      *retR0 = R0(theta, INFINITY, retVec);
    else
      *retR0 = R0(theta, -INFINITY, retVec);
    return;
  }

  void callEpidemicGrowthRate(double * theta, int * afterCbarChange, double * ret_r, double * retVec){
    // Note: this does NOT check that theta is valid parameter set, length, etc.

    if(*afterCbarChange)
      *ret_r = epidemicGrowthRate(theta, INFINITY, retVec);
    else
      *ret_r = epidemicGrowthRate(theta, -INFINITY, retVec);
    return;
  }

  void extractPi(double * theta, double * piVec){

    struct parameters param;
    SetParameters(theta, &param);
    for(size_t g = 0; g < NG; g++)
      for(size_t r = 0; r < RG; r++)
        piVec[g + r*NG] = param.pi[g][r];

    return;
  }

  SEXP callSimMod(SEXP s_theta, SEXP s_outDates, SEXP s_rateARTinit, SEXP s_stageARTelig, SEXP s_dateARTelig, SEXP s_fracScr){

    PROTECT(s_theta = coerceVector(s_theta, REALSXP));
    PROTECT(s_rateARTinit = coerceVector(s_rateARTinit, REALSXP));
    PROTECT(s_stageARTelig = coerceVector(s_stageARTelig, INTSXP));
    PROTECT(s_dateARTelig = coerceVector(s_dateARTelig, REALSXP));
    PROTECT(s_fracScr = coerceVector(s_fracScr, REALSXP));
    PROTECT(s_outDates = coerceVector(s_outDates, REALSXP));

    // prepare the output
    size_t numOutDates = length(s_outDates);

    SEXP s_out, s_out_names, s_X, s_inc, s_trans, s_art_init, s_X_dims, s_inc_dims, s_trans_dims, s_art_init_dims;
    PROTECT(s_X_dims = allocVector(INTSXP, 5));
    PROTECT(s_inc_dims = allocVector(INTSXP, 3));
    PROTECT(s_trans_dims = allocVector(INTSXP, 5));
    PROTECT(s_art_init_dims = allocVector(INTSXP, 4));
    INTEGER(s_X_dims)[0] = numOutDates; INTEGER(s_X_dims)[1] = NG; INTEGER(s_X_dims)[2] = RG+1; INTEGER(s_X_dims)[3] = DS; INTEGER(s_X_dims)[4] = ART_ST;
    INTEGER(s_inc_dims)[0] = numOutDates; INTEGER(s_inc_dims)[1] = NG; INTEGER(s_inc_dims)[2] = RG;
    INTEGER(s_trans_dims)[0] = numOutDates; INTEGER(s_trans_dims)[1] = NG; INTEGER(s_trans_dims)[2] = RG; INTEGER(s_trans_dims)[3] = DS; INTEGER(s_trans_dims)[4] = ART_ST;
    INTEGER(s_art_init_dims)[0] = numOutDates; INTEGER(s_art_init_dims)[1] = NG; INTEGER(s_art_init_dims)[2] = RG+1; INTEGER(s_art_init_dims)[3] = DS;

    PROTECT(s_X = allocVector(REALSXP, numOutDates*NG*(RG+1)*DS*ART_ST));
    PROTECT(s_inc = allocVector(REALSXP, numOutDates*NG*RG));
    PROTECT(s_trans = allocVector(REALSXP, numOutDates*NG*RG*DS*ART_ST));
    PROTECT(s_art_init = allocVector(REALSXP, numOutDates*NG*(RG+1)*DS));

    dimgets(s_X, s_X_dims);
    dimgets(s_inc, s_inc_dims);
    dimgets(s_trans, s_trans_dims);
    dimgets(s_art_init, s_art_init_dims);

    PROTECT(s_out = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(s_out, 0, s_X);
    SET_VECTOR_ELT(s_out, 1, s_inc);
    SET_VECTOR_ELT(s_out, 2, s_trans);
    SET_VECTOR_ELT(s_out, 3, s_art_init);

    PROTECT(s_out_names = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(s_out_names, 0, mkChar("X"));
    SET_VECTOR_ELT(s_out_names, 1, mkChar("inc"));
    SET_VECTOR_ELT(s_out_names, 2, mkChar("trans"));
    SET_VECTOR_ELT(s_out_names, 3, mkChar("art.init"));
    namesgets(s_out, s_out_names);

    struct modout out;
    double *Xptr = REAL(s_X);
    double *incptr = REAL(s_inc);
    double *transptr = REAL(s_trans);
    double *artptr = REAL(s_art_init);
    out.X_out = Xptr;
    out.inc_out = incptr;
    out.trans_out = transptr;
    out.art_init_out = artptr;

    SimMod(REAL(s_theta), REAL(s_outDates), numOutDates, *REAL(s_rateARTinit), *INTEGER(s_stageARTelig), *REAL(s_dateARTelig), *REAL(s_fracScr), &out);

    UNPROTECT(16);

    return s_out;
  }


  SEXP callSimPrev(SEXP s_theta){

    PROTECT(s_theta = coerceVector(s_theta, REALSXP));

    struct modprev out;
    SimPrev(REAL(s_theta), &out);

    //prepare output
    SEXP s_out, s_out_names, s_prev, s_on_art;
    PROTECT(s_prev = allocMatrix(REALSXP, nOutDates, NG));
    PROTECT(s_on_art = allocVector(REALSXP, nOutDates));

    PROTECT(s_out = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(s_out, 0, s_prev);
    SET_VECTOR_ELT(s_out, 1, s_on_art);

    PROTECT(s_out_names = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(s_out_names, 0, mkChar("prev"));
    SET_VECTOR_ELT(s_out_names, 1, mkChar("on.art"));
    namesgets(s_out, s_out_names);

    for(size_t i = 0; i < nOutDates; i++){
      REAL(s_prev)[i + 0*nOutDates] = out.male[i];
      REAL(s_prev)[i + 1*nOutDates] = out.female[i];
      REAL(s_on_art)[i] = out.artcov[i];
    }

    UNPROTECT(5);

    return s_out;
  }

} // extern "C"
