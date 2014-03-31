// AUTORIGHTS
// -------------------------------------------------------
// Copyright (C) 2013-2014 Menglong Zhu, Nikolay Atanasov
//                         Samarth Brahmbhatt
// 
// This file is part of the Active Deformable Part Models 
// code (http://cis.upenn.edu/~menglong/adpm.html)
// and is available under the terms of an MIT-like license
// provided in COPYING. Please retain this notice and
// COPYING if you use this file (or a portion of it) in
// your project.
// -------------------------------------------------------

#include "mex.h"
#include "policy.h"

// handy accessors

// return field from struct a
static inline const mxArray *F(const mxArray *a, const char *field) {
  return mxGetField(a, 0, field);
}

// return pointer to field from struct a
template <typename T>
static inline T *Fpr(const mxArray *a, const char *field) {
  return mxGetPr(F(a, field));
}

// return scalar of field from struct a
template <typename T>
static inline T Fsc(const mxArray *a, const char *field) {
  return mxGetScalar(F(a, field));
}

// return field from struct in cell of struct array a
static inline const mxArray *CF(const mxArray *a, int cell, const char *field) {
  return F(mxGetCell(a, cell), field);
}

void Policy::initpolicy(const mxArray *policy) {
  numcomponents = mxGetNumberOfElements(policy);
  pol_odp = new double**[numcomponents];
  obs_model_p = new double**[numcomponents];
  obs_model_n = new double**[numcomponents];
  interval = new double**[numcomponents];
  MAPp_min = new double[numcomponents];
  MAPp_res = new double[numcomponents];
  num_part = new int[numcomponents];
  pol_odp_rows = new int[numcomponents];
  obs_model_cols = new int[numcomponents];

  for(int c = 0; c < numcomponents; c++) {

    num_part[c] = (int)Fsc<double>(mxGetCell(policy, c), "num_part");
    MAPp_min[c] = Fsc<double>(CF(policy, c, "MAPp"), "min");
    MAPp_res[c] = Fsc<double>(CF(policy, c, "MAPp"), "res");

    pol_odp[c] = new double*[(int)mxGetN(CF(policy, c, "pol_odp"))];
    pol_odp_rows[c] = mxGetM(CF(policy, c, "pol_odp")); 
    for(unsigned int i = 0; i < mxGetN(CF(policy, c, "pol_odp")); i++)
      pol_odp[c][i] = Fpr<double>(mxGetCell(policy, c), "pol_odp") + i*pol_odp_rows[c];

    obs_model_p[c] = new double*[(int)mxGetN(CF(policy, c, "obs_model_p"))];
    obs_model_n[c] = new double*[(int)mxGetN(CF(policy, c, "obs_model_n"))];
    obs_model_cols[c] = mxGetN(CF(policy, c, "obs_model_p"));
    for(int i = 0; i < obs_model_cols[c]; i++) {
      obs_model_p[c][i] = Fpr<double>(mxGetCell(policy, c), "obs_model_p") + i*mxGetM(CF(policy, c, "obs_model_p"));
      obs_model_n[c][i] = Fpr<double>(mxGetCell(policy, c), "obs_model_n") + i*mxGetM(CF(policy, c, "obs_model_n"));
    }

    interval[c] = new double*[2];
    for(int i = 0; i < 2; i++) {
      interval[c][i] = Fpr<double>(mxGetCell(policy, c), "interval") + i*mxGetM(CF(policy, c, "interval"));
    }
  }
}

Policy::~Policy() {
  for(int c = 0; c < numcomponents; c++) {
    delete [] pol_odp[c];
    delete [] obs_model_p[c];
    delete [] obs_model_n[c];
    delete [] interval[c];
  }
  delete [] MAPp_min;
  delete [] MAPp_res;
  delete [] num_part;
  delete [] pol_odp_rows;
  delete [] obs_model_cols;
}
