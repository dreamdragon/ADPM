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

class Policy {
  public:
    double ***pol_odp, ***obs_model_p, ***obs_model_n, ***interval;
    double *MAPp_min, *MAPp_res;
    int *num_part, *pol_odp_rows, *obs_model_cols, numcomponents;
    
    Policy() {};
    Policy(const mxArray *policy) {initpolicy(policy);};
    ~Policy();

    void initpolicy(const mxArray *policy);
};
