
data {
  int<lower=0> N;              // nr of patient observations
  int<lower=0> M;              // nr of institutes
  int<lower=0> NT;             // nr of unique failure times
  int<lower=0> obs_t[N];       // observed times
  int<lower=0> ind[N];         // censor status, 0:alive, 1:dead
  int<lower=0> ctr[N];         // institute id
  real Zage[N];                // normalized patient age 
  real Zkps[N];                // normalized patient kps 
  int<lower=0> yot[N];         // year of treatment id
  int<lower=0> ind_t[N];       // t index of observed times
  int<lower=0> t[NT + 1];      // unique failure times + max followup time
  int<lower=0> dt[NT + 1];     // length per time interval
  int<lower=0> j1m;            // first time index > one month
  int<lower=0> j2y;            // first time index > two years
  int<lower=0> nperiods;       // number of periods
  int<lower=0> period[NT + 1]; // period per unique times

}
transformed data {
  int<lower=0, upper=1> Y[N, NT];    // patient i observed at beginning of time interval j, 0:no (after death or censored time), 1:yes  
  int<lower=0, upper=1> dN[N, NT];   // counting process increment for patient i at time interval j, 0:no event, 1:event
  int<lower=0, upper=1> Tind[N, NT]; // patient i having had an event (death) at end of time interval j, 0:no, 1:yes
  int<lower=0, upper=1> Tn[N, NT];   // patient i available for event count at beginning of time interval j, 0:no (after censored time), 1:yes (dead or before censored time)
  real eps;
  eps = 0.000000001;
  for(i in 1:N) {
    for(j in 1:NT) {
      // risk set = 1 if obs.t >= t
      Y[i, j]  = int_step(obs_t[i] - t[j] + eps); // patient i at risk during time interval j
      // counting process jump = 1 if obs.t in [ t[j], t[j+1] )
      //    i.e. if t[j] <= obs.t < t[j+1]
      dN[i, j] = Y[i, j] * ind[i] * int_step(t[j + 1] - obs_t[i] - eps); // patient i event during time interval j
      // patient i dead during time interval j as 0/1
      // eps, because 1 required in the time interval that the patient dies
      Tind[i, j] = int_step(t[j] - obs_t[i] + eps) * ind[i]; 
      // patient i available for count status during time interval j: 
      // either died (observed before and dead after obs time) or 
      // observed later than time interval j
      Tn[i, j] = ind[i] || Y[i, j];
    }
  }
}
parameters {
  real beta[2];
  real eta[4];
  real<lower=0> tau_e;
  real zeta[M]; 
  real<lower=0> tau_z;
  real<lower=0> dL0[nperiods]; 
} 
transformed parameters {
  real<lower=0> sigma_e; 
  real<lower=0> sigma_z; 
  sigma_e = 1 / sqrt(tau_e); 
  sigma_z = 1 / sqrt(tau_z); 
} 
model {
  // priors
  beta ~ normal(0, 1000);
  tau_e  ~ gamma(.001, .001); 
  eta ~ normal(0, sigma_e); 
  tau_z  ~ gamma(.001, .001); 
  zeta ~ normal(0, sigma_z); 
  dL0  ~ gamma(0.001, 0.001);
  for(j in 1:NT) {
    for (i in 1:N) {
      if (Y[i, j] != 0) 
        target += poisson_lpmf(dN[i, j] | 
                               Y[i, j] * 
                               dL0[period[j]] *
                               dt[j] * 
                               exp( beta[1] * Zage[i] + // normalized age
                                    beta[2] * Zkps[i] + // normalized kps
                                    eta[yot[i]] +       // random effects on year of treatment 
                                    zeta[ctr[i]] )      // random effects on center
                               ); 
    }     
  }
}
generated quantities {
    int<lower=0> Nrisk[M, NT];            // nr pts at risk per center per time interval
    int<lower=0> Tnsum[M, NT];            // nr pts with available information per center per time interval
    real<lower=0> V[M];                   // nr of pts per center
    real<lower=0> V_1m[M];                // nr of pts available for observation per center @ t[1]=1m (real to avoid int/int=int)
    real<lower=0> V_2y[M];                // nr of pts available for observation per center @ t[24]=24m (real to avoid int/int=int)
    real<lower=0, upper=1> R_1m[M];       // death rate @ 1m per center 
    real<lower=0, upper=1> R_2y[M];       // survival rate @ 2y per center 
    real<lower=0, upper=1> Rhat_1m;       // overall death rate @ 1m 
    real<lower=0, upper=1> Rhat_2y;       // overall survival rate @ 2y 
    real cumhaz_total;                    // counter for cumhaz based on average yot and average ctr
    real cumhaz_total_i;                  // counter for cumhaz based on specific yot and ctr
    real cumhaz_pred[N, NT];              // predicted cumulative hazard per patient per time interval based on average yot and average ctr
    real cumhaz_pred_i[N, NT];            // predicted cumulative hazard per patient per time interval based on specific yot and average ctr
    real<lower=0, upper=1> Surv_all[NT];  // predicted survival per center per time interval based on specific patient yot and ctr
    real<lower=0, upper=1> SurvM[M, NT];  // predicted survival per center per time interval based on specific patient yot and ctr
    real<lower=0, upper=1> SurvMG[M, NT]; // predicted survival per center per time interval based on average yot and ctr
    real<lower=0> cumsurv_all[NT];        // cumulative predicted survival per center per time interval based on specific patient characteristics
    real<lower=0> cumsurvM[M, NT];        // cumulative predicted survival per center per time interval based on specific patient characteristics
    real<lower=0> cumsurvMG[M, NT];       // cumulative predicted survival per center per time interval based on specific patient characteristics
    real<lower=0,upper=1> surv_pred1m[N]; // predicted survival @ 1m per patient
    real<lower=0,upper=1> surv_pred2y[N]; // predicted survival @ 24m per patient
    int<lower=0>  O_1m_all;               // observed nr of deaths overall at one month
    real<lower=0> E_1m_all;               // expected nr of deaths overall at one month
    real<lower=0> S_1m_all;               // standardized ratio of O/E nr of deaths overall at one month
    int<lower=0>  O_2y_all;               // observed nr of deaths overall at two years
    real<lower=0> E_2y_all;               // expected nr of deaths overall at two years
    real<lower=0> S_2y_all;               // standardized ratio of O/E nr of deaths overall at two years
    real<lower=0> O_mg[M, NT];            // observed nr of deaths per center per time interval (Martingale residual)
    real<lower=0> E_mg[M, NT];            // expected nr of deaths per center per time interval (Martingale residual)
    real<lower=0> S_mg[M, NT];            // observed minus expected nr of deaths per center per time interval
    int<lower=0>  O_1m[M];                // observed nr of pts alive per center @ t[1]=1m
    int<lower=0>  O_2y[M];                // observed nr of pts alive per center @ t[24]=24m
    real<lower=0> E_1m[M];                // expected nr of pts alive per center @ t[1]=1m
    real<lower=0> E_2y[M];                // expected nr of pts alive per center @ t[24]=24m
    real<lower=0> S_1m[M];                // standardized ratio of O/E nr of pts alive per center @ 1m
    real<lower=0> S_2y[M];                // standardized ratio of O/E nr of pts alive per center @ 24m

    // risk table per center per time, starting at 0, ending at one before list time
    Nrisk = rep_array(0, M, NT); // nr of pts at risk (still alive) at beginning of time index j
    Tnsum = rep_array(0, M, NT); // nr of pts with available information (still alive or earlier dead) at beginning of time index j
    for (i in 1:N) {
      for (j in 1:NT) {
        Nrisk[ctr[i], j] = Nrisk[ctr[i], j] + Y[i,j];
        Tnsum[ctr[i], j] = Tnsum[ctr[i], j] + Tn[i,j];
      }
    }
    // calculate center volume & 1m mortality & 2y survival
    V  = rep_array(0.0, M);
    V_1m  = rep_array(0.0, M);
    V_2y  = rep_array(0.0, M);
    R_1m  = rep_array(0.0, M);
    R_2y  = rep_array(0.0, M);
    Rhat_1m = 0.0;
    Rhat_2y = 0.0;
    for (i in 1:N) {
      V[ctr[i]] = V[ctr[i]] + 1; 
      V_1m[ctr[i]] = V_1m[ctr[i]] + Tn[i, j1m];
      V_2y[ctr[i]] = V_2y[ctr[i]] + Tn[i, j2y]; 
    }
    Rhat_1m = 1.0 * sum(Tind[,j1m]) / sum(Tn[,j1m]);
    Rhat_2y = 1 - 1.0 * sum(Tind[,j2y]) / sum(Tn[,j2y]);

    // patient-specific cumulative hazard per time interval, predicted (median) survival and survival probability @ observation time & @ 1m & @ 24m 
    cumsurv_all = rep_array(0.0, NT);
    cumsurvM = rep_array(0.0, M, NT);
    cumsurvMG = rep_array(0.0, M, NT);
    O_mg  = rep_array(0.0, M, NT);
    E_mg  = rep_array(0.0, M, NT);
    O_1m_all = 0;
    E_1m_all = 0;
    O_2y_all = 0;
    E_2y_all = 0;
    S_1m_all = 1.0;
    S_2y_all = 1.0;
    O_1m  = rep_array(0, M);
    E_1m  = rep_array(0.0, M);
    S_1m  = rep_array(1.0, M);
    O_2y  = rep_array(0, M);
    E_2y  = rep_array(0.0, M);
    S_2y  = rep_array(1.0, M);

    for (i in 1:N) {
      cumhaz_total = 0;
      cumhaz_total_i = 0;
      for (j in 1:NT) {

        // prediction based on average ctr
        cumhaz_pred[i, j] = cumhaz_total + 
                            dL0[period[j]] *
                            dt[j] *         
                            exp( beta[1] * Zage[i] + // normalized age
                                 beta[2] * Zkps[i] + // normalized kps
                                 eta[yot[i]] +       // year of treatment effect
                                 0);                 // average center effect = 0
        cumhaz_total = cumhaz_pred[i, j];
        
        // martingale residuals
          // observed nr of deaths up to t[j] in center k
          O_mg[ctr[i], j] = O_mg[ctr[i], j] + Tind[i, j];
          // expected nr of deaths up to t[j] in center k in available patients
          E_mg[ctr[i], j] = E_mg[ctr[i], j] + Tn[i, j] * (1 - exp( -cumhaz_total ));
          cumsurvMG[ctr[i], j] = cumsurvMG[ctr[i], j] + Tn[i, j] * exp( -cumhaz_total );
          cumsurv_all[j] = cumsurv_all[j] + Tn[i, j] * exp( -cumhaz_total );
        
        // prediction based on patient-specific ctr
        cumhaz_pred_i[i, j] = cumhaz_total_i + 
                            dL0[period[j]] *
                            dt[j] *         
                            exp( beta[1] * Zage[i] +    // normalized age
                                 beta[2] * Zkps[i] +    // normalized kps
                                 eta[yot[i]] +          // year of treatment effect
                                 zeta[ctr[i]]);         // center effect
        cumhaz_total_i = cumhaz_pred_i[i, j];
        cumsurvM[ctr[i], j] = cumsurvM[ctr[i], j] + Tn[i, j] * exp( -cumhaz_total_i );
      } // loop time index j

      // nr of pts ** dead ** @ 1m
      surv_pred1m[i]  = exp( -cumhaz_pred[i, j1m ] );
      if (Tn[i, j1m] == 1) { // available observation @1m
        O_1m_all = O_1m_all + Tind[i, j1m]; 
        E_1m_all = E_1m_all + (1 - surv_pred1m[i]);
        O_1m[ctr[i]] = O_1m[ctr[i]] + Tind[i, j1m]; 
        E_1m[ctr[i]] = E_1m[ctr[i]] + (1 - surv_pred1m[i]);                 
      }

     // nr of pts ** alive ** @ 24m
     surv_pred2y[i]  = exp( -cumhaz_pred[i, j2y ] );
      if (Tn[i, j2y] == 1) { // available observation @2y
        O_2y_all = O_2y_all + (1 - Tind[i, j2y]) ; //Y[i, j2y];
        E_2y_all = E_2y_all + surv_pred2y[i]; // add predicted survival in available pts
        O_2y[ctr[i]] = O_2y[ctr[i]] + (1 - Tind[i, j2y]);// Y[i, j2y];
        E_2y[ctr[i]] = E_2y[ctr[i]] + surv_pred2y[i];
      }

    } // loop patient index i

    S_1m_all = O_1m_all / E_1m_all;
    S_2y_all = O_2y_all / E_2y_all;
    Surv_all = rep_array(0.0, NT);
    SurvM = rep_array(0.0, M, NT);
    SurvMG = rep_array(0.0, M, NT);
    for (k in 1:M) {
      for (j in 1:NT) {
        if (E_1m[k] > 0) S_1m[k] = O_1m[k] / E_1m[k];
        if (E_2y[k] > 0) S_2y[k] = O_2y[k] / E_2y[k];
        if (V_1m[k] > 0) R_1m[k] = O_1m[k] / V_1m[k];
        if (V_2y[k] > 0) R_2y[k] = O_2y[k] / V_2y[k];
        if (E_mg[k, j] > 0) S_mg[k, j] = O_mg[k, j] / E_mg[k, j];
        if (Tnsum[k, j] > 0) {
          SurvM[k, j]  = cumsurvM[k, j] / Tnsum[k, j];
          SurvMG[k, j] = cumsurvMG[k, j] / Tnsum[k, j];
        }
      }
    }
    for (j in 1:NT) {
      Surv_all[j]  = cumsurv_all[j] / sum(Tn[, j]);
    }


}


