functions {
  // return (A \otimes B) v where:
  // A is n1 x n1, B = n2 x n2, V = n2 x n1 = reshape(v,n2,n1)
  matrix kron_mvprod(matrix A, matrix B, matrix V) {
    return (B * V) * transpose(A);
  }
}

data {
  /////////////////////////////////////////////////////////////////////////////
  // description of data dimensions
  /////////////////////////////////////////////////////////////////////////////

  // number of states/regions and seasons
  // corresponds to n2 in Flaxman et al. Appendix A.3
  // (except that we will really use 2*SU)
  int<lower=1> SU;

  // number of states
  int<lower=1> S;

  // largest number of seasons for each state
  int<lower=1> U;

  // for each SU, what is corresponding S and U?
  int<lower=1> SU_to_S_map[SU];
  int<lower=1> SU_to_U_map[SU];

  // number of weeks of data for region with the most weeks
  // corresponds to n1 in Flaxman et al. Appendix A.3
  int<lower=1> T_max;

  // timing of holidays
  //int<lower = 1> thanksgiving_week_ind[SU];
  int<lower = 1> christmas_week_ind[SU];
  //int<lower = 1> postchristmas_week_ind[SU];
  //int<lower = 0> max_num_holidays;
  //int<lower = 0> num_thanksgiving;
  //int<lower = 0> num_christmas;
  //int<lower = 0> num_postchristmas;

  //int<lower = 0> num_thanksgiving_by_state[S];
  //int<lower = 0> num_christmas_by_state[S];
  //int<lower = 0> num_postchristmas_by_state[S];

  //int<lower = 1> thanksgiving_week_inds[max_num_holidays, S];
  //int<lower = 1> christmas_week_inds[max_num_holidays, S];
  //int<lower = 1> postchristmas_week_inds[max_num_holidays, S];

  // ili data
  int<lower = 0> num_obs_ili[SU];
  int<lower = 0> total_num_obs_ili;

  real<lower = 0, upper = 1> ili_t_concat[total_num_obs_ili];
  int ili_obs_t[total_num_obs_ili];

  int<lower = 0> num_pred_ili[SU];
  int<lower = 0> total_num_pred_ili;
  int ili_pred_t[total_num_pred_ili];

  //real outpatient_success_adjust;
  //real outpatient_fail_adjust;

  // virology data
  int<lower = 0, upper = T_max> num_obs_virology[SU];
  int<lower = 0> total_num_obs_virology;
  real<lower = 0, upper = 1> prop_a_t_concat[total_num_obs_virology];
  int<lower = 0> p_t_concat[total_num_obs_virology];
  int<lower = 1, upper = T_max> virology_obs_t[total_num_obs_virology];

  // number of columns in design matrix
  int<lower = 1> DX;
  matrix[T_max, DX] X;
  real knots[DX];

  // number of lags to use for autoregressive spline structure
  int<lower = 1> max_lag;
}

transformed data {
  matrix[DX, DX] xd;
  real logit_ili_t_concat[total_num_obs_ili];
  real logit_prop_a_t_concat[total_num_obs_virology];

  // distances between knots in time axis
  for (i in 1:DX) {
    xd[i, i] = 0;
    for (j in (i+1):DX) {
      xd[i, j] = (knots[j]-knots[i]);
      xd[j, i] = xd[i, j];
    }
  }

  logit_ili_t_concat = logit(ili_t_concat);

  logit_prop_a_t_concat = logit(prop_a_t_concat);
}

parameters {
  vector[DX] beta_mean_raw_a;
  vector[DX] beta_mean_raw_b;
  matrix[DX, S] beta_state_raw_a;
  matrix[DX, S] beta_state_raw_b;
  matrix[DX, 1] beta_seasonal_raw_a[U];
  matrix[DX, 1] beta_seasonal_raw_b[U];
  matrix[DX, S] beta_raw_a[U];
  matrix[DX, S] beta_raw_b[U];

  real<lower = 0> sigma_beta_mean_a;
  real<lower = 0> sigma_beta_mean_b;
  real<lower = 0> sigma_beta_state_a;
  real<lower = 0> sigma_beta_state_b;
  real<lower = 0> sigma_beta_seasonal_a;
  real<lower = 0> sigma_beta_seasonal_b;
  real<lower = 0> sigma_beta_a;
  real<lower = 0> sigma_beta_b;
  real<lower=0, upper=1> gamma_a[max_lag]; // autocorrelation for splines
  real<lower=0, upper=1> gamma_b[max_lag]; // autocorrelation for splines
  real<lower=0, upper=1> gamma_state_a[max_lag]; // autocorrelation for splines
  real<lower=0, upper=1> gamma_state_b[max_lag]; // autocorrelation for splines
  real<lower=0, upper=1> gamma_seasonal_a[max_lag]; // autocorrelation for splines
  real<lower=0, upper=1> gamma_seasonal_b[max_lag]; // autocorrelation for splines

  cholesky_factor_corr[S] L_S_state_a;
  cholesky_factor_corr[S] L_S_state_b;
  cholesky_factor_corr[U] L_U_seasonal_a;
  cholesky_factor_corr[U] L_U_seasonal_b;

  cholesky_factor_corr[S] L_S_a;
  cholesky_factor_corr[S] L_S_b;
  cholesky_factor_corr[U] L_U_a;
  cholesky_factor_corr[U] L_U_b;

  real state_effect_raw[S];
  real<lower = 0> sigma_state_effect;

  real christmas_effect_raw[SU];
  real mean_christmas_effect;
  real<lower = 0> sigma_christmas_effect;

  //vector[total_num_obs_ili] outpatient_overdispersion_raw;
  real<lower = 0> sigma_outpatient_overdispersion;

  //vector[total_num_obs_virology] virology_overdispersion_raw;
  real<lower = 0> sigma_virology_overdispersion;
}

transformed parameters {
  vector[total_num_obs_ili] mean_logit_wili;
  vector[total_num_obs_virology] mean_logit_virology_a;

  {
      vector[DX] beta_mean_a;
  vector[DX] beta_mean_b;
  matrix[DX, S] beta_state_a;
  matrix[DX, S] beta_state_b;
  matrix[DX, U] beta_seasonal_a;
  matrix[DX, U] beta_seasonal_b;
  matrix[S * DX, U] beta_a;
  matrix[S * DX, U] beta_b;
  vector[T_max] log_I_a[SU];
  vector[T_max] log_I_b[SU];
  //real log_I_a_concat[total_num_obs_ili];
  //real log_I_b_concat[total_num_obs_ili];
  //real log_I_concat[total_num_obs_ili];
  //real logit_state_offset_concat[total_num_obs_ili];
  //real logit_christmas_concat[total_num_obs_ili];
  //real logit_binom_prob_outpatient[total_num_obs_ili];
  //real log_binom_prob_virology_a[total_num_obs_virology];
  vector[total_num_obs_ili] log_I_a_concat;
  vector[total_num_obs_ili] log_I_b_concat;
  vector[total_num_obs_ili] log_I_concat;
  vector[total_num_obs_ili] logit_state_offset_concat;
  vector[total_num_obs_ili] logit_christmas_concat;
  //vector[total_num_obs_ili] logit_binom_prob_outpatient;

  //vector[total_num_obs_virology] logit_binom_prob_virology_a;

  beta_mean_a = sigma_beta_mean_a * beta_mean_raw_a;
  beta_mean_b = sigma_beta_mean_b * beta_mean_raw_b;

  {
    matrix[DX, DX] Sigma1_state_a;
    matrix[DX, DX] Sigma1_state_b;
    matrix[DX, DX] L_Sigma1_state_a;
    matrix[DX, DX] L_Sigma1_state_b;
    matrix[DX, S] state_SW_blocks_a;
    matrix[DX, S] state_SW_blocks_b;

    matrix[DX, DX] Sigma1_seasonal_a;
    matrix[DX, DX] Sigma1_seasonal_b;
    matrix[DX, DX] L_Sigma1_seasonal_a;
    matrix[DX, DX] L_Sigma1_seasonal_b;
    matrix[DX, U] seasonal_SW_blocks_a;
    matrix[DX, U] seasonal_SW_blocks_b;

    matrix[DX, DX] Sigma1_a;
    matrix[DX, DX] Sigma1_b;
    matrix[DX, DX] L_Sigma1_a;
    matrix[DX, DX] L_Sigma1_b;
    matrix[S * DX, U] SW_blocks_a;
    matrix[S * DX, U] SW_blocks_b;
    matrix[DX, S] temp_a;
    matrix[DX, S] temp_b;
    real temp;

    // get betas for state and season terms
    for(i in 1:DX) {
      for(j in 1:DX) {
        Sigma1_state_a[i,j] = gamma_state_a[1]^xd[i, j];
        Sigma1_state_b[i,j] = gamma_state_b[1]^xd[i, j];
        Sigma1_seasonal_a[i,j] = gamma_seasonal_a[1]^xd[i, j];
        Sigma1_seasonal_b[i,j] = gamma_seasonal_b[1]^xd[i, j];
      }
      Sigma1_state_a[i,i] = Sigma1_state_a[i,i] + .00001;
      Sigma1_state_b[i,i] = Sigma1_state_b[i,i] + .00001;
      Sigma1_seasonal_a[i,i] = Sigma1_seasonal_a[i,i] + .00001;
      Sigma1_seasonal_b[i,i] = Sigma1_seasonal_b[i,i] + .00001;
    }
    L_Sigma1_state_a = cholesky_decompose(Sigma1_state_a);
    L_Sigma1_state_b = cholesky_decompose(Sigma1_state_b);
    L_Sigma1_seasonal_a = cholesky_decompose(Sigma1_seasonal_a);
    L_Sigma1_seasonal_b = cholesky_decompose(Sigma1_seasonal_b);

    beta_state_a = sigma_beta_state_a * kron_mvprod(L_S_state_a, L_Sigma1_state_a, beta_state_raw_a);
    beta_state_b = sigma_beta_state_b * kron_mvprod(L_S_state_b, L_Sigma1_state_b, beta_state_raw_b);

    for(u in 1:U) {
      // same as computation for beta_a and beta_b below, but we have no S for ili (no age groups)
      seasonal_SW_blocks_a[1:DX, u] = (L_Sigma1_seasonal_a * beta_seasonal_raw_a[u])[,1];
      seasonal_SW_blocks_b[1:DX, u] = (L_Sigma1_seasonal_b * beta_seasonal_raw_b[u])[,1];
    }

    beta_seasonal_a = sigma_beta_seasonal_a * seasonal_SW_blocks_a * L_U_seasonal_a; // in reverse order of seasons since we did not transpose L_U
    beta_seasonal_b = sigma_beta_seasonal_b * seasonal_SW_blocks_b * L_U_seasonal_b; // in reverse order of seasons since we did not transpose L_U

    // get betas for individual region x season terms
    for(i in 1:DX) {
      for(j in 1:DX) {
        Sigma1_a[i,j] = gamma_a[1]^xd[i, j];
        Sigma1_b[i,j] = gamma_b[1]^xd[i, j];
      }
      Sigma1_a[i,i] = Sigma1_a[i,i] + .00001;
      Sigma1_b[i,i] = Sigma1_b[i,i] + .00001;
    }
    L_Sigma1_a = cholesky_decompose(Sigma1_a);
    L_Sigma1_b = cholesky_decompose(Sigma1_b);

    for(u in 1:U) {
      temp_a = kron_mvprod(L_S_a, L_Sigma1_a, beta_raw_a[u]);
      temp_b = kron_mvprod(L_S_b, L_Sigma1_b, beta_raw_b[u]);
      for(s in 1:S) {
        beta_a[((s - 1)*DX + 1):((s - 1)*DX + DX), u] = sigma_beta_a * L_Sigma1_a * ((beta_raw_a[u])[1:DX, s]);
        beta_b[((s - 1)*DX + 1):((s - 1)*DX + DX), u] = sigma_beta_b * L_Sigma1_b * ((beta_raw_b[u])[1:DX, s]);
      }
    }

    {
      vector[DX] beta_su;

      for(su in 1:SU) {
        // Compute binomial probabilities
        real log_numeric_offset = log(10.0^(-10.0));
        for(dx in 1:DX) {
          beta_su[dx] = beta_mean_a[dx] + beta_state_a[dx, SU_to_S_map[su]] +
            beta_seasonal_a[dx, SU_to_U_map[su]] +
            beta_a[(SU_to_S_map[su] - 1) * DX + dx, U + 1 - SU_to_U_map[su]];
        }
        log_I_a[su] = -1.0 * log1p_exp(-1.0* X * beta_su); // + log1m_exp(log(2.0) + log_numeric_offset);

        for(dx in 1:DX) {
          beta_su[dx] = beta_mean_b[dx] + beta_state_b[dx, SU_to_S_map[su]] +
            beta_seasonal_b[dx, SU_to_U_map[su]] +
            beta_b[(SU_to_S_map[su] - 1) * DX + dx, U + 1 - SU_to_U_map[su]];
        }
        log_I_b[su] = -1.0 * log1p_exp(-1.0 * X * beta_su); // + log1m_exp(log(2.0) + log_numeric_offset);
      }
    }

    {
      int concat_ind_ili = 0;
      int concat_ind_vir = 0;
      for(s_ind in 1:SU) {
        for(t in 1:num_obs_ili[s_ind]) {
          concat_ind_ili += 1;
          log_I_a_concat[concat_ind_ili] = log_I_a[s_ind][ili_obs_t[concat_ind_ili]];
          log_I_b_concat[concat_ind_ili] = log_I_b[s_ind][ili_obs_t[concat_ind_ili]];
          log_I_concat[concat_ind_ili] = log_sum_exp(log_I_a_concat[concat_ind_ili], log_I_b_concat[concat_ind_ili]);

          logit_state_offset_concat[concat_ind_ili] = sigma_state_effect * state_effect_raw[SU_to_S_map[s_ind]];

          if(ili_obs_t[concat_ind_ili] == christmas_week_ind[s_ind]) {
            logit_christmas_concat[concat_ind_ili] = mean_christmas_effect + sigma_christmas_effect * christmas_effect_raw[s_ind];
            mean_logit_wili[concat_ind_ili] = log_I_concat[concat_ind_ili] - log1m_exp(log_I_concat[concat_ind_ili]) + logit_state_offset_concat[concat_ind_ili] + mean_christmas_effect + sigma_christmas_effect * christmas_effect_raw[s_ind];
            #logit_binom_prob_outpatient[concat_ind_ili] = logit_binom_prob_outpatient[concat_ind_ili] + mean_christmas_effect + sigma_christmas_effect * christmas_effect_raw[s_ind];
            // christmas week
            //real log_holiday_scale = log(0.8);
            //real log_effect_components[3];
            //log_effect_components[1] = log_I_a[s_ind][ili_obs_t[concat_ind_ili]];
            //log_effect_components[2] = log_I_b[s_ind][ili_obs_t[concat_ind_ili]];
            //log_effect_components[3] = log_holiday_scale - log1p_exp(-1.0*(mean_christmas_effect + sigma_christmas_effect * christmas_effect_raw[s_ind]));
            //log_binom_prob_outpatient[concat_ind_ili] = log_diff_exp(log_sum_exp(log_effect_components), log_holiday_scale - log(2.0));
          } else {
            logit_christmas_concat[concat_ind_ili] = 0;
            mean_logit_wili[concat_ind_ili] = log_I_concat[concat_ind_ili] - log1m_exp(log_I_concat[concat_ind_ili]) + logit_state_offset_concat[concat_ind_ili];
          }
        }

        //logit_binom_prob_outpatient[concat_ind_ili] = log_I_concat[concat_ind_ili] + log1m_exp(log_I_concat[concat_ind_ili]) + sigma_outpatient_overdispersion * outpatient_overdispersion_raw[concat_ind_ili];

        for(t in 1:num_obs_virology[s_ind]) {
          concat_ind_vir += 1;
          temp = log_I_a[s_ind][virology_obs_t[concat_ind_vir]] - log_sum_exp(log_I_a[s_ind][virology_obs_t[concat_ind_vir]], log_I_b[s_ind][virology_obs_t[concat_ind_vir]]);
          mean_logit_virology_a[concat_ind_vir] = temp - log1m_exp(temp);
        }
      }
      //logit_binom_prob_outpatient = log_I_concat + log1m_exp(log_I_concat) + logit_state_offset_concat + logit_christmas_concat + sigma_outpatient_overdispersion * outpatient_overdispersion_raw;
      //logit_binom_prob_outpatient_pred = log_I_pred_concat + log1m_exp(log_I_pred_concat) + logit_state_offset_pred_concat + logit_christmas_pred_concat + sigma_outpatient_overdispersion * outpatient_overdispersion_pred_raw;
    }
  }
  }
}

model {
  state_effect_raw ~ normal(0.0, 1.0);

  //thanksgiving_effect ~ normal(0.0, 1.0);
  christmas_effect_raw ~ normal(0.0, 1.0);
  //postchristmas_effect ~ normal(0.0, 1.0);

  //outpatient_overdispersion_raw ~ normal(0.0, 1.0);
  //virology_overdispersion_raw ~ normal(0.0, 1.0);

  // GP on spline basis function coefficients
  L_S_state_a ~ lkj_corr_cholesky(1);
  L_S_state_b ~ lkj_corr_cholesky(1);
  L_U_seasonal_a ~ lkj_corr_cholesky(1);
  L_U_seasonal_b ~ lkj_corr_cholesky(1);

  L_S_a ~ lkj_corr_cholesky(1);
  L_S_b ~ lkj_corr_cholesky(1);
  L_U_a ~ lkj_corr_cholesky(1);
  L_U_b ~ lkj_corr_cholesky(1);

  // beta ~ normal(0, L'L);
  beta_mean_raw_a ~ normal(0.0, 1.0);
  beta_mean_raw_b ~ normal(0.0, 1.0);
  for(u in 1:U) {
    for(d in 1:DX) {
      beta_raw_a[u][d] ~ normal(0.0, 1.0);
      beta_raw_b[u][d] ~ normal(0.0, 1.0);
      beta_seasonal_raw_a[u][d] ~ normal(0.0, 1.0);
      beta_seasonal_raw_b[u][d] ~ normal(0.0, 1.0);
    }
  }

  for(d in 1:DX) {
    beta_state_raw_a[d] ~ normal(0.0, 1.0);
    beta_state_raw_b[d] ~ normal(0.0, 1.0);
  }

  // TODO: add prior for gamma
  sigma_beta_mean_a ~ lognormal(0, 1);
  sigma_beta_mean_b ~ lognormal(0, 1);
  sigma_beta_state_a ~ lognormal(0, 1);
  sigma_beta_state_b ~ lognormal(0, 1);
  sigma_beta_seasonal_a ~ lognormal(0, 1);
  sigma_beta_seasonal_b ~ lognormal(0, 1);
  sigma_beta_a ~ lognormal(0,1);
  sigma_beta_b ~ lognormal(0,1);
  sigma_outpatient_overdispersion ~ lognormal(0, 1);
  sigma_virology_overdispersion ~ lognormal(0, 1);
  sigma_state_effect ~ lognormal(0, 1);

  // data model
  logit_ili_t_concat ~ normal(mean_logit_wili, sigma_outpatient_overdispersion);
  logit_prop_a_t_concat ~ normal(mean_logit_virology_a, sigma_virology_overdispersion);
}
