data{
  int N;
  int mortality[N];
  int tested[N];
  vector[N] concentration;
  int N_1;
  vector[N_1] concentration_sim;

}

parameters{
  real<lower=0> B;
  real C;
  real<lower=0> E;
  real<lower=0, upper=1> phi;
}

model{
  // likelihood
  for(i in 1:N){
    real f;

    if (concentration[i] > 0)
    f = 1 -  ((1 - phi) /((1+exp(B*(log(concentration[i])-C)))^E));
    else
    f = phi;

    mortality[i]~binomial(tested[i], f);
  }

  // priors
  B ~ normal(5,10);
  C ~ normal(0,5);
  E ~ normal(7,10);
  phi ~ normal(0,0.1);
}

generated quantities{
  vector[N_1] mean_mortality_sim;
  vector[N_1] mean_mortality_sim_bm;
  vector[N] mean_mortality;
  vector[N] mean_mortality_bm;
  vector[N] LogLikelihood;

  for(i in 1:N_1) {
    mean_mortality_sim[i] = 1 + ((-1)/((1+exp(B*(log(concentration_sim[i])-C)))^E));
    mean_mortality_sim_bm[i] = 1 -  ((1 - phi)/((1+exp(B*(log(concentration_sim[i])-C)))^E));
    }

  for(i in 1:N) {
    mean_mortality[i] = 1 + ((-1)/((1+exp(B*(log(concentration[i])-C)))^E));
    mean_mortality_bm[i] = 1 -  ((1 - phi)/((1+exp(B*(log(concentration[i])-C)))^E));
    LogLikelihood[i] = binomial_lpmf(mortality[i]|tested[i], mean_mortality_bm[i]);
  }

}
