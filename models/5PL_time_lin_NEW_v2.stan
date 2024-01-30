data{
  int N;
  int mortality[N];
  int tested[N];
  vector[N] concentration;
  int N_1;
  vector[N_1] concentration_sim;
  int Y; // year variable
  vector[N] year;

}

parameters{
  real<lower=0> B;
  real<lower=0> E;
  real F;
  real G;
  real<lower=0, upper=1> Z;

}

model{
  // likelihood
  for(i in 1:N){
    real f;
    real C_temp = F+G*year[i];

    if (concentration[i] > 0)
    f = 1 -  ((1 - Z) /((1+exp(B*(log(concentration[i])-C_temp)))^E));
    else
    f = Z;

    mortality[i]~binomial(tested[i], f);
  }

  // priors
  B ~ normal(5,10);
  E ~ normal(7,10);
  F ~ normal(0,5);
  G ~ normal(0,5);
  Z ~ normal(0,0.1);
}

generated quantities{
  matrix[N_1,Y] mean_mortality_sim;
  matrix[N_1,Y] mean_mortality_sim_bm;
  vector[N] mean_mortality;
  vector[N] mean_mortality_bm;
  vector[N] LogLikelihood;

  for(i in 1:N_1) {
    for(j in 0:(Y-1)){
      real C_temp = (F+G*j);
      mean_mortality_sim[i,(j+1)] = 1 + ((-1)/((1+exp(B*(log(concentration_sim[i])-C_temp)))^E));
      mean_mortality_sim_bm[i,(j+1)] = 1 - ((1-Z)/((1+exp(B*(log(concentration_sim[i])-C_temp)))^E));
    }}

  for(i in 1:N) {
    real C_temp = (F+G*year[i]);
    mean_mortality[i] =  1 + ((-1)/((1+exp(B*(log(concentration[i])-C_temp)))^E));
    mean_mortality_bm[i] =  1 - ((1-Z)/((1+exp(B*(log(concentration[i])-C_temp)))^E));
    LogLikelihood[i] = binomial_lpmf(mortality[i]|tested[i], mean_mortality_bm[i]);

  }

}
