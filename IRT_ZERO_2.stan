data{
  int<lower=0> n_item;
  int<lower=0> n_response;
  int<lower=0> total_length;
  int<lower=0> item_response[total_length];
}
parameters{
  real theta[n_response];  
  real beta[n_item];                   // vedio type effect
  real<lower=0,upper=1> zero;  // zero infalted
}
model{
  theta ~ cauchy(0,1);
  beta ~ cauchy(0,1);
  for (i in 1:n_response){
    for (j in 1:n_item){
      if (item_response[n_item *(i - 1) + j] == 0){
        target += log_sum_exp(bernoulli_lpmf(1 | zero), bernoulli_lpmf(0 | zero) +  poisson_log_lpmf(item_response[n_item *(i - 1) + j] | theta[i] + beta[j]));
      }else{
        target += bernoulli_lpmf(0 | zero) + poisson_log_lpmf( item_response[n_item *(i - 1) + j] |  theta[i]  +  beta[j]);
      }
    }
  }
}
generated quantities{
  int<lower=0> item_response_rep[total_length];
  real<lower=0,upper=1> zero_i_rep[total_length];
  for (i in 1:n_response){
    for (j in 1:n_item){
      zero_i_rep[n_item *(i - 1) + j] = uniform_rng(0,1);
      if (zero_i_rep[n_item *(i - 1) + j] < zero){
        item_response_rep[n_item *(i - 1) + j] = 0;
      }
    else{
      item_response_rep[n_item *(i - 1) + j] = poisson_log_rng(theta[i] +  beta[j]);
    }  
    }
  }
}