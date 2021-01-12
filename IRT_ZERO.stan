data{
  int<lower=0> n_item;                 // number of item
  int<lower=0> n_response;             // number of response 
  int<lower=0> total_length;           // number of observation (n_item * n_response)
  int<lower=0> item_response[total_length];  // observations 
}
parameters{
  real theta[n_response];                // participant parameters
  real beta[n_item];                     // item paramerters 
  real<lower=0,upper=1> lambda;          // zero infalted
}
model{
  theta ~ cauchy(0,1);                   // weekly informative prior
  for (i in 1:n_response){
    for (j in 1:n_item){
      if (item_response[n_item *(i - 1) + j] == 0){
        target += log_sum_exp(bernoulli_lpmf(1 | lambda), bernoulli_lpmf(0 | lambda) +  poisson_log_lpmf(item_response[n_item *(i - 1) + j] | theta[i] + beta[j]));
      }else{
        target += bernoulli_lpmf(0 | lambda) + poisson_log_lpmf( item_response[n_item *(i - 1) + j] | theta[i]  +  beta[j]);
      }
    }
  }
}
generated quantities{
  int<lower=0> item_response_rep[total_length];
  real<lower=0,upper=1> lambda_i_rep[total_length];
  for (i in 1:n_response){
    for (j in 1:n_item){
      lambda_i_rep[n_item *(i - 1) + j] = uniform_rng(0,1);
      if (lambda_i_rep[n_item *(i - 1) + j] < lambda){
        item_response_rep[n_item *(i - 1) + j] = 0;
      }
    else{
      item_response_rep[n_item *(i - 1) + j] = poisson_log_rng( theta[i] + beta[j]);
    }  
  }
  }
}