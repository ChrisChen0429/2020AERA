data{
  int<lower=0> n_item;
  int<lower=0> n_response;
  int<lower=0> total_length;
  int<lower=0> item_response[total_length];
}
parameters{
  real theta[3]; 
  real beta[3];
  real<lower=0,upper=1> zero;
}
model{
  theta ~ cauchy(0,1);
  int<lower=0,upper=3> theta_type[n_response];
  int<lower=0,upper=3> beta_type[n_item];
  for (i in 1:n_response){
    for (j in 1:n_item){
      if (item_response[n_item *(i - 1) + j] ==0 ){
        target += log_sum_exp(bernoulli_lpmf(1 | zero), bernoulli_lpmf(0 | zero) +  poisson_log_lpmf( item_response[n_item *(i - 1) + j] | theta[theta_type[i]] - beta[beta_type[j]]));
      }else{
        target += bernoulli_lpmf(0 | zero) + poisson_log_lpmf( item_response[n_item *(i - 1) + j] |  theta[theta_type[i]] - beta[beta_type[j]]);
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
      item_response_rep[n_item *(i - 1) + j] = poisson_log_rng( theta[theta_type[i]] - beta[beta_type[j]]);
    }  
    }
  }
}

