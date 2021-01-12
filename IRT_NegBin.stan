data{
  int<lower=0> n_item;
  int<lower=0> n_response;
  int<lower=0> total_length;
  int<lower=0> item_response[total_length];
  int<lower=0> beta_indecator[n_item];
}
parameters{
  real theta[n_response];  
  real beta;                   // vedio type effect
  real<lower=0,upper=1> zero;  // zero infalted
  real phi;
}
model{
  theta ~ normal(0,1);
  phi ~ cauchy(0,1);
  for (i in 1:n_response){
    for (j in 1:n_item){
      if (item_response[n_item *(i - 1) + j] == 0){
        target += log_sum_exp(bernoulli_lpmf(1 | zero), bernoulli_lpmf(0 | zero) +  neg_binomial_2_log_lpmf(item_response[n_item *(i - 1) + j] | theta[i] +  beta * beta_indecator[j],phi));
      }else{
        target += bernoulli_lpmf(0 | zero) + neg_binomial_2_log_lpmf( item_response[n_item *(i - 1) + j] |  theta[i]  + beta * beta_indecator[j],phi);
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
      item_response_rep[n_item *(i - 1) + j] = neg_binomial_2_log_rng(theta[i] +  beta * beta_indecator[j],phi);
    }  
    }
  }
}

