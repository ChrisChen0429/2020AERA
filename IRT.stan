data{
  int<lower=0> n_item;
  int<lower=0> n_response;
  int<lower=0> total_length;
  int<lower=0> item_response[total_length];
  //int<lower=0,upper=1> vialogue_type[n_item];
  //int<lower=0,upper=1> people_type[n_response];
  //int<lower=0,upper=1> people_department[n_response];
}
parameters{
  real theta[n_response]; 
  real beta[2];
  //real beta_vialogue_type;
  //real theta_people_type;
  //real theta_people_department;
}
model{
  theta ~ cauchy(0,1);
  beta ~ cauchy(0,1);
  for (i in 1:n_response){
    for (j in 1:n_item){
      item_response[n_item *(i - 1) + j]  ~ poisson_log((theta[i] - beta[j]));
    }
  }
}
generated quantities{
  int<lower=0> item_response_rep[total_length];
  for (i in 1:n_response){
      for (j in 1:n_item){
      item_response_rep[n_item *(i - 1) + j] = poisson_log_rng((theta[i] - beta[j]));
    }  
    }
}