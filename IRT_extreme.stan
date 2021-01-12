functions {
  real gpareto_lpdf(vector y, real ymin, real k, real sigma) {
    // generalised Pareto log pdf 
    int N = rows(y);
    real inv_k = inv(k);
    if (k<0 && max(y-ymin)/sigma > -inv_k)
      reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma)
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma)
    if (fabs(k) > 1e-15)
      return -(1+inv_k)*sum(log1p((y-ymin) * (k/sigma))) -N*log(sigma);
    else
      return -sum(y-ymin)/sigma -N*log(sigma); // limit k->0
  }
  real gpareto_cdf(vector y, real ymin, real k, real sigma) {
    // generalised Pareto cdf
    real inv_k = inv(k);
    if (k<0 && max(y-ymin)/sigma > -inv_k)
      reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma)
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma)
    if (fabs(k) > 1e-15)
      return exp(sum(log1m_exp((-inv_k)*(log1p((y-ymin) * (k/sigma))))));
    else
      return exp(sum(log1m_exp(-(y-ymin)/sigma))); // limit k->0
  }
  real gpareto_lcdf(vector y, real ymin, real k, real sigma) {
    // generalised Pareto log cdf
    real inv_k = inv(k);
    if (k<0 && max(y-ymin)/sigma > -inv_k)
      reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma)
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma)
    if (fabs(k) > 1e-15)
      return sum(log1m_exp((-inv_k)*(log1p((y-ymin) * (k/sigma)))));
    else
      return sum(log1m_exp(-(y-ymin)/sigma)); // limit k->0
  }
  real gpareto_lccdf(vector y, real ymin, real k, real sigma) {
    // generalised Pareto log ccdf
    real inv_k = inv(k);
    if (k<0 && max(y-ymin)/sigma > -inv_k)
      reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma)
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma)
    if (fabs(k) > 1e-15)
      return (-inv_k)*sum(log1p((y-ymin) * (k/sigma)));
    else
      return -sum(y-ymin)/sigma; // limit k->0
  }
  real gpareto_rng(real ymin, real k, real sigma) {
    // generalised Pareto rng
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma)
    if (fabs(k) > 1e-15)
      return ymin + (uniform_rng(0,1)^-k -1) * sigma / k;
    else
      return ymin - sigma*log(uniform_rng(0,1)); // limit k->0
  }
}

data {
  real ymin;
  int<lower=0> N;
  vector<lower=ymin>[N] y;
  int<lower=0> Nt;
  vector<lower=ymin>[Nt] yt;
}
data{
  int<lower=0> n_item;
  int<lower=0> n_response;
  int<lower=0> total_length;
  int<lower=0> item_response[total_length];
  int<lower=0> beta_indecator[n_item];
}

parameters{
  real theta[n_response]; 
  real beta[2];
  real<lower=0,upper=1> zero; 
}

model{
  theta ~ cauchy(0,1);
  beta ~  cauchy(0,1);
  for (i in 1:n_response){
    for (j in 1:n_item){
      if (item_response[n_item *(i - 1) + j] == 0){
        target += log_sum_exp(bernoulli_lpmf(1 | zero), bernoulli_lpmf(0 | zero) +  poisson_log_lpmf(item_response[n_item *(i - 1) + j] | theta[i] - beta[beta_indecator[j]]));
      }else{
        target += bernoulli_lpmf(0 | zero) + poisson_log_lpmf( item_response[n_item *(i - 1) + j] |  theta[i] - beta[beta_indecator[j]]);
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
      item_response_rep[n_item *(i - 1) + j] = poisson_log_rng(theta[i] - beta[beta_indecator[j]]);
    }  
    }
  }
}


transformed data {
  real ymax = max(y);
}
parameters {
  real<lower=0> sigma; 
  real<lower=-sigma/(ymax-ymin)> k; 
}
model {
  y ~ gpareto(ymin, k, sigma);
}
generated quantities {
  vector[N] log_lik;
  vector[N] yrep;
  vector[Nt] predccdf;
  for (n in 1:N) {
    log_lik[n] = gpareto_lpdf(rep_vector(y[n],1) | ymin, k, sigma);
    yrep[n] = gpareto_rng(ymin, k, sigma);
  }
  for (nt in 1:Nt)
    predccdf[nt] = exp(gpareto_lccdf(rep_vector(yt[nt],1) | ymin, k, sigma));
}