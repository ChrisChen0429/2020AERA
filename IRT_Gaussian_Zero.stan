functions {
  vector gp_pred_rng(real[] x2,
                     vector y1, real[] x1,
                     real alpha, real rho, real sigma, real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K =   cov_exp_quad(x1, alpha, rho)
                         + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      vector[N2] f2_mu = (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 =   cov_exp_quad(x2, alpha, rho) - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
}


data{
  int<lower=0> n_item;
  int<lower=0> n_response;
  int<lower=0> total_length;
  int<lower=0> item_response[total_length];
}
parameters{
  real theta[n_response]; 
  real beta[n_item];
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

}
model{
  theta ~ cauchy(0,1);
  matrix[N, N] cov =   cov_exp_quad(x, alpha, rho)
                     + diag_matrix(rep_vector(square(sigma), N));
  matrix[N, N] L_cov = cholesky_decompose(cov);

  for (i in 1:n_response){
    for (j in 1:n_item){
    
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
      item_response_rep[n_item *(i - 1) + j] = poisson_log_rng( theta[i] - beta[j]);
    }  
    }
  }
}