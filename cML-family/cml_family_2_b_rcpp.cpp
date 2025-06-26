// 首先是Rcpp的C++代码部分
// 保存为 mle_optimization.cpp

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// 计算t_m和f_m统计量
// [[Rcpp::export]]
List compute_statistics_cpp(const arma::mat& beta_hat, 
                           const arma::mat& gamma_current, 
                           double alpha_current,
                           const List& Sigma_beta_list) {
  int num_snps = Sigma_beta_list.size();
  arma::vec var_beta_o(num_snps);
  arma::vec var_beta_p(num_snps);
  arma::vec t_m(num_snps);
  arma::vec f_m(num_snps);
  
  for (int m = 0; m < num_snps; m++) {
    arma::mat sigma_beta = as<arma::mat>(Sigma_beta_list[m]);
    var_beta_o(m) = sigma_beta(0, 0);
    var_beta_p(m) = sigma_beta(1, 1);
    
    double diff_o = beta_hat(0, m) - alpha_current * gamma_current(0, m);
    double diff_p = beta_hat(1, m) - alpha_current * gamma_current(1, m);
    
    t_m(m) = diff_o * diff_o / var_beta_o(m);
    f_m(m) = diff_p * diff_p / var_beta_p(m);
  }
  
  return List::create(
    Named("t_m") = t_m,
    Named("f_m") = f_m,
    Named("var_beta_o") = var_beta_o,
    Named("var_beta_p") = var_beta_p
  );
}

// 更新A∩B集合中的gamma
// [[Rcpp::export]]
arma::mat update_gamma_intersect_cpp(const arma::uvec& A_intersect_B,
                                    const arma::mat& gamma_hat,
                                    const arma::mat& beta_hat,
                                    double alpha_current,
                                    const List& Sigma_gamma_inv_list,
                                    const List& Sigma_beta_inv_list,
                                    arma::mat gamma_new) {
  
  for (uword i = 0; i < A_intersect_B.n_elem; i++) {
    int m = A_intersect_B(i) - 1; // R索引转C++索引
    
    arma::mat term1_inv = as<arma::mat>(Sigma_gamma_inv_list[m]);
    arma::mat term2_inv = alpha_current * alpha_current * as<arma::mat>(Sigma_beta_inv_list[m]);
    arma::mat combined_inv = term1_inv + term2_inv;
    
    arma::vec term1_val = term1_inv * gamma_hat.col(m);
    arma::vec term2_val = alpha_current * as<arma::mat>(Sigma_beta_inv_list[m]) * beta_hat.col(m);
    
    gamma_new.col(m) = solve(combined_inv, term1_val + term2_val);
  }
  
  return gamma_new;
}

// 更新A-B集合中的gamma
// [[Rcpp::export]]
arma::mat update_gamma_diff_cpp(const arma::uvec& A_diff_B,
                               const arma::mat& gamma_hat,
                               const arma::mat& beta_hat,
                               double alpha_current,
                               const List& Sigma_gamma_list,
                               const arma::vec& var_beta_o,
                               arma::mat gamma_new) {
  
  for (uword i = 0; i < A_diff_B.n_elem; i++) {
    int m = A_diff_B(i) - 1; // R索引转C++索引
    
    arma::mat sigma_gamma = as<arma::mat>(Sigma_gamma_list[m]);
    
    double num_go = gamma_hat(0, m) * var_beta_o(m) - 
                   alpha_current * sigma_gamma(0, 0) * beta_hat(0, m);
    double den_go = var_beta_o(m) - alpha_current * alpha_current * sigma_gamma(0, 0);
    double gamma_o_new = num_go / den_go;
    
    double term_gp = sigma_gamma(0, 1) * (-alpha_current / var_beta_o(m)) * 
                    (beta_hat(0, m) - alpha_current * gamma_o_new);
    double gamma_p_new = gamma_hat(1, m) + term_gp;
    
    gamma_new(0, m) = gamma_o_new;
    gamma_new(1, m) = gamma_p_new;
  }
  
  return gamma_new;
}

// 更新alpha的分子分母
// [[Rcpp::export]]
List update_alpha_terms_cpp(const arma::uvec& A_intersect_B,
                           const arma::uvec& A_diff_B,
                           const arma::mat& gamma_current,
                           const arma::mat& beta_hat,
                           const List& Sigma_beta_inv_list,
                           const arma::vec& var_beta_o) {
  
  double num_alpha = 0.0;
  double den_alpha = 0.0;
  
  // 处理A∩B集合
  for (uword i = 0; i < A_intersect_B.n_elem; i++) {
    int m = A_intersect_B(i) - 1; // R索引转C++索引
    
    arma::mat sigma_beta_inv = as<arma::mat>(Sigma_beta_inv_list[m]);
    arma::vec gamma_m = gamma_current.col(m);
    arma::vec beta_m = beta_hat.col(m);
    
    num_alpha += as_scalar(gamma_m.t() * sigma_beta_inv * beta_m);
    den_alpha += as_scalar(gamma_m.t() * sigma_beta_inv * gamma_m);
  }
  
  // 处理A-B集合
  for (uword i = 0; i < A_diff_B.n_elem; i++) {
    int m = A_diff_B(i) - 1; // R索引转C++索引
    
    num_alpha += gamma_current(0, m) * beta_hat(0, m) / var_beta_o(m);
    den_alpha += gamma_current(0, m) * gamma_current(0, m) / var_beta_o(m);
  }
  
  return List::create(
    Named("num_alpha") = num_alpha,
    Named("den_alpha") = den_alpha
  );
}

// 计算矩阵的最大绝对变化
// [[Rcpp::export]]
double compute_max_change_cpp(const arma::mat& mat1, const arma::mat& mat2) {
  return max(max(abs(mat1 - mat2)));
}