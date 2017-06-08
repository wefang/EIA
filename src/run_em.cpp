// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List run_em(const NumericMatrix mat, const NumericMatrix bg_mean, const NumericMatrix bg_sd,
            int K, int max_iter, double tol, int num_threads,
            NumericVector p, NumericMatrix q, NumericVector theta1, NumericVector sigma1,
            const NumericVector eta, const NumericVector gamma,
	    const double lambda = 2.0, const double nu = 2.0, const double kappa = 1.0){
      
      // set number of threads
      omp_set_num_threads(num_threads);

      const int I = mat.nrow(), J = mat.ncol();
      
      bool converge_flag = FALSE;
      double loglike = -1e10, loglike_new;
      double m, d1, d2, temp_max;
      NumericVector p_new(K), theta1_new(J), sigma1_new(J), temp(K), temp_clust_sum(K), temp_post_sum(J), all_like(I), theta1_lower(J);
      NumericMatrix q_new(K, J), like1(I, J), like0(I, J), clust_like(I, K), post(I, J);
      arma::cube temp_like_sum(K, J, I), cond_like(I, K, J);

      // Setting like0 which is not being updated
#pragma omp parallel for shared(like0, mat, bg_mean, bg_sd) collapse(2)
      for (int i=0; i < I; i++){
            for (int j=0; j < J; j++){
                  like0(i, j) = R::dnorm(mat(i, j), bg_mean(i, j), bg_sd(i, j), FALSE);
            }
      }

	// computer the lower constraint for theta1      
#pragma omp parallel for shared(theta1_lower, bg_mean, bg_sd)
      for (int j = 0; j < J; j++){
	      theta1_lower[j] = mean(bg_mean(_, j)) + lambda * mean(bg_sd(_, j));
      }

      for (int iter = 0; iter < max_iter; iter++){
#pragma omp parallel for shared(like1, mat, theta1, sigma1) collapse(2)
          for (int i=0; i < I; i++){
              for (int j=0; j < J; j++){
                  like1(i, j) = R::dnorm(mat(i, j), theta1[j], sigma1[j], FALSE);
              }
          }

#pragma omp parallel for shared(temp_like_sum, cond_like, like0, like1, q) private(d1, d2) collapse(3)
          for (int i = 0; i < I; i++){
              for (int j = 0; j < J; j++){
                  for (int k = 0; k < K; k++){
                      d1 = q(k, j) * like1(i, j);
                      d2 = (1 - q(k, j)) * like0(i, j);
                      temp_like_sum(k, j, i) = d1 + d2;
                      cond_like(i, k, j) = d1 / (d1 + d2);
                  }
              }
          }

          // compute unnormalized clust_like on log scale
#pragma omp parallel for shared(clust_like, temp_like_sum, p)
          for (int i = 0; i < I; i++){
              for (int k = 0; k < K; k++){
                  clust_like(i, k) = log(p[k]) + sum(log(temp_like_sum.slice(i).row(k)));
              }
          }

          // normalize clust_like
          // #pragma omp parallel for shared(clust_like, all_like) private(temp, temp_max)
          for (int i = 0; i < I; i++){
              temp = clust_like(i, _);
              temp_max = max(temp);
              temp = clone(temp) - temp_max;
              temp = exp(temp);
              clust_like(i, _) = temp / sum(temp);
              all_like[i] = log(sum(temp)) + temp_max;
          }

          loglike_new = sum(all_like);
          loglike = loglike_new;

          // compute cond_like
#pragma omp parallel for shared(cond_like, clust_like)
          for (int j = 0; j < J; j++){
              for (int i = 0; i < I; i++){
                  for (int k = 0; k < K; k++){
                      cond_like(i ,k, j) = cond_like(i, k, j) * clust_like(i, k);
                  }
              }
          }

          for (int k =0; k < K; k++){
              temp_clust_sum[k] = sum(clust_like(_, k));
          }
          p_new = (temp_clust_sum + 1) / (I + K);

#pragma omp parallel for shared(cond_like, q_new)
          for (int j = 0; j<J; j++){
              for (int k =0; k < K; k++){
                  q_new(k, j) = (sum(cond_like.slice(j).col(k)) + 1) / (temp_clust_sum[k] + 2);
              }
          }

          if ( max( abs(p_new - p) / p ) < tol & max( abs(q_new - q) / q) < tol ) {
              converge_flag = TRUE;
          }

          p = clone(p_new);
          q = clone(q_new);

#pragma omp parallel for shared(cond_like, post)
          for (int j = 0; j < J; j++){
              for (int i = 0; i < I; i++){
                  post(i, j) = sum(cond_like.slice(j).row(i));
              }
          }

#pragma omp parallel for shared(temp_post_sum, post)
          for (int j = 0; j < J; j++){
              temp_post_sum[j] = sum(post(_, j));
          }

#pragma omp parallel for shared(theta1_new, mat, post, eta, temp_post_sum)
          for (int j = 0; j < J; j++){
              theta1_new[j] = (sum(post(_, j) * mat(_, j)) + kappa*eta[j]) / (temp_post_sum[j]+kappa);
              // restricted maximizer
              if (theta1_new[j] < theta1_lower[j]){
                  theta1_new[j] = theta1_lower[j];
              }
          }
          theta1 = clone(theta1_new);

#pragma omp parallel for shared(post, mat, theta1, nu, gamma, temp_post_sum)
          for (int j = 0; j < J; j++){
              sigma1_new[j] = std::sqrt((sum(post(_, j) * pow(mat(_, j) - theta1[j], 2)) +
				      nu*std::pow(gamma[j], 2) + kappa*std::pow(theta1[j] - eta[j], 2)) /
			      (temp_post_sum[j] + nu + 3));
          }
          sigma1 = clone(sigma1_new);

          if (converge_flag == TRUE){
              Rcout << "converged after " << iter << " iterations. " << "\n";
              break;
          }
      }

      if (converge_flag == FALSE){
	      Rcout << "max iterations reached.\n";
      }

      return List::create(
              Named("p") = p,
              Named("q") = q,
              Named("theta1") = theta1,
              Named("sigma1") = sigma1,
              Named("loglike") = loglike,
              Named("clust.like") = clust_like,
              Named("cond.like") = cond_like
              );

}
