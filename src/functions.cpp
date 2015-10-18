#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]
List SearchMode(NumericMatrix mat, NumericMatrix bg_mean, NumericMatrix bg_sd, int K, int max_iter, double tol,
        NumericVector p, NumericMatrix q, NumericVector theta1, NumericVector sigma1){
    int I = mat.nrow(), J = mat.ncol();

//    NumericVector theta1(J);
//    NumericVector sigma1(J);
    NumericVector theta1_new(J);
    NumericVector sigma1_new(J);

    // Set initial values for theta and sigma
    // for (int j =0; j < J; j++ ){
    //     theta1[j] = mean(na_omit(mat( _ , j ))) + 2 * bg_sd( _ , j);
    //     sigma1[j] = sd(na_omit(mat(_, j )));
    // }

    NumericVector mu1 = clone(theta1);
    NumericVector sigma1_0 = clone(sigma1);
    const double nu_0 = 2.0;

    // Generate initial p and q
 //   NumericVector p(K);
 //   p = rgamma(K, 1, 1);
 //   p = p / sum(p);

 //   NumericMatrix q(K, J);
  //  for (int j=0; j < J; j++){
   //     q( _ , j) = runif(K);
  //  }

    double loglike = -1e10;
    double loglike_new;
    bool converge_flag = FALSE;
    double temp_max, m, s, d1, d2;

    NumericMatrix like1(I, J);
    NumericMatrix like0(I, J);
    arma::cube temp_like_sum(I, J, K);
    arma::cube temp_like_ratio(I, J, K);
    NumericMatrix clust_like(I, K);
    NumericVector temp(K);
    NumericVector all_like(I);
    arma::cube cond_like(I, J, K);
    NumericVector p_new(K);
    NumericMatrix q_new(K, J);
    NumericMatrix post1(I, J);

    // Setting like0 which is not being updated
    for (int i=0; i < I; i++){
        for (int j=0; j < J; j++){
            like0(i, j) = R::dnorm(mat(i, j), bg_mean(i, j), bg_sd(i, j), FALSE);
        }
    }

    // beginning of iteration
    for (int iter = 0; iter < max_iter; iter++){
        for (int j=0; j < J; j++){
            like1(_, j) = dnorm(mat(_, j), theta1[j], sigma1[j]);
        }

        for (int i = 0; i < I; i++){
            for (int j = 0; j < J; j++){
                for (int k = 0; k < K; k++){
                    d1 = q(k, j) * like1(i, j);
                    d2 = (1 - q(k, j)) * like0(i, j);
                    temp_like_sum(i, j, k) = d1 + d2;
                    temp_like_ratio(i, j, k) = d1 / (d1 + d2);
                }
            }
        }

        // compute unnormalized clust_like on log scale
        for (int k = 0; k < K; k++){
            for (int i = 0; i < I; i++){
                s = 0;
                for (int j = 0; j < J; j++){
                    if (!R_IsNA(temp_like_sum(i, j, k))){
                        s += log(temp_like_sum(i, j, k));
                    }
                }
                clust_like(i, k) = log(p[k]) + s;
            }
        }

        // normalize clust_like
        for (int i = 0; i < I; i++){
            temp = clust_like(i, _);
            temp_max = max(temp);
            temp = clone(temp) - temp_max;
            temp = exp(temp);
            clust_like(i, _) = temp / sum(temp);
            all_like[i] = log(sum(temp)) + temp_max;
        }

        loglike_new = sum(all_like);
        // Rcout << (loglike_new - loglike) << std::endl;
        loglike = loglike_new;

        // compute cond_like
        for (int i = 0; i < I; i++){
            for (int j = 0; j < J; j++){
                for (int k = 0; k < K; k++){
                    cond_like(i ,j, k) = temp_like_ratio(i, j, k) * clust_like(i, k);
                }
            }
        }

        for (int k =0; k < K; k++){
            p_new[k] = ( sum(clust_like(_, k)) + 1) / (I + K);
        }

        for (int k =0; k < K; k++){
            m = sum(clust_like(_, k));
            for (int j = 0; j<J; j++){
                s = 0;
                for (int i = 0; i<I; i++){
                    if (!R_IsNA(cond_like(i, j, k))){
                        s += cond_like(i, j, k);
                    }
                }
                q_new(k, j) = (s + 1) / ( m + 2);
            }
        }

        if ( max( abs(p_new - p) / p ) < tol & max( abs(q_new - q) / q) < tol ){
            converge_flag = TRUE;
        }

        p = clone(p_new);
        q = clone(q_new);

        for (int i = 0; i < I; i++){
            for (int j = 0; j < J; j++){
                s = 0;
                for (int k = 0; k < K; k++){
                    if (!R_IsNA(cond_like(i, j, k))){
                        s += cond_like(i, j, k);
                    }
                }
                post1(i, j) = s;
            }
        }

        for (int j = 0; j < J; j++){
            s = 0;
            m = 0;
            for (int i = 0; i < I; i++){
                if (!R_IsNA(mat(i, j))){
                    s += post1(i, j) * mat(i, j);
                }
                m += post1(i, j);
            }
            theta1_new[j] = ( s + mu1[j] ) / (m + 1);
            if (theta1_new[j] < max(bg_mean( _ , j)) + max(bg_sd( _ , j))){
                theta1_new[j] = max(bg_mean( _ , j)) + max(bg_sd( _ , j));
            }
        }

        theta1 = clone(theta1_new);

        for (int j = 0; j < J; j++){
            s = 0;
            m = 0;
            for (int i = 0; i < I; i++){
                if (!R_IsNA(mat(i, j))){
                    s += post1(i, j) * std::pow( ( mat(i, j) - theta1[j] ), 2) ;
                }
                m += post1(i, j);
            }
            sigma1_new[j] = std::sqrt( (s + std::pow(nu_0, 2) * std::pow(sigma1_0[j], 2))  / ( m + nu_0 - 1 ));
        }
        sigma1 = clone(sigma1_new);

        if (converge_flag == TRUE){
            Rcout << "converged after " << iter << " iterations. " << std::endl;
            break;
        }

    }


    return List::create(
            Rcpp::Named("p") = p,
            Rcpp::Named("q") = q,
            Rcpp::Named("theta1") = theta1,
            Rcpp::Named("sigma1") = sigma1,
            Rcpp::Named("loglike") = loglike,
            Rcpp::Named("clust.like") = clust_like,
            Rcpp::Named("cond.like") = cond_like
            );

}

// [[Rcpp::export]]
NumericMatrix calLikeIso(NumericMatrix mat, NumericVector p, NumericMatrix q, NumericMatrix bg_mean, NumericMatrix bg_sd,
        NumericVector theta1, NumericVector sigma1){

    int I = mat.nrow(), J = mat.ncol(), K = q.nrow();

    NumericMatrix like1(I, J);
    NumericMatrix like0(I, J);
    arma::cube temp_like_sum(I, J, K);
    double s;

    for (int i=0; i < I; i++){
        for (int j=0; j < J; j++){
            like0(i, j) = R::dnorm(mat(i, j), bg_mean(i, j), bg_sd(i, j), FALSE);
            like1(i, j) = R::dnorm(mat(i, j), theta1[j], sigma1[j], FALSE);
        }
    }

    for (int i = 0; i < I; i++){
        for (int j = 0; j < J; j++){
            for (int k = 0; k < K; k++){
                temp_like_sum(i, j, k) = q(k, j) * like1(i, j) + (1 - q(k, j)) * like0(i, j);
            }
        }
    }

    NumericMatrix loglike_mat(I, K);
    NumericVector temp(K);

    for (int k = 0; k < K; k++){
        for (int i = 0; i < I; i++){
            s = 0;
            for (int j = 0; j < J; j++){
                s += log(temp_like_sum(i, j, k));
            }
            loglike_mat(i, k) = log(p[k]) + s;
        }
    }

//    for (int i = 0; i < I; i++){
//        temp = loglike_mat(i, _);
//        loglike_mat(i, _) = temp - max(temp);
//    }

    return loglike_mat;
}

// [[Rcpp::export]]
NumericVector count_bins(NumericVector counts, NumericVector bins){
  int i;
  for (i=0; i<bins.size(); i++){ counts[bins[i] - 1]++; }
  return counts;
}


