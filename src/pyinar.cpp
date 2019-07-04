// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace std;
using namespace Rcpp;

void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

double lfactorial(double x) {
  return lgamma(x + 1);
}

double rtbeta(double a0, double b0, double max) {
    double upper_limit = R::pbeta(max, a0, b0, TRUE, FALSE);
    double u = R::runif(0, upper_limit);
    return R::qbeta(u, a0, b0, TRUE, FALSE);
}

IntegerVector minus_m_ij(IntegerMatrix x, int i, int j) {
    int n = x.ncol();
     IntegerVector minus_matrix(n-1);
     for (int t = 0; t < n; t++) {
          if (t < j) minus_matrix[t] =  x(i, t);
          else if (t > j) minus_matrix[t - 1] = x(i, t);
     }
     return minus_matrix;
}

IntegerVector table_each_count(NumericVector x) {
    int t = x.size(); 
    IntegerVector n(t); 
  
    for (int i = 0; i < t; i++) { 
        double z = x[i]; 
        NumericVector y = x[x == z]; 
        n[i] = y.size();
    }
	
    return n;
}

int m_full_conditional(double alpha, double lambda, int yt_i, int yt, IntegerVector m_minus_t) {
    int M = m_minus_t.size(), sum_m_minus_t;

    sum_m_minus_t = 0; for(int j = 0; j < M; j++) sum_m_minus_t += m_minus_t[j];

    if (yt_i == 0 || yt - sum_m_minus_t == 0) return 0;

    NumericVector log_weights(min(yt_i, yt - sum_m_minus_t) + 1), weights(min(yt_i, yt - sum_m_minus_t) + 1);
    int m_t;
    double total = 0;

    for (m_t = 0; m_t <= min(yt_i, yt - sum_m_minus_t); m_t++) {
        log_weights[m_t] = m_t*(log(alpha) - log(lambda) - log(1 - alpha))
                           - lfactorial(m_t) - lfactorial(yt_i - m_t) - lfactorial(yt - m_t - sum_m_minus_t);
        /* if (!R_finite(log_weights[m_t])) {
            Rprintf("alpha = %f, lambda = %f\n, m_t = %f\n", alpha, lambda);
            stop("m_full_conditional: log weight is not a finite number!\n");
        } */
    }

    double max_log_weights = max(log_weights);

    for (m_t = 0; m_t <= min(yt_i, yt - sum_m_minus_t); m_t++) {
        if (R_finite(log_weights[m_t])) weights[m_t] = exp(log_weights[m_t] - max_log_weights);
        else                            weights[m_t] = 0;
        total += weights[m_t];
    }

    double u = R::runif(0, 1);

    for (m_t = 0; m_t <= min(yt_i, yt - sum_m_minus_t); m_t++) {
        if ((u -= (weights[m_t] / total)) <= 0) return m_t;
    }

    stop("lambda_full_conditional: should never get to this point!\n");
}

double lambda_full_conditional_py(NumericVector minus_lambda_t, double a0_G0, double b0_G0, double tau, double sigma, int k_minus_lambda_t, IntegerVector n_minus_t, int y_t, int sum_m) {
    int J = minus_lambda_t.size(); 
    NumericVector log_weights(J), weights(J);
    double weight_gamma, total = 0; 
  
    for(int j = 0; j < J; j++) {
        log_weights[j] = - minus_lambda_t[j] + (y_t - sum_m)*log(minus_lambda_t[j]) ;
        /* if (! R_finite(log_weights[r])) {
          Rprintf("minus_lambda_t[%d] = %f, log_weights[%d] = %f\n", r, minus_lambda_t[r], r, log_weights[r]);
          stop("lambda_full_conditional: log weight is not a finite number!\n");
        } */
    }
  
    weight_gamma = exp( log(tau + k_minus_lambda_t*sigma) + a0_G0*log(b0_G0) + lgamma(y_t - sum_m + a0_G0) - lgamma(a0_G0) - (y_t - sum_m + a0_G0)*log(b0_G0 + 1) ); 	
  
    total += weight_gamma; 
  
    for (int r = 0; r < J; r++) {
        if (R_finite(log_weights[r])) weights[r] = exp(log_weights[r] + log( 1 - ( sigma / n_minus_t[r] ) ) );
        else                          weights[r] = 0;
        /* if (! R_finite(weights[r])) stop("lambda_full_conditional: weight is not a finite number!\n"); */ 
        total += weights[r];
    }
  
    /* if (! R_finite(total)) stop("lambda_full_conditional: total sum of weights is not a finite number!\n"); */ 
  
    double u = R::runif(0, 1);
  
    if((u -= (weight_gamma / total)) <= 0) {
        return R::rgamma(a0_G0 + y_t - sum_m , 1 / (b0_G0 + 1));
    } else {
        for (int r = 0; r < J; r++) {
          if ((u -= (weights[r] / total)) <= 0) return minus_lambda_t[r];
        }
    }
  
    // Rprintf("u = %f\n", u);
    stop("lambda_full_conditional: should never get to this point!\n");
}

double lambda_star_full_conditional(int j, IntegerVector y, IntegerMatrix m, IntegerVector c, double a0_G0, double b0_G0) {
    int sum = 0, n_j = 0, p = m.ncol(), T = y.size();

    for (int t = p; t < T; t++) {
        int sum_m = 0;
        if (c[t - p] == j + 1) {
            for (int r = 0; r < p; r++) sum_m += m(t - p , r);
            sum += y[t] - sum_m;
            n_j++;
        }
    }

    return R::rgamma(a0_G0 + sum, 1 / (b0_G0 + n_j) );
}

NumericVector minus_alpha_ij(NumericMatrix x, int i, int j) {
    int n = x.ncol();
    NumericVector minus_vector(n - 1);

    for (int t = 0; t < n; t++) {
        if (t < j) minus_vector[t] =  x(i, t) ;
        else if (t > j) minus_vector[t - 1] = x(i - 1, t) ;
    }

    return minus_vector;
}

// [[Rcpp::export(.posterior)]]
List posterior(IntegerVector y,  
			   int p,
			   List prior, 
			   int burn_in, 
			   int N,
			   unsigned int random_seed, 
			   bool verbose) 
{	
    int T = y.size(), sum_m, k_minus_lambda_t;
    double sum_minus_alpha;
	NumericVector minus_lambda_t(T - p - 1), minus_alpha(p - 1), sum_1(p), sum_2(p);
    IntegerVector clusters(burn_in + N), n_minus_t(T - p - 1), minus_m(T - p - 1); 
    NumericMatrix alpha(burn_in + N, p), lambda(burn_in + N, T - p);
    IntegerMatrix m(T - p, p);
    Progress pb(burn_in + N + 1, true);
	List model;
	
	set_seed(random_seed); 
	
	NumericVector a = prior["a_alpha"]; 
	double a0_G0 = prior["a0_G0"]; 
	double b0_G0 = prior["b0_G0"]; 
	double tau = prior["tau"]; 
	double sigma = prior["sigma"]; 
	
	for (int j = 0; j < p; j++) alpha(0, j) = 0;
    for (int t = 0; t < T - p; t++)  lambda(0, t) = 1; 
  
    pb.increment();
    for (int i = 1; i < burn_in + N; i++) {
        if (Progress::check_abort()) return model;

		 for (int j = 0; j < p; j++) {
            sum_1[j] = 0; for (int t = 0; t < T - p; t++) sum_1[j] += m(t, j);
            sum_2[j] = 0; for (int t = 0; t < T - p; t++) sum_2[j] += y[t + p - j - 1] - m(t, j);
            minus_alpha = minus_alpha_ij(alpha, i, j);
            sum_minus_alpha = 0; for (int ell = 0; ell < p - 1; ell++) sum_minus_alpha += minus_alpha[ell];
            alpha(i, j) = rtbeta(a[j] + sum_1[j], 1 + sum_2[j], 1 - sum_minus_alpha);
        }
        for (int t = 0; t < T - p - 1; t++)  n_minus_t[t] = 0; 
     
        for (int t = 0; t < T - p; t++) {
            for (int s = 0; s < T - p; s++) {
                if (s < t) {
                minus_lambda_t[s] = lambda(i, s);
                } else if (s > t) {
                minus_lambda_t[s-1] = lambda(i-1, s);
                }
            }
            NumericVector minus_lambda_t_star = unique(minus_lambda_t);
            k_minus_lambda_t = minus_lambda_t_star.size();
            n_minus_t = table_each_count(minus_lambda_t); 
			sum_m = 0; for(int j = 0; j < p; j++)  sum_m += m(t, j);
            lambda(i, t) = lambda_full_conditional_py(minus_lambda_t, a0_G0, b0_G0, tau, sigma, k_minus_lambda_t, n_minus_t, y[t + p], sum_m);
        }
    
       for (int j = 0; j < p; j++) {
            for (int t = 0; t < T - p; t++) {
                minus_m = minus_m_ij(m, t, j);
                m(t, j) = m_full_conditional(alpha(i, j), lambda(i, t), y[t + p - j - 1], y[t + p], minus_m);
            }
        }
        NumericVector lambda_star = unique(lambda.row(i));
        int k = lambda_star.size();
        clusters[i] = k;
        IntegerVector c = match(lambda.row(i), lambda_star);
    
        for (int j = 0; j < k; j++) lambda_star[j] = lambda_star_full_conditional(j, y, m, c, a0_G0, b0_G0);
    
        for (int t = 0; t < T - p; t++) lambda(i, t) = lambda_star[c[t] - 1];
    
        pb.increment();
    
    }
  
    model["time_series"] = y;
	model["p"] = p;
	
	List burn_in_pars;
    burn_in_pars["alpha"] = alpha(Range(0, burn_in - 1), _);
    burn_in_pars["lambda"] = lambda(Range(0, burn_in - 1), _);
    burn_in_pars["num_clusters"] = clusters[Range(0, burn_in - 1)];
    model["burn_in"] = burn_in_pars;

	List chain_pars;
    chain_pars["alpha"] = alpha(Range(burn_in, burn_in + N - 1), _);
    chain_pars["lambda"] = lambda(Range(burn_in, burn_in + N - 1), _);
    chain_pars["num_clusters"] = clusters[Range(burn_in, burn_in + N - 1)];
    model["chain"] = chain_pars;
	
    return model;
}

double polya_blackwell_macqueen_py(double a0_G0, double b0_G0, double tau, double sigma, NumericVector lambda) {
    int t = lambda.size();
    double u = R::runif(0, 1);
    NumericVector lambda_star = unique(lambda);
    int k = lambda_star.size();
    IntegerVector n = table_each_count(lambda); 
  
    if((u -= ((tau + k*sigma) / (tau + t))) <= 0) {
      return R::rgamma(a0_G0, 1 / b0_G0);
    } else {
      for (int r = 0; r < t; r++) {
          if ((u -= ( (1 - sigma / n[r]) / (tau + t) ) )  <= 0) return lambda[r];
          }
      }
}

 double pr_y_fut_given_yt_alpha_lambda(int y_fut, int yt, double alpha, NumericVector lambda) {
    int k = lambda.size();
    double mu_k = 0;
    double sum = 0;
  
    for(int i = 0 ; i < k ; i++) mu_k += exp( (k-i-1)*log(alpha) + log(lambda[i]) );
  
    for (int m = 0; m <= min(yt, y_fut); m++)
     sum += exp( (y_fut-m)*log(mu_k) - mu_k - lfactorial(y_fut-m)
                  + R::lchoose(yt, m) + k*m*log(alpha) + (yt-m)*log(1-exp( k*log(alpha) )) );
    
    return sum;
}

// [[Rcpp::export(.generalized_median)]]
int generalized_median(NumericVector pred) {
    int median = 0;
    double cum_sum = pred[0];
    double least_abs_diff = abs(0.5 - cum_sum);
    for (int i = 1; i < pred.size(); i++) {
        cum_sum += pred[i];
        double abs_diff = abs(0.5 - cum_sum);
        if (abs_diff < least_abs_diff) {
            median = i;
            least_abs_diff = abs_diff;
        }
    }
    return median;
}

// [[Rcpp::export(.predictive_distribution_prop)]]
NumericVector predictive_distribution_prop(List model, int h) {
    IntegerVector y = model["time_series"];
    int T = y.size();
    int p = model["p"];
    List chain = model["chain"];
    NumericVector alpha = chain["alpha"];
    NumericMatrix lambda_prev = chain["lambda"];
    int N = chain["length"];
    List prior = model["prior"];
    double a0_G0 = prior["a0_G0"], b0_G0 = prior["b0_G0"], sigma = prior["sigma"], tau = prior["tau"];
	
	int count = 0;
    double cum = 0, prob = 0, mean;
    const double TOL = 0.05;
    NumericMatrix lambda(N, T - 1 + h);
    NumericVector proxy_lambda(h), pred(100);
  
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < T-1; j++) lambda(i, j) = lambda_prev(i, j);
        // Successive Polya-Blackwell-MacQueen urns
        for (int r = 1; r <= h; r++) {
            NumericVector given_lambda(T-2+r);
            for (int s = 0; s < T - 2 + r; s++) given_lambda[s] = lambda(i, s);
            lambda(i, T - 2 + r) = polya_blackwell_macqueen_py(a0_G0, b0_G0, tau, sigma, given_lambda);
        }
    }
  
    do {
        mean = 0;
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < h; j++) proxy_lambda[j] = lambda(i, T - 1 + j);
            mean += pr_y_fut_given_yt_alpha_lambda(count, y[T-1], alpha[i], proxy_lambda);
        }
        mean /= N;
        pred[count] = mean; 
        cum += mean;
        count++;
    }
	while (cum < (1-TOL));
  
    return pred;	
}

// [[Rcpp::export(.predictive_distribution_mc)]]
NumericVector predictive_distribution_mc(List model, int h, int replications) {
    IntegerVector y = model["time_series"];
    int T = y.size();
    int p = model["p"];
    List chain = model["chain"];
    NumericMatrix alpha = chain["alpha"];
    NumericMatrix lambda_prev = chain["lambda"];
    int N = chain["length"];
    List prior = model["prior"];
    double a0_G0 = prior["a0_G0"], b0_G0 = prior["b0_G0"], sigma = prior["sigma"], tau = prior["tau"];

    int y_next = 0;
    NumericMatrix lambda(N, T - p + h);
    IntegerVector y_prev(T + h);
    IntegerVector yt_plus_h(N * replications);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < T - p; j++) lambda(i, j) = lambda_prev(i, j);
        // Successive Polya-Blackwell-MacQueen urns
        for (int r = 1; r <= h; r++) {
            NumericVector given_lambda(T - p - 1 + r);
            for (int s = 0; s < T - p - 1 + r; s++) given_lambda[s] = lambda(i, s);
            lambda(i, T - p - 1 + r) = polya_blackwell_macqueen_py(a0_G0, b0_G0, tau, sigma, given_lambda);
        }
    }

    for(int j = 0; j < T; j++) y_prev[j] = y[j];

    for (int i = 0; i < N; i++) {
        for (int ell = 0; ell < replications; ell++) {
            for (int j = 1; j <= h; j++) {
                y_next = 0;
                for(int r = 0; r < p; r++) y_next += R::rbinom(y_prev[T - p - 1 + r + j], alpha(i, p - r - 1));
                y_next += R::rpois(lambda(i, T - p - 1 + j));
                y_prev[T + j - 1] = y_next;
            }
            yt_plus_h[i*replications + ell] = y_next;
        }
    }

    NumericVector count = Rcpp::as<NumericVector>(table(yt_plus_h));
    int c = count.size();
    NumericVector pred(c);
    for (int i = 0; i < c; i++) pred[i] = count[i] / (N * replications);

    return pred;
}

double log_p_y_t_given_m_t_lambda_t(int y_t, int m_t, double lambda_t) {
    return ((y_t - m_t)*log(lambda_t) - lambda_t - lfactorial(y_t - m_t));
}

double dawid_sebastiani_score(NumericVector pred, int y) { // smaller is better
    int c = pred.size();

    double mu_pred = 0;
    for (int i = 1; i < c; i++) mu_pred += (pred[i] * i);

    double mean_of_squares = 0;
    for (int i = 1; i < c; i++) mean_of_squares += (pred[i] * pow(i, 2));
    double sd_pred = sqrt(mean_of_squares - pow(mu_pred, 2));
	
    return pow(((y - mu_pred) / sd_pred), 2) + 2 * log(sd_pred);
}
