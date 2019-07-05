# pyinar.R

#' Model training
#'
#' Computes the posterior distribution
#'
#' @return A model
#' 
#' @export
pyinar <- function(time_series,
				   p = 1,
                   prior = list(a0_alpha = NULL, b0_alpha = NULL,
                                a0_G0 = NULL, b0_G0 = NULL,
                                tau = NULL, k0 = NULL,
								sigma = NULL, 
                                lambda_max = NULL),
                   burn_in = 10^3,
                   chain_length = 10^4,
                   random_seed = 1761,
                   verbose = TRUE)
{
    if (!length(time_series) > 0) stop("Time series must have positive length.")
	if (length(time_series) <= p) stop("Time series length must be bigger than p.")
    if (!burn_in >= 0) stop("Burn-in must be positive.")
    if (!chain_length >= 0) stop("Chain length must be positive.")
    if (!random_seed >= 0) stop("Random seed must be positive.")
    if (!is.logical(verbose)) stop("Verbose parameter must be TRUE or FALSE.")
    
    set.seed(random_seed)
    
    begin <- proc.time()
    
    if (is.null(prior[["a0_G0"]]) || is.null(prior[["b0_G0"]])) {
        if (is.null(prior[["lambda_max"]])) prior[["lambda_max"]] <- max(time_series) # heuristic
        opt_G0 <- .KL_base_measure(prior[["lambda_max"]])
        prior[["a0_G0"]] <- opt_G0[1]
        prior[["b0_G0"]] <- opt_G0[2]
    }	
    if (is.null(prior[["sigma"]])) prior[["sigma"]] <- 0
	if (is.null(prior[["k0"]])) prior[["k0"]] <- 1.01	
	if (is.null(prior[["tau"]])) prior[["tau"]] <- uniroot(.pitman_expected_clusters, 
														   c(-prior[["sigma"]] + 1e-6, 100), prior[["sigma"]], 
				   							               expected_number_of_clusters = prior[["k0"]], 
                   									       t = length(time_series) - p)$root
	if (is.null(prior[["a_alpha"]])) prior[["a_alpha"]] <- rep(1, p)
		                 
    post <- .posterior(time_series,
						  p,
                          prior[c("a_alpha", "a0_G0", "b0_G0", "tau", "sigma")],
                          burn_in,
                          chain_length,
                          random_seed,
                          verbose)
           
    model <- list(time_series = time_series,
				  p = p,
                  prior = prior[c("a_alpha", "a0_G0", "b0_G0", "tau", "sigma")],
                  burn_in = list(length = burn_in,
                                 alpha = post$burn_in$alpha,
                                 lambda = post$burn_in$lambda,
                                 num_clusters = post$burn_in$num_clusters),
                  chain = list(length = chain_length,
                               alpha = post$chain$alpha,
                               lambda = post$chain$lambda,
                               num_clusters = post$chain$num_clusters),
                  est = list(alpha = apply(post$chain$alpha, 2, mean),
                             lambda = apply(post$chain$lambda, 2, mean),
                             num_clusters = which.max(tabulate(post$chain$num_clusters))),
                  elapsed = proc.time() - begin)
    
    class(model) <- "pyinar"
    
    invisible(model)
}

predict.pyinar <- function(model, h = 1, replications = 10^4) {
    if(model$p == 1) {
        pred <- .predictive_distribution_prop(model, h)
        list(est = .generalized_median(pred), distr = pred)
	}
	else {
		pred <- .predictive_distribution_mc(model, h, replications)
        list(est = .generalized_median(pred), distr = pred)
	}
}

#' Model summaries
#'
#' Summarises the model
#'
#' @return A summary
#' 
#' @export
summary.pyinar <- function(model) {
    printf <- function(...) cat(sprintf(...))
    printf("\n========================\n")
    printf("PY-INAR(%d) Model Summary\n", model$p)
    printf("========================\n")
    printf("Time series length: %d\n", length(model$time_series))
    printf("Prior parameters:\n")
    printf("  a0_G0 = %.2f, b0_G0 = %.2f\n", model$prior[["a0_G0"]], model$prior[["b0_G0"]])
    cat("  ")
    for (i in 1:model$p) {
        printf("a_%d = %.2f", i, model$prior[["a_alpha"]][i])
        if (i < model$p) cat(", ")
        else cat("\n")
    }
	printf("  tau = %.2f \n", model$prior[["tau"]])
	printf("  sigma = %.2f \n", model$prior[["sigma"]])
    printf("Effective Markov chains length: %d (%d burn-in)\n", model$chain$length, model$burn_in$length)
    printf("Some posterior means with 95%% credible intervals:\n")
    for (i in 1:model$p) {
        post_qs <- unname(quantile(model$chain$alpha[, i], probs = c(0.025, 0.975)))
        printf("  alpha_%d = %.2f  [ %.2f ; %.2f ]\n", i, mean(model$chain$alpha[, i]), post_qs[1], post_qs[2])
    }
    for (i in as.integer(seq(1, length(model$time_series) - model$p, length.out = 3))) {
        post_qs <- unname(quantile(model$chain$lambda[, i], probs = c(0.025, 0.975)))
        printf("  lambda_%d = %.2f  [ %.2f ; %.2f ]\n", i, mean(model$chain$lambda[, i]), post_qs[1], post_qs[2])
    }
    printf("Posterior distribution of number of clusters:")
    print(table(model$chain$num_clusters) / model$chain$length)
    printf("Total simulation time: %.2f seconds\n\n", round(model$elapsed[3]))
}

.KL_base_measure <- function(lambda_max) {
    D_KL <- function(a0, b0, lambda_max) lgamma(a0) - a0*log(b0) - (a0 - 1)*(log(lambda_max) - 1) + b0*(lambda_max / 2) - log(lambda_max)
    optim(c(1, 1), function(x) D_KL(x[1], x[2], lambda_max), method = "L-BFGS-B", lower = c(1e-6, 1e-6))$par
}
		  
.pitman_expected_clusters <- function (tau, sigma, t, expected_number_of_clusters) {
                                 if (sigma == 0) 
							         tau * (digamma(tau + t) - digamma(tau)) - expected_number_of_clusters 	 
							     else
							         exp((lnpoch(tau + sigma, t) - log(sigma)
                                     - lnpoch(tau + 1, t - 1))) - 
                                     expected_number_of_clusters - (tau / sigma) 
}