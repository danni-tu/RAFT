# RAFT algorithm

#' RAFT Penalized Loss Function for an Individual
#'
#' This function calculates the penalized loss function for an individual based on the specified warp parameters, penalty settings, and likelihood model.
#'
#' @param bblid_i String: the identifier for the individual whose data is being transformed.
#' @param reg_type Character string: the type of warp transformation. Options are "beta_cdf" or "poly_warp".
#' @param warp_params_i Numeric vector: warp parameters for the specified individual.
#' @param pfr_obj Object: the fitted model object, if available. Default is NULL.
#' @param idvar String: the name of the variable in `dat_orig_obs` representing the IDs.
#' @param outcome_var String: the name of the variable in `dat_orig_obs` representing the outcome variable.
#' @param dat_orig_obs data.frame containing the data to be warped, with one observation per row.
#' @param argvals_orig_obs String: the name of the variable in `dat_orig_obs` representing the original time grid.
#' @param fvars Character vector: the names of the functional covariates.
#' @param penalty List: a list specifying which penalties to apply. Default is list(loglik = TRUE, dist_temp = TRUE, dist_orig = FALSE, dist_alphas = FALSE, dist_warp = TRUE).
#' @param lambdas List: a list specifying the penalty weights. Default is list(loglik = 1, dist_temp = 1, dist_orig = 1, dist_alphas = 1, dist_warp = 1).
#' @param penalize_type Character string: the type of penalization. Options are "dist" for distance-based penalization or "corr" for correlation-based penalization. Default is "dist".
#' @param template_ind Numeric scalar: the index of the template curve. Default is 5.
#' @param scale_init Numeric scalar: the initial scale parameter for the likelihood model. Default is NULL.
#'
#' @return A numeric scalar representing the penalized loss function for the individual.
#'
#' @export
loglik_pen_i <- function(bblid_i, reg_type, warp_params_i, pfr_obj = NULL,
                         idvar, outcome_var,  dat_orig_obs, argvals_orig_obs,
                         fvars,
                         penalty = list(loglik = TRUE, dist_temp = TRUE, dist_orig = FALSE, dist_alphas = FALSE, dist_warp = TRUE),
                         lambdas = list(loglik = 1, dist_temp = 1, dist_orig = 1, dist_alphas = 1, dist_warp = 1),
                         penalize_type = "dist",
                         template_ind,
                         scale_init = NULL){

  if (is.null(template_ind)){
    template_ind = 5
  }
  stopifnot(all(fvars %in% colnames(dat_orig_obs)))

  # Given warp_params_i, warp the curve for bblid_i
  dat_warped_fvars_i = warp_data_one(dat_orig_obs = dat_orig_obs, reg_type = reg_type,
                                     bblid_i = bblid_i, idvar = idvar,
                                     warp_params_i = warp_params_i, fvars = fvars)

  # Replace the data for person i with the warped data
  dat_orig_i = dat_orig_obs[which(dat_orig_obs[,idvar] == bblid_i),]
  dat_warped_i = dat_orig_i

  for (j in seq_along(fvars)){
    dat_warped_i[[fvars[j]]] <- dat_warped_fvars_i[[j]]
  }

  # Initialize values of the loss function
  loss_raw = list(loglik = 0, dist_temp = 0, dist_orig = 0, dist_alphas = 0, dist_warp = 0)

  # Maximize the likelihood of the outcome given the pfr model and warped covariates
  if (is.null(pfr_obj)|(!penalty$loglik)){
    loss_raw$loglik = 0
  } else {
    stopifnot(all(paste0("L.", fvars) %in% colnames(pfr_obj$model)))

    # Get predicted values from pfr_obj
    pred_mean = refund:::predict.pfr(object = pfr_obj, newdata = dat_warped_i, type = "response")
    pred_mean = as.vector(pred_mean)
    if (!is.null(scale_init) & is.numeric(scale_init)){ # User-supplied or estimated scale
      sigma_y_hat = scale_init
    } else {
      sigma_y_hat = pfr_obj$scale
    }

    # Log-likelihood of person i
    lik_i = dnorm(x = dat_warped_i[,outcome_var], mean = pred_mean, sd = sigma_y_hat)
    # For numeric stability, bound away from 0
    lik_i[dplyr::near(lik_i, 0)] <- .Machine$double.eps
    loss_raw$loglik =  sum(log(lik_i))
  }

  # Penalize distance between warped and template curve
  if (penalty$dist_temp){
    dist_to_template = rep(as.numeric(NA), length(fvars))
    for(j in seq_along(fvars)){

      dat_curves_ij = as.vector(dat_warped_i[[fvars[j]]])
      dat_temp_ij = dat_orig_obs[[fvars[j]]][template_ind,]

      dist_to_template[j] <- if (penalize_type == "dist"){
        # Euclidean distance between curve i and the template
        sqrt(sum((dat_curves_ij - dat_temp_ij)^2))
      } else {
        # 1 - correlation between the curve and the template
        1-cor(dat_curves_ij, dat_temp_ij)
      }
    }
    loss_raw$dist_temp = mean(dist_to_template, na.rm = TRUE)
  } else {
    loss_raw$dist_temp = 0
  }

  # Penalize distance between warped and original curve
  if (penalty$dist_orig){

    dist_to_orig = rep(as.numeric(NA), length(fvars))

    for(j in seq_along(fvars)){
      dat_curves_ij = as.vector(dat_warped_i[[fvars[j]]])
      dat_orig_ij = as.vector(dat_orig_i[[fvars[j]]])

      dist_to_orig[j] = if (penalize_type == "dist"){
        # Euclidean distance between curve i and the template
        sqrt(sum((dat_curves_ij - dat_orig_ij)^2))
        # mean((dat_curves_ij -dat_orig_ij)^2)
      } else {
        1-cor(dat_curves_ij, dat_orig_ij)
      }
    }

    loss_raw$dist_orig = mean(dist_to_orig)
  } else {
    loss_raw$dist_orig = 0
  }

  # Penalize distance between alpha parameters and (1,1)
  if (penalty$dist_alphas){
    loss_raw$dist_alphas = sqrt(sum((warp_params_i - c(1,1))^2))
  } else {
    loss_raw$dist_alphas = 0
  }

  # Penalize distance between warping function and identity warp
  if (penalty$dist_warp){
    n_t = ncol(dat_orig_obs[[argvals_orig_obs]])
    warp_identity = warp_params_to_curves_inner(warp_params_i, n_grid = n_t, reg_type = reg_type)$x
    warp_fun_i = warp_params_to_curves_inner(warp_params_i, n_grid = n_t, reg_type = reg_type)$y

    loss_raw$dist_warp = sqrt(sum((warp_identity - warp_fun_i)^2))
  } else {
    loss_raw$dist_warp = 0
  }

  # Loss function weighted by the lambdas
  loss = map2(.x = loss_raw, .y = lambdas, ~(.x)*(.y))
  loss = map2(.x = loss, .y = c(1, -1, -1, -1, -1), ~(.x)*(.y))
  loss_value = sum(unlist(loss))

  attr(loss_value, "loss") <- loss
  attr(loss_value, "loss_raw") <- loss_raw
  attr(loss_value, "penalty") <- penalty
  attr(loss_value, "lambdas") <- lambdas
  attr(loss_value, "reg_type") <- reg_type
  return(loss_value)
}


#' Maximize the loss function for an individual
#'
#' This function maximizes the penalized likelihood function \code{\link[loglik_pen_i]{loglik_pen_i}} for an individual to estimate the warp parameters that best fit the data.
#'
#' @param bblid_i String: the identifier for the individual whose data is being transformed.
#' @param reg_type Character string: the type of warp transformation. Options are "beta_cdf" or "poly_warp".
#' @param warp_params_i Numeric vector: warp parameters for the specified individual.
#' @param pfr_obj Object: the fitted `pfr` model object, if available. Default is NULL.
#' @param fvars Character vector: the names of the functional covariates.
#' @param idvar String: the name of the variable in `dat_orig_obs` representing the IDs.
#' @param outcome_var String: the name of the variable in `dat_orig_obs` representing the outcome variable.
#' @param dat_orig_obs data.frame containing the data to be warped, with one observation per row.
#' @param argvals_orig_obs String: the name of the variable in `dat_orig_obs` representing the original time grid.
#' @param penalty List: a list specifying which penalties to apply. Default is list(loglik = TRUE, dist_temp = TRUE, dist_orig = FALSE, dist_alphas = FALSE, dist_warp = TRUE).
#' @param lambdas List: a list specifying the penalty weights. Default is list(loglik = 1, dist_temp = 1, dist_orig = 1, dist_alphas = 1, dist_warp = 1).
#' @param penalize_type Character string: the type of penalization. Options are "dist" for distance-based penalization or "corr" for correlation-based penalization. Default is "dist".
#' @param template_ind Numeric scalar: the index of the template curve. Default is 5.
#' @param opt_it Numeric scalar: the number of iterations for optimization. Default is 50.
#' @param scale_init Numeric scalar: the initial scale parameter for the likelihood model. Default is NULL.
#'
#' @return A list containing the optimized warp parameters, the loss value, the raw loss value, and the runtime.
#'
#' @export
max_loglik_i <- function(bblid_i, warp_params, reg_type, pfr_obj = NULL, fvars = c("mod_fun_dmn7"),
                         idvar = "bblid", outcome_var = "y1",
                         dat_orig_obs, argvals_orig_obs = "density_grid",
                         penalty = list(loglik = TRUE, dist_temp = TRUE, dist_orig = FALSE, dist_alphas = FALSE, dist_warp = TRUE),
                         lambdas = list(loglik = 1, dist_temp = 1, dist_orig = 1, dist_alphas = 1, dist_warp = 1),
                         penalize_type,
                         template_ind,
                         opt_it = 50, scale_init = NULL){

  # Initialize at the current warp_params_i
  warp_params_i = warp_params[[which(names(warp_params) == bblid_i)]]

  t0 = Sys.time()
  lower = switch(reg_type,
                 "beta_cdf" = c(0.2, 0.2),
                 "poly" = 0.1)
  upper = switch(reg_type,
                 "beta_cdf" = c(3.5,3.5),
                 "poly" = 5)

  opt = optim(par = warp_params_i,
              fn = loglik_pen_i,
              lower = lower,
              upper = upper,
              method = "L-BFGS-B",
              control = list(fnscale = -1),
              # Other function inputs
              bblid_i = bblid_i,
              reg_type = reg_type,
              idvar = idvar,
              outcome_var = outcome_var,
              pfr_obj = pfr_obj,
              fvars = fvars,
              dat_orig_obs = dat_orig_obs,
              argvals_orig_obs = argvals_orig_obs,
              penalty = penalty,
              lambda = lambdas,
              penalize_type = penalize_type,
              template_ind = template_ind,
              scale_init = scale_init)

  if (opt$convergence == 52){
    n_deps = if (reg_type == "poly") 0.1 else  c(0.1, 0.1)
    opt = optimr(par = warp_params_i,
                 fn = loglik_pen_i,
                 lower = lower,
                 upper = upper,
                 # OPTIMR settings
                 method = "Rvmmin",
                 control = list(
                   ndeps = n_deps,
                   maxit = 500,
                   maximize = TRUE),
                 # Other function inputs
                 bblid_i = bblid_i,
                 reg_type = reg_type,
                 idvar = idvar,
                 outcome_var = outcome_var,
                 pfr_obj = pfr_obj,
                 fvars = fvars,
                 dat_orig_obs = dat_orig_obs,
                 argvals_orig_obs = argvals_orig_obs,
                 penalty = penalty,
                 lambda = lambdas,
                 penalize_type = penalize_type,
                 template_ind = template_ind,
                 scale_init = scale_init)
  }

  t1 = Sys.time()
  runtime = difftime(t0, t1, units="secs")

  # Evaluate function at the estimated value (to get the penalty values)
  loss = loglik_pen_i(bblid_i = bblid_i,
                      reg_type = reg_type,
                      warp_params_i = opt$par,
                      idvar = idvar,  outcome_var = outcome_var,
                      pfr_obj = pfr_obj, fvars = fvars,
                      dat_orig_obs = dat_orig_obs, argvals_orig_obs = argvals_orig_obs,
                      penalty = penalty,
                      lambda = lambdas,
                      penalize_type = penalize_type,
                      template_ind = template_ind,
                      scale_init = scale_init)

  # Return warp parameters and value of the loss function
  return(list(warp_params_opt_i = opt$par,
              loss_value = attributes(loss)$loss,
              loss_value_raw = attributes(loss)$loss_raw,
              runtime = runtime))
}



#' RAFT Algorithm
#'
#' The RAFT algorithm maximizes the penalized likelihood function for all individuals to estimate the optimal warp parameters.
#'
#' @param warp_params List: a list of warp parameters for all individuals.
#' @param reg_type Character string: the type of warp transformation. Options are "beta_cdf" or "poly_warp".
#' @param idvar String: the name of the variable in `dat_orig_obs` representing the IDs.
#' @param outcome_var String: the name of the variable in `dat_orig_obs` representing the outcome variable.
#' @param pfr_obj Object: the fitted `pfr` model object, if available. Default is NULL.
#' @param fvars Character vector: the names of the functional covariates.
#' @param dat_orig_obs data.frame containing the data to be warped, with one observation per row.
#' @param argvals_orig_obs String: the name of the variable in `dat_orig_obs` representing the original time grid.
#' @param penalty List: a list specifying which penalties to apply. Default is list(loglik = TRUE, dist_temp = TRUE, dist_orig = FALSE, dist_alphas = FALSE, dist_warp = TRUE).
#' @param lambdas List: a list specifying the penalty weights. Default is list(loglik = 1, dist_temp = 1, dist_orig = 1, dist_alphas = 1, dist_warp = 1).
#' @param penalize_type Character string: the type of penalization. Options are "dist" for distance-based penalization or "corr" for correlation-based penalization. Default is "dist".
#' @param template_ind Numeric scalar: the index of the template curve. Default is 5.
#' @param opt_it Numeric scalar: the number of iterations for optimization. Default is 50.
#' @param scale_init Numeric scalar: the initial scale parameter for the likelihood model. Default is NULL.
#'
#' @return A list containing the old and new warp parameters, their data frames, the distance between old and new warps, and the loss values.
#'
#' @export
max_loglik_all <- function(warp_params, reg_type = "beta_cdf", idvar = "bblid",  outcome_var = "y1",
                           pfr_obj = NULL, fvars = c("mod_fun_dmn7"),
                           dat_orig_obs, argvals_orig_obs = "density_grid",
                           penalty = list(loglik = TRUE, dist_temp = TRUE, dist_orig = FALSE, dist_alphas = FALSE, dist_warp = TRUE),
                           lambdas = list(loglik = 1, dist_temp = 1, dist_orig = 1, dist_alphas = 1, dist_warp = 1),
                           penalize_type = "dist", template_ind = NULL,
                           opt_it = 50, scale_init = NULL){

  out_warp_params = warp_params

  dat_orig_obs = as.data.frame(dat_orig_obs)
  bblids = unique(dat_orig_obs[,idvar])
  n_bblids = length(bblids)
  loss_values = vector("list", length = length(bblids))
  loss_values_raw = vector("list", length = length(bblids))

  for (i in seq_along(bblids)){
    # Print progress message for multiples of 10%
    pct = i/n_bblids
    if (pct %in% quantile(x = (1:n_bblids)/n_bblids, p = c(0.1*1:10), type = 3)) message(paste0("Progress: ", 100*round(pct, 4), "%"))


    bblid_i = bblids[i]
    # Maximize the penalized log-likelihood
    optimal_warp_i = quiet(max_loglik_i(bblid_i = bblid_i,
                                        warp_params = warp_params,
                                        reg_type = reg_type,
                                        pfr_obj = pfr_obj,
                                        idvar = idvar, outcome_var = outcome_var,
                                        dat_orig_obs = dat_orig_obs,
                                        argvals_orig_obs = argvals_orig_obs,
                                        fvars = fvars,
                                        penalty = penalty,
                                        lambdas = lambdas,
                                        penalize_type = penalize_type,
                                        template_ind = template_ind,
                                        opt_it = opt_it,
                                        scale_init = scale_init))

    # Update with the optimal value
    out_warp_params[[which(names(out_warp_params) == bblid_i)]] <- optimal_warp_i$warp_params_opt_i

    loss_values[[i]] = unlist(optimal_warp_i$loss_value)
    loss_values_raw[[i]] = unlist(optimal_warp_i$loss_value_raw)
  }

  old_warps_df = do.call(rbind,warp_params)
  new_warps_df = do.call(rbind,out_warp_params)

  # Calculate distance between old and new warps
  dist = sqrt(colMeans((old_warps_df-new_warps_df)^2))

  return(list(old_warps_df = old_warps_df,
              old_warps = warp_params,
              new_warps_df = new_warps_df,
              new_warps = out_warp_params,
              dist = dist,
              loss_values = loss_values,
              loss_values_raw = loss_values_raw))
}
