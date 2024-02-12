# Functions for alignment

#' Warping Functions based on the Beta Cumulative Distribution Function (CDF)
#'
#' This function generates a warping function based on the Beta distribution using the specified
#' shape parameters, witht the option to use truncated beta functions and/or inverse functions.
#' The resulting function can be either a quantile function (inverse cumulative distribution function)
#' or a cumulative distribution function.
#'
#' @param x Numeric vector: input values on which to evaluate the function.
#' @param shape1,shape2 Numeric scalars: shape parameters of the beta distribution.
#' @param trunc_beta Logical: indicating whether to use truncated beta functions from the \code{cascsim} package.
#' @param inverse Logical: indicating whether to use the CDF function or the inverse CDF (quantile) function.
#'
#' @return Numeric vector: warped values of x based on the Beta CDF.
#'
#' @import cascsim
#'
#' @export
warp_beta <- function(x, shape1, shape2, trunc_beta, inverse){
  warp_fun = if (trunc_beta){
    if (inverse){
      function(p) cascsim::qtbeta(p, shape1 = shape1, shape2 = shape2)
    } else {
      function(q) cascsim::ptbeta(q, shape1 = shape1, shape2 = shape2)
    }
  } else {
    if (inverse){
      function(p) qbeta(p, shape1 = shape1, shape2 = shape2)
    } else {
      function(q) pbeta(q, shape1 = shape1, shape2 = shape2)
    }
  }

  return(warp_fun(x))
}

#' Warping Functions based on simple polynomials
#'
#' This function applies a polynomial warp transformation to the input data.
#'
#' @param x Numeric vector: input values on which to evaluate the function.
#' @param exp Numeric scalar: the exponent for the warp transformation.
#' @param inverse Logical: indicating whether to use the forward transformation or the inverse transformation.
#'
#' @return NNumeric vector: warped values of x based on the polynomial warping.
#'
#' @export
warp_poly <- function(x, exp, inverse){
  out = if (!inverse){
    x^(exp)
  } else {
    x^(1/exp)
  }

  return(out)
}

#' Apply Beta CDF warp to functional data
#'
#' This function applies a Beta CDF warp function to functional data.
#'
#' @param shape1,shape2 Numeric scalar: shape parameter for the beta distribution.
#' @param trunc_beta Logical: indicating whether to use truncated beta functions from the \code{cascsim} package.
#' @param dat_person data.frame containing the functional data
#' @param argval String: the name of the variable in `dat_person` representing the density grid.
#' @param fvars String: the name of the variable in `dat_person` representing the functional covariate to warp.
#' @param return_df Logical: if TRUE, return the warped data as data frames; if FALSE, return as lists.
#' @param inverse Logical: if FALSE, apply the forward transformation; if TRUE, apply the inverse transformation.
#' @param grid_interp Character string: the type of interpolation method for grid transformation. Options are "loess" or "bspline".
#'
#' @return If return_df is TRUE, a list containing two data frames: 'warped_newgrid' with the warped data on the new grid, and 'warped_origgrid' with the warped data on the original grid. If return_df is FALSE, a list containing four elements: 'warped_curves_newgrid' with the warped data on the new grid, 'mat_density_newgrid' representing the new grid, 'warped_curves_origgrid' with the warped data on the original grid, and 'mat_density_origgrid' representing the original grid.
#'
#' @export
beta_warp_i <- function(shape1 = 1, shape2 = 1, trunc_beta = FALSE,
                        dat_person, argval = "density_grid", fvars,
                        return_df = FALSE, inverse, grid_interp){
  stopifnot(grid_interp %in% c("loess", "bspline"))

  # Create density grid if none provided
  if (!(argval %in% colnames(dat_person))){
    dat_person[[argval]] = matrix(seq(from = 0, to = 1, length.out = ncol(dat_person[,fvars[1]])), byrow = TRUE, nrow = 1)
  }

  # Apply warping function to the argvals to get the newgrid
  mat_density_origgrid = dat_person[[argval]] # Original
  mat_density_newgrid = warp_beta(mat_density_origgrid, shape1 = shape1, shape2 = shape2, trunc_beta = trunc_beta, inverse = inverse) # Warped

  # Get values for the warped curve on the old and new grid
  warped_curves_newgrid = vector("list", length(fvars))
  warped_curves_origgrid = vector("list", length(fvars))

  # For each functional covariate,
  for(i in 1:length(fvars)){

    if (grid_interp == "loess"){
      loess_i = loess(as.vector(dat_person[[fvars[i]]]) ~ as.vector(mat_density_newgrid), span = 0.1)
      pred_i = predict(object = loess_i, newdata = as.vector(mat_density_origgrid))
    }

    if (grid_interp == "bspline"){
      bspline_i = try(splines::interpSpline(obj1 = mat_density_newgrid[1,], obj2 = dat_person[[fvars[i]]][1,], bSpline = TRUE))
      # If bspline doesn't work, use LOESS interpolation
      if ("try-error" %in% class(bspline_i)){
        loess_i = loess(as.vector(dat_person[[fvars[i]]]) ~ as.vector(mat_density_newgrid), span = 0.1)
        pred_i = predict(object = loess_i, newdata = as.vector(mat_density_origgrid))
      } else {
        pred_i = splines:::predict.bSpline(object = bspline_i, x = mat_density_origgrid)$y
      }
    }

    # Name of new warped curve
    new_name = paste0(fvars[i], "_warp")
    # Warped curve on the orig grid:
    # New y values, but same grid values
    warped_curves_origgrid[[i]] = t(as.matrix(pred_i))
    names(warped_curves_origgrid)[i] <- new_name

    # Warped curve on the new grid
    # Same y values, but new grid values
    warped_curves_newgrid[[i]] = dat_person[[fvars[i]]]
    names(warped_curves_newgrid)[i] <- new_name
  }

  # Optionally, return data as data.frame()
  if (return_df){
    df_newgrid = warped_curves_newgrid %>% map(as.vector) %>% Reduce(cbind,.) %>% as.data.frame %>%
      `colnames<-`(names(warped_curves_newgrid))
    df_newgrid2 = data.frame(density_warp_ng = as.vector(mat_density_newgrid)) %>%
      cbind(df_newgrid)

    df_origgrid = warped_curves_origgrid %>% map(as.vector)  %>% Reduce(cbind,.) %>% as.data.frame %>%
      `colnames<-`(names(warped_curves_origgrid))
    df_origgrid2 = data.frame(density_warp_og = as.vector(mat_density_origgrid)) %>%
      cbind(df_origgrid)

    return(list(warped_newgrid = df_newgrid2,
                warped_origgrid = df_origgrid2))
  }

  warped_curves_all = list(warped_curves_newgrid = warped_curves_newgrid,
                           mat_density_newgrid = mat_density_newgrid,
                           warped_curves_origgrid = warped_curves_origgrid,
                           mat_density_origgrid = mat_density_origgrid)

  return(warped_curves_all)
}


#' Wrapper function to apply Beta CDF warping functions to individual observations in a dataset
#'
#' This function applies Beta CDF warping functions to individual observations in a dataset.
#'
#' @param dat Data frame containing one observation in each row.
#' @param idvar String: the name of the variable in `dat` representing the IDs.
#' @param bblid_i String: the identifier for the individual whose data is being transformed.
#' @param warp_params List: a list of warp parameters corresponding to each individual.
#' @param fvars Character vector: the names of the functional covariates.
#' @param inverse Logical: if FALSE, apply the forward transformation; if TRUE, apply the inverse transformation.
#' @param grid_interp String: the type of interpolation method for grid transformation. Options are "loess" or "bspline".
#' @param ... Additional arguments to be passed to \code{\link[beta_warp_i]{beta_warp_i}}.
#'
#' @return A list containing the warped data for the specified individual.
#'
#' @export
beta_warp_i_outer <- function(dat, idvar = "bblid", bblid_i, warp_params, fvars, inverse = FALSE, grid_interp, ...){

  stopifnot(bblid_i %in% dat[,idvar])
  stopifnot(length(bblid_i) == 1)
  stopifnot(bblid_i %in% names(warp_params))

  dat_person_i = dat[which(dat[,idvar] == bblid_i),]
  warp_params_i = warp_params[[which(names(warp_params) == bblid_i)]]
  warped_curves = beta_warp_i(shape1 = warp_params_i[1], shape2 = warp_params_i[2],
                              dat_person = dat_person_i, fvars = fvars, inverse = inverse,
                              grid_interp = grid_interp, ...)

  return(warped_curves)
}



#' Apply Polynomial warp to functional data
#'
#' This function applies a Polynomial warping function to functional data.
#'
#' @param exp Numeric scalar: the exponent for the polynomial warp transformation.
#' @param trunc_beta Logical: indicating whether to use truncated beta functions from the \code{cascsim} package.
#' @param dat_person data.frame containing the functional data
#' @param argval String: the name of the variable in `dat_person` representing the density grid.
#' @param fvars String: the name of the variable in `dat_person` representing the functional covariate to warp.
#' @param return_df Logical: if TRUE, return the warped data as data frames; if FALSE, return as lists.
#' @param inverse Logical: if FALSE, apply the forward transformation; if TRUE, apply the inverse transformation.
#' @param grid_interp Character string: the type of interpolation method for grid transformation. Options are "loess" or "bspline".
#'
#' @return If return_df is TRUE, a list containing two data frames: 'warped_newgrid' with the warped data on the new grid, and 'warped_origgrid' with the warped data on the original grid. If return_df is FALSE, a list containing four elements: 'warped_curves_newgrid' with the warped data on the new grid, 'mat_density_newgrid' representing the new grid, 'warped_curves_origgrid' with the warped data on the original grid, and 'mat_density_origgrid' representing the original grid.
#'
#' @export
poly_warp_i <- function(exp = 1,
                        trunc_beta = FALSE,
                        dat_person, argval = "density_grid", fvars,
                        return_df = FALSE, inverse, grid_interp){
  stopifnot(grid_interp %in% c("loess", "bspline"))

  # Create density grid if none provided
  if (!(argval %in% colnames(dat_person))){
    dat_person[[argval]] = matrix(seq(from = 0, to = 1, length.out = ncol(dat_person[,fvars[1]])), byrow = TRUE, nrow = 1)
  }

  # Apply warping function to the argvals to get the newgrid
  mat_density_origgrid = dat_person[[argval]] # Original
  mat_density_newgrid = warp_poly(mat_density_origgrid, exp = exp, inverse = inverse) # Warped

  # Get values for the warped curve on the old and new grid
  warped_curves_newgrid = vector("list", length(fvars))
  warped_curves_origgrid = vector("list", length(fvars))

  # For each functional covariate,
  for(i in 1:length(fvars)){

    if (grid_interp == "loess"){
      loess_i = loess(as.vector(dat_person[[fvars[i]]]) ~ as.vector(mat_density_newgrid), span = 0.1)
      pred_i = predict(object = loess_i, newdata = as.vector(mat_density_origgrid))
    }

    if (grid_interp == "bspline"){
      bspline_i = try(splines::interpSpline(obj1 = mat_density_newgrid, obj2 = dat_person[[fvars[i]]][1,], bSpline = TRUE))
      # If bspline doesn't work, just use LOESS interpolation
      if ("try-error" %in% class(bspline_i)){
        loess_i = loess(as.vector(dat_person[[fvars[i]]]) ~ as.vector(mat_density_newgrid), span = 0.1)
        pred_i = predict(object = loess_i, newdata = as.vector(mat_density_origgrid))
      } else {
        pred_i = splines:::predict.bSpline(object = bspline_i, x = mat_density_origgrid)$y
      }
    }

    # Name of new warped curve
    new_name = paste0(fvars[i], "_warp")
    # Warped curve on the orig grid:
    # New y values, but same grid values
    warped_curves_origgrid[[i]] = t(as.matrix(pred_i))
    names(warped_curves_origgrid)[i] <- new_name

    # Warped curve on the new grid
    # Same y values, but new grid values
    warped_curves_newgrid[[i]] = dat_person[[fvars[i]]]
    names(warped_curves_newgrid)[i] <- new_name
  }

  # Optionally, return data as data.frame()
  if (return_df){
    df_newgrid = warped_curves_newgrid %>% map(as.vector) %>% Reduce(cbind,.) %>% as.data.frame %>%
      `colnames<-`(names(warped_curves_newgrid))
    df_newgrid2 = data.frame(density_warp_ng = as.vector(mat_density_newgrid)) %>%
      cbind(df_newgrid)

    df_origgrid = warped_curves_origgrid %>% map(as.vector)  %>% Reduce(cbind,.) %>% as.data.frame %>%
      `colnames<-`(names(warped_curves_origgrid))
    df_origgrid2 = data.frame(density_warp_og = as.vector(mat_density_origgrid)) %>%
      cbind(df_origgrid)

    return(list(warped_newgrid = df_newgrid2,
                warped_origgrid = df_origgrid2))
  }

  warped_curves_all = list(warped_curves_newgrid = warped_curves_newgrid,
                           mat_density_newgrid = mat_density_newgrid,
                           warped_curves_origgrid = warped_curves_origgrid,
                           mat_density_origgrid = mat_density_origgrid)

  return(warped_curves_all)
}

#' Wrapper function to apply Polynomial warping functions to individual observations in a dataset
#'
#' This function applies Polynomial warping functions to individual observations in a dataset.
#'
#' @param dat Data frame containing one observation in each row.
#' @param idvar String: the name of the variable in `dat` representing the IDs.
#' @param bblid_i String: the identifier for the individual whose data is being transformed.
#' @param warp_params List: a list of warp parameters corresponding to each individual.
#' @param fvars Character vector: the names of the functional covariates.
#' @param inverse Logical: if FALSE, apply the forward transformation; if TRUE, apply the inverse transformation.
#' @param grid_interp String: the type of interpolation method for grid transformation. Options are "loess" or "bspline".
#' @param ... Additional arguments to be passed to \code{\link[poly_warp_i]{poly_warp_i}}.
#'
#' @return A list containing the warped data for the specified individual.
#'
#' @export
poly_warp_i_outer <- function(dat, idvar = "bblid", bblid_i, warp_params, fvars, inverse = FALSE, grid_interp, ...){

  stopifnot(bblid_i %in% dat[,idvar])
  stopifnot(length(bblid_i) == 1)
  stopifnot(bblid_i %in% names(warp_params))

  dat_person_i = dat[which(dat[,idvar] == bblid_i),]
  warp_params_i = warp_params[[which(names(warp_params) == bblid_i)]]
  warped_curves = poly_warp_i(exp = warp_params_i[1],
                              dat_person = dat_person_i, fvars = fvars, inverse = inverse,
                              grid_interp = grid_interp, ...)

  return(warped_curves)
}



#' Wrapper function to warp data given specified warping functions and warp parameters.
#'
#' This function warps individual observations in a dataset using specified warp parameters.
#'
#' @param dat_orig_obs data.frame containing the data to be warped, with one observation per row.
#' @param reg_type String: the type of warp transformation. Options are "beta_cdf" or "poly_warp".
#' @param bblid_i String: the ID of the individual whose data is being transformed.
#' @param idvar String: the name of the variable in `dat_orig_obs` representing the IDs.
#' @param warp_params_i Numeric vector: warp parameters for the specified individual.
#' @param fvars Character vector: the names of the functional covariates.
#' @param inverse Logical: if FALSE, apply the forward transformation; if TRUE, apply the inverse transformation.
#'
#' @return The warped data for the specified individual.
#'
#' @export
warp_data_one <- function(dat_orig_obs, reg_type = "beta_cdf", bblid_i, idvar = "bblid",
                          warp_params_i = c(1,1), fvars = c("mod_fun_dmn7"), inverse = FALSE){

  dat_orig_obs = as.data.frame(dat_orig_obs)
  stopifnot(length(bblid_i) == 1)
  stopifnot(bblid_i %in% dat_orig_obs[,idvar])

  warp_params_i2 = list(warp_params_i)
  names(warp_params_i2) = bblid_i

  if (reg_type == "beta_cdf"){
    dat_warps_i = beta_warp_i_outer(bblid_i = bblid_i,
                                    dat = dat_orig_obs,
                                    idvar = idvar, warp_params = warp_params_i2,
                                    fvars = fvars, inverse = inverse,
                                    grid_interp = "bspline")

  } else if (reg_type == "poly") {
    dat_warps_i = poly_warp_i_outer(bblid_i = bblid_i,
                                    dat = dat_orig_obs,
                                    idvar = idvar, warp_params = warp_params_i2,
                                    fvars = fvars, inverse = inverse,
                                    grid_interp = "bspline")
  }

  return(dat_warps_i$warped_curves_origgrid)

}
