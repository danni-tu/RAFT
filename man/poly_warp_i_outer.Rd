% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignment.R
\name{poly_warp_i_outer}
\alias{poly_warp_i_outer}
\title{Wrapper function to apply Polynomial warping functions to individual observations in a dataset}
\usage{
poly_warp_i_outer(
  dat,
  idvar = "bblid",
  bblid_i,
  warp_params,
  fvars,
  inverse = FALSE,
  grid_interp,
  ...
)
}
\arguments{
\item{dat}{Data frame containing one observation in each row.}

\item{idvar}{String: the name of the variable in \code{dat} representing the IDs.}

\item{bblid_i}{String: the identifier for the individual whose data is being transformed.}

\item{warp_params}{List: a list of warp parameters corresponding to each individual.}

\item{fvars}{Character vector: the names of the functional covariates.}

\item{inverse}{Logical: if FALSE, apply the forward transformation; if TRUE, apply the inverse transformation.}

\item{grid_interp}{String: the type of interpolation method for grid transformation. Options are "loess" or "bspline".}

\item{...}{Additional arguments to be passed to \code{\link[poly_warp_i]{poly_warp_i}}.}
}
\value{
A list containing the warped data for the specified individual.
}
\description{
This function applies Polynomial warping functions to individual observations in a dataset.
}
