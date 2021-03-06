% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/math_utils.R
\name{logspace_sum}
\alias{logspace_sum}
\title{Calculate log(sum(exp(logx))) in a somewhat numerically stable way.}
\usage{
logspace_sum(logx)
}
\arguments{
\item{logx}{log-scale numeric vector of values to sum}
}
\value{
log(sum(exp(logx))), but more numerically stable
}
\description{
Calculate log(sum(exp(logx))) in a somewhat numerically stable way.
}
