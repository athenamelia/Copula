% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_theta_unbounded_to_bounded.R
\name{transform_nac}
\alias{transform_nac}
\title{A helper method to transform vector of unbounded real numbers (as optim might use)
to a vector of parameters that are used for each copula in a NAC.}
\usage{
transform_nac(nac_Node)
}
\arguments{
\item{nac_Node}{a nested Archimedean copula}
}
\value{
a bounded theta
}
\description{
A helper method to transform vector of unbounded real numbers (as optim might use)
to a vector of parameters that are used for each copula in a NAC.
}
