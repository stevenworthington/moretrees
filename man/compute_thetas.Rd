% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_thetas.R
\name{compute_thetas}
\alias{compute_thetas}
\title{Here's a brief description.
\code{compute_thetas} Computes non-sparse effect estimates and credible
  intervals from spike & slba model output}
\usage{
compute_thetas(mod, ci_level, method, dsgn)
}
\arguments{
\item{mod}{List containing outputs from spike and slab VI algorithm}

\item{ci_level}{A number between 0 and 1 giving the desired credible interval.}

\item{method}{Description here}

\item{dsgn}{Description here
For example, ci_level = 0.95 (the default) returns a 95\% credible interval.}
}
\value{
A matrix containing estimated coefficients and credible intervals.
}
\description{
All the details go here!
}
\section{Details}{
NA
}

\seealso{
Other spike and slab functions: 
\code{\link{compute_betas}()},
\code{\link{ml_by_group}()},
\code{\link{spike_and_slab_logistic_moretrees}()},
\code{\link{spike_and_slab_logistic}()},
\code{\link{spike_and_slab_normal_moretrees}()},
\code{\link{spike_and_slab_normal}()},
\code{\link{spike_and_slab}()}
}
\concept{spike and slab functions}
