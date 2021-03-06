% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{define_deme_parameters}
\alias{define_deme_parameters}
\title{Define simulation deme parameters}
\usage{
define_deme_parameters(project, H = 1000, seed_infections = 100,
  M = 1000)
}
\arguments{
\item{project}{covfefe project.}

\item{H}{vector specifying human population size in each deme.}

\item{M}{vector specifying mosquito population size (strictly the number of
adult female mosquitoes) in each deme.}

\item{Eh_init}{vector specifying the initial number of infected humans in
each deme.}
}
\description{
Define simulation deme parameters
}
