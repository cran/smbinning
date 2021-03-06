% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smbinning.R
\name{smbinning.monotonic}
\alias{smbinning.monotonic}
\title{Monotonic Binning}
\usage{
smbinning.monotonic(df, y, x, p = 0.05)
}
\arguments{
\item{df}{A data frame.}

\item{y}{Binary response variable (0,1). Integer (\code{int}) is required.
Name of \code{y} must not have a dot. Name "default" is not allowed.}

\item{x}{Continuous characteristic. At least 5 different values. Value \code{Inf} is not allowed. 
Name of \code{x} must not have a dot.}

\item{p}{Percentage of records per bin. Default 5\% (0.05).}
}
\value{
The command \code{smbinning.monotonic} generates and object containing the necessary info and utilities for binning.
The user should save the output result so it can be used 
with \code{smbinning.plot}, \code{smbinning.sql}, and \code{smbinning.gen}.
}
\description{
It gives the user the ability to impose a monotonic trend for good/bad rates per bin.
}
\examples{
# Load library and its dataset
library(smbinning) # Load package and its data

# Example 1: Monotonic Binning (Increasing Good Rate per Bin)
smbinning(df=smbsimdf2,y="fgood2",x="chr2",p=0.05)$ivtable # Run regular binning
smbinning.monotonic(df=smbsimdf2,y="fgood2",x="chr2",p=0.05)$ivtable # Run monotonic binning
 
# Example 2: Monotonic Binning (Decreasing Good Rate per Bin)
smbinning(df=smbsimdf2,y="fgood3",x="chr3",p=0.05)$ivtable # Run regular binning
smbinning.monotonic(df=smbsimdf2,y="fgood3",x="chr3",p=0.05)$ivtable # Run monotonic binning
}
