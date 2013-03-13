\name{seg.profile}
\alias{seg.profile}
\title{seg profile}
\description{Run several smoothers on a profile, saving the detected breakpoint
locations to disk.}
\usage{seg.profile(profile, smooth.funs = smoothers, tosave = c("seconds", 
    "parameters", "breakpoints"), db = file.path(Sys.getenv("HOME"), 
    "seg"))}
\arguments{
  \item{profile}{Profile data.frame.}
  \item{smooth.funs}{List of smoothing functions to apply to the profile.}
  \item{tosave}{Variables to save to the db directory.}
  \item{db}{Location to save gzipped result files.}
}

\value{Nothing, the results are saved to files.}

\author{Toby Dylan Hocking}





