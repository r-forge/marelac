\name{ssd2rad}
\alias{ssd2rad}
\title{Estimate Global Radiation from Measured Sunshine Duration Time}
\description{
  The function converts values of sunshine duration (in hours) to global radiation
  (in \eqn{J/cm2/d}{J cm^{-2} d^{-1}})
}
\usage{
ssd2rad(S, doy, a = 0.25, b = 0.5, rho = 50.29)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Sunshine duration (hours)}
  \item{doy}{Julian day (for northern hemisphere only)}
  \item{a, b, rho}{site specific conversion parameters, must be fitted to
    measured data.}
}
\note{
  Don't forget to fit the function parameters to site specific values!
}


\value{
  Estimated value of global radiaton in \eqn{J/cm2/d}{J cm^{-2} d^{-1}}.
}

%% ToDo: find international reference or remove function from package
\references{Dyck, S. u. G. Peschke (1995): Grundlagen der Hydrologie. 3. Auflage.
  Verlag f{\"u}r Bauwesen, Berlin 1995, ISBN 3-345-00586-7}

\examples{
ssd2rad(8, 120)
}
\keyword{ misc }
