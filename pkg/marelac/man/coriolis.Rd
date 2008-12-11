\name{coriolis}
\alias{coriolis}
\title{The Coriolis Force as a Function of Latitude}
\description{Estimates the coriolis factor, f (in \eqn{s^{-1}}{1/second} where
   \eqn{f=2 \cdot \omega \cdot sin(lat)}{f=2*omega*sin(lat)}, where
   \eqn{\omega =7.292e-5}{omega=7.292e-5} radians/sec}

\usage{coriolis(lat)}
\arguments{
  \item{lat }{latitude in degrees north (-90:+90).}
}
\value{
  the coriolis factor (\eqn{s^{-1}}{1/second).
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{plot(-90:90,coriolis(-90:90),xlab="latitude, dg North", 
                 ylab= "/s" , main ="coriolis factor",type="l",lwd=2)}
\references{S. Pond & G. Pickard  2nd Ed. 1986. Introductory Dynamical Oceanography,
Pergamon Press, Sydney.}
\keyword{utilities}

