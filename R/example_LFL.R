#'A realization of a \emph{Gamma-GARMA} model
#'
#'A dataset contaning the realization of a \emph{Gamma-GARMA
#'(0,1)} model, with \eqn{\phi=0.1,\theta=0.1, \Beta=2,
#'\mu_{0}=10,\sigma^2=2 }
#'
#'@format A data frame with 1000 rows and 4 variables:
#'\describe{
#'  \item{indext}{indext, time index of the series}
#'  \item{mu.t}{mu.t, the observed \eqn{\mu_{t}} value}
#'  \item{yt}{\eqn{y_{t}}, the realization of the process}
#'  \item{x}{x, an intercept variable}
#'  }
#'
#'@references For the \emph{GARMA} model:
#'
#' Benjamin, Michael A., Rigby, Robert A. and Stasinopoulos, D. Mikis. 2003. Generalized Autoregressive Moving Average Models. Journal of the American Statistical Association. Mar, 2003, Vol. 98, 461, pp. 214-223.
#'
#'

"example_LFL"
