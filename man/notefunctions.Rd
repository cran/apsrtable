\name{notefunctions}
\alias{notefunctions}
\alias{se.note}
\alias{stars.note}
\title{Table notes}
\usage{se.note()
stars.note()
}
\description{Prepare notes about standard errors and statistical significance}

\details{
Table notes are part of the tabular environment and may be based on the content of the table itself. For example, the \code{stars} argument to \code{\link{apsrtable}} determines whether one or many levels of statistical significance are indicated in the output. The \code{stars.note} function creates text to place in such a note.

By default the output uses the notation \eqn{ * p <.05} and the example below shows a replacement function that states, \dQuote{significant at \code{lev} percent.}. 

Thanks to lazy evaulation, you can access variables in the call to \code{apsrtable} from functions in \code{notes}. The appropriate environment is 3 levels up the call stack (see the example).

Remember, to escape characters in Latex output, backslashes have to be doubled in R character strings.}
\value{A character string to place within the tabular environment in footnotesize beneath other output. }

\author{Michael Malecki <malecki at wustl.edu>}
\examples{
### Custom note function

signif.pct <- function() {
  env <- sys.frame(-3); 
  paste("$^*$ significant at", evalq(lev,env)*100, "percent")
}
### Continue the example from apsrtable
\dontrun{
apsrtable(lm.D90, lm.D9, glm.D9, digits=1, align="left",
          stars=1, lev=0.05, model.counter=0, order="rl",
          notes=list(se.note(), signif.pct(), 
            "Plant weight data from the lm() example" )
	 )
}
} 



