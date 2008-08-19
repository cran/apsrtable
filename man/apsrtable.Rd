\name{apsrtable}
\alias{apsrtable}
\title{APSR-style latex tables with multiple models}
\description{Produce well-formatted LaTeX tables of multiple models side-by-side.}
\usage{ apsrtable(..., 
                  se=c("robust","vcov","both"),
                  model.names	=NULL, model.counter=1,
                  digits=2, stars=1, lev=.05,
                  align=c("left","center","right"),
                  order=c("lr","rl","longest"),
                  omitcoef=NULL,
                  Sweave=FALSE, Minionfig=FALSE)}
\arguments{
  \item{...}{One or more fitted model objects of a supported class such as \code{lm} or \code{glm}. The model-object (a \code{list}) may also optionally contain an item named \code{se}: \code{model$se} may be a vector of standard errors, or a variance-covariance matrix, in which case the square root of the diagonal is used as the "robust" standard errors in the output. See the \code{se} argument.}
  \item{se}{Print the default standard errors ("vcov"), or those supplied by the user ("robust") in the model, or both? In the case of "both", user-supplied errors are printed in (parentheses) and the default are printed in [square brackets.] If any model in \code{...} contains an \code{se} element and "robust" is chosen (the default), output is labeled as "robust;" if no models have an \code{se} element (all use model vcov) but \code{se="robust"}, labeling is simply "Standard errors in parentheses." Default = "robust".}
  \item{model.names}{Optional vector of names to use as column headings in the table. If more models than names are supplied, unnamed models are numbered (starting at one more than the number of names).}
  \item{model.counter}{Change the number to start counting from when using automatically numbered models. Default = 1.}
  \item{digits}{Number of decimal places to report. Default = 2}
  \item{stars}{Show statistical significance "stars", either "1" or "default" where "default" is based on the output of \code{summary.lm}, except that a superscript dagger is used instead of a dot for \eqn{p < .10}. Default=1}
  \item{lev}{When \code{stars} is \code{1}, what level should be used for the test to reject statistical insignificance and bestow the glittering star? Default=.05.}
  \item{align}{How should columns be aligned in the output? Model summaries are always decimal-aligned using Dcolumn (and therefore also set in math mode), but Dcolumn also provides for decimal-point centering. Model names are set in \code{\multicolumn} spans with alignment given here, as are model terms (leftmost column of table). Default = "left".}
  \item{order}{Determines the order in which terms (rows) are included in the output when more than one model (column) is present. "lr" and "rl" take the order of terms from the first or last (leftmost or rightmost) model and appends new terms as they are encountered. "longest" uses the order of terms in the model with the most terms. Default = "lr".}
  \item{omitcoef}{An optional positive-integer, or character, vector of coefficients to exclude from the output.}
  \item{Sweave}{Toggle whether to include \code{\begin{table}...\end{table}}, \emph{empty} \code{\label{}} and \code{\caption{}}, or only the \code{\begin{tabular} ... \end{tabular}}. When called from within an \code{Sweave} document one would typically write such elements in the "documentation" (latex-part) rather than inside the code chunk. When called from an \code{Sweave} document, make sure to set the code chunk option \code{results=tex}. Default = FALSE}
  \item{Minionfig}{Include latex command to change the figure style to "tabular" and back to "proportional", specifically for the \code{MinionPro} latex package. Default = FALSE}
}

\details{Given one or several model objects of various types, \code{apsrtable()} will produce side-by-side output in well-formatted LaTeX using either automated numbering or user-supplied model names and \code{Dcolumn} decimal-aligned columns. Terms are matched across rows, with options for determining the order of terms. Nuisance terms (e.g. controls, or other quantities not of primary interest) may be omitted. Standard errors of parameter estimates are placed below estimates and in parentheses, with the option for the user to supply a replacement vector of standard errors or a replacement variance-covariance matrix, such as one estimated using the \code{sandwich} package. By default a single star denotes statistical significance at the .05 level, with the option to employ further decorations or specify another arbitrary level for the test. Finally, some model diagnostics are included along with a (somewhat) flexible means to program or include different items depending on model object class.}

\references{}
\author{Michael Malecki <malecki at wustl.edu>}
\keywords{}
\examples{ 
     ## Use the example from lm() to show both models:
     ## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
     ## Page 9: Plant Weight Data.
     ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
     trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
     group <- gl(2,10,20, labels=c("Ctl","Trt"))
     weight <- c(ctl, trt)
     lm.D9 <- lm(weight ~ group)
     lm.D90 <- lm(weight ~ group - 1) # omitting intercept
     apsrtable(lm.D90, lm.D9, digits=1, align="center", 
               stars="default", model.counter=0, order="rl")
}