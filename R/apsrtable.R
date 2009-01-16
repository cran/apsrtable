## New APSRtable that uses much cleverer subscripting
apsrtable <- function (..., 
                       se=c("robust","vcov","both"),
                       # model.names can be shorter, others numbered;
                       # numbers start at value of model.counter
                       model.names=NULL, model.counter=1, digits=2,
                       # stars="default" prints R default
                       # this function default is one star at .05
                       stars=1,lev=.05,
                       align=c("left","center","right"),
                       order=c("lr","rl","longest"),
                       notes=list(se.note(), stars.note() ),
                       omitcoef=NULL,coef.names=NULL,
                       Sweave=FALSE,Minionfig=FALSE) {
  x <- character(0)
  signif.stars <- TRUE
  order <- match.arg(order,c("lr","rl","longest"))
  opts <- match.call(expand.dots=FALSE)
  se <- match.arg(se,c("robust","vcov","both"))
  align <- substr(align,1,1)
  se=match.arg(se,c("robust","vcov","both"))
  align <- match.arg(align,c("l","c","r"))
  adigits <- ifelse(align=="c",
                    -1,
                    digits)
  models <- list(...)
  nmodels <- length(models)
  
  if(!Sweave){
    x <- cat("\\begin{table}[!ht]\n\\caption{}\n\\label{}\n")
  }
  if(Minionfig){
    x <- cat("%Uncomment the following line and the end one to change figure versions\n%if you are using a full-featured family such as Minion Pro.\n\\figureversion{tabular}\n")
  }


  
  
  ## get the summaries for the objects
  model.summaries <- lapply(models,
                            function(x) {
                              s <- summary(x)
                              if(!is.null(x$se) && se != "vcov") {
                                est <- coef(x)
                                if(class(x$se) == "matrix") {
                                  x$se <- sqrt(diag(x$se))
                                } 
                                s$coefficients[,3] <- tval <- est / x$se
                                s$coefficients[,4] <-
                                  2 * pt(abs(tval),
                                         length(x$residuals) - x$rank,
                                         lower.tail=FALSE)
                                s$se <- x$se }
                              return(s)
                            } )

  ## Set up the model names
  ## If there's a vector of names, use that, or as many as there are
  ## and either all or the remainder.
  ## Optionally, model.number.start allows you to resetcounter
  ## TO DO: allow model "name" attribute to be used
  ##        but overridden by vector here.
  if (is.null(model.names)) {
    m.first = model.counter; m.last=m.first+(nmodels-1)
    model.names=paste("Model", m.first:m.last)
  } else if (!is.null(model.names) && (length(model.names) < nmodels) ) {
    m.first = length(model.names)+1
    model.names=c(model.names, paste( "Model", m.first:nmodels))
  }
  
## get and order the coefficient names from all models
  coefnames <- orderCoef(model.summaries, order=order)

  ## mark those to omit from the output
  incl <- rep(TRUE,length(coefnames))
  names(incl) <- coefnames
  if(!is.null(omitcoef)) {
    incl[omitcoef] <- FALSE
  }
## now figure out position of each coef in each model
model.summaries <- coefPosition(model.summaries, coefnames)
  
  ## Now that the coef name matching is done, switch to pretty names
  ## if they are supplied. 
    if(!is.null(coef.names)) {
      if(length(coef.names) != sum(incl)) {
        warning("Supplied coef.names not the same length as output. Check automatic names before supplying 'pretty' names.\n") }
      coefnames[incl] <- coef.names}
  
  
  out.table <- lapply(model.summaries, function(x){
    var.pos <- attr(x,"var.pos")
    model.out <- model.se.out <- star.out <- rep(NA,length(coefnames))
    model.out[var.pos] <- x$coefficients[,1]
    star.out[var.pos] <- apsrStars(x$coefficients,stars=stars,lev=lev,signif.stars=TRUE)
    model.out <- ifelse(!is.na(model.out),
                        paste(formatC(model.out,digits=digits,format="f"),
                              star.out),
                        "")
    
    
    
    model.se.out[var.pos] <- x$coefficients[,2]
    if( !is.null(x$se) & se %in% c("robust","both") ) {
      model.se.out[var.pos] <- x$se
    }
    
    model.se.out <- ifelse(!is.na(model.se.out),
                           paste("(",
                                 formatC(model.se.out,
                                         digits=digits,
                                         format="f"),
                                 ")",sep=""),
                           "")
    if(se=="both" && !is.null(x$se)){      
      model.se.out[var.pos] <- ifelse(model.se.out != "",
                             paste(model.se.out," [",
                                   formatC(x$coefficients[,2],
                                           digits=digits,
                                           format="f"),
                                   "]",sep=""),
                             "")
    }

    model.out <- rep(model.out[incl], each=2)
    model.se.out <- rep(model.se.out[incl], each=2)
    pos.se <- (1:length(model.out))[(1:length(model.out) %% 2==0)]
    model.out[pos.se] <- model.se.out[pos.se]
    ## Add a new model info attribute to the model's output entry
    ## To change modelInfo for a given model, change the method for it
    ## see ?modelInfo, it is reasonably well documented.
    attr(model.out,"model.info") <- modelInfo(x)
    return(model.out)
  })
  
  out.matrix <- matrix(unlist(out.table), length(coefnames[incl])*2, nmodels)

  out.matrix <- cbind(rep(coefnames[incl],each=2), out.matrix)
  out.matrix[ (row(out.matrix)[,1] %% 2 ==0) , 1] <- ""  

  out.info <- lapply(out.table, attr, "model.info")
  info.names <- orderCoef(out.info)
  out.info <- coefPosition( out.info, orderCoef(out.info) )
  out.info <- lapply(out.info, function(x) {
    var.pos <- attr(x,"var.pos")
    model.out <- rep("",length(info.names))
    model.out[var.pos] <- coef(x)
    return(model.out)
  } )

  out.info <- matrix(unlist(out.info), length(info.names), nmodels)
  out.info <- cbind(as.character(info.names), out.info)
  out.matrix <- rbind(c("%",model.names ),
                      out.matrix)
  outrows <- nrow(out.matrix)
  
  out.matrix <- rbind(out.matrix,out.info)
  out.matrix[,-1] <- format(out.matrix[,-1])
  
  out.matrix[,1] <- format(out.matrix)[,1]
  out.matrix <- apply(out.matrix, 1, paste, collapse=" & ")

  out.info <- out.matrix[ (1+outrows) : length(out.matrix) ]
  out.matrix <- out.matrix[ 1:outrows ]

  x <- cat(paste("\\begin{tabular}{",align,
                 paste("D{.}{.}{",rep(adigits,nmodels),"}",sep="",collapse="")
                 ,"}",sep="")); cat("\\hline \n &");
  x <- cat( paste("", paste("\\multicolumn{1}{",align,"}{",
                               model.names,"}",collapse=" & ")  ))
  
  x <- cat("\\\\ \\hline\n")
  x <- cat(paste(out.matrix, collapse="\\\\ \n"))
  x <- cat("\\\\\n")
  x <- cat(paste(out.info, collapse="\\\\ \n"))
  
  ## Do notes
   ## Evaluate the notes list
    ## Switch the se to either robust or regular --
  ## Robust is the default, but if only vcov are given,
  ## quietly switch the argument.
  se <- ifelse((se != "vcov" &&
                sum(unlist(lapply(models,
                                  function(x) !is.null(x$se))) >0 ) ) ,
               "robust","vcov")

  notes <- lapply(notes,evalq)
  
  cat("\\\\ \\hline\n")
  notes <- lapply(notes, function(x) {
    paste("\\multicolumn{",nmodels+1,"}{l}{\\footnotesize{", x , "}}")
             } )
  cat(paste(notes), sep="\\\\\n")
 
  
  cat("\n\\end{tabular}\n")
  if(Minionfig) {cat("\n\\figureversion{proportional}\n") }
  if(!Sweave) { cat("\\end{table}\n") }
  return(invisible(x))
}



apsrStars <- function (x, digits = max(3, getOption("digits") - 2),
                       signif.stars = getOption("show.signif.stars"), 
                       signif.legend = signif.stars,
                       dig.tst = max(1, min(5, digits - 1)), cs.ind = 1:k,
                       tst.ind = k + 1, zap.ind = integer(0), 
                       P.values = NULL,
                       has.Pvalue = nc >= 4 && substr(colnames(x)[nc],
                                      1, 3) == "Pr(",
                       eps.Pvalue = .Machine$double.eps, na.print = "NA",
                       stars="default",lev=.05,
    ...) 
{
    if (is.null(d <- dim(x)) || length(d) != 2) 
        stop("'x' must be coefficient matrix/data frame")
    nc <- d[2]
    if (is.null(P.values)) {
        scp <- getOption("show.coef.Pvalues")
        if (!is.logical(scp) || is.na(scp)) {
            warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
            scp <- TRUE
        }
        P.values <- has.Pvalue && scp
    }
    else if (P.values && !has.Pvalue) 
        stop("'P.values' is TRUE, but 'has.Pvalue' is not")
    if (has.Pvalue && !P.values) {
        d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
        nc <- nc - 1
        has.Pvalue <- FALSE
    }
    else xm <- data.matrix(x)
    k <- nc - has.Pvalue - (if (missing(tst.ind)) 
        1
    else length(tst.ind))
    if (!missing(cs.ind) && length(cs.ind) > k) 
        stop("wrong k / cs.ind")
    Cf <- array("", dim = d, dimnames = dimnames(xm))
    ok <- !(ina <- is.na(xm))
    if (length(cs.ind) > 0) {
        acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
        if (any(is.finite(acs))) {
            digmin <- 1 + floor(log10(range(acs[acs != 0], na.rm = TRUE)))
            Cf[, cs.ind] <- format(round(coef.se, max(1, digits - 
                digmin)), digits = digits)
        }
    }
    if (length(tst.ind) > 0) 
        Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst), 
            digits = digits)
    if (length(zap.ind) > 0) 
        Cf[, zap.ind] <- format(zapsmall(xm[, zap.ind], digits = digits), 
            digits = digits)
    if (any(r.ind <- !((1:nc) %in% c(cs.ind, tst.ind, zap.ind, 
        if (has.Pvalue) nc)))) 
        Cf[, r.ind] <- format(xm[, r.ind], digits = digits)
    okP <- if (has.Pvalue) 
        ok[, -nc]
    else ok
    x1 <- Cf[okP]
    dec <- getOption("OutDec")
    if (dec != ".") 
        x1 <- chartr(dec, ".", x1)
    x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
    if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
        Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1, 
            digits - 1))
    }
    if (any(ina)) 
        Cf[ina] <- na.print
    
    if (P.values) {
        if (!is.logical(signif.stars) || is.na(signif.stars)) {
            warning("option \"show.signif.stars\" is invalid: assuming TRUE")
            signif.stars <- TRUE
        }
        
        if (any(okP <- ok[, nc])) {
          pv <- as.vector(xm[, nc])
          Cf[okP, nc] <- format.pval(pv[okP], digits = dig.tst, 
                                     eps = eps.Pvalue)
          signif.stars <- signif.stars && any(pv[okP] < 0.1)
          Signif <- ""
          if (signif.stars && stars=="default") {
            Signif <- symnum(pv, corr = FALSE, na = FALSE, 
                             cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                             symbols = c("^{***}", "^{**}", "^*", "^\\dagger\ ", " "))
            Cf <- cbind(Cf, format(Signif))
          }
          else if (signif.stars && stars==1) {
            Signif <- symnum(pv, corr = FALSE, na = FALSE, 
                             cutpoints = c(0,lev,1), 
                             symbols = c("^*"," "))
          }
          return(Signif)
        }
      }

    return()
  }


setGeneric("modelInfo", function(x) standardGeneric("modelInfo") )


modelInfo.summary.lm <- function(x) {
  env <- sys.parent()
  digits <- evalq(digits, env)
  model.info <- list(
                     "$N$"=formatC(sum(x$df[1:2]),format="d"),
                     "$R^2$"=formatC(x$r.squared,format="f",digits=digits),
                     "adj. $R^2$"=formatC(x$adj.r.squared,format="f",digits=digits),
                     "Resid. sd" = formatC(x$sigma,format="f",digits=digits))
  class(model.info) <- "model.info"
  invisible(model.info) 
}

modelInfo.summary.glm <- function(x) {
  env <- sys.parent()
  digits <- evalq(digits, env)
  model.info <- list(
                       "$N$"=formatC(sum(x$df[1:2]),format="d"),
                       
                       AIC=formatC(x$aic,format="f",digits=digits),
                       BIC= formatC(
                         ( (x$aic - 2*(length(x$coef)) ) +
                             log(sum(x$df[1:2]))*
                             length(coef(x)) ),
                         format="f",digits=digits),
                       "$\\log L$"=formatC( ((x$aic - 2*(length(x$coef))) / -2),
                         format="f",digits=digits))
  class(model.info) <- "model.info"
  invisible(model.info)
}
setOldClass(c("summary.lm","summary.glm"))

setMethod("modelInfo", "summary.lm", modelInfo.summary.lm )
setMethod("modelInfo","summary.glm", modelInfo.summary.glm )

"coef.model.info" <- function(x) { x <- as.matrix(unlist(x)); invisible(x)} 


## RULES: All according to longest model,
##        then left to right
## RESULT: union of all models' coefficient names in requested order.
orderCoef <- function(model.summaries,order="lr") {
  nmodels <- length(model.summaries)
  mlength <- sapply(model.summaries, function(x) length(coef(x)) )
  longest <- which.max(mlength) # longest model
  if(order=="rl") {
    modelorder <- nmodels:1 } else {
      modelorder <- 1:nmodels }
  if(order=="longest") {
    coefnames <-  rownames(coef(model.summaries[[longest]]))
  } else {
    coefnames <- rownames(coef(model.summaries[[modelorder[1]]])) }
  
  for(i in seq_along(model.summaries)) {
    matched <- match(rownames(coef(model.summaries[[i]])), coefnames, nomatch=0)
    unmatched <- which(is.na(matched) | matched==0)
    coefnames <- c(coefnames,
                   rownames(coef(model.summaries[[i]]))[unmatched]
                   )
  }
  return(coefnames)
}
## Given a list of model summaries (or anything with a coef method),
## and a master (unioned) list of coef names,
##
## Append an attribute to each element containing its coefs' position in the
## master coefficient list
"coefPosition" <- function(model.summaries, coefnames) {
  model.summaries <- lapply(model.summaries, function(x) {
    pos <- match(rownames(coef(x)), coefnames)
    attr(x,"var.pos") <- pos
    return(x)
  })
return(model.summaries)
}

"se.note" <- function() {
  env <- sys.frame(-3)
  note <- paste(ifelse( evalq(se,env) != "vcov","Robust s","S"),
      "tandard errors in parentheses",
      ifelse(evalq(se,env)=="both",
             paste("\\\\\n\\multicolumn{",
                   evalq(nmodels,env)+1,"}{l}{",
                   'Na\\"ive standard errors in brackets',
                   collapse="",sep=""),
             "" ) ,sep="")
  return(note)
}
"stars.note" <- function() {
  env <- sys.frame(-3)
  paste(ifelse(evalq(stars,env)=="default",
               paste("$^\\dagger$ significant at $p<.10$; $^* p<.05$; $^{**} p<.01$; $^{***} p<.001$"),
               paste("$^*$ indicates significance at $p<",evalq(lev,env),"$")))
}
