#' @title Firestar - Computerized Adaptive Testing (CAT) simulation program
#'
#' @author Seung W Choi, \email{schoi@@austin.utexas.edu}
#'
#' @description
#' \code{Firestar} simulates CAT with dichotomous and polytomous IRT models and generates results in various tables and plots
#'
#' @details
#' \code{Firestar} is designed for simulating CAT with dichotomous and polytomous items.
#' The item response theory models supported by the program include the dichotomous models (Birnbaum, 1968), Samejima's (1969) graded response model (GRM) and
#' Muraki's (1992) generalized partial credit model (GPCM). Both Masters' (1982) partial credit model (PCM) and
#' Andrich's (1978) rating scale model are also supported as special cases of the GPCM.
#'
#' @param filename.ipar Name of a required item parameter file (comma separated, no headers, columns in the order of id, model, a, cb1, cb2,...,cbk, blank for NA; model: 1=1PL, 2=2PL, 3=3PL, 4=PC, 5=GPC, 6=GR)
#' @param item.pool Object of item.pool class
#' @param filename.resp Name of an optional item response file (comma separated, no headers, item responses, base 1, blank for missing)
#' @param filename.content Name of an optional content specification file
#' @param ncc Number of Content Categories, effective only if content balancing is invoked by providing filename.content
#' @param filename.theta Name of an optional true or external theta file
#' @param true.theta True theta values (default: NULL)
#' @param min.score.0 TRUE if the minimum item score is 0 not 1 (default: FALSE)
#' @param simulate.theta TRUE to simulate item responses or FALSE to read in from an external file (filename.theta)
#' @param pop.dist Population distribution type for simulated theta: NORMAL, UNIFORM, or GRID
#' @param pop.par Population distribution parameters: For example, pop.par=c(M,SD) if pop.dist="NORMAL", pop.par=c(LL,UL) if pop.dist="UNIFORM", or pop.par=c(-3,-2,...3) if pop.dist="GRID"
#' @param n.simulee Toral number of simulees to generate if pop.dist in c("NORMAL","UNIFORM") or the number per theta point if pop.dist="GRID"
#' @param eap.full.length TRUE to generate EAP theta estimates based on all items or FALSE to supress
#' @param max.cat Maximum number of response categories across items
#' @param min.theta Minimum theta value
#' @param max.theta Maximum theta value
#' @param inc Theta increment value to generate a grid between min.theta and max.theta
#' @param min.NI Minimum number of items to administer (default: 4)
#' @param max.NI Maximum number of items to administer (default: 12)
#' @param max.SE Maximum SE for stopping
#' @param exposure.control TRUE to invoke exposure control or FALSE to supress (default: FALSE)
#' @param exposure.control.method Exposure control method: RD, PR, SH (defaul: "RD")
#' @param top.N Top N items from which a next item is selected randomly; effective when exposure.control.method == "Randomesque" (default: 1)
#' @param PAS A vector of the Probability of Administration given Selection, P(A|S), for each item; effective when exposure.control.method == "SH" (default: 1)
#' @param r.max Maximumum target exposure rate; effective when exposure.control.method == "SH" (default = 0.25)
#' @param stop.SE Minimum reduction in predicted SE to override continuing and stop under PSER (default: 0.01)
#' @param continue.SE Minimum reduction in predicted SE to override stopping and continue under PSER (default: 0.03)
#' @param min.SE.change Minimum reduction in SE to continue beyond satisfying min.NI (default: 0.0); not effective under PSER
#' @param extreme.response.check Check for repeated extreme responses: L for checking in the left side (low) only, R for right (high) only, E for either, or N for neither (default: N)
#' @param max.extreme.response Maximum number of responses allowed before stopping (default: 4)
#' @param selection.method Item selection method: MFI, MKL, MLWI, MPWI, MPWKL, MEI, MEPV, MEPWI, RND, KET, LOC, SEQ, TSB, PSER, MI, or AMC (default: MPWI)
#' @param info.AMC Information method for AMC: KL, MI, PWKL, or FI (default: KL)
#' @param stop.AMC Test statistic for AMC to determine whether to stop: SE, Z, LR, or ST (default: SE)
#' @param alpha.AMC Type-I error rate for AMC test statistic (default: 0.05)
#' @param BH TRUE to apply Benjamini-Hotchberg correction (default: FALSE)
#' @param interim.theta Interim theta estimator: EAP or MLE
#' @param Fisher.scoring TRUE to use Fisher's method of scoring for MLE
#' @param shrinkage.correction TRUE to correct for the bias of EAP (default: FALSE)
#' @param se.method SE estimation method: 1 = Posterior Standard Deviation or 2 = Inverse of Square Root of Information
#' @param first.item.selection Alternative first item selection method: 1 = Prior Mean, 2 = At a fixed value specified by first.at.theta, 3 = Use a specific item identified by first.item, or 4 = At external or theta values specified by filename.theta
#' @param first.at.theta Specific theta location at which the first item is optimized
#' @param first.item Specific item number to be selected as the first item
#' @param show.theta.audit.trail TRUE to generate CAT audit trail plots or FALSE to suppress
#' @param plot.usage TRUE to generate item usage plot or FALSE to suppress
#' @param plot.info TRUE to generate item intormation plots or FALSE to suppress
#' @param plot.prob TRUE to generate item response probability plots or FALSE to suppress
#' @param add.final.theta TRUE to append three additional final theta estimates (MLE, MAP, and WLE) to file.other.thetas or FALSE to supress
#' @param bank.diagnosis TRUE to generate item bank diagnostic plots or FALSE to suppress
#' @param prior.dist Type of prior distribution: 1 = Normal or 2 = Losgistic
#' @param prior.mean Prior distribution mean (default: 0.0)
#' @param prior.sd Prior distribution standard deviation (default: 1.0)
#' @param file.items.used Name of the file to contain information on items administered
#' @param file.theta.history Name of the file to contain information on history of theta estimates
#' @param file.se.history Name of the file to contain information on history of SE estimates
#' @param file.final.theta.se Name of the file to contain final theta and SE estimates
#' @param file.other.thetas Name of the file to contain other theta estimates (MLE, MAP, and WLE)
#' @param file.likelihood.dist Name of the file to contain likelihood functions
#' @param file.posterior.dist Name of the file to contain posterior distributions
#' @param file.matrix.info Name of the file to contain the item information matrix
#' @param file.full.length.theta Name of the file to contain theta estimates based on all items in the bank
#' @param file.selected.item.resp Name of the file to contain item responses for the selected items only
#' @param output.previous List object from Firestar for the previous test
#'
#' @return List of summary statistics and output results:
#' \itemize{
#'   \item{\code{call}} Call with all of the specified arguments
#'   \item{\code{nia}} Total number of items administered
#'   \item{\code{mean.nia}} Mean of the number of items administered
#'   \item{\code{cor.theta}} Correlation between true theta and theta from CAT
#'   \item{\code{rmsd.theta}} RMSE based on true theta and theta from CAT
#'   \item{\code{true.theta}} True theta
#'   \item{\code{mean.SE}} Mean standard error
#'   \item{\code{item.pool}} Item pool object
#'   \item{\code{resp}} Item response matrixc
#'   \item{\code{items.used}} Items used by examinee
#'   \item{\code{theta.history}} Theta history by examinee
#'   \item{\code{se.history}} Standard error history by examinee
#'   \item{\code{selected.item.resp}} Selected item responses by examinee
#'   \item{\code{final.theta.se}} Final theta and standard error by examinee
#'   \item{\code{likelihood.dist}} Final likelihood distribution by examinee
#'   \item{\code{posterior.dist}} Final posterior distribution by examinee
#'   \item{\code{matrix.info}} Matrix of item information
#'   \item{\code{ni.administered}} Number of items administered by examinee
#'   \item{\code{Z}} Z-test statistic if selection.method == 'AMC'
#'   \item{\code{LR}} Likelihood-ratio test statistic if selection.method == 'AMC'
#'   \item{\code{ST}} Score test statistic if selection.method == 'AMC'
#' }
#'
#' @references{
#' Andrich, D. (1978). A rating formulation for ordered response categories. Psychometrika, 43, 561-573.
#' Birnbaum, A. (1968). Some latent trait models and their use in inferring an examinee’s ability. In F. M. Lord & M. R. Novick (Eds.), Statistical theories of mental test scores (pp. 395-479). Reading, MA: Addison-Wesley.
#' Choi, S. W. (2009). Computerized Adaptive Testing Simulation Program for Polytomous IRT Models. Applied Psychological Measurement. 33, 644-645.
#' Choi, S. W., & Swartz, J. R. (2009). Comparison of CAT Item Selection Criteria for Polytomous Items. Applied Psychological Measurement. 33, 419-440.
#' Choi, S. W., Grady, M., & Dodd, B. G. (2011). A new stopping rule for computerized adaptive testing. Educational and Psychological Measurement. 71, 37-53.
#' Choi, S. W., Podrabsky, T., & McKinney, N. (2012). Firestar-D: Computerized Adaptive Testing Simulation Program for Dichotomous Item Response Theory Models. Applied Psychological Measurement, 36, 67-68.
#' Choi, S. W. (2018). Firestar: Simulating Computerized Adaptive Testing. In W. J. van der Linden (Ed.), Handbook of Item Response Theory. Chapman and Hall/CRC.
#' Finkelman, M. D., Weiss, D. J., Kim-Kang, G. (2010). Item selection and hypothesis testing for the adaptive measurement of change. Applied Psychological Measurement, 34, 238-254.
#' Masters, G. N. (1982). A Rasch model for partial credit scoring. Psychometrika, 47, 149-174.
#' Muraki, E. (1992). A generalized partial credit model: Application of an EM algorithm. Applied Psychological Measurement, 16, 159-176.
#' Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. Psychometrika Monograph Supplement, No. 17.
#' }
#'
#' @import TestDesign
#' @import stats
#' @import utils
#' @import graphics
#'
#' @export

Firestar <- function(filename.ipar = "", item.pool = NULL, filename.resp = "", filename.content = "", ncc = 1, filename.theta = "", true.theta = NULL, min.score.0 = FALSE,
                     simulate.theta = FALSE, pop.dist = "NORMAL", pop.par = c(0,1), n.simulee = 1000, eap.full.length = TRUE, max.cat = 5, min.theta = -4.0, max.theta = 4.0, inc = 0.1,
                     min.NI = 4, max.NI = 12, max.SE = 0.3, exposure.control = FALSE, exposure.control.method = "RD", top.N = 1, PAS = 1, r.max = 0.25, stop.SE = 0.01, continue.SE = 0.03, min.SE.change = 0.0, extreme.response.check = "N", max.extreme.response = 4,
                     selection.method = "MPWI", info.AMC = "KL", stop.AMC = "SE", alpha.AMC = 0.05, BH = FALSE, interim.theta = "EAP", Fisher.scoring = TRUE, shrinkage.correction = FALSE, se.method = 1,
                     first.item.selection = 1, first.at.theta = 0.0, first.item = 1, show.theta.audit.trail = FALSE, plot.usage = FALSE, plot.info = FALSE, plot.prob = FALSE, add.final.theta = FALSE, bank.diagnosis = FALSE,
                     prior.dist = 1, prior.mean = 0.0, prior.sd = 1.0, file.items.used = "", file.theta.history = "", file.se.history = "", file.final.theta.se = "", file.other.thetas = "", file.likelihood.dist = "",
                     file.posterior.dist = "", file.matrix.info = "", file.full.length.theta = "", file.selected.item.resp = "", output.previous = NULL) {

  call <- match.call()

  if (filename.ipar != "") {
    item.pool <- TestDesign::loadItemPool(filename.ipar)
  } else if (is.null(item.pool)) {
    stop("filename.ipar or item.pool must be supplied")
  }

  NCAT <- item.pool@NCAT
  ni <- item.pool@ni
  minScore <- as.numeric(!min.score.0)
  exposure.rate <- numeric(ni)

  if (exposure.control) {
    if (exposure.control.method %in% c("SH", "SYMPSON-HETTER")) {
      if (any(PAS > 1 | PAS <= 0)) {
        stop("invalid value(s) found in PAS")
      }
      if (length(PAS) == 1) {
        PAS = rep(PAS, ni)
      } else if (length(PAS) != ni) {
        stop("PAS must be a vector of length ni")
      }
      selection.rate <- numeric(ni)
    }
  }

  content.balancing <- FALSE

  if (ncc > 1 && filename.content != "" && !(toupper(selection.method) %in% c("SEQ", "TSB"))){
    target.content.dist <- as.numeric(read.csv(filename.content, header = FALSE, nrows = 1))
    if (all(target.content.dist == 0)) {
      warning("WARNING: all values in target content distribution are zero\n:content balancing not used")
    }
    content.cat <- read.csv(filename.content, header = FALSE, skip = 1)[[2]]
    if (abs(sum(target.content.dist) -1) > .1) {
      warning("WARNING: the sum of content proportions should add up to 1.0\n:content balancing not used")
    } else if (length(target.content.dist)!=ncc) {
      warning("WARNING: the number of content categories (ncc) does not match the number of target proportions in the content control file\n:content balancing not used")
    } else if (length(content.cat)!=ni) {
      warning("WARNING: the number of records in the content control file does not match the number of items in the bank\n:content balancing not used")
    } else {
      if (max.NI > sum(content.cat %in% which(target.content.dist > 0))) {
        warning("WARNING: max.NI cannot be larger than the sum of items where target.content.dist > 0")
      }
      overall.content.freq <- numeric(ncc)
      content.balancing <- TRUE
    }
  }

  .GetNextContent <- function() {
    available.content <- which(target.content.dist > 0 & as.numeric(tapply(items.available, content.cat, sum) > 0))
    idx <- which.max(target.content.dist[available.content] - current.content.dist[available.content])
    return(available.content[idx])
  }

  .UpdateContentDist <- function() {
    idx <- content.cat[item.selected]
    current.content.freq[idx] <<- current.content.freq[idx] + 1
    overall.content.freq[idx] <<- overall.content.freq[idx] + 1
    current.content.dist <<- current.content.freq / ni.given
  }

  if (filename.resp != "" && simulate.theta == FALSE) {
    resp.data <- read.csv(filename.resp, sep = ",", header = FALSE, col.names = paste("R", 1:ni, sep = ""))
    resp.matrix <- data.matrix(resp.data)
    true.theta <- NULL
  } else if (!is.null(true.theta)) {
      resp.matrix <- simResp(item.pool, true.theta)
      n.simulee <- length(true.theta)
      if (!min.score.0) resp.matrix <- resp.matrix + 1
  } else if (!is.na(n.simulee) && n.simulee > 0) {
    if (toupper(pop.dist) == "NORMAL") {
      true.theta <- rnorm(n.simulee) * pop.par[2] + pop.par[1]
      resp.matrix <- simResp(item.pool, true.theta)
    } else if (toupper(pop.dist) == "UNIFORM") {
      true.theta <- runif(n.simulee, pop.par[1], pop.par[2])
      resp.matrix <- simResp(item.pool, true.theta)
    } else if (toupper(pop.dist) == "GRID") {
      true.theta <- rep(pop.par, each = n.simulee)
      resp.matrix <- simResp(item.pool, true.theta)
    } else {
      stop("invalid option specified for pop.dist")
    }
    if (!min.score.0) resp.matrix <- resp.matrix + 1
  }

  theta <- seq(min.theta, max.theta, inc)
  nq <- length(theta)

  if (first.item.selection == 2 && first.at.theta >= min.theta && first.at.theta <= max.theta) {
    start.theta <- first.at.theta
  } else {
    start.theta <- prior.mean
  }

  if (prior.dist == 1) {
    prior <- dnorm((theta - prior.mean) / prior.sd)
  } else if (prior.dist == 2) {
    prior <- exp((theta - prior.mean) / prior.sd) / (1 + exp((theta - prior.mean) / prior.sd))^2
  } else {
    prior <- dnorm(theta)
  }

  nExaminees <- dim(resp.matrix)[1]
  items.used <- matrix(NA, nExaminees, max.NI)
  selected.item.resp <- matrix(NA, nExaminees, max.NI)
  ni.administered <- numeric(nExaminees)
  theta.CAT <-rep(NA, nExaminees)
  sem.CAT <- rep(NA, nExaminees)
  theta.history <- matrix(NA, nExaminees, max.NI)
  se.history <- matrix(NA, nExaminees, max.NI)
  posterior.matrix <- matrix(NA, nExaminees, nq)
  LH.matrix <- matrix(NA, nExaminees, nq)
  ppp <- TestDesign::calcProb(item.pool, theta)
  pp <- array(dim = c(nq, ni, max.cat))

  for (i in 1:ni) {
    pp[, i, ] <- ppp[[i]]
  }

  matrix.info <- TestDesign::calcFisher(item.pool, theta) #nq x ni

  .CalcFullLengthEAP <- function() {
    posterior <- matrix(rep(prior, nExaminees), nExaminees, nq, byrow=TRUE)
    for (i in 1:ni) {
      resp <- matrix(resp.matrix[, i], nExaminees, 1) + min.score.0
      if (!all(is.na(resp))) {
        prob <- t(pp[, i, resp])
        prob[is.na(prob)] <- 1.0
        posterior <- posterior * prob
      }
    }
    EAP <- posterior %*% theta / rowSums(posterior)
    SEM <- sqrt(rowSums(posterior * (matrix(theta, nExaminees, nq, byrow = TRUE) - matrix(EAP ,nExaminees, nq))^2) / rowSums(posterior))
    return(data.frame(theta = EAP, SE = SEM))
  }

  if (!(filename.theta == "") & eap.full.length == FALSE) {
    ext.theta <- read.csv(filename.theta,sep = ",", header = F, col.names = "theta")
  } else {
    ext.theta <- .CalcFullLengthEAP()
  }

  .CalcInfo <- function(th) {
    info <- numeric(ni)
    available <- items.available
    if (content.balancing) {
      available <- items.available & (content.cat == .GetNextContent())
    }
    if (sum(available) == 0) {
      return(info)
    }
    info <- as.vector(TestDesign::calcFisher(item.pool, th))
    info[!available] <- 0
    return(info)
  }

  .CalcLocInfo <- function(th) {
    info <- numeric(ni)
    loc <- unlist(TestDesign::calcLocation(item.pool))
    available <- items.available
    if (content.balancing) {
      available <- items.available & (content.cat == .GetNextContent())
    }
    for (i in 1:ni) {
      if (available[i] == TRUE) {
        p <- 1 / (1 + exp(-1 * (th - loc[i])))
        info[i] <- p * (1 - p)
      }
    }
    return(info)
  }

  .CalcLWInfo <- function(lk) {
    info <- numeric(ni)
    info <- apply(matrix.info * lk, 2, sum)
    info[items.available == FALSE] <- 0
    if (content.balancing) {
      info[content.cat != .GetNextContent()] <- 0
    }
    return(info)
  }

  .CalcPWInfo <- function(pos) {
    info <- numeric(ni)
    info <- apply(matrix.info * pos, 2, sum)
    info[items.available == FALSE] <- 0
    if (content.balancing) {
      info[content.cat != .GetNextContent()] <- 0
    }
    return(info)
  }

  .CalcExpectedInfo <- function(pos, current.theta) {
    info <- numeric(ni)
    available <- items.available
    if (content.balancing) {
      available <- items.available & (content.cat == .GetNextContent())
    }
    for (i in 1:ni) {
      if (available[i]) {
        ncat <- NCAT[i]
        EAP.k <- numeric(ncat)
        wt <- numeric(ncat)
        for (k in 1:ncat) {
          posterior.k <- pos * pp[, i, k]
          wt[k] <- sum(posterior.k)
          EAP.k[k] <- sum(posterior.k * theta) / wt[k]
        }
        wt <- wt / sum(wt)
        for (r in 1:ncat) {
          info.r = TestDesign::calcFisher(item.pool@parms[[i]], EAP.k[r])
          info[i] <- info[i] + wt[r] * info.r
        }
      }
    }
    return(info)
  }

  .CalcExpectedVar <- function(pos, current.theta) {
    epv <- numeric(ni)
    available <- items.available
    if (content.balancing) {
      available <- items.available & (content.cat == .GetNextContent())
    }
    for (i in 1:ni) {
      if (available[i]) {
        ncat <- NCAT[i]
        wt <- numeric(ncat)
        EAP.k <- numeric(ncat)
        Var.k <- numeric(ncat)
        for (k in 1:ncat) {
          posterior.k <- pos * pp[, i, k]
          wt[k] <- sum(posterior.k)
          EAP.k[k] <- sum(posterior.k * theta) / wt[k]
          Var.k[k] <- sum(posterior.k * (theta - EAP.k[k])^2) / wt[k]
        }
        wt <- wt / sum(wt)
        for (r in 1:ncat) {
          epv[i] <- epv[i] + wt[r] * Var.k[r]
        }
        epv[i] <- 1 / epv[i]
      }
    }
    return(epv)
  }

  .CalcMI <- function(pos) {
    MI <- numeric(ni)
    available <- items.available
    if (content.balancing) {
      available <- items.available & (content.cat == .GetNextContent())
    }
    for (i in 1:ni) {
      if (available[i]) {
        ncat <- NCAT[i]
        p <- numeric(ncat)
        posterior.k <- matrix(NA, ncat, length(pos))
        for (k in 1:ncat) {
          posterior.k[k,] <- pos * pp[, i, k]
          p[k] <- sum(posterior.k[k, ])
        }
        p <- p / sum(p)
        for (k in 1:ncat) {
          MI[i] <- MI[i] + sum(posterior.k[k, ] * log(posterior.k[k, ] / (pos * p[k])))
        }
      }
    }
    return(MI)
  }

  .CalcMI.AMC <- function(pos.H1, pos.H0) {
    MI <- numeric(ni)
    available <- items.available
    if (content.balancing) {
      available <- items.available & (content.cat == .GetNextContent())
    }
    for (i in 1:ni) {
      if (available[i]) {
        ncat <- NCAT[i]
        p <- numeric(ncat)
        posterior.k <- matrix(NA, ncat, length(pos.H1))
        for (k in 1:ncat) {
          posterior.k[k,] <- pos.H1 * pp[, i, k]
          p[k] <- sum(posterior.k[k, ])
        }
        p <- p / sum(p)
        for (k in 1:ncat) {
          MI[i] <- MI[i] + sum(posterior.k[k, ] * log(posterior.k[k, ] / (pos.H0 * p[k])))
        }
      }
    }
    return(MI)
  }

  .CalcKL <- function(current.theta, d) {
    KL <- numeric(ni)
    interval <- seq(current.theta - d, current.theta + d, length.out = 10)
    available <- items.available
    if (content.balancing) {
      available <- items.available & (content.cat == .GetNextContent())
    }
    for (i in 1:ni) {
      if (available[i]) {
        ncat <- NCAT[i]
        p <- as.vector(TestDesign::calcProb(item.pool@parms[[i]], current.theta))
        p.interval <- TestDesign::calcProb(item.pool@parms[[i]], interval)
        for (k in 1:ncat) {
          KL[i] <- KL[i] + sum(p[k] * log(p[k] / p.interval[, k]))
        }
      }
    }
    return(KL)
  }

  .CalcKL.AMC <- function(theta.H1, theta.H0) {
    KL <- numeric(ni)
    available <- items.available
    if (content.balancing) {
      available <- items.available & (content.cat == .GetNextContent())
    }
    for (i in 1:ni) {
      if (available[i]) {
        ncat <- NCAT[i]
        p.H1 <- as.vector(TestDesign::calcProb(item.pool@parms[[i]], theta.H1))
        p.H0 <- as.vector(TestDesign::calcProb(item.pool@parms[[i]], theta.H0))
        for (k in 1:ncat) {
          KL[i] <- KL[i] + p.H1[k] * log(p.H1[k] / p.H0[k])
        }
      }
    }
    return(KL)
  }

  .CalcPW.KL <- function(pos, current.theta) {
    KL <- numeric(ni)
    available <- items.available
    if (content.balancing) {
      available <- items.available & (content.cat == .GetNextContent())
    }
    for (i in 1:ni) {
      if (available[i]) {
        ncat <- NCAT[i]
        p <- as.vector(TestDesign::calcProb(item.pool@parms[[i]], current.theta))
        KL.k <- matrix(NA, ncat, length(pos))
        for (k in 1:ncat) {
          KL.k[k, ] <- p[k] * log(p[k] / pp[, i, k])
        }
        KL[i] <- sum(colSums(KL.k) * pos)
      }
    }
    return(KL)
  }

  .CalcPW.KL.AMC <- function(pos.H1, pos.H0, current.theta) {
    KL <- numeric(ni)
    available <- items.available
    if (content.balancing) {
      available <- items.available & (content.cat == .GetNextContent())
    }
    for (i in 1:ni) {
      if (available[i]) {
        ncat <- NCAT[i]
        p <- as.vector(TestDesign::calcProb(item.pool@parms[[i]], current.theta))
        KL.k <- matrix(NA, ncat, length(pos.H1))
        for (k in 1:ncat) {
          KL.k[k, ] <- p[k] * log(pos.H1 / pos.H0)
        }
        KL[i] <- sum(colSums(KL.k) * pos.H1)
      }
    }
    return(KL)
  }

  .CalcPredictedPSD <- function(pos, current.theta) {
    ppsd <- rep(NA, ni)
    available <- items.available
    if (content.balancing) {
      available <- items.available & (content.cat == .GetNextContent())
    }
    for (i in 1:ni) {
      if (available[i]) {
        ncat <- NCAT[i]
        wt <- numeric(ncat)
        EAP.k <- numeric(ncat)
        posterior.k <- matrix(NA, nq, ncat)
        for (k in 1:ncat) {
          posterior.k[,k] <- pos * pp[, i, k]
          wt[k] <- sum(pp[, i, k] * pos / sum(pos))
        }
        wt <- wt / sum(wt)
        ppsd[i] <- 0
        for (k in 1:ncat) {
          EAP.k[k] <- sum(posterior.k[, k] * theta) / sum(posterior.k[, k])
          ppsd[i] <- ppsd[i] + wt[k] * sqrt(sum(pos * pp[, i, k] * (theta - EAP.k[k])^2) / sum(pos * pp[, i, k]))
        }
      }
    }
    return(ppsd)
  }

  .CalcExpectedPWInfo <- function(pos, current.theta) {
    info <- numeric(ni)
    available <- items.available
    if (content.balancing) {
      available <- items.available & (content.cat == .GetNextContent())
    }
    for (i in 1:ni) {
      if (available[i]) {
        ncat <- NCAT[i]
        wt <- numeric(ncat)
        info.i <- matrix.info[, i]
        info.k <- numeric(ncat)
        for (k in 1:ncat) {
          posterior.k <- pos * pp[, i, k] / sum(pos * pp[, i, k])
          info.k[k] <- sum(info.i * posterior.k)
          wt[k] <- sum(pos * pp[, i, k])
        }
        wt <- wt / sum(wt)
        info[i] <- sum(info.k * wt)
      }
    }
    return(info)
  }

  .SelectMaxInfo <- function () {
    if (!exposure.control) {
      if (ni.available > 0) {
        item.selected <- info.index[1]
      }
    } else {
      if (toupper(exposure.control.method) %in% c("RANDOMESQUE", "RD")) {
        if (ni.available >= top.N) {
          item.selected <- info.index[sample(top.N, 1)]
        } else if (ni.available > 0) {
          item.selected <- info.index[sample(ni.available, 1)]
        }
      } else if (toupper(exposure.control.method) %in% c("PROGRESSIVE-RESTRICTED", "PR")) {
        rc <- runif(ni, max = max(array.info))
        rc[!items.available] <- 0
        rel <- 1 - se.history[j, ni.given]^2
        w.array.info <- rel * (1 - exposure.rate / j) * array.info + (1 - rel) * rc
        item.selected <- order(w.array.info, decreasing = TRUE)[1]
      } else if (toupper(exposure.control.method) %in% c("SYMPSON-HETTER", "SH")) {
        found = FALSE
        for (i in 1:ni.available) {
          if (PAS[info.index[i]] == 1) {
            item.selected <- info.index[i]
            selection.rate[item.selected] <<- selection.rate[item.selected] + 1
            found = TRUE
            break
          } else {
            random <- runif(1)
            selection.rate[info.index[i]] <<- selection.rate[info.index[i]] + 1
            if (random <= PAS[info.index[i]]) {
              item.selected <- info.index[i]
              found = TRUE
              break
            } else {
              items.available[info.index[i]] <<- FALSE
            }
          }
        }
        if (!found) {
          item.selected <- info.index[ni.available]
        }
      }
    }
    exposure.rate[item.selected] <<- exposure.rate[item.selected] + 1
    return (item.selected)
  }

  .CalcSE.d <- function(examinee, th, ngiven) {
    info.1 <- 0
    info.2 <- 0
    for (i in 1:output.previous$ni.administered[examinee]) {
      item <- output.previous$items.used[examinee, i]
      info.1 <- info.1 + TestDesign::calcFisher(item.pool@parms[[item]], th)
    }
    for (i in 1:ngiven) {
      item <- items.used[examinee, i]
      info.2 <- info.2 + TestDesign::calcFisher(item.pool@parms[[item]], th)
    }
    return(sqrt(1 / info.1 + 1 / info.2))
  }

  .CalcSE <- function(examinee, ngiven, th, pooled = FALSE) {
    info <- 0
    if (pooled){
      for (i in 1:output.previous$ni.administered[examinee]) {
        item <- output.previous$items.used[examinee, i]
        info <- info + TestDesign::calcFisher(item.pool@parms[[item]], th)
      }
    }
    for (i in 1:ngiven) {
      item <- items.used[examinee, i]
      info <- info + TestDesign::calcFisher(item.pool@parms[[item]], th)
    }
    SEM <- 1 / sqrt(info)
    return(SEM)
  }

  .CalcLH <- function(th, items, resps) {
    LH <- 1
    resps <- resps + min.score.0
    for (i in 1:length(items)) {
      p <- TestDesign::calcProb(item.pool@parms[[items[i]]], th)
      LH <- LH * p[resps[i]]
    }
    return(LH)
  }

  .CalcEAP <- function (examinee, ngiven, pooled = FALSE) {
    if (pooled) {
      LH <- output.previous$likelihood.dist[j, ]
    } else {
      LH <- rep(1, nq)
    }
    for (i in 1:ngiven) {
      item <- items.used[examinee, i]
      resp <- resp.matrix[examinee, item] + min.score.0
      prob <- pp[, item, resp]
      LH <- LH * prob
    }
    posterior <- prior * LH
    EAP <- sum(posterior * theta) / sum(posterior)
    if (se.method == 1) {
      SEM <- sqrt(sum(posterior * (theta - EAP)^2) / sum(posterior))
    } else if (se.method == 2) {
      SEM <- .CalcSE(examinee, ngiven, EAP, pooled = pooled)
    }
    if (shrinkage.correction) {
      EAP <- EAP * (1 + SEM^2)
      if (se.method == 1) {
        SEM <- 1 / (sqrt(1 / SEM^2 - 1 / prior.sd^2))
      } else if (se.method == 2) {
        SEM <- .CalcSE(examinee, ngiven, EAP, pooled = pooled)
      }
    }
    return(list(THETA = EAP, SEM = SEM, LH = LH, posterior = posterior))
  }

  .CalcST <- function(examinee, ngiven, th){
    ngiven.previous <- output.previous$ni.administered[examinee]
    sf1 <- 0
    sf2 <- 0
    info1 <- 0
    info2 <- 0
    for (i in 1:ngiven.previous) {
      item <- output.previous$items.used[examinee, i]
      sf1 <- sf1 + TestDesign::calcJacobian(item.pool@parms[[item]], th, output.previous$selected.item.resp[examinee, i] - !min.score.0)
      info1 <- info1 + TestDesign::calcFisher(item.pool@parms[[item]], th)
    }
    for (i in 1:ngiven) {
      item <- items.used[examinee, i]
      sf2 <- sf2 + TestDesign::calcJacobian(item.pool@parms[[item]], th, resp.matrix[examinee, item] - !min.score.0)
      info2 <- info2 + TestDesign::calcFisher(item.pool@parms[[item]], th)
    }
    ST <- sf1^2 / info1 + sf2^2 / info2
    return(ST)
  }

  .CheckScore <- function(resp, ncat) {
    not.missing <- !is.na(resp)
    sum.score <- sum(resp[not.missing])
    min.score <- sum(not.missing) * !min.score.0
    max.score <- sum(ncat[not.missing] - min.score.0)
    return(sum.score > min.score && sum.score < max.score)
  }

  .CalcMLE <- function(examinee, ngiven, maxIter = 50, crit = 0.0001, pooled = FALSE) {
    EAP.estimates <- .CalcEAP(examinee, ngiven, pooled = pooled)
    if (pooled) {
      ngiven.previous <- output.previous$ni.administered[examinee]
    } else {
      ngiven.previous <- 0
    }
    resp <- numeric(ngiven.previous + ngiven)
    ncat <- numeric(ngiven.previous + ngiven)
    if (pooled) {
      for (i in 1:ngiven.previous) {
        item <- output.previous$items.used[examinee, i]
        resp[i] <- output.previous$resp[examinee, item]
        ncat[i] <- NCAT[i]
      }
    }
    for (i in 1:ngiven) {
      item <- items.used[examinee, i]
      resp[ngiven.previous + i] <- resp.matrix[examinee, item]
      ncat[ngiven.previous + i] <- NCAT[i]
    }
    if (!.CheckScore(resp, ncat)) {
      MLE <- EAP.estimates$THETA
      SEM <- EAP.estimates$SEM
    } else {
      change <- 1000
      nIter <- 0
      post.theta <- EAP.estimates$THETA
      while (nIter <= maxIter && change > crit) {
        pre.theta <- post.theta
        deriv1 <- 0
        deriv2 <- 0
        if (pooled) {
          for (i in 1:ngiven.previous) {
            item <- output.previous$items.used[examinee, i]
            deriv1 <- deriv1 + TestDesign::calcJacobian(item.pool@parms[[item]], pre.theta, resp[i] - !min.score.0)
            if (Fisher.scoring) deriv2 <- deriv2 + TestDesign::calcFisher(item.pool@parms[[item]], pre.theta)
            else deriv2 <- deriv2 + TestDesign::calcHessian(item.pool@parms[[item]], pre.theta, resp[i] - !min.score.0)
          }
        }
        for (i in 1:ngiven) {
          item <- items.used[examinee, i]
          deriv1 <- deriv1 + TestDesign::calcJacobian(item.pool@parms[[item]], pre.theta, resp[ngiven.previous + i] - !min.score.0)
          if (Fisher.scoring) {
            deriv2 <- deriv2 + TestDesign::calcFisher(item.pool@parms[[item]], pre.theta)
          } else {
            deriv2 <- deriv2 + TestDesign::calcHessian(item.pool@parms[[item]], pre.theta, resp[ngiven.previous + i] - !min.score.0)
          }
        }
        SEM <- 1/sqrt(abs(deriv2))
        if (Fisher.scoring) {
          post.theta <- pre.theta + deriv1 / deriv2
        } else {
          post.theta <- pre.theta - deriv1 / deriv2
        }
        change <- abs(post.theta - pre.theta)
        nIter <- nIter + 1
      }
      if (post.theta < min.theta) {
        MLE <- min.theta
      } else if (post.theta > max.theta) {
        MLE <- max.theta
      } else {
        MLE <- post.theta
      }
    }
    return(list(THETA = MLE, SEM = SEM, LH = EAP.estimates$LH, posterior = EAP.estimates$posterior))
  }

  .PlotThetaAuditTrail <- function() {
    par(mfrow = c(2, 1))
    plot(1:max.NI, seq(min.theta, max.theta, length = max.NI), main = paste("CAT Audit Trail - Examinee ", j, sep = ""), xlab = "Items Administered", ylab = "Theta", type = "n", las = 1)
    points(1:ni.given, theta.history[j, 1:ni.given], type = "b", pch = 9, col = "blue")
    abline(h = theta.CAT[j], lty = 2, col = "red")
    item.string <- paste(items.used[j, 1:ni.given], collapse = ",")
    text(1, max.theta, paste("Items: ", item.string, sep = ""), cex = 0.7, adj = 0)
    text(1, min.theta + 0.3, paste("Theta: ", round(estimates$THETA, digits = 2), " SE: ", round(estimates$SEM, digits = 2)), cex = 0.8, adj = 0)
    if (content.balancing) text(1, max.theta-0.7, paste("Content Distribution:", paste(round(current.content.dist, digits = 2), collapse = ",")), cex = 0.7, adj = 0)
    for (i in 1:ni.given) {
      lines(rep(i, 2), c(theta.history[j, i]-1.96*se.history[j, i], theta.history[j, i] + 1.96*se.history[j, i]))
    }
    resp.string <- paste(resp.matrix[j, items.used[j, 1:ni.given]], collapse = ",")
    plot(theta, posterior.matrix[j, ], main = "Final Posterior Distribution", xlab = "Theta", ylab = "Posterior", type = "l", col = "blue", yaxt = "n")
    text(min.theta, max(posterior.matrix[j, ]), paste("Responses: ", resp.string, sep = ""), cex = 0.7, adj = 0)
  }

  .PlotItemUsage <- function () {
    par(mfrow = c(1, 1))
    if (!is.null(true.theta)) {
      if (toupper(pop.dist) == "GRID") {
        boxplot(rowSums(!is.na(items.used)) ~ true.theta, col = "skyblue", boxwex = 0.5, ylim = c(min.NI, max.NI), names = (format(pop.par, digits = 1)), xlim = c(1, length(pop.par)), xlab = "True Theta", ylab = "Number of Items Administered")
        for (i in min.NI:max.NI) {
          abline(h = i, lty = 3, col = "light grey")
        }
      } else {
        plot(true.theta, jitter(rowSums(!is.na(items.used)), amount = 0.01), ylim = c(min.NI, max.NI), xlab = "True Theta", ylab = "Number of Items Administered", col = 4, las = 1)
        grid()
      }
    } else  {
      plot(theta.CAT, jitter(rowSums(!is.na(items.used)), amount = 0.01), ylim = c(min.NI, max.NI), xlab = "CAT Theta", ylab = "Number of Items Administered", col = 4, las = 1)
      grid()
    }
    pct.items.used <- numeric(ni)
    tot.ni.used <- sum(!is.na(items.used))
    for (i in 1:ni) {
      pct.items.used[i] <- sum(items.used == i, na.rm = T)*100/tot.ni.used
    }
    plot(c(1, ni), c(0, max(pct.items.used)), type = "n", xlab = "Items", ylab = "Percent Used", las = 1)
    for (i in 1:ni) {
      lines(rep(i, 2), c(0, pct.items.used[i]), col = "blue")
    }
  }

  .PlotItemInfo <- function () {
    par(mfrow = c(1, 1))
    bank.info <- rowSums(matrix.info)
    bank.se <- 1/sqrt(bank.info)
    bticks <- pretty(theta,  ceiling(max.theta - min.theta))
    rticks <- pretty(1/sqrt(bank.info))
    scale.factor <- (max(bank.info)-min(bank.info)) / (max(bank.se) - min(bank.se))
    plot(theta, bank.info, type = "l", xaxt = "n", las = 2, xlab = "Theta", ylab = "Total Information", lty = 1)
    axis(4, at = rticks*scale.factor, labels = rticks, tck = 0.01, lty = 2)
    mtext("Standard Error", side = 4, line = 3)
    axis(1, at = bticks, labels = bticks, tck = 0.01, lty = 2)
    lines(theta, scale.factor * bank.se, lty = 3)
    legend(min(theta), max(bank.info), legend = c("Information", "SE"), lty = c(1, 3), bg = "white")
    par(mfrow = c(3, 4))
    max.info <- max(matrix.info)
    for (i in 1:ni){
      plot(theta, seq(0, max.info, length = length(theta)), type = "n", xlab = "Theta", ylab = "Information", main = paste("Item", i, sep = " "))
      lines(theta, matrix.info[, i], lty = 1, col = 4)
      theta.at.max <- theta[which(matrix.info[, i] == max(matrix.info[, i]))]
      points(theta.at.max, 0, pch = "|", col = 6)
      text(mean(c(min.theta, max.theta)), max.info, paste("Max at Theta =", round(theta.at.max, digits = 1), sep = ""), cex = 0.7)
    }
  }

  .PlotItemProb <- function () {
    par(mfrow = c(3, 4))
    for (i in 1:ni){
      ncat <- NCAT[i]
      plot(theta, seq(0, 1, length = length(theta)), type = "n", xlab = "Theta", ylab = "Probability", main = paste("Item", i, sep = " "))
      for (k in 1:ncat){
        lines(theta, pp[, i, k], lty = k, col = k)
      }
      legend(min(theta), 1, legend = 1:ncat, lty = 1:ncat, cex = 0.5, col = 1:ncat, bg = "white")
    }
  }

  .RunFinalThetaEstimators <- function () {
    info.wt <- sqrt(rowSums(matrix.info))
    WLH.matrix <- LH.matrix * info.wt
    find.max <- function(vec) {
      theta[which(vec == max(vec))]
    }
    theta.MAP <- apply(posterior.matrix, 1, find.max)
    theta.MLE <- apply(LH.matrix, 1, find.max)
    theta.WLE <- apply(WLH.matrix, 1, find.max)
    final.theta.estimators <- data.frame(EAP = theta.CAT, MAP = theta.MAP, MLE = theta.MLE, WLE = theta.WLE)
    return (final.theta.estimators)
  }

  .PlotMaxInfo <- function() {
    par(mfrow = c(1, 1))
    matrix.order <- t(apply(matrix.info, 1, order, decreasing = TRUE))
    sorted.info <- matrix(0, nrow(matrix.order), ncol(matrix.order))
    for (j in 1:nrow(matrix.order)) {
      sorted.info[j, ] <- matrix.info[j, matrix.order[j, ]]
    }
    cum.info <- t(apply(sorted.info, 1, cumsum))
    bticks <- pretty(theta,  ceiling(max.theta - min.theta))
    plot(theta, cum.info[, 1], xlab = "Theta", ylab = "Max Attainable Information", ylim = c(0, max(cum.info)), type = "n", cex.lab = 1.0, las = 1, xaxt = "n")
    axis(1, at = bticks, labels = bticks, tck = 1, lty = 2, col = "grey")
    box()
    abline(h = 1 / max.SE^2, col = "blue", lty = 2)
    text(min.theta, 1 / max.SE^2, paste("SE=", max.SE, sep = ""), adj = c(0, 0), cex = 0.8)
    for (i in 1:ni) {
      color <- ifelse(i%%5 == 0, "black", "grey")
      lines(theta, cum.info[, i], col = color)
      text(theta[which(cum.info[, i] == max(cum.info[, i]))], max(cum.info[, i]), i, cex = 0.5, adj = c(0, 0), col = "red")
    }
  }

  .PlotQ3 <- function(theta) {
    par(mfrow = c(1, 1))
    es <- matrix(NA, length(theta), ni)
    for (i in 1:ni){
      es[, i] <- TestDesign::calcEscore(item.pool@parms[[i]], theta)
    }
    res <- resp.matrix-es
    Q3 <- cor(res, use = "pairwise.complete.obs")
    plot(1:ni, 1:ni, type = "n", ylab = "Item", xlab = "Item", xaxt = "n", yaxt = "n")
    axis(1, at = 1:ni, cex = 0.5)
    axis(2, at = 1:ni, cex = 0.5)
    axis(3, at = 1:ni, cex = 0.5)
    axis(4, at = 1:ni, cex = 0.5)
    for (r in 1:ni) {
      for (c in 1:ni) {
        if (r>c) {
          if (abs(Q3[r, c]) >= 0.4) {
            points(r, c, pch = 16, col = "red")
            abline(h = c, lty = 2, col = "dark grey")
            abline(v = r, lty = 2, col = "dark grey")
          }
          else if (abs(Q3[r, c]) >= 0.3) points(r, c, pch = 10)
          else if (abs(Q3[r, c]) >= 0.2) points(r, c, pch = 21)
        }
      }
    }
    legend(1, ni, c("| r | >= 0.4", "| r | >= 0.3", "| r | >= 0.2"), pch = c(16, 10, 21), bg = "white", col = c("red", "black", "black"))
  }

  .CheckExtremeResponse <- function() {
    flag <- FALSE
    if (toupper(extreme.response.check) %in% c("L", "R","H", "E")) {
      if (ni.given == max.extreme.response) {
        resp.string <- paste0(selected.item.resp[j, 1:ni.given], collapse = "")
        if (toupper(extreme.response.check) == "L") {
          if (resp.string == paste0(rep(as.numeric(!min.score.0), ni.given), collapse = "")) {
            flag <- TRUE
          }
        } else if (toupper(extreme.response.check) %in% c("R",  "H")) {
          if (resp.string == paste0(NCAT[items.used[j, 1:ni.given]] - min.score.0, collapse = "")) {
            flag <- TRUE
          }
        } else {
          if (resp.string == paste0(rep(1, ni.given), collapse = "") || resp.string == paste0(NCAT[items.used[j, 1:ni.given]] - min.score.0, collapse = "")) {
            flag <- TRUE
          }
        }
      }
    }
    return (flag)
  }

  .CheckSEChange <- function() {
    flag <- FALSE
    if (min.SE.change>0 && ni.given >= min.NI && ni.given >= 2) {
      if ((se.history[j, ni.given-1] - se.history[j, ni.given]) < min.SE.change) {
        flag <- TRUE
      }
    }
    return (flag)
  }

  .CheckSE <- function() {
    return(estimates$SEM <= max.SE && ni.given >= min.NI)
  }

  .StopCAT <- function(){
    switch(stop.AMC,
           SE = estimates$SEM <= max.SE && ni.given >= min.NI,
           Z = Z[j, ni.given] >= qnorm(1 - FDR / 2, mean = 0, sd = 1) && ni.given >= min.NI,
           LR = LR[j, ni.given] >= qchisq(1 - FDR / 2, df = 1) && ni.given >= min.NI,
           ST = ST[j, ni.given] >= qchisq(1 - FDR / 2, df = 1) && ni.given >= min.NI
           )
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  pb = txtProgressBar(0, nExaminees, char = "|", style = 3)

  if (toupper(selection.method) == "MFI") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      if (first.item.selection == 4) {
        theta.current <- ext.theta$theta[j]
      } else {
        theta.current <- start.theta
      }
      while (critMet == FALSE && ni.given<max.to.administer) {
        array.info <- .CalcInfo(theta.current)
        ni.available <- sum(array.info > 0)
        info.index <- order(array.info, decreasing = TRUE)
        item.selected <- .SelectMaxInfo()
        if (ni.given == 0) {
          if (first.item.selection == 3 && first.item >= 1 && first.item <= ni) {
            if (items.available[first.item] == TRUE) {
              item.selected <- first.item
            }
          } else item.selected <- .SelectMaxInfo()
        }
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) .UpdateContentDist()
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        if (toupper(interim.theta) == "EAP") {
          estimates <- .CalcEAP(j, ni.given)
        } else if (toupper(interim.theta) == "MLE") {
          estimates <- .CalcMLE(j, ni.given)
        }
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "MKL") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      if (first.item.selection == 4) theta.current <- ext.theta$theta[j]
      else theta.current <- start.theta
      while (critMet == FALSE && ni.given<max.to.administer) {
        array.info <- .CalcKL(theta.current, 0.75)
        ni.available <- sum(array.info>0)
        info.index <- order(array.info, decreasing = TRUE)
        item.selected <- .SelectMaxInfo()
        if (ni.given == 0) {
          if (first.item.selection == 3 && first.item >= 1 && first.item <= ni) {
            if (items.available[first.item] == TRUE) {
              item.selected <- first.item
            }
          } else item.selected <- .SelectMaxInfo()
        }
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) .UpdateContentDist()
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        if (toupper(interim.theta) == "EAP") {
          estimates <- .CalcEAP(j, ni.given)
        } else if (toupper(interim.theta) == "MLE") {
          estimates <- .CalcMLE(j, ni.given)
        }
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "MLWI") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      if (first.item.selection == 4) theta.current <- ext.theta$theta[j]
      else theta.current <- start.theta
      likelihood <- rep(1, length(theta))
      while (critMet == FALSE && ni.given<max.to.administer) {
        array.info <- .CalcLWInfo(likelihood)
        ni.available <- sum(array.info>0)
        info.index <- order(array.info, decreasing = TRUE)
        item.selected <- .SelectMaxInfo()
        if (ni.given == 0) {
          if (first.item.selection == 3 && first.item >= 1 && first.item <= ni) {
            if (items.available[first.item] == TRUE) {
              item.selected <- first.item
            }
          } else if (first.item.selection == 2 || first.item.selection == 4) {
            array.info <- .CalcInfo(theta.current)
            info.index <- order(array.info, decreasing = TRUE)
            item.selected <- .SelectMaxInfo()
          }
        }
        resp <- resp.matrix[j, item.selected]
        prob <- pp[, item.selected, resp]
        likelihood <- likelihood*prob
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) .UpdateContentDist()
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        if (toupper(interim.theta) == "EAP") {
          estimates <- .CalcEAP(j, ni.given)
        } else if (toupper(interim.theta) == "MLS") {
          estimates <- .CalcMLE(j, ni.given)
        }
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "MPWI") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      if (first.item.selection == 4) {
        theta.current <- ext.theta$theta[j]
      } else {
        theta.current <- start.theta
      }
      posterior <- prior
      while (critMet == FALSE && ni.given<max.to.administer) {
        array.info <- .CalcPWInfo(posterior)
        ni.available <- sum(array.info>0)
        info.index <- order(array.info, decreasing = TRUE)
        item.selected <- .SelectMaxInfo()
        if (ni.given == 0) {
          if (first.item.selection == 3 && first.item >= 1 && first.item <= ni) {
            if (items.available[first.item] == TRUE) {
              item.selected <- first.item
            }
          } else if (first.item.selection == 2 || first.item.selection == 4) {
            array.info <- .CalcInfo(theta.current)
            info.index <- order(array.info, decreasing = TRUE)
            item.selected <- .SelectMaxInfo()
          }
        }
        resp <- resp.matrix[j, item.selected]
        prob <- pp[, item.selected, resp]
        posterior <- posterior*prob
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) .UpdateContentDist()
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        estimates <- .CalcEAP(j, ni.given)
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "MPWKL") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      if (first.item.selection == 4) theta.current <- ext.theta$theta[j]
      else theta.current <- start.theta
      posterior <- prior
      while (critMet == FALSE && ni.given<max.to.administer) {
        array.info <- .CalcPW.KL(posterior, theta.current)
        ni.available <- sum(array.info>0)
        info.index <- order(array.info, decreasing = TRUE)
        item.selected <- .SelectMaxInfo()
        if (ni.given == 0) {
          if (first.item.selection == 3 && first.item >= 1 && first.item <= ni) {
            if (items.available[first.item] == TRUE) {
              item.selected <- first.item
            }
          } else if (first.item.selection == 2 || first.item.selection == 4) {
            array.info <- .CalcInfo(theta.current)
            info.index <- order(array.info, decreasing = TRUE)
            item.selected <- .SelectMaxInfo()
          }
        }
        resp <- resp.matrix[j, item.selected]
        prob <- pp[, item.selected, resp]
        posterior <- posterior*prob
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) .UpdateContentDist()
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        estimates <- .CalcEAP(j, ni.given)
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "MEI") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      if (first.item.selection == 4) theta.current <- ext.theta$theta[j]
      else theta.current <- start.theta
      posterior <- prior
      while (critMet == FALSE && ni.given<max.to.administer) {
        array.info <- .CalcExpectedInfo(posterior, theta.current)
        ni.available <- sum(array.info>0)
        info.index <- order(array.info, decreasing = TRUE)
        item.selected <- .SelectMaxInfo()
        if (ni.given == 0) {
          if (first.item.selection == 3 && first.item >= 1 && first.item <= ni) {
            if (items.available[first.item] == TRUE) {
              item.selected <- first.item
            }
          } else if (first.item.selection == 2 || first.item.selection == 4) {
            array.info <- .CalcInfo(theta.current)
            info.index <- order(array.info, decreasing = TRUE)
            item.selected <- .SelectMaxInfo()
          }
        }
        resp <- resp.matrix[j, item.selected]
        prob <- pp[, item.selected, resp]
        posterior <- posterior*prob
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) .UpdateContentDist()
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp
        estimates <- .CalcEAP(j, ni.given)
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "MEPV") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      if (first.item.selection == 4) theta.current <- ext.theta$theta[j]
      else theta.current <- start.theta
      posterior <- prior
      while (critMet == FALSE && ni.given<max.to.administer) {
        array.info <- .CalcExpectedVar(posterior, theta.current)
        ni.available <- sum(array.info>0)
        info.index <- order(array.info, decreasing = TRUE)
        item.selected <- .SelectMaxInfo()
        if (ni.given == 0) {
          if (first.item.selection == 3 && first.item >= 1 && first.item <= ni) {
            if (items.available[first.item] == TRUE) {
              item.selected <- first.item
            }
          } else if (first.item.selection == 2 || first.item.selection == 4) {
            array.info <- .CalcInfo(theta.current)
            info.index <- order(array.info, decreasing = TRUE)
            item.selected <- .SelectMaxInfo()
          }
        }
        resp <- resp.matrix[j, item.selected]
        prob <- pp[, item.selected, resp]
        posterior <- posterior*prob
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) .UpdateContentDist()
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp
        estimates <- .CalcEAP(j, ni.given)
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "MEPWI") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      if (first.item.selection == 4) {
        theta.current <- ext.theta$theta[j]
      } else {
        theta.current <- start.theta
      }
      posterior <- prior
      while (critMet == FALSE && ni.given<max.to.administer) {
        array.info <- .CalcExpectedPWInfo(posterior, theta.current)
        ni.available <- sum(array.info>0)
        info.index <- order(array.info, decreasing = TRUE)
        item.selected <- .SelectMaxInfo()
        if (ni.given == 0) {
          if (first.item.selection == 3 && first.item >= 1 && first.item <= ni) {
            if (items.available[first.item] == TRUE) {
              item.selected <- first.item
            }
          } else if (first.item.selection == 2 || first.item.selection == 4) {
            array.info <- .CalcInfo(theta.current)
            info.index <- order(array.info, decreasing = TRUE)
            item.selected <- .SelectMaxInfo()
          }
        }
        resp <- resp.matrix[j, item.selected]
        prob <- pp[, item.selected, resp]
        posterior <- posterior*prob
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) .UpdateContentDist()
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp
        estimates <- .CalcEAP(j, ni.given)
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "RND") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) {
          items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
        }
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      while (critMet == FALSE && ni.given<max.to.administer) {
        random <- runif(ni)
        random[!items.available] <- 0
        if (content.balancing) {
          random[content.cat != .GetNextContent()] <- 0
        }
        item.selected <- order(random, decreasing = TRUE)[1]
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) {
          .UpdateContentDist()
        }
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        estimates <- .CalcEAP(j, ni.given)
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "KET") {
    info.table <- .CalcInfo(ext.theta[[1]])
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) {
          items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
        }
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      while (critMet == FALSE && ni.given<max.to.administer) {
        array.info <- info.table[j, ]
        array.info[!items.available] <- 0
        if (content.balancing) {
          array.info[content.cat != .GetNextContent()] <- 0
        }
        item.selected <- order(array.info, decreasing = TRUE)[1]
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) {
          .UpdateContentDist()
        }
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        estimates <- .CalcEAP(j, ni.given)
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "LOC") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) {
          items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
        }
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      if (first.item.selection == 4) theta.current <- ext.theta$theta[j]
      else theta.current <- start.theta
      while (critMet == FALSE  && ni.given<max.to.administer) {
        array.info <- .CalcLocInfo(theta.current)
        ni.available <- sum(array.info>0)
        info.index <- order(array.info, decreasing = TRUE)
        item.selected <- .SelectMaxInfo()
        if (ni.given == 0) {
          if (first.item.selection == 3 && first.item >= 1 && first.item <= ni) {
            if (items.available[first.item] == TRUE) {
              item.selected <- first.item
            }
          }
        }
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) {
          .UpdateContentDist()
        }
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        estimates <- .CalcEAP(j, ni.given)
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "SEQ") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      item.order <- 1:ni
      if (sum(items.available) < ni) {
        item.order <- item.order[-which(items.available == FALSE)]
      }
      while (critMet == FALSE && ni.given<length(item.order)) {
        item.selected <- item.order[ni.given + 1]
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        estimates <- .CalcEAP(j, ni.given)
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "TSB") {
    items.available <- rep(TRUE, ni)
    array.info <- .CalcPWInfo(prior)
    info.index <- order(array.info, decreasing = TRUE)
    locator.item <- info.index[1]
    ncat <- NCAT[locator.item]
    sequence.matrix <- matrix(NA, ncat, ni)
    posterior.k <- matrix(rep(prior, ncat), ncat, length(prior), byrow = TRUE)
    for (k in 1:ncat) {
      posterior.k[k, ] <- prior*pp[, locator.item, k]
      array.info <- .CalcPWInfo(posterior.k[k, ])
      info.index <- order(array.info, decreasing = TRUE)
      info.index <- info.index[-which(info.index == locator.item)]
      sequence.matrix[k, ] <- c(locator.item, info.index)
    }
    par(mfrow = c(1, 1))
    plot(1:15, seq(1, (ncat + 4), length = 15), xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "")
    text(2, round(median(1:ncat)) + 2, paste("Locator: ", locator.item, sep = ""))
    for (k in 1:ncat) {
      text(5, ncat + 3 - k, paste("(", k, ")", sep = ""))
      text(6, ncat + 3 - k, paste(sequence.matrix[k, 2:max.NI], collapse = "-"), adj = c(0), cex = 0.7, col = "blue")
    }
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      locator.response <- resp.matrix[j, locator.item]
      if (is.na(locator.response) == TRUE) {
        locator.response <- (round(median(1:ncat)))
      }
      item.order <- sequence.matrix[locator.response, ]
      if (sum(items.available)<ni) {
        item.order <- setdiff(item.order, which(items.available == FALSE))
      }
      while (critMet == FALSE && ni.given < max.to.administer) {
        item.selected <- item.order[ni.given + 1]
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        estimates <- .CalcEAP(j, ni.given)
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "PSER") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      old.SE <- 1
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) {
          items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
        }
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      if (first.item.selection == 4) {
        theta.current <- ext.theta$theta[j]
      } else {
        theta.current <- start.theta
      }
      posterior <- prior
      while (ni.given<max.to.administer) {
        array.se <- .CalcPredictedPSD(posterior, theta.current)
        array.info <- 1/array.se^2
        array.info[is.na(array.info)] <- 0
        ni.available <- sum(array.info>0)
        info.index <- order(array.info, decreasing = TRUE)
        item.selected <- .SelectMaxInfo()
        if (ni.given == 0) {
          if (first.item.selection == 3 && first.item >= 1 && first.item <= ni) {
            if (items.available[first.item] == TRUE) {
              item.selected <- first.item
            }
          } else if (first.item.selection == 2 || first.item.selection == 4) {
            array.info <- .CalcInfo(theta.current)
            info.index <- order(array.info, decreasing = TRUE)
            item.selected <- .SelectMaxInfo()
          } else item.selected <- .SelectMaxInfo()
        }
        se.change <- old.SE - array.se[item.selected]
        if (ni.given >= min.NI) {
          if ((old.SE <= max.SE) && (se.change < continue.SE)) {
            break
          } else if ((old.SE > max.SE) && (se.change < stop.SE)) {
            break
          }
        }
        if (.CheckExtremeResponse()) {
          break
        }
        resp <- resp.matrix[j, item.selected]
        prob <- pp[, item.selected, resp]
        posterior <- posterior*prob
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) {
          .UpdateContentDist()
        }
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        estimates <- .CalcEAP(j, ni.given)
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        old.SE <- estimates$SEM
        theta.current <- estimates$THETA
        theta.CAT[j] <- estimates$THETA
        sem.CAT[j] <- estimates$SEM
        LH.matrix[j, ] <- estimates$LH
        posterior.matrix[j, ] <- estimates$posterior
        ni.administered[j] <- ni.given
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "MI") {
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) {
          items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
        }
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      if (first.item.selection == 4) {
        theta.current <- ext.theta$theta[j]
      } else {
        theta.current <- start.theta
      }
      posterior <- prior
      while (critMet == FALSE && ni.given < max.to.administer) {
        array.info <- .CalcMI(posterior)
        ni.available <- sum(array.info>0)
        info.index <- order(array.info, decreasing = TRUE)
        item.selected <- .SelectMaxInfo()
        if (ni.given == 0) {
          if (first.item.selection == 3 && first.item >= 1 && first.item <= ni) {
            if (items.available[first.item] == TRUE) {
              item.selected <- first.item
            }
          } else if (first.item.selection == 2 || first.item.selection == 4) {
            array.info <- .CalcInfo(theta.current)
            info.index <- order(array.info, decreasing = TRUE)
            item.selected <- .SelectMaxInfo()
          }
        }
        resp <- resp.matrix[j, item.selected]
        prob <- pp[, item.selected, resp]
        posterior <- posterior*prob
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) {
          .UpdateContentDist()
        }
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        estimates <- .CalcEAP(j, ni.given)
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        if (ni.given >= max.to.administer || .CheckSE() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  if (toupper(selection.method) == "AMC") {
    if (is.null(output.previous)) {
      stop("for selection.method = \'AMC\' output.previous is required")
    }
    posterior.matrix.previous <- output.previous$posterior.dist
    LH.matrix.previous <- output.previous$likelihood.dist
    if (nrow(posterior.matrix.previous)!=nExaminees) {
      stop("the number of simulees in output.previous must be equal to nExaminees")
    }
    Z <- matrix(NA, nExaminees, max.NI)
    LR <- matrix(NA, nExaminees, max.NI)
    ST <- matrix(NA, nExaminees, max.NI)
    if (BH) {
      n.tests <- max.NI - min.NI + 1
      if (n.tests > 1) {
        alpha.AMC.BH <- (1:n.tests)*rep(alpha.AMC, max.NI - min.NI + 1)/n.tests
      }
    }
    for (j in 1:nExaminees) {
      critMet <- FALSE
      items.available <- rep(TRUE, ni)
      items.available[is.na(resp.matrix[j, 1:ni])] <- FALSE
      if (content.balancing) {
        current.content.dist <- numeric(ncc)
        current.content.freq <- numeric(ncc)
        if (any(target.content.dist == 0)) items.available[content.cat %in% which(target.content.dist == 0)] <- FALSE
      }
      max.to.administer <- ifelse(sum(items.available) <= max.NI, sum(items.available), max.NI)
      ni.given <- 0
      FDR <- alpha.AMC
      if (first.item.selection == 4) {
        theta.current <- output.previous$final.theta.se[j, 1]
      } else {
        theta.current <- start.theta
      }
      while (critMet == FALSE && ni.given<max.to.administer) {
        if (ni.given == 0) {
          array.info <- .CalcInfo(theta.current)
        } else {
          if (toupper(info.AMC) == "KL") {
            array.info <- .CalcKL.AMC(theta.current, theta.pooled)
          } else if (toupper(info.AMC) == "PWKL") {
            array.info <- .CalcPW.KL.AMC(estimates$posterior, estimates.pooled$posterior, theta.current)
          } else if (toupper(info.AMC) == "MI") {
            array.info <- .CalcMI.AMC(estimates$posterior, estimates.pooled$posterior)
          } else if (toupper(info.AMC) == "FI") {
            array.info <- .CalcInfo(theta.current)
          }
        }
        ni.available <- sum(array.info > 0)
        info.index <- order(array.info, decreasing = TRUE)
        item.selected <- .SelectMaxInfo()
        if (ni.given == 0) {
          if (first.item.selection == 3 && first.item >= 1 && first.item <= ni) {
            if (items.available[first.item] == TRUE) {
              item.selected <- first.item
            }
          }
        }
        ni.given <- ni.given + 1
        items.used[j, ni.given] <- item.selected
        if (content.balancing) .UpdateContentDist()
        items.available[item.selected] <- FALSE
        selected.item.resp[j, ni.given] <- resp.matrix[j, item.selected]
        if (toupper(interim.theta) == "EAP") {
          estimates <- .CalcEAP(j, ni.given)
          estimates.pooled <- .CalcEAP(j, ni.given, pooled = TRUE)
        } else if (toupper(interim.theta) == "MLE") {
          estimates <- .CalcMLE(j, ni.given)
          estimates.pooled <- .CalcMLE(j, ni.given, pooled = TRUE)
        }
        theta.history[j, ni.given] <- estimates$THETA
        se.history[j, ni.given] <- estimates$SEM
        theta.current <- estimates$THETA
        theta.pooled <- estimates.pooled$THETA
        Z[j, ni.given] <- abs(theta.current-output.previous$final.theta.se[j, 1])/.CalcSE.d(j, theta.pooled, ni.given)
        LR.0 <- .CalcLH(theta.pooled, c(output.previous$items.used[j, 1:output.previous$ni.administered[j]], items.used[j, 1:ni.given]), c(output.previous$selected.item.resp[j, 1:output.previous$ni.administered[j]], selected.item.resp[j, 1:ni.given]))
        LR.1 <- .CalcLH(theta.current, items.used[j, 1:ni.given], selected.item.resp[j, 1:ni.given]) * .CalcLH(output.previous$final.theta.se[j, 1], output.previous$items.used[j, 1:output.previous$ni.administered[j]], output.previous$selected.item.resp[j, 1:output.previous$ni.administered[j]])
        LR[j, ni.given] <- -2*(log(LR.0/LR.1))
        ST[j, ni.given] <- .CalcST(j, ni.given, theta.pooled)
        if (BH && ni.given >= min.NI) {
          FDR <- alpha.AMC.BH[ni.given - min.NI + 1]
        }
        if (ni.given >= max.to.administer || .StopCAT() || .CheckExtremeResponse() || .CheckSEChange()) {
          critMet <- TRUE
          theta.CAT[j] <- estimates$THETA
          sem.CAT[j] <- estimates$SEM
          LH.matrix[j, ] <- estimates$LH
          posterior.matrix[j, ] <- estimates$posterior
          ni.administered[j] <- ni.given
        }
      }
      if (show.theta.audit.trail) {
        .PlotThetaAuditTrail()
      }
      setTxtProgressBar(pb, j)
    }
  }

  cat("\n")
  par(mfrow = c(1, 1))

  if (!is.null(true.theta)) {
    cor.theta <- round(cor(true.theta, theta.CAT), 3)
    rmsd.theta <- round(sqrt(mean((true.theta - theta.CAT)^2)), 3)
    if (toupper(pop.dist) == "GRID") {
      boxplot(theta.CAT~true.theta, col = "yellow", boxwex = 0.5, ylim = c(min.theta, max.theta), names = (format(pop.par, digits = 1)), xlim = c(1, length(pop.par)), xlab = "True Theta", ylab = "CAT Theta", main = "CAT vs. True Theta")
      text(pop.par[1], max.theta, adj = 0, paste("r = ", cor.theta, sep = ""))
      segments(1, min.theta, length(pop.par), max.theta, col = "red", lwd = 2, lty = 2)
      for (i in min.theta:max.theta) {
        abline(h = i, lty = 3, col = "light grey")
      }
    } else {
      plot(min.theta:max.theta, min.theta:max.theta, xlab = "True Theta", ylab = "CAT Theta", main = "CAT vs. True Theta", type = "n", las = 1)
      points(true.theta, theta.CAT, col = "blue")
      text(min.theta, max.theta, adj = 0, paste("r = ", cor.theta, sep = ""))
      abline(0, 1)
      grid()
    }
  } else if (eap.full.length) {
    cor.theta <- round(cor(ext.theta$theta, theta.CAT), 3)
    rmsd.theta <- round(sqrt(mean((ext.theta$theta - theta.CAT)^2)), 3)
    plot(min.theta:max.theta, min.theta:max.theta, xlab = "Full-Bank Theta", ylab = "CAT Theta", main = "CAT vs. Full-Bank Theta Estimates", type = "n", las = 1)
    points(ext.theta$theta, theta.CAT, col = "blue")
    text(min.theta, max.theta, adj = 0, paste("r = ", cor.theta, sep = ""))
    abline(0, 1)
    grid()
  } else {
    cor.theta <- round(cor(ext.theta$theta, theta.CAT), 3)
    rmsd.theta <- round(sqrt(mean((ext.theta$theta - theta.CAT)^2)), 3)
    plot(min.theta:max.theta, min.theta:max.theta, xlab = "External Theta", ylab = "CAT Theta", main = "CAT vs. External Theta Estimates", type = "n", las = 1)
    points(ext.theta$theta, theta.CAT, col = "blue")
    text(min.theta, max.theta, adj = 0, paste("r = ", cor.theta, sep = ""))
    abline(0, 1)
    grid()
  }

  if (plot.usage) {
    .PlotItemUsage()
  }

  if (content.balancing) {
    overall.content.dist <- overall.content.freq/sum(overall.content.freq)
    content.dist <- rbind(target.content.dist, overall.content.dist)
    par.mar <- par()$mar
    par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 5))
    barplot(content.dist, ylab = "Proportion", xlab = "Content Category", las = 1, names.arg = paste(1:ncc), col = c("black", "grey"), beside = TRUE)
    legend(ncc*3 + 0.5, max(target.content.dist, current.content.dist)/2, c("Target", "Current"), fill = c("black", "grey"))
    par(xpd = FALSE, mar = par.mar)
  }

  if (bank.diagnosis) {
    .PlotMaxInfo()
    .PlotQ3(ext.theta$theta)
  }

  if (plot.info) {
    .PlotItemInfo()
  }

  if (plot.prob) {
    .PlotItemProb()
  }

  if (add.final.theta || !(file.other.thetas == "")) {
    final.thetas <- .RunFinalThetaEstimators()
    if (add.final.theta) {
      pairs(final.thetas,
            panel = function(x, y) {
              points(x, y, col = 4)
              abline(0, 1, lwd = 2)
            },
            diag.panel = function(x) {
              par(new = TRUE)
              hist(x, main = "", axes = FALSE, nclass = 12, probability = TRUE)
              lines(density(x))
              points(x, rep(0, length(x)), pch = "|")
            }
      )
    }
    if (!(file.other.thetas == "")) {
      write.table(final.thetas, file = file.other.thetas, sep = ",", na = " ", row.names = FALSE, col.names = TRUE)
    }
  }

  if (!(file.items.used == "")) {
    colnames(items.used) <- paste("Item", seq(1:max.NI), sep = "")
    write.table(items.used, file = file.items.used, sep = ",", na = " ", row.names = FALSE, col.names = TRUE)
  }

  if (!(file.theta.history == "")) {
    colnames(theta.history) <- paste("Item", seq(1:max.NI), sep = "")
    write.table(theta.history, file = file.theta.history, sep = ",", na = " ", row.names = FALSE, col.names = TRUE)
  }

  if (!(file.se.history == "")) {
    colnames(se.history) <- paste("Item", seq(1:max.NI), sep = "")
    write.table(se.history, file = file.se.history, sep = ",", na = " ", row.names = FALSE, col.names = TRUE)
  }

  final.theta.se <- as.data.frame(cbind(theta.CAT, sem.CAT))
  if (!(file.final.theta.se == "")) {
    colnames(final.theta.se) <- c("Theta", "SEM")
    write.table(final.theta.se, file = file.final.theta.se, sep = ",", na = " ", row.names = FALSE, col.names = TRUE)
  }

  if (!(file.likelihood.dist == "")) {
    colnames(LH.matrix) <- paste("Theta=", theta, sep = "")
    write.table(LH.matrix, file = file.likelihood.dist, sep = ",", na = " ", row.names = FALSE, col.names = TRUE)
  }

  if (!(file.posterior.dist == "")) {
    colnames(posterior.matrix) <- paste("Theta=", theta, sep = "")
    write.table(posterior.matrix, file = file.posterior.dist, sep = ",", na = " ", row.names = FALSE, col.names = TRUE)
  }

  if (!(file.matrix.info == "")) {
    colnames(matrix.info) <- paste("Item", 1:ni, sep = "")
    write.table(matrix.info, file = file.matrix.info, sep = ",", na = " ", row.names = FALSE, col.names = TRUE)
  }

  if (!(file.full.length.theta == "")) {
    write.table(ext.theta, file = file.full.length.theta, sep = ",", na = " ", row.names = FALSE, col.names = TRUE)
  }

  if (!(file.selected.item.resp == "")) {
    colnames(selected.item.resp) <- paste("Item", seq(1:max.NI), sep = "")
    write.table(selected.item.resp, file = file.selected.item.resp, sep = ",", na = " ", row.names = FALSE, col.names = TRUE)
  }

  nia <- rowSums(!is.na(items.used))
  mean.nia <- mean(nia)
  mean.SE <- mean(sem.CAT)
  exposure.rate <- exposure.rate / j
  out <- list(call = call, nia, mean.nia = mean.nia, cor.theta = cor.theta, rmsd.theta = rmsd.theta, exposure.rate = exposure.rate, true.theta = true.theta, mean.SE = mean.SE, item.pool = item.pool, resp = resp.matrix, items.used = items.used, theta.history = theta.history, se.history = se.history, selected.item.resp = selected.item.resp, final.theta.se = final.theta.se, likelihood.dist = LH.matrix, posterior.dist = posterior.matrix, matrix.info = matrix.info, ni.administered = ni.administered)

  if (toupper(selection.method) == "AMC") {
    out[['Z']] <- Z
    out[['LR']] <- LR
    out[['ST']] <- ST
  }

  if (exposure.control.method %in% c("SH", "SYMPSON-HETTER")) {
    selection.rate <- selection.rate / j
    K <- rep(1, ni)
    K[selection.rate > r.max] <- r.max / selection.rate[selection.rate > r.max]
    out[['K']] <- K
  }

  return(out)
}

