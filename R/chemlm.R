#' Regression models for multiple chemical exposures
#'
#' \code{chemlm} Run regression model for chemical data
#'
#' This is a function to run multiple regression models for chemical data,
#' where each model is fitted to a single chemical at a time
#'
#' @title chemlm
#' @param data dataset.  data should be in "long" format with the chemical names in one column and the values in another
#' @param outcome character outcome variable
#' @param chem column name for chemical type (string format)
#' @param value column name for value (string format)
#' @param adjust function to modify chemicals (default is NULL)
#' @param confound character vector of confound variables
#' @param family regression family (default is linear model)
#' @export
#' @examples
#' data(simchemdat)
#' res <- chemlm(simchemdat, outcome = "out", chem = "chem")
#' print(res)
#' plot(res)
chemlm <- function(x, ...) UseMethod("chemlm")



#' @rdname chemlm
#' @export
chemlm.default <- function(data, outcome, chem, value = "value", adjust = NULL, confound = NULL, family = "gaussian") {
  # adjust data
  if(!is.null(adjust)) {
     stop("Sorry!  This is not yet programmed")
  } else {
    adjust1 <- "none"
  }

  # confounders
  if(is.null(confound)) {
    confound1 <- "none"
  } else {
    confound1 <- paste(confound, collapse = "-")
  }

  # group by each chemical
  data <- nest(data, -chem, .key = "chemdat")  %>%
    # run regression, getting relevant output
    mutate(fit = map(chemdat, ~ innerchemlm(data = ., outcome = outcome, value = value,
                                            confound = confound, family = family))) %>%
    # reformat
    unnest(fit) %>%
    select(., -chemdat) %>%
    mutate(adjust = adjust1, confound = confound1)


  chem <- list()
  chem$results <- data
  chem$outcome <- outcome

  class(chem) <- "chemlm"

  return(chem)
}




#' Running one regression model for chemical data
#'
#' \code{innerchemlm} Run regression model for chemical data
#'
#' This is a function to run a single regression model for chemical data
#'
#' @title innerchemlm
#' @param data dataset
#' @param outcome character outcome variable (string format)
#' @param value column name for value (string format)
#' @param confound character vector of confound variable (default is null)
#' @param family regression family (default is linear model)
#' @export
#' @examples
#' data(simchemdat)
#' # limit to one chemical
#' dat1 <- filter(simchemdat, chem == "chem1")
#' res <- innerchemlm(dat1, outcome = "out")
innerchemlm <- function(data, outcome, value = "value", confound = NULL, family = "gaussian") {


  # get linear predictor
  confound1 <- paste(c(value, confound), collapse = "+")
  # get equation
  eqn <- paste0(outcome, "~", confound1)

  # run model, get 95% CI
  glm1 <- glm(eval(eqn), family = family, data = data)
  output <- summary(glm1)$coef
  output <- data.frame(output, confint(glm1))

  colnames(output) <- c("est", "SE", "z", "pvalue", "lb", "ub")
  output <- output[value, ]


  return(output)
}


#' @rdname chemlm
#' @export
print.chemlm <- function(x) {

  cat("Outcome:\n")
  print(x$outcome)
  cat("Regression results:\n")

  print(x$results)
}


#' @rdname chemlm
#' @export
plot.chemlm <- function(x, scales = "free", ncol = 3, facetchem = F) {
  res <- x$results
  g1 <- ggplot(res) + geom_pointrange(aes(x = chem, y = est, ymin = lb, ymax = ub,
                                    colour = confound)) +
    ylab("Estimate") +
    xlab("") +
    geom_hline(yintercept = 0, colour = "grey50", linetype = 2) +
    theme_bw() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5))

  if(facetchem) {
    g1 <- g1 + facet_wrap(~chem, scales = scales, ncol = ncol)
  }

  g1
}
