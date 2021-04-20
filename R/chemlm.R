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
chemlm.default <- function(data, outcome, chem, value = "value", adjust = NULL,
                           confound = NULL, family = "gaussian", type = "filter") {
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
    confound1 <- paste(confound, collapse = ";")
  }

  warning("Be sure categorical variables are coded as factors.  Currently cannot handle mix of categorical chemicals and numeric chemicals")

  # group by each chemical
  nest_vars <- colnames(data)[colnames(data) != chem]
  dataC <- nest(data, chemdat = one_of(nest_vars)) %>%
    # run regression, getting relevant output
    mutate(fit = map(chemdat, ~ innerchemlm(data = ., outcome = outcome, value = value,
                                            confound = confound, family = family, type = type))) %>%
    # reformat
    unnest(fit) %>%
    select(., -chemdat) %>%
    mutate(adjust = adjust1, confound = confound1)

  # if confounders included, run null model
  if(!is.null(confound)) {
    nest_vars <- colnames(data)[colnames(data) != chem]
    dataU <- nest(data, chemdat = one_of(nest_vars))  %>%
      # run regression, getting relevant output
      mutate(fit = map(chemdat, ~ innerchemlm(data = ., outcome = outcome, value = value,
                                              confound = NULL, family = family, type = type))) %>%
      # reformat
      unnest(fit) %>%
      select(., -chemdat) %>%
      mutate(adjust = adjust1, confound = "none")

    dataC <- full_join(dataC, dataU)
  }



  chem <- list()
  chem$results <- dataC
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
#' @param type Whether to return only chemical effect (default is filter)
#' @export
#' @examples
#' data(simchemdat)
#' # limit to one chemical
#' dat1 <- filter(simchemdat, chem == "chem1")
#' res <- innerchemlm(dat1, outcome = "out")
innerchemlm <- function(data, outcome, value = "value", confound = NULL, family = "gaussian", type = "filter") {


  # get linear predictor
  confound1 <- paste(c(value, confound), collapse = "+")
  # get equation
  eqn <- paste0(outcome, "~", confound1)

  data <- data[, c(value, confound, outcome)]
  data <- data %>%
    filter_all(all_vars(!is.infinite(.))) %>%
                  filter_all(all_vars(!is.na(.)))


  # run model, get 95% CI
  glm1 <- glm(eval(eqn), family = family, data = data)
  output <- summary(glm1)$coef
  output <- try(data.frame(names = rownames(output), output, suppressMessages(confint(glm1))))

  if(class(output) == "try-error") {browser()}

  colnames(output) <- c("names", "est", "SE", "z", "pvalue", "lb", "ub")
  #return all
  if(type == "filter") {
    output <- output[value, ]
  }

  if(family == "binomial") {
    output <- mutate(output, OR = exp(est), ORlb = exp(lb), ORub = exp(ub))
  }

  return(output)
}






#' Running one regression model for chemical data
#'
#' \code{innerchemlmer} Run linear mixed effects model with weights for chemical data
#'
#' This is a function to run a single regression model for chemical data
#'
#' @title innerchemlmer
#' @param data dataset
#' @param outcome character outcome variable (string format)
#' @param id ID variable
#' @param weight weight variable
#' @param value column name for value (string format)
#' @param confound character vector of confound variable (default is null)
#' @param type Whether to return only chemical effect (default is filter)
#' @export
#' @examples
#' data(simchemdatmixed)
#' # limit to one chemical
#' dat1 <- dplyr::filter(simchemdatmixed, chem == "chem1")
#' res <- innerchemlmer(dat1, outcome = "out", id = "id", weights = "weights")
innerchemlmer <- function(data, outcome, id, weights, value = "value",
                           confound = NULL, type = "filter") {


  # get linear predictor
  confound1 <- paste(paste(c(value, confound), collapse = "+"), " + (1 |", id, ")")
  # get equation
  eqn <- paste0(outcome, "~", confound1)

  data <- data[, c(value, confound, outcome, id, weights)]
  data <- data %>%
    filter_all(all_vars(!is.infinite(.))) %>%
    filter_all(all_vars(!is.na(.)))


  # run model, get 95% CI
  lmer1 <- lmer(eval(eqn), data = data, weights = weights) %>%
    tidy(conf.int = T)

  if(class(lmer1) == "try-error") {browser()}

  #return all
  if(type == "filter") {
    lmer1 <- dplyr::filter(lmer1, effect == "fixed", term == "value") %>%
      select(term, estimate : conf.high)
  }

  lmer1 <- rename(lmer1, names = term, est = estimate, SE = std.error,
                  z = statistic, lb = conf.low, ub = conf.high)
  return(lmer1)
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

  # if logistic model
  if("OR" %in% colnames(res)) {
    res <- select(res, -est, -lb, -ub) %>%
      rename(., est = OR, lb = ORlb, ub = ORub)
    ylab1 <- "OR"
    yint <- 1
  } else {
    ylab1 <- "Estimate"
    yint <- 0
  }

  g1 <- ggplot(res) + geom_pointrange(aes(x = chem, y = est, ymin = lb, ymax = ub,
                                    colour = confound), position = position_dodge(0.5)) +
    ylab(ylab1) +
    xlab("") +
    geom_hline(yintercept = yint, colour = "grey50", linetype = 2) +
    theme_bw() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "top")

  if(facetchem) {
    g1 <- g1 + facet_wrap(~chem, scales = scales, ncol = ncol)
  }

  g1
}
