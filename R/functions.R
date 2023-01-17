#' Perform a Bayes Factor \eqn{t}-test on two named groups in a dataframe
#'
#' @description
#'
#' Performs a Bayes factor $t$-test on a named test and control group in the provided
#' dataframe. A thin wrapper around the `BayesFactor::ttestBF()` function. Always set up so that effect size is test - control.
#'
#' @param df dataframe
#' @param group_col name of column with group names
#' @param data_col name of column with measurements in
#' @param control name of control group
#' @param test name of test group
#' @param h_1 the h_1 hypothesis, one of `c("test_greater_than_control", "test_smaller_than_control","test_not_equal_to_control")` default = "test_greater_than_control"
#' @param rscale scale for the prior, one of `c("medium", "wide", "ultrawide")` default = "medium"
#' @return dataframe with one row and 7 columns
#' \enumerate{
#' \item `control_group` - the name of the group considered control
#' \item. `test_group` - the group considered test
#' \item `h_0` - statement of the $H_0$ hypothesis
#' \item `h_1` - statement of the $H_1$ hypothesis
#' \item `BayesFactor` - the resulting Bayes Factor from \eqn{H_{1}/H_{0}}
#' \item `odds_h_1` - statement of the Bayes Factor in terms of odds $H_0:H_1$
#' \item `summary` - summary statement of relative evidence for the hypothesis
#' }
#' @examples
#' named_pair_ttestbf(PlantGrowth, group_col='group', data_col='weight', control='crtl', test='trt2')
#' @export
named_pair_ttestbf <- function(df, group_col=NA, data_col=NA, control=NA, test=NA, h_1 = "test_greater_than_control", rscale="medium") {

  #argcheck(df,group_col,data_col,control,test,h_1,rscale)

  if (! "tbl_df" %in% class(df)){
    df <- dplyr::mutate_if(df, is.factor, as.character) %>%
      tibble::as_tibble()
  }

  c_grp_data <- dplyr::filter(df, .data[[group_col]] == control)[[data_col]]
  t_grp_data <- dplyr::filter(df, .data[[group_col]] == test)[[data_col]]

  null_hyp <- glue::glue("{test} equal to {control}")
  alt_hyp <- glue::glue("{test} not equal to {control}")
  null_interval = NULL

  if (h_1 == "test_greater_than_control"){
    alt_hyp <- glue::glue("{test} greater than {control}")
    null_interval = c(0,Inf)
  }
  else if (h_1 == "test_smaller_than_control"){
    alt_hyp <- glue::glue("{test} less than {control}")
    nullInterval = c(-Inf,0)
  }

  bfres <- BayesFactor::ttestBF(x = t_grp_data,
                             y = c_grp_data,
                             rscale = rscale,
                             nullInterval = null_interval)

  bf <- BayesFactor::extractBF(bfres)$bf[1]
  summ <- get_summ(bf)
  data.frame(
         control_group = c(control),
         test_group = c(test),
         h_0 = c(null_hyp),
         h_1 = c(alt_hyp),
         BayesFactor = c(bf),
         odds_h_1 = paste0("1:",round(bf, 4)),
         summary = c(get_summ(bf))
  )

}



#' Automates Bayes Factor \eqn{t}-test on all pairs of groups in a dataframe
#' @description
#'
#' Performs a Bayes factor $t$-test for all pairs of groups in a named group column in the provided
#' dataframe. Group pairings are made in both directions 'grp1 - grp2' and 'grp2 - grp1'.
#' A thin wrapper around the `BayesFactor::ttestBF()` function. Always set up so that effect size is test - control.
#'
#' @param df dataframe
#' @param group_col name of column with group names
#' @param data_col name of column with measurements in
#' @param h_1 the h_1 hypothesis, one of `c("test_greater_than_control", "test_smaller_than_control","test_not_equal_to_control")` default = "test_greater_than_control"
#' @param rscale scale for the prior, one of `c("medium", "wide", "ultrawide")` default = "medium"
#' @return dataframe with 7 columns
#' \enumerate{
#' \item `control_group` - the name of the group considered control
#' \item. `test_group` - the group considered test
#' \item `h_0` - statement of the $H_0$ hypothesis
#' \item `h_1` - statement of the $H_1$ hypothesis
#' \item `BayesFactor` - the resulting Bayes Factor from \eqn{H_{1}/H_{0}}
#' \item `odds_h_1` - statement of the Bayes Factor in terms of odds $H_0:H_1$
#' \item `summary` - summary statement of relative evidence for the hypothesis
#' }
#' @examples
#' allpairs_ttestbf(PlantGrowth, group_col='group', data_col='weight')
#' @export
allpairs_ttestbf <- function(df, group_col=NA, data_col=NA, h_1 = "test_greater_than_control", rscale="medium"){

  # check <- ArgumentCheck::newArgCheck()
  # if (is.na(group_col) | is.na(data_col)) {
  #   ArgumentCheck::addError(
  #     msg = "group_col or data_col not defined.",
  #     check = check
  #   )
  # }

  #ArgumentCheck::finishArgCheck(check)

  if (! "tbl_df" %in% class(df)){
  df <- dplyr::mutate_if(df, is.factor, as.character) %>%
    tibble::as_tibble()
  }

  pairs <- expand.grid(unique(df[[group_col]]), unique(df[[group_col]]) ) %>%
            dplyr::filter(Var1 != Var2) %>%
            dplyr::rename(control = Var1, test = Var2) %>%
            dplyr::mutate_if( is.factor, as.character) %>%
            tibble::as_tibble()

  purrr::map2(
    pairs$control, pairs$test,
    function(control,test){
      named_pair_ttestbf(df, group_col = group_col, data_col = data_col, control = control, test = test, h_1 = h_1, rscale = rscale)
    }
  ) %>%
    dplyr::bind_rows()
}

#' Automates Bayes Factor proportion test on all pairs of groups in a dataframe
#' @description
#'
#' Performs a Bayes factor proportoin test for all pairs of groups in a named group column in the provided
#' dataframe. Group pairings are made in one direction only since order makes no difference.
#' A thin wrapper around the `BayesFactor::contingencyTableBF()` function. Always set up so that effect size is test - control.
#'
#' @param df dataframe
#' @param group_col name of column with group names
#' @param data_col name of column with measurements in
#' @param sample_type sampling scheme, one of `c("indepMulti", "jointMulti", "poisson")` default = "indepMulti"
#' @return dataframe with 7 columns
#' \enumerate{
#' \item `control_group` - the name of the group considered control
#' \item. `test_group` - the group considered test
#' \item `h_0` - statement of the $H_0$ hypothesis
#' \item `h_1` - statement of the $H_1$ hypothesis
#' \item `BayesFactor` - the resulting Bayes Factor from \eqn{H_{1}/H_{0}}
#' \item `odds_h_1` - statement of the Bayes Factor in terms of odds $H_0:H_1$
#' \item `summary` - summary statement of relative evidence for the hypothesis
#' }
#' @examples
#' allpairs_ttestbf(PlantGrowth, group_col='group', data_col='weight')
#' @export
allpairs_proportionbf <- function(df, group_col="NA", count_col="NA", sample_type="indepMulti"){
  if (! "tbl_df" %in% class(df)){
    df <- dplyr::mutate_if(df, is.factor, as.character) %>%
      tibble::as_tibble()
  }

  pairs <- suppressMessages(tibble::as_tibble(t(combn(unique(df[[group_col]]),2)), .name_repair="unique")) %>%
    dplyr::rename(control = ...1, test = ...2)

  purrr::map2(
    pairs$control, pairs$test,
    function(control,test){
      named_pair_proportionbf(df, group_col = group_col, count_col = count_col, control = control, test = test, sample_type = sample_type)
    }
  ) %>%
    dplyr::bind_rows()

}

#' Perform a Bayes Factor contingency table test on two named groups in a dataframe
#'
#' @description
#'
#' Performs a Bayes factor contingency table test on a named test and control group in the provided
#' dataframe. A thin wrapper around the `BayesFactor::contingencyTableBF()` function. Always set up so that effect size is test - control.
#'
#' @param df dataframe
#' @param group_col name of column with group names
#' @param data_col name of column with measurements in
#' @param control name of control group
#' @param test name of test group
#' @param sample_type sampling scheme, one of `c("indepMulti", "jointMulti", "poisson")` default = "indepMulti"
#' @return dataframe with one row and 7 columns
#' \enumerate{
#' \item `control_group` - the name of the group considered control
#' \item. `test_group` - the group considered test
#' \item `h_0` - statement of the $H_0$ hypothesis
#' \item `h_1` - statement of the $H_1$ hypothesis
#' \item `BayesFactor` - the resulting Bayes Factor from \eqn{H_{1}/H_{0}}
#' \item `odds_h_1` - statement of the Bayes Factor in terms of odds $H_0:H_1$
#' \item `summary` - summary statement of relative evidence for the hypothesis
#' }
#' @examples
#' named_pair_proportionbf(PlantGrowth, group_col='group', data_col='weight', control='crtl', test='trt2')
#' @export
named_pair_proportionbf <- function(df, group_col=NA, count_col=NA, control=NA, test=NA, sample_type="indepMulti") {

  #argcheck_proportion(df,group_col,count_col,control,test,sample_type)

  if (! "tbl_df" %in% class(df)){
    df <- dplyr::mutate_if(df, is.factor, as.character) %>%
      tibble::as_tibble()
  }
  df <- df %>% dplyr::filter(.data[[group_col]] %in% c(control, test))
  ctable <- table(df[[count_col]], df[[group_col]] )
  null_hyp <- glue::glue("{test} proportions equal to {control} proportions")
  alt_hyp <- glue::glue("{test} proportions not equal to {control} proportions")


  bfres <- BayesFactor::contingencyTableBF(ctable,
                                sampleType = sample_type,
                                fixedMargin = "cols")

  bf <- BayesFactor::extractBF(bfres)$bf
  summ <- get_summ(bf)
  data.frame(
    control_group = c(control),
    test_group = c(test),
    h_0 = c(null_hyp),
    h_1 = c(alt_hyp),
    BayesFactor = c(bf),
    odds_h_1 = paste0("1:",round(bf, 4)),
    summary = c(get_summ(bf))
  )

}

get_summ <- function(bf){
  if (bf > 100) "Extreme evidence for H_1 compared to H_0"
  else if (bf > 30 & bf < 100 ) "Very strong evidence for H_1 compared to H_0"
  else if (bf > 10 & bf < 30 ) "Strong evidence for H_1 compared to H_0"
  else if (bf > 3 & bf < 10) "Substantial evidence for H_1 compared to H_0"
  else if (bf > 1 & bf  < 3) "Anecdotal evidence for H_1 compared to H_0"
  else if (bf > 1/3 & bf < 1) "Anecdotal evidence for H_0 compared to H_1"
  else if (bf > 1/10 & bf < 1/3) "Substantial evidence for H_0 compared to H_1"
  else if (bf > 1/30 & bf < 1/10) "Strong evidence for H_0 compared to H_1"
  else if (bf > 1/100 & bf < 1/30) "Very strong evidence for H_0 compared to H_1"
  else if (bf < 1/100) "Extreme evidence for H_0 compared to H_1"
}
# argcheck <- function(df, group_col, data_col, control, test, h_1,rscale){
#   hypotheses <- c(
#     "test_greater_than_control",
#     "test_smaller_than_control",
#     "test_not_equal_to_control"
#   )
#   check <- ArgumentCheck::newArgCheck()
#   if (!h_1 %in% hypotheses){
#     ArgumentCheck::addError(
#       msg = glue::glue(
#         "Unknown H1 specified - must be one of:
#       {glue::glue_collapse(hypotheses, sep='\\n') }"),
#       argcheck = check
#     )
#   }
#   if (is.na(group_col) ){
#     ArgumentCheck::addError(
#       msg = "No value for argument group_col specified",
#       argcheck = check
#     )
#   }
#   if (is.na(data_col)){
#     ArgumentCheck::addError(
#       msg = "No value for argument data_col specified",
#       argcheck = check
#     )
#   }
#   if (is.na(control)){
#     ArgumentCheck::addError(
#       msg = "No value for argument control (control group factor level) specified",
#       argcheck = check
#     )
#   }
#   if (is.na(test)){
#     ArgumentCheck::addError(
#       msg = "No value for argument test (test group factor level) specified",
#       argcheck = check
#     )
#   }
#   if (is.null(df)){
#     ArgumentCheck::addError(
#       msg = "No data frame provided",
#       argcheck = check
#     )
#   }
#   if (!group_col %in% colnames(df)){
#     ArgumentCheck::addError(
#       msg = glue::glue("{group_col} not found in dataframe"),
#       argcheck = check
#     )
#   }
#   if (!data_col %in% colnames(df)){
#     ArgumentCheck::addError(
#       msg = glue::glue("{data_col} not found in dataframe"),
#       argcheck = check
#     )
#   }
#   if (!control %in% df[[group_col]]){
#     ArgumentCheck::addError(
#       msg = glue::glue("{control} not found in {group_col}"),
#       argcheck = check
#     )
#   }
#   if (!test %in% df[[group_col]]){
#     ArgumentCheck::addError(
#       msg = glue::glue("{test} not found in {group_col}"),
#       argcheck = check
#     )
#   }
#   if (!rscale %in% c("medium", "wide", "ultrawide")){
#     ArgumentCheck::addError(
#       msg = glue::glue("rscale value incorrect - must be one of:
#                        {glue::glue_collapse( c('medium', 'wide', 'ultrawide'), sep='\\n') }"),
#       argcheck = check
#     )
#   }
#   ArgumentCheck::finishArgCheck(check)
# }
# argcheck_proportion <- function(df, group_col, count_col, control, test, sample_type){
#   check <- ArgumentCheck::newArgCheck()
#   if (is.na(group_col) ){
#     ArgumentCheck::addError(
#       msg = "No value for argument group_col specified",
#       argcheck = check
#     )
#   }
#   if (is.na(count_col)){
#     ArgumentCheck::addError(
#       msg = "No value for argument count_col specified",
#       argcheck = check
#     )
#   }
#   if (is.na(control)){
#     ArgumentCheck::addError(
#       msg = "No value for argument control (control group factor level) specified",
#       argcheck = check
#     )
#   }
#   if (is.na(test)){
#     ArgumentCheck::addError(
#       msg = "No value for argument test (test group factor level) specified",
#       argcheck = check
#     )
#   }
#   if (is.null(df)){
#     ArgumentCheck::addError(
#       msg = "No data frame provided",
#       argcheck = check
#     )
#   }
#   if (!group_col %in% colnames(df)){
#     ArgumentCheck::addError(
#       msg = glue::glue("{group_col} not found in dataframe"),
#       argcheck = check
#     )
#   }
#   if (!count_col %in% colnames(df)){
#     ArgumentCheck::addError(
#       msg = glue::glue("{data_col} not found in dataframe"),
#       argcheck = check
#     )
#   }
#   if (!control %in% df[[group_col]]){
#     ArgumentCheck::addError(
#       msg = glue::glue("{control} not found in {group_col}"),
#       argcheck = check
#     )
#   }
#   if (!test %in% df[[group_col]]){
#     ArgumentCheck::addError(
#       msg = glue::glue("{test} not found in {group_col}"),
#       argcheck = check
#     )
#   }
#   if (!sample_type %in% c("indepMulti", "jointMulti", "poisson")){
#     ArgumentCheck::addError(
#       msg = glue::glue("sample_type value incorrect - must be one of:
#                        {glue::glue_collapse( c('indepMulti', 'jointMulti', 'poisson'), sep='\\n') }"),
#       argcheck = check
#     )
#   }
#   ArgumentCheck::finishArgCheck(check)
# }

test_func <- function(){
  
  "This function exists only to make a PR

}
