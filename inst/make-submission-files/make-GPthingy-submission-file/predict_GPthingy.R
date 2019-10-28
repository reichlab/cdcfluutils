# R code for SIR fit

library(readr)
library(dplyr)
library(ggplot2)
library(hforecast)
library(here)
library(grid)
library(tidyr)

# Set working directory to package root if you started off anywhere within the package directory structure
setwd(here())

# Command line arguments specify season and season week of "last observation"
# We predict forward from that point.
# args <- c("2019/2020", "22", "11", "8")
args <- commandArgs(trailingOnly=TRUE)
predict_season <- args[1]
predict_last_obs_season_week <- as.integer(args[2])
tree_depth <- as.integer(args[3])
chain_num <- as.integer(args[4])

set.seed(as.integer(paste0(
  substr(predict_season, 1, 4),
  predict_last_obs_season_week,
  chain_num)))

# Cluster
cmdstan_root <- "~/cmdstan"
output_path <- "/project/uma_nicholas_reich/flu-prediction/GPthingy/"

# Evan's machine
cmdstan_root <- "~/Documents/research/cmdstan/"
output_path <- "~/Documents/research/epi/flu/cdcfluutils/inst/make-submission-files/make-GPthingy-submission-file/"

# Some rstan internals here just because I was unable to install the rstan
# package on the cluster.

stan_kw1 <- c("for", "in", "while", "repeat", "until", "if", "then", "else", "true", "false")
stan_kw2 <- c("int", "real", "vector", "simplex", "ordered", "positive_ordered", "row_vector", "matrix", "corr_matrix", "cov_matrix", "lower", "upper")
stan_kw3 <- c("model", "data", "parameters", "quantities", "transformed", "generated")
cpp_kw <- c("alignas", "alignof", "and", "and_eq", "asm", "auto", "bitand",
"bitor", "bool", "break", "case", "catch", "char", "char16_t",
"char32_t", "class", "compl", "const", "constexpr", "const_cast",
"continue", "decltype", "default", "delete", "do", "double",
"dynamic_cast", "else", "enum", "explicit", "export", "extern",
"false", "float", "for", "friend", "goto", "if", "inline", "int",
"long", "mutable", "namespace", "new", "noexcept", "not", "not_eq",
"nullptr", "operator", "or", "or_eq", "private", "protected",
"public", "register", "reinterpret_cast", "return", "short",
"signed", "sizeof", "static", "static_assert", "static_cast",
"struct", "switch", "template", "this", "thread_local", "throw",
"true", "try", "typedef", "typeid", "typename", "union", "unsigned",
"using", "virtual", "void", "volatile", "wchar_t", "while", "xor",
"xor_eq")

data_list2array <- function (x) {
    len <- length(x)
    if (len == 0L)
        return(NULL)
    dimx1 <- dim(x[[1]])
    if (any(sapply(x, function(xi) !is.numeric(xi))))
        stop("all elements of the list should be numeric")
    if (is.null(dimx1))
        dimx1 <- length(x[[1]])
    lendimx1 <- length(dimx1)
    if (len > 1) {
        d <- sapply(x[-1], function(xi) {
            dimxi <- dim(xi)
            if (is.null(dimxi))
                dimxi <- length(xi)
            identical(dimxi, dimx1)
        })
        if (!all(d))
            stop("the dimensions for all elements (array) of the list are not same")
    }
    x <- do.call(c, x)
    dim(x) <- c(dimx1, len)
    aperm(x, c(lendimx1 + 1L, seq_len(lendimx1)))
}

is_legal_stan_vname <- function (name) {
    if (grepl("\\.", name))
        return(FALSE)
    if (grepl("^\\d", name))
        return(FALSE)
    if (grepl("__$", name))
        return(FALSE)
    if (name %in% stan_kw1)
        return(FALSE)
    if (name %in% stan_kw2)
        return(FALSE)
    if (name %in% stan_kw3)
        return(FALSE)
    !name %in% cpp_kw
}

real_is_integer <- function (x) {
    if (length(x) < 1L)
        return(TRUE)
    if (any(is.infinite(x)) || any(is.nan(x)))
        return(FALSE)
    all(floor(x) == x)
}

stan_rdump <- function (list, file = "", append = FALSE, envir = parent.frame(),
    width = options("width")$width, quiet = FALSE) {
    if (is.character(file)) {
        ex <- sapply(list, exists, envir = envir)
        if (!all(ex)) {
            notfound_list <- list[!ex]
            if (!quiet)
                warning(paste("objects not found: ", paste(notfound_list,
                  collapse = ", "), sep = ""))
        }
        list <- list[ex]
        if (!any(ex))
            return(invisible(character()))
        if (nzchar(file)) {
            file <- file(file, ifelse(append, "a", "w"))
            on.exit(close(file), add = TRUE)
        }
        else {
            file <- stdout()
        }
    }
    for (x in list) {
        if (!is_legal_stan_vname(x) & !quiet)
            warning(paste("variable name ", x, " is not allowed in Stan",
                sep = ""))
    }
    l2 <- NULL
    addnlpat <- paste0("(.{1,", width, "})(\\s|$)")
    for (v in list) {
        vv <- get(v, envir)
        if (is.data.frame(vv)) {
            vv <- data.matrix(vv)
        }
        else if (is.list(vv)) {
            vv <- data_list2array(vv)
        }
        else if (is.logical(vv)) {
            mode(vv) <- "integer"
        }
        else if (is.factor(vv)) {
            vv <- as.integer(vv)
        }
        if (!is.numeric(vv)) {
            if (!quiet)
                warning(paste0("variable ", v, " is not supported for dumping."))
            next
        }
        if (!is.integer(vv) && max(abs(vv)) < .Machine$integer.max &&
            real_is_integer(vv))
            storage.mode(vv) <- "integer"
        if (is.vector(vv)) {
            if (length(vv) == 0) {
                cat(v, " <- integer(0)\n", file = file, sep = "")
                next
            }
            if (length(vv) == 1) {
                cat(v, " <- ", as.character(vv), "\n", file = file,
                  sep = "")
                next
            }
            str <- paste0(v, " <- \nc(", paste(vv, collapse = ", "),
                ")")
            str <- gsub(addnlpat, "\\1\n", str)
            cat(str, file = file)
            l2 <- c(l2, v)
            next
        }
        if (is.matrix(vv) || is.array(vv)) {
            l2 <- c(l2, v)
            vvdim <- dim(vv)
            cat(v, " <- \n", file = file, sep = "")
            if (length(vv) == 0) {
                str <- paste0("structure(integer(0), ")
            }
            else {
                str <- paste0("structure(c(", paste(as.vector(vv),
                  collapse = ", "), "),")
            }
            str <- gsub(addnlpat, "\\1\n", str)
            cat(str, ".Dim = c(", paste(vvdim, collapse = ", "),
                "))\n", file = file, sep = "")
            next
        }
    }
    invisible(l2)
}


# Set up stan model for estimation

# Seasons to be used as part of training or prediction data set
seasons <- paste0(1997:(as.integer(substr(predict_season, 1, 4))), "/", 1998:(as.integer(substr(predict_season, 6, 9))))

stan_data <- cdcfluutils::preprocess_reg_data_for_GPthingy_stan(
  us_flu = readRDS("inst/make-submission-files/make-GPthingy-submission-file/us_flu.rds"),
  seasons = seasons,
  predict_season = predict_season,
  predict_last_obs_season_week = predict_last_obs_season_week,
  get_all_region_seasons = TRUE,
  pred_all_t = FALSE
)
stan_data$U <- stan_data$U_max
stan_data$ili_t_concat[stan_data$ili_t_concat == 0] <- 0.001
stan_data$ili_t_concat[stan_data$ili_t_concat == 1] <- 0.999
stan_data$prop_a_t_concat[stan_data$prop_a_t_concat == 0] <- 0.001
stan_data$prop_a_t_concat[stan_data$prop_a_t_concat == 1] <- 0.999

temp <- matrix(0.9, nrow = stan_data$S, ncol = stan_data$S)
diag(temp) <- 1
L_S <- t(chol(temp))

temp <- matrix(0.9, nrow = stan_data$U_max, ncol = stan_data$U_max)
diag(temp) <- 1
L_U <- t(chol(temp))

init_param_list <- list(
  #beta0 = rep(-3.0, stan_data$SU),
  beta_mean_raw_a = rnorm(stan_data$DX, mean = log(0.05), sd = 0.01),
  beta_mean_raw_b = rnorm(stan_data$DX, mean = log(0.05), sd = 0.01),
  beta_state_raw_a = rnorm(stan_data$DX * stan_data$S, mean = log(0.05), sd = 0.01),
  beta_state_raw_b = rnorm(stan_data$DX * stan_data$S, mean = log(0.05), sd = 0.01),
  beta_seasonal_raw_a = array(rnorm(stan_data$U * stan_data$DX, mean = log(0.05), sd = 0.01), dim = c(stan_data$U, stan_data$DX, 1)),
  beta_seasonal_raw_b = array(rnorm(stan_data$U * stan_data$DX, mean = log(0.05), sd = 0.01), dim = c(stan_data$U, stan_data$DX, 1)),
  beta_raw_a = rnorm(stan_data$SU * stan_data$DX, mean = log(0.05), sd = 0.01),
  beta_raw_b = rnorm(stan_data$SU * stan_data$DX, mean = log(0.05), sd = 0.01),
  sigma_beta_mean_a = 1.0,
  sigma_beta_mean_b = 1.0,
  sigma_beta_a = 1.0,
  sigma_beta_b = 1.0,
  sigma_beta_state_a = 1.0,
  sigma_beta_state_b = 1.0,
  sigma_beta_seasonal_a = 1.0,
  sigma_beta_seasonal_b = 1.0,
  #sigma_ili = 1.0,
  gamma_a = 0.8,
  gamma_b = 0.8,
  gamma_state_a = 0.8,
  gamma_state_b = 0.8,
  gamma_seasonal_a = 0.8,
  gamma_seasonal_b = 0.8,
  L_S_state_a = L_S,
  L_S_state_b = L_S,
  L_U_seasonal_a = L_U,
  L_U_seasonal_b = L_U,
  L_S_a = L_S,
  L_S_b = L_S,
  L_U_a = L_U,
  L_U_b = L_U,
  state_effect_raw = rnorm(stan_data$S, sd = 0.1),
  sigma_state_effect = 1.0,
  christmas_effect_raw = rnorm(stan_data$SU, sd = 0.1),
  mean_christmas_effect = 0.0,
  sigma_christmas_effect = 1.0,
  outpatient_overdispersion_raw = rnorm(stan_data$total_num_obs_ili, sd = 0.1),
  sigma_outpatient_overdispersion = 1.0,
  outpatient_kappa = 0.1,
  sigma_virology_overdispersion = 1.0,
  virology_kappa = 0.1
  #sigma1_a = 1.0,
  #sigma1_b = 1.0
  #ili_t_raw = rep(0.0, stan_data$total_num_obs_ili)
)
dim(init_param_list$gamma_a) <- 1L
dim(init_param_list$gamma_b) <- 1L
dim(init_param_list$gamma_state_a) <- 1L
dim(init_param_list$gamma_state_b) <- 1L
dim(init_param_list$gamma_seasonal_a) <- 1L
dim(init_param_list$gamma_seasonal_b) <- 1L
dim(init_param_list$beta_raw_a) <- c(stan_data$U_max, stan_data$DX, stan_data$S)
dim(init_param_list$beta_raw_b) <- c(stan_data$U_max, stan_data$DX, stan_data$S)
dim(init_param_list$beta_state_raw_a) <- c(stan_data$DX, stan_data$S)
dim(init_param_list$beta_state_raw_b) <- c(stan_data$DX, stan_data$S)
dim(init_param_list$beta_raw_a) <- c(stan_data$U_max, stan_data$DX, stan_data$S)
dim(init_param_list$beta_raw_b) <- c(stan_data$U_max, stan_data$DX, stan_data$S)

predict_season_filename <- gsub("/", "_", predict_season)

data_dump_file_base <- paste0(
    "data_dump_",
    predict_season_filename,
    "_",
    predict_last_obs_season_week,
    "_depth", tree_depth,
    "_chain", chain_num,
    ".R")
data_dump_file <- paste0(cmdstan_root, "/models/GPthingy/", data_dump_file_base)


stan_rdump(
  names(stan_data),
  file = data_dump_file,
  envir = as.environment(stan_data)
)

init_dump_file_base <- paste0(
    "init_dump_",
    predict_season_filename,
    "_",
    predict_last_obs_season_week,
    "_depth", tree_depth,
    "_chain", chain_num,
    ".R")
init_dump_file <- paste0(cmdstan_root, "/models/GPthingy/", init_dump_file_base)

stan_rdump(
  names(init_param_list),
  file = init_dump_file,
  envir = as.environment(init_param_list)
)

setwd(cmdstan_root)

setwd("models/GPthingy/")
output_file <- paste0(output_path, "samples_depth", tree_depth, "_chain", chain_num, "_last_obs_", predict_season_filename, "_", predict_last_obs_season_week, ".csv")

system(paste0("./GPthingy sample algorithm=hmc engine=nuts max_depth=", tree_depth, " num_samples=1000 num_warmup=500 adapt delta=0.9 data file=", data_dump_file_base, " init=", init_dump_file_base, " output file=", output_file, " refresh=1"))
