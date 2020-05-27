
#' Create a quiet function
#'
#' @param .f A function to be silenced
#' @param silenceMessages Whether messages should be silenced. Defaults to TRUE.
#' @param silenceWarnings Whether warnings should be silenced Defaults to true.
#'
#' @return A wrapped function with  messages and/or warnings suppressed. 
#' If both messages and warnings are suppressed, the function is returned unmodified. 
#' @export
#'
#' @examples
#' quiet_lmer <- quiet(lme4::lmer)
#' lme4::lmer(formula = Petal.Width ~ Petal.Length + (Petal.Length | Species), 
#'            data = iris) # Returns singularity warning
quiet <- function(.f, silenceMessages = TRUE, silenceWarnings = TRUE) {
	.f <- purrr::as_mapper(.f)
	
	warnings <- character()
	wHandler <- function(w) {
		warnings <<- c(warnings, w$message)
		invokeRestart("muffleWarning")
	}	
	
	messages <- character()
	mHandler <- function(m) {
		messages <<- c(messages, m$message)
		invokeRestart("muffleMessage")
	}
	
	if(silenceMessages & silenceWarnings) {
		function(...) withCallingHandlers(.f(...), warning = wHandler, message = mHandler)
	} else if(silenceMessages & !silenceWarnings) {
		function(...) withCallingHandlers(.f(...), message = mHandler)
	} else if(!silenceMessages & silenceWarnings) {
		function(...) withCallingHandlers(.f(...), warning = wHandler)
	} else {
		warning("Nothing was silenced. This function was useless.")
		.f(...)
	}
}

#' R-squared between or within groups
#' 
#' Function for calculating \eqn{R^2} either within or between groups.
#' 
#' By specifying a some group in the \code{group} argument of this function, \eqn{R^2} is calculated within groups. 
#' Note that this is only meaningful for groups with more than one observation and where variance exists.
#' 
#' @param pred a vector containing numeric predictions (i.e. estimated values).
#' @param obs a vector containing numeric observations (i.e. truth values).
#' @param group An optional group by which to calculate \eqn{R^2}. If left empty. \eqn{R^2} is calculated between groups.
#'
#' @return calculated \eqn{R^2} as a single numeric value.
#' @export
r2 <- function(pred, obs, group = NULL) {
	r2_fun <- rlang::expr(1-(mean((pred - obs)^2) / var(obs, na.rm = T)))
	
	if(is.null(group)) {
		return(rlang::eval_tidy(r2_fun))
	}
	else {
		data <- data.frame(group = group,
											 pred = pred,
											 obs = obs)
		data %>% 
			dplyr::group_by(group) %>% 
			dplyr::filter(n() > 1) %>% 
			dplyr::summarise(r2.inner = rlang::eval_tidy(r2_fun), n = n()) %>% 
			dplyr::filter(!is.infinite(r2.inner)) %>% 
			dplyr::summarise(r2 = weighted.mean(r2.inner, n, na.rm = T)) %>% 
			dplyr::pull(r2) %>% 
			return()
	}
}

#' Standard error of predictions between or within groups
#' 
#' Function for calculating the standard error for predictions either within or between groups.
#' 
#' By specifying a some group in the \code{group} argument of this function, the standard error is calculated within groups. 
#' Note that this is only meaningful for groups with more than one observation and where variance exists.
#' 
#' @param pred a vector containing numeric predictions (i.e. estimated values).
#' @param obs a vector containing numeric observations (i.e. truth values).
#' @param group An optional group by which to calculate the standard error. If left empty. the standard error is calculated between groups.
#'
#' @return calculated standard error as a single numeric value.
#' @export
se <- function(pred, obs, group = NULL) {
	se_fun <- rlang::expr(sd((pred - obs)^2) / sqrt(length(obs)) / var(obs))
	
	if(is.null(group)) {
		return(rlang::eval_tidy(se_fun))
	}
	else {
		data <- data.frame(group = group,
											 pred = pred,
											 obs = obs)
		data %>% 
			dplyr::group_by(group) %>% 
			dplyr::filter(n() > 1) %>% 
			dplyr::summarise(se.inner = rlang::eval_tidy(se_fun), n = n()) %>% 
			dplyr::summarise(se = weighted.mean(r2.inner, n, na.rm = T)) %>% 
			dplyr::pull(se) %>% 
			return()
	}
}

#' Intraclass correlation coefficient (ICC) using mixed effects models
#' 
#' Calculate the intraclass correlation coefficient (ICC) using mixed effects (also known as multilevel) models.
#'
#' @param data a data frame containing the data needed for calculated the ICC.
#' @param ... unquoted variables in the \code{data} argument for which the ICC needs to be calculated.
#' @param group an unquoted variable specifying the group by which to calculate the ICC.
#'
#' @return a data frame containing two columns: A \code{Variable} column for which the ICC was calculated and a \code{ICC} which contains the ICC value.
#' @export
icc <- function(data, ..., group) {
	if(rlang::dots_n(...) == 0) {
		stop("No variables selected to be evaluated")
	}
	
	group <- rlang::enquo(group)
	vars <- rlang::enquos(...)
	
	data <- data %>% 
		dplyr::select(!!group, !!!vars)
	
	group <- colnames(data)[1]
	cols <- colnames(data)[2:ncol(data)]
	res <- data.frame(Variable = cols, ICC = NA)
	
	quiet_lmer <- quiet(purrr::possibly(lme4::lmer, NA))
	for(i in 1:length(cols)) {
		lmer.res <- quiet_lmer(formula = as.formula(paste0("`", cols[i], "` ~ 1 + (1 | ", group, " )")), 
										 data = data)
		intraccc <- as.data.frame(lme4::VarCorr(lmer.res))$sdcor^2    
		res$ICC[i] <- intraccc[1] / sum(intraccc)
	}
	res
}
