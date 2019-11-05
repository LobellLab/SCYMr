#' Generic accuracy assessment function
#'
#' Calculates a suite of agreement statistics for two sets of paired observations.
#' @param x_pred vector of predictions
#' @param y_obs vector of validation observations
#' @keywords accuracy metrics
#' @export
#' @examples
#' # get statistics for all observations
#' agreement <- calcStats(df$x_pred,df$y_obs)
#'
#' # get statistics for groups of observations
#' library(dplyr)
#' agreementByGroup <- df %>%
#'   group_by(group) %>%
#'   group_map(~calcStats(.$x_pred,.$y_obs))

calcStats <- function(x_pred1, y_obs1){
  df <- data.frame(x_pred = x_pred1,
                   y_obs = y_obs1)
  # run a lm
  lm1 <- lm(y_obs ~ x_pred, data = df)

  # stats
  r2 <- summary(lm1)$r.squared
  rmse <- sqrt(mean(lm1$residuals^2))
  r <- cor(df$y_obs, df$x_pred)
  m <- coef(lm1)[2]
  int <- coef(lm1)[1]
  rmspe <- sqrt(mean(((df$y_obs - df$x_pred)/df$y_obs)^2))*100
  vecv <- (1 - sum((df$y_obs - df$x_pred)^2)/sum((df$y_obs - mean(df$y_obs))^2)) * 100
  mae <- sum(abs(df$y_obs-df$x_pred))/nrow(df)
  me <- mean(df$x_pred-df$y_obs)

  # package into df
  outdf <- data.frame(r2 = r2, rmse = rmse, r = r, m = m, int = int,
                      rmspe = rmspe, vecv = vecv, mae = mae, me = me)
  return(outdf)
}







