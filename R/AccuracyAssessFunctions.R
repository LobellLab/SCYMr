library(here)
library(tidyverse)
library(MLmetrics)

# load scym and derive mean county yield in t/ha
getSCYMTable <- function(scymTablePath, scaled = TRUE, scaleBy = NULL){
  if (scaled) {
    scym = read_csv(scymTablePath) %>%
      filter(masterid != '0') %>%
      rename(scym_biomass = yield_mean) %>%
      mutate(fips5 = str_pad(masterid, 5, pad='0'),
             yield_scym_mt_ha = scym_biomass) %>%
      select(c(fips5, year, yield_scym_mt_ha, scym_biomass))
    return(scym)
  } else {
    scym = read_csv(scymTablePath) %>%
      filter(masterid != '0') %>%
      rename(scym_biomass = yield_mean) %>%
      mutate(fips5 = str_pad(masterid, 5, pad='0'),
             yield_scym_mt_ha = scym_biomass * scaleBy) %>%
      select(c(fips5, year, yield_scym_mt_ha, scym_biomass))
    return(scym)
  }
}

# load nass
getNASSTable <- function(nassTablePath) {
  nass <- read_csv(nassTablePath) %>%
    mutate(fips5 = str_pad(fips5, 5, pad='0'))
  return(nass)
}

#get the combo df
getCombinedDF <- function(scymTablePath, scymTableScaled = TRUE, scaleBy = NULL, nassTablePath) {
  scym_table <- getSCYMTable(scymTablePath, scymTableScaled, scaleBy)
  nass_table <- getNASSTable(nassTablePath)
  combo <- scym_table %>%
    inner_join(nass_table, by=c('fips5','year'))
  return(combo)
}

# agreement stats by state (all years)
calcStats_state = function(combined_df){
  stats = combined_df %>%
    na.omit() %>%
    select(c(state_name, yield_scym_mt_ha, nassyield_mt_ha)) %>%
    group_by(state_name) %>%
    mutate(R2 = MLmetrics::R2_Score(yield_scym_mt_ha,nassyield_mt_ha),
           RMSE = MLmetrics::RMSE(yield_scym_mt_ha,nassyield_mt_ha),
           r = cor(yield_scym_mt_ha,nassyield_mt_ha),
           RMSPE = MLmetrics::RMSPE(yield_scym_mt_ha,nassyield_mt_ha)*100,
           vecv = spm::vecv(obs = nassyield_mt_ha, yield_scym_mt_ha),
           m = coef(lm(nassyield_mt_ha ~ yield_scym_mt_ha))[2],
           int = coef(lm(nassyield_mt_ha ~ yield_scym_mt_ha))[1]) %>%
    summarize(R2 = min(R2),
              RMSE = min(RMSE),
              r = min(r),
              RMSPE = min(RMSPE),
              vecv = min(vecv),
              m = min(m),
              int = min(int))
}


#make agreement plots by state

# create a manual grid for 9 corn belt states. the 'code' should match
# whatever variable you facet by (so yours might be full state names)
mygrid <- data.frame(
  row = c( 1, 1, 1, 1, 2, 2, 2, 2, 3),
  col = c( 4, 2, 1, 3, 3, 4, 2,  5,   2),
  code = c( "MICHIGAN", "MINNESOTA", "SOUTH DAKOTA", "WISCONSIN", "ILLINOIS", "INDIANA", "IOWA",  "OHIO",   "MISSOURI"),
  name = c( "Michigan", "Minnesota", "South Dakota", "Wisconsin", "Illinois", "Indiana", "Iowa",  "Ohio",   "Missouri"),
  stringsAsFactors = FALSE
)

stateStatsPlotter <- function(combined_df, statestatdf, ncol=4, title, axisScale = 15){
  ggplot(combined_df %>% na.omit(),
         aes(x = yield_scym_mt_ha, y = nassyield_mt_ha)) +
    geom_point(cex = .7, color = 'skyblue3') +
    geom_text(data = statestatdf, x = (10.2/15)*axisScale, y = (2.2/15)*axisScale, size = 3,
              aes(label = paste0('          r = ',base::round(r, digits = 2), "\n",
                                 '  RMSE = ',round(RMSE, digits = 2), "\n",
                                 'RMSPE = ', round(RMSPE), "%"))) +
    geom_text(data = statestatdf, x = (3.5/15)*axisScale, y = (12.2/15)*axisScale, size = 3,
              aes(label = paste0('vecv = ',base::round(vecv, digits = 2), '\n',
                                 'yint = ',base::round(int, digits = 2), '  \n',
                                 'm = ', base::round(m, digits = 2), '   ')) )+
    geom_abline(slope = 1, intercept = 0, linetype='dashed') +
    #facet_wrap(~state_name, ncol=ncol) +
    facet_geo(~state_name, grid = mygrid) +
    coord_equal(xlim=c(0,axisScale), ylim=c(0,axisScale)) +
    xlab('Scym yield estimates (t/ha)') + ylab('NASS county-level yield (t/ha)') +
    theme_bw() + ggtitle(title)
}

# agreement stats by year
calcStats_year = function(combined_df){
  yearstats = combined_df %>%
    na.omit() %>%
    #filter(state_name %in% c('Indiana','Illinois','Iowa')) %>%
    select(c(state_name, year, yield_scym_mt_ha, nassyield_mt_ha)) %>%
    group_by(state_name, year) %>%
    mutate(R2 = MLmetrics::R2_Score(yield_scym_mt_ha,nassyield_mt_ha),
           RMSE = MLmetrics::RMSE(yield_scym_mt_ha,nassyield_mt_ha),
           r = cor(yield_scym_mt_ha,nassyield_mt_ha),
           RMSPE = MLmetrics::RMSPE(yield_scym_mt_ha,nassyield_mt_ha)*100,
           vecv = spm::vecv(obs = nassyield_mt_ha, yield_scym_mt_ha),
           m = coef(lm(nassyield_mt_ha ~ yield_scym_mt_ha))[2],
           int = coef(lm(nassyield_mt_ha ~ yield_scym_mt_ha))[1]) %>%
    summarize(R2 = min(R2),
              RMSE = min(RMSE),
              r = min(r),
              RMSPE = min(RMSPE),
              vecv = min(vecv),
              m = min(m),
              int = min(int))
}


#make agreement plots by year
yearStatsPlotter <- function(statename, dotColor, combined_df, yearstatdf, ncol=5, title, axisScale = 15){
  ggplot(combined_df %>% na.omit() %>% filter(state_name == statename),
         aes(x = yield_scym_mt_ha, y = nassyield_mt_ha)) +
    geom_point(cex = .7, color = dotColor) +
    geom_text(data = yearstatdf %>% filter(state_name == statename),
              x = (10.2/15)*axisScale, y = (2.2/15)*axisScale, size = 3,
              aes(label = paste0('          r = ',base::round(r, digits = 2), "\n",
                                 '  RMSE = ',round(RMSE, digits = 2), "\n",
                                 'RMSPE = ', round(RMSPE), "%"))) +
    geom_text(data = yearstatdf %>% filter(state_name == statename),
              x = (3.5/15)*axisScale, y = (12.5/15)*axisScale, size = 3,
              aes(label = paste0('vecv = ',base::round(vecv, digits = 2), '\n',
                                 'yint = ',base::round(int, digits = 2), '  \n',
                                 'm = ', base::round(m, digits = 2), '   ')) )+
    geom_abline(slope = 1, intercept = 0, linetype='dashed') +
    facet_wrap(~year, ncol=ncol) +
    coord_equal(xlim=c(0,axisScale), ylim=c(0,axisScale)) +
    xlab('Scym yield estimates (t/ha)') + ylab('NASS county-level yield (t/ha)') +
    theme_bw() + ggtitle(title)
}

overallR2Calc <- function(combined_df) {
  stats = combined_df %>%
    na.omit() %>%
    select(c(state_name, yield_scym_mt_ha, nassyield_mt_ha)) %>%
    mutate(R2 = MLmetrics::R2_Score(yield_scym_mt_ha,nassyield_mt_ha))
  stats$R2[1]
}

overallRMSECalc <- function(combined_df) {
  stats = combined_df %>%
    na.omit() %>%
    select(c(state_name, yield_scym_mt_ha, nassyield_mt_ha)) %>%
    mutate(RMSE = MLmetrics::RMSE(yield_scym_mt_ha,nassyield_mt_ha))
  stats$RMSE[1]
}
