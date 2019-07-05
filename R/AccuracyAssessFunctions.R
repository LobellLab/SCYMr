library(here)
library(tidyverse)
library(MLmetrics)

# load scym and derive mean county yield in t/ha
getSCYMTable <- function(scymTablePath, scaled){
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
             yield_scym_mt_ha = scym_biomass * .45 * .01) %>%
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
getCombinedDF <- function(scymTablePath, scymTableScaled, nassTablePath) {
  scym_table <- getSCYMTable(scymTablePath, scymTableScaled)
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
    mutate(RMSE = MLmetrics::RMSE(yield_scym_mt_ha,nassyield_mt_ha),
           r = cor(yield_scym_mt_ha,nassyield_mt_ha),
           RMSPE = MLmetrics::RMSPE(yield_scym_mt_ha,nassyield_mt_ha)*100,
           vecv = spm::vecv(obs = nassyield_mt_ha, yield_scym_mt_ha),
           m = coef(lm(nassyield_mt_ha ~ yield_scym_mt_ha))[2],
           int = coef(lm(nassyield_mt_ha ~ yield_scym_mt_ha))[1]) %>%
    summarize(RMSE = min(RMSE),
              r = min(r),
              RMSPE = min(RMSPE),
              vecv = min(vecv),
              m = min(m),
              int = min(int))
}


#make agreement plots by state
stateStatsPlotter <- function(combined_df, statestatdf, ncol=4, title){
  ggplot(combined_df %>% na.omit(),
         aes(x = yield_scym_mt_ha, y = nassyield_mt_ha)) +
    geom_point(cex = .7, color = 'skyblue3') +
    geom_text(data = statestatdf, x = 10.2, y = 2.2, size = 3,
              aes(label = paste0('          r = ',base::round(r, digits = 2), "\n",
                                 '  RMSE = ',round(RMSE, digits = 2), "\n",
                                 'RMSPE = ', round(RMSPE), "%"))) +
    geom_text(data = statestatdf, x = 3.5, y = 12.2, size = 3,
              aes(label = paste0('vecv = ',base::round(vecv, digits = 2), '\n',
                                 'yint = ',base::round(int, digits = 2), '  \n',
                                 'm = ', base::round(m, digits = 2), '   ')) )+
    geom_abline(slope = 1, intercept = 0, linetype='dashed') +
    facet_wrap(~state_name, ncol=ncol) +
    coord_equal(xlim=c(0,15), ylim=c(0,15)) +
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
    mutate(RMSE = MLmetrics::RMSE(yield_scym_mt_ha,nassyield_mt_ha),
           r = cor(yield_scym_mt_ha,nassyield_mt_ha),
           RMSPE = MLmetrics::RMSPE(yield_scym_mt_ha,nassyield_mt_ha)*100,
           vecv = spm::vecv(obs = nassyield_mt_ha, yield_scym_mt_ha),
           m = coef(lm(nassyield_mt_ha ~ yield_scym_mt_ha))[2],
           int = coef(lm(nassyield_mt_ha ~ yield_scym_mt_ha))[1]) %>%
    summarize(RMSE = min(RMSE),
              r = min(r),
              RMSPE = min(RMSPE),
              vecv = min(vecv),
              m = min(m),
              int = min(int))
}


#make agreement plots by year
yearStatsPlotter <- function(statename, dotColor, combined_df, yearstatdf, ncol=5, title){
  ggplot(combined_df %>% na.omit() %>% filter(state_name == statename),
         aes(x = yield_scym_mt_ha, y = nassyield_mt_ha)) +
    geom_point(cex = .7, color = dotColor) +
    geom_text(data = yearstatdf %>% filter(state_name == statename),
              x = 10.2, y = 2.2, size = 3,
              aes(label = paste0('          r = ',base::round(r, digits = 2), "\n",
                                 '  RMSE = ',round(RMSE, digits = 2), "\n",
                                 'RMSPE = ', round(RMSPE), "%"))) +
    geom_text(data = yearstatdf %>% filter(state_name == statename),
              x = 3.5, y = 12.5, size = 3,
              aes(label = paste0('vecv = ',base::round(vecv, digits = 2), '\n',
                                 'yint = ',base::round(int, digits = 2), '  \n',
                                 'm = ', base::round(m, digits = 2), '   ')) )+
    geom_abline(slope = 1, intercept = 0, linetype='dashed') +
    facet_wrap(~year, ncol=ncol) +
    coord_equal(xlim=c(0,15), ylim=c(0,15)) +
    xlab('Scym yield estimates (t/ha)') + ylab('NASS county-level yield (t/ha)') +
    theme_bw() + ggtitle(title)
}



