#function to analyse spatial pattern, returns list with 3 tibbles: densities, ripley'sK values, mark corr values
data_params <- function(scenario_run, file_path){
  
  #open all parameters files for current scenario
  Parameters_files <- list.files(path=file_path, pattern=sprintf('Parameters %s *', scenario_run))
  #initialising dataframe for current scenario_run
  params_scen <- tibble(pf=factor(), scenario=factor(), active=factor(), CA_w_env=factor(),
                        CA_day=factor(), response_func=factor(),choosers_on_followers=numeric(),
                        bi_CA=factor(), weight_day=numeric())
  #for every scenario, parameters are extracted from parameters-files
  for (file_title in Parameters_files){
    title_full = paste(file_path, file_title, sep="/")
    d = read.table(title_full, header=T)
    d[[6]]=ifelse(d[[6]]=='None', NA, d[[6]])
    params_scen <- add_row(params_scen,
                           pf=file_title, scenario=d[[1]],
                           active=ifelse(d[[1]] %in% c('random', 'bottom-up', 'prev-exp-BU'),NA, as.character(d[[2]])),
                           CA_w_env=ifelse(d[[1]] %in% c('random', 'bottom-up', 'prev-exp-BU'),NA, as.character(d[[3]])),
                           CA_day=ifelse(d[[1]] %in% c('random', 'bottom-up', 'prev-exp-BU'),NA, as.character(d[[4]])),
                           response_func=ifelse(d[[1]] %in% c('random', 'bottom-up', 'prev-exp-BU'), NA, as.character(d[[5]])),
                           choosers_on_followers=ifelse(d[[1]]!='CF', NA,d[[6]]),
                           bi_CA=ifelse(d[[1]]!='bimodal', NA, as.character(d[[7]])),
                           weight_day = ifelse(d[[1]] %in% c('random', 'bottom-up', 'prev-exp-BU'), NA, d[[8]]))
  }
  return(params_scen)
}