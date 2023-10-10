set_FP_data <- function(outcome_names , 
                        study_names, 
                        treatment_names) {
  
  for(outcome_name in outcome_names) {
    FP_data[[outcome_name]] <- list()
    for(study_name in study_names[[outcome_name]]) {
      FP_data[[outcome_name]][[study_name]] <- list()
      for(treatment_name in treatment_names[[study_name]]) {
        natrisk <- read_excel(paste0(km_data_directory, "/", study_name, " ", treatment_name, " ", outcome_name, ".xlsx"), sheet = 2)
        ipd <- read_csv(paste0(ipd_directory, "/", outcome_name, "_", study_name, "_", treatment_name, "_ipd", ".csv"))
        FP_data_file <- paste0(FP_data_directory, "/", outcome_name, "_", study_name, "_", treatment_name, "_FPdata.csv")
        
        s <- rep(1, nrow(natrisk)-1)  # index of the study by time interval
        a <- rep(1, nrow(natrisk)-1)  # index of arm by time interval
        dt <- diff(natrisk$Time)      # time intervals length
        time_max <- max(natrisk$Time) # maximum time of follow-up
        cutpoints <- natrisk$Time
        ipd <- ipd %>% filter(time <= time_max)
        r <- ipd %>% mutate(ranges = cut2(time, cutpoints)) %>% # number of events by time interval 
          group_by(ranges) %>% 
          summarise(r = sum(event))# %>% select(r)
        r <- r$r
        
        z <- natrisk$Natrisk[-length(natrisk$Natrisk)] # number of patients at risk at each interval
        time_j <- natrisk$Time[-1] # start of every time interval
        
        FP_temp <- qpcR:::cbind.na(s, a, dt, r, z, time_j)
        
        FP_data[[outcome_name]][[study_name]][[treatment_name]] <- FP_temp
        write.csv(FP_temp, FP_data_file)
      }
    }
  }
}