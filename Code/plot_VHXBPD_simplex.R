plot_VHXBPD_simplex <- function(Uniform_xy, TimeToEvent_xy){
  
  # Get treatment info from Combined_Time data because MS_final doesn't contain
  # episode identifiers for 6 recurrences that could not be estimated
  load("../RData/Combined_Time_Event.RData")
  Combined_Time_Data <- as.data.frame(Combined_Time_Data)
  rownames(Combined_Time_Data) <- paste(Combined_Time_Data$patientid, 
                                        Combined_Time_Data$episode, sep = "_")
  if (!all(colnames(Uniform_xy) %in% rownames(Combined_Time_Data))) stop()
  if (!all(colnames(TimeToEvent_xy) %in% rownames(Combined_Time_Data))) stop()
  
  
  par(mfrow = c(2,2), mar = c(1,0,1,0))
  PMQ_uniform <- grepl("PMQ", Combined_Time_Data[colnames(Uniform_xy), "arm_num"])
  PMQ_timetoevent <- grepl("PMQ", Combined_Time_Data[colnames(TimeToEvent_xy), "arm_num"])
  
  # Plot PMQ Uniform
  plot_simplex(v.labels = c("Recrudescence", "Relapse", "Reinfection"))
  title(main = sprintf("PMQ+ uniform: %s estimates", sum(PMQ_uniform)), 
        line = -0.5)
  for(i in which(PMQ_uniform)){
    points(x = Uniform_xy["x",i], y = Uniform_xy["y",i], 
           col = adjustcolor("black", alpha.f = 0.3), 
           pch = 17-Uniform_xy["joint", i])}
  if(any(Uniform_xy["joint",which(PMQ_uniform)] == 0)){
    legend(x = 0.1, y = 0.5, bty = "n", pch = 16:17, 
           col = adjustcolor("black", alpha.f = 0.3), 
           legend = c("Jointly modelled data","Pairwise modelled data"))}
   
  # Plot BPD TimeToEvent
  plot_simplex(v.labels = c("Recrudescence", "Relapse", "Reinfection"))
  title(main = sprintf("PMQ+ time-to-event: %s estimates", sum(PMQ_timetoevent)), 
        line = -0.5)
  for(i in which(PMQ_timetoevent)){
    points(x = TimeToEvent_xy["x",i], y = TimeToEvent_xy["y",i], 
           col = adjustcolor("black", alpha.f = 0.3), 
           pch = 17-TimeToEvent_xy["joint",i])}
  if(any(TimeToEvent_xy["joint",which(PMQ_timetoevent)] == 0)){
  legend(x = 0.1, y = 0.5, bty = "n", pch = 16:17, 
         col = adjustcolor("black", alpha.f = 0.3), 
         legend = c("Jointly modelled data","Pairwise modelled data"))}
  
  # Plot VHX Uniform
  plot_simplex(v.labels = c("Recrudescence", "Relapse", "Reinfection"))
  title(main = sprintf("no PMQ uniform: %s estimates", sum(!PMQ_uniform)), 
        line = -0.5)
  for(i in which(!PMQ_uniform)){
    points(x = Uniform_xy["x",i], y = Uniform_xy["y",i], 
           col = adjustcolor("black", alpha.f = 0.3), 
           pch = 17-Uniform_xy["joint", i])}
  if(any(Uniform_xy["joint",which(!PMQ_uniform)] == 0)){
  legend(x = 0.1, y = 0.5, bty = "n", pch = 16:17, 
         col = adjustcolor("black", alpha.f = 0.3), 
         legend = c("Jointly modelled data","Pairwise modelled data"))}
  
  # Plot No PMQ TimeToEvent
  plot_simplex(v.labels = c("Recrudescence", "Relapse", "Reinfection"))
  title(main = sprintf("No PMQ time-to-event: %s estimates", sum(!PMQ_timetoevent)), 
        line = -0.5)
  for(i in which(!PMQ_timetoevent)){
    points(x = TimeToEvent_xy["x",i], y = TimeToEvent_xy["y",i], 
           col = adjustcolor("black", alpha.f = 0.3), 
           pch = 17-TimeToEvent_xy["joint", i])}
  if(any(TimeToEvent_xy["joint", which(!PMQ_timetoevent)] == 0)){
    legend(x = 0.1, y = 0.5, bty = "n", pch = 16:17, 
           col = adjustcolor("black", alpha.f = 0.3), 
           legend = c("Jointly modelled data","Pairwise modelled data"))}

}




