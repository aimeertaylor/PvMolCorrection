plot_VHXBPD_simplex <- function(Uniform_xy, TimeToEvent_xy){
  
  pt_col <- adjustcolor("black", alpha.f = 0.3)
  legend_arg <- list(x = 0.1, y = 0.5, bty = "n", pch = 16, col = pt_col)
  
  # Get treatment info from Combined_Time data because MS_final doesn't contain
  # episode identifiers for 6 recurrences that could not be estimated
  load("../RData/Combined_Time_Event.RData")
  Combined_Time_Data <- as.data.frame(Combined_Time_Data)
  rownames(Combined_Time_Data) <- paste(Combined_Time_Data$patientid, 
                                        Combined_Time_Data$episode, sep = "_")
  if (!all(colnames(Uniform_xy) %in% rownames(Combined_Time_Data))) stop()
  if (!all(colnames(TimeToEvent_xy) %in% rownames(Combined_Time_Data))) stop()
  
  
  par(mfrow = c(2,2), mar = c(1,0,1,0))
  PQ_uniform <- grepl("PMQ", Combined_Time_Data[colnames(Uniform_xy), "arm_num"])
  PQ_timetoevent <- grepl("PMQ", Combined_Time_Data[colnames(TimeToEvent_xy), "arm_num"])
  
  # Plot PQ Uniform
  plot_simplex(v.labels = c("Recrudescence", "Relapse", "Reinfection"))
  title(main = sprintf("PQ-treated uniform: %s estimates", sum(PQ_uniform)), 
        line = -0.5)
  for(i in which(PQ_uniform)){
    points(x = Uniform_xy["x",i], y = Uniform_xy["y",i], col = pt_col, 
           pch = 17-Uniform_xy["joint", i])}
  
  if(all(Uniform_xy["joint", which(PQ_uniform)] == 1)){
    legend(x = 0.1, y = 0.5, bty = "n", pch = 16, col = pt_col, legend = "All joint")
  } else {
    legend(x = 0.1, y = 0.5, bty = "n", pch = 16:17, col = pt_col, 
           legend = c("Joint","Approximate"))      
  }
  
  
  # Plot BPD TimeToEvent
  plot_simplex(v.labels = c("Recrudescence", "Relapse", "Reinfection"))
  title(main = sprintf("PQ-treated time-to-event: %s estimates", sum(PQ_timetoevent)), 
        line = -0.5)
  for(i in which(PQ_timetoevent)){
    points(x = TimeToEvent_xy["x",i], y = TimeToEvent_xy["y",i], col = pt_col, 
           pch = 17-TimeToEvent_xy["joint",i])}
  
  if(all(TimeToEvent_xy["joint",which(PQ_timetoevent)] == 1)){
    legend(x = 0.1, y = 0.5, bty = "n", pch = 16, col = pt_col, 
           legend = "All joint")
  } else {
    legend(x = 0.1, y = 0.5, bty = "n", pch = 16:17, col = pt_col, 
           legend = c("Joint","Approximate"))    
  }
  
  
  # Plot VHX Uniform
  plot_simplex(v.labels = c("Recrudescence", "Relapse", "Reinfection"))
  title(main = sprintf("PQ-untreated uniform: %s estimates", sum(!PQ_uniform)), 
        line = -0.5)
  for(i in which(!PQ_uniform)){
    points(x = Uniform_xy["x",i], y = Uniform_xy["y",i], col = pt_col, 
           pch = 17-Uniform_xy["joint", i])}
  
  if(all(Uniform_xy["joint",which(!PQ_uniform)] == 1)){
    legend(x = 0.1, y = 0.5, bty = "n", pch = 16, col = pt_col, legend = "All joint")
  } else {
    legend(x = 0.1, y = 0.5, bty = "n", pch = 16:17, col = pt_col, 
           legend = c("Joint","Approximate"))
  }
  
  # Plot PQ-untreated TimeToEvent
  plot_simplex(v.labels = c("Recrudescence", "Relapse", "Reinfection"))
  title(main = sprintf("PQ-untreated time-to-event: %s estimates", sum(!PQ_timetoevent)), 
        line = -0.5)
  for(i in which(!PQ_timetoevent)){
    points(x = TimeToEvent_xy["x",i], y = TimeToEvent_xy["y",i], col = pt_col, 
           pch = 17-TimeToEvent_xy["joint", i])}
  if(all(TimeToEvent_xy["joint", which(!PQ_timetoevent)] == 1)){
    legend(x = 0.1, y = 0.5, bty = "n", pch = 16, col = pt_col, 
           legend = "All joint") 
  } else {
    legend(x = 0.1, y = 0.5, bty = "n", pch = 16:17, col = pt_col, 
           legend = c("Joint","Approximate"))
  }
  
}




