plot_VHXBPD_simplex <- function(Uniform_xy, TimeToEvent_xy){
  
  par(mfrow = c(2,2), mar = c(0,0,0,0))
  BPD_uniform <- grepl("BPD", colnames(Uniform_xy))
  BPD_timetoevent <- grepl("BPD", colnames(TimeToEvent_xy))
  
  # Plot BPD Uniform
  plot_simplex()
  title(main = sprintf("BPD uniform: %s estimates", sum(BPD_uniform)), 
        line = -1.5)
  for(i in which(BPD_uniform)){
    points(x = Uniform_xy["x",i], y = Uniform_xy["y",i], 
           col = adjustcolor("black", alpha.f = 0.3), 
           pch = 17-Uniform_xy["joint", i])}
  if(any(Uniform_xy["joint",which(BPD_uniform)] == 0)){
    legend(x = 0.1, y = 0.5, bty = "n", pch = 16:17, 
           col = adjustcolor("black", alpha.f = 0.3), 
           legend = c("Jointly modelled data","Pairwise modelled data"))}
   
  # Plot BPD TimeToEvent
  plot_simplex()
  title(main = sprintf("BPD time-to-event: %s estimates", sum(BPD_timetoevent)), 
        line = -1.5)
  for(i in which(BPD_timetoevent)){
    points(x = TimeToEvent_xy["x",i], y = TimeToEvent_xy["y",i], 
           col = adjustcolor("black", alpha.f = 0.3), 
           pch = 17-TimeToEvent_xy["joint",i])}
  if(any(TimeToEvent_xy["joint",which(BPD_timetoevent)] == 0)){
  legend(x = 0.1, y = 0.5, bty = "n", pch = 16:17, 
         col = adjustcolor("black", alpha.f = 0.3), 
         legend = c("Jointly modelled data","Pairwise modelled data"))}
  
  # Plot VHX Uniform
  plot_simplex()
  title(main = sprintf("VHX uniform: %s estimates", sum(!BPD_uniform)), 
        line = -1.5)
  for(i in which(!BPD_uniform)){
    points(x = Uniform_xy["x",i], y = Uniform_xy["y",i], 
           col = adjustcolor("black", alpha.f = 0.3), 
           pch = 17-Uniform_xy["joint", i])}
  if(any(Uniform_xy["joint",which(!BPD_uniform)] == 0)){
  legend(x = 0.1, y = 0.5, bty = "n", pch = 16:17, 
         col = adjustcolor("black", alpha.f = 0.3), 
         legend = c("Jointly modelled data","Pairwise modelled data"))}
  
  # Plot VHX TimeToEvent
  plot_simplex()
  title(main = sprintf("VHX time-to-event: %s estimates", sum(!BPD_timetoevent)), 
        line = -1.5)
  for(i in which(!BPD_timetoevent)){
    points(x = TimeToEvent_xy["x",i], y = TimeToEvent_xy["y",i], 
           col = adjustcolor("black", alpha.f = 0.3), 
           pch = 17-TimeToEvent_xy["joint", i])}
  if(any(TimeToEvent_xy["joint", which(!BPD_timetoevent)] == 0)){
    legend(x = 0.1, y = 0.5, bty = "n", pch = 16:17, 
           col = adjustcolor("black", alpha.f = 0.3), 
           legend = c("Jointly modelled data","Pairwise modelled data"))}

}




