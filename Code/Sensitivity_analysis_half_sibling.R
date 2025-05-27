# Needs sorting 

# Plot data from cherry picked pids
cherries <- c("VHX_225", "VHX_457", "VHX_475", "VHX_554", "VHX_551") # "VHX_541","VHX_650"
plot_data(ys = ys_VHX_BPD[cherries], fs = fs_VHX_BPD)
dev.off()

Uniform_Pv3Rs[grepl("VHX_551_", rownames(Uniform_Pv3Rs)), ]
Uniform_Pv3Rs[grepl("VHX_457_", rownames(Uniform_Pv3Rs)), ]

#VHX_532 (an example of half sibs with a high posterior relapse probability)
#BPD_45 (a rare example of a possible reinfection with data on nine markers)
#VHX_622 Example of alternating relapses




