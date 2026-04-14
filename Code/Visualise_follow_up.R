rm(list = ls())
par(mfrow = c(1,2))

# Load VHX and BPD pooled non-genetic data set and convert to data frame
load('../jwatowatson-RecurrentVivax-4870715/RData/TimingModel/Combined_Time_Event.RData')
x <- split(Combined_Time_Data, Combined_Time_Data$patientid)
z <- as.data.frame(do.call(rbind, lapply(x, function(x) x[2,])))

# BPD follow-up: week 1, 2, 4 weeks and monthly thereafter
inds <- grepl("VHX", z)
plot(Combined_Time_Data$Time_since_enrolment[inds], 
     ylab = "VHX time since enrolment", xlab = "Participant index", 
     panel.first = abline(h = c(seq(0, 7*8, 7), seq(7*8, 200, 28)), col = "grey"))

# BPD follow-up: fornightly until 28 days and monthly thereafter
inds <- grepl("BPD", z) 
plot(Combined_Time_Data$Time_since_enrolment[inds], 
     ylab = "BPD time since enrolment",  xlab = "Participant index", 
     panel.first = abline(h = c(seq(0, 28, 7), seq(28, 200, 28)), col = "grey"))



