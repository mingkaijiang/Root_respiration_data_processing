#### read in root respiration raw data and process the data
####
#### created by Mingkai Jiang (m.jiang@westernsydney.edu.au)
####

############################# set up #################################
#### clear wk space
rm(list=ls(all=TRUE))

#### read in necessary stuffs
source("prepare.R")

############################# data processing #####################################
#### read in data
myDF <- read.csv("data/CR3000_rootresp_1_chamber1.dat", skip=3)
colnames(myDF) <- c("DateTime", "Record", "Label", "CurrentLevel", "TPel_1", "Troot_1", "Tbot_1", "CO2")

myDF$Temp <- rowMeans(cbind(myDF$TPel_1, myDF$Troot_1))

#### ignore the Null and CO2 = 0 checks
myDF <- subset(myDF, Label!="Null")
myDF <- subset(myDF, CO2 > 100)

#### get the label list
label.list <- as.character(unique(myDF$Label))

#### create the output DF
outDF <- unique(myDF[,c("Label", "CurrentLevel")])
outDF <- cbind(outDF, NA, NA, NA, NA,
               NA, NA, NA)
colnames(outDF) <- c("Label", "CurrentLevel","Tmean", "Tmin", "Tmax",
                     "linear_slope", "linear_intercept",
                     "log_slope", "log_intercept")


### create pdf
pdf("output/fit_curves.pdf")

#### perform fitting
for (i in label.list) {
    ### subset DF 
    subDF1 <- subset(myDF, Label == i)
    
    ### get the current list within each label list
    current.list <- unique(subDF1$CurrentLevel)
    
    ### loop through the current list
    for (j in current.list) {
        
        ### subset DF
        subDF2 <- subset(subDF1, CurrentLevel==j)
        
        ### get the record and force to starting value of 1
        length <- nrow(subDF2)
        subDF2$Record2 <- c(1:length)
        
        ### check if you have sufficient data to calculate slope
        if (length > 61) {
            ### select the time range
            subDF3 <- subset(subDF2, Record2 >= 61 & Record2 <= 420)
            
            ### get temperature range
            tmin <- min(subDF3$Temp)
            tmax <- max(subDF3$Temp)
            tmean <- mean(subDF3$Temp)
            
            ### assign temperature range
            outDF$Tmean[outDF$Label==i&outDF$CurrentLevel==j] <- tmean
            outDF$Tmin[outDF$Label==i&outDF$CurrentLevel==j] <- tmin
            outDF$Tmax[outDF$Label==i&outDF$CurrentLevel==j] <- tmax
            
            ### linear fit
            fit1 <- lm(CO2~Record2,data=subDF3)
            
            ### assign fit coefficients
            outDF$linear_slope[outDF$Label==i&outDF$CurrentLevel==j] <- coefficients(fit1)[2]
            outDF$linear_intercept[outDF$Label==i&outDF$CurrentLevel==j] <- coefficients(fit1)[1]
            
            ### log fit
            fit2 <- lm(CO2~log(Record2), data=subDF3)
            
            ### assign fit coefficients
            outDF$log_slope[outDF$Label==i&outDF$CurrentLevel==j] <- coefficients(fit2)[2]
            outDF$log_intercept[outDF$Label==i&outDF$CurrentLevel==j] <- coefficients(fit2)[1]
            
            ### make some plots
            with(subDF3, plot(CO2~Record2, xlab = "Record", ylab = expression(CO[2] * " conc.")))
            lines(subDF3$Record2,predict(fit1),type='l',col='blue')
            lines(subDF3$Record2,predict(fit2),col='red')
            legend("topleft", legend=c("linear fit", "log fit"),
                   col=c("blue", "red"), lty=1, cex=1.2)
            title(paste0(i, "_Record_", j))
        }
    }
}

dev.off()

#### define parameters
# get flux (refer to R package 'flux')
v <- 0.25 # chamber volume - sensor, m3.
M = 44.01 #if you want g CO2
p = 101325 # pressure
t = 1  #if transform seconds into hours 1/3600
R <- 8.314 # gas constant for the used units
m <-1 # root mass

#### additional flux calculations
#### return in unit of nmol CO2 g-2 s-1
outDF$flux.lin <- ((outDF$linear_slope * v * M * p) / (t * R * (outDF$Tmean + 273.15) * m)) / M*1000 #nmol CO2 g-2 s-1

### save file
write.csv(outDF, "output/Rroot_coefficients_summary.csv", row.names=F)


#### fit exponential slope to get final estimates
outDF2 <- data.frame(unique(outDF$Label), NA, NA, NA, NA)
colnames(outDF2) <- c("Label", "a", "b", "r2", "p")


pdf("output/exp_fitting_to_temperature.pdf")

for (i in unique(outDF2$Label)) {
    subDF <- subset(outDF, Label == i)
    
    ### check number of data points
    length <- nrow(subDF)
    
    if (length > 5) {

        fit_nls <- lm(log(flux.lin) ~ Tmean,
                      data = subDF)
        
        outDF2$a[outDF2$Label==i] <- coefficients(fit_nls)[1]
        outDF2$b[outDF2$Label==i] <- coefficients(fit_nls)[2]
        outDF2$r2[outDF2$Label==i] <- summary(fit_nls)$adj.r.squared
        outDF2$p[outDF2$Label==i] <- summary(fit_nls)$coefficients[4]
        
        with(subDF, plot(linear_slope~Tmean))
        lines(subDF$Tmean,predict(fit_nls),lty=2,col="red",lwd=3)
        title(paste0(i))
    }
}

dev.off()

write.csv(outDF2, "output/temperature_function.csv", row.names=F)




