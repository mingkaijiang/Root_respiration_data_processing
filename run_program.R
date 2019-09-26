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
        subDF2 <- subset(subDF1, CurrentLevel=j)
        
        ### get the record and force to starting value of 1
        m <- min(subDF2$Record)
        subDF2$Record2 <- subDF2$Record - m + 1
        
        ### get temperature range
        tmin <- min(subDF2$TPel_1)
        tmax <- max(subDF2$TPel_1)
        tmean <- mean(subDF2$TPel_1)
        
        ### assign temperature range
        outDF$Tmean[outDF$Label==i&outDF$CurrentLevel==j] <- tmean
        outDF$Tmin[outDF$Label==i&outDF$CurrentLevel==j] <- tmin
        outDF$Tmax[outDF$Label==i&outDF$CurrentLevel==j] <- tmax
        
        ### linear fit
        fit1 <- lm(CO2~Record2,data=subDF2)
        
        ### assign fit coefficients
        outDF$linear_slope[outDF$Label==i&outDF$CurrentLevel==j] <- coefficients(fit1)[2]
        outDF$linear_intercept[outDF$Label==i&outDF$CurrentLevel==j] <- coefficients(fit1)[1]
        
        ### log fit
        fit2 <- lm(CO2~log(Record2), data=subDF2)
        
        ### assign fit coefficients
        outDF$log_slope[outDF$Label==i&outDF$CurrentLevel==j] <- coefficients(fit2)[2]
        outDF$log_intercept[outDF$Label==i&outDF$CurrentLevel==j] <- coefficients(fit2)[1]
        
        ### make some plots
        with(subDF2, plot(CO2~Record2, xlab = "Record", ylab = expression(CO[2] * " conc.")))
        lines(subDF2$Record2,predict(fit1),type='l',col='blue')
        lines(subDF2$Record2,predict(fit2),col='red')
        legend("topleft", legend=c("linear fit", "log fit"),
               col=c("blue", "red"), lty=1, cex=1.2)
        title(paste0(i, "_Record_", j))
    }
    
}

dev.off()


### save file
write.csv(outDF, "output/Rroot_coefficients_summary.csv", row.names=F)
