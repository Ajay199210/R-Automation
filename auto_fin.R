
# blueprint.R: A script that automates financial statistical tasks.

##########################
#-- Function definitions #
##########################

#-- Prepare environment
prep.env <- function(req.lib,dir.loc) {
  # Load required libraries
  lapply(req.lib, require, character.only = TRUE)
  # Set working directory location
  setwd(dir.loc)
  # Get all CSV files in the current directory
  print(paste("Listing all CSV files in", getwd()))
  csv.files <- dir()[grep("*.csv", dir())]
  csv.files
}

#-- Import data into zoo format 
import <- function(data.file){
  read.zoo(data.file, header=TRUE, sep=",", format="%m/%d/%Y", index.column=1)
}

#-- Clean the data from unwanted variables (columns) IF NECESSARY
clean.data <- function(mydata,cols) {
  col.ind <- list()
  counter <- 1
  for (i in 1:length(cols)) {
    col.ind[[counter]] <- grep(c(cols)[i], colnames(mydata))
    counter <- counter + 1
  }
  # print(col.ind)
  counter <- 1
  subset.vec <- c()
  for (j in col.ind) {
    subset.vec[counter] <- length(j)
    counter <- counter + 1
  }
  subset.vec <- rep(NA,sum(subset.vec))
  # print(subset.vec)
  counter <- 1
  for (k in col.ind) {
    for (l in k) {
      subset.vec[counter] <- l
      counter <- counter + 1
    }
  }
  
  # print(subset.vec)
  return(mydata[,-subset.vec])
}

#-- Transform the columns to numeric variables IF NECESSARY
col.to.num <- function(mydata) {
  # Converting the columns to numeric and storing them in a list
  num.list <- list()
  for (i in 1:ncol(mydata)) {
    num.list[[i]] <- as.numeric(mydata[,i])
    names(num.list)[i] <- colnames(mydata)[i]
  }
  # print(num.list)
  
  # Creating again the zoo series in a matrix
  data.mat <- matrix(NA,nrow=nrow(mydata),ncol=ncol(mydata))
  ind <- 1
  for (j in num.list) {
    data.mat[,ind] <- j
    ind <- ind + 1
  }
  colnames(data.mat) <- colnames(mydata)
  return(zoo(data.mat, order.by = index(mydata)))
}

#-- Descriptive stats
des.info <- function(mydata, write=F) {
  m <- apply(mydata,2,mean)
  M. <- apply(mydata,2,max)
  m. <- apply(mydata,2,min)
  std.dev <- apply(mydata,2,sd)
  skw <- apply(mydata,2,skewness)
  data.des.info <- cbind(m,M.,m.,std.dev,skw)
  colnames(data.des.info) <- c("Mean", "Max", "Min", "Standard Deviation",
                               "skewness")
  if (write==T) {
    write.csv(data.des.info, "descriptiveInfo.csv", row.names=T)
  }
  
  return(data.des.info)
}

#-- ADF test to check if the data is stationary
aug.df.test <- function(mydata,write=F) {
  
  # Constructing the data frame for the ADF test data frame with row names
  adfTest.df <- as.data.frame(matrix(NA,nrow=ncol(mydata),ncol=4))
  colnames(adfTest.df) <- c("ADF Test Stat.", "P-Value", "Integration Order",
                            "Stationary")
  rownames(adfTest.df) <- colnames(mydata)
  
  # Creating & filling the 'ADF Test Statistic' & 'P-value' columns - col 1 & 2
  statistic <- c()
  p.values <- c()
  counter <- 1
  for (i in 1:ncol(mydata)) {
    statistic[counter] <- adf.test(mydata[,i])$statistic
    p.values[counter] <- adf.test(mydata[,i])$p.value
    
    # Filling the 'Stationary' column - col 4
    if (p.values[counter] < 0.05) {
      adfTest.df[,4][counter] <- TRUE
    } else {
      adfTest.df[,4][counter] <- FALSE
    }
    counter <- counter + 1
  }
  
  adfTest.df[,1] <- round(statistic,3) # col 1
  adfTest.df[,2] <- round(p.values,3) # col 2
  adfTest.df[,3] <- "I(0)" # col 3
  
  # print(statistic)
  # print(p.values)
  
  if (write==T) {
    write.csv(adfTest.df, "adfTest-Table.csv", row.names=T)
  }
  
  return(adfTest.df)
}

#-- Calculate the 1st differences for the logarithmic values of each variable
# (compounded returns)
diff.log <- function(mydata,adfTable) {
  if (all(adfTable$Stationary)) {
    return(mydata)
  } else if (any(adfTable$Stationary)) {
    stat.var.name <- rownames(adfTable[adfTable$Stationary==T,])
    non.stat.var.name <- rownames(adfTable[adfTable$Stationary==F,])
    log.transform <- diff(log(clean.data(mydata,stat.var.name)))
    # clean.data() is a nested function which was defined previously
    no.log.transform <- clean.data(mydata,non.stat.var.name) # stationary var.
    data.diff <- merge(log.transform,no.log.transform)
    return(data.diff)
    
  } else if (!any(adfTable$Stationary)) {
    return(diff(log(mydata)))
  }
}

#-- Granger test for Reverse Causality
grang.caus <- function(mydata,dep.var,write=F) {
  
  # Constructing the granger test table
  grang.df <-  as.data.frame(matrix(NA,nrow=2*(ncol(mydata)-1),ncol=7))
  colnames(grang.df) <- c("Variables","Obs","F-Stat.", "P-Value","Stat. Sig.",
                          "Decision", "Causality Type")
  
  # Determining the dependent & independent variables
  y <- which(colnames(mydata)==dep.var) # dep. var
  x <- which(colnames(mydata)!=dep.var) # indep. var
  
  # Data arranged columns with the dependent variable in the 1st col
  data.arr.cols <- mydata[, c(y,x)]
  
  # Filling row names with the corresponding variables - col 1
  arr.cols <- colnames(data.arr.cols)
  ind <- 2
  counter <- 1
  for (k in 1:nrow(grang.df)/2) {
    # print(paste(arr.cols[1],"-->", arr.cols[ind]))
    # print(paste(arr.cols[ind],"-->", arr.cols[1]))
    grang.df[counter,1] <- paste(arr.cols[1],"-->", arr.cols[ind])
    grang.df[counter+1,1] <- paste(arr.cols[ind],"-->", arr.cols[1])
    
    counter <- counter + 2
    ind <- ind + 1
    
    if (ind > (nrow(grang.df)/2)+1) {
      break
    }
  } 
  
  # Filling the 'Obs' (observations) column  - col 2
  grang.df[,2] <- nrow(mydata)
  
  # Filling the 'F-Stat.' & 'P-Value' columns - col 3 & 4
  val.list <- list(f.stat=0, p.values=0)
  ind <- 1
  for (i in 2:ncol(data.arr.cols)) {
    val.list[[1]][ind] <- grangertest(data.arr.cols[,1],
                                      data.arr.cols[,i])$`F`[2]
    val.list[[1]][ind+1] <- grangertest(data.arr.cols[,i],
                                        data.arr.cols[,1])$`F`[2]
    
    val.list[[2]][ind] <- grangertest(data.arr.cols[,1],
                                      data.arr.cols[,i])$`Pr(>F)`[2]
    val.list[[2]][ind+1] <- grangertest(data.arr.cols[,i],
                                        data.arr.cols[,1])$`Pr(>F)`[2]
    
    ind <- ind + 2
  }
  
  grang.df[,3] <- round(val.list$f.stat, 4)
  grang.df[,4] <- round(val.list$p.values, 4)
  
  # Filling the 'Stat. Sig.', 'Decision' & 'Causality Type' columns  - col 5,6
  # & 7 respectively
  ind <- 1
  
  for (j in grang.df[,4]) {
    
    # Filling the statistical significance - col 5
    if (j <= 1 && j > 0.1) {
      grang.df[ind,5] <- "No sig"
    } else if (j <= 0.1 && j > 0.05) {
      grang.df[ind,5] <- "."
    } else if (j <= 0.05 && j > 0.01) {
      grang.df[ind,5] <- "*" 
    } else if (j <= 0.01 && j > 0.001) {
      grang.df[ind,5] <- "**"
    } else if (j <= 0.001 && j >= 0) {
      grang.df[ind,5] <- "***"
    }
    
    # Filling the decision - col 6
    if (j < 0.05) {
      grang.df[ind,6] <- "Reject H0"
    } else {
      grang.df[ind,6] <- "DNR H0" # Do Not Reject H0
    }
    
    ind <- ind + 1
  }
  
  # Filling the causality - col 7
  if (any(grang.df[,6] == "DNR H0")) {
    grang.df[,7][which(grang.df[,6]=="DNR H0")] <- "No causality"
  } 
  if (any(grang.df[,6] == "Reject H0")) {
    caus.var <- grang.df[,1][which(grang.df[,6]=="Reject H0")]
    # print(caus.var)
    indep.var <- arr.cols[which(arr.cols!=dep.var)]
    # print(indep.var)
    ind <- 1
    caus.var.dir <- c() # stands for causality.variable.direction
    for (k in 1:length(indep.var)) {
      caus.var.dir[ind] <- length(grep(indep.var[ind],caus.var))
      ind <- ind + 1
    }
    # print(caus.var.dir)
    caus.var.dir <- caus.var.dir[which(caus.var.dir != 0)]
    # print(caus.var.dir)
    caus.var.bidir <- which(caus.var.dir==2) # causality.variable.bi-directional
    dummy.caus.var.dir <- c(caus.var.dir,rep(0,length(caus.var.bidir)))
    # print(dummy.caus.var.dir)
    ind <- 1
    for (m in 1:length(dummy.caus.var.dir)) {
      if (dummy.caus.var.dir[m] == 1) {
        grang.df[which(grang.df[,6]=="Reject H0"),][ind,7] <- "Uni-directional"
      } else if (dummy.caus.var.dir[m] == 2) {
        grang.df[which(grang.df[,6]=="Reject H0"),][ind,7] <- "Bi-directional"
        ind <- ind + 1
        grang.df[which(grang.df[,6]=="Reject H0"),][ind,7] <- "Bi-directional"
      }
      ind <- ind + 1
    } 
  }
  
  if (write==T) {
    write.csv(grang.df, "grangerTest-table.csv", row.names=T)
  }
  
  # print(val.list)
  return(grang.df)
} 

#-- Simple Linear Model (OLS) 
simp.lm.list <- function(mydata,dep.var,write=F) {
  
  # Determining the dependent & independent variables
  y <- which(colnames(mydata)==dep.var) # dep. var
  x <- which(colnames(mydata)!=dep.var) # indep. var
  
  # Data arranged columns with the dependent variable in the 1st col
  data.arr.cols <- mydata[,c(y,x)]
  
  # Storing the simple regression summaries in a list
  sumr.list <- list()
  coef.list <- list()
  ind <- 1
  for (i in 2:ncol(data.arr.cols)) {
    sumr.list[[ind]] <- summary(lm(data.arr.cols[,1] ~ data.arr.cols[,i]))
    coef.list[[ind]] <- as.data.frame(coef(sumr.list[[ind]]))
    rownames(coef.list[[ind]])[2] <- colnames(mydata[,c(x)])[ind]
    coef.list[[ind]]$Stat.Sig <- NA
    ind <- ind + 1
  }
  names(sumr.list) <- colnames(mydata[,c(x)])
  names(coef.list) <- colnames(mydata[,c(x)])
  
  # Filling the Statistical significance
  counter <- 1
  for (j in coef.list) {
    ind <- 1
    for (k in j[,4]) {
      # print(k)
      if (k <= 1 && k > 0.1) {
        coef.list[[counter]][ind,5] <- "No sig"
      } else if (k <= 0.1 && k > 0.05) {
        coef.list[[counter]][ind,5] <- "."
      } else if (k <= 0.05 && k > 0.01) {
        coef.list[[counter]][ind,5] <- "*"
      } else if (k <= 0.01 && k > 0.001) {
        coef.list[[counter]][ind,5] <- "**"
      } else if (k <= 0.001 && k >= 0) {
        coef.list[[counter]][ind,5] <- "***"
      }
      ind <- ind + 1
    }
    counter <- counter + 1
  }
  
  # Writing the SIMPLE linear results to the hard drive (optional)
  if (write==T) {
    ind <- 1
    for (l in coef.list) {
      file.name <- paste(names(coef.list)[ind], "-LMcoef",".csv", sep="")
      # print(file.name)
      write.csv(l, file.name, row.names=T)
      ind <- ind + 1
    }
  }
  return(coef.list) # return the simple lm results in a form of a list
}

###############################
#-- Functions definitions end #
###############################

#-- Prepare your environment (change it to your own folder)
prep.env(c("corrplot", "ggcorrplot","xts", "zoo", "tseries", "readxl",
              "moments", "vars"),"")

#-- Import the data
usd.ind <- import("US Dollar Index Historical Data.csv")
sp500 <- import("S&P 500 Historical Data.csv")
wti <- import("Crude Oil WTI Futures Historical Data.csv")
ftse <- import("FTSE 100 Historical Data.csv")
stox600 <- import("STOXX 600 Historical Data.csv")

#-- Merge the data into a matrix
merged <- merge(usd.ind,sp500,wti,ftse,stox600, all=F)
is.matrix(merged) # True

merged <- clean.data(merged, c("Open","High","Low","Change"))
head(merged)
View(merged)

# Data summary
summary(merged)

#-- Overview about the data
desInfo <- des.info(merged, write=F)
desInfo
 
#-- ADF test
adf.table <- aug.df.test(merged, write=F)
adf.table

# Compounded returns
merged.ret <- diff.log(merged, adf.table)
colnames(merged.ret)
colnames(merged.ret) <- paste(substr(colnames(merged.ret),7,15),
                              ".ret", sep="")
head(merged.ret)

#-- Granger reverse causality test
granger.table <- grang.caus(merged, "Price.wti", write=F)
granger.table
granger.table.ret <- grang.caus(merged.ret, "wti.ret", write=F)
granger.table.ret

#-- Graphics
# Overview of the data prices & returns
par(mfrow=c(2,1))
plot(merged$Price.wti, xlab="Year", ylab="Crude Oil WTI Price",
     cex.lab=1.15, cex.axis=1.25, main="WTI")
plot(merged.ret$wti.ret, xlab="Year", ylab="Crude Oil WTI Return",
     cex.lab=1.15, cex.axis=1.25, main="WTI")

par(mfrow=c(2,2))
plot(merged$Price.sp500, xlab="Year", ylab="S&P500 Price", col=2,
     cex.lab=1.15, cex.axis=1.25, main="S&P500")
plot(merged.ret$sp500.ret, xlab="Year", ylab="S&P500 Return", col=2,
     cex.lab=1.15, cex.axis=1.25, main="S&P500")
plot(merged$Price.usd.ind, xlab="Year", ylab="USD Index Price", col=3,
     cex.lab=1.15, cex.axis=1.25, main="USD Index")
plot(merged.ret$usd.ind.ret, xlab="Year", ylab="USD Return", col=3,
     cex.lab=1.15, cex.axis=1.25, main="USD Index")

plot(merged$Price.ftse, xlab="Year", ylab="FTSE Price", col=4,
     cex.lab=1.15, cex.axis=1.25, main="FTSE")
plot(merged.ret$ftse.ret, xlab="Year", ylab="FTSE Return", col=4,
     cex.lab=1.15, cex.axis=1.25, main="FTSE")
plot(merged$Price.stox600, xlab="Year", ylab="Stox600 Price", col=6,
     cex.lab=1.15, cex.axis=1.25, main="Stox600")
plot(merged.ret$stox600.ret, xlab="Year", ylab="Stox600 Return", col=6,
     cex.lab=1.15, cex.axis=1.25, main="Stox600")

# Correlation Matrices for prices & returns
ggcorrplot(cor(merged), hc.order = TRUE, type = "lower",
           outline.col = "white", lab = TRUE, lab_size = 4)

ggcorrplot(cor(na.omit(merged.ret)), hc.order = TRUE, type = "lower",
           outline.col = "white", lab = TRUE, lab_size = 4)

#-- Building statistical models

## Linear models (OLS)
# Simple linear regressions
sr.usd.ind <- lm(merged.ret$wti.ret ~ merged.ret$usd.ind.ret)
sumr.usd.ind <- summary(sr.usd.ind)

sr.sp500 <- lm(merged.ret$wti.ret ~ merged.ret$sp500.ret)
sumr.sp500 <- summary(sr.sp500)

sr.ftse <- lm(merged.ret$wti.ret ~ merged.ret$ftse.ret)
sumr.ftse <- summary(sr.ftse)

sr.stox600 <- lm(merged.ret$wti.ret ~ merged.ret$stox600.ret)
sumr.stox600 <- summary(sr.stox600)

sumr.usd.ind
sumr.sp500
sumr.ftse
sumr.stox600

# Multiple linear regressions
mul.lm1 <- lm(merged.ret$wti.ret ~ merged.ret$usd.ind.ret +
                merged.ret$sp500.ret)
sumr.lm1 <- summary(mul.lm1)

mul.lm2 <- lm(merged.ret$wti.ret ~ merged.ret$usd.ind.ret +
                merged.ret$sp500.ret + merged.ret$ftse.ret)
sumr.lm2 <- summary(mul.lm2)

mul.lm3 <- lm(merged.ret$wti.ret ~ merged.ret$usd.ind.ret +
                merged.ret$sp500.ret + merged.ret$ftse.ret +
                merged.ret$stox600.ret)
sumr.lm3 <- summary(mul.lm3)

sumr.lm1
sumr.lm2
sumr.lm3

#-- Exporting the regression results

# Simple linear models results table
lm.results.table <- simp.lm.list(merged.ret, "wti.ret", write=F)
lm.results.table

# Multiple linear models results
mult.reg <- ls()[grep("sumr.lm*",ls())]
mult.reg
mult.reg.list <- list(sumr.lm1, sumr.lm2, sumr.lm3)
names(mult.reg.list) <- mult.reg
mult.reg.list

count <- 1
for (i in mult.reg.list) {
  file.name <- paste(mult.reg[count], ".csv", sep="")
  # Un-comment the below line to save the results on the hard drive
   write.csv(i, file.name, row.names=T)
  count <- count + 1
}


