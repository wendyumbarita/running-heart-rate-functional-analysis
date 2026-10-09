# Functional Data Analysis of Running Heart Rates
# Based on the original coursework implementation.
# NOTE: K-means cluster labels are arbitrary and may change between runs.

#----------------------------Different type of runs----------------------------
#------------------Easy runs, Interval runs, and race runs---------------------
require(fda)
require(refund)
library(ggplot2)

#This one is to remove background. 
#par(bg=NA) #comment this if showing plots during the presentation. 

data_path <- file.path("data", "Diff_runs.txt")
if (!file.exists(data_path)) {
  stop("Missing data/Diff_runs.txt. See data/README.md for the expected file format.")
}
data <- read.delim(data_path, header = TRUE, stringsAsFactors = FALSE)
run.names <- names(data)[-1]
group <- factor(c("Aerobic","Aerobic","Aerobic","Aerobic","Aerobic","Aerobic",
                  "Intervals","Intervals","Intervals","Intervals","Intervals","Race","Race",
                  "Race")) #Adding the category for each of the runs. This will be used
#later for the function on scalar model, this will be the scalar predictor

names(group) <- run.names


lap.index <- seq(0.4,9.2,0.4)

hr.matrix <- as.matrix(data[, -1])
#Create a cubic b splines with 10 basis functions
#penalty of second derivative, smooth the data so that now we can represent
#each run as a curve
dat.basis <- create.bspline.basis(rangeval = c(0.4,9.2), nbasis = 10)
basis.sm <- smooth.basisPar(argvals = lap.index, hr.matrix, dat.basis)
hr.fd <- basis.sm$fd


#-------------------Function on scalar----------------------------

#Fit the model using the Generalized Least Squares estimator
mod.mat <- cbind(1, model.matrix(~group - 1)) #This is the matrix containing the
#intercept and the indicator variable. 1 or 0 depending on the type of run
constraints=matrix(c(0,1,1,1),1)
gls.mod=fosr(fdobj=hr.fd,X=mod.mat, con=constraints, method="GLS")
par(mfrow=c(1,4))
plot(gls.mod, split=1, set.mfrow=FALSE, titles=c("GLS: Intercept", levels(factor(group))), ylab="", xlab="Kilometer")




#-----------------------------------------------------------------
par(mfrow=c(1,3))

plot(lap.index, data$Dec17, type = "l", ylim = c(100, 200),
     xlab = "Distance (km)", ylab = "Heart rate (bpm)", main = "Example Aerobic Run")

plot(lap.index, data$Oct_25, type = "l", ylim = c(100, 200),
     xlab = "Distance (km)", ylab = "Heart rate (bpm)", main = "Example Interval Run")

plot(lap.index, data$Jun_Race, type = "l",  ylim = c(100, 200),
     xlab = "Distance (km)", ylab = "Heart rate (bpm)", main = "Example Race Run")



#-------------------------Perform PCA-----------------------------
#Create a cubic b splines with 25 basis. Add second derivative penalty
dat.basis <- create.bspline.basis(rangeval = c(0.4,9.2), nbasis = 25)
harmfdpar <- fdPar(dat.basis, int2Lfd(2), lambda = 1e-2)
numharms = 4
hrpcacobj <- pca.fd(hr.fd, nharm = numharms, harmfdpar)
hrpcacobjVM <- varmx.pca.fd(hrpcacobj)


par(mfrow=c(2,2), pty="m")
plot.pca.fd(hrpcacobjVM)

scores <- hrpcacobjVM$scores

par(mfrow=c(1,2))
#Plot each component one against the other one 

plot(scores[,1], scores[,2], pch=16,
     xlab="PC1 score", ylab="PC2 score",
     main="Runs in PC1–PC2", xlim  = c(-40, 45))
text(scores[,1], scores[,2], labels = run.names,
     pos = 4, cex = 0.7)

plot(scores[,1], scores[,3], pch=16,
     xlab="PC1 score", ylab="PC3 score",
     main="Runs in PC1–PC3", xlim  = c(-40, 45))
text(scores[,1], scores[,3], labels = run.names,
     pos = 4, cex = 0.7)

plot(scores[,1], scores[,4], pch=16,
     xlab="PC1 score", ylab="PC4 score",
     main="Runs in PC1–PC4", xlim  = c(-40, 45))
text(scores[,1], scores[,4], labels = run.names,
     pos = 4, cex = 0.7)

plot(scores[,2], scores[,3], pch=16,
     xlab="PC2 score", ylab="PC3 score",
     main="Runs in PC2–PC3", xlim  = c(-15, 32))
text(scores[,2], scores[,3], labels = run.names,
     pos = 4, cex = 0.7)

plot(scores[,2], scores[,4], pch=16,
     xlab="PC2 score", ylab="PC4 score",
     main="Runs in PC2–PC4", xlim  = c(-15, 32))
text(scores[,2], scores[,4], labels = run.names,
     pos = 4, cex = 0.7)

plot(scores[,3], scores[,4], pch=16,
     xlab="PC3 score", ylab="PC4 score",
     main="Runs in PC3–PC4", xlim  = c(-25, 25))
text(scores[,3], scores[,4], labels = run.names,
     pos = 2, cex = 0.8)


#Let's do kmeans on the scores. 
kmeans.res <- kmeans(scores[, c(1,2,3,4)], centers = 3, nstart = 25)
cluster.labels <- factor(kmeans.res$cluster, levels=c(1,2,3))#, labels = c("Interval Runs","Aerobic Runs","Race Runs"))
scores.df <- as.data.frame(scores)
colnames(scores.df) <- c("PC1", "PC2", "PC3", "PC4")
scores.df$Cluster <- cluster.labels
scores.df$Names <- run.names

p <- ggplot(scores.df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "PCA Score Plot with K-means Clustering (k=3)",
       x = "PC1", y = "PC2") +
  theme_minimal()

#Let's get back the representing function for each the groups. 
#Take the mean of the cluster in four dim,

c.k <- kmeans.res$centers 
f.bar <- hrpcacobjVM$meanfd 
Psi <- hrpcacobjVM$harmonics[1:4] 
f.k <- vector("list", nrow(c.k))
dist.grid <- seq(0.4, 9.2, length = 200)

for(k in 1:nrow(c.k)){
  temp.sum <- f.bar
  for(l in 1:4){
    temp.sum <- temp.sum + c.k[k,l] * Psi[l]
  }
  f.k[[k]] <- temp.sum
}

#This can change
func.group1 <- eval.fd(dist.grid, f.k[[3]])
func.group2 <- eval.fd(dist.grid, f.k[[2]])
func.group3 <- eval.fd(dist.grid, f.k[[1]])

par(mfrow=c(1,3))

plot(dist.grid, func.group1, type = "l", xlab = "Kilometer", ylab = "Heart Rate", main = "Function for group 1 (Aerobic)",
     ylim = c(120, 190), col = "black")
plot(dist.grid, func.group2, type = "l", xlab = "Kilometer", ylab = "Heart Rate", main = "Function for group 2 (Intervals)",
     ylim = c(120, 190), col = "black")
plot(dist.grid, func.group3, type = "l", xlab = "Kilometer", ylab = "Heart Rate", main = "Function for group 3 (Race)",
     ylim = c(120, 190), col = "black")

data.intervals <- data[,8:12]#Interval runs
matplot(lap.index, data.intervals, type = "l", xlab = "Distance km", ylab = "Heart rate")

