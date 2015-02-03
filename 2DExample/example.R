# EXAMPLE.R
#
# Example of a proteomic data analysis workflow in R
# (everything written after a "hash sign" is either a comment or "deactivated code")



# PART 1: SETUP

# Load required R packages

require(caret)
require(klaR)
require(e1071)
require(effsize)
require(pcaPP)
require(mixOmics)
require(fastICA)
require(pheatmap)
require(vsn)
require(qvalue)
require(limma)

# Import example dataset and corresponding metadata

# (Note: at least one of these three commands should work)
input.data <- read.delim("http://files.figshare.com/1860747/example_dataset.txt",dec=".")
#input.data <- read.delim("http://tomessilva.github.io/2DExample/example_dataset.txt",dec=".")
#input.data <- read.delim("http://tomesilva.com/2DExample/example_dataset.txt",dec=".")

head(input.data[,1:6]) # check if data was correctly imported by looking at the first few rows/columns

# (Note: at least one of these three commands should work)
input.metadata <- read.delim("http://files.figshare.com/1860748/example_sample_metadata.txt",dec=".")
#input.metadata <- read.delim("http://tomessilva.github.io/2DExample/example_sample_metadata.txt",dec=".")
#input.metadata <- read.delim("http://tomesilva.com/2DExample/example_sample_metadata.txt",dec=".")

print(input.metadata) # check if metadata was correctly imported
treatment <- input.metadata$treatment



# PART 2: DATA NORMALIZATION AND TRANSFORMATION

boxplot(t(input.data)) # show distribution of spot abundances for each gel before normalization
meanSdPlot(t(input.data)) # show how standard deviation of each spot depends on its mean value before normalization

norm.data <- t(justvsn(t(as.matrix(input.data)))) # normalize/transform dataset using VSN method

boxplot(t(norm.data)) # show distribution of spot abundances for each gel after normalization
meanSdPlot(t(norm.data)) # show how standard deviation of each spot depends on its mean value after normalization



# PART 3: DATA VISUALIZATION

# Principal Component Analysis
results.pca <- prcomp(norm.data,scale=TRUE,retx=TRUE) # perform Principal Component Analysis with mean-centering and scaling
plot(results.pca) # scree plot
biplot(results.pca) # biplot of scores and loadings (not very clear)
plot(results.pca$x[,1],results.pca$x[,2],col=1+treatment,xlab="PC1",ylab="PC2",pch=15,cex=2) # alternative scores plot
legend("top",c("class 0","class 1"),col=1:2,pch=15,cex=2) # add a legend to the plot

# Multidimensional Scaling
inter.gel.distances <- as.dist(1 - cor.fk(t(norm.data))^2) # dissimilarities/distances based on Kendall's tau correlation coefficient
results.mds <- cmdscale(inter.gel.distances) # perform metric Multidimensional Scaling using the calculated distances
plot(results.mds,col=1+treatment,xlab="Dimension 1",ylab="Dimension 2",pch=15) # plot samples/gels according to the 2D embedding generated in the previous step
legend("top",c("class 0","class 1"),col=1:2,pch=15) # add a legend to the plot

# Independent Component Analysis
results.ica <- fastICA(norm.data,2) # perform Independent Component Analysis using the FastICA algorithm
plot(results.ica$S,col=1+treatment,xlab="IC1",ylab="IC2",pch=15) # scores plot
legend("top",c("class 0","class 1"),col=1:2,pch=15) # add a legend to the plot

# Heatmap
pheatmap(t(norm.data),scale="row",show_rownames=FALSE) # plot heatmap with mean-centering and scaling by row



# PART 4: DIFFERENTIAL ANALYSIS

# Hypothesis testing (Smyth's moderated t-test) with multiple comparison correction (Storey's q-value method)
design.mtrx <- cbind(rep(1,length(treatment)),treatment) # define design matrix (Intercept + treatment)
fit <- lmFit(t(norm.data),design=design.mtrx,method="robust",maxit=1024) # apply robust linear fit to each spot
fit <- eBayes(fit,robust=TRUE) # compute moderated t-statistics
qval <- (qvalue(fit$p.value[,2],pi0.method="bootstrap",gui=FALSE))$qvalues # calculate q-values
fx.size <- apply(norm.data,2,function(d,f) cohen.d(d,f)$estimate,f=factor(treatment)) # calculate Cohen's d for each spot
sig.spots <- names(qval[qval < 0.1]) # assuming FDR < 10%
n.sig.spots <- length(sig.spots) # number of significant spots
spot.class <- as.numeric(colnames(norm.data) %in% sig.spots) # spot significance (binary variable)

# Plot volcano plot
spot.colours <- 1 + spot.class + as.numeric(colnames(norm.data) %in% names(qval[qval < 0.05])) # calculate colors for each spot class
plot(fx.size,-log(qval)/log(10),col=spot.colours,xlab="Effect size (Cohen's d)",ylab="Significance (-log10(q-value))",pch=15,cex=1) # volcano plot
legend("bottomleft",c("FDR > 10%","FDR < 10%","FDR < 5%"),col=1:3,pch=15) # add a legend to the plot

# Export results of the univariate analysis
univariate.results <- data.frame(spot.name=colnames(norm.data),p.value=fit$p.value[,2],
                                 q.value=qval,effect.size=fx.size,significant=spot.class)
write.table(univariate.results,"univariate_results.txt",row.names=FALSE,sep="\t") # output results to a text file, using TAB as field separator

# Sparse Partial Least Squares Discriminant Analysis
results.splsda <- mixOmics::splsda(norm.data, as.factor(treatment), ncomp=2, keepX=c(n.sig.spots,n.sig.spots)) # perform sPLS-DA keeping only "n.sig.spots" per component
sig.spots.splsda <- names((((results.splsda$loadings)$X)[,1])[abs(((results.splsda$loadings)$X)[,1]) > (.Machine$double.eps)]) # check which spots have zero weight for the first component

# Recursive Feature Elimination (with Naive Bayes classifiers)
rfeCtrl <- rfeControl(functions = nbFuncs,method = "LOOCV") # define RFE parameters (use Naive Bayes, use leave-one-out cross-validation)
results.rfe <- rfe(norm.data,factor(treatment),sizes=n.sig.spots,rfeControl=rfeCtrl) # perform RFE (slow)
sig.spots.nbrfe <- rownames(varImp(results.rfe))[(1:n.sig.spots)] # choose the "n.sig.spots" most relevant spots

# Comparison of the different feature selection approaches
all.spots <- colnames(norm.data)
tmp.counts <- matrix(0, nrow=length(all.spots), ncol=3)
for (i in 1:length(all.spots)) { # check which spots have been selected by each method
  tmp.counts[i,1] <- all.spots[i] %in% sig.spots
  tmp.counts[i,2] <- all.spots[i] %in% sig.spots.splsda
  tmp.counts[i,3] <- all.spots[i] %in% sig.spots.nbrfe
}
colnames(tmp.counts) <- c("moderated t-test","sPLS-DA","RFE (Naive Bayes)") # define labels
vennDiagram(vennCounts(tmp.counts)) # plot Venn diagram

# Heatmap with the consensus set
consensus.set <- intersect(sig.spots,intersect(sig.spots.splsda,sig.spots.nbrfe)) # calculate consensus set
pheatmap(t(norm.data[,consensus.set]),scale="row") # plot heatmap only with the consensus set
