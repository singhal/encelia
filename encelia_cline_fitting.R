library(hzar)
library(scales)

##############################
# read in data
##############################

dir = '/Users/sonal/encelia/cline_fitting/'
setwd(dir)
d = read.csv("encelia_pca_scores.csv", stringsAsFactors=F)
s = read.csv("encelia_soil_data_17Apr16.csv", stringsAsFactors=F)

s$distance2 = rep(NA, nrow(s))
s[s$habitat == "desert", "distance2"] = s[s$habitat == "desert", "distance"] * -1 + 100
s[s$habitat == "ecotone", "distance2"] = s[s$habitat == "ecotone", "distance"] + 100
s[s$habitat == "dune", "distance2"] = s[s$habitat == "dune", "distance"] + 100

##############################
# the functions
##############################

fitCline <- function(df, distance, param, min, max, to_fix = NULL, val = NULL) {
	dist = unique(df[,distance])
	names(dist) = dist

	cl <- hzar.doNormalData1DRaw(dist, df[,distance], df[,param])
	cl_model = hzar.makeCline1DNormal(cl, tails="none")
	cl_model <- hzar.model.addBoxReq(cl_model, min, max)
	
	if (!missing(to_fix)) {
		cl_model[['parameterTypes']][[to_fix]]$val <- val
		attr(cl_model[['parameterTypes']][[to_fix]], 'fixed') <- TRUE
		}
		
	cl_model_fitR = hzar.first.fitRequest.gC(cl_model, cl, verbose=TRUE)
	cl_model_fitR$mcmcParam$chainLength = 1e5
	cl_model_fitR$mcmcParam$burnin = 1e4
	cl_model_fit = hzar.doFit(cl_model_fitR)

	cl_model_data <- hzar.dataGroup.add(cl_model_fit);
	cl_model_data <- hzar.dataGroup.add(cl_model_data, hzar.chain.doSeq(hzar.next.fitRequest(cl_model_fit)));
	
	return(cl_model_data)
	}
	
plot_clines1 <- function(cl, color) {
	# plot the 2LL cline uncertainty, transformed for min max values of the cline
	min = cl$ML.cline$param.free$muL
	max = cl$ML.cline$param.free$muR
	fzCoor <- hzar.getCredParamRed(cl)$fzCline(xSeries)
	polygon(x = c(fzCoor$x, rev(fzCoor$x)), y = c(transform(fzCoor$yMin, min, max), transform(rev(fzCoor$yMax), min, max)), col = alpha(color, 0.3), border=NA)
	lines(x = xSeries, y = transform(cl$ML.cline$clineFunc(xSeries), min, max), col = color, lwd=2)
	points(cl$obsData$frame$dist, transform(cl$obsData$frame$mu, min, max), col= color, pch=16, cex=0.8)
	}

plot_clines2 <- function(cl, color) {
	# plot the cline & original points, transformed for min max values of the cline
	min = cl$ML.cline$param.free$muL
	max = cl$ML.cline$param.free$muR
	lines(x = xSeries, y = transform(cl$ML.cline$clineFunc(xSeries), min, max), col = color, lwd=2)
	points(cl$obsData$frame$dist, transform(cl$obsData$frame$mu, min, max), col= color, pch=16, cex=0.8)
	}
	
transform <- function(vector, min, max) {
	# transform values to be standardized from 0 to 1
	return( (vector - min) / (max - min) )
}

testing_clines <- function(cl1, cl2,  cl1_con, cl2_con, type) {
	# test for concordance and coincidence
	# by looking at AIC values
	aic1 = hzar.AIC.hzar.dataGroup(cl1)
	aic1_con = hzar.AIC.hzar.dataGroup(cl1_con)
	aic2 = hzar.AIC.hzar.dataGroup(cl2)
	aic2_con = hzar.AIC.hzar.dataGroup(cl2_con)
	print(paste("AIC for original cline model 1: ", aic1, sep=""))
	print(paste("AIC for ", type, " cline model 1: ", aic1_con, sep=""))
	aics = c(aic1, aic1_con)
	diff = exp((min(aics) - max(aics))/2)
	print(paste("***The constrained model explains ", diff, " of the information.", sep=""))
	print(paste("AIC for original cline model 2: ", aic2, sep=""))
	print(paste("AIC for ", type, " cline model 2: ", aic2_con, sep=""))
	aics = c(aic2, aic2_con)
	diff = exp((min(aics) - max(aics))/2)
	print(paste("***The constrained model explains ", diff, " of the information.", sep=""))
	}
	
##############################
# run it all!
##############################

cl_morph = fitCline(d, 'distance2', 'PC1', -20, 220)
cl_soil = fitCline(s, 'distance2', 'pctmedsand', -20, 220)

# set up the plot
plot(NULL, xlim=c(-5, 205), ylim=c(-0.2, 1.2), xlab="distance (m)", ylab="standarized trait value", yaxt="n")
axis(side=2, at=c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), labels=T)
# get the range for x values
xSeries <- seq(from = par("usr")[1], to = par("usr")[2], length.out = 109)
plot_clines1(cl_morph, "#d8b365")
plot_clines1(cl_soil, "#5ab4ac")
plot_clines2(cl_morph, "#d8b365")
plot_clines2(cl_soil, "#5ab4ac")
legend("topleft", c("morphology, PC1", "% fine sand"), col=c("#d8b365", "#5ab4ac"), bty="n", lwd=2)

# test for coincidence (same center)
center_soil = as.numeric(cl_soil$ML.cline$param.free['center'])
center_morph = as.numeric(cl_morph$ML.cline$param.free['center'])
cl_morph_coi = fitCline(d, 'distance2', 'PC1', -20, 220, "center", center_soil)
cl_soil_coi = fitCline(s, 'distance2', 'pctmedsand', -20, 220, "center", center_morph) 
testing_clines(cl_morph, cl_soil, cl_morph_coi, cl_soil_coi, "coincidence")

# test for concordance (same width)
width_soil = as.numeric(cl_soil$ML.cline$param.free['width'])
width_morph = as.numeric(cl_morph$ML.cline$param.free['width'])
cl_morph_con = fitCline(d, 'distance2', 'PC1', -20, 220, "width", width_soil)
cl_soil_con = fitCline(s, 'distance2', 'pctmedsand', -20, 220, "width", width_morph)
testing_clines(cl_morph, cl_soil, cl_morph_con, cl_soil_con, "concordance")
