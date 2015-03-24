
library(ggplot2)
library(ks)
library(gplots)
library(GenomicRanges)
require(ggbio)


setwd("~/Desktop/Dist_matrix_TE_div/")



reuben.biplot <- function (x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2), 
    xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL, 
    arrow.len = 0.1, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, x.col = 1,
    ...) 
{
    n <- nrow(x)
    p <- nrow(y)
    if (missing(xlabs)) {
        xlabs <- dimnames(x)[[1L]]
        if (is.null(xlabs)) 
            xlabs <- 1L:n
    }
    xlabs <- as.character(xlabs)
    dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
    if (missing(ylabs)) {
        ylabs <- dimnames(y)[[1L]]
        if (is.null(ylabs)) 
            ylabs <- paste("Var", 1L:p)
    }
    ylabs <- as.character(ylabs)
    dimnames(y) <- list(ylabs, dimnames(y)[[2L]])
    if (length(cex) == 1L) 
        cex <- c(cex, cex)
    if (missing(col)) {
        col <- par("col")
        if (!is.numeric(col)) 
            col <- match(col, palette(), nomatch = 1L)
        col <- c(col, col + 1L)
    }
    else if (length(col) == 1L) 
        col <- c(col, col)
    unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)), 
        abs(max(x, na.rm = TRUE)))
    rangx1 <- unsigned.range(x[, 1L])
    rangx2 <- unsigned.range(x[, 2L])
    rangy1 <- unsigned.range(y[, 1L])
    rangy2 <- unsigned.range(y[, 2L])
    if (missing(xlim) && missing(ylim)) 
        xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
    else if (missing(xlim)) 
        xlim <- rangx1
    else if (missing(ylim)) 
        ylim <- rangx2
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    on.exit(par(op))
    op <- par(pty = "s")
    if (!is.null(main)) 
        op <- c(op, par(mar = par("mar") + c(0, 0, 1, 0)))
    plot(
    	x, 
    	type = "p", 
    	xlim = xlim, 
    	ylim = ylim, 
    	col = x.col, 
        xlab = xlab, 
        ylab = ylab, 
        sub = sub, 
        main = main,
        pch = 16,
        cex = .5, 
        ...)
#    text(
#    	x, 
#    	xlabs, 
#    	cex = cex[1L], 
#    	col = x.col, 
#    	...)
    par(new = TRUE)
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    plot(
    	y, 
    	axes = FALSE, 
    	type = "n", 
    	xlim = xlim * ratio, 
    	ylim = ylim * ratio, 
    	xlab = "", 
    	ylab = "", 
    	col = col[1L], 
    	...)
    axis(3, col = col[2L], ...)
    axis(4, col = col[2L], ...)
    box(col = col[1L])
    text(y, 
    	labels = ylabs, 
    	cex = cex[2L], 
    	col = col[2L], 
    	...)
    if (var.axes) 
        arrows(0, 0, y[, 1L] * 0.8, y[, 2L] * 0.8, col = col[2L], 
            length = arrow.len)
    invisible()
}



# how about remapping two species to human and comparing them to each other



# read in evert


spec1 <- "Human"
spec2 <- "Horse"
ucsc1 <- "hg19"
rem.un <- "yes"
no.xy <- F

keep.limit <- 1350000
mb <- 1000000
# bin.size is set so that the indexing stays consitent across pipeline
bin.size = 500000

# trying out the unscaled method
keep.NGenes = "yes"
keep.NG4s = "no"
keep.NCpGI = "no"
keep.CpGBP = "no"
keep.GC = "no"
SCALE = "no"

#create objects into whcih i will store the binsizes
s1name <- paste("count_tables/",spec1, "_AllBinCounts.txt", sep = "")
s2name <- paste("count_tables/",spec2, "_AllBinCounts.txt", sep = "")
s1 <- read.table(s1name, header = TRUE)
s2 <- read.table(s2name, header = TRUE)
slist <- list(s1,s2)


for(i in seq(along=slist)){
	count <- slist[[i]]
	count <- count[count$Known >= bin.size,]
	count$GC <- count$GC/count$Known*100
    count[,5:(length(count)-2)] <- ((count[,5:(length(count)-2)]/count$Known) * mb) 
    if(keep.NGenes == "no"){count <- count[,!(colnames(count) == "NGenes")]}
    if(keep.NG4s ==	"no"){count <- count[,!(colnames(count) == "NG4s")]}
    if(keep.NCpGI ==	"no"){count <- count[,!(colnames(count) == "NCpGI")]}
    if(keep.CpGBP  ==	"no"){count <- count[,!(colnames(count) == "CpGBP")]}
    if(keep.GC ==	"no"){count <- count[,!(colnames(count) == "GC")]}
	if(SCALE == "yes"){count[,5:(length(count)-1)] <- scale(count[,5:(length(count)-1)])}

    #count <- count[,!(colnames(count) == "Known")]
    colnames(count)[1:4] <- c("chr", "binID", "start", "end")
    count$binID <- 1:dim(count)[1]
      slist[[i]] <- count
}

KnownS1 <- data.frame(slist[[1]]$binID, slist[[1]]$Known)
KnownS2	<- data.frame(slist[[2]]$binID,	slist[[2]]$Known)

s1 <- slist[[1]][,!(colnames(slist[[1]]) == "Known")]
s2 <- slist[[2]][,!(colnames(slist[[2]]) == "Known")]

if(rem.un == "yes"){
	if(length(grep("U", s1$chr)) > 0){s1 <- s1[-(grep( "U", s1$chr)),]}
	if(length(grep("_", s1$chr)) > 0){s1 <- s1[-(grep("_", s1$chr)),]}
	if(length(grep("M", s1$chr)) > 0){s1 <- s1[-(grep("M", s1$chr)),]}
	}
if(rem.un == "yes"){
	if(length(grep("U", s2$chr)) > 0){s2 <- s2[-(grep( "U", s2$chr)),]}
	if(length(grep("_", s2$chr)) > 0){s2 <- s2[-(grep("_", s2$chr)),]}
	if(length(grep("M", s2$chr)) > 0){s2 <- s2[-(grep("M", s2$chr)),]}
	}
	

human <- s1
bovine <- s2


keep.human <- KnownS1[KnownS1[,2] > keep.limit,1]
keep.bovine <- KnownS2[KnownS2[,2] > keep.limit,1]

human <- human[human$binID %in% keep.human,]
bovine <- bovine[bovine$binID %in% keep.bovine,]

chr_main = "whole_genome"
if(no.xy){
	human <- human[!(human$chr == "chrX"),]
	human <- human[!(human$chr == "chrY"),]
	bovine <- bovine[!(bovine$chr == "chrX"),]
	bovine <- bovine[!(bovine$chr == "chrY"),]
	chr_main = "autosome_only"
}




# we will try the change
# remeber to sort the colums correctly so the whole thing works
b.sbin <- read.table(paste("S_bin/Human_sbin/",spec1,"_aligning_",spec2, "_select", sep = ""), header = TRUE)

#b.sbin <- b.sbin[,c(3,4,1,2)]
#colnames(b.sbin) <- c("S1.bin", "S1.Proportion", "S2.bin", "S2.Proportion")




# after weve made the decisions we can set the species and the cutoff here


cutoff <- 0.2







	

S.bin_s1_s2 <- b.sbin
s1 <- human
s2 <- bovine


S.bin_s1_s2 <- S.bin_s1_s2[S.bin_s1_s2[,4] > cutoff,]
# make sure each bin has the df has the same bins
S.bin_s1_s2 <- S.bin_s1_s2[S.bin_s1_s2$S1.bin %in% s1$binID ,]
S.bin_s1_s2 <- S.bin_s1_s2[S.bin_s1_s2$S2.bin %in% s2$binID ,]
s1 <- s1[s1$binID %in% S.bin_s1_s2$S1.bin,]



s2TE <- s2[,c(2, 5:length(s2))]
merg.s1_s2 <- merge(S.bin_s1_s2, s2TE, by.x=3, by.y=1)	
	
S2_TE_name <- colnames(s2TE)[2:length(s2TE)] 
for( i in S2_TE_name){
		merg.s1_s2[,i] <- merg.s1_s2[,i]*merg.s1_s2$S2.Proportion
	}
# now we sum them
# then we sum the S1.proportion and divide em
newS2.coord <- matrix(data=rep(0, (length(s2TE)-1)*length(unique(merg.s1_s2$S1.bin)) ), ncol = (length(s2TE)-1), nrow =	length(unique(merg.s1_s2$S1.bin)) )
colnames(newS2.coord) <- S2_TE_name
replot_s1_s2 <- data.frame(S1_bin = unique(merg.s1_s2$S1.bin),newS2.coord)
for(i in seq(dim(replot_s1_s2)[1])){	
	divide <- sum(merg.s1_s2[replot_s1_s2[i,'S1_bin'] == merg.s1_s2[,'S1.bin'],'S1.Proportion'])
	replot_s1_s2[i,2:length(replot_s1_s2)] <- colSums(merg.s1_s2[replot_s1_s2[i,'S1_bin'] == merg.s1_s2[,'S1.bin'],	5:length(merg.s1_s2)]) / divide
}
replot_s1_s2 <- replot_s1_s2[order(replot_s1_s2[,'S1_bin']),]



# the guts of this need a rewrite because there are things occuring that don't resemble what we saw before 

all(replot_s1_s2$S1_bin == s1$binID)

rownames(replot_s1_s2) <- 1:dim(replot_s1_s2)[1]
rownames(s1) <- 1:dim(s1)[1]

	
	
TE.replot <- (replot_s1_s2[,2:length(replot_s1_s2)])
TE.s2 <- (s2TE[s2TE$binID %in% unique(S.bin_s1_s2$S2.bin),2:length(s2TE)])
#TE.s2 <- s2TE[,2:(length(s2TE)-1)]



plot(scale(s1$LINE_L1+s1$SINE1_7SL), scale(s1$SINE2_MIR))
points(scale(replot_s1_s2$LINE_L1+replot_s1_s2$SINE2_CanSINE), scale(replot_s1_s2$SINE2_MIR), col = 2)




#### assign a new name here to the s2


if(spec2 == "Dog"){
	assign("Dog.s2", replot_s1_s2)
}


if(spec2 == "Bovine"){
	assign("Dog.s2", replot_s1_s2)
}

if(spec2 == "Horse"){
	assign("Horse.s2", replot_s1_s2)
}







Horse.s2.keep <- Horse.s2[Horse.s2$S1_bin %in% Dog.s2$S1_bin,]
Dog.s2.keep <- Dog.s2[Dog.s2$S1_bin %in% Horse.s2.keep$S1_bin,]

all(Dog.s2.keep$S1_bin == Horse.s2.keep$S1_bin)


# some sort of distance function to see where TE contnet has changed
# so we probably have to scale becasue the counts have changed while density has remained the same
# what we are comparing is whather a bin is still the mean for a given TE distribution 




	s1.comp <- data.frame(
		new_LINE = Dog.s2.keep$LINE_L1 + Dog.s2.keep$LINE_RTE_BovB,  
		new_SINE = Dog.s2.keep$SINE2_BOV.tA +  Dog.s2.keep$SINE_BOVA2, 
		MIR = Dog.s2.keep$SINE2_MIR,
		L2 = Dog.s2.keep$LINE_L2
	)


	s2.comp <- data.frame(
		new_LINE = Horse.s2.keep$LINE_L1 , 
		new_SINE = Horse.s2.keep$SINE2_ERE1 + Horse.s2.keep$SINE2_ERE2 + Horse.s2.keep$SINE2_ERE3 + Horse.s2.keep$SINE2_ERE4, 
		MIR = Horse.s2.keep$SINE2_MIR,
		L2 = Horse.s2.keep$LINE_L2
	)







#s1.comp <- data.frame(scale(s1.comp, center = FALSE))
#s2.comp <- data.frame(scale(s2.comp, center = FALSE))



pca.s1 <- prcomp(s1.comp, scale. = F)
pca.s2 <- prcomp(s2.comp, scale. = F)

reuben.biplot(pca.s1$x, pca.s1$rotation, x.col = "grey")
reuben.biplot(pca.s2$x, pca.s2$rotation, x.col = "grey")



all_2_all = F

if(all_2_all){
	A.dist = as.matrix(dist(s1.comp))
	B.dist = as.matrix(dist(s2.comp))
	d.te = colMeans(abs(A.dist - B.dist))
}else{
	d.te <- rep(NA, dim(s1.comp)[1])
	for(i in 1:dim(s1.comp)[1]){
		d.te[i] = dist(rbind(s2.comp[i,], s1.comp[i,]))
	}
}





d.te <- log(d.te)

dense.dist <- density((d.te))
plot(dense.dist)



qplot(s1.comp$new_SINE, s2.comp$new_SINE, col = d.te) + scale_color_gradient(low = "yellow", high = "blue") 
qplot(s1.comp$new_LINE, s2.comp$new_LINE, col = d.te) + scale_color_gradient(low = "yellow", high = "blue") + stat_smooth(method = "glm", level= 0.9)
qplot(s1.comp$MIR, s2.comp$MIR, col = d.te) + scale_color_gradient(low = "yellow", high = "blue") + stat_smooth(method = "glm", level= 0.9)
qplot(s1.comp$L2, s2.comp$L2, col = d.te) + scale_color_gradient(low = "yellow", high = "blue")

