

library(ggplot2)
library(ks)
library(gplots)
library(GenomicRanges)
require(ggbio)


setwd("~/Desktop/mapping_stats_analysis/")


rm(list = ls())






# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#







multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

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



spec1 <- "Human"
spec2 <- "Dog"
ucsc1 <- "hg19"
rem.un <- "yes"
no.xy <- F

keep.limit <- 1350000
mb <- 1000000
# bin.size is set so that the indexing stays consitent across pipeline
bin.size = 500000

# trying out the unscaled method
keep.NGenes = "no"
keep.NG4s = "no"
keep.NCpGI = "no"
keep.CpGBP = "no"
keep.GC = "no"
SCALE = "no"

#create objects into whcih i will store the binsizes
s1name <- paste("raw_data/bin_cov/",spec1, "_bin_coverage", sep = "")
s2name <- paste("raw_data/bin_cov/",spec2, "_bin_coverage", sep = "")
s1 <- read.table(s1name, header = TRUE)
s2 <- read.table(s2name, header = TRUE)
slist <- list(s1,s2)


for(i in seq(along=slist)){
	count <- slist[[i]]
#	count <- count[count$Known >= bin.size,]
#	count$GC <- count$GC/count$Known*100
    count[,5:(length(count)-2)] <- ((count[,5:(length(count)-2)]/count$Known) * mb) 
#    if(keep.NGenes == "no"){count <- count[,!(colnames(count) == "NGenes")]}
#    if(keep.NG4s ==	"no"){count <- count[,!(colnames(count) == "NG4s")]}
#    if(keep.NCpGI ==	"no"){count <- count[,!(colnames(count) == "NCpGI")]}
#    if(keep.CpGBP  ==	"no"){count <- count[,!(colnames(count) == "CpGBP")]}
#    if(keep.GC ==	"no"){count <- count[,!(colnames(count) == "GC")]}
#	if(SCALE == "yes"){count[,5:(length(count)-1)] <- scale(count[,5:(length(count)-1)])}

    #count <- count[,!(colnames(count) == "Known")]
    colnames(count)[1:4] <- c("chr", "binID", "start", "end")
 #   count$binID <- 1:dim(count)[1]
      slist[[i]] <- count
}

KnownS1 <- data.frame(slist[[1]]$binID, slist[[1]]$Known)
KnownS2	<- data.frame(slist[[2]]$binID,	slist[[2]]$Known)

s1 <- slist[[1]][,!(colnames(slist[[1]]) == "Known")]
s2 <- slist[[2]][,!(colnames(slist[[2]]) == "Known")]

s1 <- s1[,!(colnames(s1) == "binID.1")]
s2 <- s2[,!(colnames(s2) == "binID.1")]




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
b.sbin <- read.table(paste("~/Desktop/Dist_matrix_TE_div/S_bin/Human_sbin/",spec1,"_aligning_",spec2, "_select", sep = ""), header = TRUE)


cutoff <- 0.15







	

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






# some sort of distance function to see where TE contnet has changed
# so we probably have to scale becasue the counts have changed while density has remained the same
# what we are comparing is whather a bin is still the mean for a given TE distribution 

s1.comp <- data.frame(
		new_LINE = s1$new_LINE, 
		new_SINE = s1$new_SINE, 
		MIR = s1$MIR,
		L2 = s1$L2
)


	s2.comp <- data.frame(
		new_LINE = replot_s1_s2$new_LINE, 
		new_SINE = replot_s1_s2$new_SINE, 
		MIR = replot_s1_s2$MIR,
		L2 = replot_s1_s2$L2
	)



#s1.comp.s <- data.frame(scale(s1.comp, center = T))
#s2.comp.s <- data.frame(scale(s2.comp, center = T))


all_2_all = F

if(all_2_all){
	A.dist = as.matrix(dist(s1.comp), mehtod = "manhatten")
	B.dist = as.matrix(dist(s2.comp), mehtod = "manhatten")
	F.d.te = (colMeans(abs(A.dist - B.dist)))
}else{
	d.te <- rep(NA, dim(s1.comp)[1])
	for(i in 1:dim(s1.comp)[1]){
		d.te[i] = dist(rbind(s2.comp[i,], s1.comp[i,]), method = "manhattan")
	}
}

A.dist = mahalanobis(s1.comp, c(0,0,0,0), cov(s1.comp))
B.dist = mahalanobis(s2.comp, c(0,0,0,0), cov(s2.comp))
F.d.te = A.dist - B.dist

d.te.n <- d.te
d.te <- scale((d.te))

dense.dist <- density((d.te))
plot((dense.dist))


groups = rep("medium", length(d.te))
groups[d.te > 1.96] = "high"
groups[d.te < - 1.96] = "low"





t1 <- qplot(s1.comp$new_LINE, s2.comp$new_LINE, col = d.te.n) + scale_color_gradient(low = "lightblue", high = "red") + labs(title = colnames(s1.comp)[1]) + guides(col = FALSE)
t2 <- qplot(s1.comp$new_SINE, s2.comp$new_SINE, col = d.te.n) + scale_color_gradient(low = "lightblue", high = "red") + labs(title = colnames(s1.comp)[2]) + guides(col = FALSE)
t3 <- qplot(s1.comp$MIR, s2.comp$MIR, col = d.te.n) + scale_color_gradient(low = "lightblue", high = "red") + labs(title = colnames(s1.comp)[3]) + guides(col = FALSE)
t4 <- qplot(s1.comp$L2, s2.comp$L2, col = d.te.n) + scale_color_gradient(low = "lightblue", high = "red")+ labs(title = colnames(s1.comp)[4]) + guides(col = FALSE)

multiplot(t1,t2,t3,t4, cols = 2)


# lets try plotting distribtuions in 2d
i = 1


for(i in 1:4){
	dat <- data.frame(S1 = s1.comp[,i], S2 = s2.comp[,i], groups = groups)
	
	assign(
		paste("p", i , sep =""), 
		ggplot(dat, aes(x = S1, y = S2, col = groups)) +  geom_point() + labs(title = colnames(s1.comp)[i]) + guides(col = FALSE) + xlim(min(dat$S1), max(dat$S1)) + ylim(min(dat$S2), max(dat$S2))
	)
	
	assign(
		paste("q", i , sep =""),
		ggplot(dat, aes(x = S1, y = S2, col = groups,) ) + geom_density2d()  + labs(title = colnames(s1.comp)[i]) + guides(col = FALSE) + xlim(min(dat$S1), max(dat$S1)) + ylim(min(dat$S2), max(dat$S2))	
	)
		
	assign(
		paste("r", i , sep =""),
		ggplot(dat, aes(x = S1, y = S2, col = groups) ) + geom_density2d()  + labs(title = colnames(s1.comp)[i]) + guides(col = FALSE) + xlim(min(dat$S1), max(dat$S1)) + ylim(min(dat$S2), max(dat$S2))	 + geom_point()
	)

}


multiplot(p1, p2, p3, p4, cols = 2)

multiplot(q1, q2, q3, q4, cols = 2)

multiplot(r1, r2, r3, r4, cols = 2)





i = 1
dat <- data.frame(S1 = s1.comp[,i], S2 = s2.comp[,i], groups = groups)
pca.line <- prcomp(dat[,1:2])
dat <- data.frame(PC1 = pca.line$x[,1], PC2 = pca.line$x[,2], groups = groups)
p = ggplot(dat, aes( scale(PC1), fill = groups)) 
p = p +  geom_density(alpha = .7) 
p = p + labs(title = colnames(s1.comp)[i]) 
p




i = 2
dat <- data.frame(S1 = s1.comp[,i], S2 = s2.comp[,i], groups = groups)
pca.line <- prcomp(dat[,1:2], scale. = T)
dat <- data.frame(PC1 = pca.line$x[,1], PC2 = pca.line$x[,2], groups = groups)

q = ggplot(dat, aes(scale(PC1), scale(PC2), col = groups)) 
#q = q +  geom_density2d(alpha = .7) 
q = q + geom_point()
q = q + labs(title = colnames(s1.comp)[i]) 
q = q + geom_hline(yintercept= + 1.96)
q = q + geom_hline(yintercept= - 1.96)
q = q + geom_vline(xintercept= + 1.96)
q = q + geom_vline(xintercept= - 1.96)
q



multiplot(r2, q)


# for bins that a differnet we identify wheather there is an over or under abundance of a particular TE fmily 

# lets see if we can get the binIDS



# module comparisons


big.df <- read.table("~/Desktop/big.df.txt", header = T)

cols <- as.character(unique(big.df$mod.col))
pca <- prcomp(s1.comp)

for(c in 1:length(cols)){
	magenta <- big.df[big.df$mod.col == cols[c],]
	bins.col <- unique(magenta$binID)
	colors <- rep("grey", dim(pca$x)[1])
	colors[s1$binID %in% bins.col] <- "darkblue"
	reuben.biplot(pca$x, pca$rotation, x.col = colors, main = cols[c])
}




# so the plan here is to make a data frame of colors and bin IDS and a data frame of groups and binIDs \
# we then do a merge and see what happens

stats = data.frame(color = cols, odds_ratio = NA, p_val = NA)
for( c in 1:length(cols)){
	magenta <- big.df[big.df$mod.col == cols[c],]
	bins.col <- unique(magenta$binID)
	colors <- rep("grey", dim(pca$x)[1])
	colors[s1$binID %in% bins.col] <- "darkblue"
	colors_groups.df <- data.frame(groups = groups, colors = colors)
	con.tab <- table(colors_groups.df)
	f.test <- fisher.test(con.tab)
	stats[c, "odds_ratio"] <- f.test$estimate
	stats[c, "p_val"] <- f.test$p.value
}


# if we tighten our inteteresting classification we should expect a tightening of the odds ratio

barplot(stats[,"odds_ratio"], col = cols, names = stats$color, ylab = "odds ratio", ylim = c(0,5))

# maybe the thing isn't needed because both are made form the same table 
p_vals <- data.frame(pos = rep(0, length(cols)), sym = rep("", length(cols)), col = cols)
p_vals[stats$p_val < .001,"pos"] = 4.5
p_vals[stats$p_val < .001,"sym"] = "*"

text(p_vals$pos, p_vals$sym, col = cols)




# groups1 will be the high alu changes 







i = 2
dat <- data.frame(S1 = s1.comp[,i], S2 = s2.comp[,i], groups = groups)
pca.line <- prcomp(dat[,1:2])
dat <- data.frame(PC1 = pca.line$x[,1], PC2 = pca.line$x[,2], groups = groups)

groups2 <- rep("boring2", length(groups))
groups2[scale(dat$PC2) < -1.96  & groups == "high"] = "interesting2"
groups2[scale(dat$PC2) > 1.96  & groups == "high"] = "interesting2"







# this one needs some work
# before we can whole heartedly trust it !!





i = 2
dat <- data.frame(S1 = s1.comp[,i], S2 = s2.comp[,i], groups = groups)
pca.line <- prcomp(dat[,1:2])
dat <- data.frame(PC1 = pca.line$x[,1], PC2 = pca.line$x[,2], groups = groups)

groups1 <- rep("boring", length(groups))
groups1[scale(dat$PC2) < -1.96 & groups == "high"] = "interesting"


# so the interesting LINE1 stuff is not the same as the ALu stuff
# also we have only 4 bins that do not have interesting new line or sine content, maybe they have both ?



stats = data.frame(color = cols, odds_ratio = NA, p_val = NA)
for( c in 1:length(cols)){
	magenta <- big.df[big.df$mod.col == cols[c],]
	bins.col <- unique(magenta$binID)
	colors <- rep("grey", dim(pca$x)[1])
	colors[s1$binID %in% bins.col] <- "darkblue"
	colors_groups.df <- data.frame(groups = groups1, colors = colors)
	con.tab <- con.tab1 <- table(colors_groups.df)
	con.tab[1,] <- con.tab1["interesting",]
	con.tab[2,] <- con.tab1["boring",]
#	con.tab[,1] <- con.tab1[,"darkblue"]
#	con.tab[,2] <- con.tab1[,"grey"]
	rownames(con.tab) <- c("interesting", "boring")
#	colnames(con.tab) <- c("darkblue", "grey")
	f.test <- fisher.test(con.tab)
	stats[c, "odds_ratio"] <- f.test$estimate
	stats[c, "p_val"] <- f.test$p.value
}


# if we tighten our inteteresting classification we should expect a tightening of the odds ratio

barplot(stats[,"odds_ratio"], col = cols, names = stats$color, ylab = "odds ratio")

# maybe the thing isn't needed because both are made form the same table 
p_vals <- data.frame(pos = rep(0, length(cols)), sym = rep("", length(cols)), col = cols)
p_vals[stats$p_val < .001,"pos"] = 4.5
p_vals[stats$p_val < .001,"sym"] = "*"

text(p_vals$pos, p_vals$sym, col = cols)

# if we use our secondary criteria we are finding that the odds ratio goes up

# it also seems some of the depletion of old stuff in human is because of the influx of ALU. 




