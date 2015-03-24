# what's going on with comparing all the species on the same accies 

# or at least too speices 

# we can just do the same approach as before but do it to more than just MIR 

# Run it through and get the cutoffs at which enough things work 

# then we can plot it


# do we have a mass option?




library(ggplot2)
library(ks)
library(gplots)
library(GenomicRanges)
require(ggbio)


setwd("~/Desktop/Dist_matrix_TE_div/")



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



# how about remapping two species to human and comparing them to each other

# to pick coverage or not is all about getting the right dir
# adding in coverage has confused the remapping 

cov = TRUE
species2 <- c("Horse", "Bovine", "Dog")

for (s in species2){
spec1 <- "Human"
spec2 <- s
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

if(cov){
  s1name <- paste("~/Desktop/mapping_stats_analysis/raw_data/bin_cov/",spec1, "_bin_coverage", sep = "")
  s2name <- paste("~/Desktop/mapping_stats_analysis/raw_data/bin_cov/",spec2, "_bin_coverage", sep = "")
}else{
  s1name <- paste("~/Desktop/Dist_matrix_TE_div/count_tables/",spec1, "_AllBinCounts.txt", sep = "")
  s2name <- paste("~/Desktop/Dist_matrix_TE_div/count_tables/",spec2, "_AllBinCounts.txt", sep = "")
}




s1 <- read.table(s1name, header = TRUE)
s2 <- read.table(s2name, header = TRUE)
slist <- list(s1,s2)

if(!cov){
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
}else{ 
  for(i in seq(along=slist)){
    count <- slist[[i]]
    count[,5:(length(count)-2)] <- ((count[,5:(length(count)-2)]/count$Known) * mb) 
    colnames(count)[1:4] <- c("chr", "binID", "start", "end")
    slist[[i]] <- count
  }  
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



# element keeping, so keep them then remodeling


if(!cov){

human <- data.frame(
		chr = human$chr,
		binID = human$binID,
		start = human$start,
		end = human$end,
		new_LINE = human$LINE_L1, 
		new_SINE = human$SINE1_7SL, 
		MIR = human$SINE2_MIR,
		L2 = human$LINE_L2
)


if(spec2 == "Bovine"){
	bovine <- data.frame(
		chr = bovine$chr,
		binID = bovine$binID,
		start = bovine$start,
		end = bovine$end,
#		new_LINE = bovine$LINE_L1 +  bovine$LINE_RTE_BovB,
    new_LINE = bovine$LINE_L1,
		new_SINE = bovine$SINE2_BOV.tA + bovine$SINE_BOVA2, 
		MIR = bovine$SINE2_MIR,
		L2 = bovine$LINE_L2
	)
}

if(spec2 == "Dog"){
	bovine <- data.frame(
		chr = bovine$chr,
		binID = bovine$binID,
		start = bovine$start,
		end = bovine$end,
		new_LINE = bovine$LINE_L1,  
		new_SINE = bovine$SINE2_CanSINE, 
		MIR = bovine$SINE2_MIR,
		L2 = bovine$LINE_L2
	)
}

if(spec2 == "Horse"){
	bovine <- data.frame(
		chr = bovine$chr,
		binID = bovine$binID,
		start = bovine$start,
		end = bovine$end,
		new_LINE = bovine$LINE_L1 , 
		new_SINE = bovine$SINE2_ERE1 + bovine$SINE2_ERE2 + bovine$SINE2_ERE3 + bovine$SINE2_ERE4, 
		MIR = bovine$SINE2_MIR,
		L2 = bovine$LINE_L2
	)
}

}else{
  human <- human[,1:(length(human) - 1)]
  bovine <- bovine[,1:(length(bovine) - 1)]
}





# we will try the change
# remeber to sort the colums correctly so the whole thing works
b.sbin <- read.table(paste("S_bin/Human_sbin/",spec1,"_aligning_",spec2, "_select", sep = ""), header = TRUE)

#b.sbin <- b.sbin[,c(3,4,1,2)]
#colnames(b.sbin) <- c("S1.bin", "S1.Proportion", "S2.bin", "S2.Proportion")




# after weve made the decisions we can set the species and the cutoff here


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


# Get each species there and save all the human ones 


assign(paste(spec1, spec2, "comp", sep= "_"), s1)
assign(paste(spec2, spec1, "comp", sep= "_"), replot_s1_s2)
}

# now we select a core set of features 
Core.ID <- NULL
NAMES <- NULL
for(s in species2){
  core.a <- get(paste(s, spec1,"comp", sep = "_"))[,"S1_bin"]
#  core.b <- get(paste(spec1, s,"comp", sep = "_"))[,"binID"]
  Core.ID <- c(Core.ID,
               list(core.a)
            #   list(core.b)
               )
  NAMES <- c(NAMES, 
             paste(s, spec1,"comp", sep = "_")
            # paste(spec1, s,"comp", sep = "_")
             )
}
names(Core.ID) <- NAMES
VENN <- venn(Core.ID)

core.id <- Core.ID[[1]][Core.ID[[1]] %in% Core.ID[[2]]]
core.id <- core.id[core.id %in% Core.ID[[3]]]
#core.id <- core.id[core.id %in% Core.ID[[4]]]


# here we can select the names we are interested in 







spec1 <- "Mouse_Human"
spec2 <- "Human_Mouse"

s1 <- get(paste(spec1, "comp", sep = "_"))
replot_s1_s2 <- get(paste(spec2, "comp", sep = "_"))



if(colnames(s1)[1] == "chr"){
  s1 <- s1[s1$binID %in% core.id,]
  s1.comp <- s1[,5:length(s1)]
}else{
  s1 <- s1[s1$S1_bin %in% core.id,]
  s1.comp <- s1[,2:length(s1)]
}


if(colnames(replot_s1_s2)[1] == "chr"){
  replot_s1_s2 <- replot_s1_s2[replot_s1_s2$binID %in% core.id,]
  s2.comp <- replot_s1_s2[,5:length(replot_s1_s2)]
}else{
  replot_s1_s2 <- replot_s1_s2[replot_s1_s2$S1_bin %in% core.id,]
  s2.comp <- replot_s1_s2[,2:length(replot_s1_s2)]
}

# so lets transform the s1.comp
#colMax <- function(X) apply(X, 2, max)
#colM <- colMax(rbind(colMax(s1.comp), colMax(s2.comp)))
#for(i in 1:4){
#	s1.comp[i] <- s1.comp[i]/colM[i]
#	s2.comp[i] <- s2.comp[i]/colM[i]
#}


pca.s1 <- prcomp(s1.comp, scale. = F)
pca.s2 <- prcomp(s2.comp, scale. = F)

reuben.biplot(pca.s1$x, pca.s1$rotation, x.col = "grey")
reuben.biplot(pca.s2$x, pca.s2$rotation, x.col = "grey")


d.te <- rep(NA, dim(s1.comp)[1])
for(i in 1:dim(s1.comp)[1]){
	d.te[i] = dist(rbind(s2.comp[i,], s1.comp[i,]))
}


d.te.n <- d.te
d.te <- scale((d.te))
dense.dist <- density((d.te))
plot((dense.dist))


groups = rep("medium", length(d.te))
groups[d.te > 1.96] = "high"
#groups[d.te < - 1.96] = "low"

t1 <- qplot(s1.comp$new_LINE, s2.comp$new_LINE, col = d.te.n) + scale_color_gradient(low = "lightblue", high = "red") + labs(title = paste(spec1, "vs.", spec2,"for" ,colnames(s1.comp)[1])) + guides(col = FALSE)
t2 <- qplot(s1.comp$new_SINE, s2.comp$new_SINE, col = d.te.n) + scale_color_gradient(low = "lightblue", high = "red") + labs(title = paste(spec1, "vs.", spec2,"for" , colnames(s1.comp)[2])) + guides(col = FALSE)
t3 <- qplot(s1.comp$MIR, s2.comp$MIR, col = d.te.n) + scale_color_gradient(low = "lightblue", high = "red") + labs(title = paste(spec1, "vs.", spec2,"for" , colnames(s1.comp)[3])) + guides(col = FALSE)
t4 <- qplot(s1.comp$L2, s2.comp$L2, col = d.te.n) + scale_color_gradient(low = "lightblue", high = "red")+ labs(title = paste(spec1, "vs.", spec2,"for" , colnames(s1.comp)[4])) + guides(col = FALSE)

multiplot(t1,t2,t3,t4, cols = 2)


# lets try plotting distribtuions in 2d


for(i in 1:4){
	dat <- data.frame(S1 = s1.comp[,i], S2 = s2.comp[,i], groups = groups)
	
	assign(
		paste("p", i , sep =""), 
		ggplot(dat, aes(x = S1, y = S2, col = groups)) +  geom_point() + labs(title = paste(spec1, "vs.", spec2,"for" , colnames(s1.comp)[i])) + guides(col = FALSE) + xlim(min(dat$S1), max(dat$S1)) + ylim(min(dat$S2), max(dat$S2))
	)
	
	assign(
		paste("q", i , sep =""),
		ggplot(dat, aes(x = S1, y = S2, col = groups) ) + geom_density2d()  + labs(title = paste(spec1, "vs.", spec2,"for" , colnames(s1.comp)[i])) + guides(col = FALSE) + xlim(min(dat$S1), max(dat$S1)) + ylim(min(dat$S2), max(dat$S2))	
	)
		
	assign(
		paste("r", i , sep =""),
		ggplot(dat, aes(x = S1, y = S2, col = groups) ) + geom_density2d()  + labs(title = paste(spec1, "vs.", spec2,"for" , colnames(s1.comp)[i])) + guides(col = FALSE) + xlim(min(dat$S1), max(dat$S1)) + ylim(min(dat$S2), max(dat$S2))	 + geom_point()
	)

}



multiplot(p1, p2, p3, p4, cols = 2)

multiplot(q1, q2, q3, q4, cols = 2)

multiplot(r1, r2, r3, r4, cols = 2)





# is there a cooler way to show these results? # we pr obably should get it so we can either put the lines in or go one dimensional
# or put the specie sname on th epca plots


# if we are going to look at differences it may be worht doing a pooled difference spectrum
# this way species that look like each other wont be too differnt 



spec.list <- paste(species2, "Human", sep = "_")
spec.list.a <- t(combn(spec.list, 2))
spec.list.human <- cbind(paste("Human", species2, sep = "_"), spec.list)
spec.list <- rbind(spec.list.human, spec.list.a)


# we'll just look at the other species for a moment
spec.list <- spec.list.a


distances <- matrix(ncol = nrow(spec.list), nrow = length(core.id))
d.all <- NULL
for(s in 1:nrow(spec.list)){
  
  
  spec1 <- spec.list[s,1]
  spec2 <- spec.list[s,2]
  
  s1 <- get(paste(spec1, "comp", sep = "_"))
  replot_s1_s2 <- get(paste(spec2, "comp", sep = "_"))
  
  if(colnames(s1)[1] == "chr"){
    s1 <- s1[s1$binID %in% core.id,]
    s1.comp <- s1[,5:length(s1)]
  }else{
    s1 <- s1[s1$S1_bin %in% core.id,]
    s1.comp <- s1[,2:length(s1)]
  }
  
  
  if(colnames(replot_s1_s2)[1] == "chr"){
    replot_s1_s2 <- replot_s1_s2[replot_s1_s2$binID %in% core.id,]
    s2.comp <- replot_s1_s2[,5:length(replot_s1_s2)]
  }else{
    replot_s1_s2 <- replot_s1_s2[replot_s1_s2$S1_bin %in% core.id,]
    s2.comp <- replot_s1_s2[,2:length(replot_s1_s2)]
  }
  
  d.te <- rep(NA, dim(s1.comp)[1])
  for(i in 1:dim(s1.comp)[1]){
    d.te[i] = dist(rbind(s2.comp[i,], s1.comp[i,]))
  }
  distances[,s] <- d.te
  d.all <- rbind(d.all, d.te)
}


dd2 <- distances
distances <- log(distances)
dd <- d.all
d.all <- log(d.all)
plot(density((d.all)))
points(mean(d.all) + 1.96*sd(d.all), 0, col = 2)
points(mean(d.all) - 1.96*sd(d.all), 0, col = 2)
points(mean(d.all), 0, col= 3)
# so now we have a real spectrum of differences between all the genomes

colnames(distances) <- paste(spec.list[,1], spec.list[,2], sep="-")
groups <- distances
groups[groups > mean(d.all) + 1.96*sd(d.all)] <- "high"
groups[groups != "high"] <- "medium"
groups[distances < mean(d.all) - 1.96*sd(d.all)] <- "low"
# now we have each group sorted,
# we also have a spectrum of scores that encompases all species comparisons
# now we don't have to worry about the context 


# here we have to reintroduce our groups like before 
s= 3

for(s in 1:nrow(spec.list)){

spec1 <- spec.list[s,1]
spec2 <- spec.list[s,2]

s1 <- get(paste(spec1, "comp", sep = "_"))
replot_s1_s2 <- get(paste(spec2, "comp", sep = "_"))
if(colnames(s1)[1] == "chr"){
  s1 <- s1[s1$binID %in% core.id,]
  s1.comp <- s1[,5:length(s1)]
}else{
  s1 <- s1[s1$S1_bin %in% core.id,]
  s1.comp <- s1[,2:length(s1)]
}
if(colnames(replot_s1_s2)[1] == "chr"){
  replot_s1_s2 <- replot_s1_s2[replot_s1_s2$binID %in% core.id,]
  s2.comp <- replot_s1_s2[,5:length(replot_s1_s2)]
}else{
  replot_s1_s2 <- replot_s1_s2[replot_s1_s2$S1_bin %in% core.id,]
  s2.comp <- replot_s1_s2[,2:length(replot_s1_s2)]
}

# also we need a table of matches so as to pick the right column
col.pick <- paste(spec.list[s,1], spec.list[s,2], sep="-")





for(i in 1:4){
  dat <- data.frame(S1 = s1.comp[,i], S2 = s2.comp[,i], groups = groups[,col.pick], col = distances[,col.pick])
  
  t <- ggplot(dat,aes(x = S1, y = S2, col = col)) + geom_point() + scale_color_gradient(limits = c(min(d.all), max(d.all)), low = "lightblue", high = "red") + labs(title = paste(col.pick ,colnames(s1.comp)[i])) + guides(col = FALSE)
  assign(paste("t", i , sep =""),t)


  print(i)
  dat <- data.frame(S1 = s1.comp[,i], S2 = s2.comp[,i], groups = groups[,col.pick])

  assign(
    paste("p", i , sep =""), 
    ggplot(dat, aes(x = S1, y = S2, col = groups)) +  geom_point() + labs(title = paste(col.pick , colnames(s1.comp)[i])) + guides(col = FALSE) + xlim(min(dat$S1), max(dat$S1)) + ylim(min(dat$S2), max(dat$S2))
  )
  
  assign(
    paste("q", i , sep =""),
    ggplot(dat, aes(x = S1, y = S2, col = groups)) + geom_density2d()  + labs(title = paste(col.pick , colnames(s1.comp)[i])) + guides(col = FALSE) + xlim(min(dat$S1), max(dat$S1)) + ylim(min(dat$S2), max(dat$S2))	
  )
  
  assign(
    paste("r", i , sep =""),
    ggplot(dat, aes(x = S1, y = S2, col = groups)) + geom_density2d()  + labs(title = paste(col.pick , colnames(s1.comp)[i])) + guides(col = FALSE) + xlim(min(dat$S1), max(dat$S1)) + ylim(min(dat$S2), max(dat$S2))	 + geom_point()
  )
  
}


pdf(file = paste("~/Desktop/mapping_stats_analysis/plots_cov/pairwise_distribution/", col.pick, "_no_human.pdf", sep = ""), onefile = TRUE)
multiplot(t1, t2, t3, t4, cols = 2)
multiplot(p1, p2, p3, p4, cols = 2)
multiplot(q1, q2, q3, q4, cols = 2)
multiplot(r1, r2, r3, r4, cols = 2)
dev.off()

}

# maybe it is worht showing the differences attributed to human and the differences attributed to everything else 
# basicly there is two peaks when we analyse human, and i trhing this is becasue there is greater differnce over all 


# but if we do the comparison in whihc TE families get there own axis this way we can compare all the families we are intersted in
# it works a bit like a kmer spectrum
# this way we can do  a PCA and plot a species in the same PC space. 
# with the rotation













i = 4
dat <- data.frame(S1 = s1.comp[,i], S2 = s2.comp[,i], groups = groups)
pca.line <- prcomp(dat[,1:2])
dat <- data.frame(PC1 = pca.line$x[,1], PC2 = pca.line$x[,2], groups = groups)
p = ggplot(dat, aes( scale(PC2), fill = groups)) 
p = p +  geom_density(alpha = .7) 
p = p + labs(title = colnames(s1.comp)[i]) 
p




i = 2
dat <- data.frame(S1 = s1.comp[,i], S2 = s2.comp[,i], groups = groups)
pca.line <- prcomp(dat[,1:2], scale. = F)
dat <- data.frame(PC1 = pca.line$x[,1], PC2 = pca.line$x[,2], groups = groups)

q = ggplot(dat, aes(scale(PC2), scale(PC1), col = groups)) 
#q = q +  geom_density2d(alpha = .7) 
q = q + geom_point()
q = q + labs(title = colnames(s1.comp)[i]) 
q = q + geom_hline(yintercept= + 1.96)
q = q + geom_hline(yintercept= - 1.96)
q = q + geom_vline(xintercept= + 1.96)
q = q + geom_vline(xintercept= - 1.96)
q






# so we do the rotation and have a look

# how do we get the rotation lines in though

# what ever pca does to calculate the axis




colors <- rep("grey", length(groups))
colors[groups == "high"] = "orange"


reuben.biplot(pca.s1$x, pca.s1$rotation, x.col = colors)
reuben.biplot(pca.s2$x, pca.s2$rotation, x.col = colors)

biplot(prcomp(t(s2.comp)))



enrich_samples = NULL
enrich_samples_e = NULL

for(i in 1:4){
	# getting the TE we're interested in
	dat <- data.frame(S1 = s1.comp[,i], S2 = s2.comp[,i], groups = groups, binID = s1$binID)
	dat[,1:2] <- scale(dat[,1:2])
	pca.line <- prcomp(dat[dat$group == "medium",1:2])	
	
	# working out the pca rotations and how to rotate the data
	
	
	# selecting the right PC has been a probelem so far, we obviously need a method that doesn't sometimes pick the wrong one
	
	
	if(cor(pca.line$x[,2], d.te.n[dat$group == "medium"]) < cor(pca.line$x[,1], d.te.n[dat$group == "medium"])){
		PC = 1
	}else{
		PC = 2
	}
	# we have an abline and rotate it
	line.try = NULL
	line.try <- matrix(data = c(rep(pca.line$sdev[PC] * -2.58, 10), seq(from = 1 , to = 10)), ncol = 2)
	for(p in 1:nrow(line.try)){
		line.try[p,1:2] <- line.try[p,1:2] %*% pca.line$rotation
	}
	if(cor(line.try[,1], line.try[,2]) < 0 ){
		flip = TRUE
	}else{
		flip = FALSE
	}
	
	# drawing the upper and lower bound of interesting regions
	line1 <- matrix(data = c(rep(pca.line$sdev[PC] * -2.58, 10), seq(1: 10)), ncol = 2)
	if(flip){
		line1 <- line1[,c(2,1)]
	}
	for(p in 1:nrow(line1)){
		line1[p,1:2] <- line1[p,1:2] %*% pca.line$rotation
	}
	line1 <- lm(line1[,2] ~ line1[,1])
	
	line2 <- matrix(data = c(rep(pca.line$sdev[PC] * 2.58, 10), seq(1: 10)), ncol = 2)
	if(flip){
		line2 <- line2[,c(2,1)]
	}
	for(p in 1:nrow(line2)){
		line2[p,1:2] <- line2[p,1:2] %*% pca.line$rotation
	}
	line2 <- lm(line2[,2] ~ line2[,1])
	
	# plotting the findings
	plot(dat[,1:2], col = colors, pch = 16, main = colnames(s1.comp)[i])
	abline(line1, col = 2)
	abline(line2, col = 2)
	
	# have to make sure the actual y is greaterr than the expected y
	# makeing sure we have the right direcgtion is important with our lines
	# because sometimes line1 is the lowerbound and line2 is the upperbound
	
	
	# getting the bins with over and under TE amounts
	e.y.s1 <- line2$coefficients[2] * dat$S1 - line2$coefficients[1]
	e.y.s2 <- line1$coefficients[2] * dat$S1 - line1$coefficients[1]
	if(nrow(dat[dat$S2 > e.y.s1,]) < nrow(dat[dat$S2 < e.y.s1,])){
		H.s2 <- dat[dat$S2 > e.y.s1,]
		H.s1 <- dat[dat$S2 < e.y.s2,]
	}else{
		H.s1 <- dat[dat$S2 < e.y.s1,]
		H.s2 <- dat[dat$S2 > e.y.s2,]
	}
	
	assign(paste("H.s1", colnames(s1.comp)[i], sep = "_"), H.s1$binID)
	assign(paste("H.s1", colnames(s1.comp)[i],"e", sep = "_"), H.s1$binID[as.character(H.s1$group) == "high"])
	
	assign(paste("H.s2", colnames(s1.comp)[i], sep = "_"), H.s2$binID)
	assign(paste("H.s2", colnames(s1.comp)[i],"e", sep = "_"), H.s2$binID[as.character(H.s2$group) == "high"])
	
	
	
	enrich_samples = c(enrich_samples, paste("H.s1", colnames(s1.comp)[i], sep = "_"))	
	enrich_samples = c(enrich_samples, paste("H.s2", colnames(s1.comp)[i], sep = "_"))
	
	if(nrow(H.s1[H.s1$groups == "high",]) > 0){
		enrich_samples_e = c(enrich_samples_e, paste("H.s1", colnames(s1.comp)[i], "e", sep = "_"))
	} 
	
	if(nrow(H.s2[H.s2$groups == "high",]) > 0){
		enrich_samples_e = c(enrich_samples_e, paste("H.s2", colnames(s1.comp)[i], "e", sep = "_"))
	} 

	# from this we need to gather lists of things with interesting TE content 
	# so make two lists for each think, the H.s1 and the H.s1.e
	# and also get the binID
	
	
	# then we can do enrichment type stuff
	
}


# now i have bin IDs and I have plots i can start to print fpr






# can now get confidence intervals on rotations
# make plots 




# so we have PC2 and PC1, we want the PC that eiter correlates more or less with d





# so we have cut offs now for interesting TE.content 
# we just need a way of getting those bins


# so lets run these IDs, through a gene enrichment analysis

big.df <- read.table("~/Desktop/big.df.txt", header = T)



# so we have all the genes and the binIDs, we need to cycle through the module colors and regions of interest
# it will pretty much be a nested loop thing , but at the end we'll get a list of gnenes in a module that is interesting.

for(i in 1:length(enrich_samples)){
	H.s <- get(enrich_samples[i])
	
	ID.list <- rep("B", nrow(big.df))
	ID.list[big.df$binID %in% H.s] <- "A"
	
	assign(paste(enrich_samples[i], "IDs", sep = "_"), ID.list)
}


for(i in 1:length(enrich_samples_e)){
	
	ID.list <- rep("B", nrow(big.df))
	ID.list[big.df$binID %in% H.s] <- "A"
	
	assign(paste(enrich_samples_e[i], "IDs", sep = "_"), ID.list)
}

enrich_samples_ID <- paste(enrich_samples, "IDs", sep = "_")

enrich_samples_e_ID <- paste(enrich_samples_e, "IDs", sep = "_")

# with TE content "A" is ther interesting one
# it's because R puts things in alphabetical order


colors <- as.character(unique(big.df$mod.col))

for(i in 1:length(enrich_samples_ID)){
	stats = data.frame(color = colors, odds_ratio = NA, p_val = NA)
	for(c in 1:length(colors)){
		col <- rep("grey", nrow(big.df))
		col[big.df$mod.col == colors[c]] = "color"
		pre.tab <- data.frame(TE.content = get(enrich_samples_ID[i]), module = col)
		con.tab <- table(pre.tab)
		fisher <- fisher.test(con.tab)
		stats[c, "odds_ratio"] <- fisher$estimate
		stats[c, "p_val"] <- fisher$p.value

	}
	assign(paste("stats",enrich_samples_ID[i], sep ="_"), stats)
}




for(i in 1:length(enrich_samples_e_ID)){
	stats = data.frame(color = colors, odds_ratio = NA, p_val = NA)
	for(c in 1:length(colors)){
		col <- rep("grey", nrow(big.df))
		col[big.df$mod.col == colors[c]] = "color"
		pre.tab <- data.frame(TE.content = get(enrich_samples_ID[i]), module = col)
		con.tab <- table(pre.tab)
		fisher <- fisher.test(con.tab)
		stats[c, "odds_ratio"] <- fisher$estimate
		stats[c, "p_val"] <- fisher$p.value

	}
	assign(paste("stats",enrich_samples_e_ID[i], sep ="_"), stats)
}


# it likly seems that none of the genes from bins with interesting TE content accosiate with any modules, still nees to chack the ones that are just generally different!!


highID <- s1$binID[groups == "high"]
highID_2 <- rep("B", nrow(big.df))
highID_2[big.df$binID %in% highID] <- "A"

stats = data.frame(color = colors, odds_ratio = NA, p_val = NA)
for(c in 1:length(colors)){
	col <- rep("grey", nrow(big.df))
	col[big.df$mod.col == colors[c]] = "color"
	pre.tab <- data.frame(TE.content = highID_2, module = col)
	con.tab <- table(pre.tab)
	fisher <- fisher.test(con.tab)
	stats[c, "odds_ratio"] <- fisher$estimate
	stats[c, "p_val"] <- fisher$p.value

}

# using the gene wide approach we aren't really seeing much, I'm not sure if our gene list has too many genes that just aren't involved in the analysis
# also there is probably some inherent bias when going binwise instead of gene wise

# we are looking at wheater or not genes that come from a bin of interest, accosiate with a particular module. 
# so at lwast for human and dog, genes from a particular bin don't tend to accosiate with any modules




par(mar=c(10,5,5,3))

barplot(stats[,2], col = as.character(stats[,1]), ylab = "odds ratio")


	
labs = stats[,1]

end_point = 0.5 + nrow(stats) +3


text(seq(.7,end_point,by=1.2), par("usr")[3] - .05, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = labs, cex=1)

pv1 <- as.character(round(stats[,3], digit = 2))
pv <- rep("", length(pv1))
pv[stats[,3] < .05] <- "*"


text(seq(.7,end_point,by=1.2), par("usr")[3] + 2.3, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = pv, cex=1)




# what is going on with this blue module, why is it so depleted of things with anything new



par(mar=c(10,5,5,3))

barplot(stats_H.s1_new_SINE_e_IDs[,2], col = as.character(stats_H.s1_new_SINE_e_IDs[,1]), ylab = "odds ratio", ylim = c(0,2))


	
labs = stats_H.s1_new_SINE_e_IDs[,1]

end_point = 0.5 + nrow(stats_H.s1_new_SINE_e_IDs) +3


text(seq(.7,end_point,by=1.2), par("usr")[3] - .05, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = labs, cex=1)

pv1 <- as.character(round(stats_H.s1_new_SINE_e_IDs[,3], digit = 2))
pv <- rep("", length(pv1))
pv[stats_H.s1_new_SINE_e_IDs[,3] < .05] <- "*"


text(seq(.7,end_point,by=1.2), par("usr")[3] + 2.3, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = pv, cex=1)







points(pca.line$rotation, col = 2, pch = 16)
text(pca.line$rotation, c("PC1", "PC2"))
points(0,0, col = 4, pch = 16)

#points(pca.line$sdev %*% pca.line$rotation, col = 3, pch = 16)

points(c(1,pca.line$sdev[2]) %*% pca.line$rotation, col = 3)
pca.line



# take these points and multiply them by the rotation 
points(c(0,pca.line$sdev[1] * -2)  %*% pca.line$rotation, col = 5)
points(c(pca.line$sdev[1] * 2, 2) %*% pca.line$rotation, col = 5)




# this is harder than first thought




# do the pca plot with the blue module highlighted 






B <- unique(big.df[big.df$mod.col == "blue", "binID"])

colours <- s1$binID
colours[colours %in% B] <- "blue"
colours[colours != "blue"] <- "grey"

reuben.biplot(pca.s1$x, pca.s1$rotation, x.col=colours)

