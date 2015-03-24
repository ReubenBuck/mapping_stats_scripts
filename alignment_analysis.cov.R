


library(ggplot2)
library(ks)
library(gplots)
library(GenomicRanges)
require(ggbio)


setwd("~/Desktop/mapping_stats_analysis/")


rm(list = ls())


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
spec2 <- "Mouse"
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
	

human <- s1[]
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



cutoff <- c(0, .05,.1, .15,.2)
#cutoff <- .1

results.p <- data.frame(row.names = as.character(1:(length(s2)-4)))
results.p.s1 <- data.frame(row.names = as.character(1:(length(s1)-4)))
results.p.s2 <- data.frame(row.names = as.character(1:(length(s2)-4)))
results.s <- data.frame(row.names = as.character(1:(length(s2)-4)))
results.s.s1 <- data.frame(row.names = as.character(1:(length(s1)-4)))
results.s.s2 <- data.frame(row.names = as.character(1:(length(s2)-4)))
results.m <- data.frame(cut_off = cutoff, cor = 0, pval = 0, ref_bin_no = 0, query_bin_no = 0)
for(c in 1:length(cutoff)){
	
	

S.bin_s1_s2 <- b.sbin
s1 <- human
s2 <- bovine


	S.bin_s1_s2 <- S.bin_s1_s2[S.bin_s1_s2[,4] > cutoff[c],]
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
#	TE.s2 <- s2TE[,2:(length(s2TE)-1)]

	
	
	# so in here somewhere we can look 
	
	
	# paste the cutoff value to the resulting dataframes to do multiple heatmaps
	assign(paste("s2_replot", cutoff[c], sep = "_"), replot_s1_s2)
	assign(paste("s1", cutoff[c], sep = "_"), s1)
	
	
	tests <- data.frame(TEs = colnames(TE.s2), pval = 0, stat = 0)
	for(i in 1:dim(TE.s2)[2]){
		Ttest <- ks.test(TE.s2[,i], TE.replot[,i])
		tests[i,2] <- Ttest$p.value
		tests[i,3] <- Ttest$statistic
	}
	results.p[,c] <- tests[,2]
	results.s[,c] <- tests[,3]
	results.m[c,2:5] <- c(0 ,0, length(unique(S.bin_s1_s2$S1.bin)), length(unique(S.bin_s1_s2$S2.bin))) 
	
	
	# look at TE.s1 and do a test for that 
	
	tests <- data.frame(TEs = colnames(s1[,5:length(s1)]), pval = 0, stat = 0)
	for(i in 5:dim(s1)[2]){
		Ttest <- ks.test(human[,i], s1[,i])
		tests[i - 4,2] <- Ttest$p.value
		tests[i - 4,3] <- Ttest$statistic
	}
	results.p.s1[,c] <- tests[,2]
	results.s.s1[,c] <- tests[,3]

	# look at TE.s2 compared to the unsampled distribution and do a test for that 
	
	tests <- data.frame(TEs = colnames(TE.s2), pval = 0, stat = 0)
	for(i in 1:dim(TE.s2)[2]){
		Ttest <- ks.test(bovine[,i + 4], TE.replot[,i])
		tests[i,2] <- Ttest$p.value
		tests[i,3] <- Ttest$statistic
	}
	results.p.s2[,c] <- tests[,2]
	results.s.s2[,c] <- tests[,3]
	
	
	
	# pictures at each cutoff, we do the pca of cow and human and we colour in the dots 
	# use the reuben biplot function to get the heatmaps all we really need to save is the colours,
	
	cutoff.col = rep(3, length(human$binID))
	cutoff.col[human$binID %in% s1$binID] <- 1	
	assign(paste("cutoff.col.spec1_", cutoff[c], sep = ""), cutoff.col)
	
	cutoff.col = rep(3, length(bovine$binID))
	cutoff.col[bovine$binID %in% unique(S.bin_s1_s2$S2.bin)] <- 1	
	assign(paste("cutoff.col.spec2_", cutoff[c], sep = ""), cutoff.col)

}
colnames(results.p) <- paste("cutoff", cutoff, sep="_")
rownames(results.p) = colnames(TE.s2)
rownames(results.s) = colnames(TE.s2)

colnames(results.p.s1) <- paste("cutoff", cutoff, sep="_")
rownames(results.p.s1) = colnames(s1)[5:length(s1)]

colnames(results.p.s2) <- paste("cutoff", cutoff, sep="_")
rownames(results.p.s2) = colnames(s2)[5:length(s2)]


pdf(file = paste("~/Desktop/mapping_stats_analysis/plots_cov/bar_line_plot/", spec2,"_stats_",chr_main,".pdf", sep = ""),onefile = TRUE)

barplot(t((results.p)), beside = T, legend = T , main = paste("remodeled", spec2, "sampling compared to sample", chr_main), cex.names = .6, xlab = "TE family", ylab = "p-value", las = 2, ylim = c(0,1))
abline(h = (.05), col = 2)

barplot(t((results.p.s1)), beside = T, legend = T , main = paste(spec1, "sampling compared to total", chr_main), cex.names = .6, xlab = "TE family", ylab = "p-value", las = 2, ylim = c(0,1))
abline(h = (.05), col = 2)

barplot(t((results.p.s2)), beside = T, legend = T , main = paste("remodeled", spec2, "sampling compared to total", chr_main), cex.names = .6, xlab = "TE family", ylab = "p-value", las = 2, ylim = c(0,1))
abline(h = (.05), col = 2)

plot(results.m[,1], results.m[,4],  ylim = c(0,2000), main = paste("1:1 mappings across" ,spec1, "and", spec2, chr_main), type = "l", ylab = "number of bins", col = 2, xlab = "cutoff")
lines(results.m[,1], results.m[,5], col = 3)
legend("topright",c( paste(spec1, "bins"), paste(spec2, "bins")),fill=c(2:3))

dev.off()

# now we get a sense of the TE content of the bins we are cutting out
pca.spec1 <- prcomp(human[,5:length(human)], scale. = FALSE)
pca.spec2 <- prcomp(bovine[,5:length(bovine)], scale. = FALSE)
pdf(file = paste("~/Desktop/mapping_stats_analysis/plots_cov/pca_plot/", spec2,"_biplot_",chr_main,"_no_scale",".pdf", sep = ""),onefile = TRUE)

for(c in 1:length(cutoff)){
	
	reuben.biplot(
	pca.spec1$x, 
	pca.spec1$rotation, 
	x.col = get(paste("cutoff.col.spec1_",cutoff[c], sep = "")),
	main = paste(spec1, "biplot compared to",spec2,"at cutoff", cutoff[c],"for", chr_main,"no scale")
	)
}

for(c in 1:length(cutoff)){
	reuben.biplot(pca.spec2$x, 
	pca.spec2$rotation, 
	x.col = get(paste("cutoff.col.spec2_",cutoff[c], sep = "")),
	main = paste(spec2, "biplot at cutoff", cutoff[c],"for", chr_main,"no scale")
	)
}

dev.off()

pca.spec1 <- prcomp(human[,5:length(human)], scale. = TRUE)
pca.spec2 <- prcomp(bovine[,5:length(bovine)], scale. = TRUE)
pdf(file = paste("~/Desktop/mapping_stats_analysis/plots_cov/pca_plot/", spec2,"_biplot_",chr_main,"_scaled",".pdf", sep = ""),onefile = TRUE)

for(c in 1:length(cutoff)){
	
	reuben.biplot(
	pca.spec1$x, 
	pca.spec1$rotation, 
	x.col = get(paste("cutoff.col.spec1_",cutoff[c], sep = "")),
	main = paste(spec1, "biplot compared to",spec2,"at cutoff", cutoff[c],"for", chr_main,"scaled")
	)
}

for(c in 1:length(cutoff)){
	reuben.biplot(pca.spec2$x, 
	pca.spec2$rotation, 
	x.col = get(paste("cutoff.col.spec2_",cutoff[c], sep = "")),
	main = paste(spec2, "biplot at cutoff", cutoff[c],"for", chr_main,"scaled")
	)
}

dev.off()





## so the next part is to do the karyogram maybe we could use thecolours from before to do it



web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", ucsc1, "/database/chromInfo.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
chrom_info <- read.delim(textConnection(txt), header = FALSE)
chrom_info_1 <- chrom_info[,2]
names(chrom_info_1) <- chrom_info[,1]
chrom_info_1 <- chrom_info_1[names(chrom_info_1) %in% unique(human$chr)]

#web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", ucsc2, "/database/chromInfo.txt.gz", sep = "")
#con <- gzcon(url(web))
#txt <- readLines(con)
#chrom_info <- read.delim(textConnection(txt), header = FALSE)
#chrom_info_2 <- chrom_info[,2]
#names(chrom_info_2) <- chrom_info[,1]
#chrom_info_2 <- chrom_info_2[names(chrom_info_2) %in% unique(bovine$chr)]

#IDs <- NULL
#for(i in 1:length(chrom_info_2)){
#	if(max(bovine[bovine$chr == names(chrom_info_2[i]),"end"]) > chrom_info_2[i]){
#		ID <- bovine[bovine$chr == names(chrom_info_2[i]),][(bovine[bovine$chr == names(chrom_info_2[i]),"end"]) == max(bovine[bovine$chr == names(chrom_info_2[i]),"end"]),"binID"]
#	}
#	IDs <- c(ID,IDs)
#}
#bovine <- bovine[!(bovine$binID %in% IDs),]

spec1.gr <- GRanges(
	seqnames = Rle(human$chr),
	ranges = IRanges(
		start = human$start,
		end = human$end
	),
	seqlengths = chrom_info_1
)

# species and cutoff level 


for(c in 1:length(cutoff)){
	q <- autoplot(spec1.gr,layout = "karyogram",fill = "#AEAEFF", col = "#4141FF", xlab = paste(spec1, "karyogram compared to",spec2,"at cutoff", cutoff[c],"for", chr_main))
	
	if(any(get(paste("cutoff.col.spec1_",cutoff[c],sep = "")) == 3)){
		q <- q + layout_karyogram(spec1.gr[get(paste("cutoff.col.spec1_",cutoff[c],sep = "")) == 3], cgeom = "rect", ylim = c(0,10), fill = "#AEAEFF", col = "#FF2C2C")
		}
	
	ggsave(filename = paste("~/Desktop/mapping_stats_analysis/plots_cov/karyogram/",spec2,"/", spec1, "_karyogram_compared_to_",spec2,"_at_cutoff_", cutoff[c],"_for_", chr_main,".pdf", sep = ""))
}




spec2.gr <- GRanges(
	seqnames = Rle(bovine$chr),
	ranges = IRanges(
		start = bovine$start,
		end = bovine$end
	),
#	seqlengths = chrom_info_2 +1
)

for(c in 1:length(cutoff)){
	q <- autoplot(spec2.gr,layout = "karyogram",fill = "#AEAEFF", col = "#4141FF", xlab = paste(spec2, "karyogram at cutoff", cutoff[c],"for", chr_main))
	if(any(get(paste("cutoff.col.spec2_",cutoff[c],sep = "")) == 3)){
		q <- q + layout_karyogram(spec2.gr[get(paste("cutoff.col.spec2_",cutoff[c],sep = "")) == 3], cgeom = "rect", ylim = c(0,10), fill = "#FF8A8A", col = "#FF2C2C")
		}
	
	
	ggsave(filename = paste("~/Desktop/mapping_stats_analysis/plots_cov/karyogram/", spec2,"/", spec2,"_karyogram_at_",cutoff[c],"_for_",chr_main,".pdf", sep = ""))
    
}


# this works for mouse over these TE families
# ks.test(TE.s2[,c(1,2,4,5,8,11,12,13,14,10)], TE.replot[,c(1,2,4,5,8,11,12,13,14,10)])


# this could be automated and we could have a cutoff
# so tomorrow I'll do a couple of species.   

# run it through on differnt cutoffs


my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 19)
my_p_val_palette <- colorRampPalette(c( "white", "black"))(n = 101)
pv <- seq(1,0, by = -.1)
pv.text <- rep("", length(pv))
pv.text[as.integer(quantile(1:(length(pv))))] <- as.character(pv[as.integer(quantile(1:(length(pv))))])
# now i can get the row side colours to give me a p,val




pdf(file = paste("~/Desktop/mapping_stats_analysis/plots_cov/heatmap/", spec2,"_heatmap_",chr_main,".pdf", sep = ""),onefile = TRUE) 
for(c in 1:length(cutoff)){
	A <- get(paste("s1", cutoff[c], sep = "_"))
	B <- get(paste("s2_replot", cutoff[c], sep = "_"))
	heatmap.2(cor(A[,5:length(s1)], B[,2:length(replot_s1_s2)], method = "spearman"), 
			scale = "none",
			col = my_palette,
			symkey = TRUE,
			breaks = seq(from = -1.1, to = 1.1, length = 20),
			trace = "none",
			symbreaks = TRUE,
			margins = c(8,8),
			xlab = spec2,
			ylab = spec1,
			main = paste("spatial correlations \nof TE families at cutoff", cutoff[c],"\n", chr_main),
			density.info = "histogram",
			denscol=1,
			RowSideColors = my_p_val_palette[as.integer(results.p.s1[,c] * 100) + 1],
			ColSideColors = my_p_val_palette[as.integer(results.p.s2[,c] * 100) + 1]
			)
	legend("topright", pv.text,fill = my_p_val_palette[pv*100 +1], y.intersp = .5, border = 0, box.col = 0, title = "P-value", cex = .8)
}
dev.off()



# spec1 self correlation cut off

heatmap.2(cor(A[,5:length(s1)], method = "spearman"), 
			scale = "none",
			col = my_palette,
			symkey = TRUE,
			breaks = seq(from = -1, to = 1, length = 20),
			trace = "none",
			symbreaks = TRUE,
			margins = c(8,8),
#			xlab = spec2,
#			ylab = spec1,
			main = paste(spec1,"spatial correlations \nof TE families at cutoff", cutoff[c]),
			denscol=1,
			RowSideColors = my_p_val_palette[as.integer(results.p.s1[,c] * 100) + 1],
			ColSideColors = my_p_val_palette[as.integer(results.p.s1[,c] * 100) + 1]
			)


# re plotted genome self correlation of sample 

heatmap.2(cor(B[,2:length(replot_s1_s2)], method = "spearman"), 
			scale = "none",
			col = my_palette,
			symkey = TRUE,
			breaks = seq(from = -1, to = 1, length = 20),
			trace = "none",
			symbreaks = TRUE,
			margins = c(8,8),
#			xlab = spec2,
#			ylab = spec1,
			main = paste(spec2, "spatial correlations \nof TE families at cutoff", cutoff[c]),
			denscol=1,
			RowSideColors = my_p_val_palette[as.integer(results.p.s2[,c] * 100) + 1],
			ColSideColors = my_p_val_palette[as.integer(results.p.s2[,c] * 100) + 1]
			)






# original species 2 data self correlation


heatmap.2(cor(bovine[,5:length(bovine)], method = "spearman"), 
			scale = "none",
			col = my_palette,
			symkey = TRUE,
			breaks = seq(from = -1, to = 1, length = 20),
			trace = "none",
			symbreaks = TRUE,
			margins = c(8,8),
#			xlab = spec2,
#			ylab = spec1,
			main = paste(spec2, "spatial correlations \n no processing"),
			denscol=1,
#			RowSideColors = my_p_val_palette[as.integer(results.p.s2[,c] * 100) + 1],
#			ColSideColors = my_p_val_palette[as.integer(results.p.s2[,c] * 100) + 1]
			)









# this one is sample to sample p values


heatmap.2(cor(A[,5:length(s1)], B[,2:length(replot_s1_s2)], method = "spearman"), 
			scale = "none",
			col = my_palette,
			symkey = TRUE,
			breaks = seq(from = -1.1, to = 1.1, length = 20),
			trace = "none",
			symbreaks = TRUE,
			margins = c(8,8),
			xlab = spec2,
			ylab = spec1,
			main = paste("spatial correlations \nof TE families at cutoff", cutoff[c]),
			density.info = "histogram",
			denscol=1,
			RowSideColors = my_p_val_palette[as.integer(results.p.s1[,c] * 100) + 1],
			ColSideColors = my_p_val_palette[as.integer(results.p[,c] * 100) + 1]
			)
legend("topright", pv.text,fill = my_p_val_palette[pv*100 +1], y.intersp = .5, border = 0, box.col = 0, title = "P-value", cex = .8)


