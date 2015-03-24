# what's going on with comparing all the species on the same accies 

# or at least too speices 

# we can just do the same approach as before but do it to more than just MIR 

# Run it through and get the cutoffs at which enough things work 

# then we can plot it



library(ggplot2)
library(ks)
library(gplots)
library(GenomicRanges)
require(ggbio)


setwd("~/Desktop/Dist_matrix_TE_div/")



# how about remapping two species to human and comparing them to each other






spec <- "Bovine"
rem.un <- "yes"
no.xy <- F

#keep.limit <- 1350000
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
s1name <- paste("count_tables/",spec, "_AllBinCounts.txt", sep = "")
s1 <- read.table(s1name, header = TRUE)
slist <- list(s1)


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


s1 <- slist[[1]][,!(colnames(slist[[1]]) == "Known")]

if(rem.un == "yes"){
	if(length(grep("U", s1$chr)) > 0){s1 <- s1[-(grep( "U", s1$chr)),]}
	if(length(grep("_", s1$chr)) > 0){s1 <- s1[-(grep("_", s1$chr)),]}
	if(length(grep("M", s1$chr)) > 0){s1 <- s1[-(grep("M", s1$chr)),]}
	}
	

# just end up with s1, keep all the same IDs

KnownS1 <- KnownS1[KnownS1[,1] %in% s1$binID,]
KnownS1 <- KnownS1[order(KnownS1[,1]),]
all(KnownS1[,1] == s1$binID)


s1.bins.gr <- GRanges(
		seqnames = Rle(s1$chr),
		ranges = IRanges(start = s1$start, end = s1$end - 1),
		binID = s1$binID
	)



# we read in two species and get there bin information. 


# then we readin the repeat files and sumerise the repeatfiles to only include the repeats we want




# run through the gRanges stuff to get the overlaps, or do the intersects, then an overlap and then aggregate over the file

# maybe we could do a pipe command rather than read the whole thing in
# the system command

# so run awk on the think and read in the data from that 
# so this runs a bit lower on memmory and we don't have to handle big objects


if(spec == "Human"){
	TE.1 <- c("new.LINE",  "new.SINE", "L2", "MIR")
	TE.2 <- c("L1", "ALU" , "L2", "MIR")
}


if(spec == "Bovine"){
  
  # removing bov b for this calculation
	TE.1 <- c(
 #           "new.LINE", 
            "new.LINE", 
            "new.SINE",  
            "new.SINE", 
            "L2", 
            "MIR"
            )
	TE.2 <- c(
#            "BovB", 
            "L1", 
            "BOVTA" , 
            "BOVA2", 
            "L2", 
            "MIR"
            )
}

if(spec =="Horse"){
	TE.1 <- c( "new.LINE",  "new.SINE", "L2", "MIR")
	TE.2 <- c( "L1", "ERE", "L2", "MIR")
}

if(spec == "Dog"){
	TE.1 <- c( "new.LINE",  "new.SINE", "L2", "MIR")
	TE.2 <- c( "L1", "SINEC", "L2", "MIR")
}


if(spec == "Mouse"){
  TE.1 <- c("new.LINE",  "new.SINE", "L2", "MIR")
  TE.2 <- c("L1", "B1" , "L2", "MIR")
}




for(i in 1:length(TE.2)){
	system(paste("grep ", TE.2[i]," ~/Desktop/TE_div/Rep_files/all_chr_rep_",spec ," | cat > ~/Desktop/TE_div/Rep_files/tmp", sep = ""))
	
	TE.map <- read.table("~/Desktop/TE_div/Rep_files/tmp")
	TE.map.gr <- GRanges(
		seqnames = Rle(TE.map[,1]),
		ranges = IRanges(start = TE.map[,2], end = TE.map[,3])
	)
	int.gr <- intersect(s1.bins.gr,TE.map.gr)
	OL <- as.matrix(findOverlaps(s1.bins.gr, int.gr))
	bin.cov <- data.frame(binID = elementMetadata(s1.bins.gr)[[1]][OL[,1]], rep.cov = width(int.gr)[OL[,2]])
	ID.cov <- aggregate(bin.cov$rep.cov, by = list(bin.cov$binID), FUN = sum)
	bin.cov <- data.frame(binID = ID.cov[,1], rep.cov = ID.cov[,2])
	TE.cov <- rep(NA, dim(s1)[1])
	for(e in 1:length(TE.cov)){
		if(s1$binID[e] %in% bin.cov$binID){
			TE.cov[e] <- bin.cov[s1$binID[e] == bin.cov$binID, "rep.cov"]
		}else{
			TE.cov[e] <- 0
		}		
	}
	assign(TE.2[i], TE.cov)
}


for(i in 1:4){
	TE.name <- unique(TE.1)[i]
	TE.cov <- matrix(ncol = length(TE.1[TE.1 == TE.name]), nrow = length(s1$binID))
	for(e in 1:dim(TE.cov)[2]){
		TE.cov[,e] <- get(TE.2[TE.1 == TE.name][e])
	}	
	assign(TE.name, rowSums(TE.cov))
}

new.s1 <- data.frame(new_LINE = new.LINE, new_SINE = new.SINE, L2 = L2, MIR = MIR )

biplot(prcomp(new.s1, scale.=T))

plot(sqrt(new.s1[,3:4]))



read_out <- data.frame(s1[,1:4], new.s1, KnownS1)
colnames(read_out) <- c("Chromosome", "Bin", "Bin_Start", "Bin_End", colnames(new.s1), "binID", "Known")

write.table(read_out, file = paste("~/Desktop/mapping_stats_analysis/raw_data/bin_cov/", spec, "_bin_coverage", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t") 

# all we need to do is read it out and readjust our alignmnet analysis to consider the new data format
# in these cases the data had been filtered to the s.bin format
# so we only need to include our rate calculation and secondary filtrations




