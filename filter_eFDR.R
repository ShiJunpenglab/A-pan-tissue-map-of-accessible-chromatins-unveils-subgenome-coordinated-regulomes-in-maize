# load data
args <- commandArgs(T)
true.acrs <- as.character(args[1])
fake.acrs <- as.character(args[2])
fdr <- as.numeric(args[3])
output <- as.character(args[4])

# import
a <- read.table(true.acrs)
b <- read.table(fake.acrs)

# estimate Tn5 density per 1000bp
a$tn5_density <- a[,ncol(a)]/((a$V3-a$V2)/1000)
b$tn5_density <- b[,ncol(b)]/((b$V3-b$V2)/1000)

# find emprical threshold
threshold <- quantile(b$tn5_density, 1-fdr)

# filter ACRs
a.filtered <- subset(a, a$tn5_density > threshold)
a.filtered$tn5_density <- NULL
a.filtered[,c(1:(ncol(a.filtered)-1))]
write.table(a.filtered, file=output, quote=F, row.names=F, col.names=F, sep="\t")

##usage:    >Rscript filter_eFDR.R <input_ACRs.bed> <controls.bed> <FDR> <output_filename.bed>
