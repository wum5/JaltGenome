setwd("/Users/mengwu/Documents/Research/JaltGenome/github")
library(genoPlotR)
library(ggplot2)
library(phytools)


### below script is plot the synteny graph of the SEUSS locus
## dataframe 1 show the physical locations of target genes on jaltomata assembly scaffolds
df1 <- data.frame(name=c("scf29960g00010", "scf29960g00020", "scf29960g00030", "scf29960g00040", 
                         "scf29960g00050", "scf29960g00060", "scf29960g00070", "scf29960g00080",
                         "scf29960g00090", "scf29960g00100", "scf29960g00110", 
                         "scf29960g00120", "scf29960g00130", "scf29960g00140",
                         "Jalsin_scf29960g00150 (seu1)", "seu2", "Jalsin_scf29960g00160 (seu3)", 
                         "seu4", "seu5", "Jalsin_scf29960g00170 (seu6)", "Jalsin_scf29960g00180 (seu7)", 
                         "Jalsin_scf31961g00030 (seu8)", "Jalsin_scf31961g00020 (seu9)", 
                         "Jalsin_scf31961g00010 (seu10)"),
                  start=c(11934, 19252, 19591, 20023, 21396, 49761, 49187, 55932,
                          61334, 76213, 83871, 127932, 134122, 153608, 
                          168404, 177503, 203913, 216631, 235948, 249109,
                          261959, 302429, 340522, 347719),
                  end=c(12860, 19563, 19957, 20300, 43968, 54800, 54394, 56161,
                        74667, 76931, 86396, 132387, 135314, 160912, 
                        173482, 183490, 209775, 222712, 238067, 255182,
                        263132, 307998, 344035, 353369),
                  strand=c(-1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
                  col=c(rep("black", 24)),
                  fill=c(rep("blue", 14), rep("red", 10)))
dna_seg1 <- dna_seg(df1)

## dataframe 2 show the relative locations of target genes on tomato chromosome 11 (start from 700000)
df2 <- data.frame(name=c("Solyc11g005900", "Solyc11g005910", "Solyc11g005920", "Solyc11g005930", 
                         "Solyc11g005940", "Solyc11g005950", "Solyc11g005960", "Solyc11g005970",
                         "Solyc11g005980", "Solyc11g005985", "Solyc11g005990"),
                  start=c(14133, 18690, 43300, 51381, 57836, 62434, 67365, 71804, 76841, 83400, 91177),
                  end=c(15047, 39675, 47156, 52853, 62293, 66224, 69644, 73765, 83348, 90570, 99767),
                  strand=c(-1, 1, 1, 1, -1, 1, -1, -1, 1, 1, -1),
                  col=c(rep("black", 11)), gene_type="arrows",
                  fill=c(rep("blue", 10), "red"))
dna_seg2 <- dna_seg(df2)
dna_segs <- list(dna_seg1, dna_seg2)

mycol <- rgb(162, 162, 162, max = 255, alpha = 75, names = "lightgrey")

## dataframe 3 show the synteny relationship between orthologous genes
df3 <- data.frame(start1=c(11934,21396,49761,61334,83871,127932,134122,153608,168404,177503,203913,216631,
                           235948,249109,261959,302429,340522,347719),
                  end1=c(12860,43968,54800,74667,86396,132387,135314,160912,173482,183490,209775,222712,
                         238067,255182,263132,307998,344035,353369),
                  start2=c(14133,18690,43300,51381,57836,62434,67365,83400,91177,91177,91177,91177,
                           91177,91177,91177,91177,91177,91177),
                  end2=c(15047,39675,47156,52853,62293,66224,69644,90570,99767,99767,99767,99767,
                         99767,99767,99767,99767,99767,99767),
                  col=c(rep(mycol, 18)))

comparisons <- list(comparison(df3))
## there are two scaffolds, here I splitted them by their positions (0-263938; 300000-394325)
xlims <- list(c( 0, 263938, 300000, 394325),
              c( 0, 110000))
 
pdf("synteny_plotA.pdf", width=5, height=1.5) 
plot_gene_map(dna_segs=dna_segs, xlims=xlims, comparisons=comparisons)
dev.off()



### Below the script is to plot the read depth across the investigated scaffolds
### using 500-bp sliding windows
data <- read.csv("seuss_depth.csv")
window <- c()
depth <- c()
n <- 0
depth_sum <- 0
win_size <- 500
total_length <- 394325

for (n in 1:total_length) {
  index = which(data$position==n)
  cur_depth  <- data[index, "read_depth"]
  if (length(cur_depth) == 0) cur_depth = 0

  if (n %% win_size == 0) {
    window <- c(window, (n-win_size))
    depth_avg <- depth_sum/win_size
    depth <- c(depth, depth_avg)
    print (n-win_size) 
    print (depth_sum/win_size)
    depth_sum <- cur_depth
  }
  else depth_sum <- depth_sum + cur_depth
}


## below script is to highlight the gene regions in the plot (target/other genes)
data2 <- data.frame(window=window, depth=depth)
write.csv(data2, "windows_read_depth.csv", row.names=FALSE)

pdf("synteny_plotB.pdf", width=5, height=3) 
par(mfrow=c(2,1))
data2$category <- "Background"

# label the target genes 
data2$category[ (data2$window+250)>168404 & (data2$window+250)<173482 ] <- "SEUSS"
data2$category[ (data2$window+250)>177503 & (data2$window+250)<183490 ] <- "SEUSS"
data2$category[ (data2$window+250)>203913 & (data2$window+250)<209775 ] <- "SEUSS"
data2$category[ (data2$window+250)>216631 & (data2$window+250)<222712 ] <- "SEUSS"
data2$category[ (data2$window+250)>235948 & (data2$window+250)<238067 ] <- "SEUSS"
data2$category[ (data2$window+250)>249109 & (data2$window+250)<255182 ] <- "SEUSS"
data2$category[ (data2$window+250)>261959 & (data2$window+250)<263132 ] <- "SEUSS"
data2$category[ (data2$window+250)>302429 & (data2$window+250)<307998 ] <- "SEUSS"
data2$category[ (data2$window+250)>340522 & (data2$window+250)<344035 ] <- "SEUSS"
data2$category[ (data2$window+250)>347719 & (data2$window+250)<353369 ] <- "SEUSS"

# label the other genes
data2$category[ (data2$window+250)>11934 & (data2$window+250)<12860 ] <- "OTHER"
data2$category[ (data2$window+250)>19252 & (data2$window+250)<19563 ] <- "OTHER"
data2$category[ (data2$window+250)>19591 & (data2$window+250)<19957 ] <- "OTHER"
data2$category[ (data2$window+250)>20023 & (data2$window+250)<20300 ] <- "OTHER"
data2$category[ (data2$window+250)>21396 & (data2$window+250)<43968 ] <- "OTHER"
data2$category[ (data2$window+250)>49761 & (data2$window+250)<54800 ] <- "OTHER"
data2$category[ (data2$window+250)>49187 & (data2$window+250)<54394 ] <- "OTHER"
data2$category[ (data2$window+250)>55932 & (data2$window+250)<56161 ] <- "OTHER"
data2$category[ (data2$window+250)>61334 & (data2$window+250)<74667 ] <- "OTHER"
data2$category[ (data2$window+250)>76213 & (data2$window+250)<76931 ] <- "OTHER"
data2$category[ (data2$window+250)>83871 & (data2$window+250)<86396 ] <- "OTHER"
data2$category[ (data2$window+250)>127932 & (data2$window+250)<132387 ] <- "OTHER"
data2$category[ (data2$window+250)>134122 & (data2$window+250)<135314 ] <- "OTHER"
data2$category[ (data2$window+250)>153608 & (data2$window+250)<160912 ] <- "OTHER"

Palette <- c('dark grey',"blue", 'red')
p <- ggplot(data2[which(data2$depth>0),]) + 
  geom_point(aes(x=window, y=depth, color=category), size=2, alpha=0.8)
p + scale_colour_manual(values=Palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(legend.title=element_blank()) + theme(legend.position="none") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank())

dev.off()


## plot the SEUSS gene tree 
pdf("synteny_plotC.pdf", width=10, height=10) 
raxml.tree <- read.tree("RAxML_bestTree.seuss")
plot(raxml.tree)
dev.off()

