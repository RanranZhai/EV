

load('/Users/ranran/Documents/GIANT/cancer120_proteins.RData')



fl <- list.files('/Users/ranran/PASCAL/output/GIANT')
fl <- fl[grep(fl, pattern = 'sum.genescores')]

traits <- unlist(strsplit(fl, '.sum.genescores.txt'))
names <- c('WC (2015)', 'WCadjBMI (2015)', 'WHR (2015)', 'WHRadjBMI (2015)', 'WHRadjBMI (Additive)', 'WHRadjBMI (Recessive)',
			'WAISTadjBMI (Active)', 'WAISTadjBMI (Inactive)', 'WAISTadjBMI SNPadjPA', 'WCadjBMI (NonSMK)', 'WCadjBMI (SMK)', 'WCadjBMI SNPadjSMK', 
			'WHRadjBMI (Active)', 'WHRadjBMI (Inactive)', 'WHRadjBMI SNPadjPA', 'WHRadjBMI (NonSMK)', 'WHRadjBMI (SMK)', 'WHRadjBMI SNPadjSMK', 
			'WHRadjBMI (GIANT+UKBB)')



giant <- matrix(NA, ncol = length(proteins), nrow = length(fl))

for (i in 1: 113) {
	for (j in 1:length(fl)){
		data <- read.table(paste0('/Users/ranran/PASCAL/output/GIANT/', fl[j]), header = T, stringsAsFactors=F)
		if(proteins[i] %in% data$gene_symbol){
			giant[j,i] <- data[data$gene_symbol == proteins[i],'pvalue']
		} else {
			giant[j,i] <- NA
		}
	}
	cat(i, ' ')
}

colnames(giant) <- proteins
rownames(giant) <- names

new_giant <- giant

new_giant[is.na(new_giant)] <- runif(sum(is.na(new_giant)))

fdr_giant <- matrix(p.adjust(new_giant, method = 'fdr'), ncol = ncol(new_giant), nrow = nrow(new_giant))

dimnames(fdr_giant) <- dimnames(new_giant)

giant_.05 <- fdr_giant[,colSums(fdr_giant < .05) > 0]

save(giant, giant_.05, new_giant, file = '/Users/ranran/Documents/GIANT/pascal_results.RData')





load('/Users/lanez/Documents/Projects/GIANT/pascal_results.RData')

require(corrplot)

nejm.colors <- function (n = 8, alpha = 1){
  pal <- c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1",
           "#6F99AD", "#FFDC91", "#EE4C97")
  acode <- substr(rgb(0, 0, 0, alpha = alpha), 8, 9)
  return(paste0(pal, acode)[1:n])
}


##giant_.05 <- giant_.05[, order(colSums(giant_.05 < .01), decreasing = T)]
#
#col <- colorRampPalette(c("#FFFFFF", "#0072B5"))
#
#
#pdf('/Users/lanez/Documents/Projects/GIANT/Figure1_20201006.pdf', 12,6)
#corrplot(-log10(new_giant[rownames(giant_.05), colnames(giant_.05)]), p.mat = giant_.05,
#	     is.corr = F, method = 'square', 
#	     addgrid.col = 'grey80', 
#	     #col = col(200),
#	     insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 1, pch.col = "white", 
#         tl.srt = 65, tl.col = 'black', 
#         cl.ratio = .1, cl.cex = 1)
#dev.off()






wc <- grep(rownames(giant_.05), pattern = 'WC')
wc <- c(wc, grep(rownames(giant_.05), pattern = 'WAIST'))

whr <- grep(rownames(giant_.05), pattern = 'WHR')


col <- colorRampPalette(c("#FFFFFF", "#FFFFFF", "#0072B5"))
pdf('/Users/lanez/Documents/Projects/GIANT/Figure1_20201008.pdf', 12,6)
corrplot(-log10(new_giant[rownames(giant_.05), colnames(giant_.05)][c(wc,whr),]), p.mat = giant_.05[c(wc,whr),],
	     is.corr = F, method = 'square', 
	     addgrid.col = 'grey80', 
	     col = col(200),
	     insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 1, pch.col = 'white', 
         tl.srt = 65, tl.col = 'black', 
         cl.ratio = .1, cl.cex = 1)
dev.off()









