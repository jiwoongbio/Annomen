args <- commandArgs(TRUE)

table <- read.table(args[1], header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "")
height <- table[, 2]
names(height) <- table[, 1]

color <- 1:length(unique(substr(table[, 1], 2, 6))) + 1
names(color) <- unique(substr(table[, 1], 2, 6))
color <- color[substr(table[, 1], 2, 6)]

pdf(file = args[2], width = as.numeric(args[3]), height = as.numeric(args[4]))
if(length(args) == 5) {
	barplot(height, las = 2, col = color, main = args[5])
} else if(length(args) == 6) {
	barplot(height, las = 2, col = color, main = args[5], ylim = c(0, as.numeric(args[6])))
}
dev.off()
