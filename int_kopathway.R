#intw1
pathwayintw1 <- read.csv("pathwayintw1.csv", header = TRUE, row.names=1)
meta <- read.csv("kopathintw1meta.csv", header = TRUE, row.names=1)

all(colnames(orthologsintw1) %in% rownames(meta))
all(colnames(orthologsintw1) == rownames(meta))

pathwayintw1dds <- DESeqDataSetFromMatrix(countData = pathwayintw1, colData = meta, design = ~Dose)

pathwaydds1 <- DESeq(pathwayintw1dds)
resultpathwayw1 <- results(pathwaydds1)
write.csv(resultpathwayw1, file = "resultpathwayw1.csv")


#intw4
pathwayintw4_1 <- read.csv("pathwayintw4_1.csv", header = TRUE, row.names=1)
meta_1 <- read.csv("kopathintw1meta_1.csv", header = TRUE, row.names=1)

pathwayintw4dds_1 <- DESeqDataSetFromMatrix(countData = pathwayintw4_1, colData = meta_1, design = ~Dose)

pathwaydds4_1 <- DESeq(pathwayintw4dds_1)
resultpathwayw4_1 <- results(pathwaydds4_1)
write.csv(resultpathwayw4_1, file = "resultpathwayw4_1.csv")


#intw1s
pathwayintw1s <- read.csv("pathwayintw1s.csv", header = TRUE, row.names=1)


pathwayintw1sdds <- DESeqDataSetFromMatrix(countData = pathwayintw1s, colData = meta, design = ~Dose)

pathwaydds1s <- DESeq(pathwayintw1sdds)
resultpathwayw1s <- results(pathwaydds1s)
write.csv(resultpathwayw1s, file = "resultpathwayw1s.csv")