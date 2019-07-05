library("limma")

#venn
intestine <- read.csv("intestinevenn.csv")
intestine
attach(intestine)

#w1
a <- (a >= 1)
b <- (b >= 1)
intw1 <- cbind(a, b)
intw1 <- vennCounts(intw1)
intw1
vennDiagram(intw1, include = "both", 
            names = c("Control W1", "5-ALA W1"), 
            cex=2,circle.col=c(3,2), counts.col = "blue")


#w4
c <- (c >= 1)
d <- (d >= 1)
intw4 <- cbind(c, d)
intw4 <- vennCounts(intw4)
intw4
vennDiagram(intw4, include = "both", 
            names = c("Control W4", "5-ALA W4"), 
            cex=2,circle.col=c(3,2), counts.col = "blue")


#s1w
e <- (e >= 1)
f <- (f >= 1)
ints1w <- cbind(e, f)
ints1w <- vennCounts(ints1w)
ints1w
vennDiagram(ints1w, include = "both", 
            names = c("Control s1W", "5-ALA s1W"), 
            cex=2,circle.col=c(3,2), counts.col = "blue")




intw1and4 <- cbind(a, b, c, d)
intw1and4 <- vennCounts(intw1and4)
intw1and4
vennDiagram(intw1and4, include = "both", 
            names = c("Control W1", "5-ALA W1", "Control W4", "5-ALA W4"), 
            cex=2,circle.col=c(3,2), counts.col = "blue")


