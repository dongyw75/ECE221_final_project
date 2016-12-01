# -----------------------------------------------------------------------------
# Functions to apply intensity dependent filter.
# Typically applied to differential expression output (edgeR, DEseq) where
# large fold changes are seen at low intensity genes.
#
# See SEQanswers, post #34 in http://seqanswers.com/forums/showthread.php?t=9998
#
## Example usage:
## Assuming a table of differentially expressed genes (logFC) and intensities (logCPM)
## has been generated with edgeR using topTags():
#
# detable<- topTags(lrt, n= nrow(d))$table
# 
# EXAMPLE
#source('~/svn_checkout/bioinformatics-misc/intensityFilter.R')
#set.seed(1234)
#logCPM<- sort(rnorm(n= 10000, mean= 5, sd= 1))
#set.seed(12345)
#logFC<- sapply(length(logCPM):1, function(i) rnorm(n= 1, mean= 0, sd= i^(1/1.5)))
#z<- localZ(logCPM, logFC, nbins= 10)
#
#smoothScatter(x= logCPM, y= logFC, nrpoints= 1000)
#points(logCPM, logFC, col= ifelse(abs(z) > 1.5, 'red', NA), pch= 19, cex= 0.5)
#
# -----------------------------------------------------------------------------

z.score<- function(x){
    ## Return z-scores for vector x
    z<- (x - mean(x))/sd(x)
    return(z)
}

localZ<- function(x, y, nbins= 10){
    # Return z-score of each y for each window around x
    # x:
    #   Vector of intensities. Typically logCPM (x-axis in MA plot). 
    # y:
    #   Vector to be z-score'd. Typically logFC (y-axis in MA plot).
    # nbins:
    #   Divide the x vector into this many bins, each of them containing the
    #   same number of datapoints. z-scores are calculated within each bin.
    # -------------------------------------------------------------------------
    ## Number of datapoints surrounding the target point
    n<- round(length(x) / nbins, 0) 

    ## Record the order of input x (and y)
    xorder<- order(x)
    
    ## Sort y by x
    yorder<- y[xorder]

    ## Initialize vector of z-scores zvec
    zvec<- vector(length= length(x))
    ## For each point p in x,
    xleft<- floor(n/2)
    xright<- ceiling(n/2)
    for (i in 1:length(y)){
        if (i <= xleft){
            slice<- yorder[1:n]
        } else if ((i + xright) >= length(x)) {
            slice<- yorder[((length(x) - n)+1) : length(x)]
        } else {
            slice<- yorder[((i - xleft)+1) : (i + xright)]
        }
        zvec[i]<- (yorder[i] - mean(slice)) / sd(slice)
    }
    ## Return zvec in the original order of x (and y)
    zvec_ordered<- zvec[order(xorder)]
    return(zvec_ordered)
}
