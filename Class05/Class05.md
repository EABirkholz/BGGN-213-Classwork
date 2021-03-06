Graphics and Plots
================
Erica Birkholz
January 2019

``` r
#Class 5 Graphics and Plots
```

This is some narrative text that I can style **bold** or *italic* and add links [webpages](https://rmarkdown.rstudio.com/articles_report_from_r_script.html)

``` r
#Section 2A: line plot

weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
plot(weight, pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="Age (months)", ylab="Weight (kg)", main="Baby Weight")
```

![](Class05_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
plot(weight, typ="o", pch=15, cex=1.5, lwd=2, ylim=c(2,10), col="magenta", xlab="Age (months)", ylab="Weight (kg)", main="Baby Weight")
```

![](Class05_files/figure-markdown_github/unnamed-chunk-2-2.png)

``` r
#Section 2B: Barplot

feat <- read.table("bimm143_05_rstats/feature_counts.txt", header=T, sep="\t")
par(mar=c(3, 11, 4, 2))
barplot(feat$Count, horiz=T, ylab="Features", names.arg = feat$Feature, main="Mouse Features", las=1, xlim=c(0,80000))
```

![](Class05_files/figure-markdown_github/unnamed-chunk-2-3.png)

``` r
#Section 2C: Histograms

tenthou <- c(rnorm(10000), rnorm(10000)+4)
par(mar=c(3, 5, 4, 2))
hist(tenthou, breaks=30)
```

![](Class05_files/figure-markdown_github/unnamed-chunk-2-4.png)

``` r
#Section 3A: Color
mf <- read.table("bimm143_05_rstats/male_female_counts.txt", header=T, sep="\t")
par(mar=c(6, 4, 4, 2))
barplot(mf$Count, col=c("blue2", "red2"), main="Male and Female", ylim=c(0,20), names.arg=mf$Sample, las=2)        
```

![](Class05_files/figure-markdown_github/unnamed-chunk-2-5.png)

``` r
#Section 3B: Color by Value

genes <- read.table("bimm143_05_rstats/up_down_expression.txt", header=T, sep="\t")
table(genes$State)
```

    ## 
    ##       down unchanging         up 
    ##         72       4997        127

``` r
palette(c("magenta", "cyan", "yellow"))
par(mar=c(5, 5, 4, 2))
plot(genes$Condition1, genes$Condition2, col=genes$State, xlab="Condition 1", ylab="Condition 2", main="Gene Expression")
```

![](Class05_files/figure-markdown_github/unnamed-chunk-2-6.png)
