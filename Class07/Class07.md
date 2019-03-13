Class 7: Packages from CRAN and BioConductor
================

Functions Revisited
-------------------

Load (i.e. **source**) our rescale() function from last class.

``` r
source("http://tinyurl.com/rescale-R")
```

Test this function

``` r
rescale(1:5)
```

    ## [1] 0.00 0.25 0.50 0.75 1.00

``` r
#"Non-numeric argument to binary operator"
#rescale( c(1:5, "string"))
```

We want to make this function more robust to these types of errors.

``` r
#"Input x should be numeric"
#rescale2( c(1:5, "string"))
```

``` r
is.numeric(1:5)
```

    ## [1] TRUE

``` r
is.numeric("string")
```

    ## [1] FALSE

``` r
!is.numeric("string")
```

    ## [1] TRUE

``` r
#Define example x and y
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
#How many TRUEs (matches) are there?
sum(is.na(x) & is.na(y))
```

    ## [1] 1

``` r
#Where is vector does the TRUE lay?
which(is.na(x) & is.na(y))
```

    ## [1] 3

Now we take our working snippet and make a first function.

``` r
both_na <- function(x,y){
  #Check for NA elements in both input vectors
  sum(is.na(x) & is.na(y))
}
```

``` r
both_na(x,y)
```

    ## [1] 1

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
```

``` r
both_na(x, y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
#vectors not the same length
#both_na2(x, y2)
```

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

both_na3(x,y)
```

    ## Found 1 NA's at position(s):3

    ## $number
    ## [1] 1
    ## 
    ## $which
    ## [1] 3

Exploration of rentrez.

``` r
#install.packages("rentrez")
library(rentrez)
```

    ## Warning: package 'rentrez' was built under R version 3.4.4

``` r
entrez_dbs()
```

    ##  [1] "pubmed"          "protein"         "nuccore"        
    ##  [4] "ipg"             "nucleotide"      "nucgss"         
    ##  [7] "nucest"          "structure"       "sparcle"        
    ## [10] "genome"          "annotinfo"       "assembly"       
    ## [13] "bioproject"      "biosample"       "blastdbinfo"    
    ## [16] "books"           "cdd"             "clinvar"        
    ## [19] "clone"           "gap"             "gapplus"        
    ## [22] "grasp"           "dbvar"           "gene"           
    ## [25] "gds"             "geoprofiles"     "homologene"     
    ## [28] "medgen"          "mesh"            "ncbisearch"     
    ## [31] "nlmcatalog"      "omim"            "orgtrack"       
    ## [34] "pmc"             "popset"          "probe"          
    ## [37] "proteinclusters" "pcassay"         "biosystems"     
    ## [40] "pccompound"      "pcsubstance"     "seqannot"       
    ## [43] "snp"             "sra"             "taxonomy"       
    ## [46] "biocollections"  "unigene"         "gencoll"        
    ## [49] "gtr"

``` r
entrez_db_searchable("protein")
```

    ## Searchable fields for database 'protein'
    ##   ALL     All terms from all searchable fields 
    ##   UID     Unique number assigned to each sequence 
    ##   FILT    Limits the records 
    ##   WORD    Free text associated with record 
    ##   TITL    Words in definition line 
    ##   KYWD    Nonstandardized terms provided by submitter 
    ##   AUTH    Author(s) of publication 
    ##   JOUR    Journal abbreviation of publication 
    ##   VOL     Volume number of publication 
    ##   ISS     Issue number of publication 
    ##   PAGE    Page number(s) of publication 
    ##   ORGN    Scientific and common names of organism, and all higher levels of taxonomy 
    ##   ACCN    Accession number of sequence 
    ##   PACC    Does not include retired secondary accessions 
    ##   GENE    Name of gene associated with sequence 
    ##   PROT    Name of protein associated with sequence 
    ##   ECNO    EC number for enzyme or CAS registry number 
    ##   PDAT    Date sequence added to GenBank 
    ##   MDAT    Date of last update 
    ##   SUBS    CAS chemical name or MEDLINE Substance Name 
    ##   PROP    Classification by source qualifiers and molecule type 
    ##   SQID    String identifier for sequence 
    ##   GPRJ    BioProject 
    ##   SLEN    Length of sequence 
    ##   MLWT    Molecular Weight 
    ##   FKEY    Feature annotated on sequence 
    ##   PORG    Scientific and common names of primary organism, and all higher levels of taxonomy 
    ##   ASSM    Assembly 
    ##   DIV     Division 
    ##   STRN    Strain 
    ##   ISOL    Isolate 
    ##   CULT    Cultivar 
    ##   BRD     Breed

``` r
search_protein <- entrez_search(db = "protein", term = "(HNH endonuclease[protein name]) AND pseudomonas [ORGN]")
print(search_protein)
```

    ## Entrez search result with 9841 hits (object contains 20 IDs and no web_history object)
    ##  Search term (as translated):  HNH endonuclease[protein name] AND "Pseudomonas"[O ...

Now I'll plot the trend in how much the word *synergy* was used in Pubmed over 2010-2018.

``` r
search_year <- function(year, term){
    query <- paste(term, "AND (", year, "[PDAT])")
    entrez_search(db="pubmed", term=query, retmax=0)$count
}
#

year <- 2010:2018
papers <- sapply(year, search_year, term="synergy", USE.NAMES=FALSE)

plot(year, papers, type='b', main="Abuse of Synergy")
```

![](Class07_files/figure-markdown_github/unnamed-chunk-20-1.png)
