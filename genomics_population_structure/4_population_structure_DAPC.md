Population structure estimation with DAPC
================
Olivia Angelin-Bonnet
June 14, 2022

``` r
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(legend.position = "bottom")
```

## Data path

``` r
adegenet_input = here("genomics_data_preprocessing/processed_data/adegenet_input.csv")
samples_info_file = here("data/samples_info.csv")
gtf_file = paste0("/home/oangelin/Documents/Potato_data/Reference/PGSC_DM_V403_fixed_representative_genes.gtf")
```

``` r
samples_info = read_csv(samples_info_file, col_types = cols())
```

## Creating the `adegenet` data object

``` r
data_raw = read_csv(adegenet_input, show_col_types = FALSE)

temp = data_raw %>% 
  select(-Chrom, -Position) %>% 
  pivot_longer(cols = -Marker,
               names_to = "Sample",
               values_to = "Dosage") %>% 
  pivot_wider(names_from = Marker,
              values_from = Dosage)
mat = as.matrix(temp[, -1])
rownames(mat) = temp$Sample

data = new('genlight', mat, n.cores = 3)

adegenet::chromosome(data) <- str_remove(colnames(mat), "_.+$")
adegenet::position(data)   <- as.numeric(str_extract(colnames(mat), "(?<=_)\\d+$"))
adegenet::locNames(data)   <- colnames(mat)

save(data, file = here("genomics_population_structure/processed_data/data_genlight.RData"))
```

``` r
load(here("genomics_population_structure/processed_data/data_genlight.RData"))
```

### Loading the samples information (and pedigree)

``` r
samples_info = samples_info %>% 
  filter(VCF_sample_RGID %in% indNames(data))

dim(samples_info)
```

    ## [1] 176   8

``` r
samples_fam = samples_info$Family
names(samples_fam) = samples_info$VCF_sample_RGID

## for the parents, the family will be "parent"
samples_fam[samples_info$role == "parent"] = "parent"

samples_fam = samples_fam[indNames(data)]
```

``` r
families = list()
for(fam in setdiff(unique(samples_fam), c("parent"))){
  families[[fam]] = list("progeny" = samples_info$VCF_sample_RGID[samples_info$Family == fam])
  
  temp = lapply(c("Mum", "Dad"), function(x){
    parent = unique(samples_info[[x]][samples_info$Family == fam])
    parent = str_remove(parent, "mum")
    parent = str_remove(parent, "dad")
    parent[str_detect(parent, "DolCiVita")] = "Dolcevita"
    parent[str_detect(parent, "RedRascal")] = "Red"
    parent[str_detect(parent, "1335_8")] = "Crop_52"
    if("Crop_52" %in% parent) parent = c(parent, "Crop52")
    return(samples_info$VCF_sample_RGID[str_detect(samples_info$VCF_sample_RGID, parent)])
  })
  
  families[[fam]]$Mum = temp[[1]]
  families[[fam]]$Dad = temp[[2]]
  families[[fam]]$parents = c(families[[fam]]$Mum, families[[fam]]$Dad)
  families[[fam]]$all = c(families[[fam]]$parents, families[[fam]]$progeny)
}
```

``` r
## Add other parent column for PCA plots
samples_info <- samples_info %>% 
  mutate(other_parent = case_when(role == "parent" ~ Sampl_ID,
                                  Dad %in% c("Crop20dad", "1335_8dad") ~ Mum,
                                  Mum %in% c("Crop20", "1335_8") ~ Dad),
         other_parent = case_when(str_detect(other_parent, "(1335_8)|(Crop20)") ~ "(Crop20/Crop52)",
                                  TRUE ~ other_parent))
```

### k-means clustering

The first step for the DAPC analysis is to estimate the optimal number
of groups in the population of samples. In order to speed up this
clustering step, the data is first transformed using PCA. The problem is
that the PCA implementation used in the adegenet package is very slow.
We used a faster function that was proposed by a `adegenet` user on
GitHub, that is very fast.

``` r
## Copied from https://github.com/thibautjombart/adegenet/pull/150

glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){

  if(!inherits(x, "genlight")) stop("x is not a genlight object")

  # keep the original mean / var code, as it's used further down
  # and has some NA checks..
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }

  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }

  # convert to full data, try to keep the NA handling as similar
  # to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)

  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }

  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)

  # all dot products at once using underlying BLAS
  # to support thousands of samples, this could be
  # replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights

  ## PERFORM THE ANALYSIS ##
  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]

  ## scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }

  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging

  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")


  ## GET LOADINGS ##
  ## need to decompose X^TDV into a sum of n matrices of dim p*r
  ## but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }

      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }

    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }


  ## FORMAT OUTPUT ##
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }

  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }

  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }

  res$call <- match.call()

  class(res) <- "glPca"

  return(res)
}
```

``` r
## Performing the PCA
pcaX = glPcaFast(data,
                 center = TRUE,
                 scale = FALSE,
                 nf=min(c(nInd(data), nLoc(data))), ## number of retained PCs
                 loadings=TRUE,
                 returnDotProd = FALSE)
```

``` r
tibble(pc = 1:length(pcaX$eig),
       cumVar = 100 * cumsum(pcaX$eig)/sum(pcaX$eig)) %>% 
  ggplot(aes(x = pc, y = cumVar)) +
  geom_col(fill = "cornflowerblue", colour = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Number of retained PCs",
       y = "Cumulative variance (%)",
       title = paste("Variance explained by PCA - ", length(pcaX$eig), "PCs computed"))
```

<img src="4_population_structure_DAPC_files/figure-gfm/unnamed-chunk-13-1.png" width="\linewidth" style="display: block; margin: auto;" />

``` r
set.seed(1)
grp_all_k = find.clusters(data,
                          max.n.clust=10,
                          glPca = pcaX,
                          n.pca = min(c(nInd(data), nLoc(data))),
                          choose.n.clust = FALSE,
                          criterion = "smoothNgoesup") ## recommended criterion 
```

``` r
tibble(K = names(grp_all_k$Kstat),
       stat = unname(grp_all_k$Kstat)) %>% 
  mutate(K = as.numeric(str_extract(K, "\\d+")),
         colour = K == 5) %>% 
  ggplot(aes(x = K, y = stat)) +
  geom_line(colour = "gray") +
  geom_point(aes(colour = colour), size = 3) +
  scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'dodgerblue'), guide = 'none') +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Number of clusters K",
       y = "BIC",
       title = 'Number of group estimation')
```

<img src="4_population_structure_DAPC_files/figure-gfm/unnamed-chunk-15-1.png" width="\linewidth" style="display: block; margin: auto;" />

The k-means clustering requires as input the number of clusters K. As we
want to estimate this number, the clustering is repeated for different
K, and we can compare the BIC score resulting from the clustering with
the different values of K. We want to choose the K that minimises the
BIC score, while keeping the results of the clustering biologically
relevant. Here we’ll retain 5 clusters.

``` r
set.seed(12)
grp = find.clusters(data,
                    max.n.clust=10,
                    glPca = pcaX,
                    n.pca = min(c(nInd(data), nLoc(data))),
                    choose.n.clust = FALSE,
                    n.clust = 5)

save(pcaX, grp_all_k, grp, file = here("genomics_population_structure/processed_data/clusters.RData"))
```

### Discriminant Analysis on the clusters

The last step of the DAPC analysis is to perform a discriminant analysis
on the principal components from the PCA, with each sampled being
labelled as member of one of the K clusters, from the k-means clustering
step. The discriminant analysis aims at finding linear combinations of
the variables (here PCs from the PCA) that maximises the ratio of
variance between clusters over the ratio of variance within clusters.
These linear combinations are called discriminant functions.

Contrary to the clustering step, we don’t want to keep too many of the
PCA principal components for the discriminant analysis, as it would lead
to overfitting. With the package `adegenet`, we can estimate the optimal
number of PCs to retain via with cross-validation.

#### Optimal number of PCs to retain - cross-validation

This cross-validation method consists in retaining a certain number of
PCs and performing DAPC on a subset of the samples (training set), and
assess the performance of the resulting model in assigning the rest of
the samples (validation set) to their clusters. This is repeated for
different numbers of PCs retained, in order to find the trade off
between retaining too few and too many PCs for the analysis.

``` r
set.seed(15)
xval = xvalDapc(tab(data, NA.method="mean"),
                grp$grp,
                n.pca.max = length(pcaX$eig),
                training.set = 0.9,
                result = "groupMean",
                center = TRUE,
                scale = FALSE,
                n.pca = NULL,
                n.rep = 30,
                xval.plot = TRUE,
                parallel = "multicore",
                ncpus = 3)
save(xval, file = here("genomics_population_structure/processed_data/xval.RData"))
```

``` r
smoothScatter(xval[[1]]$n.pca, 
              as.vector(xval[[1]]$success),
              nrpoints=Inf,
              pch=20,
              col=transp("black"),
              ylim=c(0,1),
              xlab="Number of PCA axes retained",
              ylab="Proportion of successful outcome prediction", 
              main="DAPC Cross-Validation")
abline(h=xval[[2]], lty=c(2,1,2))
```

<img src="4_population_structure_DAPC_files/figure-gfm/unnamed-chunk-19-1.png" width="\linewidth" style="display: block; margin: auto;" />

The plot above shows the proportion of successful outcome predictions
for an increasing number of PCs retained. The solid and dotted lines
represent the median and confidence interval of the proportion of
successful outcome obtained by random chance. We can see that as we
increase the number of PCs that are retained, the performance decreases,
indicating a problem of overfitting. We can also plot the RMSE (root
mean squared error) score for increasing number of PCs retained:

``` r
tibble(npc = as.numeric(names(xval$`Root Mean Squared Error by Number of PCs of PCA`)),
       RMSE = xval$`Root Mean Squared Error by Number of PCs of PCA`) %>% 
  ggplot(aes(x = npc, y = RMSE)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = as.numeric(names(xval$`Root Mean Squared Error by Number of PCs of PCA`))) + 
  labs(x = "Number of PCs retained",
       y = "RMSE")
```

<img src="4_population_structure_DAPC_files/figure-gfm/unnamed-chunk-20-1.png" width="\linewidth" style="display: block; margin: auto;" />

    ## Number of PCs Achieving Lowest MSE: 20

The number of PCs achieving the lowest MSE is 20. However this is a bit
low, so we would rather retain 60 PCs, which amounts to retaining 66.8%
of the variance from our dataset.

``` r
dapc_res = dapc(data, 
                pop = grp$grp, 
                n.pca = 60, 
                glPca = pcaX, 
                n.da = 4, 
                var.contrib = TRUE, 
                var.loadings = TRUE)

save(dapc_res, file = here("genomics_population_structure/processed_data/dapc_res.RData"))
```

``` r
dapc_res
```

    ##  #################################################
    ##  # Discriminant Analysis of Principal Components #
    ##  #################################################
    ## class: dapc
    ## $call: dapc.genlight(x = data, pop = grp$grp, n.pca = 60, n.da = 4, 
    ##     var.contrib = TRUE, var.loadings = TRUE, glPca = pcaX)
    ## 
    ## $n.pca: 60 first PCs of PCA used
    ## $n.da: 4 discriminant functions saved
    ## $var (proportion of conserved variance): 0.668
    ## 
    ## $eig (eigenvalues): 1182 743.7 512.2 241  vector    length content                   
    ## 1 $eig      4      eigenvalues               
    ## 2 $grp      176    prior group assignment    
    ## 3 $prior    5      prior group probabilities 
    ## 4 $assign   176    posterior group assignment
    ## 5 $pca.cent 72847  centring vector of PCA    
    ## 6 $pca.norm 72847  scaling vector of PCA     
    ## 7 $pca.eig  175    eigenvalues of PCA        
    ## 
    ##   data.frame    nrow  ncol content                                          
    ## 1 $tab          176   60   retained PCs of PCA                              
    ## 2 $means        5     60   group means                                      
    ## 3 $loadings     60    4    loadings of variables                            
    ## 4 $ind.coord    176   4    coordinates of individuals (principal components)
    ## 5 $grp.coord    5     4    coordinates of groups                            
    ## 6 $posterior    176   5    posterior membership probabilities               
    ## 7 $pca.loadings 72847 60   PCA loadings of original variables               
    ## 8 $var.contr    72847 4    contribution of original variables               
    ## 9 $var.load     72847 4    loadings of original variables

## Saving the results

We are going to use the results of the DAPC analysis in our GWAS
analysis. Here we are saving the samples posterior group membership
probabilities as well as their coordinates in the distriminant functions
space.

``` r
temp = round(dapc_res$posterior, 3)
colnames(temp) = paste0("Grp", colnames(temp))
towrite = as_tibble(temp) %>% 
  mutate(Sample = rownames(temp)) %>% 
  select(Sample, everything())

write_csv(towrite, here("genomics_population_structure/processed_data/dapc_grp_membership.csv"))
```
