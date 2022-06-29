Transcriptomics data pre-processing
================
Olivia Angelin-Bonnet
June 29, 2022

``` r
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             legend.position = "bottom",
             text = element_text(size = 16))
```

## Data path

``` r
counts_path <- here("data/RNAseq/")
summary_path <- here("data/RNAseq/Alignment_stats_PGSCv4.03.csv")

counts_files_list <- list.files(path = counts_path, pattern = "*out.tab$")
cat(length(counts_files_list), "counts files.")
```

    ## 100 counts files.

``` r
gtf_path <- here("data/PGSC_DM_V403_fixed_representative_genes.gtf")
```

## Reading gene read counts

``` r
counts_df = lapply(counts_files_list, function(f){
  sample_name = str_remove(f, "_S\\d+_STARReadsPerGene.out.tab")
  counts_raw = read.table(paste0(counts_path, f), sep = "\t", header = FALSE)
  
  ## Getting the gene read counts (column 4 is for second_read_strand)
  counts = counts_raw[5:nrow(counts_raw), c(1, 4)]
  colnames(counts) = c("gene_id", sample_name)
  
  return(counts)
}) %>% 
  reduce(full_join, by = "gene_id")

dim(counts_df)
```

    ## [1] 39028   101

``` r
write_csv(counts_df, here("transcriptomics_differential_expression/processed_data/rna_seq_read_counts_raw.csv"))
```

## Genes filtering

We will remove from the dataset any gene for which more more than 5% of
the samples have a raw read count \< 5.

``` r
genes_low_expr = rowSums(counts_df[, -1] < 5) / (ncol(counts_df) - 1)
names(genes_low_expr) = counts_df$gene_id

sum(genes_low_expr > 0.95)
```

    ## [1] 13865

This filtering leads to removing 13865 genes, i.e. 35.53% of all genes.

``` r
retained_genes = names(genes_low_expr)[genes_low_expr <= 0.95]

counts_df_filtered = counts_df %>% 
  filter(gene_id %in% retained_genes)

write_csv(counts_df_filtered,
          here("transcriptomics_differential_expression/processed_data/rna_seq_read_counts_raw_filtered.csv"))
```

## Conversion of genes ID to GO terms

In this section, we query the Ensembl database, through the BiomaRt
package, to assign GO terms to the genes in the dataset.

``` r
genes_list = counts_df_filtered$gene_id
```

``` r
## Select a database
ensembl_plant = useMart("plants_mart", host="plants.ensembl.org")
```

    ## Warning: Ensembl will soon enforce the use of https.
    ## Ensure the 'host' argument includes "https://"

``` r
## Select a dataset
sol_ensembl = useMart("plants_mart", host="plants.ensembl.org", dataset = "stuberosum_eg_gene")
```

    ## Warning: Ensembl will soon enforce the use of https.
    ## Ensure the 'host' argument includes "https://"

We can then query the dataset using a list of filters. In the
transcriptomic dataset, the genes are identified with an ID starting
with `PGSC` and with `DMG` in it (which is different from the proteins
which have `DMT` in their PGSC ID):

``` r
filters = listFilters(sol_ensembl)
filters[str_which(filters$description, "PGSC\\d+DMG"), ]
```

    ##         name                            description
    ## 42 pgsc_gene PGSC ID(s) [e.g. PGSC0003DMG400000001]

Note that using `ensembl_gene_id` as filter seems to yield exactly the
same results.

We can extract different attributes from the dataset. For example, we
can get the GO term accession of the queries, their GO name, or the KEGG
ID.The attribute definition_1006 was removed from the query because one
of the query’s formatting for this attribute is problematic and blocks
the search:

``` r
attributes = listAttributes(sol_ensembl)
attributes[str_which(attributes$description, "(GO)|(KEGG)"), ]
```

    ##                      name             description         page
    ## 27                  go_id       GO term accession feature_page
    ## 28              name_1006            GO term name feature_page
    ## 29        definition_1006      GO term definition feature_page
    ## 30        go_linkage_type   GO term evidence code feature_page
    ## 31         namespace_1003               GO domain feature_page
    ## 32   goslim_goa_accession GOSlim GOA Accession(s) feature_page
    ## 33 goslim_goa_description  GOSlim GOA Description feature_page

``` r
attrs = c("go_id", "name_1006", "namespace_1003")
```

Notes about the different GO attributes:

the GO term evidence code refers to the evidence (e.g. inferred from
experiment, publication, etc) used to associate a gene/transcript to a
particular GO term.

A GO term can refer to 3 domains (namespace_1003):

-   Molecular function: the molecular activites of the individual gene
    products;
-   Cellular Component: where the gene products are active;
-   Biological Process: the pathways and larger processes to which that
    gene product’s activity contributes.

We use the function getBM() to build the query. I perform the query in
several batches, otherwise it generally fails.

``` r
batches = seq(1, length(genes_list), by = 200)
if(batches[length(batches)] != length(genes_list)) batches = c(batches, length(genes_list))

go_annot = lapply(2:length(batches), function(i){
  getBM(filters = "pgsc_gene",
        values = genes_list[batches[i-1]:batches[i]],
        attributes = c("pgsc_gene",
                       attrs),
        mart = sol_ensembl)
}) %>% 
  reduce(bind_rows) 

go_annot = go_annot %>% 
  filter(namespace_1003 %in% c("biological_process", "molecular_function", "cellular_component")) %>% 
  mutate(GO_domain = case_when(namespace_1003 == "biological_process" ~ 'Biological process',
                               namespace_1003 == "cellular_component" ~ 'Cellular component',
                               namespace_1003 == "molecular_function" ~ 'Molecular function')) %>% 
  rename(gene_id = pgsc_gene)

write_csv(go_annot, here("transcriptomics_differential_expression/processed_data/go_annot.csv"))
```

Out of the 25163 genes in the filtered transcriptomics dataset, 19646 of
them, or 78.1%, are associated with at least one GO term.

For future GO enrichment analyses, we will gather the genes assigned to
each GO term in the annotation retrieved for the dataset. They will
constitute the background set in enrichment analyses.

``` r
go_terms = unique(go_annot$go_id)

batches = seq(1, length(go_terms), by = 200)
if(batches[length(batches)] != length(go_terms)) batches = c(batches, length(go_terms))

go_background = lapply(2:length(batches), function(i){
  getBM(filters = "go",
        values = go_terms[batches[i-1]:batches[i]],
        attributes = c("go_id", "namespace_1003", "name_1006", "pgsc_gene"),
        mart = sol_ensembl)
}) %>% 
  reduce(bind_rows) %>% 
  distinct() %>% 
  filter(go_id %in% go_terms) %>% 
  mutate(GO_domain = case_when(namespace_1003 == "biological_process" ~ 'Biological process',
                               namespace_1003 == "cellular_component" ~ 'Cellular component',
                               namespace_1003 == "molecular_function" ~ 'Molecular function')) %>% 
  rename(gene_id = pgsc_gene) %>% 
  arrange(go_id)

## Checking how many background genes are not in the filtered dataset
# sum(!(unique(go_background$gene_id) %in% genes_list))

write_csv(go_background, here("transcriptomics_differential_expression/processed_data/go_background.csv"))
```

## Retrieving the genomic position of the genes

``` r
annotation = makeTxDbFromGFF(gtf_path, format = c("gtf"))
```

    ## Import genomic features from the file as a GRanges object ... OK
    ## Prepare the 'metadata' data frame ... OK
    ## Make the TxDb object ... OK

``` r
annot_list = as.list(annotation)

genes_info_df = left_join(
  annot_list$transcripts,
  annot_list$genes,
  by = "tx_id"
) %>% 
  select(gene_id, tx_chrom, tx_start, tx_end) %>% 
  distinct() %>% 
  rename(chrom = tx_chrom) %>% 
  mutate(pos = ceiling((tx_start + tx_end) / 2),
         posmb = pos / 1e6)
  # select(gene_id, chrom, posmb)
```

## Querying gene description

We will also use the biomaRt package to query the Ensembl database in
order to obtain a name and description for the different genes.

``` r
## Run the query by batch
batches = seq(1, length(genes_list), by = 200)
if(batches[length(batches)] != length(genes_list)) batches = c(batches, length(genes_list))

descr_annot = lapply(2:length(batches), function(i){
  getBM(filters = "pgsc_gene",
        values = genes_list[batches[i-1]:batches[i]],
        attributes = c("pgsc_gene", "ensembl_peptide_id", "description"),
        mart = sol_ensembl)
}) %>% 
  purrr::reduce(bind_rows) %>% 
  mutate(description = str_remove(description, "\\[Source:.+$")) %>% 
  rename(gene = pgsc_gene,
         transcript = ensembl_peptide_id) %>% 
  group_by(gene, description) %>% 
  summarise(transcript = paste0(transcript, collapse = ", ")) %>% 
  filter(description != "")

write_csv(descr_annot, here("transcriptomics_differential_expression/processed_data/descr_annot_nopos.csv"))
```

Out of the 25163 genes in the filtered transcriptomics dataset, 25157 of
them, or 99.98% have a description.

``` r
all(descr_annot$gene %in% genes_info_df$gene_id)
```

    ## [1] TRUE

``` r
genes_info_df %>% 
  full_join(descr_annot, by = c("gene_id" = "gene")) %>% 
  write_csv(here("transcriptomics_differential_expression/processed_data/descr_annot.csv"))
```
