library(tidyverse)
library(readxl)
library(here)
library(janitor)
library(Gviz)
library(biomaRt)
library(patchwork)

#data(geneModels)





listMarts()
r=useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org")#host="grch37.ensembl.org")
#r=useMart("ENSEMBL_MART_FUNCGEN")

listDatasets(r)
d=useDataset(dataset = "hsapiens_gene_ensembl", mart=r)

listFilters(d) #d@filters with_refseq_mrna
listAttributes(d)
listFilterValues(d, filter = "chromosome_name")

#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#listAttributes(ensembl)
genes <- getBM(attributes = c("ensembl_gene_id",
                              "hgnc_symbol",
                              "chromosome_name",
                              "start_position",
                              "end_position",
                              "5_utr_start",
                              "5_utr_end",
                              "3_utr_start",
                              "3_utr_end",
                              "exon_chrom_start",
                              "exon_chrom_end"),
               values = "15", #"*", 
               filters = "chromosome_name", #listFilters(ensembl)
               mart = d)

genes <- getBM(attributes = c("ensembl_gene_id",
                              "hgnc_symbol"),
                              #"chromosome_name",
                              #"start_position",
                              #"end_position",
                              #"5_utr_start",
                              #"5_utr_end",
                              #"3_utr_start",
                              #"3_utr_end",
                              #"exon_chrom_start",
                              #"exon_chrom_end"),
               values = "15", #"*", 
               filters = "chromosome_name", #listFilters(ensembl)
               mart = d)

#convert to Granges
genes$chromosome_name = paste0("chr", genes$chromosome_name)
genesGr <- with(genes, GRanges(chromosome_name, IRanges(start_position+1, end_position), ensg=ensembl_gene_id))
#remove genes not on std chromosome
genesGr <- genesGr[genesGr@seqnames %in% chrs,]



file_path_to_data <- here::here("data/PWSAtypDeletions.xlsx")
d_tan <- read_excel(file_path_to_data, sheet = "Sheet2", range = "C2:L16")

d_lit <- read_excel(file_path_to_data, sheet = "Sheet2", range = "C23:L32") %>% clean_names() %>%
  mutate(y_axis_idx = nrow(.):1,
         y_axis_idx = ifelse(author == "Cao", 8, y_axis_idx)) #+
  #geom_line(start_of_)

# 1 or 2 common deletions for PWS
# list author names instead of case 1, 2, etc.
# present case at the top
# zoom in other area of interest
# add in a blow up section
# overlap enhancer / encode data

min_coord0 <- min(d_lit$liftover_hg19_start) #- round(min(d_lit$liftover_hg19_start) * percent_expansion, 0)  #- 1000 #+100000#
max_coord0 <- max(d_lit$liftover_hg19_end) #+ round(max(d_lit$liftover_hg19_end) * percent_expansion, 0)#+ 1000
range_coords <- max_coord0 - min_coord0
percent_expansion <- 1#0.30
min_coord <- min(d_lit$liftover_hg19_start) - round(range_coords * percent_expansion, 0)  #- 1000 #+100000#
max_coord <- max(d_lit$liftover_hg19_end) + round(range_coords * percent_expansion, 0)#+ 1000

chromos <- "chr15"

cren_start <- d_lit %>% filter(author == "Crenshaw") %>% pull(liftover_hg19_start)
cren_end <- d_lit %>% filter(author == "Crenshaw") %>% pull(liftover_hg19_end)

p1 <- ggplot(d_lit, aes(x = liftover_hg19_start, 
                  xend = liftover_hg19_end,
                  #y= y_axis_idx-0.2, yend = y_axis_idx+0.2)) + 
                  y= y_axis_idx, 
                  yend = y_axis_idx)) + 
  geom_segment(size=2) + 
  geom_text(aes(label=paste0(author," ", date)), vjust=1.35, hjust=-0.05) +
  theme_classic() +
 # xlim(min(d_lit$start_of_deletion_breakpoint), max(d_lit$end_of_deletion_breakpoint))
  xlim(min_coord, max_coord) +
  #xlim(c(min(d_lit$liftover_hg19_start) - 30000, max(d_lit$liftover_hg19_end))) +
  xlab("Genomic Coordinates (hg19)") + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank() 
  ) +
  geom_vline(xintercept=cren_start, linetype="dotted", color="red", size=0.5) +
  geom_vline(xintercept=cren_end, linetype="dotted", color="red", size=0.5)
p1
#now, overlay genomic tracks from hg19 over this region
#ggbio: an R package for extending the grammar of graphics for genomic data
#https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-8-r77
#https://compgenomr.github.io/book/visualizing-and-summarizing-genomic-intervals.html

#Look on table browser - download annotations from here, may be more up to date - not sure why the snoU13 is present in this region, looks to be elsewhere when looking in UCSC browser
#https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1427561773_6nyof963MfmORbcpyPsYmfDKqeh4
gene.track <- BiomartGeneRegionTrack(genome = "hg19",
                                     GRCh = 37,
                                     chromosome = chromos, 
                                     start = min_coord, 
                                     end = max_coord,
                                     name = "ENSEMBL")
                                     #filter = list(with_refseq_mrna=TRUE))#"ENSEMBL") #"REFSEQ" 
gene.track2 <- gene.track
table(gene.track@range$feature)
#gene.track2@range <- gene.track2@range[gene.track2@range$feature %in% c("protein_coding", "snoRNA", "lincRNA")]#== "protein_coding"]
test_range1 <- gene.track2@range[gene.track2@range$feature == "protein_coding"]
test_range2 <- gene.track2@range[gene.track2@range$feature %in% c("protein_coding", "snoRNA", "lincRNA", "utr5")]

#convert this into a ggplot2 object
tbl_obj_genes <- test_range2 %>% as_tibble()

tbl_features <- tbl_gene.track2$feature %>% unique()

tbl_gene.track2 <- gene.track2@range %>% as_tibble()

tbl_gene.track2_filt <- tbl_gene.track2 %>% filter(feature == "processed_transcript")

lapply(tbl_features, function(x) tbl_gene.track2 %>% filter(feature == x))
#lapply(tbl_features, function(x) tbl_gene.track2 %>% filter(feature == x) %>% pull(symbol) %>% table())


feat_to_remove <- c("pseudogene", "sense_intronic", "non_coding")
#convert "sense_intronic" -> "lincRNA"
#processed_transcript only has 2 genes: RP11-701H24.10 and SNHG14
tbl_gene.track3 <- tbl_gene.track2 %>% filter(!feature %in% feat_to_remove) %>% 
  mutate(feature = ifelse(substr(symbol,1,5) == "SNORD", "snoRNA", feature)) %>%
  mutate(ylab_dat=as.factor(feature) %>% as.numeric())%>%
  mutate(ylab_dat = ifelse(symbol == "SNRPN", 6, ylab_dat),
         ylab_dat = ifelse(symbol == "SNURF", 4, ylab_dat),
         feature = ifelse(feature == "processed_transcript", "lincRNA", feature),
         ylab_dat = ifelse(feature == "snoRNA", 5, ylab_dat),
         ylab_dat = ifelse(substr(symbol,1,4) == "RP11", 3, ylab_dat),
         symbol_display = ifelse(symbol == "SNHG14", "", symbol))

avg_start_snhg14 <- tbl_gene.track3 %>% filter(symbol == "SNHG14") %>% pull(start) %>% mean()
height_snhg14 <- tbl_gene.track3 %>% filter(symbol == "SNHG14") %>% head(1) %>% pull(ylab_dat) #%>% mean()


p3 <- ggplot(tbl_gene.track3, aes(xmin=start, xmax=end,ymin=ylab_dat-0.5, ymax=ylab_dat, color=feature, fill=feature)) + 
  geom_line(aes(x = start, y = ylab_dat-0.25, group=paste0(tbl_gene.track3$feature, tbl_gene.track3$symbol)), color="black", size=0.25) +
  geom_rect() + 
  geom_text(aes(x=(start+end) / 2, y=ylab_dat+0.05, label=symbol_display), angle=90, hjust=0, size=3) +
  geom_text(aes(x=avg_start_snhg14, y=height_snhg14+0.05,  label="SNHG14"), color="black", hjust=0, size=3) +
  #coord_fixed(ratio=100000) +
  xlim(min_coord, max_coord) + #+
  ylim(0,max(tbl_gene.track3$ylab_dat)+2) +
  theme_classic() +
  # xlim(min(d_lit$start_of_deletion_breakpoint), max(d_lit$end_of_deletion_breakpoint))
  #xlim(c(min(d_lit$liftover_hg19_start) - 30000, max(d_lit$liftover_hg19_end))) +
  #xlab("Genomic Coordinates (hg19)") + 
  theme(axis.text=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  geom_vline(xintercept=cren_start, linetype="dotted", color="red", size=0.5) +
  geom_vline(xintercept=cren_end, linetype="dotted", color="red", size=0.5) +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  scale_fill_brewer(palette = "Dark2", direction = -1) +
  guides(color = guide_legend(reverse=TRUE),
         fill = guide_legend(reverse=TRUE))

  

p3/p1
#does it overlap with the ranges

min_max_range <- GRanges(
  seqnames = chromos,
  ranges = IRanges(start = min_coord, end = max_coord))

intersect_granges <- intersect(test_range2, min_max_range, ignore.strand=TRUE)
intersect_granges <- subsetByOverlaps(min_max_range, test_range2, ignore.strand=TRUE)

#findoverlaps
#subsetbyoverlaps

test_range2 %>% as_tibble() %>% View()

test_range2$symbol %>% table()
gene.track2@range$symbol %>% table()
p2 <- plotTracks(gene.track2,from=min_coord,to=max_coord,chromosome=chromos)
p2

p1 / plotTracks(gene.track2,from=min_coord,to=max_coord,chromosome=chromos)

  cpgi.track=AnnotationTrack(cpgi.gr,
                           name = "CpG")

CpGiFile=filePath=system.file("extdata",
                              "CpGi.hg19.table.txt",
                              package="compGenomRData")
cpgi.gr=genomation::readGeneric(CpGiFile, 
                                chr = 1, start = 2, end = 3,header=TRUE, 
                                keep.all.metadata =TRUE,remove.unusual=TRUE )


ggplot(d_lit, aes(x = start_of_deletion_breakpoint, 
           xend = end_of_deletion_breakpoint,
           #y= y_axis_idx-0.2, yend = y_axis_idx+0.2)) + 
          y= y_axis_idx, yend = y_axis_idx)) + 
  geom_segment()

geom_rect(aes(xmin = 1, xmax = 3, ymin = 10, ymax = 15))


ens_version <- biomaRt::useEnsembl(biomart = 'genes', 
                                   GRCh = 37,
                                   dataset = 'hsapiens_gene_ensembl'
)


genes2 <- getBM(attributes = c("ensembl_gene_id",
                              "hgnc_symbol",
               "chromosome_name",
               "start_position",
               "end_position"),
               #"5_utr_start",
               #"5_utr_end"),
               #"3_utr_start",
               #"3_utr_end",
               #"exon_chrom_start",
               #"exon_chrom_end"),
               values = "15", #"*", 
               filters = "chromosome_name", #listFilters(ensembl)
               mart = d)

genes2_filt <- genes2 %>% filter(start_position >=min_coord & end_position <=max_coord)

p4 <- ggplot(tbl_gene.track3, aes(xmin=start, xmax=end,ymin=ylab_dat-0.5, ymax=ylab_dat, color=feature)) + 
  geom_line(aes(x = start, y = ylab_dat-0.25, group=paste0(tbl_gene.track3$feature, tbl_gene.track3$symbol)), color="black", size=0.25) +
  geom_rect() + 
  geom_text(aes(x=start, y=ylab_dat+0.05, label=symbol), angle=90, hjust=0, size=2) +
  #coord_fixed(ratio=100000) +
  xlim(min_coord, max_coord) + #+
  ylim(0,max(tbl_gene.track3$ylab_dat)+2) +
  theme_classic() +
  # xlim(min(d_lit$start_of_deletion_breakpoint), max(d_lit$end_of_deletion_breakpoint))
  #xlim(c(min(d_lit$liftover_hg19_start) - 30000, max(d_lit$liftover_hg19_end))) +
  #xlab("Genomic Coordinates (hg19)") + 
  theme(axis.text=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  geom_vline(xintercept=cren_start, linetype="dotted", color="red", size=0.5) +
  geom_vline(xintercept=cren_end, linetype="dotted", color="red", size=0.5) +
  scale_color_brewer(palette = "Dark2")





genes3 <- getBM(attributes = c("ensembl_gene_id",
                               "hgnc_symbol",
                               "chromosome_name",
                               "start_position",
                               "end_position",
                #"5_utr_start",
                #"5_utr_end"),
                #"3_utr_start",
                #"3_utr_end",
                #"exon_chrom_start",
                #"exon_chrom_end"),
                values = "15", #"*", 
                filters = "chromosome_name", #listFilters(ensembl)
                mart = ens_version))

gene.track3 <- BiomartGeneRegionTrack(biomart = ens_version,#genome = "hg19",
                                     #GRCh = 37,
                                     chromosome = chromos, 
                                     start = min_coord, 
                                     end = max_coord,
                                     name = "ENSEMBL")
