QCB_455_final
================
Sena Cetin
2024-11-13

## Figure 1 (f-g), Supplementary Figure 1 (d-e)

Import data of significant Foxp3 and Foxp1 ChIP-seq peaks in Treg cells
and Tconv cells

``` r
# foxp3 peaks
foxp3_peaks <- read.csv(
  "/Users/senacetin/Documents/qcb455final/foxp3_peaks.csv"
)
# foxp1 peaks in treg
foxp1_peaks_treg <- read.csv(
  "/Users/senacetin/Documents/qcb455final/foxp1_peaks_Treg.csv"
)
# foxp1 peaks in tconv
foxp1_peaks_tconv <- read.csv(
  "/Users/senacetin/Documents/qcb455final/foxp1_peaks_Tconv.csv"
)

#foxp1_peaks_treg[!is.na(foxp1_peaks_treg$deseq.treg.q)&foxp1_peaks_treg$deseq.treg.q < 0.01, ]

#foxp3_peaks[!is.na(foxp3_peaks$deseq.foxp1wt.q) & foxp3_peaks$deseq.foxp1wt.q < 0.01, ]
```

Convert data from csv files to gff format so they can be used to extract
the corresponding sequences from the mouse genome

``` r
# form new df with the required columns for gff format
foxp3 <- foxp3_peaks %>% 
  select(seqnames, start, end, score, strand, name) %>% 
  mutate(
    source = "source", # placeholder for source
    phase = ".", # placeholder for phase
    feature = "peak", 
    attributes = name
  )

# reorganize cols into correct order
foxp3 <- foxp3 %>% 
  select(seqnames, source, feature, start, end, score, strand, phase, attributes)

head(foxp3, n = 10)
```

    ##    seqnames source feature    start      end score strand phase      attributes
    ## 1      chr1 source    peak  6214224  6215056   558      *     .  foxp3_peak_119
    ## 2      chr1 source    peak  7088416  7089588   210      *     . foxp3_peak_150b
    ## 3      chr1 source    peak  7396710  7399082  1052      *     . foxp3_peak_181a
    ## 4      chr1 source    peak  9824044  9824400    87      *     .  foxp3_peak_294
    ## 5      chr1 source    peak 10232654 10233234   292      *     .  foxp3_peak_334
    ## 6      chr1 source    peak 13267215 13268307   808      *     .  foxp3_peak_431
    ## 7      chr1 source    peak 13371875 13374740   266      *     . foxp3_peak_443g
    ## 8      chr1 source    peak 13383111 13384004   292      *     . foxp3_peak_446b
    ## 9      chr1 source    peak 13588367 13590207   170      *     . foxp3_peak_473d
    ## 10     chr1 source    peak 13593142 13593932   438      *     .  foxp3_peak_474

Repeat for each data set

``` r
# form new df with the required columns for gff format
foxp1_tconv <- foxp1_peaks_tconv %>% 
  select(seqnames, start, end, score, strand, name) %>% 
  mutate(
    source = "source", # placeholder for source
    phase = ".", # placeholder for phase
    feature = "peak", 
    attributes = name
  )

# reorganize cols into correct order
foxp1_tconv <- foxp1_tconv %>% 
  select(seqnames, source, feature, start, end, score, strand, phase, attributes)

head(foxp1_tconv, n = 10)
```

    ##    seqnames source feature    start      end score strand phase      attributes
    ## 1      chr1 source    peak  9772758  9773068   125      *     .   foxp1_peak_82
    ## 2      chr1 source    peak 20969275 20969700   320      *     .  foxp1_peak_230
    ## 3      chr1 source    peak 30873439 30874160   401      *     .  foxp1_peak_387
    ## 4      chr1 source    peak 30949552 30950417   211      *     . foxp1_peak_391a
    ## 5      chr1 source    peak 36307389 36308248   325      *     . foxp1_peak_490a
    ## 6      chr1 source    peak 36511567 36511958   760      *     .  foxp1_peak_498
    ## 7      chr1 source    peak 36761501 36762076   191      *     .  foxp1_peak_507
    ## 8      chr1 source    peak 39367608 39368099   642      *     .  foxp1_peak_613
    ## 9      chr1 source    peak 39936302 39937006   100      *     . foxp1_peak_641b
    ## 10     chr1 source    peak 43933414 43934160   428      *     .  foxp1_peak_740

``` r
# form new df with the required columns for gff format
foxp1_treg <- foxp1_peaks_treg %>% 
  select(seqnames, start, end, score, strand, name) %>% 
  mutate(
    source = "source", # placeholder for source
    phase = ".", # placeholder for phase
    feature = "peak", 
    attributes = name
  )

# reorganize cols into correct order
foxp1_treg <- foxp1_treg %>% 
  select(seqnames, source, feature, start, end, score, strand, phase, attributes)

head(foxp1_treg, n = 10)
```

    ##    seqnames source feature    start      end score strand phase      attributes
    ## 1      chr1 source    peak  6214552  6215093   140      *     .   foxp1_peak_30
    ## 2      chr1 source    peak  9772758  9773068   125      *     .   foxp1_peak_82
    ## 3      chr1 source    peak  9847984  9848564   242      *     .   foxp1_peak_84
    ## 4      chr1 source    peak 13373658 13374573   106      *     . foxp1_peak_137a
    ## 5      chr1 source    peak 16619088 16619682   143      *     .  foxp1_peak_178
    ## 6      chr1 source    peak 16687997 16689009    94      *     . foxp1_peak_182b
    ## 7      chr1 source    peak 21079130 21079627    58      *     .  foxp1_peak_235
    ## 8      chr1 source    peak 21333198 21334102   116      *     . foxp1_peak_245a
    ## 9      chr1 source    peak 24678391 24678907   125      *     .  foxp1_peak_302
    ## 10     chr1 source    peak 30873439 30874160   401      *     .  foxp1_peak_387

Export gff files that can be inputted into “Extract Genomic DNA” tool on
Galaxy

``` r
# foxp3 peaks
write.table(
  foxp3, 
  file = "/Users/senacetin/Documents/foxp3_peaks.gff", 
  sep = "\t", 
  row.names = FALSE, 
  col.names = FALSE, 
  quote = FALSE
  )

# foxp1 tconv peaks
write.table(
  foxp1_tconv, 
  file = "/Users/senacetin/Documents/foxp1_tconv_peaks.gff", 
  sep = "\t", 
  row.names = FALSE, 
  col.names = FALSE, 
  quote = FALSE
  )

# foxp1 tconv peaks
write.table(
  foxp1_treg, 
  file = "/Users/senacetin/Documents/foxp1_treg_peaks.gff", 
  sep = "\t", 
  row.names = FALSE, 
  col.names = FALSE, 
  quote = FALSE
  )
```

In Galaxy, import gff files and use Genome Extract tool with the
following settings: - Interpret features when possible: No - Souce of
refernce genome: locally cached - Using reference genome: mm10 (as used
in paper) - output format: fasta

``` r
# import fasta files as DNA String Set 
library(Biostrings)
```

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     intersect, setdiff, union

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     second, second<-

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:lubridate':
    ## 
    ##     %within%

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'XVector'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
foxp1_tconv_seq <- readDNAStringSet("/Users/senacetin/Documents/Galaxy37-[foxp1_tconv_seq].fasta")
foxp1_treg_seq <- readDNAStringSet("/Users/senacetin/Documents/Galaxy38-[foxp1_treg_seq].fasta")
foxp3_seq <- readDNAStringSet("/Users/senacetin/Documents/Galaxy39-[foxp3_seq].fasta")

# extract cycle length and mutate sequences back to tables 
foxp1_peaks_tconv <- foxp1_peaks_tconv %>% 
  mutate(data.frame(
    length = width(foxp1_tconv_seq),
    sequence = as.character(foxp1_tconv_seq),
    stringsAsFactors = FALSE
    )
  )
foxp1_peaks_treg <- foxp1_peaks_treg %>% 
  mutate(data.frame(
    length = width(foxp1_treg_seq),
    sequence = as.character(foxp1_treg_seq),
    stringsAsFactors = FALSE
    )
  )
foxp3_peaks <- foxp3_peaks %>% 
  mutate(data.frame(
    length = width(foxp3_seq),
    sequence = as.character(foxp3_seq),
    stringsAsFactors = FALSE
    )
  )
```

In the paper, sequence motifs were identified in 150 bp windows around
peak summits Although we do not have access to the data of all of the
summits in each peak, we have the position of the highest summit in each
peak.

For our MEME-ChIP analysis, we will sample 150 bp regions centered at
the position of the highest summit of each peak (confusingly called
‘peak’ in the data).

``` r
# x is peak data, size is desired region size
subset_sequence <- function(x, size) {
  center = as.integer(size / 2) # set center position of sub sequence
  sub_seqs <- c() # intialize list of sub-sequences
  
  for (i in 1:nrow(x)) {
    # if summit position is less than half of desired size
    if (x$peak[i] <= center) { 
      # subset sequence of desired size starting from 0 position
      sub_seq <- str_sub(x$sequence[i], 1, size) 
      
      # if summit position closer to end than half of desired size
    } else if ((length(x$sequence[i]) - x$peak[i]) >= center) { 
      
      # subset seq of desired size starting from end
      sub_seq <- str_sub(x$sequence[i], -size) 
      
    } else {
      # subset seq centering at summit position
      sub_seq <- str_sub(
        x$sequence[i], 
        (x$peak[i] - center + 1), 
        (x$peak[i] + center)) 
    }
    sub_seqs <- c(sub_seqs, sub_seq) # concatenate sub-sequence to list
  }
  return(sub_seqs)
}
```

Apply function to each data set

``` r
foxp1_tconv_subseqs <- subset_sequence(foxp1_peaks_tconv, 150)

print(foxp1_tconv_subseqs[1:10])
```

    ##  [1] "TGAAGTATGCTATTTTATTAAATCTTTTGTGAGAATTTCATACACACACACTGTACTGTGTTTACAGCATGGAACTTTCTAACGCTAGCCTCCTGGAGGACTGTAAACATGAGCGCCCTGGCCTCTAGGGGGCGGGCTTTCCCTTTGAAT"
    ##  [2] "AAACAAAGAGGGGGAAATGACTTAAGAAACTGAAAAGGGCAGGTTTGGTTTTTGTTGTTGCTGTTGCTGTTGTTGTTGTTGTTGTTGTTGCTGCTGCTGCTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTATTTCTCAAGATC"
    ##  [3] "GAGCCCGCAGCCATCGGCTTGCGGCGGCGGCGGCGCGGGAGCAGGAGCGGGGCCGCGTCACGTGACGCCACTGTTTACCGTTCCGGGGGCTGGGGCAGCGTGCGCCGCGCCTGCGCGCTAGAGCCTGCGCGCAGCCGCCCGCGAGGCTCT"
    ##  [4] "AGGCGCCCCTGCGTCACAAGTCCGCCCCTTTATAGCCGGCGCGACGTGCGCACAGGCGGCAGCCAATCACGGCGCGGGACGCCTGCGGCTCGAGGCGGGGCCCGGGACGCAGACGTCACTCGGGGCAGCCGCCACCGGTGCGAACCGGTT"
    ##  [5] "CCTCGGGTCTAAGGGGGCGGAGTCTGCTCCGCTCTAGAAAATAGGATTCACCAATGGGCGCCCAGAGCGCTGGGGGGGGGCGGGTCCATGGGCGGGTCCGCGGGGGAGGAGGCAGTCGAGCGGCGCCGGTAGTCAGTCGCGCCGAGAGCG"
    ##  [6] "CGGAGCCCCACGCCCCTCCGGCTCGCCCCGCCCCGCGCCGCTTGTTTACTCCCCGGCGCAGCCTAGTTCTCCCCAGAGCCCGCCCCCGCCAGCCGCTTATTGGTCAAAGTGACGCGGAGGCCCGCGCCGATTGGCTGCAACAACAGAGGC"
    ##  [7] "GGCTGGTACTGCGTGACCCCGGAAATGTGCTCTGCCCCGCCAGGAGCATGTGGCCTGTGATTTTCCTTCCACACAGCATTTTGTTGCTTTGACGTCAACTTCAAACACGGCTCTCTTTCCTGTCTTCGCATCAAGGAAATAGGCTGAGCT"
    ##  [8] "CCTGGGCCCGCCCTGGGTCCGCCCCATGCCGCTCGCTCCGCCCCTGACGCGCGCTGTATTTACGCCCCAGCCTAACTCCTAACTTCCCTCGCCACACCTCTTCCTCCTCGCCTCCCCGTGACCCGGAAGTTGTACGGCTACGCGACTTTC"
    ##  [9] "GTGAGTGCTTCGTGGCGATTTTGTTGGAAAGATGCTGTCACCCTTTGTTTCTGTTCTGTCAGCCCCACCTCAGATGCGGTTACCGTTGTTTTGTTTTGTTTTGTTTTTGTTTTTTCCCAAGAAGAAAACTACTACAATGCCCTCTACTGC"
    ## [10] "ACAGGAGGGAGCTCATTTCTGTCACGTGTTTTGCTTTGGGGACTCTGGCGACAAAACAGCCGACAACCCACTTCCTGTTGACTAAACAGCACCTCACCCTGCCAGAAAGGCCGGGCCGCTAAATCTGAGACTTCCGCGCTCGACTTTACA"

``` r
foxp1_treg_subseqs <- subset_sequence(foxp1_peaks_treg, 150)

print(foxp1_treg_subseqs[1:10])
```

    ##  [1] "CCCCCGGCTGGCTGTTTCCCGGCGTGCTCTGCGGCGGAAAGGCCGAGTCGACAATAACAAACCCCACGGCGGCCGCGACCCAGCCCTGCCAAGCTCTCAGTGCCTCGGCCGGCGGACTCGGGTCCCCGCGCGGAGCCGAGGGGCCGGAGC"
    ##  [2] "TGAAGTATGCTATTTTATTAAATCTTTTGTGAGAATTTCATACACACACACTGTACTGTGTTTACAGCATGGAACTTTCTAACGCTAGCCTCCTGGAGGACTGTAAACATGAGCGCCCTGGCCTCTAGGGGGCGGGCTTTCCCTTTGAAT"
    ##  [3] "ACCTCCCTGGATGCGGTCTGCGGAAACCTTTGTCAGGCGCCAGGTGACTGGCGCGGGGAGGGTGGGGGTGGGGGGTGGGGAGGCGGCCACGGATGTGGGCGGGACGCCGATGCGGGCGGGGCGCCCTCAGGCAGGCGCCCTAGGGGCAGG"
    ##  [4] "CGTGTTAATAACTGCTGCCTGGTTGTTTATTTCAATGGAGATCCTCCCCCAGCTCCCTCCTCCTCCTCCCCCTCCTCCTCCGCGGCTCTGCACTTGCGGAGAGAGACACAGACCGCGCGAGCTCGCGAGCAGGGGAGGGGGCTCGGGCAC"
    ##  [5] "CAAGGCCAGCGAGGCCCGAGACGCAGGGGGCGGAGCCGCCGTGCTCGCGTCACGCCGCGGGGGCCCAACATCCCGGGCACGGGGGTGGAGTCGCCGATGGGAAGGGGCCCGCGTGTCCGATTGGCTAGCGCGCGGGTGCCCGGGGCGGGG"
    ##  [6] "GAGCTGCGGCGCCAACTATTTCGAGAGAAAGTTTCCACTGAACCTGAAGCTGGAGGGCTGGCAGGAACGCTGCTCCAGCAGTGAGCCCGCTGCGGTGTGGGGTGGGGGTGGAGCAATGCCGACCACTTCCGGGTGCGGATCATCCTTTCC"
    ##  [7] "AACCGCCCCCGGCCACAGCCCTCGCCCCGCAGCGGCTGCACAGGCGCGGGCCTCTAGCCACACCCCCTTCCTCCGTCCGCCCCGCCCCCAGGCCGCAAAGGCGCGCGCCTCTGGCCACACCCCCTTCCTCCATCCGCCCCGCCCTCAAGC"
    ##  [8] "TACCAATCAGAATTAACTGGGGGCAGGTGCGTAGAAGCTATGTGCAGACTCTCTCGTTGTAAACAATTTTGGGGTGGGGGTGGGGCGAAATTATCATCAAAATACAAGCAGCTACAGCATACAATACAAATAAATCTTAAGAAAATTAAC"
    ##  [9] "CCAGCAAGATCAAGTGGCACGCAGGCCGCGCGTGGGCGGGCAGCGTGGTGACGTCACGGTGCGTCGGGTGCGCCCCGCCTCGGCTCGGCTGGTCGCCAGGGGACGCTCGCGGGCCGGGTCTCCTTGGCAACCCTCCGCTCTGCCCCTCCG"
    ## [10] "GAGCCCGCAGCCATCGGCTTGCGGCGGCGGCGGCGCGGGAGCAGGAGCGGGGCCGCGTCACGTGACGCCACTGTTTACCGTTCCGGGGGCTGGGGCAGCGTGCGCCGCGCCTGCGCGCTAGAGCCTGCGCGCAGCCGCCCGCGAGGCTCT"

``` r
foxp3_subseqs <- subset_sequence(foxp3_peaks, 150)

print(foxp3_subseqs[1:10])
```

    ##  [1] "AGACACCTGAGTCACCCCTCCCCCGGCTGGCTGTTTCCCGGCGTGCTCTGCGGCGGAAAGGCCGAGTCGACAATAACAAACCCCACGGCGGCCGCGACCCAGCCCTGCCAAGCTCTCAGTGCCTCGGCCGGCGGACTCGGGTCCCCGCGC"
    ##  [2] "TCCCGGGTTGCGGGTCGGGCTCCGCCCTCCCGCGTCGGCGCTGGCGGACGCTGTGCCGGGCCAGTCGGCTGTCTCTGCGGCGCTCCGCCCGCCGCCGGTCCGCCGCTGCTTTGGCCTGGAGGGCTCCCGCCTACGCAGGCTCGTGACGCG"
    ##  [3] "GTGCGAGAGCAACTGTCTCAGGGTTAAATCCAGAGAAACCAAGTCTTGGGGAAAAACCAAACAAACAAACAAAAACCAAACAAACAAACAAAAAAACCAAAACCAGGCAGGTGGTAGTGGCGGACGCCTTTAATCCCAGCACTGGGGAAG"
    ##  [4] "TCTTTTGGGAAAAATATAAAAGCCTTTTGCAGTGTTTTTAAAAATATAGATCACAATGCATTGGTAAGTCATGAAATCGACTTAGGAGTTCTGACCTGAATTAAATAGGCCAAAGCAAAAATAATGGAAATTTGCTGTGTATTAAAAGGT"
    ##  [5] "CCGTGGGGACCAGCATCCGCCGCAGCCGCCGACTGGCCTCCATCACCGCCGACCTGCCCCGGACGCCATGTTGAGAGGAAACGACGCGTCACGCAGCGGCCGAACGGCTGCGGAGAGGAAGCGGCGGCGGCGGCGGTGCCAGCGCCAGCC"
    ##  [6] "CAACGCTGCTGCTGGTGCCTTCGCTCTTCACTGTCAGTTTCCTGTTTACATCATCACTTCCCCCCTTTGATAAGGCAAATTTTATAAGGCGTAAAGCATTTTTTTTCCCTGACCACAGATTTAACAGCAGTCTCCTCTGGGACCGACCTT"
    ##  [7] "TGTTCCCGTGTTAATAACTGCTGCCTGGTTGTTTATTTCAATGGAGATCCTCCCCCAGCTCCCTCCTCCTCCTCCCCCTCCTCCTCCGCGGCTCTGCACTTGCGGAGAGAGACACAGACCGCGCGAGCTCGCGAGCAGGGGAGGGGGCTC"
    ##  [8] "TGGCATTTTAAAAAATGAAAGAGTCGGGAGTTCCCCAATATGATGACGTGTGTTTATGAAGAGGCTGTAACAGGTGTCATACTAATAGTGAAACATTGCAACTTCCTTTTGGGGTTGGGAAAGAAAATTAGGGACTCTATGCCTAGTGTT"
    ##  [9] "AGCGGTGGCGAGCATGCGCACCAGGAAGAAACCGCGCGCTCGCGCCGTGGCCGCACCTGGTGGCCCACAGCGCCCCCTGTGGGCTGGAGGTCGCGGAACGCGAGTGGCGCGCGGAGGTTTCGCGCTGGGAGTTGGGAGTCGTTTGGGGCG"
    ## [10] "GGTCCTTGGTTGCTGGCATTGTTTGGATTGGTTTAGGAGATATGCGTTGCTGGAAGAAGTGTGTCAGGGAACAGGCTTTGAGATTGTACAGTCATCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTGTTTTGCATTT"

Export fasta files for the Foxp1 and Foxp3 bound summits that can be
analyzed for motifs in MEME-ChIP

``` r
# 4,159 sequences total for foxp1 (tconv and treg)
foxp1_summits <- DNAStringSet(c(foxp1_tconv_subseqs, foxp1_treg_subseqs))
names(foxp1_summits) <- c(
  paste(foxp1_peaks_tconv$name, "_conv", sep = ""), 
  paste(foxp1_peaks_treg$name, "_reg", sep = "")
  )

# 7,174 sequences total for foxp3 
foxp3_summits <- DNAStringSet(foxp3_subseqs)
names(foxp3_summits) <- foxp3_peaks$name

# Export as fasta files
writeXStringSet(foxp1_summits, "/Users/senacetin/Documents/foxp1_summits.fasta")
writeXStringSet(foxp3_summits, "/Users/senacetin/Documents/foxp3_summits.fasta")
```
