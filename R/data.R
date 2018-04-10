#' HM27 probe ordering
#'
#' Data frame of 27578 rows each describing the chip address of a probe.
#' Column names include Probe ID, methylated allele address, unmethylated
#' allele address, probe design type, color channel and abbreviations.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#'
#' HM27 probe chip address
"HM27.ordering"

#' HM450 probe ordering
#'
#' Data frame of 485577 rows each describing the chip address of a probe.
#' Column names include Probe ID, methylated allele address, unmethylated
#' allele address, probe design type, color channel and abbreviations.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
#' HM450 probe chip address
"HM450.ordering"

#' EPIC probe ordering
#'
#' Data frame of 865918 rows each describing the chip address of a probe.
#' Column names include Probe ID, methylated allele address, unmethylated
#' allele address, probe design type, color channel and abbreviations.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
#' EPIC probe chip address
"EPIC.ordering"

#' Cytoband information of hg19
#'
#' Cytoband coordinates for genome build hg19. Data frame with columns
#' chrom, chromStart, chromEnd, name and geiStain.
#'
#' Source: UCSC genome browser
#' 
"cytoBand.hg19"

#' Cytoband information of hg38
#' 
#' Cytoband coordinates for genome build hg38. Data frame with columns
#' chrom, chromStart, chromEnd, name and geiStain.
#'
#' Source: UCSC genome browser
#' 
"cytoBand.hg38"

#' EPIC.controls
#' 
#' Data frame of 635 rows each describing chip address of a control probe,
#' its color channel and the control probe type in EPIC.
#' 
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"EPIC.controls"

#' HM27.controls
#' 
#' Data frame of 144 rows each describing chip address of a control probe,
#' its color channel and the control probe type in HM27.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"HM27.controls"

#' HM450.controls
#'
#' Data frame of 850 rows each describing chip address of a control probe,
#' its color channel and the control probe type in HM450.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"HM450.controls"

#' EPIC.female.clean.chrY.probes
#'
#' Vector of 314 probe IDs for chromosome Y probes in EPIC array excluding
#' pseudoautosomal regions. Pseudoautosomal regions was identified based on
#' cross-hybridization signal in TCGA normal samples.
#'
"EPIC.female.clean.chrY.probes"

#' M450.female.clean.chrY.probes
#'
#' Vector of 337 probe IDs for chromosome Y probes in HM450 array excluding
#' pseudoautosomal regions. Pseudoautosomal regions was identified based on
#' cross-hybridization signal in TCGA normal samples.
#' 
"HM450.female.clean.chrY.probes"

#' EPIC.female.xlinked.chrX.probes
#'
#' Vector of 3433 probe IDs for X-linked probes in EPIC array. X-linkage is
#' based on intermediate DNA methylation signal in TCGA normal female samples.
#' 
"EPIC.female.xlinked.chrX.probes"

#' HM450.female.xlinked.chrX.probes
#'
#' Vector of 3797 probe IDs for X-linked probes in HM450 array.  X-linkage is
#' based on intermediate DNA methylation signal in TCGA normal female samples.
#' 
"HM450.female.xlinked.chrX.probes"

#' EPIC.hg19.probe2chr
#'
#' Named vector that maps 866895 EPIC probe to chromosomes based on genome
#' build hg19.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"EPIC.hg19.probe2chr"

#' HM27.hg19.probe2chr
#'
#' Named vector that maps 27588 HM27 probe to chromosomes based on genome
#' build hg19.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"HM27.hg19.probe2chr"

#' HM450.hg19.probe2chr
#'
#' Named vector that maps 485577 HM450 probe to chromosomes based on genome
#' build hg19.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"HM450.hg19.probe2chr"

#' EPIC.mapped.probes.hg19
#'
#' GenomicRange object for 840750 CpG probes included in EPIC array based
#' on genome build hg19.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"EPIC.mapped.probes.hg19"

#' EPIC.mapped.probes.hg38
#'
#' GenomicRange object for 838881 CpG probes included in EPIC array based
#' on genome build hg38.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"EPIC.mapped.probes.hg38"

#' HM450.mapped.probes.hg19
#'
#' GenomicRange object for 467097 CpG probes included in EPIC array based
#' on genome build hg19.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"HM450.mapped.probes.hg19"

#' HM450.mapped.probes.hg38
#'
#' GenomicRange object for 466007 CpG probes included in EPIC array based
#' on genome build hg38.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"HM450.mapped.probes.hg38"

#' EPIC.mask
#'
#' Vector of 82108 probe IDs recommended for masking in EPIC array based on
#' Zhou et al. 2017 Nucleic Acids Research. MASK.general column of the
#' original annotation is used here.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"EPIC.mask"

#' HM27.mask
#'
#' Vector of 2857 probe IDs recommended for masking in HM27 array based on
#' Zhou et al. 2017 Nucleic Acids Research. MASK.general column of the
#' original annotation is used here.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"HM27.mask"

#' HM450.mask
#'
#' Vector of 50186 probe IDs recommended for masking in HM450 array based on
#' Zhou et al. 2017 Nucleic Acids Research. MASK.general column of the
#' original annotation is used here.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"HM450.mask"

#' HM450.mask.tcga
#'
#' Vector of 89512 probe IDs recommended for masking in HM450 array based on
#' TCGA legacy.
#'
#' Source: Cancer Genome Atlas Network 2013 New England Journal of Medicine
#' 
"HM450.mask.tcga"

#' EPIC.typeI.extC
#'
#' Vector of 46733 C-extension Type-I probes included in EPIC array.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"EPIC.typeI.extC"

#' EPIC.typeI.extT
#'
#' Vector of 15638 T-extension Type-I probes included in EPIC array.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"EPIC.typeI.extT"

#' HM450.typeI.extC
#'
#' Vector of 45427 C-extension Type-I probes included in HM450 array.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"HM450.typeI.extC"

#' HM450.typeI.extT
#'
#' Vector of 15148 T-extension Type-I probes included in HM450 array.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"HM450.typeI.extT"

#' ethnicity.ccs.probes
#'
#' Vector of 332 color-channel switching probes for ethnicity inference.
#' Details in Zhou et al. 2017 Nucleic Acids Research
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"ethnicity.ccs.probes"

#' ethnicity.rs.probes
#'
#' Vector of 59 built-in SNP probes for ethnicity inference.
#'
#' Source: http://zwdzwd.github.io/InfiniumAnnotation
#' 
"ethnicity.rs.probes"

#' ethnicity.model
#'
#' Random forest model trained for ethnicity inference. An object of class
#' "randomForest". Details in Zhou et al. 2017 Nucleic Acids Research
#'
#' See scripts/make-data.R for details
#'
"ethnicity.model"

#' sex.inference.model
#'
#' Random forest model trained for sex inference. An object of class
#' "randomForest"
#'
#' See scripts/make-data.R for details
#' 
"sex.inference.model"

#' hg19.chrominfo
#'
#' List of two GenomicRange objects. The first represents the genomic
#' coordinates for gap, the second for sequence lengths. The information
#' is intended for genome build hg19.
#'
#' See scripts/make-data.R for building details.
#' 
"hg19.chrominfo"

#' hg38.chrominfo
#'
#' List of two GenomicRange objects. The first represents the genomic
#' coordinates for gap, the second for sequence lengths. The information
#' is intended for genome build hg38.
#'
#' See scripts/make-data.R for building details.
#' 
"hg38.chrominfo"

#' UCSC.refGene.gene2txn.hg19
#'
#' List of 27129 elements each contains the transcript ID for a gene.
#' The data was adapted from UCSC refGene (2016) for genome build hg19.
#'
#' See scripts/make-data.R for building details.
#' 
"UCSC.refGene.gene2txn.hg19"

#' UCSC.refGene.gene2txn.hg38
#'
#' List of 27221 elements each contains the transcript ID for a gene.
#' The data was adapted from UCSC refGene (2016) for genome build hg38.
#'
#' See scripts/make-data.R for building details.
#' 
"UCSC.refGene.gene2txn.hg38"

#' UCSC.refGene.txn2gene.hg19
#'
#' List of 57956 elements each contains the gene names for a transcript.
#' The data was adapted from UCSC refGene (2016) for genome build hg19.
#' 
"UCSC.refGene.txn2gene.hg19"

#' UCSC.refGene.txn2gene.hg38
#'
#' List of 58056 elements each contains the gene names for a transcript.
#' The data was adapted from UCSC refGene (2016) for genome build hg38.
#'
#' See scripts/make-data.R for building details.
#' 
"UCSC.refGene.txn2gene.hg38"

#' UCSC.refGene.txns.hg19
#'
#' GRangesList object of length 62083. Each contains a transcript and
#' the coordinates of its exons. The data was adapted from UCSC refGene
#' (2016) for genome build hg19.
#'
#' See scripts/make-data.R for building details.
#' 
"UCSC.refGene.txns.hg19"

#' UCSC.refGene.txns.hg38
#'
#' GRangesList object of length 66553. Each contains a transcript and
#' the coordinates of its exons. The data was adapted from UCSC refGene
#' (2016) for genome build hg38.
#'
#' See scripts/make-data.R for building details.
#' 
"UCSC.refGene.txns.hg38"

#' EPIC.betas.leuko.whole
#'
#' Matrix of 485577 probes x 3 samples of whole blood DNA methylation
#' beta values. Data obtained from Reinus et al 2013 PLoS One.
#'
"EPIC.betas.leuko.whole"

#' HM27.betas.leuko.whole
#'
#' Matrix of 27578 probes x 2 samples of whole blood DNA methylation
#' beta values. Data obtained from TCGA pilot study.
#' 
"HM27.betas.leuko.whole"

#' HM450.betas.leuko.whole
#'
#' Matrix of 485577 probes x 3 samples of whole blood DNA methylation
#' beta values. Data obtained from Reinus et al 2013 PLoS One.
#' 
"HM450.betas.leuko.whole"

#' Age Predictor based on Horvath 2013
#'
#' Data frame of 353 rows each corresponding a probe. Columns include
#' probe IDs and regression coefficients.
#' The data is adapted from Horvath 2013 Genome Biology.
#' 
"agePredHorvath353"

#' cellref.buccal
#'
#' Vector of 299850 binarized DNA methylation status for buccal cells.
#' 
'cellref.buccal'

#' cellref.CD14Monocytes
#'
#' Vector of 351705 binarized DNA methylation status for CD14Monocytes cells.
#' 
'cellref.CD14Monocytes'

#' cellref.CD19B
#'
#' Vector of 342174 binarized DNA methylation status for CD19B cells.
#'
'cellref.CD19B'

#' cellref.CD4T
#'
#' Vector of 343658 binarized DNA methylation status for CD4T cells.
#'
'cellref.CD4T'

#' cellref.CD56NK
#'
#' Vector of 336647 binarized DNA methylation status for CD56NK cells.
#'
'cellref.CD56NK'

#' cellref.CD8T
#'
#' Vector of 332973 binarized DNA methylation status for CD8T cells.
#'
'cellref.CD8T'

#' cellref.eosinophil
#'
#' Vector of 342377 binarized DNA methylation status for eosinophil cells.
#'
'cellref.eosinophil'

#' cellref.granulocytes
#'
#' Vector of 347903 binarized DNA methylation status for granulocytes cells.
#'
'cellref.granulocytes'

#' cellref.hair
#' 
#' Vector of 317041 binarized DNA methylation status for hair cells.
#'
'cellref.hair'

#' cellref.liver
#' 
#' Vector of 292658 binarized DNA methylation status for liver cells.
#'
'cellref.liver'

#' cellref.muscle
#' 
#' Vector of 292577 binarized DNA methylation status for muscle cells.
#'
'cellref.muscle'

#' cellref.neutrophil
#' 
#' Vector of 343822 binarized DNA methylation status for neutrophil cells.
#'
'cellref.neutrophil'

#' cellref.omentum
#' 
#' Vector of 277709 binarized DNA methylation status for omentum cells.
#'
'cellref.omentum'

#' cellref.saliva
#' 
#' Vector of 323797 binarized DNA methylation status for saliva cells.
#'
'cellref.saliva'

#' cellref.scFat
#' 
#' Vector of 283869 binarized DNA methylation status for scFat cells.
#'
'cellref.scFat'

#' cellref.skin
#' 
#' Vector of 280657 binarized DNA methylation status for skin cells.
#'
'cellref.skin'

#' cellref.spleen
#' 
#' Vector of 302593 binarized DNA methylation status for spleen cells.
#'
'cellref.spleen'
