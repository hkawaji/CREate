# CREate

__Cis-Regulatory Element identification based on Architecture of Transcription
Initiation Events__
CREate is a software tool designed to identify cis-regulatory elements by
analyzing transcription initiation activity measured by CAGE (Cap Analysis of
Gene Expression). It also predicts the target genes of the CREs, based on the
ABC (Activity-By-Contact) model, adapted to take the level of transcription
as the activity term.


Cis-regulatory elements, such as promoters and enhancers, are characterized by
active transcription and typically exhibit a divergent architecture. For
example, promoters that produce mRNA and long noncoding RNA often generate
upstream antisense RNA (uaRNA) or PROMoter uPstream Transcripts (PROMPTs),
while enhancers are frequently associated with bidirectional transcription,
although unidirectional patterns can also occur. Divergent transcription is
thus a common hallmark of both promoters and enhancers.

Based on transcription initiation profiles, CREate detects genomic regions
exhibiting non-convergent transcription and classifies them as candidate
cis-regulatory elements. Given a genomic position N and the transcription
initiation profiles in both the sense and antisense directions, N is classified
as divergent if the following condition is satisfied:


    [ left_rev + right_fwd ] > [ left_fwd + right_rev ]

where

                     left        right

      forward     (left_fwd)     (right_fwd)
               -------------[N]----------------
      reverse    (left_rev)     (right_rev)






Requirements

* gawk (GNU awk; CREate uses --assign and 3-argument match(), tested in v5.x)
* samtools (https://github.com/samtools/samtools , tested in version 1.9)
* bedtools (https://github.com/arq5x/bedtools2 , tested in v2.29.2)
* jksrc (http://hgdownload.cse.ucsc.edu/admin/ , tested in v425), which provides
  bedGraphToBigWig, bigWigAverageOverBed, bigWigMerge and bigWigToBedGraph
* R (https://www.r-project.org/ , tested in v4.0.2)
* tidyverse (https://www.tidyverse.org/ , tested in v1.3.0)
* (STAR, to generate the input alignment)



## Usage

    CREate.sh <subcommand> <args>



Subcommands

* bam2ctss : Generate a CTSS (CAGE tag starting site) count file from a BAM
* run      : Run the peak call steps at once (total, call, filter, classify,
             eachcount)
* total    : Generate a total CTSS count file from multiple CTSS BED files
* call     : Call divergently transcribed regions as candidate cis-regulatory
             elements
* filter   : Filter the regions
* classify : Classify the regions as PLA (promoter level activity) or ELA
             (enhancer level activity)
* eachcount: Count reads belonging to the divergently transcribed regions per
             sample
* makeprom : Build promoter regions of known genes (strand-aware)
* promact  : Quantify promoter activity (TPM, sense strand) from one CTSS file
* creact   : Quantify cis-element (CRE) activity (TPM) from one CTSS file
* hicprep  : Extract the contact values link needs from a Hi-C contact file
* link     : Predict regulatory interactions between CREs and known promoters
* version  : Print the CREate version


## Tests

    ./test/run_tests.sh

Runs in a few seconds against the fixtures in test/. Each check regenerates an
output and compares it to a committed reference, so any change in a result fails
until the reference is updated on purpose. Covered: total, call, filter,
classify, makeprom, promact, creact, hicprep and link, with and without measured
contact, plus the documented defaults.

Fixtures are named by where they come from. EGR1_* are real, from the EGR1 locus
reads in this directory, and each is what the subcommand that names it produces.
synth_* are generated, because two things cannot be checked with real data
alone:

* synth_contacts.bed.gz comes from a known power law (gamma = 1.0, scale = 5.0,
  5kb bins, one bin in seven omitted), so hicprep is checked for recovering
  those numbers rather than for reproducing its own past output. Dropping
  whole bins leaves the decay intact and lowers the intercept predictably.
* synth_promact.bed.gz adds two promoters to the real one. A single gene never
  exercises what ABC does, since there is only one denominator; with three the
  denominators differ, and one promoter lies outside every CRE, which is what
  puts an uncovered promoter into the candidate set.

synth_link.bed.gz and synth_link_hic.bed.gz are link run on those, so they carry
the synth_ prefix too. synth_link_Ocre.bed.gz is the same run under -O cre -B
0.2, the candidate rule and exponent that were the defaults before 2026-08; it
is kept so that results recorded under them stay reproducible, and a test holds
link to it.

## bam2ctss

Generate a CTSS (CAGE tag starting site) count file in BED format, from a BAM
alignment file

    CREate.sh bam2ctss -i infile.star.bam -g genome.fa \
        [-q mapQ (default:20)] \
        [-p parallel (default:8)] \
        [-w filter_window (default:20)] \
        [-b filter_target_base (default:G; otherwise N)] \
        [-r ratio_Gaddition (default:0.89)] \
        [-f full_set_of_intermediate_files (default: null)]

This subcommand counts capped RNA 5-ends at single-nucleotide resolution based
on CAGE (Cap Analysis of Gene Expression) read alignments. It selectively
considers reads that start with a genomically unencoded G (guanine), a hallmark
of high-confidence capped 5-ends (Science. 385(6704):eadd8394, 2024; Genome Res.
24(4):708-17, 2014). However, accurately identifying such unencoded Gs involves
technical challenges, especially in regions where the genome itself also
contains Gs near transcription start sites.

Determining whether a read carries an unencoded G can be difficult when the
genomic sequence immediately upstream (one nucleotide before the true TSS) also
harbors G. To resolve this ambiguity, the subcommand examines the surrounding
sequence within a defined window to determine whether such reads should be
retained or discarded.

Accurate identification of transcription start sites (TSSs) within G-stretch
regions is further complicated by ambiguity over whether the 5-end guanine is
genomically encoded or the result of post-transcriptional G-addition. The
subcommand applies a correction strategy by shifting a proportion of reads based
on an assumed G-addition rate.

Note: The input file must be a BAM file generated by STAR
(https://github.com/alexdobin/STAR) with the --alignEndsType Local option. This
alignment mode ensures that mismatched ends, particularly 5-ends, are treated as
soft clips through local alignment.


## run

Run the peak call steps (total, call, filter, classify, eachcount) at once

    CREate.sh run  \
        -i infiles(comma separated) \
        -c chrom_sizes \
        -o out_prefix \
        [-w window_size(default:200)] \
        [-l min_length] \
        [-n min_counts] \
        [-t min_tpm] \
        [-e PLA_or_ELA_threshold_in_TPM(default: 2)] \
        [-p parallel(default:20)] \
        [-z compress_decompress_program(default:gzip; zstd is recommended)]

Note that each of the input files should be a BED-formatted CTSS profile with
gzip (*.ctss.bed.gz). The resulting files are:

* out_prefix_param.txt
* out_prefix_{total|max}.ctss.{bed.gz|fwd.bw|rev.bw}
* out_prefix_each.ctss.bw.tar
* out_prefix_each.ctssTpm.bw.tar
* out_prefix_region.bed.gz
* out_prefix_region_filter_classify.bed.gz
* out_prefix_region_filter_classify_eachcounts.bed9pls3.gz
* out_prefix_region_filter_classify_eachcounts_totals.txt
* out_prefix_region_filter_classify_eachcounts.txt.gz


## total

Generate a total CTSS count file from multiple CTSS BED files

    CREate.sh total \
        -i infiles(comma separated) \
        -c chrom_sizes \
        -o out_prefix \
        [-p parallel(default:20)]

which produces the following files:

* PREFIX_total.ctss.{bed.gz|fwd.bw|rev.bw}
* PREFIX_max.ctss.{bed.gz|fwd.bw|rev.bw}


## call

Call divergently transcribed regions as candidate cis-regulatory elements

    CREate.sh call \
        -i infile.ctss.bed.gz \
        -c chrom_sizes \
        [-w window_size(default:200)] \
        [-m outfileMax.ctss.bed.gz] \
        [-p parallel(default:20)] \
        [-z compress_decompress_program(default:gzip; zstd is recommended)] \
    | gzip -c > output.bed.gz


The output file is formatted in BED9, where the thickStart/thickEnd specify the
'core' region which is predominantly divergent ('#') and the start/end specify
the region where signals are considered ('+'). The 'core' regions do not overlap
each other, while the signal considered region may. The 5th (score) column
indicates CRE activity based on TPM (transcripts per million).

            ffffffffffff
       +++++########++++
       rrrrrrrrrrrrr

Forward signals are indicated in 'f' and reverse signals in 'r'.


## filter

Filter the regions

    gunzip -c output.bed.gz \
    | CREate.sh filter \
        [-D max_directionality(default:1)] \
        [-d min_directionality(default:0)] \
        [-C max_counts(default:-1)] \
        [-c min_counts(default:4)] \
        [-T max_tpm(default:-1)] \
        [-t min_tpm(default:0)] \
        [-e min_counts_in_each_strand(default:0)] \
        [-x min_ctssMax(default:0)] \
        [-y min_ctssMax_in_each_strand(default:0)] \
        [-L max_length(default:-1)] \
        [-l min_length(default:5)] \
        [-v] (for invert match, such as "grep -v")

Filtering based on expression (-t TPM) is recommended, where 0.05 and 0.5 are
reasonable for CAGE.


## classify

Classify the regions as PLA (promoter level activity) or ELA (enhancer level
activity)

    gunzip -c output.bed.gz \
    | CREate.sh classify \
        [-e PLA_or_ELA_threshold_in_TPM (default: 2)]


## eachcount

Count reads belonging to the divergently transcribed regions per sample

    CREate.sh eachcount \
        -i divergently_transcribed_region.bed.gz  (BED9 format) \
        -e infile_each.ctss.bw.tar (archive of bigWig files for each strand) \
        -o out_prefix \
        [-p parallel(default:20)]


## makeprom

Build promoter regions of known genes (strand-aware)

    CREate.sh makeprom \
        -g gene.bed12.gz (BED12 gene models) \
        -c chrom_sizes \
        [-w promoterWindowHalf(default:500)] \
    | gzip -c > promoter.bed.gz

The TSS of each gene (strand-aware) is extended by promoterWindowHalf on both
sides, and overlapping windows on the SAME strand are merged (opposite strands
are kept separate). Output is BED6: chrom, start, end, gene(s) (with '|'
replaced by ','), 0, strand. No activity is assigned here; use promact to
quantify.


## promact

Quantify promoter activity (TPM, sense strand) from a single-sample CTSS

    CREate.sh promact \
        -i promoter.bed.gz (output of makeprom; gzip or plain) \
        -r ctss.bed.gz (single-sample CTSS; gzip or plain) \
        [-p parallel(default:20)]

For each promoter the same-strand (sense) CAGE is summed over [start,end] and
normalized to TPM by the total of all CTSS counts (both strands) in the input.
The output is the input BED6 with column 5 replaced by the activity (TPM). The
input needs at least 6 columns, since the strand decides which strand is sense.


## creact

Quantify cis-element (CRE) activity (TPM) from a single-sample CTSS

    CREate.sh creact \
        -i cre (call/classify output, BED9 with core in thickStart/thickEnd) \
        -r ctss.bed.gz (single-sample CTSS; gzip or plain) \
        [-m mode(default:divergent; divergent or whole)] \
        [-p parallel(default:20)]

For each CRE the sense/antisense CAGE is summed over its divergent extents
(forward over [thickStart,end], reverse over [start,thickEnd]) and normalized to
TPM by the total of all CTSS counts (both strands) in the input. The output is
the input BED9 with column 5 replaced by the activity (TPM).

-m whole sums both strands over the whole feature [start,end] and ignores
thickStart/thickEnd, for regions defined outside CREate that carry no called
core. It therefore takes as few as 4 columns, whereas the default (-m divergent)
needs the 8 that carry thickStart/thickEnd. A 4-column input gains the activity
as a new column 5.


## hicprep

Extract the contact values required by 'link' from a Hi-C contact file

    CREate.sh hicprep \
        -i contacts.bed.gz (chrom, bin1, bin2, contact; e.g. ENCFF134PUN) \
        -k promoter.bed (output of makeprom; gzip or plain) \
        [-r resolution(default:auto, detected from the bin coordinates)] \
        [-w maxDistance(default:5000000)] \
        [-o distanceProfile.tsv] \
    | gzip -c > contacts.prom.tsv.gz

The Hi-C file must already be coverage-normalized (KR or SCALE).
As only a small chunk of the Hi-C data is needed to predict links between
CRE and promoters, this subcommand strips a large Hi-C dataset into a small
table containing pairs that touch a promoter bin, and records in the header
the additional information required downstream (so that link needs nothing
but the file itself):

    #resolution=5000
    #gamma=0.986946
    #scale=5.029419
    #fitDistances=1000
    #fitR2=0.996486
    #recordedFraction=0.8060

gamma and scale are the slope and intercept of contact against distance in
log-log space, fitted on this file. fitDistances and fitR2 report how many
distances entered the fit and how well it held. recordedFraction is the share
of possible bin pairs the file actually records, which link uses to decide
whether dropping the distance expectation (-N) is safe. With -o, the
distance profile is written as well: distance, sum_contact, n_observed,
n_expected (all bin pairs possible at that distance, recorded or not).


## link

Predict regulatory interactions between CREs and known promoters

    CREate.sh link \
        -k promoterActivity.bed.gz (output of promact) \
        -r creActivity.bed.gz (output of creact) \
        [-w windowSizeForLinking(default:1000000)] \
        [-s scoreType(default:ABC; ABC or ACI)] \
        [-S scoreCutoff(default:0.001 for ABC, 1 for ACI)] \
        [-p promTpmCutoffForLinking(default:1)] \
        [-d pseudoCpm(default:0.001)] \
        [-b creActivityExponent(default:0.2)] \
        [-B promoterActivityExponent(default:1.0)] \
        [-a powerLawAlpha(default:1.0)] \
        [-f distanceFloor(default:5000)] \
        [-i contacts.tsv.gz (output of hicprep); \
            estimated contact based on power-law is used if not specified] \
        [-N] (use the measured contact alone, without the distance expectation \
            included in the contacts.tsv.gz) \
        [-O categoryTypeIfOverlapping(default:promoter; promoter or cre)] \
        [-W windowSizeForDenominator(default:5000000)] \
        [-t linkTargetRegions.bed.gz(default: the same file to -r)]

For each pair of CRE and promoter, ACI (activity contact index) is calculated:

    ACI = (creActivity + pseudoCpm)^creActivityExponent x (dist_eff/100kb)^-alpha

where dist_eff = max(|CRE midpoint - promoter midpoint|, distanceFloor). ACI is
scaled so that a 1cpm CRE at 100kb equals 1. ABC is ACI normalized per target
promoter. The denominator sums ACI of neighboring promoters (including the
target promoter itself) and the CREs overlapping none of them. -O decides which
of the two is used when a promoter and a CRE overlap.

The output is BEDPE-like (10 columns): CRE (chrom, start, end), promoter
(chrom, start, end), name, score, CRE strand, promoter strand. Both scores
(ABC and ACI) are always kept in the name column. The name (7th column) is
a flat 'key:value|key:value|...' string in three namespaces:

  src_*  : source CRE fields inherited from its name (src_id, src_tpm, ...)
  tgt_id : target promoter/gene name
  link_* : link_dist(Mb; raw, before flooring), link_estCnt (uses dist_eff),
           link_expCre, link_expProm, link_ACI, link_ABC

With -i the contact term becomes

    contact = [ powerlaw(dist_eff) + measured ] / ( 2 * powerlaw(100kb) )

with the bin size, gamma and scale read from the header hicprep wrote, so the
power law is the one fitted to that same file. It is the expected contact at a
distance, so a pair with nothing recorded keeps its distance ranking instead of
collapsing, and a pair carrying the contact typical for its distance reads 1, as
it does without -i. -N drops that expectation, leaving

    contact = measured / powerlaw(100kb).

What decides whether to use it is coverage, not the kind of Hi-C: two thirds
of pairs go unrecorded in a single intact Hi-C experiment and need the floor,
while an average over many leaves only a few percent. hicprep reports
recordedFraction for exactly this decision, and link warns if -N is used below
0.9.


## Author

KAWAJI, Hideya <h.kawaji@gmail.com>


