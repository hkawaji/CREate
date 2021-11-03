# CREate 

__Cis-Regulatory Element indentification based on Architecture of Transcription initiation Events__

Cis-regulatory elements, in particular promoters and enhancers, have a feature
of being transcribed. The architecture of transcription in these regions are
predominantly divergent, sometime bidirectional and unidirectional, but not convergent.
Based on such transcription initiation profiles (such as CAGE), `CREate.sh` identifies
active genomic regions with non-convergent transcription as cis-regulatory elements.

An accompanied script, `selectCorrectUnencodedG2ctss.sh`, produces a data file 
of 5'-end frequencies (CTSS bed file) from CAGE read alignments with a
reference genome. Genuine transcription start site (TSS) is further enriched
by selecting reads with "G" at their 5'-ends that are not derived from the genome. 
Such "unencoded G" can be identified from their alignments in case that immediate
upstream of TSS is not G. In case that its immediate upstream is G, "unencoded G"
is estimated from their neighbors, and the frequency of starting sites are adjusted
if possible. Note that it is not applicable to CAGE data for HeliScope sequencing.
It is tested with [STAR](https://github.com/alexdobin/STAR) alignments of CAGE for
Illumina seuqneicng, and its use is highly recommended as a pre-processing of `CREate.sh`.

Author
------
KAWAJI, Hideya
