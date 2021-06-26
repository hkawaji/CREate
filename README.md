# CREate 

__Cis-Regulatory Element indentification based on Architecture of Transcription initiation Events__

Cis-regulatory elements, in particular promoters and enhancers, has a feature
of being transcribed. The architecture of transcription in these regions are
predominantly divergent, often unidirectional, but not convergent. Based on
transcription initiation profiles (such as CAGE), `CREate.sh` identifies active
genomic regions with non-convergent transcription as cis-regulatory elements.

An accompanied script, `selectCorrectUnencodedG2ctss.sh`, produces a data file 
of 5'-end frequencies (CTSS bed file) from CAGE read alignments with a
reference genome. Authentic transcription start site (TSS) is further enriched
by selecting reads with "G" at their 5'-ends that are not derived from the genome. 
Such "unencoded G" can be identified from their alignments unambiguously in case
that immediate upstream of TSS is not G, and can be estimated from their neighbors
otherwise. The frequency of starting sites are adjusted when such unencoded G
matched to the upstream of TSSs. Its use is recommended in a pre-processing of
`CREate.sh`.

Author
------
KAWAJI, Hideya
