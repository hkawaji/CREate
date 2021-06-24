# CREate 

__Identification of Cis-Regulatory Elements based on Architecture of Transcription initiation Events__

Cis-regulatory elements, in particular promoters and enhancers, has a feature
of being transcribed. The architecture of such transcription is predominantly
divergent, sometime unidirectional, but not convergent. Based on CAGE (Cap Analysis
of Gene Expression) profiles, `CREate.sh` identifies genomic regions with non-convergent
transcription, as cis-regulatory elements.

`selectCorrectUnencodedG2ctss.sh` produces 5'-end frequencies based on CAGE
read alignments with a reference genome through enriching 5'-capped ends.
It takes reads starting with G unencoded in the genome, and reads with G
matching to the genome in case its neighbrs harbor unencoded G. It use
is recommended as a pre-processor of `CREate.sh`.

Author
------
KAWAJI, Hideya
