# dtCREs: Divergently Transcribed CREs (cis-regulatory elements)

Many class of cis-regulatory elements (such as promoters and enhancers) are transcribed
divergently. This repo contains two tools to process CAGE data, to identify divergently
transcribed regions regardless the level of directionality. 

* Count CAGE read 5'ends per CTSS by using only ones with unencoded "G" as possible (even if it matches to the genome), which is added by terminal transferase activity of the reverse transcriptase.
* Identify divergently transcribed regions, not convergent, through scanning genomic bases without TSS clustering nor peak identification.


Author
------
KAWAJI, Hideya
