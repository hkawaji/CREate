# dtCREs: Divergently Transcribed CREs (cis-regulatory elements)

Some classes of cis-regulatory elements (in particular, promoters and enhancers) are characterized
with divergent transcription. This repo contains two scripts to process CAGE data, to identify divergently
transcribed regions unifromly regardless the level of directionality. 

* __selectCorrectUnencodedG2ctss.sh__  Count CAGE read 5'ends per CTSS by using only ones with unencoded "G" as possible (in case the 5'-end "G" matches to the genome, it make an educated guess from the neighboring region), which is added by terminal transferase activity of the reverse transcriptase.
* __dtCREs.sh__ Identify divergently transcribed regions, not convergent, through scanning genomic bases without TSS clustering nor peak identification.


Author
------
KAWAJI, Hideya
