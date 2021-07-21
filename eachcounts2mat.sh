#!/bin/sh

usage()
{
  cat <<EOF

  $0 -i infileBed -t infileTotals

EOF
  exit
}


### handle options
while getopts b:t: opt
do
  case ${opt} in
  b) infileBed=${OPTARG};;
  t) infileTotals=${OPTARG};;
  *) usage;;
  esac
done

if [ ! -n "${infileBed-}" ]; then usage; fi
if [ ! -n "${infileTotals-}" ]; then usage; fi


cat <<EOF | R --slave 
  library(tidyverse)

  buf = read_tsv("${infileTotals}", col_names=F, comment = "#")
  totals = unlist( buf[,2] )
  names(totals) = unlist( buf[,1] )

  bed = read_tsv("${infileBed}", col_names=F, comment = "#")
  colnames(bed)[c(4,10,11,12)] = c("00Annotation","counts","countsFwd","countsRev")
  counts = bed %>% 
    select("00Annotation", counts) %>%
    separate( counts, sep=",", into=names(totals) )

  counts = rbind(
    c( "01STAT:TOTAL", totals ),
    counts
  )

  write.table(counts, sep="\t", quote = F,row.names=F)
EOF
