#!/bin/bash

set -ue

### important for making 'sort' and 'join' to be compatible
export LC_ALL=C


### define error message
usage()
{
  cat <<EOF

usage: $0 <command> <args>

  #------------------------------------------------
  # run all processes at once, except for filtering
  #------------------------------------------------
  $0 run 
    -i infilesList
    -c chrom_sizes
    -o out_prefix
    [-w window_size(default:200)]
    [-p parallel(default:20)]
    [-z compress_decompress_program(default:gzip; zstd is recommended if available)]

  Note that each of the input files should be
  BED-formatted CTSS profiles with gzip (*.ctss.bed.gz),
  and the file 'infilesList' should have their paths.


  #---------------------------------
  # run individual processes each by each
  #---------------------------------

  # produce total (and max) CTSS counts by merging CTSS profiles
  $0 total
          -i infileList
          -c chrom_sizes
          -o out_prefix
          [-p parallel(default:20)]

    which produces the following files:
      * PREFIX_total.ctss.{bed.gz|fwd.bw|rev.bw}
      * PREFIX_max.ctss.{bed.gz|fwd.bw|rev.bw}


  # call (identify) divergent regions
  $0 call
          -i infile.ctss.bed.gz
          -c chrom_sizes
          [-w window_size(default:200)]
          [-m outfileMax.ctss.bed.gz]
          [-p parallel(default:20)]
          [-z compress_decompress_program(default:gzip; zstd is recommended if available)]
  | gzip -c > output.bed.gz


  # count reads per samples
  $0 eachcount
          -i divergently_transcribed_region.bed.gz  (BED9 format)
          -e infile_each.ctss.bw.tar (archive of bigWig files for each strand)
          -o out_prefix
          [-p parallel(default:20)]


  # filter the divergently transcribed regions
  gunzip -c output.bed.gz
  | $0 filter
          [-D max_directionality(default:1)]
          [-d min_directionality(default:0)]
          [-C max_counts(default:-1)]
          [-c min_counts(default:0)]
          [-T max_tpm(default:-1)]
          [-t min_tpm(default:0)]
          [-e min_counts_in_each_strand(default:0)]
          [-x min_ctssMax(default:0)]
          [-y min_ctssMax_in_each_strand(default:0)]
          [-W max_width(default:-1)]
          [-w min_width(default:0)]
          [-v] (for invert match, such as "grep -v")

  recommended parameters for both cis-regulatory elements of promoters and enhancers:
    * -c 4 -t 0.05  (for (NET-)CAGE data with G selection and correction )
    * -c 4 -t 0.5   (for (NET-)CAGE data without G selection and correction )


(version. 2021.6.6)


Overview
--------
Cis-regulatory elements, promoter and enhancer, are known
to be transcribed in both orientation (foward and reverse)
divergently. This nature is evident in particular for enhancer,
as candidate regions are often defined through finding bidirectionally
transcribed regions, while many promoters are known to have PROMPTs
(promoter upstream transcripts). 

In those analysis, identification of divergently transcribed
regions is the first step. The pioneering work to identify
transcribed enhancers (Andersson et al. 2014) takes TSS clusters
or peaks, and finds divergent pairs of them. However, existing
methods of TSS clustering (or peak calling), such as single-likage
clustering (Carninci et al. 2016; or 'mergeBed' in bedtools),
DPI (Forrest et al. 2014), Paraclu (Frith et al. 2008) are not
necessarily designed for the purpose.

This program is designed to identify divergently transcribed regions,
not convergent, through scanning genomic bases without TSS clustering
nor peak identification.

'convergency' is defined as below, and any regions with convergency
is less than 1 are divergent.


  convergency = ( left_fwd + right_rev ) / ( left_rev + right_fwd )

where

                 left        right

  foward     (left_fwd)     (right_fwd)
           -------------[N]----------------
  reverse    (left_rev)     (right_rev)


'filter' subcommand can be used for finding bidirectionally
transcribed regions.


Output
-------
BED9 format, where the thickStart/thickEnd specify 'core' region which is divergent('#').
The start/end positions indicate the covered region where the signals are aggregated('+').

        ffffffffffff
   +++++########++++
   rrrrrrrrrrrrr

forward signals in 'f' and reverse signals in 'r' are accumulated for counts of 
foward (F) and  reverse (R). The directionality is computed as (F - R)/(F + R).

Note that divergent 'core' regions do not overlap each other, while the covered
regions may overlap.



Requirements
-------------
* bedtools
* jksrc



Author
------
KAWAJI, Hideya <h.kawaji@gmail.com>

EOF
  exit 1;
}



#-------------------------------------
# functions
#-------------------------------------

ctssbed_to_bg ()
{
  local strand=$1

  awk --assign strand=$strand 'BEGIN{OFS="\t"}{
    if($6 == strand )
    {
      for (i=$2; i < $3; i++) 
      {
        print $1,i,i+1,$5
      }
    }
  }' \
  | sort $SORT_OPT_BED
}



bg_to_bed4 ()
{
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"|"$2"|"$3}'
}


ctss_density ()
{
  local infile_bed4=$1
  local infile_bw=$2
  local density_width=$3
  local outfile_bg=$4

  cat ${chrom_sizes} \
  | awk --assign density_width=$density_width 'BEGIN{OFS="\t"}{print $1, density_width, $2 - density_width}' \
  > ${infile_bed4}.tmp.scope

  intersectBed -u -wa -a ${infile_bed4} -b ${infile_bed4}.tmp.scope \
  | split -l ${splitLines} - ${infile_bed4}.tmp.non_boundary.split.

  ls ${infile_bed4}.tmp.non_boundary.split.* \
  | xargs -n 1 -P ${parallel} -L 1 -I % bigWigAverageOverBed -sampleAroundCenter=$density_width $infile_bw % %.OUT

  cat ${infile_bed4}.tmp.non_boundary.split.*.OUT \
  | cut -f 1,5 \
  | sed -e 's/|/\t/g' \
  | sort $SORT_OPT_BED \
  > ${outfile_bg}

  rm -f ${infile_bed4}.tmp.non_boundary.split.*
}



# --- prepare the following files ---
# ${tmpdir}/infile.fwd.bg
# ${tmpdir}/infile.rev.bg
# ${tmpdir}/infile.fwd.bw
# ${tmpdir}/infile.rev.bw
# ${tmpdir}/infile.density.fwd.bg
# ${tmpdir}/infile.density.rev.bg
# ${tmpdir}/infile.density.fwd.bw
# ${tmpdir}/infile.density.rev.bw
# ${tmpdir}/potential_center.bg.COMP
prep_input ()
{
  local infile=$1
  local window_double_size=$2
  local density_width=20

  gunzip -c ${infile} | ctssbed_to_bg + > ${tmpdir}/infile.fwd.bg
  gunzip -c ${infile} | ctssbed_to_bg - > ${tmpdir}/infile.rev.bg
  bedGraphToBigWig ${tmpdir}/infile.fwd.bg $chrom_sizes ${tmpdir}/infile.fwd.bw
  bedGraphToBigWig ${tmpdir}/infile.rev.bg $chrom_sizes ${tmpdir}/infile.rev.bw

  cat ${tmpdir}/infile.fwd.bg | bg_to_bed4 > ${tmpdir}/infile.fwd.bed4
  cat ${tmpdir}/infile.rev.bg | bg_to_bed4 > ${tmpdir}/infile.rev.bed4

  ctss_density \
    ${tmpdir}/infile.fwd.bed4 \
    ${tmpdir}/infile.fwd.bw \
    ${density_width} \
    ${tmpdir}/infile.density.fwd.bg

  ctss_density \
    ${tmpdir}/infile.rev.bed4 \
    ${tmpdir}/infile.rev.bw \
    ${density_width} \
    ${tmpdir}/infile.density.rev.bg

  bedGraphToBigWig ${tmpdir}/infile.density.fwd.bg $chrom_sizes ${tmpdir}/infile.density.fwd.bw
  bedGraphToBigWig ${tmpdir}/infile.density.rev.bg $chrom_sizes ${tmpdir}/infile.density.rev.bw

  sort $SORT_OPT_BED ${tmpdir}/infile.fwd.bg ${tmpdir}/infile.rev.bg \
  | mergeBed  -i stdin -d $window_double_size \
  | awk 'BEGIN{OFS="\t"}{ for (i=$2;i<$3;i++){print $1,i,i+1,$4} }' \
  | $PROG_COMP \
  > ${tmpdir}/potential_center.bg.COMP
}



acc_left ()
{
  local infile_bw=$1
  local chrom_sizes=$2
  local window_size=$3
  local outfile=$4
  awk 'BEGIN{OFS="\t"}{name=$1"|"$2"|"$3; print $1,$2,$3,name}' \
  | slopBed -i - -g $chrom_sizes -l $window_size -r 0 \
  | bigWigAverageOverBed $infile_bw /dev/stdin $outfile
}



acc_right ()
{
  local infile_bw=$1
  local chrom_sizes=$2
  local window_size=$3
  local outfile=$4
  awk 'BEGIN{OFS="\t"}{name=$1"|"$2"|"$3; print $1,$2,$3,name}' \
  | slopBed -i - -g $chrom_sizes -l 0 -r $window_size \
  | bigWigAverageOverBed $infile_bw /dev/stdin $outfile
}



acc_lr ()
{
  local lr=$1
  local infile_bw=$2
  local chrom_sizes=$3
  local window_size=$4
  local outfile=$5

  if [ "${lr}" = "left" ]; then
    awk 'BEGIN{OFS="\t"}{name=$1"|"$2"|"$3; print $1,$2,$3,name}' \
    | slopBed -i - -g $chrom_sizes -l $window_size -r 0 \
    | bigWigAverageOverBed $infile_bw /dev/stdin $outfile
  elif [ "${lr}" = "right" ]; then
    awk 'BEGIN{OFS="\t"}{name=$1"|"$2"|"$3; print $1,$2,$3,name}' \
    | slopBed -i - -g $chrom_sizes -l 0 -r $window_size \
    | bigWigAverageOverBed $infile_bw /dev/stdin $outfile
  else
    printf "acc_lr: no operation" >&2
  fi
}




# accumulate signals for (left | right) x (forward | reverse)
accumulate_neiboring_signals () 
{
  local infile_center=$1
  local infile_fwd_bw=$2
  local infile_rev_bw=$3
  local window_size=$4

  declare -A infile_bw
  infile_bw["fwd"]=$infile_fwd_bw
  infile_bw["rev"]=$infile_rev_bw

  for lr in left right
  do
    for strand in fwd rev
    do

      cat $infile_center | $PROG_DECOMP > ${tmpdir}/accumulate_neiboring_signals_${lr}_${strand}.txt.DECOMP

      chunkN=$(expr $( cat ${tmpdir}/accumulate_neiboring_signals_${lr}_${strand}.txt.DECOMP | wc -l ) / $splitLines + 1 )
      split -n "r/${chunkN}" \
        ${tmpdir}/accumulate_neiboring_signals_${lr}_${strand}.txt.DECOMP  \
        ${tmpdir}/accumulate_neiboring_signals_${lr}_${strand}.txt.split.
      rm -f ${tmpdir}/accumulate_neiboring_signals_${lr}_${strand}.txt.DECOMP &

      export -f acc_lr
      infbw=${infile_bw[${strand}]}
      ls ${tmpdir}/accumulate_neiboring_signals_${lr}_${strand}.txt.split.* \
      | xargs -n 1 -P ${parallel} -L 1 -I % bash -c "cat % | acc_lr ${lr} $infbw $chrom_sizes $window_size %.OUT"

      cat ${tmpdir}/accumulate_neiboring_signals_${lr}_${strand}.txt.split.*.OUT \
      | cut -f 1,4 \
      > ${tmpdir}/accumulate_neiboring_signals_${lr}_${strand}.txt.split.OUT.CONCAT

      sort $SORT_OPT_NAME ${tmpdir}/accumulate_neiboring_signals_${lr}_${strand}.txt.split.OUT.CONCAT \
      | $PROG_COMP \
      > ${tmpdir}/accumulate_neiboring_signals_${lr}_${strand}.txt.COMP

      rm -f ${tmpdir}/accumulate_neiboring_signals_${lr}_${strand}.txt.split.* &
    done
  done

  wait

  cat ${tmpdir}/accumulate_neiboring_signals_left_fwd.txt.COMP \
  | $PROG_DECOMP \
  | join -t "	" /dev/stdin <( cat ${tmpdir}/accumulate_neiboring_signals_left_rev.txt.COMP  | $PROG_DECOMP ) \
  | join -t "	" /dev/stdin <( cat ${tmpdir}/accumulate_neiboring_signals_right_fwd.txt.COMP | $PROG_DECOMP ) \
  | join -t "	" /dev/stdin <( cat ${tmpdir}/accumulate_neiboring_signals_right_rev.txt.COMP | $PROG_DECOMP )

  #cat ${tmpdir}/accumulate_neiboring_signals_left_fwd.txt \
  #| join -t "	" /dev/stdin ${tmpdir}/accumulate_neiboring_signals_left_rev.txt \
  #| join -t "	" /dev/stdin ${tmpdir}/accumulate_neiboring_signals_right_fwd.txt \
  #| join -t "	" /dev/stdin ${tmpdir}/accumulate_neiboring_signals_right_rev.txt
}



classify_convergent_divergent ()
{
  local outfile_convergent=$1
  local outfile_divergent=$2

  awk \
    --assign outfile_convergent=$outfile_convergent \
    --assign outfile_divergent=$outfile_divergent \
  '{
    left_fwd = $2; left_rev = $3; right_fwd = $4; right_rev = $5;

    if ( (left_rev + right_fwd)  == 0 ) { next }
    convergency = (left_fwd + right_rev) / (left_rev + right_fwd) 
    if ( convergency > 1 ) { print > outfile_convergent".txt"; next}
    print > outfile_divergent".txt"
  }'

  cat ${outfile_convergent}".txt" \
  | sed -e 's/|/\t/g' | cut -f 1-3 | sort $SORT_OPT_BED \
  | mergeBed -i stdin \
  > ${outfile_convergent}

  cat ${outfile_divergent}".txt" \
  | sed -e 's/|/\t/g' | cut -f 1-3 | sort $SORT_OPT_BED \
  | mergeBed -i stdin \
  > ${outfile_divergent}
}



find_max_score_in_left () {
  local infile_bg=$1
  local extend_size=$2

  awk --assign extend_size=$extend_size 'BEGIN{OFS="\t"}{
    name=$1"|"$2"|"$3;
    start = $2 - extend_size;
    if (start < 0) {start = 0}
    print $1, start, $3, name
  }' \
  | sort $SORT_OPT_BED \
  | intersectBed -sorted -wa -wb -a stdin -b ${infile_bg} \
  | groupBy -grp 4 -opCols 5,6,7,8 -ops collapse \
  | awk 'BEGIN{OFS="\t"}{
      orig_n = split($1, orig_a, "|")
      chr_n = split($2, chr_a, ",")
      start_g = split($3, start_a, ",")
      stop_n = split($4, stop_a, ",")
      score_n = split($5, score_a, ",")

      chr = orig_a[1]; score = 0;
      start = orig_a[2]; stop = start + 1;
      leftmost_start = start;
      for (i=1; i<=score_n; i++)
      { 
        ### max score
        if ( score < score_a[i] )
        {
          chr = chr_a[i]; score = score_a[i];
          start = start_a[i]; stop = stop_a[i];
        } else if (( score == score_a[i]) && (start < start_a[i])) {
          chr = chr_a[i]; score = score_a[i];
          start = start_a[i]; stop = stop_a[i];
        }
        ### leftmost
        if (( score_a[i] > 0 ) && ( start_a[i] < leftmost_start ))
        {
          leftmost_start = start_a[i]
        }
      }
      print $1, chr"|"start"|"stop"|"score"|leftmost_start|"leftmost_start
    }'
}



find_max_score_in_right () {
  local infile_bg=$1
  local extend_size=$2

  awk --assign extend_size=$extend_size 'BEGIN{OFS="\t"}{
    name=$1"|"$2"|"$3;
    print $1,$2,$3 + extend_size,name
  }' \
  | sort $SORT_OPT_BED \
  | intersectBed -sorted -wa -wb -a stdin -b ${infile_bg} \
  | groupBy -grp 4 -opCols 5,6,7,8 -ops collapse \
  | awk 'BEGIN{OFS="\t"}{
      orig_n = split($1, orig_a, "|")
      chr_n = split($2, chr_a, ",")
      start_n = split($3, start_a, ",")
      stop_n = split($4, stop_a, ",")
      score_n = split($5, score_a, ",")

      chr = orig_a[1]; score = 0;
      stop = orig_a[3]; start = stop - 1;
      rightmost_stop = stop;
      for (i=1; i<=score_n; i++)
      {
        ### max score
        if ( score_a[i] > score )
        {
          chr = chr_a[i]; score = score_a[i];
          start = start_a[i]; stop = stop_a[i];
        } else if (( score == score_a[i]) && (start > start_a[i])) {
          chr = chr_a[i]; score = score_a[i];
          start = start_a[i]; stop = stop_a[i];
        }
        ### rightmost
        if (( score_a[i] > 0 ) && ( stop_a[i] > rightmost_stop ))
        {
          rightmost_stop = stop_a[i]
        }
      }
      print $1, chr"|"start"|"stop"|"score"|rightmost_stop|"rightmost_stop
  }'
}



find_max_scores () {
  local infile=$1
  local infile_density_fwd_bg=$2
  local infile_density_rev_bg=$3
  local window_size=$4

  cat $infile \
  | find_max_score_in_left ${infile_density_rev_bg} ${window_size} \
  | sort $SORT_OPT_NAME \
  > ${tmpdir}/find_max_scores_left.txt &

  cat $infile \
  | find_max_score_in_right ${infile_density_fwd_bg} ${window_size} \
  | sort $SORT_OPT_NAME \
  > ${tmpdir}/find_max_scores_right.txt &

  wait


  cat ${tmpdir}/find_max_scores_left.txt ${tmpdir}/find_max_scores_right.txt \
  | cut -f 1 \
  | sort $SORT_OPT_NAME \
  | uniq \
  | join -j 1 -t "	" -a 1 -o 1.1,2.2 -e "NA|NA|NA|NA|NA|NA" - ${tmpdir}/find_max_scores_left.txt \
  | join -j 1 -t "	" -a 1 -o 1.1,1.2,2.2 -e "NA|NA|NA|NA|NA|NA" - ${tmpdir}/find_max_scores_right.txt \
  > ${tmpdir}/find_max_scores.txt

  cat ${tmpdir}/find_max_scores.txt \
  | sed -e 's/|/\t/g' \
  | awk 'BEGIN{OFS="\t"}{
      chr = $1; start = $2; stop = $3;
      name = chr"|"start"|"stop;

      chrL = $4; startL = $5; stopL = $6;
      scoreL = $7; leftmost = $9;
      if  (chrL == "NA")
      {
        chrL = chr; startL = start; stopL = stop; 
        scoreL = 0; leftmost = start;
      }

      chrR = $10; startR = $11; stopR = $12;
      scoreR = $13; rightmost = $15;
      if ( chrR == "NA")
      {
        chrR = chr; startR = start; stopR = stop;
        scoreR = 0; rightmost = stop;
      }

      if ( ( leftmost < rightmost ) && (startL < stopR) )
      {
        print chr, leftmost, rightmost, name, 0, ".",
              startL, stopR,"0,0,0"
      } else {
        print chr, start,    stop,      name, 0, ".",
              start,      stop,     "0,0,0"
      }
    }' \
  > ${tmpdir}/find_max_scores.bed

  cat ${tmpdir}/find_max_scores.bed \
  | awk 'BEGIN{OFS="\t"}{print $1,$7,$8,$2","$3}' \
  | sort $SORT_OPT_BED \
  | mergeBed -i stdin -c 4 -o collapse \
  | awk 'BEGIN{OFS="\t"}{
      chr = $1
      pos_n = split($4,pos_a,",")
      start = pos_a[1]
      stop = pos_a[1]
      for (i = 2; i<= pos_n; i++)
      {
        if (start > pos_a[i] ) { start = pos_a[i] }
        if (stop < pos_a[i] ) { stop = pos_a[i] }
      }
      name = chr"|"start"|"stop
        print chr, start, stop, name, 0, ".", $2, $3, "0,0,0"
    }'
}


tighten_unidirectional ()
{
  local infile=$1

  cat $infile \
  | cut -f 1-4 \
  | sort $SORT_OPT_BED \
  | uniq \
  > ${tmpdir}/trim_unidirectional.outbound.bed

  bigWigAverageOverBed \
    ${tmpdir}/infile.fwd.bw \
    ${tmpdir}/trim_unidirectional.outbound.bed \
    /dev/stdout \
  | awk '{if($4 == 0){print $1}}' \
  | sed -e 's/|/\t/g' \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"|"$2"|"$3}' \
  | sort $SORT_OPT_BED \
  | intersectBed -sorted -wa -wb -a - -b ${tmpdir}/infile.rev.bg \
  | groupBy -g 1,2,3,4 -c 7 -o collapse \
  | cut -f 4,5 \
  | awk 'BEGIN{OFS="\t"}{
      stop_n = split($2,stop_a,",")
      rightmost = stop_a[1]
      for (i=2;i<=stop_n;i++)
      {
        if (rightmost < stop_a[i]) { rightmost = stop_a[i]}
      }
      print $1,"rightmost",rightmost
    }' \
  > ${tmpdir}/trim_unidirectional.outbound.adjust.txt

  bigWigAverageOverBed \
    ${tmpdir}/infile.rev.bw \
    ${tmpdir}/trim_unidirectional.outbound.bed \
    /dev/stdout \
  | awk '{if($4 == 0){print $1}}' \
  | sed -e 's/|/\t/g' \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"|"$2"|"$3}' \
  | sort $SORT_OPT_BED \
  | intersectBed -sorted -wa -wb -a - -b ${tmpdir}/infile.fwd.bg \
  | groupBy -g 1,2,3,4 -c 6 -o collapse \
  | cut -f 4,5 \
  | awk 'BEGIN{OFS="\t"}{
      start_n = split($2,start_a,",")
      leftmost = start_a[1]
      for (i=2;i<=stop_n;i++)
      {
        if (leftmost < start_a[i]) { leftmost = start_a[i]}
      }
      print $1,"leftmost",leftmost
    }' \
  >> ${tmpdir}/trim_unidirectional.outbound.adjust.txt

  cat ${tmpdir}/trim_unidirectional.outbound.adjust.txt \
  | sort $SORT_OPT_NAME \
  > ${tmpdir}/trim_unidirectional.outbound.adjust.txt.sort

  cat $infile \
  | sort $SORT_OPT_BED_NAME \
  | join -1 4 -2 1 -a 1 -t "	" - ${tmpdir}/trim_unidirectional.outbound.adjust.txt.sort \
  | awk 'BEGIN{OFS="\t"}{
      name = $1; chr = $2; start = $3; stop = $4;
      score = $5 ; strand = $6; thickStart = $7; thickStop = $8;
      color = $9 ; flag = $10; pos = $11;
      if (flag == "leftmost") { start = pos; thickStart = pos}
      if (flag == "rightmost") { stop = pos; thickStop = pos}
      name = chr"|"start"|"stop
      print chr, start, stop, name, score, strand, thickStart, thickStop, color
    }'
}


counts_fr ()
{
  local infile=$1
  local infileMax=
  if [ $# -ge 2 ]; then
    infileMax=$2
  fi

  cat $infile \
  | sed -e 's/\t/;/g' \
  | awk -F ";" 'BEGIN{OFS="\t"}{
      print $1, $7 ,$3, $0
  }' \
  | sort $SORT_OPT_BED \
  > ${tmpdir}/counts_fr.region_for_fwd.bed


  cat $infile \
  | sed -e 's/\t/;/g' \
  | awk -F ";" 'BEGIN{OFS="\t"}{
      print $1, $2 ,$8, $0
  }' \
  | sort $SORT_OPT_BED \
  > ${tmpdir}/counts_fr.region_for_rev.bed


  for strand in fwd rev
  do
    chunkN=$( expr  $( cat ${tmpdir}/counts_fr.region_for_${strand}.bed | wc -l ) / $splitLines + 1 )
    split -n "r/${chunkN}" ${tmpdir}/counts_fr.region_for_${strand}.bed ${tmpdir}/counts_fr.region_for_${strand}.bed.split.

    ls ${tmpdir}/counts_fr.region_for_${strand}.bed.split.* \
    | xargs -n 1 -P ${parallel} -L 1 -I % bigWigAverageOverBed ${tmpdir}/infile.${strand}.bw  % %.OUT

    cat ${tmpdir}/counts_fr.region_for_${strand}.bed.split.*.OUT \
    | cut -f 1,4 \
    | sort $SORT_OPT_NAME \
    >  ${tmpdir}/counts_fr.${strand}.txt

    rm -f cat ${tmpdir}/counts_fr.region_for_${strand}.bed.split.* &

    if [ -n "${infileMax}" ]; then
      gunzip -c ${infileMax} \
      | awk --assign strand=$strand '{ 
          if((strand =="fwd")&&($6 == "+")){print};
          if((strand =="rev")&&($6 == "-")){print}; 
        }' \
      | sort $SORT_OPT_BED \
      | intersectBed -sorted -wa -wb -a ${tmpdir}/counts_fr.region_for_${strand}.bed -b stdin \
      | sort $SORT_OPT_NAME \
      | groupBy -grp 4 -opCols 9 -ops max \
      | sort $SORT_OPT_NAME \
      > ${tmpdir}/max_fr.${strand}.txt &
    fi
  done

  wait

  #chr22   17141513        17142052        chr22,17141513,17142052 0       .       17141713        17142052        0,0,0   1       23
  join -t "	" ${tmpdir}/counts_fr.fwd.txt ${tmpdir}/counts_fr.rev.txt \
  > ${tmpdir}/counts_fr.joined.txt 

  if [ -n "${infileMax}" ]; then
    cat ${tmpdir}/counts_fr.joined.txt \
    | sort $SORT_OPT_NAME \
    | join -t "	" -a 1  -e 0 -o 0,1.2,1.3,2.2 /dev/stdin ${tmpdir}/max_fr.fwd.txt \
    | sort $SORT_OPT_NAME \
    | join -t "	" -a 1  -e 0 -o 0,1.2,1.3,1.4,2.2 /dev/stdin ${tmpdir}/max_fr.rev.txt \
    > ${tmpdir}/counts_fr.joined.txt.tmp
    mv --force ${tmpdir}/counts_fr.joined.txt.tmp ${tmpdir}/counts_fr.joined.txt
  fi


  cat ${tmpdir}/counts_fr.joined.txt \
  | sed -e 's/;/\t/g' \
  | awk 'BEGIN{OFS="\t"}{
      chr = $1; start = $2; stop = $3; name = $4; score = $5; strand = $6;
      thickStart = $7; thickStop = $8; color = $9; counts_fwd = $10; counts_rev = $11; max_fwd = $12 ; max_rev = $13 ;
      counts = counts_fwd + counts_rev
      if (counts == 0){next}
      directionality = (counts_fwd - counts_rev) / counts
      directionality_abs = directionality
      if (directionality_abs < 0) { directionality_abs = -1 * directionality_abs}
      name = name"|counts:"counts"|countsFwd:"counts_fwd"|countsRev:"counts_rev"|directionality:"directionality
      if ( max_fwd != "" ) {
        max_ctss = max_fwd;
        if ( max_ctss < max_rev ){max_ctss = max_rev}
        name = name"|ctssMax:"max_ctss"|ctssMaxFwd:"max_fwd"|ctssMaxRev:"max_rev
      }

      color = sprintf( "%d,0,%d", 127 + 127 * directionality , 127 - 127 * directionality )
      print chr, start, stop, name, directionality_abs, strand, thickStart, thickStop, color
    }'
}


counts_each ()
{
  local infileEach=$1
}



#--------------------------------------------------------
# main
#--------------------------------------------------------

cmd_total ()
{
  ### handle options
  parallel=20
  while getopts i:c:p:o: opt
  do
    case ${opt} in
    i) infilesList=${OPTARG};;
    c) chrom_sizes=${OPTARG};;
    p) parallel=${OPTARG};;
    o) out_prefix=${OPTARG};;
    *) usage;;
    esac
  done
  if [ ! -n "${infilesList-}" ]; then usage; fi
  if [ ! -n "${chrom_sizes-}" ]; then usage; fi
  if [ ! -n "${out_prefix-}" ]; then usage; fi

  ### setup for later
  tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
  trap "test -d $tmpdir && rm -rf $tmpdir" 0 1 2 3 15
  mkdir -p ${tmpdir}/each
  mkdir -p ${tmpdir}/merge
  SORT_OPT_BASE="--batch-size=100"
  export SORT_OPT_BED="${SORT_OPT_BASE} -k1,1 -k2,2n -k3,3n"
  export SORT_OPT_BED_NAME="${SORT_OPT_BASE} -k4,4"
  export SORT_OPT_NAME="${SORT_OPT_BASE} -k1,1"

  export -f ctssbed_to_bg
  cat $infilesList | xargs  -n 1 -P ${parallel} -I {} bash -c \
    "gunzip -c {} | ctssbed_to_bg + > ${tmpdir}/each/\$(basename {}).fwd.bg"
  cat $infilesList | xargs  -n 1 -P ${parallel} -I {} bash -c \
    "gunzip -c {} | ctssbed_to_bg - > ${tmpdir}/each/\$(basename {}).rev.bg"
  cat $infilesList | xargs  -n 1 -P ${parallel} -I {} bash -c \
    "bedGraphToBigWig ${tmpdir}/each/\$(basename {}).fwd.bg $chrom_sizes ${tmpdir}/each/\$(basename {}).fwd.bw"
  cat $infilesList | xargs  -n 1 -P ${parallel} -I {} bash -c \
    "bedGraphToBigWig ${tmpdir}/each/\$(basename {}).rev.bg $chrom_sizes ${tmpdir}/each/\$(basename {}).rev.bw"
 
  find  ${tmpdir}/each -name '*.fwd.bw' -type f -print > ${tmpdir}/merge/infileF.list
  find  ${tmpdir}/each -name '*.rev.bw' -type f -print > ${tmpdir}/merge/infileR.list

  bigWigMerge -inList ${tmpdir}/merge/infileF.list /dev/stdout \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",$4,"+"}' \
  > ${tmpdir}/merge/totalF.bed &

  bigWigMerge -inList ${tmpdir}/merge/infileR.list /dev/stdout \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",$4,"-"}' \
  > ${tmpdir}/merge/totalR.bed &

  bigWigMerge -max -inList ${tmpdir}/merge/infileF.list /dev/stdout \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",$4,"+"}' \
  > ${tmpdir}/merge/maxF.bed &

  bigWigMerge -max -inList ${tmpdir}/merge/infileR.list /dev/stdout \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",$4,"-"}' \
  > ${tmpdir}/merge/maxR.bed &

  wait

  for operation in total max
  do
    cat ${tmpdir}/merge/${operation}F.bed ${tmpdir}/merge/${operation}R.bed \
    | sort $SORT_OPT_BED \
    | awk 'BEGIN{OFS="\t"}
      {
        for (start=$2; start < $3; start++ ) 
        {
          print $1, start, start + 1, $4, $5, $6
        }
      }' \
    | gzip -c > ${tmpdir}/merge/${operation}.ctss.bed.gz

    gunzip -c ${tmpdir}/merge/${operation}.ctss.bed.gz \
    | ctssbed_to_bg + > ${tmpdir}/merge/${operation}.ctss.fwd.bg
    gunzip -c ${tmpdir}/merge/${operation}.ctss.bed.gz \
    | ctssbed_to_bg - > ${tmpdir}/merge/${operation}.ctss.rev.bg
    bedGraphToBigWig ${tmpdir}/merge/${operation}.ctss.fwd.bg $chrom_sizes ${tmpdir}/merge/${operation}.ctss.fwd.bw
    bedGraphToBigWig ${tmpdir}/merge/${operation}.ctss.rev.bg $chrom_sizes ${tmpdir}/merge/${operation}.ctss.rev.bw

    for suffix in .bed.gz .fwd.bw .rev.bw
    do
      mv -f ${tmpdir}/merge/${operation}.ctss${suffix} ${out_prefix}_${operation}.ctss${suffix}
    done
  done

  (cd ${tmpdir}/each/;  ls *.bw | sort | xargs tar cavf - ) > ${out_prefix}_each.ctss.bw.tar

  #mv ${tmpdir}/merge/total.ctss.bed.gz ${out_prefix}_total.ctss.bed.gz
  #mv ${tmpdir}/merge/max.ctss.bed.gz ${out_prefix}_max.ctss.bed.gz
  #(cd ${tmpdir}/each/;  ls *.bw | sort | xargs tar cavf - ) > ${out_prefix}_each.ctss.bw.tar
}


cmd_call ()
{
  ### handle options
  parallel=20
  window_size=200
  export splitLines=10000000
  infileMax=''
  infileEach=''
  prog_compression=gzip
  while getopts i:m:c:d:w:p:m:z: opt
  do
    case ${opt} in
    i) infile=${OPTARG};;
    m) infileMax=${OPTARG};;
    c) chrom_sizes=${OPTARG};;
    w) window_size=${OPTARG};;
    p) parallel=${OPTARG};;
    z) prog_compression=${OPTARG};;
    *) usage;;
    esac
  done
  if [ ! -n "${infile-}" ]; then usage; fi
  if [ ! -n "${chrom_sizes-}" ]; then usage; fi

  ### setup for later
  tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
  trap "test -d $tmpdir && rm -rf $tmpdir" 0 1 2 3 15
  window_half_size=$(( $window_size / 2 ))
  window_double_size=$(( $window_size * 2 ))

  SORT_OPT_BASE="--batch-size=100"
  SORT_OPT_BED="${SORT_OPT_BASE} -k1,1 -k2,2n -k3,3n"
  SORT_OPT_BED_NAME="${SORT_OPT_BASE} -k4,4"
  SORT_OPT_NAME="${SORT_OPT_BASE} -k1,1"

  PROG_COMP="$prog_compression -c"
  PROG_DECOMP="$prog_compression -d"
  if [ "${prog_compression}" = "zstd"  ]; then
    PROG_COMP="$PROG_COMP -T${parallel} "
  fi


  printf "### prepare target positions\n" >&2 
  prep_input $infile $window_double_size


  printf "### compute neighbouring signals\n" >&2
  accumulate_neiboring_signals \
    ${tmpdir}/potential_center.bg.COMP \
    ${tmpdir}/infile.fwd.bw \
    ${tmpdir}/infile.rev.bw \
    $window_size \
  | $PROG_COMP \
  > ${tmpdir}/potential_center.txt.COMP


  printf "### select bidirectional regions\n" >&2
  cat ${tmpdir}/potential_center.txt.COMP \
  | $PROG_DECOMP \
  | classify_convergent_divergent \
      ${tmpdir}/potential_center_conv.bed \
      ${tmpdir}/potential_center_divergent.bed

  printf "### find cores with left/right boundaries\n"  >&2
  find_max_scores \
    ${tmpdir}/potential_center_divergent.bed \
    ${tmpdir}/infile.density.fwd.bg \
    ${tmpdir}/infile.density.rev.bg \
    $window_size \
  > ${tmpdir}/potential_center_divergent_with_core.bed

  printf "### tighten\n"  >&2
  tighten_unidirectional ${tmpdir}/potential_center_divergent_with_core.bed \
  > ${tmpdir}/potential_center_divergent_with_core_tight.bed

  printf "### counts_fr\n"  >&2
  total=$( gunzip -c ${infile} | awk '{sum += $5 * ($3-$2)}END{print sum}' )
  counts_fr ${tmpdir}/potential_center_divergent_with_core_tight.bed ${infileMax} \
  | awk --assign total=$total 'BEGIN{OFS="\t"}{
      if (match( $4, /counts:[-.0-9,]+/)) {
        c = substr($4, RSTART, RLENGTH)
        sub("counts:", "", c)
        c += 0
      }
      tpm = 1e6 * c / total
      $4 = $4"|total:"total"|tpm:"tpm
      print
    }'
}



cmd_eachcount ()
{
  ### handle options
  parallel=20
  infile=
  infileEach=
  while getopts i:e:p:o: opt
  do
    case ${opt} in
    i) infile=${OPTARG};;
    e) infileEach=${OPTARG};;
    p) parallel=${OPTARG};;
    o) out_prefix=${OPTARG};;
    *) usage;;
    esac
  done
  if [ ! -n "${infile-}" ]; then usage; fi
  if [ ! -n "${infileEach-}" ]; then usage; fi
  if [ ! -n "${out_prefix-}" ]; then usage; fi

  ### setup for later
  tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
  trap "test -d $tmpdir && rm -rf $tmpdir" 0 1 2 3 15

  mkdir -p ${tmpdir}/each
  tar xvf ${infileEach} -C ${tmpdir}/each/


  ###
  ### prep regions (core + extended for forward or reverse)
  ###
  gunzip -c ${infile} \
  | sed -e 's/\t/#/g' \
  | awk 'BEGIN{OFS="\t";FS="#"}{
    chr=$1;start=$2;stop=$3;core_start=$7;core_stop=$8;name=$4
    print chr, core_start, stop, $_
  }' > ${tmpdir}/fwd.bed

  gunzip -c ${infile} \
  | sed -e 's/\t/#/g' \
  | awk 'BEGIN{OFS="\t";FS="#"}{
    chr=$1;start=$2;stop=$3;core_start=$7;core_stop=$8;name=$4
    print chr, start, core_stop, $_
  }' > ${tmpdir}/rev.bed


  ###
  ### prep totals
  ###
  find ${tmpdir}/each/ -name '*.fwd.bw' \
  | xargs -L 1 -I % bash -c \
     "bigWigToBedGraph % /dev/stdout | awk '{sum += \$4 * (\$3 - \$2)}END{print \"%,\"sum}'" \
  | sed 's/.fwd.bw,/\t/'  \
  | sed 's/^.*\///'  \
  | sort -k1,1 \
  > ${tmpdir}/totalsFwd.txt

  find ${tmpdir}/each/ -name '*.rev.bw' \
  | xargs -L 1 -I % bash -c \
     "bigWigToBedGraph % /dev/stdout | awk '{sum += \$4 * (\$3 - \$2)}END{print \"%,\"sum}'" \
  | sed 's/.rev.bw,/\t/'  \
  | sed 's/^.*\///'  \
  | sort -k1,1 \
  > ${tmpdir}/totalsRev.txt

  join -t "	" ${tmpdir}/totalsFwd.txt ${tmpdir}/totalsRev.txt \
  | awk 'BEGIN{OFS="\t"}{print $1,$2+$3}' \
  > ${tmpdir}/totals.txt


  ###
  ### counts in parallel
  ###
  cut -f 1 ${tmpdir}/totals.txt \
  | xargs -L 1 -P ${parallel} -I % bigWigAverageOverBed ${tmpdir}/each/%.fwd.bw ${tmpdir}/fwd.bed ${tmpdir}/each/%.fwd.bw.txt

  cut -f 1 ${tmpdir}/totals.txt \
  | xargs -L 1 -P ${parallel} -I % bigWigAverageOverBed ${tmpdir}/each/%.rev.bw ${tmpdir}/rev.bed ${tmpdir}/each/%.rev.bw.txt


  ###
  ### merge each counts 
  ###
  cut -f 1 ${tmpdir}/each/$( head -1 ${tmpdir}/totals.txt | cut -f 1 ).fwd.bw.txt \
  | tee ${tmpdir}/count.fwd.txt \
  > ${tmpdir}/count.rev.txt 

  for X in $( cat ${tmpdir}/totals.txt | cut -f 1  )
  do
    X=${tmpdir}/each/${X}.fwd.bw.txt
    join -t "	" ${tmpdir}/count.fwd.txt <(cut -f 1,4 $X) > ${tmpdir}/count.fwd.txt.tmp
    mv -f ${tmpdir}/count.fwd.txt.tmp ${tmpdir}/count.fwd.txt
  done

  for X in $( cat ${tmpdir}/totals.txt | cut -f 1  )
  do
    X=${tmpdir}/each/${X}.rev.bw.txt
    join -t "	" ${tmpdir}/count.rev.txt <(cut -f 1,4 $X) > ${tmpdir}/count.rev.txt.tmp
    mv -f ${tmpdir}/count.rev.txt.tmp ${tmpdir}/count.rev.txt
  done


  ###
  ### merge fwd & rev
  ###
  cat ${tmpdir}/count.fwd.txt \
  | awk '{buf="eachCountsFwd:";for(i=2;i <=NF; i++){buf=buf","$i};print $1"\t"buf}' \
  | sed 's/Fwd:,/Fwd:/' \
  > ${tmpdir}/count.fwd.txt.tmp
  mv -f ${tmpdir}/count.fwd.txt.tmp ${tmpdir}/count.fwd.txt

  cat ${tmpdir}/count.rev.txt \
  | awk '{buf="eachCountsRev:";for(i=2;i <=NF; i++){buf=buf","$i};print $1"\t"buf}' \
  | sed 's/Rev:,/Rev:/' \
  > ${tmpdir}/count.rev.txt.tmp
  mv -f ${tmpdir}/count.rev.txt.tmp ${tmpdir}/count.rev.txt

  join -t "	" \
    ${tmpdir}/count.fwd.txt \
    ${tmpdir}/count.rev.txt \
  > ${tmpdir}/count.fwdrev.txt

  eachTotals=$( cut -f 2 ${tmpdir}/totals.txt | xargs | sed -e 's/ /,/g' )
  cat ${tmpdir}/count.fwdrev.txt \
  | awk --assign eachTotals=$eachTotals 'BEGIN{OFS="\t"}{
    split($2,keyValue,":");fwdN=split(keyValue[2],eachCountsFwd,",");
    split($3,keyValue,":");revN=split(keyValue[2],eachCountsRev,",");
    buf="eachCounts:"
    bufD="eachDirectionalities:"
    for (i=1; i<=fwdN; i++) {
      tmpTotal = eachCountsFwd[i] + eachCountsRev[i]
      tmpDir = "NA"
      if (tmpTotal > 0) {
        tmpDir = (eachCountsFwd[i] - eachCountsRev[i]) / tmpTotal
        if (tmpDir < 0) {tmpDir = tmpDir * -1}
      }
      buf = buf","tmpTotal
      bufD = bufD","tmpDir
    }
    print $1,$2,$3,buf,bufD,"eachTotals:"eachTotals
  }' \
  | sed -e 's/eachCounts:,/eachCounts:/'  \
  | sed -e 's/eachDirectionalities:,/eachDirectionalities:/'  \
  | sed -e 's/#/\t/g' \
  | awk 'BEGIN{OFS="\t"}{$4=$4"|"$10"|"$11"|"$12"|"$13"|"$14;print}' \
  | cut -f 1-9 \
  | awk 'BEGIN{OFS="\t"}{
      if (match( $4, /eachCounts:[-.0-9,]+/)) {
        ec_str = substr($4, RSTART, RLENGTH)
        sub("eachCounts:", "", ec_str)
        ec_n = split(ec_str,ec,",")
      }
      if (match( $4, /eachTotals:[-.0-9,]+/)) {
        et_str = substr($4, RSTART, RLENGTH)
        sub("eachTotals:", "", et_str)
        et_n = split(et_str,et,",")
      }

      countsMax = 0
      tpmMax = 0
      for (i=1;i<=ec_n;i++)
      {
        tpm = 1e6 * (ec[i]+0)/(et[i]+0)
        if ( tpm > tpmMax){tpmMax = tpm}
        if ( ec[i] > countsMax){countsMax = ec[i]}
      }
      $4 = $4"|countsMax:"countsMax"|tpmMax:"tpmMax
      print
    }' \
  | gzip -c > ${out_prefix}_eachcounts.bed.gz

  mv -f ${tmpdir}/totals.txt ${out_prefix}_eachtotals.txt
}



cmd_filter ()
{
  ### handle options
  max_directionality=1
  min_directionality=0
  max_counts=-1
  min_counts=0
  max_tpm=-1
  min_tpm=0
  min_counts_in_each_strand=0
  min_ctssMax=0
  min_ctssMax_in_each_strand=0
  max_width=-1
  min_width=0
  invert_match="false"

  while getopts D:d:C:c:T:t:e:x:y:w:W:v opt
  do
    case ${opt} in
    D) max_directionality=${OPTARG};;
    d) min_directionality=${OPTARG};;
    C) max_counts=${OPTARG};;
    c) min_counts=${OPTARG};;
    T) max_tpm=${OPTARG};;
    t) min_tpm=${OPTARG};;
    e) min_counts_in_each_strand=${OPTARG};;
    x) min_ctssMax=${OPTARG};;
    y) min_ctssMax_in_each_strand=${OPTARG};;
    W) max_width=${OPTARG};;
    w) min_width=${OPTARG};;
    v) invert_match="true";;
    *) usage;;
    esac
  done


  awk \
    --assign max_directionality=$max_directionality \
    --assign min_directionality=$min_directionality \
    --assign max_counts=$max_counts \
    --assign min_counts=$min_counts \
    --assign max_tpm=$max_tpm \
    --assign min_tpm=$min_tpm \
    --assign min_counts_in_each_strand=$min_counts_in_each_strand \
    --assign min_ctssMax=$min_ctssMax \
    --assign min_ctssMax_in_each_strand=$min_ctssMax_in_each_strand \
    --assign max_width=$max_width \
    --assign min_width=$min_width \
    --assign invert_match=$invert_match \
  'BEGIN{OFS="\t"}{

    if (match( $0, /directionality:[-.0-9]+/)) {
      s = substr($0, RSTART, RLENGTH)
      sub("directionality:", "", s)
      d = s + 0;
      if (d < 0) { d = d * -1}
    }

    if (match( $0, /counts:[-.0-9]+/)) {
      c = substr($0, RSTART, RLENGTH)
      sub("counts:", "", c)
      c += 0;
    }
    if (match( $0, /countsFwd:[-.0-9]+/)) {
      f = substr($0, RSTART, RLENGTH)
      sub("countsFwd:", "", f)
      f += 0;
    }
    if (match( $0, /countsRev:[-.0-9]+/)) {
      r = substr($0, RSTART, RLENGTH)
      sub("countsRev:", "", r)
      r += 0;
    }

    cm = 0; cmf = 0 ; cmr = 0;
    if (match( $0, /ctssMax:[-.0-9]+/)) {
      cm = substr($0, RSTART, RLENGTH)
      sub("ctssMax:", "", cm)
      cm += 0;
    }
    if (match( $0, /ctssMaxFwd:[-.0-9]+/)) {
      cmf = substr($0, RSTART, RLENGTH)
      sub("ctssMaxFwd:", "", cmf)
      cmf += 0;
    }
    if (match( $0, /ctssMaxRev:[-.0-9]+/)) {
      cmr = substr($0, RSTART, RLENGTH)
      sub("ctssMaxRev:", "", cmr)
      cmr += 0;
    }
    if (match( $0, /tpm:[-.0-9]+/)) {
      tpm = substr($0, RSTART, RLENGTH)
      sub("tpm:", "", tpm)
      tpm += 0;
    }

    if ( ( cm == 0 ) && ( min_ctssMax > 0) ) {
      print "Error: No ctssMax found. Related filters cannot be used.\n" > "/dev/stderr"
      exit 1
    }
    if ( ( cm == 0 ) && ( min_ctssMax_in_each_strand > 0) ) {
      print "Error: No ctssMax found. Related filters cannot be used.\n" > "/dev/stderr"
      exit 1
    }

    width = $3 - $2

    flagM = "true"
    if      ( d > max_directionality )           { flagM = "false" }
    else if ( d < min_directionality )           { flagM = "false" }

    else if ( ( max_counts >= 0) &&
              ( c > max_counts ) )               { flagM = "false" }
    else if ( c < min_counts )                   { flagM = "false" }

    else if ( ( max_tpm >= 0  ) &&
              ( tpm > max_tpm ) )                { flagM = "false" }
    else if ( tpm < min_tpm )                    { flagM = "false" }

    else if ( f < min_counts_in_each_strand )    { flagM = "false" }
    else if ( r < min_counts_in_each_strand )    { flagM = "false" }

    else if ( cm < min_ctssMax )                 { flagM = "false" }
    else if ( cmf < min_ctssMax_in_each_strand ) { flagM = "false" }
    else if ( cmr < min_ctssMax_in_each_strand ) { flagM = "false" }

    else if ( ( max_width >= 0  ) &&
              ( width > max_width ) )                { flagM = "false" }
    else if ( width < min_width )                    { flagM = "false" }


    # invert
    if ( invert_match == "true" ) {
      if ( flagM == "true" ) { flagM = "false" } else { flagM = "true" }
    }
    #color = sprintf( "%d,0,%d", 127 + 127 * s , 127 - 127 * s )
    #$9 = color

    if ( flagM == "true" ) { print }
  }'
}




cmd_classify ()
{
  ### handle options
  tpm_threshold=2
  while getopts t: opt
  do
    case ${opt} in
    t) tpm_threshold=${OPTARG};;
    *) usage;;
   esac
  done

  awk --assign tpm_threshold=$tpm_threshold 'BEGIN{OFS="\t"}{
    if (match( $0, /tpm:[-.0-9]+/)) {
      tpm = substr($0, RSTART, RLENGTH)
      sub("tpm:", "", tpm)
      tpm += 0;
    }
    if (match( $0, /directionality:[-.0-9]+/)) {
      s = substr($0, RSTART, RLENGTH)
      sub("directionality:", "", s)
      d = s + 0;
    }

    if (tpm >= tpm_threshold ) { $4 = "class:PrmL|"$4}
    else { $4 = "class:EnhL|"$4; $9 = "253,167,1"}
    print
  }'
}


cmd_run ()
{

  ### handle options
  window_size=200
  parallel=20
  prog_compression=gzip

  while getopts i:c:o:w:p:z: opt
  do
    case ${opt} in
    i) infilesList=${OPTARG};;
    c) chrom_sizes=${OPTARG};;
    o) out_prefix=${OPTARG};;
    w) window_size=${OPTARG};;
    p) parallel=${OPTARG};;
    z) prog_compression=${OPTARG};;
    *) usage;;
   esac
  done

  if [ ! -n "${infilesList-}" ]; then usage; fi
  if [ ! -n "${chrom_sizes-}" ]; then usage; fi
  if [ ! -n "${out_prefix-}" ]; then usage; fi

  $0 total \
    -i $infilesList \
    -c $chrom_sizes \
    -o ${out_prefix} \
    -p ${parallel}

  $0 call \
    -i ${out_prefix}_total.ctss.bed.gz \
    -c ${chrom_sizes} \
    -w ${window_size} \
    -m ${out_prefix}_max.ctss.bed.gz \
    -p ${parallel} \
    -z ${prog_compression} \
  | gzip -c > ${out_prefix}_region.bed.gz

  $0 eachcount \
    -i ${out_prefix}_region.bed.gz \
    -e ${out_prefix}_each.ctss.bw.tar \
    -o ${out_prefix}_region \
    -p ${parallel}

  #gunzip -c ${out_prefix}_region_eachcounts.bed.gz \
  #| $0 filter \
  #| gzip -c > ${out_prefix}_region_eachcounts_filtered.bed.gz
}


#--------------------------------------------------------
# main
#--------------------------------------------------------

case "${1:-}" in
  run)
    shift;
    cmd_run "$@"
    ;;
  total)
    shift;
    cmd_total "$@"
    ;;
  call)
    shift;
    cmd_call "$@"
    ;;
  eachcount)
    shift;
    cmd_eachcount "$@"
    ;;
  filter)
    shift;
    cmd_filter "$@"
    ;;
  classify)
    shift;
    cmd_classify "$@"
    ;;
  *)
    usage
    ;;
esac




