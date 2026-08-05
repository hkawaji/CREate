#!/bin/bash
#
# Regression tests for CREate.
#
# Everything here runs off the EGR1 locus fixtures in this directory and takes a
# few seconds. Each test regenerates an output and compares it to a committed
# reference, so a change in any result is a failure until the reference is
# deliberately updated.
#
#   ./test/run_tests.sh          from the repository root
#
# Options are pinned rather than left to the defaults, so that changing a default
# does not silently rewrite what these tests check. Defaults are covered
# separately, at the end.

set -u

cd "$(dirname "$0")/.." || exit 1
CREATE=./CREate.sh
TESTDIR=test
tmpdir=$(mktemp -d -p "${TMPDIR:-/tmp}")
trap 'test -d "$tmpdir" && rm -rf "$tmpdir"' 0 1 2 3 15

pass=0
fail=0

report () {  # name, ok(0/1), detail
  if [ "$2" -eq 0 ]; then
    printf 'ok   %s\n' "$1"
    pass=$((pass + 1))
  else
    printf 'FAIL %s\n' "$1"
    [ -n "${3-}" ] && printf '     %s\n' "$3"
    fail=$((fail + 1))
  fi
}

same () {  # name, produced, reference
  if diff <(zcat -f "$2") <(zcat -f "$3") > "$tmpdir/diff.txt" 2>&1; then
    report "$1" 0
  else
    report "$1" 1 "$(head -3 "$tmpdir/diff.txt" | tr '\n' ' ')"
  fi
}


### peak calling: total -> call -> filter -> classify

$CREATE total -i $TESTDIR/K562_EGR1.bed.gz,$TESTDIR/MCF7_EGR1.bed.gz \
              -c $TESTDIR/hg38.chrom_sizes -o "$tmpdir/egr1" > /dev/null 2>&1
$CREATE call -i "$tmpdir/egr1_total.ctss.bed.gz" -c $TESTDIR/hg38.chrom_sizes \
             -m "$tmpdir/egr1_max.ctss.bed.gz" 2>/dev/null \
| $CREATE filter 2>/dev/null \
| $CREATE classify -e 2 2>/dev/null \
| gzip -c > "$tmpdir/classify.bed.gz"
same "total/call/filter/classify" "$tmpdir/classify.bed.gz" $TESTDIR/EGR1_create.bed.gz


### makeprom + promact: promoter regions and their activity

$CREATE makeprom -g $TESTDIR/gene.bed12.gz -c $TESTDIR/hg38.chrom_sizes 2>/dev/null \
| gzip -c > "$tmpdir/promoter.bed.gz"
same "makeprom" "$tmpdir/promoter.bed.gz" $TESTDIR/promoter.bed.gz

$CREATE promact -i $TESTDIR/promoter.bed.gz -r $TESTDIR/K562_EGR1.bed.gz 2>/dev/null \
| gzip -c > "$tmpdir/promact.bed.gz"
same "promact" "$tmpdir/promact.bed.gz" $TESTDIR/EGR1_promact.bed.gz


### creact: CRE activity from a single-sample CTSS

$CREATE creact -i $TESTDIR/EGR1_create.bed.gz -r $TESTDIR/K562_EGR1.bed.gz 2>/dev/null \
| gzip -c > "$tmpdir/creact.bed.gz"
same "creact" "$tmpdir/creact.bed.gz" $TESTDIR/EGR1_creact.bed.gz


### hicprep: bin size and the power law it fits
#
# synth_contacts.bed.gz is built from gamma = 1.0, scale = 5.0 at 5kb bins, with
# one bin in seven left out entirely. Dropping whole bins keeps the decay intact
# and lowers the intercept by about log((6/7)^2) = -0.31, so gamma should come
# back near 1.0 and scale near 4.7.

$CREATE hicprep -i $TESTDIR/synth_contacts.bed.gz -k $TESTDIR/synth_promact.bed.gz \
        2>/dev/null > "$tmpdir/contacts.tsv"

hp_res=$(  sed -n 's/^#resolution=//p' "$tmpdir/contacts.tsv")
hp_gamma=$(sed -n 's/^#gamma=//p'      "$tmpdir/contacts.tsv")
hp_scale=$(sed -n 's/^#scale=//p'      "$tmpdir/contacts.tsv")

[ "$hp_res" = "5000" ]
report "hicprep detects the 5kb bin size" $? "got '$hp_res'"

awk -v g="$hp_gamma" 'BEGIN{ exit !(g > 0.95 && g < 1.05) }'
report "hicprep recovers gamma near 1.0" $? "got '$hp_gamma'"

awk -v s="$hp_scale" 'BEGIN{ exit !(s > 4.4 && s < 5.0) }'
report "hicprep recovers scale near 4.7" $? "got '$hp_scale'"


### link, without and with measured contact
#
# synth_promact.bed.gz holds the real EGR1 promoter plus two synthetic ones, so
# that more than one gene is normalised and the two differ. One of the synthetic
# promoters sits outside every CRE, which is what puts an uncovered promoter into
# the ABC candidate set.

$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
             -S 0 -p 0 -W 200000 2>/dev/null | gzip -c > "$tmpdir/link.bed.gz"
same "link" "$tmpdir/link.bed.gz" $TESTDIR/synth_link.bed.gz

$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
             -S 0 -p 0 -W 200000 -i "$tmpdir/contacts.tsv" 2>/dev/null \
| gzip -c > "$tmpdir/link_hic.bed.gz"
same "link -i (measured contact)" "$tmpdir/link_hic.bed.gz" $TESTDIR/synth_link_hic.bed.gz

if diff -q <(zcat "$tmpdir/link.bed.gz") <(zcat "$tmpdir/link_hic.bed.gz") > /dev/null 2>&1
then report "link -i changes the score" 1 "output identical with and without -i"
else report "link -i changes the score" 0
fi

### -i must also take a stream, since it is read for both the header and the data
$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
             -S 0 -p 0 -W 200000 -i <(cat "$tmpdir/contacts.tsv") 2>/dev/null \
| gzip -c > "$tmpdir/link_hic_pipe.bed.gz"
same "link -i from a process substitution" "$tmpdir/link_hic_pipe.bed.gz" $TESTDIR/synth_link_hic.bed.gz


### -N drops the pseudocount, which the ABC model also does for average Hi-C

hp_frac=$(sed -n 's/^#recordedFraction=//p' "$tmpdir/contacts.tsv")
awk -v f="$hp_frac" 'BEGIN{ exit !(f > 0.70 && f < 0.76) }'
report "hicprep reports the recorded fraction" $? "got '$hp_frac', expected ~0.735 = (6/7)^2"

$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
             -S 0 -p 0 -W 200000 -i "$tmpdir/contacts.tsv" -N 2>/dev/null \
| gzip -c > "$tmpdir/link_nopseudo.bed.gz"
if diff -q <(zcat "$tmpdir/link_hic.bed.gz") <(zcat "$tmpdir/link_nopseudo.bed.gz") > /dev/null 2>&1
then report "-N changes the contact term" 1 "output identical with and without -N"
else report "-N changes the contact term" 0
fi


### defaults
#
# Not a reference comparison: these only check that the documented defaults are
# what the code uses, so that a change to one has to be made on purpose.

# Run against -O cre -B 0.2 rather than the defaults, because on this three-promoter
# fixture the default denominator is large enough that every ABC lands below 1.8e-05:
# 0.001 and 0.016 would both keep nothing and the check would pass vacuously. The
# legacy rule spreads ABC over 8.9e-04 to 0.29, which separates the two candidates.
leg="-O cre -B 0.2 -p 0 -W 200000"
n_abc=$($CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
                $leg 2>/dev/null | wc -l)
n_cut=$($CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
                $leg -s ABC -S 0.001 2>/dev/null | wc -l)
n_old=$($CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
                $leg -s ABC -S 0.016 2>/dev/null | wc -l)
if [ "$n_abc" = "$n_cut" ] && [ "$n_abc" != "$n_old" ]; then r=0; else r=1; fi
report "default score and cutoff are ABC / 0.001" $r \
       "default kept $n_abc, -S 0.001 kept $n_cut, -S 0.016 kept $n_old"

n_aci=$($CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
                -p 0 -W 200000 -s ACI 2>/dev/null | wc -l)
n_aci1=$($CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
                 -p 0 -W 200000 -s ACI -S 1 2>/dev/null | wc -l)
[ "$n_aci" = "$n_aci1" ]
report "-s ACI defaults its cutoff to 1" $? "default kept $n_aci, -S 1 kept $n_aci1"




### creact -m: divergent (default) against whole-region

# -m divergent is the default, so naming it must change nothing.
$CREATE creact -m divergent -i $TESTDIR/EGR1_create.bed.gz -r $TESTDIR/K562_EGR1.bed.gz 2>/dev/null \
| gzip -c > "$tmpdir/creact_div.bed.gz"
same "creact -m divergent equals the default" "$tmpdir/creact_div.bed.gz" $TESTDIR/EGR1_creact.bed.gz

# -m whole ignores thickStart/thickEnd, so it must agree with the same regions fed
# in with the core widened to the whole feature -- which is how this was done
# before the option existed.
zcat -f $TESTDIR/EGR1_create.bed.gz \
| gawk 'BEGIN{OFS="\t"}{ $7 = $2; $8 = $3; print }' | gzip -c > "$tmpdir/wide.bed.gz"
$CREATE creact -i "$tmpdir/wide.bed.gz" -r $TESTDIR/K562_EGR1.bed.gz 2>/dev/null \
| gawk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $5 }' > "$tmpdir/wide.act.txt"
$CREATE creact -m whole -i $TESTDIR/EGR1_create.bed.gz -r $TESTDIR/K562_EGR1.bed.gz 2>/dev/null \
| gawk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $5 }' > "$tmpdir/whole.act.txt"
same "creact -m whole equals widening the core by hand" "$tmpdir/whole.act.txt" "$tmpdir/wide.act.txt"

# Measuring the whole feature cannot give the same answer as measuring the core,
# or the option would be doing nothing.
if cmp -s "$tmpdir/creact_div.bed.gz" "$tmpdir/whole.act.txt"; then r=1; else r=0; fi
report "creact -m whole changes the activity" $r

$CREATE creact -m sideways -i $TESTDIR/EGR1_create.bed.gz -r $TESTDIR/K562_EGR1.bed.gz \
        > /dev/null 2> "$tmpdir/mode_err.txt"
grep -q "must be divergent or whole" "$tmpdir/mode_err.txt"
report "creact rejects an unknown -m" $?


### link -t: score one set of regions, normalise by another

# Omitting -t must reproduce the run that has no -t at all.
$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
             -t $TESTDIR/EGR1_creact.bed.gz -S 0 -p 0 -W 200000 2>/dev/null \
| gzip -c > "$tmpdir/link_t_self.bed.gz"
same "link -t pointed at -r equals no -t" "$tmpdir/link_t_self.bed.gz" $TESTDIR/synth_link.bed.gz

# With a different -t the reported elements are the target's, not the candidates'.
zcat -f $TESTDIR/EGR1_creact.bed.gz \
| gawk 'BEGIN{OFS="\t"} NR % 3 == 1 { $2 = $2 - 200; $3 = $3 + 200; $7 = $2; $8 = $3; print }' \
| gzip -c > "$tmpdir/target.bed.gz"
n_target=$(zcat -f "$tmpdir/target.bed.gz" | wc -l)
$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
             -t "$tmpdir/target.bed.gz" -S 0 -p 0 -W 200000 2>/dev/null > "$tmpdir/link_t.bedpe"
n_out=$(cut -f1-3 "$tmpdir/link_t.bedpe" | sort -u | wc -l)
n_in=$(zcat -f "$tmpdir/target.bed.gz" | cut -f1-3 | sort -u | wc -l)
[ "$n_out" -gt 0 ] && [ "$n_out" -le "$n_in" ]
report "link -t reports the target regions" $? "target had $n_in distinct regions, output had $n_out"

# The denominator must not move: it is built from -r, which did not change. ABC is
# ACI over a per-gene constant, so that constant has to match the no-t run.
den () {  # first pair of a given gene -> ACI / ABC
  gawk -v g="$1" 'BEGIN{FS="\t"}
    { name = $7
      if (name !~ ("[|]tgt_id:" g "[|]")) next
      aci = name; sub(/.*[|]link_ACI:/, "", aci); sub(/[|].*/, "", aci)
      abc = name; sub(/.*[|]link_ABC:/, "", abc)
      if (abc + 0 > 0) { printf "%.4g\n", aci / abc; exit } }' "$2"
}
gene=$(zcat -f $TESTDIR/synth_promact.bed.gz | head -1 | cut -f4)
zcat -f $TESTDIR/synth_link.bed.gz > "$tmpdir/link_plain.bedpe"
d_plain=$(den "$gene" "$tmpdir/link_plain.bedpe")
d_t=$(den "$gene" "$tmpdir/link_t.bedpe")
[ -n "$d_plain" ] && [ "$d_plain" = "$d_t" ]
report "link -t leaves the denominator alone" $? "gene $gene: without -t $d_plain, with -t $d_t"



### column requirements

# -m whole only reads chrom/start/end, so the thick columns are not needed and a
# plain BED6 -- a DNase peak file, say -- must give the same activities as the
# BED9 the regions came from.
$CREATE creact -m whole -i $TESTDIR/EGR1_create.bed.gz -r $TESTDIR/K562_EGR1.bed.gz 2>/dev/null \
| gawk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $5 }' > "$tmpdir/act9.txt"
zcat -f $TESTDIR/EGR1_create.bed.gz | cut -f1-6 | gzip -c > "$tmpdir/b6.bed.gz"
$CREATE creact -m whole -i "$tmpdir/b6.bed.gz" -r $TESTDIR/K562_EGR1.bed.gz 2>/dev/null \
| gawk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $5 }' > "$tmpdir/act6.txt"
same "creact -m whole takes BED6" "$tmpdir/act6.txt" "$tmpdir/act9.txt"

# BED4 has no score column, so the activity has to be appended as a new column 5
# rather than overwriting the name or falling off the end of the record.
zcat -f $TESTDIR/EGR1_create.bed.gz | cut -f1-4 | gzip -c > "$tmpdir/b4.bed.gz"
$CREATE creact -m whole -i "$tmpdir/b4.bed.gz" -r $TESTDIR/K562_EGR1.bed.gz 2>/dev/null \
> "$tmpdir/b4.out"
n_col=$(head -1 "$tmpdir/b4.out" | gawk '{print NF}')
gawk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $5 }' "$tmpdir/b4.out" > "$tmpdir/act4.txt"
if [ "$n_col" = "5" ] && cmp -s "$tmpdir/act4.txt" "$tmpdir/act9.txt"; then r=0; else r=1; fi
report "creact -m whole puts BED4 activity in column 5" $r "output had $n_col columns"

# A shorter input would be measured over an empty extent and silently reported as
# zero activity, so it has to be refused instead.
zcat -f $TESTDIR/EGR1_create.bed.gz | cut -f1-3 | gzip -c > "$tmpdir/b3.bed.gz"
$CREATE creact -m whole -i "$tmpdir/b3.bed.gz" -r $TESTDIR/K562_EGR1.bed.gz \
        > /dev/null 2> "$tmpdir/b3.err" && r=1 || r=0
grep -q "needs at least 4" "$tmpdir/b3.err" || r=1
report "creact -m whole rejects BED3" $r "$(head -1 "$tmpdir/b3.err")"

# The default mode reads thickStart/thickEnd, which BED6 does not carry.
$CREATE creact -i "$tmpdir/b6.bed.gz" -r $TESTDIR/K562_EGR1.bed.gz \
        > /dev/null 2> "$tmpdir/div.err" && r=1 || r=0
grep -q "divergent needs at least 8" "$tmpdir/div.err" || r=1
report "creact -m divergent rejects BED6" $r "$(head -1 "$tmpdir/div.err")"

# promact reads the strand, which BED4 does not carry.
$CREATE promact -i "$tmpdir/b4.bed.gz" -r $TESTDIR/K562_EGR1.bed.gz \
        > /dev/null 2> "$tmpdir/prom.err" && r=1 || r=0
grep -q "sense needs at least 6" "$tmpdir/prom.err" || r=1
report "promact rejects BED4" $r "$(head -1 "$tmpdir/prom.err")"


### link -B: a separate exponent for the promoter half of the denominator

# Naming the default has to change nothing.
$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
             -B 1.0 -S 0 -p 0 -W 200000 2>/dev/null | gzip -c > "$tmpdir/link_B_same.bed.gz"
same "link -B 1.0 equals the default" "$tmpdir/link_B_same.bed.gz" $TESTDIR/synth_link.bed.gz

# -B is a denominator-only knob, so ABC must move and ACI must not.
$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
             -B 0.6 -S 0 -p 0 -W 200000 2>/dev/null > "$tmpdir/link_B06.bedpe"
zcat -f $TESTDIR/synth_link.bed.gz > "$tmpdir/link_B_ref.bedpe"
aci () { grep -o 'link_ACI:[^|]*' "$1" | md5sum | cut -d' ' -f1; }
abc () { grep -o 'link_ABC:[^|]*' "$1" | md5sum | cut -d' ' -f1; }
if [ "$(aci "$tmpdir/link_B06.bedpe")" = "$(aci "$tmpdir/link_B_ref.bedpe")" ]; then r=0; else r=1; fi
report "link -B leaves ACI alone" $r
if [ "$(abc "$tmpdir/link_B06.bedpe")" != "$(abc "$tmpdir/link_B_ref.bedpe")" ]; then r=0; else r=1; fi
report "link -B changes ABC" $r

# Under -s ACI the reported score is ACI, which -B cannot touch, so -B must change
# neither the values nor which pairs are reported. The name field still carries
# link_ABC, which -B does move, so that one field is stripped before comparing --
# and the rows are sorted, since the ranking is over a column that did not change
# and ties may fall either way.
strip_abc () { sed 's/|link_ABC:[^|\t]*//' "$1" | sort; }
$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
             -s ACI -S 0 -p 0 -W 200000 2>/dev/null > "$tmpdir/link_aci_plain.bedpe"
$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
             -s ACI -B 0.6 -S 0 -p 0 -W 200000 2>/dev/null > "$tmpdir/link_aci_B.bedpe"
strip_abc "$tmpdir/link_aci_plain.bedpe" > "$tmpdir/aci_plain.txt"
strip_abc "$tmpdir/link_aci_B.bedpe"     > "$tmpdir/aci_B.txt"
same "link -B does not affect what -s ACI reports" "$tmpdir/aci_B.txt" "$tmpdir/aci_plain.txt"


### link -O: which of a promoter and an overlapping CRE goes into the denominator

# -O cre is the pre-2026-08 rule. Together with the exponent and cutoff defaults of
# that time it has to reproduce results recorded then, byte for byte -- that is the
# whole reason the option is kept, so it gets its own committed reference.
$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
             -O cre -B 0.2 -S 0 -p 0 -W 200000 2>/dev/null | gzip -c > "$tmpdir/link_Ocre.bed.gz"
same "link -O cre -B 0.2 reproduces the pre-2026-08 output" \
     "$tmpdir/link_Ocre.bed.gz" $TESTDIR/synth_link_Ocre.bed.gz

# Naming the default has to change nothing.
$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
             -O promoter -S 0 -p 0 -W 200000 2>/dev/null | gzip -c > "$tmpdir/link_Oprom.bed.gz"
same "link -O promoter equals the default" "$tmpdir/link_Oprom.bed.gz" $TESTDIR/synth_link.bed.gz

$CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
        -O both > /dev/null 2> "$tmpdir/O_err.txt" && r=1 || r=0
grep -q "must be promoter or cre" "$tmpdir/O_err.txt" || r=1
report "link rejects an unknown -O" $r "$(head -1 "$tmpdir/O_err.txt")"

# -O feeds only the candidate set, which feeds only the denominator. So the two
# rules must return the same pairs with the same ACI and a different ABC.
o_aci () { grep -o 'link_ACI:[^|]*' "$1" | md5sum | cut -d' ' -f1; }
o_abc () { grep -o 'link_ABC:[^|]*' "$1" | md5sum | cut -d' ' -f1; }
zcat -f "$tmpdir/link_Ocre.bed.gz" > "$tmpdir/Ocre.bedpe"
zcat -f $TESTDIR/synth_link.bed.gz > "$tmpdir/Oprom.bedpe"
if [ "$(o_aci "$tmpdir/Ocre.bedpe")" = "$(o_aci "$tmpdir/Oprom.bedpe")" ]; then r=0; else r=1; fi
report "link -O leaves ACI alone" $r
if [ "$(o_abc "$tmpdir/Ocre.bedpe")" != "$(o_abc "$tmpdir/Oprom.bedpe")" ]; then r=0; else r=1; fi
report "link -O changes ABC" $r

cut -f1-3,7 "$tmpdir/Ocre.bedpe"  | sed 's/|link_ABC:[^|]*//' | sort > "$tmpdir/O_pairs_cre.txt"
cut -f1-3,7 "$tmpdir/Oprom.bedpe" | sed 's/|link_ABC:[^|]*//' | sort > "$tmpdir/O_pairs_prom.txt"
same "link -O reports the same pairs" "$tmpdir/O_pairs_cre.txt" "$tmpdir/O_pairs_prom.txt"

# Under -O cre the promoter half of the candidate set is a few percent of the
# denominator, so -B barely registers; under -O promoter it is most of it. Both
# still have to respond, or -B is not reaching the right terms.
for o in cre promoter; do
  $CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
               -O $o -B 0.6 -S 0 -p 0 -W 200000 2>/dev/null > "$tmpdir/link_O${o}_B06.bedpe"
  $CREATE link -k $TESTDIR/synth_promact.bed.gz -r $TESTDIR/EGR1_creact.bed.gz \
               -O $o -B 0.3 -S 0 -p 0 -W 200000 2>/dev/null > "$tmpdir/link_O${o}_B03.bedpe"
  if [ "$(o_abc "$tmpdir/link_O${o}_B06.bedpe")" != "$(o_abc "$tmpdir/link_O${o}_B03.bedpe")" ]
  then r=0; else r=1; fi
  report "link -B acts under -O $o" $r
done
printf '\n%d passed, %d failed\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
