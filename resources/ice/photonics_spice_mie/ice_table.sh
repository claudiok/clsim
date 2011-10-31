#!/bin/sh

model=mie
wvs=`seq 305 10 600`
DATE=`date +'%B %d, %G'`
export PPCTABLESDIR=dat/$model

if ! test -d $PPCTABLESDIR; then
cat << EOF
Run this script within the ppc directory with the compiled "ppc"
executable. ppc version must be v50, from Apr 08, 2011 or later.
Unpack the ice table archive dat.tgz within this same directory.
EOF
exit; fi

echr() { echo -ne "$*" 1>&2; }

cos=`awk '/^[ ]*[-.0-9]/ {n++; if(n==4) print $1}' $PPCTABLESDIR/cfg.txt`
echr "model=$model cos=$cos "

for w in $wvs; do
WFLA=$w ./ppc - > w.$w 2>n.$w
done

file=Ice_table.$model.i3coords.cos`printf %3.2f 0.9 | sed 's-\.--'`.`date +%d%b%G`.txt
echr "file=$file "

(
cat << EOF
# Ice table for the "SPICE $model" model (made on $DATE) for the depth range 1100-2810 m
NLAYER 171
NWVL 30 300 10
EOF

x=-855.4

for i in `seq 1 171`; do

y=`echo $x | awk '{print $1+10}'`
echo LAYER $x $y

echo -n ABS
for w in $wvs; do
awk $i'==NR {printf " "$2}' w.$w
done
echo

echo -n SCAT
for w in $wvs; do
awk $i'==NR {printf " "$3}' w.$w
done
echo

echo -n COS
for w in $wvs; do
echo -n " $cos"
done
echo

echo -n N_GROUP
for w in $wvs; do
awk 'BEGIN {FS="[ =]*"} /cm=/ {printf " "(0.299792458/$8)}' n.$w
done
echo

echo -n N_PHASE
for w in $wvs; do
awk 'BEGIN {FS="[ =]*"} /np=/ {printf " "$6}' n.$w
done
echo

x=$y
if test $[$i%40] = 0; then
if test $i -lt 150; then echr .; else echr " "; fi
fi

done
) > $file

for w in $wvs; do rm {w,n}.$w; done
echr "($DATE)\n"
