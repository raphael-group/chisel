# Checks
: ex: set ft=markdown ;:<<'```shell' #

This script runs all the tests to check that the current CHISEL implementation is correct and behaves as expected.

## Set up

```shell
set -e
set -o xtrace
PS4='[\t]'
cd $( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )
rm -rf X/ complete/ callingE/ cloningE/ plottingE/ pseudonormal/  
:<<'```shell' # Ignore this line
```

## Check function

```shell
check () {
      if cmp $1 $2
      then
          echo "CHECK $3: TEST $1 SUCCESS!"
      else
          echo "CHECK $3: TEST $1 FAILED!"
          exit 1
      fi
}
:<<'```shell' # Ignore this line
```

## Check complete

```shell
mkdir X/
sed 's/^chisel/python2.7 ..\/..\/bin\/chisel.py/g' ../demos/complete/demo-complete.sh > X/demo-complete.sh
curl -L https://github.com/raphael-group/chisel-data/raw/master/tests/complete.tar.gz | tar -xvz
check complete.chk <(bash X/demo-complete.sh |& grep -v -e "Progress:" -e "UserWarning" -e "--:--:--" -e "chisel" -e "curl" -e "Speed" -e "gzip" -e "samtools" -e "rundir" |& sed 's/\x1b\[[0-9;]*m//g' |& sed -u 's/\[[^]]*\]//g' |& tee log) "complete"
for F in complete/*.png; do check ${F} X/plots/$(basename ${F}) "complete"; done
check complete/calls.tsv X/calls/calls.tsv "complete"
check complete/mapping.tsv X/clones/mapping.tsv "complete"
rm -rf X/ complete/
:<<'```shell' # Ignore this line
```

## Check callingE

```shell
mkdir X/
sed 's/^chisel_calling/python2.7 ..\/..\/bin\/chisel_calling.py/g' ../demos/callingE/demo-callingE.sh > X/demo-callingE.sh
curl -L https://github.com/raphael-group/chisel-data/raw/master/tests/callingE.tar.gz | tar -xvz
check callingE.chk <(bash X/demo-callingE.sh |& grep -v -e "Progress:" -e "UserWarning" -e "--:--:--" -e "chisel" -e "curl" -e "Speed" -e "gzip" -e "samtools" -e "rundir" |& sed 's/\x1b\[[0-9;]*m//g' |& sed -u 's/\[[^]]*\]//g') "callingE"
for F in callingE/*.png; do check ${F} X/plots/$(basename ${F}) "callingE"; done
check callingE/calls.tsv X/calls/calls.tsv "callingE"
check callingE/mapping.tsv X/clones/mapping.tsv "callingE"
rm -rf X/ callingE/
:<<'```shell' # Ignore this line
```

## Check cloningE

```shell
mkdir X/
sed 's/^chisel_cloning/python2.7 ..\/..\/bin\/chisel_cloning.py/g' ../demos/cloningE/demo-cloningE.sh > X/demo-cloningE.sh
curl -L https://github.com/raphael-group/chisel-data/raw/master/tests/cloningE.tar.gz | tar -xvz
check cloningE.chk <(bash X/demo-cloningE.sh |& grep -v -e "Progress:" -e "UserWarning" -e "--:--:--" -e "chisel" -e "curl" -e "Speed" -e "gzip" -e "samtools" -e "rundir" |& sed 's/\x1b\[[0-9;]*m//g' |& sed -u 's/\[[^]]*\]//g') "cloningE"
for F in cloningE/*.png; do check ${F} X/plots/$(basename ${F}) "cloningE"; done
check cloningE/mapping.tsv X/clones/mapping.tsv "cloningE"
rm -rf X/ cloningE/
:<<'```shell' # Ignore this line
```

## Check plottingE

```shell
mkdir X/
sed 's/^chisel_plotting/python2.7 ..\/..\/bin\/chisel_plotting.py/g' ../demos/plottingE/demo-plottingE.sh > X/demo-plottingE.sh
curl -L https://github.com/raphael-group/chisel-data/raw/master/tests/plottingE.tar.gz | tar -xvz
check plottingE.chk <(bash X/demo-plottingE.sh |& grep -v -e "Progress:" -e "UserWarning" -e "--:--:--" -e "chisel" -e "curl" -e "Speed" -e "gzip" -e "samtools" -e "rundir" |& sed 's/\x1b\[[0-9;]*m//g' |& sed -u 's/\[[^]]*\]//g') "plottingE"
for F in plottingE/*.png; do check ${F} X/plots/$(basename ${F}) "plottingE"; done
rm -rf X/ plottingE/
:<<'```shell' # Ignore this line
```

## Checking pseudonormal

```shell
mkdir X/
sed 's/^chisel_pseudonormal/python2.7 ..\/..\/bin\/chisel_pseudonormal.py/g' ../demos/pseudonormal/demo-pseudonormal.sh > X/demo-pseudonormal.sh
curl -L https://github.com/raphael-group/chisel-data/raw/master/tests/pseudonormal.tar.gz | tar -xvz
check pseudonormal.chk <(bash X/demo-pseudonormal.sh |& grep -v -e "Progress:" -e "UserWarning" -e "--:--:--" -e "chisel" -e "curl" -e "Speed" -e "gzip" -e "samtools" -e "rundir" |& sed 's/\x1b\[[0-9;]*m//g' |& sed -u 's/\[[^]]*\]//g') "pseudonormal"
check pseudonormal/diploid.tsv X/diploid.tsv "pseudonormal"
rm -rf X/ pseudonormal/
:<<'```shell' # Ignore this line
```

## Successful checks

```shell
echo "ALL CHECKS PASSED SUCCESSFULLY!"
exit $?
```
