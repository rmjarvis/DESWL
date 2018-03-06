#!/bin/bash
function run {
    echo "host: $(hostname)"
    echo "tag = $tag"
    echo "exp = $exp"

    command=/direct/astro+u/mjarvis/rmjarvis/DESWL/psf/run_piff.py
    echo ${command} --mag_cut=1 --max_mag=21 --work=$work --tag=$tag --rm_files=1 --reserve=0.2 --use_ngmix --exps $exp --scratch=$TMPDIR --blacklist=0 --clear_output=0
    ${command} --mag_cut=1 --max_mag=21 --work=$work --tag=$tag --rm_files=1 --reserve=0.2 --use_ngmix --exps $exp --scratch=$TMPDIR --blacklist=0 --clear_output=0
    status=$?

    echo "time: $SECONDS"

    if [[ $status != "0" ]]; then
        echo "error running command: $status"
    fi

    return $status
}

tag=$1
exp=$2
work=$3
logfile=$4

cd /direct/astro+u/mjarvis/rmjarvis/DESWL/psf
source /direct/astro+u/mjarvis/.bashrc

if [[ -n $_CONDOR_SCRATCH_DIR ]]; then
    tmpdir=$_CONDOR_SCRATCH_DIR
    export TMPDIR=$tmpdir
else
    tmpdir='/data/mjarvis/y3_piff'
    export TMPDIR=$tmpdir
    mkdir -p $tmpdir
fi

pushd $tmpdir

tmplog=$(basename $logfile)

run &> ${tmplog}
status=$?

echo "moving log file ${tmplog} -> ${logfile}" >> ${tmplog}

mv -fv "${tmplog}" "${logfile}" 1>&2

popd
