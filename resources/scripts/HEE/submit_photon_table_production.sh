#!/bin/bash

# possible production script for photon tables
# has to be adjusted to fit different purposes
# paths need to be changed before usage!


STDOUT="/lustre/fs20/group/icecube/tkittler/stdout/"
PRESCALE=1000
NEVENTS=7500     # number of events per table
NTABLES=5     # number of tables per zenith angle
SOURCE="infinite-muon"     # cascade, infinite-muon, flasher
PREFIX="photonTables__mDOM_${SOURCE}_N${NEVENTS}_p${PRESCALE}"
#PREFIX="photon_table_photon_flux_cascade_N${NEVENTS}"
FOLDER="/lustre/fs20/group/icecube/tkittler/photonTables/tableOutput"
POSTFIX=".fits"
SCRIPT="/afs/ifh.de/group/amanda/scratch/tkittler/simulation/iceTray/combo/trunk/src/clsim/resources/tablemaker/photomc.py"
RANDOM=$(date +%s)
NOFFSET=0
ICEMODEL="homogeneous_ice"


#for ZENITH in $(seq 0 10 180)
for ZENITH in 180
do
    for TN_ in $(seq ${NTABLES})
    do
	TN=$(($TN_ + $NOFFSET))
        SEED=$(/afs/ifh.de/user/t/tkittler/bin/random_int.py --high 1000000)
        JOBNAME="PT_tracks_tn${TN}_n${NEVENTS}_p${PRESCALE}_Z${ZENITH}"
        FILENAME="${FOLDER}/${PREFIX}_Z${ZENITH}_TN${TN}${POSTFIX}"
        COMMAND="qsub -N ${JOBNAME} -j y -o ${STDOUT} -l h_rss=4G ${SCRIPT} --seed ${SEED} --zenith ${ZENITH} --prescale ${PRESCALE} --ice-model ${ICEMODEL} --nevents ${NEVENTS} --record-errors --light-source ${SOURCE} ${FILENAME}"
        #echo ${COMMAND}
        ${COMMAND}
    done
done

