#Variables
USER='mschlaepfer'
WORKDIR="/data/users/${USER}/rnaseq_course"
OUTDIR="${WORKDIR}/FEATURECOUNTS"
FCTABLE="${OUTDIR}/merged_counts.txt"
OUTFILE="${OUTDIR}/feature_counts_modified.txt"

#Take off Specifically, the first line and the columns containing Chr, Start, End, Strand and Length.
tail -n +2 ${FCTABLE} | cut --complement -f2,3,4,5,6 > ${OUTFILE}


#complement = "keep everything but 2, 3, 4, 5, 6"










