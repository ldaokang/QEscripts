#norm.sh for F-G group
#$ ./norm.sh (NATOMS) (NELECTRONS) < [projout] &> [norm-projout]
#(1) determines normalisation constant
#(2) outputs a normalised projwfc file
NATOMS=$1
NELECTRONS=$2
tail -n $((4*NATOMS+10)) | head -n $((4*NATOMS)) | sed -n '1~4p' > epop.tmp
cut -c35-40 epop.tmp > epop-charge.tmp
cut -c48-53 epop.tmp > epop-s.tmp
cut -c61-66 epop.tmp > epop-p.tmp
CHARGE_EXP=$(paste -sd+ epop-charge.tmp)
CHARGE_TOT=$(echo "scale = 4; $CHARGE_EXP"|bc)
NORM_CONST=$(echo "scale = 6; $NELECTRONS / $CHARGE_TOT" | bc)
echo "Normalisation completed"
echo "TOTAL CHARGE:" $CHARGE_TOT
echo "NORM CONST:" $NORM_CONST
echo "Ran with "$NATOMS " atoms and " $NELECTRONS "electrons."
printf "ATOM\ttot_q\ts_orb\tp_orb\n"
for ((i=1;i <= NATOMS ; i++)) 
do 
	CHARGE=$(sed -n "${i}p" epop-charge.tmp)
	NCHARGE=$(echo "scale = 4;$CHARGE * $NORM_CONST" | bc)
	S_ORB=$(sed -n "${i}p" epop-s.tmp)
	NS_ORB=$(echo "scale = 4;$S_ORB * $NORM_CONST" | bc)
	P_ORB=$(sed -n "${i}p" epop-p.tmp)
	NP_ORB=$(echo "scale = 4;$P_ORB * $NORM_CONST" | bc)
	printf "%u\t%.4f\t%.4f\t%.4f\n" $i $NCHARGE $NS_ORB $NP_ORB
done
rm epop.tmp epop-charge.tmp epop-s.tmp epop-p.tmp
