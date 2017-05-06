#QERUN.SH V0 by Michael Lee
#INPUT: *.xyz
#OUPUT: *.in & *pbs (with option to qsub)

update_elem () {
	case $1 in
		'H')
			pseudopot=""
			val=1;;
		'C')
			pseudopot="C 12.0107 C.pbe-n-rrkjus_psl.0.1.UPF"
			val=4;;
		'O')
			pseudopot=""
			val=6;;
		'F')
			pseudopot="F 18.99840325 F.pbe-n-rrkjus_psl.0.1.UPF"
			val=7;;
		'Pt')
			pseudopot=""
			val=10;;
		*)
			echo "Invalid Element" $1
	esac
}

printparams () {
	echo "Writing $prefix.in with parameters:"
	printf "NATOM: $NATOM \n"
	printf "NTYP: $NTYP \n"
	printf "NBND: $NBND \n"
	printf "ecutwfc = %.1f \n" $ecutwfc
	printf "ecutrho = %.1f \n" $ecutrho
	printf "ecutrho = %.1g \n" $degauss
	printf "KPOINTS: %s \n" "${kpoints[@]}"
}

rm_filelist (){ #rm .tmp files
	rm NAT_TYP_ARR.tmp
}

GLOBAL_ECUTWFC=40.0
GLOBAL_ECUTRHO=480.0
GLOBAL_DEGAUSS=0.001
GLOBAL_KPOINTS=(18 18 1 0 0 0)
while read -r -u "$pf" prefix
do
	ecutwfc=$GLOBAL_ECUTWFC
	ecutrho=$GLOBAL_ECUTRHO
	degauss=$GLOBAL_DEGAUSS
	kpoints=${GLOBAL_KPOINTS[@]}
	pseudopot=""
	val=-1
	INPUT_FILE=$(printf $prefix".in")
	PBS_FILE=$(printf "run-"$prefix".pbs")
	XYZ=$(printf $prefix".xyz")
	
	if [ ! -f $XYZ ] #Check if XYZ file exists
	then
		echo "$XYZ not found. Terminating job."
		continue
	else
		echo "Analysing $XYZ..."
	fi

	if [ ! -d $prefix ] #Check if prefix directory exists
	then
		echo "Subdirectory $prefix not found"
		read -r -p "mkdir $prefix [m], run here[h]? Press any other key to skip run" reply
		case "$reply" in
			[mM])
				echo "Making subdirectory $prefix"
				mkdir $prefix
				cd $prefix
				cp ../$XYZ .
				;;
			[hH])
				echo "Running $prefix in this folder";;
			*)
				echo "Skipping $prefix"
				continue;;
		esac
	else
		cd $prefix
		if [ -f $INPUT_FILE ] #Check if input file exists
		then
			read -r -p "Replace $INPUT_FILE?[y/n]" reply
			case $reply in
				[yY])
				        echo "Overwriting $INPUT_FILE"	
					rm $INPUT_FILE;;
				*) 
					echo "Skipping $prefix"
					cd ..
					continue;;
			esac
		fi
		cp ../$XYZ .
	fi
#START XYZ ANALYSIS
	tail -n+3 $XYZ | cut -c1-2 | sort | uniq -c > NAT_TYP_ARR.tmp #array of elements
	TYP_ARR=($(cut -c9-10 NAT_TYP_ARR.tmp))  #array of elements
	NATOM_ARR=($(cut -c6-7 NAT_TYP_ARR.tmp))  #array of number of atoms
	NTYP=${#TYP_ARR[@]} #number of unique elements
	NATOM=$(head -n1 $XYZ | tr -d '\r') #can check with sum atoms
	NBND=0 
	PSEUDO_FILE="pseudo-template"
	NATOM_CHK=0
	elem_index=0
	for i in "${TYP_ARR[@]}"
	do
		#PSEUDO_POTS[$elem_index]=$(grep -w $i $PSEUDO_FILE) #match whole words (not 'F' in 'UPF')
		update_elem $i
		NBND=$(($NBND+(${NATOM_ARR[$elem_index]}*val+1)/2))
		NATOM_CHK=$(($NATOM_CHK+${NATOM_ARR[$elem_index]}))
		elem_index=$(($elem_index+1))
		PSEUDO_POTS[$elem_index]=$pseudopot #match whole words (not 'F' in 'UPF')
	done
	NBND=$(($NBND+10)) #buffer NBND by 10
#END XYZ ANLAYSIS
	if [ ! $NATOM_CHK -eq $NATOM ] #Check no. of atoms
	then
		echo "Number of atoms in $XYZ is wrong. Skipping $prefix..."
		continue
	fi
	EDIT_PARAMS="y"
	while true
	do
		printparams
		read -r -p "Change: ecutwfc/ecutrho[e] K-points[k] degauss[d]? Press any other key to continue" reply
	       	case "$reply" in
			[eE]) 
				read -r -p "Enter new value for ecutwfc: " ecutwfc
				read -r -p "Enter new value for ecutrho: " ecutrho
				read -r -p "Apply to all files?[y/n]" yn
				if [ yn = "y" ]
				then GLOBAL_ECUTWFC=$ecutwfc;GLOBAL_ECUTRHO=$ecutrho
				fi;;
			[dD]) 
				read -r -p "Enter new value for degauss: " degauss
				read -r -p "Apply to all files?[y/n]" yn
				if [ yn = "y" ]
				then GLOBAL_DEGAUSS=$degauss
				fi;;
			[kK]) 
				read -r -p "Enter new values for K Points: " kpoints;;
			*) break;;
		esac	
	done
printf "&control
    calculation = 'relax'
    restart_mode='from_scratch'
    title='%s'
    prefix='%s'
    nstep = 1000
 /
 &system
   ibrav = 12, celldm(1) = 27.8923581741617, celldm(2) = 1.0 , celldm(3) = 1.0
    celldm(4) = -0.5, celldm(5) = 0.0, celldm(6) = 0.0
    ecutwfc = %.1f
    ecutrho = %.1f
    nat = %i
    ntyp = %i 
    nbnd = %i
    occupations = 'smearing'
    smearing = 'cold'
    degauss = %.1g
    nspin = 2
    starting_magnetization = 1.0
 /
 &electrons
    electron_maxstep = 1000
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr =  1.0d-6
 /
 &ions
 ion_dynamics = 'damp'
 wfc_extrapolation = 'second order'
 pot_extrapolation = 'second order'
 /
 ATOMIC_SPECIES
" $prefix $prefix $ecutwfc $ecutrho $NATOM $NTYP $NBND $degauss >> $INPUT_FILE #print main parameters
	printf "%s\n" "${PSEUDO_POTS[@]}" >> $INPUT_FILE #print pseudopotentials
	printf "ATOMIC_POSITIONS (angstrom)\n" >> $INPUT_FILE
	tail -n+3 $XYZ | tr -d '\r' >> $INPUT_FILE #print coordinates, fixes non-unix endline escape key
	printf "\nK_POINTS (automatic) \n%s\n" "${kpoints[@]}" >> $INPUT_FILE

#PRINT TO PBS_FILE
echo "PRINTING $PBS_FILE"
printf "#!/bin/sh
#PBS -N %s
# ask for 3 chunks of 24 cores, flat MPI (48 processes in total)
#PBS -l select=3:ncpus=24:mpiprocs=24:ompthreads=1:mem=96gb
# modify for required runtime
#PBS -l walltime=20:00:00
# Send standard output and standard error to same output file
#PBS -j oe

inputfile=%s

echo 'PBS_O_WORKDIR is : \$PBS_O_WORKDIR'
echo 'PBS JOB DIR is: \$PBS_JOBDIR'
# Notice that the output of pwd will be in lustre scratch space
echo 'PWD is : \`pwd\`\'

# load default module environment for Quantum Espresso
# The default build is a flat MPI build
module load quantum-espresso

# change directory into the directory which has been copied into sandbox
cd \${PBS_O_WORKDIR}

# run the executable
mpirun pw.x -i \$inputfile > output-%s.\$PBS_JOBID

# all files will now be copied back into location specified by stageout setting
" $prefix $INPUT_FILE $prefix >> $PBS_FILE
read -r -p "qsub $PBS_FILE?[y/n]" yn
if [ yn = "y" ]
then
	echo "Submitting $PBS_FILE"
	qsub $PBS_FILE > jobid
fi
rm_filelist
cd ..
printf "\n"
done {pf}<prefixes
