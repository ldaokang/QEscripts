#QERUN.SH V1
#MICHAEL LEE
#10 May 2017

#INPUT: *.xyz, prefixes
#OUTPUT: *.in & *.pbs (with option to submit job)
#updates: easier for manual edit of default parameters, pseudopotentials
#EDITABLE:

#DEFAULT PARMETERS:
DEF_ECUTWFC=40.0
DEF_ECUTRHO=480.0
DEF_DEGAUSS=0.001
DEF_KPOINTS=(18 18 1 0 0 0)

get_elem (){ #gets element & updates pseudopotential string & valence number
	case $1 in #only populated for ORR & F-G applications
		'H') val=1; pseudopot="" ;;
		'C') val=4; pseudopot="C 12.0107 C.pbe-n-rrkjus_psl.0.1.UPF" ;;
		'Si') val=4; pseudopot="" ;;
		'O') val=6; pseudopot="" ;;
		'F') val=7; pseudopot="F 18.99840325 F.pbe-n-rrkjus_psl.0.1.UPF" ;;
		'Pt') val=10; pseudopot="" ;;
		*) val=-1; pseudopot=""; echo "Invalid Element $1";;
	esac
}

INPUT_TEMPLATE="\
&control
    calculation = 'relax'
    restart_mode = 'from_scratch'
    title='%s'
    prefix='%s'
    nstep=1000
 /
 &system
 ibrav=12, celldm(1) = 27.8923581741617, celldm(2) = 1.0 , celldm(3) = 1.0
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
    conv_thr = 1.0d-6
 /
 &ions
 ion_dynamics = 'damp'
 wfc_extrapolation = 'second order'
 pot_extrapolation = 'second order'
 /
"

make_input (){ #parameters: [prefix,prefix,ecutwfc,ecutrho,nat,ntyp,nbnd,degauss]
	printf "$INPUT_TEMPLATE" $1 $1 $2 $3 $4 $5 $6 $7 $8
}

#DON'T EDIT:
origin=$PWD

terminate() {
	echo "ERROR: $1 terminated"
	cd $origin
}

print_params (){ #prints current input file parameters
	echo "$INPUT_FILE PARAMETERS:"
	printf "NATOM: %i \n" $NATOM
	printf "NTYP: %i \n" $NTYP
	printf "NBND: %i \n" $NBND
	printf "ecutwfc = %.1f \n" $ecutwfc
	printf "ecutrho = %.1f \n" $ecutrho
	printf "degauss = %.1g \n" $degauss
	printf "KPOINTS: %s \n" "${kpoints[*]}"
}

check_xyz_dir() { #Takes in XYZ_FILE and prefix

	if [ ! -f $1 ] #Check if XYZ file exists
	then
		echo "ERROR: $1 not found"
		terminate $2
		continue
	elif [ ! -d $2 ] #Check if prefix directory exists
	then
		echo "Directory $2 not found"
		read -r -p "mkdir $2 [m] or write here [h]? Terminate $2 [any other key]: " reply
		case "$reply" in
			[mM])
				echo "MAKING: subdir $2"; mkdir $2;;
			[hH])	echo "staying HERE $PWD";;
			*)	terminate $2; continue;;
		esac
	fi
}

rm_ifExist() { #Remove existing files
	if [ -f $1 ] 
		then
			echo "WARNING: $1 already exists"
			read -r -p "Replace $1?[y/n] " reply
			case $reply in
				[yY]) rm $1;;
				*) terminate $prefix; continue;;
			esac
	fi
}

analyse_xyz() { #Returns NATOM, NBND, NTYP, PSEUDOPOT_ARR
	tail -n+3 $1 | cut -c1-2 | sort | uniq -c > NAT_TYP_ARR.tmp #rough cut: NATOM lines of NTYP
	NATOM_ARR=($(tr -s ' ' < NAT_TYP_ARR.tmp | cut -d ' ' -f2))  #array of number of atoms
	TYP_ARR=($(tr -s ' ' < NAT_TYP_ARR.tmp | cut -d ' ' -f3))  #array of types of atoms
	NTYP=${#TYP_ARR[@]}
	NATOM=$(head -n1 $1 | tr -d '\r')
	elem_index=0
	NATOM_CHK=0

	for i in "${TYP_ARR[@]}"
	do
		pseudopot=""
		val=-1
		get_elem $i
		NBND=$(($NBND+(${NATOM_ARR[$elem_index]}*val+1)/2))
		NATOM_CHK=$(($NATOM_CHK+${NATOM_ARR[$elem_index]}))
		PSEUDOPOT_ARR[$elem_index]=$pseudopot
		elem_index=$(($elem_index+1))
	done
	NBND=$(($NBND+10)) #buffer NBND by 10

	if [ ! $NATOM_CHK -eq $NATOM ] #Check no. of atoms
	then
		echo "Incorrect NATOMS in $XYZ_FILE : $NATOM_CHK != $NATOM"
		terminate $prefix
		continue
	fi
	rm NAT_TYP_ARR.tmp
}

confirm_params(){
	while true
	do
		print_params
		read -r -p "CHANGE: ecutwfc/ecutrho[e] | degauss[d] | K-points[k] ? ACCEPT parameters [any other key]: " reply
	       	case "$reply" in
			[eE]) 
				read -r -p "Enter new value for ecutwfc: " ecutwfc
				read -r -p "Enter new value for ecutrho: " ecutrho
				;;
			[dD]) 
				read -r -p "Enter new value for degauss: " degauss
				;;
			[kK]) 

				read -r -p "New KPOINTS:" kp
				kpoints=($(printf "%i %i %i %i %i %i" $kp))
				;;

			*) break;;
		esac	
	done
}

while read -r -u "$pf" prefix
do
	XYZ_FILE=$(printf "$prefix.xyz")
	INPUT_FILE=$(printf "$prefix.in")
	PBS_FILE=$(printf "run-$prefix.pbs")
	
	check_xyz_dir $XYZ_FILE $prefix
	rm_ifExist $INPUT_FILE
	rm_ifExist $PBS_FILE
	cd $prefix
	cp ../$XYZ_FILE .
	
	NTYP=0;NATOM=0;NBND=0;PSEUDOPOT_ARR=("")
	analyse_xyz $XYZ_FILE

	ecutwfc=$DEF_ECUTWFC
	ecutrho=$DEF_ECUTRHO
	degauss=$DEF_DEGAUSS
	kpoints=("")
	kpoints="${DEF_KPOINTS[@]}"
	confirm_params
	
	make_input $prefix $ecutwfc $ecutrho $NATOM $NTYP $NBND $degauss > $INPUT_FILE
	printf " ATOMIC_SPECIES\n" >> $INPUT_FILE
        printf " %s\n" "${PSEUDOPOT_ARR[@]}" >> $INPUT_FILE
	printf " ATOMIC_POSITIONS (angstrom)\n" >> $INPUT_FILE
	tail -n+3 $XYZ_FILE | tr -d '\r' >> $INPUT_FILE
	printf "K_POINTS (automatic)\n%s\n" "${kpoints[*]}" >> $INPUT_FILE
	less $INPUT_FILE
	cd ..
done {pf}<prefixes

