#QERUN.SH 2017 by Michael Lee
#TURNS XYZ FILES INTO INPUT FILES & EXECUTES (OPTIONAL)

update_elem () { #update this function depending on elements used
	case $1 in
		'H')
			val=1;;
		'C')
			pseudopot="C 12.0107 C.pbe-n-rrkjus_psl.0.1.UPF"
			val=4;;
		'O')
			val=6;;
		'F')
			pseudopot="F 18.99840325 F.pbe-n-rrkjus_psl.0.1.UPF"
			val=7;;
		'Pt')
			val=10;;
		*)
			echo "Invalid Element" $1
	esac
}

rm_filelist (){
	rm NAT_TYP_ARR.tmp
}

while read prefix
do
	INPUT_FILE=$(printf $prefix".in") #pipe all to this file
	>$INPUT_FILE
	XYZ=$(printf $prefix".xyz")
	tail -n+3 $XYZ | cut -c1-2 | sort | uniq -c > NAT_TYP_ARR.tmp #array of elements
	TYP_ARR=($(cut -c9-10 NAT_TYP_ARR.tmp))  #array of elements
	NATOM_ARR=($(cut -c6-7 NAT_TYP_ARR.tmp))  #array of number of atoms
	NTYP=${#TYP_ARR[@]} #number of unique elements
	NATOM=$(head -n1 $XYZ | tr -d '\r') #can check with sum atoms
	NBND=0 
	PSEUDO_FILE="pseudo-template"
	elem_index=0
	for i in "${TYP_ARR[@]}"
	do
		#PSEUDO_POTS[$elem_index]=$(grep -w $i $PSEUDO_FILE) #match whole words (not 'F' in 'UPF')
		update_elem $i
		NBND=$(($NBND+(${NATOM_ARR[$elem_index]}*val+1)/2))
		elem_index=$(($elem_index+1))
		PSEUDO_POTS[$elem_index]=$pseudopot #match whole words (not 'F' in 'UPF')
	done
	NBND=$(($NBND+10)) #buffer NBND by 10

#SUBSTITUTION OF VARIABLES INTO DEFAULT INPUT FILE
#TO INSERT NEW VARIABLE, BREAK PRINTF STREAM & PRINTF AT DESIRED LINE
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
    ecutwfc = 40.0
    ecutrho = 480.0
    nat = %i
    ntyp = %i 
    nbnd = %i
    occupations = 'smearing'
    smearing = 'cold'
    degauss = 0.001
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
" $prefix $prefix $NATOM $NTYP $NBND >> $INPUT_FILE #print main parameters
	printf "%s\n" "${PSEUDO_POTS[@]}" >> $INPUT_FILE #print pseudopotentials
	printf "ATOMIC_POSITIONS (angstrom)\n" >> $INPUT_FILE
	tail -n+3 $XYZ | tr -d '\r' >> $INPUT_FILE #print coordinates, fixes non-unix endline escape key
	printf "K_POINTS (automatic) \n18 18 1 0 0 0" >> $INPUT_FILE

rm_filelist
done < prefixes
