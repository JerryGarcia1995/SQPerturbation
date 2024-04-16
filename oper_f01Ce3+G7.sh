#!/bin/sh

#SBATCH -p i8cpu
#SBATCH -N 2
#SBATCH -n 182
#SBATCH -c 1

PrefixName=f01Ce3+G7
DirNameList=(-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1.0)
HtSHList=(HtSH_Oh.pfpi.-0.1.txt HtSH_Oh.pfpi.-0.2.txt HtSH_Oh.pfpi.-0.3.txt HtSH_Oh.pfpi.-0.4.txt HtSH_Oh.pfpi.-0.5.txt HtSH_Oh.pfpi.-0.6.txt HtSH_Oh.pfpi.-0.7.txt HtSH_Oh.pfpi.-0.8.txt HtSH_Oh.pfpi.-0.9.txt HtSH_Oh.pfpi.-1.0.txt)
StartIndex=1
EndIndex=${#HtSHList[@]}

copyFILES(){

ManualMPIFilename=ManualMPI_f01.tsv
StartFilename=f01Ce3+_JH0.1.dat
fupFilename=f02Pr3+_JH0.6.dat
fdownFilename=VAC
ProgramName=SQPerturbation.x

cp $ManualMPIFilename $DirName/.
cp $StartFilename $DirName/.
cp $fupFilename $DirName/.
cp $fdownFilename $DirName/.
cp $HtSHFilename $DirName/.
cp $ProgramName $DirName/.

}

writeINPUT(){

cat > basic_properties_SQPerturbation.txt << EOF
OutputFilename	f01Ce3+G7										
OutputPrecision	20										
ChopLevel	12										
ReadPrtb(Yes1No0)	1										
WritePrtb(Yes1No0)	1										
NOrbitalRooms	28										
HoppintHamiltonianFilename(IfVaried,"VA")	$HtSHFilename										
ManualMPI(Yes1No0)	1										
ManualMPIFilename(WorksForManualMPI=1Only;TSV)	ManualMPI_f01.tsv										
SaveMemory(Yes1No0)	1										
GetLabel(Yes1No0;WorksForSaveMemory=1Only)	1										
ReadSHorCH(SH1CH0)	1										
MonitoringIndex	1										
FullMonitoring(Super2Yes1No0)	2										
*****************************************************INITIAL***************************************************											
NInitialStateFile	2										
HoppintHamiltonianFilename	$HtSHFilename										
===========================================Filename&NOrbital&NOccupied=========================================											
f01Ce3+_JH0.1.dat	14	1									
-------------------------------------------NKets,NKets*(Re/Im/L/S/J/Jz)----------------------------------------											
NPickupState	2										
NKets=	2										
Re=	0	Im=	-0.912870929175277	L=	3	S=	0.5	J=	2.5	Jz=	-1.5
Re=	0	Im=	0.408248290463863	L=	3	S=	0.5	J=	2.5	Jz=	2.5
NKets=	2										
Re=	0	Im=	-0.912870929175277	L=	3	S=	0.5	J=	2.5	Jz=	1.5
Re=	0	Im=	0.408248290463863	L=	3	S=	0.5	J=	2.5	Jz=	-2.5
---------------------------------------------------------------------------------------------------------------											
f01Ce3+_JH0.1.dat	14	1									
-------------------------------------------NKets,NKets*(Re/Im/L/S/J/Jz)----------------------------------------											
NPickupState	2										
NKets=	2										
Re=	0	Im=	-0.912870929175277	L=	3	S=	0.5	J=	2.5	Jz=	-1.5
Re=	0	Im=	0.408248290463863	L=	3	S=	0.5	J=	2.5	Jz=	2.5
NKets=	2										
Re=	0	Im=	-0.912870929175277	L=	3	S=	0.5	J=	2.5	Jz=	1.5
Re=	0	Im=	0.408248290463863	L=	3	S=	0.5	J=	2.5	Jz=	-2.5
---------------------------------------------------------------------------------------------------------------											
***************************************************************************************************************											
***************************************************INTERMEDIATE************************************************											
NIntermediates	1										
===================================================INTERMEDIATE1===============================================											
NIntermediateStateGroup	2										
HoppintHamiltonianFilename	$HtSHFilename										
-------------Filename(orVAC/FULL)&NOrbital&NOccupied&ReadStart(IfAll,-1)&ReadEnd(IfAll,-1)&EnergyOffset-------------											
NIntermediateStateFile	2										
f02Pr3+_JH0.6.dat	14	2	-1	-1	4.59871471159652						
VAC	14	0	-1	-1	1.80000000000000						
-------------Filename(orVAC/FULL)&NOrbital&NOccupied&ReadStart(IfAll,-1)&ReadEnd(IfAll,-1)&EnergyOffset-------------											
NIntermediateStateFile	2										
VAC	14	0	-1	-1	1.80000000000000						
f02Pr3+_JH0.6.dat	14	2	-1	-1	4.59871471159652						
***************************************************************************************************************											
*****************************************************FINAL*****************************************************											
RE											
***************************************************************************************************************											
NUserdefinedTerms	4										
TermName	J										
NTermComposites	2										
////Index_Coefficient_i_j_reim////											
1	2	1	4	re							
2	2	2	3	re							
TermName	K										
NTermComposites	4										
////Index_Coefficient_i_j_reim////											
1	2	1	1	re							
2	-2	2	2	re							
3	-2	1	4	re							
4	-2	2	3	re							
TermName	G										
NTermComposites	2										
////Index_Coefficient_i_j_reim////											
1	2	2	3	im							
2	-2	1	4	im							
TermName	Gp										
NTermComposites	1										
////Index_Coefficient_i_j_reim////											
1	4	1	2	re							
											
EOF

}

runPROGRAM(){

mkdir $DirName
copyFILES
cd $DirName
writeINPUT
srun SQPerturbation.x
cd ..

}

for RunningIndex in $(seq $((StartIndex-1)) 1 $((EndIndex-1)))
 do
  DirName=$PrefixName"_"${DirNameList[$RunningIndex]}
  HtSHFilename=${HtSHList[$RunningIndex]}
  runPROGRAM
 done