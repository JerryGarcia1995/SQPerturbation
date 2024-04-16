#!/bin/sh

#SBATCH -p i8cpu
#SBATCH -N 2
#SBATCH -n 182
#SBATCH -c 1

PrefixName=f13Yb3+G6
DirNameList=(-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1.0)
HtSHList=(HtSH_Oh.pfpi.-0.1.txt HtSH_Oh.pfpi.-0.2.txt HtSH_Oh.pfpi.-0.3.txt HtSH_Oh.pfpi.-0.4.txt HtSH_Oh.pfpi.-0.5.txt HtSH_Oh.pfpi.-0.6.txt HtSH_Oh.pfpi.-0.7.txt HtSH_Oh.pfpi.-0.8.txt HtSH_Oh.pfpi.-0.9.txt HtSH_Oh.pfpi.-1.0.txt)
StartIndex=1
EndIndex=${#HtSHList[@]}

copyFILES(){

ManualMPIFilename=ManualMPI_f13.tsv
StartFilename=f13Yb3+_JH0.1.dat
fupFilename=FULL
fdownFilename=f12Tm3+_JH0.5.dat
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
OutputFilename	f13Yb3+G6										
OutputPrecision	20										
ChopLevel	12										
ReadPrtb(Yes1No0)	1										
WritePrtb(Yes1No0)	1										
NOrbitalRooms	28										
HoppintHamiltonianFilename(IfVaried,"VA")	$HtSHFilename										
ManualMPI(Yes1No0)	1										
ManualMPIFilename(WorksForManualMPI=1Only;TSV)	ManualMPI_f13.tsv										
SaveMemory(Yes1No0)	1										
GetLabel(Yes1No0;WorksForSaveMemory=1Only)	1										
ReadSHorCH(SH1CH0)	1										
MonitoringIndex	1										
FullMonitoring(Super2Yes1No0)	2										
*****************************************************INITIAL***************************************************											
NInitialStateFile	2										
HoppintHamiltonianFilename	$HtSHFilename										
===========================================Filename&NOrbital&NOccupied=========================================											
f13Yb3+_JH0.1.dat	14	13									
-------------------------------------------NKets,NKets*(Re/Im/L/S/J/Jz)----------------------------------------											
NPickupState	2										
NKets=	2										
Re=	0	Im=	-0.645497224367903	L=	3	S=	0.5	J=	3.5	Jz=	-3.5
Re=	0	Im=	-0.763762615825973	L=	3	S=	0.5	J=	3.5	Jz=	0.5
NKets=	2										
Re=	0	Im=	0.645497224367903	L=	3	S=	0.5	J=	3.5	Jz=	3.5
Re=	0	Im=	0.763762615825973	L=	3	S=	0.5	J=	3.5	Jz=	-0.5
---------------------------------------------------------------------------------------------------------------											
f13Yb3+_JH0.1.dat	14	13									
-------------------------------------------NKets,NKets*(Re/Im/L/S/J/Jz)----------------------------------------											
NPickupState	2										
NKets=	2										
Re=	0	Im=	-0.645497224367903	L=	3	S=	0.5	J=	3.5	Jz=	-3.5
Re=	0	Im=	-0.763762615825973	L=	3	S=	0.5	J=	3.5	Jz=	0.5
NKets=	2										
Re=	0	Im=	0.645497224367903	L=	3	S=	0.5	J=	3.5	Jz=	3.5
Re=	0	Im=	0.763762615825973	L=	3	S=	0.5	J=	3.5	Jz=	-0.5
---------------------------------------------------------------------------------------------------------------											
***************************************************************************************************************											
***************************************************INTERMEDIATE************************************************											
NIntermediates	1										
===================================================INTERMEDIATE1===============================================											
NIntermediateStateGroup	2										
HoppintHamiltonianFilename	$HtSHFilename										
-------------Filename(orVAC/FULL)&NOrbital&NOccupied&ReadStart(IfAll,-1)&ReadEnd(IfAll,-1)&EnergyOffset-------------											
NIntermediateStateFile	2										
FULL	14	14	-1	-1	6.40000000000000						
f12Tm3+_JH0.5.dat	14	12	-1	-1	24.20329259035310						
-------------Filename(orVAC/FULL)&NOrbital&NOccupied&ReadStart(IfAll,-1)&ReadEnd(IfAll,-1)&EnergyOffset-------------											
NIntermediateStateFile	2										
f12Tm3+_JH0.5.dat	14	12	-1	-1	24.20329259035310						
FULL	14	14	-1	-1	6.40000000000000						
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