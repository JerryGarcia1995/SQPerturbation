#!/bin/sh

#SBATCH -p F16cpu
#SBATCH -N 8
#SBATCH -n 1010
#SBATCH -c 1

PrefixName=f11Er3+G7
DirNameList=(-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1.0)
HtSHList=(HtSH_Oh.pfpi.-0.1.txt HtSH_Oh.pfpi.-0.2.txt HtSH_Oh.pfpi.-0.3.txt HtSH_Oh.pfpi.-0.4.txt HtSH_Oh.pfpi.-0.5.txt HtSH_Oh.pfpi.-0.6.txt HtSH_Oh.pfpi.-0.7.txt HtSH_Oh.pfpi.-0.8.txt HtSH_Oh.pfpi.-0.9.txt HtSH_Oh.pfpi.-1.0.txt)
StartIndex=1
EndIndex=${#HtSHList[@]}

copyFILES(){

ManualMPIFilename=ManualMPI_f11.tsv
StartFilename=f11Er3+_JH0.9.dat
fupFilename=f12Tm3+_JH0.5.dat
fdownFilename=f10Ho3+_JH0.8.dat
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
OutputFilename	f11Er3+G7										
OutputPrecision	20										
ChopLevel	12										
ReadPrtb(Yes1No0)	1										
WritePrtb(Yes1No0)	1										
NOrbitalRooms	28										
HoppintHamiltonianFilename(IfVaried,"VA")	$HtSHFilename										
ManualMPI(Yes1No0)	1										
ManualMPIFilename(WorksForManualMPI=1Only;TSV)	ManualMPI_f11.tsv										
SaveMemory(Yes1No0)	1										
GetLabel(Yes1No0;WorksForSaveMemory=1Only)	1										
ReadSHorCH(SH1CH0)	1										
MonitoringIndex	1										
FullMonitoring(Super2Yes1No0)	2										
*****************************************************INITIAL***************************************************											
NInitialStateFile	2										
HoppintHamiltonianFilename	$HtSHFilename										
===========================================Filename&NOrbital&NOccupied=========================================											
f11Er3+_JH0.9.dat	14	11									
-------------------------------------------NKets,NKets*(Re/Im/L/S/J/Jz)----------------------------------------											
NPickupState	2										
NKets=	4										
Re=	0	Im=	0.239356776939085	L=	6	S=	1.5	J=	7.5	Jz=	-5.5
Re=	0	Im=	0.450693909432999	L=	6	S=	1.5	J=	7.5	Jz=	-1.5
Re=	0	Im=	-0.581843335157039	L=	6	S=	1.5	J=	7.5	Jz=	2.5
Re=	0	Im=	-0.633278506398778	L=	6	S=	1.5	J=	7.5	Jz=	6.5
NKets=	4										
Re=	0	Im=	-0.239356776939085	L=	6	S=	1.5	J=	7.5	Jz=	5.5
Re=	0	Im=	-0.450693909432999	L=	6	S=	1.5	J=	7.5	Jz=	1.5
Re=	0	Im=	0.581843335157039	L=	6	S=	1.5	J=	7.5	Jz=	-2.5
Re=	0	Im=	0.633278506398778	L=	6	S=	1.5	J=	7.5	Jz=	-6.5
---------------------------------------------------------------------------------------------------------------											
f11Er3+_JH0.9.dat	14	11									
-------------------------------------------NKets,NKets*(Re/Im/L/S/J/Jz)----------------------------------------											
NPickupState	2										
NKets=	4										
Re=	0	Im=	0.239356776939085	L=	6	S=	1.5	J=	7.5	Jz=	-5.5
Re=	0	Im=	0.450693909432999	L=	6	S=	1.5	J=	7.5	Jz=	-1.5
Re=	0	Im=	-0.581843335157039	L=	6	S=	1.5	J=	7.5	Jz=	2.5
Re=	0	Im=	-0.633278506398778	L=	6	S=	1.5	J=	7.5	Jz=	6.5
NKets=	4										
Re=	0	Im=	-0.239356776939085	L=	6	S=	1.5	J=	7.5	Jz=	5.5
Re=	0	Im=	-0.450693909432999	L=	6	S=	1.5	J=	7.5	Jz=	1.5
Re=	0	Im=	0.581843335157039	L=	6	S=	1.5	J=	7.5	Jz=	-2.5
Re=	0	Im=	0.633278506398778	L=	6	S=	1.5	J=	7.5	Jz=	-6.5
---------------------------------------------------------------------------------------------------------------											
***************************************************************************************************************											
***************************************************INTERMEDIATE************************************************											
NIntermediates	1										
===================================================INTERMEDIATE1===============================================											
NIntermediateStateGroup	2										
HoppintHamiltonianFilename	$HtSHFilename										
-------------Filename(orVAC/FULL)&NOrbital&NOccupied&ReadStart(IfAll,-1)&ReadEnd(IfAll,-1)&EnergyOffset-------------											
NIntermediateStateFile	2										
f12Tm3+_JH0.5.dat	14	12	-1	-1	18.60329259035310						
f10Ho3+_JH0.8.dat	14	10	-1	-1	27.41081246392870						
-------------Filename(orVAC/FULL)&NOrbital&NOccupied&ReadStart(IfAll,-1)&ReadEnd(IfAll,-1)&EnergyOffset-------------											
NIntermediateStateFile	2										
f10Ho3+_JH0.8.dat	14	10	-1	-1	27.41081246392870						
f12Tm3+_JH0.5.dat	14	12	-1	-1	18.60329259035310						
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