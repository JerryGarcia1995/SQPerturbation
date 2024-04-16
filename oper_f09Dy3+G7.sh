#!/bin/sh

#SBATCH -p F144cpu
#SBATCH -N 144
#SBATCH -n 18018
#SBATCH -c 1

PrefixName=f09Dy3+G7
#DirNameList=(-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1.0)
#HtSHList=(HtSH_Oh.pfpi.-0.1.txt HtSH_Oh.pfpi.-0.2.txt HtSH_Oh.pfpi.-0.3.txt HtSH_Oh.pfpi.-0.4.txt HtSH_Oh.pfpi.-0.5.txt HtSH_Oh.pfpi.-0.6.txt HtSH_Oh.pfpi.-0.7.txt HtSH_Oh.pfpi.-0.8.txt HtSH_Oh.pfpi.-0.9.txt HtSH_Oh.pfpi.-1.0.txt)
DirNameList=(-0.5)
HtSHList=(HtSH_Oh.pfpi.-0.5.txt)
StartIndex=1
EndIndex=${#HtSHList[@]}

copyFILES(){

ManualMPIFilename=ManualMPI_f09.tsv
StartFilename=f09Dy3+_JH1.1.dat
fupFilename=f10Ho3+_JH0.8.dat
fdownFilename=f08Tb3+_JH1.0.dat
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
OutputFilename	f09Dy3+G7										
OutputPrecision	20										
ChopLevel	12										
ReadPrtb(Yes1No0)	1										
WritePrtb(Yes1No0)	1										
NOrbitalRooms	28										
HoppintHamiltonianFilename(IfVaried,"VA")	$HtSHFilename										
ManualMPI(Yes1No0)	1										
ManualMPIFilename(WorksForManualMPI=1Only;TSV)	ManualMPI_f09.tsv										
SaveMemory(Yes1No0)	1										
GetLabel(Yes1No0;WorksForSaveMemory=1Only)	1										
ReadSHorCH(SH1CH0)	1										
MonitoringIndex	1										
FullMonitoring(Super2Yes1No0)	2										
*****************************************************INITIAL***************************************************											
NInitialStateFile	2										
HoppintHamiltonianFilename	$HtSHFilename										
===========================================Filename&NOrbital&NOccupied=========================================											
f09Dy3+_JH1.1.dat	14	9									
-------------------------------------------NKets,NKets*(Re/Im/L/S/J/Jz)----------------------------------------											
NPickupState	2										
NKets=	4										
Re=	0	Im=	0.239356776939085	L=	5	S=	2.5	J=	7.5	Jz=	-5.5
Re=	0	Im=	0.450693909432999	L=	5	S=	2.5	J=	7.5	Jz=	-1.5
Re=	0	Im=	-0.581843335157039	L=	5	S=	2.5	J=	7.5	Jz=	2.5
Re=	0	Im=	-0.633278506398778	L=	5	S=	2.5	J=	7.5	Jz=	6.5
NKets=	4										
Re=	0	Im=	-0.239356776939085	L=	5	S=	2.5	J=	7.5	Jz=	5.5
Re=	0	Im=	-0.450693909432999	L=	5	S=	2.5	J=	7.5	Jz=	1.5
Re=	0	Im=	0.581843335157039	L=	5	S=	2.5	J=	7.5	Jz=	-2.5
Re=	0	Im=	0.633278506398778	L=	5	S=	2.5	J=	7.5	Jz=	-6.5
---------------------------------------------------------------------------------------------------------------											
f09Dy3+_JH1.1.dat	14	9									
-------------------------------------------NKets,NKets*(Re/Im/L/S/J/Jz)----------------------------------------											
NPickupState	2										
NKets=	4										
Re=	0	Im=	0.239356776939085	L=	5	S=	2.5	J=	7.5	Jz=	-5.5
Re=	0	Im=	0.450693909432999	L=	5	S=	2.5	J=	7.5	Jz=	-1.5
Re=	0	Im=	-0.581843335157039	L=	5	S=	2.5	J=	7.5	Jz=	2.5
Re=	0	Im=	-0.633278506398778	L=	5	S=	2.5	J=	7.5	Jz=	6.5
NKets=	4										
Re=	0	Im=	-0.239356776939085	L=	5	S=	2.5	J=	7.5	Jz=	5.5
Re=	0	Im=	-0.450693909432999	L=	5	S=	2.5	J=	7.5	Jz=	1.5
Re=	0	Im=	0.581843335157039	L=	5	S=	2.5	J=	7.5	Jz=	-2.5
Re=	0	Im=	0.633278506398778	L=	5	S=	2.5	J=	7.5	Jz=	-6.5
---------------------------------------------------------------------------------------------------------------											
***************************************************************************************************************											
***************************************************INTERMEDIATE************************************************											
NIntermediates	1										
===================================================INTERMEDIATE1===============================================											
NIntermediateStateGroup	2										
HoppintHamiltonianFilename	$HtSHFilename										
-------------Filename(orVAC/FULL)&NOrbital&NOccupied&ReadStart(IfAll,-1)&ReadEnd(IfAll,-1)&EnergyOffset-------------											
NIntermediateStateFile	2										
f10Ho3+_JH0.8.dat	14	10	-1	-1	23.01081246392870						
f08Tb3+_JH1.0.dat	14	8	-1	-1	26.27596786450340						
-------------Filename(orVAC/FULL)&NOrbital&NOccupied&ReadStart(IfAll,-1)&ReadEnd(IfAll,-1)&EnergyOffset-------------											
NIntermediateStateFile	2										
f08Tb3+_JH1.0.dat	14	8	-1	-1	26.27596786450340						
f10Ho3+_JH0.8.dat	14	10	-1	-1	23.01081246392870						
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