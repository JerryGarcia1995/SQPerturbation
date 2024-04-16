#!/bin/sh

#SBATCH -p F16cpu
#SBATCH -N 8
#SBATCH -n 1010
#SBATCH -c 1

PrefixName=f03Nd3+G6
DirNameList=(-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1.0)
HtSHList=(HtSH_Oh.pfpi.-0.1.txt HtSH_Oh.pfpi.-0.2.txt HtSH_Oh.pfpi.-0.3.txt HtSH_Oh.pfpi.-0.4.txt HtSH_Oh.pfpi.-0.5.txt HtSH_Oh.pfpi.-0.6.txt HtSH_Oh.pfpi.-0.7.txt HtSH_Oh.pfpi.-0.8.txt HtSH_Oh.pfpi.-0.9.txt HtSH_Oh.pfpi.-1.0.txt)
StartIndex=1
EndIndex=${#HtSHList[@]}

copyFILES(){

ManualMPIFilename=ManualMPI_f03.tsv
StartFilename=f03Nd3+_JH0.7.dat
fupFilename=f04Pm3+_JH0.8.dat
fdownFilename=f02Pr3+_JH0.6.dat
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
OutputFilename	f03Nd3+G6										
OutputPrecision	20										
ChopLevel	12										
ReadPrtb(Yes1No0)	1										
WritePrtb(Yes1No0)	1										
NOrbitalRooms	28										
HoppintHamiltonianFilename(IfVaried,"VA")	$HtSHFilename										
ManualMPI(Yes1No0)	1										
ManualMPIFilename(WorksForManualMPI=1Only;TSV)	ManualMPI_f03.tsv										
SaveMemory(Yes1No0)	1										
GetLabel(Yes1No0;WorksForSaveMemory=1Only)	1										
ReadSHorCH(SH1CH0)	1										
MonitoringIndex	1										
FullMonitoring(Super2Yes1No0)	2										
*****************************************************INITIAL***************************************************											
NInitialStateFile	2										
HoppintHamiltonianFilename	$HtSHFilename										
===========================================Filename&NOrbital&NOccupied=========================================											
f03Nd3+_JH0.7.dat	14	3									
-------------------------------------------NKets,NKets*(Re/Im/L/S/J/Jz)----------------------------------------											
NPickupState	2										
NKets=	3										
Re=	0	Im=	0.204124145231931	L=	6	S=	1.5	J=	4.5	Jz=	-3.5
Re=	0	Im=	0.763762615825973	L=	6	S=	1.5	J=	4.5	Jz=	0.5
Re=	0	Im=	0.612372435695794	L=	6	S=	1.5	J=	4.5	Jz=	4.5
NKets=	3										
Re=	0	Im=	0.204124145231931	L=	6	S=	1.5	J=	4.5	Jz=	3.5
Re=	0	Im=	0.763762615825973	L=	6	S=	1.5	J=	4.5	Jz=	-0.5
Re=	0	Im=	0.612372435695794	L=	6	S=	1.5	J=	4.5	Jz=	-4.5
---------------------------------------------------------------------------------------------------------------											
f03Nd3+_JH0.7.dat	14	3									
-------------------------------------------NKets,NKets*(Re/Im/L/S/J/Jz)----------------------------------------											
NPickupState	2										
NKets=	3										
Re=	0	Im=	0.204124145231931	L=	6	S=	1.5	J=	4.5	Jz=	-3.5
Re=	0	Im=	0.763762615825973	L=	6	S=	1.5	J=	4.5	Jz=	0.5
Re=	0	Im=	0.612372435695794	L=	6	S=	1.5	J=	4.5	Jz=	4.5
NKets=	3										
Re=	0	Im=	0.204124145231931	L=	6	S=	1.5	J=	4.5	Jz=	3.5
Re=	0	Im=	0.763762615825973	L=	6	S=	1.5	J=	4.5	Jz=	-0.5
Re=	0	Im=	0.612372435695794	L=	6	S=	1.5	J=	4.5	Jz=	-4.5
---------------------------------------------------------------------------------------------------------------											
***************************************************************************************************************											
***************************************************INTERMEDIATE************************************************											
NIntermediates	1										
===================================================INTERMEDIATE1===============================================											
NIntermediateStateGroup	2										
HoppintHamiltonianFilename	$HtSHFilename										
-------------Filename(orVAC/FULL)&NOrbital&NOccupied&ReadStart(IfAll,-1)&ReadEnd(IfAll,-1)&EnergyOffset-------------											
NIntermediateStateFile	2										
f04Pm3+_JH0.8.dat	14	4	-1	-1	8.22792700703035						
f02Pr3+_JH0.6.dat	14	2	-1	-1	6.59871471159652						
-------------Filename(orVAC/FULL)&NOrbital&NOccupied&ReadStart(IfAll,-1)&ReadEnd(IfAll,-1)&EnergyOffset-------------											
NIntermediateStateFile	2										
f02Pr3+_JH0.6.dat	14	2	-1	-1	6.59871471159652						
f04Pm3+_JH0.8.dat	14	4	-1	-1	8.22792700703035						
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