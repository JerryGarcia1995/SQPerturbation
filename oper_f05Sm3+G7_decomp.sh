#!/bin/sh

#SBATCH -p i8cpu
#SBATCH -N 8
#SBATCH -n 1024
#SBATCH -c 1

PrefixName=f05Sm3+G7
DirNameList=(-0.2 -0.5 -0.8)
DStartIndex=1
DEndIndex=${#DirNameList[@]}
CouplingConstList=(K J)
CStartIndex=1
CEndIndex=${#CouplingConstList[@]}
fupFilename=f06Eu3+_JH0.9.txt
fdownFilename=f04Pm3+_JH0.8.txt
ProgramName=SQPerturbationDecompositionAnalyzer.x

writeINPUT(){

cat > basic_properties_SQPerturbationDecompositionAnalyzer.txt << EOF
CouplingConstFilename[+SizeMinimumLimit(-1ForNoLimit)][+SizeMinimumLimit(-1ForNoLimit)] $DirName////$PrefixName.$CouplingConst.decomposition.txt	-1
IntmFilename1[orVACorFULL;+MaximumBlocks(-1ForNoLimit)] $fupFilename	5
IntmFilename2[orVACorFULL;+MaximumBlocks(-1ForNoLimit)] $fdownFilename	5
OutputFilename $DirName.$CouplingConst.decomposition.analyzation.txt
VerboseMode(Yes1No0)	0
MonitoringIndex	1
EOF

}

runPROGRAM(){

writeINPUT
srun $ProgramName

}

for DRunningIndex in $(seq $((DStartIndex-1)) 1 $((DEndIndex-1)))
 do
  DirName=$PrefixName"_"${DirNameList[$DRunningIndex]}
  for CRunningIndex in $(seq $((CStartIndex-1)) 1 $((CEndIndex-1)))
   do
    CouplingConst=${CouplingConstList[$CRunningIndex]}
    runPROGRAM
   done
 done