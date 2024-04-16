#!/bin/sh

#SBATCH -p i8cpu
#SBATCH -N 8
#SBATCH -n 1024
#SBATCH -c 1

PrefixName=f13Yb3+G6
DirNameList=(-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1.0)
DStartIndex=1
DEndIndex=${#DirNameList[@]}
CouplingConstList=(K J)
CStartIndex=1
CEndIndex=${#CouplingConstList[@]}
fupFilename=FULL
fdownFilename=f12Tm3+_JH0.5.txt
ProgramName=SQPerturbationDecompositionAnalyzer.x

writeINPUT(){

cat > basic_properties_SQPerturbationDecompositionAnalyzer.txt << EOF
CouplingConstFilename[+SizeMinimumLimit(-1ForNoLimit)] $DirName////$PrefixName.$CouplingConst.decomposition.txt	-1
IntmFilename1[orVACorFULL;+MaximumBlocks(-1ForNoLimit)] $fupFilename	-1
IntmFilename2[orVACorFULL;+MaximumBlocks(-1ForNoLimit)] $fdownFilename	-1
OutputFilename $DirName.$CouplingConst.decomposition.analyzation.txt
VerboseMode(Yes1No0)	0
MonitoringIndex	1000
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