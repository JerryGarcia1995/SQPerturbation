#!/bin/sh

#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -c 1

#nfList=(1 2 3 4 5 6 7 8 9 10 11 12 13)
nfList=(12)
StartIndex=1
EndIndex=${#nfList[@]}

LSjob(){

cat > basic_properties_LSCoupling.txt << EOF
NOrbitalRoom(e.g.,d:10,f:14)	14
NOirbtalOccupied	$nf
CubicHarmoicOperators	operatorfzetaDO.dat	operatorfzetaUP.dat	operatorfetaDO.dat	operatorfetaUP.dat	operatorfxiDO.dat	operatorfxiUP.dat	operatorfADO.dat	operatorfAUP.dat	operatorfgammaDO.dat	operatorfgammaUP.dat	operatorfbetaDO.dat	operatorfbetaUP.dat	operatorfalphaDO.dat	operatorfalphaUP.dat
ChopLevel	12
OutputFilePrefix(IfWantAutomode,"Auto";Then,OutputFilePrefix,SLlist,andConditions_will_be_automated)	Auto
EigenvectorOutputMode(HumanReadable0Binary1)	1
SLlist(Recipe)	SLlist1.ini
operatorSx	operatorSx.dat
operatorSy	operatorSy.dat
operatorSz	operatorSz.dat
operatorSm	operatorSm.dat
operatorSp	operatorSp.dat
operatorLx	operatorLx.dat
operatorLy	operatorLy.dat
operatorLz	operatorLz.dat
operatorLm	operatorLm.dat
operatorLp	operatorLp.dat
operatorSOCtraverse	operatorSOCtraverse.dat
operatorSOCvertical	operatorSOCvertical.dat
operatorSOC	operatorSOC.dat
NCoulombHamiltonian	4
F0	CoulombHamiltonianF0.dat
F2	CoulombHamiltonianF2.dat
F4	CoulombHamiltonianF4.dat
F6	CoulombHamiltonianF6.dat
SOC_coefficient_name	Lambda
NConditions 16
Index	ConditionName	F0	F2	F4	F6	Lambda (UnderAutomode,"InputIndex&ConditionName&JH"Only;F0_will_be_0)
1 JH0.1	0.1
2 JH0.2	0.2
3 JH0.3	0.3
4 JH0.4	0.4
5 JH0.5	0.5
6 JH0.6	0.6
7 JH0.7	0.7
8 JH0.8	0.8
9 JH0.9	0.9
10 JH1.0	1.0
11 JH1.1	1.1
12 JH1.2	1.2
13 JH1.3	1.3
14 JH1.4	1.4
15 JH1.5	1.5
16 JH1.6	1.6

****Example(s) for Automode****
1 JH0.1	0.1
2 JH0.2	0.2
3 JH0.3	0.3
4 JH0.4	0.4
5 JH0.5	0.5
6 JH0.6	0.6
7 JH0.7	0.7
8 JH0.8	0.8
9 JH0.9	0.9
10 JH1.0	1.0
11 JH1.1	1.1
12 JH1.2	1.2
13 JH1.3	1.3
14 JH1.4	1.4
15 JH1.5	1.5
16 JH1.6	1.6

***Recommended JH***
f01Ce3+ Any
f02Pr3+	0.7
f03Nd3+	0.7
f04Pm3+	0.8
f05Sm3+	0.7
f06Eu3+	1.0
f07Gd3+	0.9
f08Tb3+	1.1
f09Dy3+	0.9
f10Ho3+	0.8
f11Er3+	0.8
f12Tm3+	1.6
f13Yb3+ Any

****Example(s) for not-Automode****
17	f01U01JH0.16	0	1.97473538	1.238023947	0.893644046	0.0996
EOF

srun LSCoupling.x

}

for RunningIndex in $(seq $((StartIndex-1)) 1 $((EndIndex-1)))
 do
  nf=${nfList[$RunningIndex]}
  LSjob
 done
 
