#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include <set>
#include "Eigen/Dense"
#include <mpi.h>

using namespace std;
using namespace Eigen;

typedef struct {

	string condition_name;
	vector<long double> CoulombCoeff;
	long double SOCCoeff;

}CONDITIONS;

vector<CONDITIONS> Conditions;

typedef struct {

	long double coeffcientRe;
	long double coeffcientIm;
	vector<int> CRList;
	vector<int> ANList;

}OPERATORLINE;

typedef struct {

	vector<OPERATORLINE> OperatorLine;

}OPERATOR;

typedef struct {

	long double S;
	vector<long double> L;

}SLLISTLINE;

typedef struct {

	vector<SLLISTLINE> SLListLine;

}SLLIST;

SLLIST _SLlist;

typedef struct {

	long double coeffcientRe;
	long double coeffcientIm;
	signed long long int base;

}KETELEMENT;

typedef struct {

	vector<KETELEMENT> KetElement;

}KET;

typedef struct {

	string CoulombName, CoulombFile;
	long double CoulombCoeff;
	OPERATOR CoulombHamiltonian;

}COULOMB;

vector<COULOMB> Coulomb;

vector<OPERATOR> operatorCubicHarmonic;

typedef struct {

	long double Lx, Ly, Lz;
	long double Sx, Sy, Sz;
	signed long long int ket;

}LISTREPETITIVENC;

vector<LISTREPETITIVENC> listRepetitiveNC;

typedef struct {

	signed long long int starket;
	string starket_rep;
	KET ket;

}LISTREPETITIVENCCUBICHARMONIC;

vector<LISTREPETITIVENCCUBICHARMONIC> listRepetitiveNCCubicHarmonic;

typedef struct {

	long double L, Lx, Ly, Lz;
	long double S, Sx, Sy, Sz;
	long double *CoulombEnerg;
	long double CoulombEnergySum;
	KET ket;

}LSKETLIST;

vector<LSKETLIST> LSKetList;

typedef struct {

	long double L, Lx, Ly, Lz;
	long double S, Sx, Sy, Sz;
	long double J, Jx, Jy, Jz;
	long double *CoulombEnerg;
	long double CoulombEnergySum;
	long double SOCEnerg;
	long double CoulombAndSOCEnerg;
	KET ket;

}LSJKETLIST;

vector<LSJKETLIST> LSJKetList;

typedef struct {

	signed long long int istart;
	signed long long int iend;
	int multiplicity;
	long double L, S;

}LSBLOCK;

vector<LSBLOCK> LSBlock;

typedef struct {

	signed long long int istart;
	signed long long int iend;
	long double L, S, J;

}LSJBLOCK;

vector<LSJBLOCK> LSJBlock;

int pr_size, ierr, pr_id;
int *pr_job_count;
int *pr_iconditionstart, *pr_iconditionend;

OPERATOR OperatorSx;
OPERATOR OperatorSy;
OPERATOR OperatorSz;
OPERATOR OperatorSm;
OPERATOR OperatorSp;
OPERATOR OperatorLx;
OPERATOR OperatorLy;
OPERATOR OperatorLz;
OPERATOR OperatorLm;
OPERATOR OperatorLp;
OPERATOR OperatorJm;
OPERATOR OperatorJp;
OPERATOR OperatorSOCtraverse;
OPERATOR OperatorSOCvertical;
OPERATOR OperatorSOC;

int NOrbitalRoom;
int NOirbtalOccupied;
vector<string> CubicHarmoicOperators;
long double ChopLevel, Chop;
string OutputFilePrefix;
int EigenvectorOutputMode;
string SLlist;
string operatorSx;
string operatorSy;
string operatorSz;
string operatorSm;
string operatorSp;
string operatorLx;
string operatorLy;
string operatorLz;
string operatorLm;
string operatorLp;
string operatorSOCtraverse;
string operatorSOCvertical;
string operatorSOC;
int NCoulombHamiltonian;
string SOC_coefficient_name;
int NConditions;
long double SOC_coefficient;

ifstream readpara;
ofstream writeoutput;
FILE *_writeoutput;
string _writeoutput_filename;

string sbuf;
int CoulombHamiltonian_count;
int CubicHarmonic_count;
int *SecondQuantizedKet;
int *ConvertedQuantizedKet;

KET vac;

KET InputKet;
KET OutputKet;
KET SimplifiedOutputKet;
KET SSimplifiedOutputKet;
KET NormalizedSimplifiedOutputKet;
KET InputBra;
vector<KET> GramSchmidtKet;

bool Automode;
long double operate_scalar_re, operate_scalar_im;
signed long long int initial_mLindex;
signed long long int LSKET_size_buf;
int iGramSchmidt;
vector<int> same_L_series;
int i_same_L_series;
vector<signed long long int> catched;
signed long long int global_istart, global_iend;

struct sort_class {

	bool operator() (long double i, long double j) {

		return (i > j);

	}

} sort_object;

signed long long int ConvertedKetNumber() {

	signed long long int Bin = 0;

	for (int i = 0; i < NOrbitalRoom; i++) {
		Bin += ConvertedQuantizedKet[NOrbitalRoom - (i + 1)] * (int)pow(2, i);
	}

	return Bin;

}

void Secondization(signed long long int ketNumber) {

	signed long long int ketNumberBin = ketNumber;

	for (int i = 0; i < NOrbitalRoom; i++) {
		SecondQuantizedKet[NOrbitalRoom - (i + 1)] = ketNumberBin % 2;
		ketNumberBin = ketNumberBin / 2;
	}

}

long double ComplexMTwoRe(double a1, double a2, double b1, double b2) {
	return a1 * b1 - a2 * b2;
}

long double ComplexMTwoIm(double a1, double a2, double b1, double b2) {
	return a2 * b1 + a1 * b2;
}

long double ComplexMThreeRe(double a1, double a2, double b1, double b2, double c1, double c2) {
	return a1 * b1 * c1 - a2 * b2 * c1 - a2 * b1 * c2 - a1 * b2 * c2;
}

long double ComplexMThreeIm(double a1, double a2, double b1, double b2, double c1, double c2) {
	return a2 * b1 * c1 + a1 * b2 * c1 + a1 * b1 * c2 - a2 * b2 * c2;
}

bool ChopBool(long double v) {

	if (abs(v) < Chop) {
		return true;
	}
	else {
		return false;
	}

}

signed long long int nCr(signed long long int n, signed long long int r) {

	signed long long int nf, nrf, rf;

	nf = 1;
	for (signed long long int inf = 1; inf <= n; inf++) {
		nf *= inf;
	}

	nrf = 1;
	for (signed long long int inrf = 1; inrf <= n - r; inrf++) {
		nrf *= inrf;
	}

	rf = 1;
	for (signed long long int irf = 1; irf <= r; irf++) {
		rf *= irf;
	}

	return nf / (nrf * rf);

}

void read_para() {

	readpara.open("basic_properties_LSCoupling.txt");

	readpara >> sbuf >> NOrbitalRoom;

	SecondQuantizedKet = new int[NOrbitalRoom];
	ConvertedQuantizedKet = new int[NOrbitalRoom];

	readpara >> sbuf >> NOirbtalOccupied;
	readpara >> sbuf;
	for (int i = 0; i < NOrbitalRoom; i++) {
		readpara >> sbuf;
		CubicHarmoicOperators.push_back(sbuf);
	}

	readpara >> sbuf >> ChopLevel;

	Chop = pow(10.0, -ChopLevel);

	readpara >> sbuf >> OutputFilePrefix;

	if (OutputFilePrefix == "Auto") {
		Automode = true;
		if (NOirbtalOccupied == 1) {
			OutputFilePrefix = "f01Ce3+";
		}
		else if (NOirbtalOccupied == 2) {
			OutputFilePrefix = "f02Pr3+";
		}
		else if (NOirbtalOccupied == 3) {
			OutputFilePrefix = "f03Nd3+";
		}
		else if (NOirbtalOccupied == 4) {
			OutputFilePrefix = "f04Pm3+";
		}
		else if (NOirbtalOccupied == 5) {
			OutputFilePrefix = "f05Sm3+";
		}
		else if (NOirbtalOccupied == 6) {
			OutputFilePrefix = "f06Eu3+";
		}
		else if (NOirbtalOccupied == 7) {
			OutputFilePrefix = "f07Gd3+";
		}
		else if (NOirbtalOccupied == 8) {
			OutputFilePrefix = "f08Tb3+";
		}
		else if (NOirbtalOccupied == 9) {
			OutputFilePrefix = "f09Dy3+";
		}
		else if (NOirbtalOccupied == 10) {
			OutputFilePrefix = "f10Ho3+";
		}
		else if (NOirbtalOccupied == 11) {
			OutputFilePrefix = "f11Er3+";
		}
		else if (NOirbtalOccupied == 12) {
			OutputFilePrefix = "f12Tm3+";
		}
		else if (NOirbtalOccupied == 13) {
			OutputFilePrefix = "f13Yb3+";
		}
	}

	readpara >> sbuf >> EigenvectorOutputMode;
	readpara >> sbuf >> SLlist;

	if (Automode == true) {
		if (NOirbtalOccupied == 1) {
			SLlist = "SLlist1.ini";
		}
		else if (NOirbtalOccupied == 2) {
			SLlist = "SLlist2.ini";
		}
		else if (NOirbtalOccupied == 3) {
			SLlist = "SLlist3.ini";
		}
		else if (NOirbtalOccupied == 4) {
			SLlist = "SLlist4.ini";
		}
		else if (NOirbtalOccupied == 5) {
			SLlist = "SLlist5.ini";
		}
		else if (NOirbtalOccupied == 6) {
			SLlist = "SLlist6.ini";
		}
		else if (NOirbtalOccupied == 7) {
			SLlist = "SLlist7.ini";
		}
		else if (NOirbtalOccupied == 8) {
			SLlist = "SLlist8.ini";
		}
		else if (NOirbtalOccupied == 9) {
			SLlist = "SLlist9.ini";
		}
		else if (NOirbtalOccupied == 10) {
			SLlist = "SLlist10.ini";
		}
		else if (NOirbtalOccupied == 11) {
			SLlist = "SLlist11.ini";
		}
		else if (NOirbtalOccupied == 12) {
			SLlist = "SLlist12.ini";
		}
		else if (NOirbtalOccupied == 13) {
			SLlist = "SLlist13.ini";
		}
	}

	readpara >> sbuf >> operatorSx;
	readpara >> sbuf >> operatorSy;
	readpara >> sbuf >> operatorSz;
	readpara >> sbuf >> operatorSm;
	readpara >> sbuf >> operatorSp;
	readpara >> sbuf >> operatorLx;
	readpara >> sbuf >> operatorLy;
	readpara >> sbuf >> operatorLz;
	readpara >> sbuf >> operatorLm;
	readpara >> sbuf >> operatorLp;
	readpara >> sbuf >> operatorSOCtraverse;
	readpara >> sbuf >> operatorSOCvertical;
	readpara >> sbuf >> operatorSOC;
	readpara >> sbuf >> NCoulombHamiltonian;

	for (int i = 0; i < NCoulombHamiltonian; i++) {
		COULOMB obj;
		readpara >> obj.CoulombName >> obj.CoulombFile;
		Coulomb.push_back(obj);
	}

	readpara >> sbuf >> SOC_coefficient_name;
	readpara >> sbuf >> NConditions;
	for (int i = 0; i < NCoulombHamiltonian + 4; i++) {
		readpara >> sbuf;
	}

	if (Automode == false) {

		for (int i = 0; i < NConditions; i++) {
			CONDITIONS obj;
			readpara >> sbuf;
			readpara >> obj.condition_name;
			for (int j = 0; j < NCoulombHamiltonian; j++) {
				long double ldbuf;
				readpara >> ldbuf;
				obj.CoulombCoeff.push_back(ldbuf);
			}
			readpara >> obj.SOCCoeff;
			Conditions.push_back(obj);
		}

	}
	else {

		long double F2whenJH1, F4whenJH1, F6whenJH1, SOCrecommended;
		//PR127(6),2058
		if (NOirbtalOccupied == 1) {
			F2whenJH1 = 12.34209613; F4whenJH1 = 7.73764967; F6whenJH1 = 5.58527529; SOCrecommended = 0.099600000;
		}
		else if (NOirbtalOccupied == 2) {
			F2whenJH1 = 12.33115193; F4whenJH1 = 7.75544146; F6whenJH1 = 5.58391785; SOCrecommended = 0.117600000;
		}
		else if (NOirbtalOccupied == 3) {
			F2whenJH1 = 12.31801444; F4whenJH1 = 7.76678019; F6whenJH1 = 5.59010294; SOCrecommended = 0.135600000;
		}
		else if (NOirbtalOccupied == 4) {
			F2whenJH1 = 12.32117776; F4whenJH1 = 7.75727499; F6whenJH1 = 5.59389816; SOCrecommended = 0.155003521;
		}
		else if (NOirbtalOccupied == 5) {
			F2whenJH1 = 12.31547934; F4whenJH1 = 7.75504808; F6whenJH1 = 5.60215414; SOCrecommended = 0.177600000;
		}
		else if (NOirbtalOccupied == 6) {
			F2whenJH1 = 12.31761977; F4whenJH1 = 7.75986412; F6whenJH1 = 5.59594898; SOCrecommended = 0.199518400;
		}
		else if (NOirbtalOccupied == 7) {
			F2whenJH1 = 12.31380817; F4whenJH1 = 7.75878886; F6whenJH1 = 5.60114814; SOCrecommended = 0.224255325;
		}
		else if (NOirbtalOccupied == 8) {
			F2whenJH1 = 12.31806858; F4whenJH1 = 7.76011200; F6whenJH1 = 5.59524219; SOCrecommended = 0.250645240;
		}
		else if (NOirbtalOccupied == 9) {
			F2whenJH1 = 12.32393862; F4whenJH1 = 7.76261669; F6whenJH1 = 5.58657320; SOCrecommended = 0.277200000;
		}
		else if (NOirbtalOccupied == 10) {
			F2whenJH1 = 12.32192444; F4whenJH1 = 7.75837829; F6whenJH1 = 5.59218337; SOCrecommended = 0.308384040;
		}
		else if (NOirbtalOccupied == 11) {
			F2whenJH1 = 12.33376168; F4whenJH1 = 7.74864573; F6whenJH1 = 5.58623297; SOCrecommended = 0.339600000;
		}
		else if (NOirbtalOccupied == 12) {
			F2whenJH1 = 12.32881596; F4whenJH1 = 7.75488784; F6whenJH1 = 5.58702202; SOCrecommended = 0.372734800;
		}
		else if (NOirbtalOccupied == 13) {
			F2whenJH1 = 12.32692411; F4whenJH1 = 7.75712981; F6whenJH1 = 5.58743756; SOCrecommended = 0.408000000;
		}
	
		for (int i = 0; i < NConditions; i++) {
			CONDITIONS obj;
			readpara >> sbuf;
			readpara >> obj.condition_name;
			long double vJH;
			readpara >> vJH;
			obj.CoulombCoeff.push_back(0.0);
			obj.CoulombCoeff.push_back(F2whenJH1 * vJH);
			obj.CoulombCoeff.push_back(F4whenJH1 * vJH);
			obj.CoulombCoeff.push_back(F6whenJH1 * vJH);
			obj.SOCCoeff = SOCrecommended;
			Conditions.push_back(obj);
		}

	}
	
	readpara.close();

}

void read_SLlist() {

	ifstream readSLlist;
	readSLlist.open(SLlist);

	while (true) {

		SLLISTLINE obj_SLLine;
		getline(readSLlist, sbuf);
		stringstream line_reader;
		line_reader << sbuf;
		line_reader >> obj_SLLine.S;

		while (true) {

			long double ldbuf;
			line_reader >> ldbuf;
			obj_SLLine.L.push_back(ldbuf);
			if (line_reader.eof()) {
				break;
			}

		}

		_SLlist.SLListLine.push_back(obj_SLLine);

		if (readSLlist.eof()) {
			break;
		}

	}

	readSLlist.close();

}

void read_operator(string operator_type) {

	string operator_filename;
	if (operator_type == "operatorSx") {
		operator_filename = operatorSx;
	}
	else if (operator_type == "operatorSy") {
		operator_filename = operatorSy;
	}
	else if (operator_type == "operatorSz") {
		operator_filename = operatorSz;
	}
	else if (operator_type == "operatorSm") {
		operator_filename = operatorSm;
	}
	else if (operator_type == "operatorSp") {
		operator_filename = operatorSp;
	}
	else if (operator_type == "operatorLx") {
		operator_filename = operatorLx;
	}
	else if (operator_type == "operatorLy") {
		operator_filename = operatorLy;
	}
	else if (operator_type == "operatorLz") {
		operator_filename = operatorLz;
	}
	else if (operator_type == "operatorLm") {
		operator_filename = operatorLm;
	}
	else if (operator_type == "operatorLp") {
		operator_filename = operatorLp;
	}
	else if (operator_type == "operatorSOCtraverse") {
		operator_filename = operatorSOCtraverse;
	}
	else if (operator_type == "operatorSOCvertical") {
		operator_filename = operatorSOCvertical;
	}
	else if (operator_type == "operatorSOC") {
		operator_filename = operatorSOC;
	}
	else if (operator_type == "Coulomb") {
		operator_filename = Coulomb[CoulombHamiltonian_count].CoulombFile;
	}
	else {
		operator_filename = CubicHarmoicOperators[CubicHarmonic_count];
	}

	ifstream readoperator;
	readoperator.open(operator_filename);

	OPERATOR obj;
	while (true) {

		OPERATORLINE obj_line;
		getline(readoperator, sbuf);
		stringstream line_reader;
		line_reader << sbuf;
		line_reader >> obj_line.coeffcientRe >> obj_line.coeffcientIm;
		for (int i = 0; i < NOrbitalRoom; i++) {
			int ibuf;
			line_reader >> ibuf;
			obj_line.CRList.push_back(ibuf);
		}
		for (int i = 0; i < NOrbitalRoom; i++) {
			int ibuf;
			line_reader >> ibuf;
			obj_line.ANList.push_back(ibuf);
		}
		obj.OperatorLine.push_back(obj_line);

		if (readoperator.eof()) {
			break;
		}

	}

	if (operator_type == "operatorSx") {
		OperatorSx = obj;
	}
	else if (operator_type == "operatorSy") {
		OperatorSy = obj;
	}
	else if (operator_type == "operatorSz") {
		OperatorSz = obj;
	}
	else if (operator_type == "operatorSm") {
		OperatorSm = obj;
	}
	else if (operator_type == "operatorSp") {
		OperatorSp = obj;
	}
	else if (operator_type == "operatorLx") {
		OperatorLx = obj;
	}
	else if (operator_type == "operatorLy") {
		OperatorLy = obj;
	}
	else if (operator_type == "operatorLz") {
		OperatorLz = obj;
	}
	else if (operator_type == "operatorLm") {
		OperatorLm = obj;
	}
	else if (operator_type == "operatorLp") {
		OperatorLp = obj;
	}
	else if (operator_type == "operatorSOCtraverse") {
		OperatorSOCtraverse = obj;
	}
	else if (operator_type == "operatorSOCvertical") {
		OperatorSOCvertical = obj;
	}
	else if (operator_type == "operatorSOC") {
		OperatorSOC = obj;
	}
	else if (operator_type == "Coulomb") {
		Coulomb[CoulombHamiltonian_count].CoulombHamiltonian = obj;
	}
	else {
		operatorCubicHarmonic.push_back(obj);
	}

	readoperator.close();

}

void read_operators() {

	read_operator("operatorSx");
	read_operator("operatorSy");
	read_operator("operatorSz");
	read_operator("operatorSm");
	read_operator("operatorSp");
	read_operator("operatorLx");
	read_operator("operatorLy");
	read_operator("operatorLz");
	read_operator("operatorLm");
	read_operator("operatorLp");
	read_operator("operatorSOCtraverse");
	read_operator("operatorSOCvertical");
	read_operator("operatorSOC");

	for (CoulombHamiltonian_count = 0; CoulombHamiltonian_count < NCoulombHamiltonian; CoulombHamiltonian_count++) {
		read_operator("Coulomb");
	}

	for (int i = 0; i < (signed)OperatorSm.OperatorLine.size(); i++) {
		OperatorJm.OperatorLine.push_back(OperatorSm.OperatorLine[i]);
	}
	for (int i = 0; i < (signed)OperatorLm.OperatorLine.size(); i++) {
		OperatorJm.OperatorLine.push_back(OperatorLm.OperatorLine[i]);
	}

	for (int i = 0; i < (signed)OperatorSp.OperatorLine.size(); i++) {
		OperatorJp.OperatorLine.push_back(OperatorSp.OperatorLine[i]);
	}
	for (int i = 0; i < (signed)OperatorLp.OperatorLine.size(); i++) {
		OperatorJp.OperatorLine.push_back(OperatorLp.OperatorLine[i]);
	}

	for (CubicHarmonic_count = 0; CubicHarmonic_count < NOrbitalRoom; CubicHarmonic_count++) {
		read_operator("CubicHarmoic");
	}

}

void create_vac() {

	KETELEMENT obj;
	obj.coeffcientRe = 1.0;
	obj.coeffcientIm = 0.0;
	obj.base = 0;
	vac.KetElement.push_back(obj);

}

void clear_kets() {

	InputKet.KetElement.clear();
	OutputKet.KetElement.clear();
	SimplifiedOutputKet.KetElement.clear();
	SSimplifiedOutputKet.KetElement.clear();
	NormalizedSimplifiedOutputKet.KetElement.clear();
	InputBra.KetElement.clear();

}

void normalized_simplifiedoutputket() {

	long double norm = 0.0;
	for (int i = 0; i < (signed)SSimplifiedOutputKet.KetElement.size(); i++) {
		norm += pow(SSimplifiedOutputKet.KetElement[i].coeffcientRe, 2.0) + pow(SSimplifiedOutputKet.KetElement[i].coeffcientIm, 2.0);
	}
	norm = sqrt(norm);

	for (int i = 0; i < (signed)SSimplifiedOutputKet.KetElement.size(); i++) {
		KETELEMENT obj;
		obj.coeffcientRe = SSimplifiedOutputKet.KetElement[i].coeffcientRe / norm;
		obj.coeffcientIm = SSimplifiedOutputKet.KetElement[i].coeffcientIm / norm;
		obj.base = SSimplifiedOutputKet.KetElement[i].base;
		NormalizedSimplifiedOutputKet.KetElement.push_back(obj);
	}

}

void operate(OPERATOR _operator) {

	for (int i = 0; i < (signed)InputKet.KetElement.size(); i++) {

		Secondization(InputKet.KetElement[i].base);
		/// It might be not so significant to discuss, but for those who are confused...
		/// Definition goes like Array(SecondQuantizedKet/ConvertedQuantizedKet){[0], [1], .... , [NSite-1]} = Exists/Nope at {Site_(NSite-1), Site_(NSite-2), ... , Site_0}
		/// Exists = 1; Nope = 0 

		/// OperatorFiles also should be prepared like...
		/// Coeff. CR{Site_(NSite-1), Site_(NSite-2), ... , Site_0} AN{Site_(NSite-1), Site_(NSite-2), ... , Site_0}

		for (int j = 0; j < (signed)_operator.OperatorLine.size(); j++) {

			bool IfAn = true;
			for (int iSite = NOrbitalRoom - 1; iSite >= 0; iSite--) {
				if (_operator.OperatorLine[j].ANList[iSite] != -1 && _operator.OperatorLine[j].ANList[iSite] != SecondQuantizedKet[iSite]) {
					IfAn = false;
				}
				if (IfAn == false) {
					break;
				}
			}

			if (IfAn == true) {

				int MeaningfulAn = 0;
				for (int iSite = NOrbitalRoom - 1; iSite >= 0; iSite--) {
					if (_operator.OperatorLine[j].ANList[iSite] != -1) {
						MeaningfulAn++;
					}
				}

				vector<int> omitted;
				bool invalid_operate = false;

				for (int iSite = NOrbitalRoom - 1; iSite >= 0; iSite--) {
					if (_operator.OperatorLine[j].ANList[iSite] == -1 && SecondQuantizedKet[iSite] == 1) {
						if (_operator.OperatorLine[j].CRList[iSite] != SecondQuantizedKet[iSite]) {
							omitted.push_back(iSite);
						}
						else {
							invalid_operate = true;
						}
					}
					if (invalid_operate == true) {
						break;
					}
				}

				int count1 = 0;
				int count2 = 0;
				int total_count = 0;

				if (invalid_operate == false) {

					for (int iOmitted = (signed)omitted.size() - 1; iOmitted >= 0; iOmitted--) {
						for (int iSite = omitted[iOmitted] - 1; iSite >= 0; iSite--) {
							if (_operator.OperatorLine[j].ANList[iSite] != -1) {
								count1++;
							}
						}
					}

					for (int iOmitted = 0; iOmitted < (signed)omitted.size(); iOmitted++) {
						for (int iSite = 0; iSite <= omitted[iOmitted] - 1; iSite++) {
							if (_operator.OperatorLine[j].CRList[iSite] != -1) {
								count2++;
							}
						}
					}

					total_count = MeaningfulAn * (MeaningfulAn - 1) / 2 + count1 + count2;

					KETELEMENT obj_element;
					long double buf_re = ComplexMTwoRe(InputKet.KetElement[i].coeffcientRe, InputKet.KetElement[i].coeffcientIm, _operator.OperatorLine[j].coeffcientRe, _operator.OperatorLine[j].coeffcientIm);
					long double buf_im = ComplexMTwoIm(InputKet.KetElement[i].coeffcientRe, InputKet.KetElement[i].coeffcientIm, _operator.OperatorLine[j].coeffcientRe, _operator.OperatorLine[j].coeffcientIm);
					obj_element.coeffcientRe = pow(-1, total_count) * buf_re;
					obj_element.coeffcientIm = pow(-1, total_count) * buf_im;
					for (int iSite = NOrbitalRoom - 1; iSite >= 0; iSite--) {
						ConvertedQuantizedKet[iSite] = 0;
						if (_operator.OperatorLine[j].CRList[iSite] != -1) {
							ConvertedQuantizedKet[iSite] = 1;
						}
						for (int iOmitted = 0; iOmitted < (signed)omitted.size(); iOmitted++) {
							if (omitted[iOmitted] == iSite) {
								ConvertedQuantizedKet[iSite] = 1;
								break;
							}
						}
					}

					obj_element.base = ConvertedKetNumber();
					OutputKet.KetElement.push_back(obj_element);

				}

			}

		}

	}

	bool *ElementChecked;
	ElementChecked = new bool[(signed)OutputKet.KetElement.size()];
	for (int i = 0; i < (signed)OutputKet.KetElement.size(); i++) {
		ElementChecked[i] = false;
	}

	for (int i = 0; i < (signed)OutputKet.KetElement.size(); i++) {
		if (ElementChecked[i] == false) {
			KETELEMENT obj_ketelement;
			obj_ketelement.coeffcientRe = 0.0;
			obj_ketelement.coeffcientIm = 0.0;
			obj_ketelement.base = OutputKet.KetElement[i].base;
			for (int j = i; j < (signed)OutputKet.KetElement.size(); j++) {
				if (ElementChecked[j] == false) {
					if (OutputKet.KetElement[j].base == OutputKet.KetElement[i].base) {
						obj_ketelement.coeffcientRe += OutputKet.KetElement[j].coeffcientRe;
						obj_ketelement.coeffcientIm += OutputKet.KetElement[j].coeffcientIm;
						ElementChecked[j] = true;
					}
				}
			}
			SimplifiedOutputKet.KetElement.push_back(obj_ketelement);
		}
	}

	for (int i = 0; i < (signed)SimplifiedOutputKet.KetElement.size(); i++) {
		if ((ChopBool(SimplifiedOutputKet.KetElement[i].coeffcientRe) && ChopBool(SimplifiedOutputKet.KetElement[i].coeffcientIm)) == false) {
			SSimplifiedOutputKet.KetElement.push_back(SimplifiedOutputKet.KetElement[i]);
		}
	}

}

void _operate_scalar() {

	operate_scalar_re = 0.0;
	operate_scalar_im = 0.0;

	for (int i = 0; i < (signed)InputKet.KetElement.size(); i++) {
		for (int j = 0; j < (signed)InputBra.KetElement.size(); j++) {
			if (InputKet.KetElement[i].base == InputBra.KetElement[j].base) {
				long double re_buf = ComplexMTwoRe(InputBra.KetElement[j].coeffcientRe, -InputBra.KetElement[j].coeffcientIm, InputKet.KetElement[i].coeffcientRe, InputKet.KetElement[i].coeffcientIm);
				long double im_buf = ComplexMTwoIm(InputBra.KetElement[j].coeffcientRe, -InputBra.KetElement[j].coeffcientIm, InputKet.KetElement[i].coeffcientRe, InputKet.KetElement[i].coeffcientIm);
				operate_scalar_re += re_buf;
				operate_scalar_im += im_buf;
			}
		}
	}

}

void operate_scalar() {

	operate_scalar_re = 0.0;
	operate_scalar_im = 0.0;

	for (int i = 0; i < (signed)SSimplifiedOutputKet.KetElement.size(); i++) {
		for (int j = 0; j < (signed)InputBra.KetElement.size(); j++) {
			if (SSimplifiedOutputKet.KetElement[i].base == InputBra.KetElement[j].base) {
				long double re_buf = ComplexMTwoRe(InputBra.KetElement[j].coeffcientRe, -InputBra.KetElement[j].coeffcientIm, SSimplifiedOutputKet.KetElement[i].coeffcientRe, SSimplifiedOutputKet.KetElement[i].coeffcientIm);
				long double im_buf = ComplexMTwoIm(InputBra.KetElement[j].coeffcientRe, -InputBra.KetElement[j].coeffcientIm, SSimplifiedOutputKet.KetElement[i].coeffcientRe, SSimplifiedOutputKet.KetElement[i].coeffcientIm);
				operate_scalar_re += re_buf;
				operate_scalar_im += im_buf;
			}
		}
	}

}

void write_subset() {

	writeoutput << "ket_count\ttestket\tSecondQuantizedKet\tLx.re\tLx.im\tLy.re\tLy.im\tLz.re\tLz.im\tSx.re\tSx.im\tSy.re\tSy.im\tSz.re\tSz.im\n";
	writeoutput.flush();

	signed long long int ket_count = 0;
	for (signed long long int testket = 0; testket < (signed long long int)pow(2, NOrbitalRoom); testket++) {
		Secondization(testket);
		int count_occupied = 0;
		for (int i = 0; i < NOrbitalRoom; i++) {
			if (SecondQuantizedKet[i] == 1) {
				count_occupied++;
			}
		}
		if (count_occupied == NOirbtalOccupied) {

			LISTREPETITIVENC obj;
			obj.ket = testket;

			writeoutput << ket_count << "\t" << testket << "\t";
			writeoutput.flush();
			for (int i = 0; i < NOrbitalRoom; i++) {
				writeoutput << SecondQuantizedKet[i];
				writeoutput.flush();
			}
			writeoutput << "\t";
			writeoutput.flush();

			///////////////////////////////////////////////////
			KETELEMENT obj_ketelement;
			obj_ketelement.coeffcientRe = 1.0;
			obj_ketelement.coeffcientIm = 0.0;
			obj_ketelement.base = testket;

			clear_kets();

			InputKet.KetElement.push_back(obj_ketelement);
			InputBra.KetElement.push_back(obj_ketelement);

			operate(OperatorLx);
			operate_scalar();
			writeoutput << fixed << setprecision(1) << operate_scalar_re << "\t" << operate_scalar_im << "\t";
			writeoutput.flush();
			obj.Lx = operate_scalar_re;

			clear_kets();

			InputKet.KetElement.push_back(obj_ketelement);
			InputBra.KetElement.push_back(obj_ketelement);

			operate(OperatorLy);
			operate_scalar();
			writeoutput << operate_scalar_re << "\t" << operate_scalar_im << "\t";
			writeoutput.flush();
			obj.Ly = operate_scalar_re;

			clear_kets();

			InputKet.KetElement.push_back(obj_ketelement);
			InputBra.KetElement.push_back(obj_ketelement);

			operate(OperatorLz);
			operate_scalar();
			writeoutput << operate_scalar_re << "\t" << operate_scalar_im << "\t";
			writeoutput.flush();
			obj.Lz = operate_scalar_re;

			clear_kets();

			InputKet.KetElement.push_back(obj_ketelement);
			InputBra.KetElement.push_back(obj_ketelement);

			operate(OperatorSx);
			operate_scalar();
			writeoutput << operate_scalar_re << "\t" << operate_scalar_im << "\t";
			writeoutput.flush();
			obj.Sx = operate_scalar_re;

			clear_kets();

			InputKet.KetElement.push_back(obj_ketelement);
			InputBra.KetElement.push_back(obj_ketelement);

			operate(OperatorSy);
			operate_scalar();
			writeoutput << operate_scalar_re << "\t" << operate_scalar_im << "\t";
			writeoutput.flush();
			obj.Sy = operate_scalar_re;

			clear_kets();

			InputKet.KetElement.push_back(obj_ketelement);
			InputBra.KetElement.push_back(obj_ketelement);

			operate(OperatorSz);
			operate_scalar();
			writeoutput << operate_scalar_re << "\t" << operate_scalar_im;
			writeoutput.flush();
			obj.Sz = operate_scalar_re;

			clear_kets();
			///////////////////////////////////////////////////

			writeoutput << "\n";
			writeoutput.flush();
			listRepetitiveNC.push_back(obj);
			ket_count++;

			///////////////////////////////////////////////////
			clear_kets();

			Secondization(testket);
			int *SecondQuantizedKet_buf;
			SecondQuantizedKet_buf = new int[NOrbitalRoom];
			for (int i = 0; i < NOrbitalRoom; i++) {
				SecondQuantizedKet_buf[i] = SecondQuantizedKet[i];
			}

			KET InputKet_buf = vac;
			for (int i = 0; i < NOrbitalRoom; i++) {
				if (SecondQuantizedKet_buf[i] == 1) {

					///////////////////////////////////////////////////
					clear_kets();

					InputKet = InputKet_buf;
					operate(operatorCubicHarmonic[i]);
					InputKet_buf = SSimplifiedOutputKet;

					clear_kets();
					///////////////////////////////////////////////////

				}
			}
			
			LISTREPETITIVENCCUBICHARMONIC obj_ch;
			obj_ch.starket = testket;
			obj_ch.ket = InputKet_buf;
			listRepetitiveNCCubicHarmonic.push_back(obj_ch);

			clear_kets();
			///////////////////////////////////////////////////

		}
	}

	writeoutput << listRepetitiveNC.size() << " bases (all in a ket form) are created; ";
	writeoutput.flush();
	if ((signed long long int)listRepetitiveNC.size() == nCr((signed long long int)NOrbitalRoom, (signed long long int)NOirbtalOccupied)) {
		writeoutput << "as expected.";
		writeoutput.flush();
	}
	else {
		writeoutput << "something wrong! (signed)listRepetitiveNC.size() != nCr(NOrbitalRoom, NOirbtalOccupied).";
		writeoutput.flush();
	}

	writeoutput << "\n";
	writeoutput.flush();
	for (int i = 0; i < (signed)listRepetitiveNCCubicHarmonic.size(); i++) {
		Secondization(listRepetitiveNCCubicHarmonic[i].starket);
		listRepetitiveNCCubicHarmonic[i].starket_rep = "*";
		for (int j = 0; j < NOrbitalRoom; j++) {
			listRepetitiveNCCubicHarmonic[i].starket_rep += to_string(SecondQuantizedKet[j]);
		}
		writeoutput << i << "\t" << listRepetitiveNCCubicHarmonic[i].starket_rep << "\t=\t||\tNbases=\t" << (signed)listRepetitiveNCCubicHarmonic[i].ket.KetElement.size() << "\t||\t";
		writeoutput.flush();
		for (int j = 0; j < (signed)listRepetitiveNCCubicHarmonic[i].ket.KetElement.size(); j++) {
			writeoutput << fixed << setprecision(20) << listRepetitiveNCCubicHarmonic[i].ket.KetElement[j].coeffcientRe << "\t" << listRepetitiveNCCubicHarmonic[i].ket.KetElement[j].coeffcientIm << "\t";
			writeoutput.flush();
			Secondization(listRepetitiveNCCubicHarmonic[i].ket.KetElement[j].base);
			for (int k = 0; k < NOrbitalRoom; k++) {
				writeoutput << SecondQuantizedKet[k];
				writeoutput.flush();
			}
			if (j != (signed)listRepetitiveNCCubicHarmonic[i].ket.KetElement.size() - 1) {
				writeoutput << "\t|\t";
				writeoutput.flush();
			}
			else {
				writeoutput << "\t||\n";
				writeoutput.flush();
			}
		}
	}

	writeoutput << listRepetitiveNCCubicHarmonic.size() << " bases (all in a cubic harmonic ket form) are created; ";
	writeoutput.flush();
	if ((signed long long int)listRepetitiveNCCubicHarmonic.size() == nCr((signed long long int)NOrbitalRoom, (signed long long int)NOirbtalOccupied)) {
		writeoutput << "as expected.";
		writeoutput.flush();
	}
	else {
		writeoutput << "something wrong! (signed)listRepetitiveNCCubicHarmonic.size() != nCr(NOrbitalRoom, NOirbtalOccupied).";
		writeoutput.flush();
	}

}

void eitherLmorSm(int _LmSm_, signed long long int iListRepetitiveNC, int iGramSchmidt, signed long long int iLSKetList) {

	///////////////////////////////////////////////////
	clear_kets();

	if (_LmSm_ == 0) {
		KETELEMENT obj_ketelement;
		obj_ketelement.coeffcientRe = 1.0;
		obj_ketelement.coeffcientIm = 0.0;
		obj_ketelement.base = listRepetitiveNC[iListRepetitiveNC].ket;
		InputKet.KetElement.push_back(obj_ketelement);
	}
	else if (_LmSm_ == 1) {
		InputKet = GramSchmidtKet[iGramSchmidt];
	}
	else if (_LmSm_ == 2 || _LmSm_ == 3) {
		InputKet = LSKetList[iLSKetList].ket;
	}

	if (_LmSm_ == 0 || _LmSm_ == 1) {
		LSKETLIST obj_LSKet;
		obj_LSKet.ket = InputKet;
		LSKetList.push_back(obj_LSKet);
	}
	else if (_LmSm_ == 2) {
		operate(OperatorSm);
		normalized_simplifiedoutputket();
		LSKETLIST obj_LSKet;
		obj_LSKet.ket = NormalizedSimplifiedOutputKet;
		LSKetList.push_back(obj_LSKet);
	}
	else {
		operate(OperatorLm);
		normalized_simplifiedoutputket();
		LSKETLIST obj_LSKet;
		obj_LSKet.ket = NormalizedSimplifiedOutputKet;
		LSKetList.push_back(obj_LSKet);
	}

	if (_LmSm_ == 0 || _LmSm_ == 1 || _LmSm_ == 2) {
		initial_mLindex = (signed long long int)LSKetList.size() - 1;
	}

	clear_kets();
	///////////////////////////////////////////////////

}

void invalid_eitherLmorSm_exit() {

	writeoutput << "\n***************\n***************\nListRepetitiveNCLocked != true\n***************\n***************\n";
	writeoutput.flush();
	exit(EXIT_FAILURE);

}

void invalid_catchedlocked_exit() {

	writeoutput << "\n***************\n***************\ncatchedlocked != true\n***************\n***************\n";
	writeoutput.flush();
	exit(EXIT_FAILURE);

}

void _LmSm(long double S, long double L) {

	signed long long int iListRepetitiveNC;
	bool ListRepetitiveNCLocked = false;
	for (signed long long int i = 0; i < (signed)listRepetitiveNC.size(); i++) {
		if (ChopBool(listRepetitiveNC[i].Sx) && ChopBool(listRepetitiveNC[i].Sy) && ChopBool(listRepetitiveNC[i].Sz - S) && ChopBool(listRepetitiveNC[i].Lx) && ChopBool(listRepetitiveNC[i].Ly) && ChopBool(listRepetitiveNC[i].Lz - L)) {
			iListRepetitiveNC = i;
			ListRepetitiveNCLocked = true;
		}
	}

	if (ListRepetitiveNCLocked == true) {

		for (long double mS = S; mS + S >= -Chop; mS -= 1.0) {
			for (long double mL = L; mL + L >= -Chop; mL -= 1.0) {
				if (ChopBool(mS - S)) {
					if (ChopBool(mL - L)) {
						eitherLmorSm(0, iListRepetitiveNC, 0, 0);
					}
					else {
						eitherLmorSm(3, 0, 0, (signed long long int)LSKetList.size() - 1);
					}
				}
				else {
					if (ChopBool(mL - L)) {
						eitherLmorSm(2, 0, 0, initial_mLindex);
					}
					else {
						eitherLmorSm(3, 0, 0, (signed long long int)LSKetList.size() - 1);
					}
				}
			}
		}

	}
	else {
		invalid_eitherLmorSm_exit();
	}

}

void GramSchmidt(int LS0orLSJ1) {
	
	bool *common_bases;
	vector<signed long long int> _common_bases;
	int NVectors;

	if (LS0orLSJ1 == 0) {

		common_bases = new bool[(signed)listRepetitiveNC.size()];
		for (signed long long int ilistRepetitiveNC = 0; ilistRepetitiveNC < (signed)listRepetitiveNC.size(); ilistRepetitiveNC++) {
			common_bases[ilistRepetitiveNC] = false;
		}

		for (signed long long int icatched = 0; icatched < (signed)catched.size(); icatched++) {
			for (signed long long int icatchedbases = 0; icatchedbases < (signed)LSKetList[catched[icatched]].ket.KetElement.size(); icatchedbases++) {
				for (signed long long int ilistRepetitiveNC = 0; ilistRepetitiveNC < (signed)listRepetitiveNC.size(); ilistRepetitiveNC++) {
					if (LSKetList[catched[icatched]].ket.KetElement[icatchedbases].base == listRepetitiveNC[ilistRepetitiveNC].ket) {
						common_bases[ilistRepetitiveNC] = true;
					}
				}
			}
		}

		for (signed long long int ilistRepetitiveNC = 0; ilistRepetitiveNC < (signed)listRepetitiveNC.size(); ilistRepetitiveNC++) {
			if (common_bases[ilistRepetitiveNC] == true) {
				_common_bases.push_back(listRepetitiveNC[ilistRepetitiveNC].ket);
			}
		}

		NVectors = (signed)catched.size() + (signed)same_L_series.size();
	
	}
	else {

		common_bases = new bool[global_iend - global_istart + 1];
		for (signed long long int global_i = 0; global_i < global_iend - global_istart + 1; global_i++) {
			common_bases[global_i] = false;
		}

		for (signed long long int icatched = 0; icatched < (signed)catched.size(); icatched++) {
			for (signed long long int global_i = global_istart; global_i <= global_iend; global_i++) {
				///////////////////////////////////////////////////
				clear_kets();
				InputKet = LSJKetList[catched[icatched]].ket;
				InputBra = LSKetList[global_i].ket;
				_operate_scalar();
				if (!ChopBool(operate_scalar_re) || !ChopBool(operate_scalar_im)) {
					common_bases[global_i - global_istart] = true;
				}
				clear_kets();
				///////////////////////////////////////////////////
			}
		}

		for (signed long long int global_i = global_istart; global_i <= global_iend; global_i++) {
			if (common_bases[global_i - global_istart] == true) {
				_common_bases.push_back(global_i);
			}
		}

		NVectors = (signed)catched.size() + 1;
	
	}
		
	long double **rearrangedRe, **rearrangedIm;
	rearrangedRe = new long double *[NVectors];
	rearrangedIm = new long double *[NVectors];
	for (signed long long int icatched = 0; icatched < NVectors; icatched++) {
		rearrangedRe[icatched] = new long double[(signed)_common_bases.size()];
		rearrangedIm[icatched] = new long double[(signed)_common_bases.size()];
		for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
			rearrangedRe[icatched][i_common_bases] = 0.0;
			rearrangedIm[icatched][i_common_bases] = 0.0;
		}
	}

	for (signed long long int icatched = 0; icatched < NVectors; icatched++) {
		if (icatched < (signed)catched.size()) {
			if (LS0orLSJ1 == 0) {
				for (signed long long int icatchedbases = 0; icatchedbases < (signed)LSKetList[catched[icatched]].ket.KetElement.size(); icatchedbases++) {
					for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
						if (LSKetList[catched[icatched]].ket.KetElement[icatchedbases].base == _common_bases[i_common_bases]) {
							rearrangedRe[icatched][i_common_bases] = LSKetList[catched[icatched]].ket.KetElement[icatchedbases].coeffcientRe;
							rearrangedIm[icatched][i_common_bases] = LSKetList[catched[icatched]].ket.KetElement[icatchedbases].coeffcientIm;
						}
					}
				}
			}
			else {
				for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
					///////////////////////////////////////////////////
					clear_kets();
					InputKet = LSJKetList[catched[icatched]].ket;
					InputBra = LSKetList[_common_bases[i_common_bases]].ket;
					_operate_scalar();
					rearrangedRe[icatched][i_common_bases] = operate_scalar_re;
					rearrangedIm[icatched][i_common_bases] = operate_scalar_im;
					clear_kets();
					///////////////////////////////////////////////////
				}
			}
		}
		else {
			for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
				rearrangedRe[icatched][i_common_bases] = icatched + i_common_bases + 1.0;
				rearrangedIm[icatched][i_common_bases] = 0.0;
			}
			long double renorm = 0.0;
			for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
				renorm += pow(rearrangedRe[icatched][i_common_bases], 2.0) + pow(rearrangedIm[icatched][i_common_bases], 2.0);
			}
			renorm = sqrt(renorm);
			for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
				rearrangedRe[icatched][i_common_bases] /= renorm;
				rearrangedIm[icatched][i_common_bases] /= renorm;
			}
		}	
	}

	int GS_count = 0;
	while (true) {

		for (signed long long int icatched = (signed)catched.size(); icatched < NVectors; icatched++) {

			long double *orthonormalizedRe;
			long double *orthonormalizedIm;
			orthonormalizedRe = new long double[(signed)_common_bases.size()];
			orthonormalizedIm = new long double[(signed)_common_bases.size()];
			for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
				orthonormalizedRe[i_common_bases] = rearrangedRe[icatched][i_common_bases];
				orthonormalizedIm[i_common_bases] = rearrangedIm[icatched][i_common_bases];
			}

			for (signed long long int iicatched = 0; iicatched < icatched; iicatched++) {
				long double vdotuRe = 0.0;
				long double vdotuIm = 0.0;
				for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
					vdotuRe += ComplexMTwoRe(orthonormalizedRe[i_common_bases], -orthonormalizedIm[i_common_bases], rearrangedRe[iicatched][i_common_bases], rearrangedIm[iicatched][i_common_bases]);
					vdotuIm += ComplexMTwoIm(orthonormalizedRe[i_common_bases], -orthonormalizedIm[i_common_bases], rearrangedRe[iicatched][i_common_bases], rearrangedIm[iicatched][i_common_bases]);
				}
				for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
					orthonormalizedRe[i_common_bases] -= ComplexMTwoRe(vdotuRe, vdotuIm, rearrangedRe[iicatched][i_common_bases], rearrangedIm[iicatched][i_common_bases]);
					orthonormalizedIm[i_common_bases] -= ComplexMTwoIm(vdotuRe, vdotuIm, rearrangedRe[iicatched][i_common_bases], rearrangedIm[iicatched][i_common_bases]);
				}
			}

			long double normnorm = 0.0;
			for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
				normnorm += pow(orthonormalizedRe[i_common_bases], 2.0) + pow(orthonormalizedIm[i_common_bases], 2.0);
			}
			normnorm = sqrt(normnorm);
			for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
				orthonormalizedRe[i_common_bases] /= normnorm;
				orthonormalizedIm[i_common_bases] /= normnorm;
			}
			for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
				rearrangedRe[icatched][i_common_bases] = orthonormalizedRe[i_common_bases];
				rearrangedIm[icatched][i_common_bases] = orthonormalizedIm[i_common_bases];
			}

		}

		// ORTHO?
		long double **ortho_check_re, **ortho_check_im;
		ortho_check_re = new long double *[NVectors];
		ortho_check_im = new long double *[NVectors];
		for (int i = 0; i < NVectors; i++) {
			ortho_check_re[i] = new long double[NVectors];
			ortho_check_im[i] = new long double[NVectors];
		}

		for (int i = 0; i < NVectors; i++) {
			for (int j = i; j < NVectors; j++) {
				ortho_check_re[i][j] = 0.0;
				ortho_check_im[i][j] = 0.0;
				for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
					ortho_check_re[i][j] += ComplexMTwoRe(rearrangedRe[i][i_common_bases], -rearrangedIm[i][i_common_bases], rearrangedRe[j][i_common_bases], rearrangedIm[j][i_common_bases]);
					ortho_check_im[i][j] += ComplexMTwoIm(rearrangedRe[i][i_common_bases], -rearrangedIm[i][i_common_bases], rearrangedRe[j][i_common_bases], rearrangedIm[j][i_common_bases]);
				}
			}
		}

		for (int i = 0; i < NVectors; i++) {
			for (int j = 0; j < i; j++) {
				ortho_check_re[i][j] = ortho_check_re[j][i];
				ortho_check_im[i][j] = -ortho_check_im[j][i];
			}
		}

		bool ortho_check = true;
		for (int i = 0; i < NVectors; i++) {
			for (int j = 0; j < NVectors; j++) {
				if (i == j) {
					if (!ChopBool(ortho_check_re[i][j] - 1.0) || !ChopBool(ortho_check_im[i][j])) {
						ortho_check = false;
						break;
					}
				}
				else {
					if (!ChopBool(ortho_check_re[i][j]) || !ChopBool(ortho_check_im[i][j])) {
						ortho_check = false;
						break;
					}
				}
			}
		}

		GS_count++;
		
		if (ortho_check == true) {
			break;
		}
		else {
			if (GS_count == 100) {
				writeoutput << "\n-100\t=============ortho_check_re_warning=============\n";
				writeoutput.flush();
				for (int i = 0; i < NVectors; i++) {
					for (int j = 0; j < NVectors; j++) {
						writeoutput << ortho_check_re[i][j];
						writeoutput.flush();
						if (j != NVectors - 1) {
							writeoutput << "\t";
							writeoutput.flush();
						}
						else {
							writeoutput << "\n";
							writeoutput.flush();
						}
					}
				}
				writeoutput << "-100\t=============ortho_check_re_warning=============";
				writeoutput.flush();
				break;
			}
		}

	}

	for (signed long long int icatched = (signed)catched.size(); icatched < NVectors; icatched++) {
		
		if (LS0orLSJ1 == 0) {

			KET obj_ket;
			for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
				KETELEMENT obj_ket_ele;
				obj_ket_ele.coeffcientRe = rearrangedRe[icatched][i_common_bases];
				obj_ket_ele.coeffcientIm = rearrangedIm[icatched][i_common_bases];
				obj_ket_ele.base = _common_bases[i_common_bases];
				obj_ket.KetElement.push_back(obj_ket_ele);
			}
			GramSchmidtKet.push_back(obj_ket);

		}
		else {

			KET obj_ket;
			KET obj_ket_re;

			for (signed long long int i_common_bases = 0; i_common_bases < (signed)_common_bases.size(); i_common_bases++) {
				for (int iket = 0; iket < (signed)LSKetList[_common_bases[i_common_bases]].ket.KetElement.size(); iket++) {
					KETELEMENT obj_ket_ele;
					obj_ket_ele.coeffcientRe = rearrangedRe[icatched][i_common_bases] * LSKetList[_common_bases[i_common_bases]].ket.KetElement[iket].coeffcientRe;
					obj_ket_ele.coeffcientIm = rearrangedIm[icatched][i_common_bases] * LSKetList[_common_bases[i_common_bases]].ket.KetElement[iket].coeffcientIm;
					obj_ket_ele.base = LSKetList[_common_bases[i_common_bases]].ket.KetElement[iket].base;
					obj_ket.KetElement.push_back(obj_ket_ele);
				}
			}

			bool *ElementChecked;
			ElementChecked = new bool[(signed)obj_ket.KetElement.size()];
			for (int i = 0; i < (signed)obj_ket.KetElement.size(); i++) {
				ElementChecked[i] = false;
			}

			for (int i = 0; i < (signed)obj_ket.KetElement.size(); i++) {
				if (ElementChecked[i] == false) {
					KETELEMENT obj_ketelement;
					obj_ketelement.coeffcientRe = 0.0;
					obj_ketelement.coeffcientIm = 0.0;
					obj_ketelement.base = obj_ket.KetElement[i].base;
					for (int j = i; j < (signed)obj_ket.KetElement.size(); j++) {
						if (ElementChecked[j] == false) {
							if (obj_ket.KetElement[j].base == obj_ket.KetElement[i].base) {
								obj_ketelement.coeffcientRe += obj_ket.KetElement[j].coeffcientRe;
								obj_ketelement.coeffcientIm += obj_ket.KetElement[j].coeffcientIm;
								ElementChecked[j] = true;
							}
						}
					}
					obj_ket_re.KetElement.push_back(obj_ketelement);
				}
			}

			long double _norm = 0.0;
			for (int i = 0; i < (signed)obj_ket_re.KetElement.size(); i++) {
				_norm += pow(obj_ket_re.KetElement[i].coeffcientRe, 2.0) + pow(obj_ket_re.KetElement[i].coeffcientIm, 2.0);
			}
			_norm = sqrt(_norm);
			for (int i = 0; i < (signed)obj_ket_re.KetElement.size(); i++) {
				obj_ket_re.KetElement[i].coeffcientRe /= _norm;
				obj_ket_re.KetElement[i].coeffcientIm /= _norm;
			}

			GramSchmidtKet.push_back(obj_ket_re);
			
		}
				
	}

	if (LS0orLSJ1 == 0) {

		MatrixXcd DiagonalizeGramSchimdtHamiltonian;
		DiagonalizeGramSchimdtHamiltonian = MatrixXcd::Random((signed)GramSchmidtKet.size(), (signed)GramSchmidtKet.size());
		for (int i = 0; i < (signed)GramSchmidtKet.size(); i++) {
			for (int j = 0; j < (signed)GramSchmidtKet.size(); j++) {

				long double u_re = 0.0;
				long double u_im = 0.0;
				///////////////////////////////////////////////////
				for (int iCoulomb = 0; iCoulomb < NCoulombHamiltonian; iCoulomb++) {
					clear_kets();
					InputKet = GramSchmidtKet[j];
					InputBra = GramSchmidtKet[i];
					operate(Coulomb[iCoulomb].CoulombHamiltonian);
					operate_scalar();
					clear_kets();
					u_re += operate_scalar_re * Coulomb[iCoulomb].CoulombCoeff;
					u_im += operate_scalar_im * Coulomb[iCoulomb].CoulombCoeff;
				}
				///////////////////////////////////////////////////
				DiagonalizeGramSchimdtHamiltonian(i, j) = complex<double>(u_re, u_im);

			}
		}

		ComplexEigenSolver<MatrixXcd> HamKsolver(DiagonalizeGramSchimdtHamiltonian);

		long double ***eigenvector;
		eigenvector = new long double **[(signed)GramSchmidtKet.size()];
		for (int i = 0; i < (signed)GramSchmidtKet.size(); i++) {
			eigenvector[i] = new long double *[(signed)GramSchmidtKet.size()];
			for (int j = 0; j < (signed)GramSchmidtKet.size(); j++) {
				eigenvector[i][j] = new long double[2];
			}
		}

		int element_count = 0;
		for (int i = 0; i < (signed)GramSchmidtKet.size(); i++) {
			for (int j = 0; j < (signed)GramSchmidtKet.size(); j++) {
				if (!ChopBool(HamKsolver.eigenvectors().data()[element_count].real())) {
					eigenvector[i][j][0] = HamKsolver.eigenvectors().data()[element_count].real();
				}
				else {
					eigenvector[i][j][0] = 0.0;
				}
				if (!ChopBool(HamKsolver.eigenvectors().data()[element_count].imag())) {
					eigenvector[i][j][1] = HamKsolver.eigenvectors().data()[element_count].imag();
				}
				else {
					eigenvector[i][j][1] = 0.0;
				}
				element_count++;
			}
		}

		vector<KET> GramSchmidt_buf;
		for (int i = 0; i < (signed)GramSchmidtKet.size(); i++) {
			KET obj_GSbuf;
			for (int k = 0; k < (signed)_common_bases.size(); k++) {
				KETELEMENT obj_GSbufEle;
				obj_GSbuf.KetElement.push_back(obj_GSbufEle);
			}
			for (int k = 0; k < (signed)_common_bases.size(); k++) {
				obj_GSbuf.KetElement[k].coeffcientRe = 0.0;
				obj_GSbuf.KetElement[k].coeffcientIm = 0.0;
				obj_GSbuf.KetElement[k].base = _common_bases[k];
			}
			for (int j = 0; j < (signed)GramSchmidtKet.size(); j++) {
				for (int k = 0; k < (signed)_common_bases.size(); k++) {
					obj_GSbuf.KetElement[k].coeffcientRe += ComplexMTwoRe(eigenvector[i][j][0], eigenvector[i][j][1], GramSchmidtKet[j].KetElement[k].coeffcientRe, GramSchmidtKet[j].KetElement[k].coeffcientIm);
					obj_GSbuf.KetElement[k].coeffcientIm += ComplexMTwoIm(eigenvector[i][j][0], eigenvector[i][j][1], GramSchmidtKet[j].KetElement[k].coeffcientRe, GramSchmidtKet[j].KetElement[k].coeffcientIm);
				}
			}
			GramSchmidt_buf.push_back(obj_GSbuf);
		}

		// Final check for U
		for (int i = 0; i < (signed)GramSchmidt_buf.size(); i++) {
			for (int j = i + 1; j < (signed)GramSchmidt_buf.size(); j++) {

				long double u_re = 0.0;
				long double u_im = 0.0;
				///////////////////////////////////////////////////
				for (int iCoulomb = 0; iCoulomb < NCoulombHamiltonian; iCoulomb++) {
					clear_kets();
					InputKet = GramSchmidt_buf[j];
					InputBra = GramSchmidt_buf[i];
					operate(Coulomb[iCoulomb].CoulombHamiltonian);
					operate_scalar();
					clear_kets();
					u_re += operate_scalar_re * Coulomb[iCoulomb].CoulombCoeff;
					u_im += operate_scalar_im * Coulomb[iCoulomb].CoulombCoeff;
				}
				///////////////////////////////////////////////////
				if (!ChopBool(u_re) || !ChopBool(u_im)) {
					writeoutput << "\n******************NOT FULLY DIAGONALIZED!!!!******************\n";
					writeoutput.flush();
				}

			}
		}

		GramSchmidtKet.clear();
		for (int i = 0; i < (signed)GramSchmidt_buf.size(); i++) {
			KET obj_GS;
			long double norm_calc = 0.0;
			for (int k = 0; k < (signed)_common_bases.size(); k++) {
				norm_calc += pow(GramSchmidt_buf[i].KetElement[k].coeffcientRe, 2.0) + pow(GramSchmidt_buf[i].KetElement[k].coeffcientIm, 2.0);
			}
			norm_calc = sqrt(norm_calc);
			for (int k = 0; k < (signed)_common_bases.size(); k++) {
				if (!ChopBool(GramSchmidt_buf[i].KetElement[k].coeffcientRe) || !ChopBool(GramSchmidt_buf[i].KetElement[k].coeffcientIm)) {
					KETELEMENT obj_GSKE;
					obj_GSKE.base = _common_bases[k];
					obj_GSKE.coeffcientRe = GramSchmidt_buf[i].KetElement[k].coeffcientRe / norm_calc;
					obj_GSKE.coeffcientIm = GramSchmidt_buf[i].KetElement[k].coeffcientIm / norm_calc;
					obj_GS.KetElement.push_back(obj_GSKE);
				}
			}
			GramSchmidtKet.push_back(obj_GS);
		}

	}

}

void LmSm(long double S, long double L) {

	if (iGramSchmidt == 0) {

		GramSchmidtKet.clear();
		catched.clear();
		bool catchedlocked = false;
		for (signed long long int i = 0; i < (signed)LSKetList.size(); i++) {
			if (ChopBool(LSKetList[i].Sz - S) && ChopBool(LSKetList[i].Lz - L)) {
				catched.push_back(i);
				catchedlocked = true;
			}
		}

		if (catchedlocked == true) {
			GramSchmidt(0);
		}
		else {
			invalid_catchedlocked_exit();
		}

	}
	
	for (long double mS = S; mS + S >= -Chop; mS -= 1.0) {
		for (long double mL = L; mL + L >= -Chop; mL -= 1.0) {
			if (ChopBool(mS - S)) {
				if (ChopBool(mL - L)) {
					eitherLmorSm(1, 0, iGramSchmidt, 0);
					iGramSchmidt++;
				}
				else {
					eitherLmorSm(3, 0, 0, (signed long long int)LSKetList.size() - 1);
				}
			}
			else {
				if (ChopBool(mL - L)) {
					eitherLmorSm(2, 0, 0, initial_mLindex);
				}
				else {
					eitherLmorSm(3, 0, 0, (signed long long int)LSKetList.size() - 1);
				}
			}
		}
	}

}

void write_LSKetList(signed long long int LSKetList_start, signed long long int LSKetList_end, long double S, long double L) {

	if (EigenvectorOutputMode == 0) {
		writeoutput << "\n";
		writeoutput.flush();
	}
	
	for (signed long long int i = LSKetList_start; i <= LSKetList_end; i++) {

		LSKetList[i].S = S;
		LSKetList[i].L = L;

		///////////////////////////////////////////////////
		clear_kets();

		InputKet = LSKetList[i].ket;
		InputBra = LSKetList[i].ket;
		operate(OperatorSx);
		operate_scalar();
		LSKetList[i].Sx = operate_scalar_re;

		clear_kets();

		InputKet = LSKetList[i].ket;
		InputBra = LSKetList[i].ket;
		operate(OperatorSy);
		operate_scalar();
		LSKetList[i].Sy = operate_scalar_re;

		clear_kets();

		InputKet = LSKetList[i].ket;
		InputBra = LSKetList[i].ket;
		operate(OperatorSz);
		operate_scalar();
		LSKetList[i].Sz = operate_scalar_re;

		clear_kets();

		InputKet = LSKetList[i].ket;
		InputBra = LSKetList[i].ket;
		operate(OperatorLx);
		operate_scalar();
		LSKetList[i].Lx = operate_scalar_re;

		clear_kets();

		InputKet = LSKetList[i].ket;
		InputBra = LSKetList[i].ket;
		operate(OperatorLy);
		operate_scalar();
		LSKetList[i].Ly = operate_scalar_re;

		clear_kets();

		InputKet = LSKetList[i].ket;
		InputBra = LSKetList[i].ket;
		operate(OperatorLz);
		operate_scalar();
		LSKetList[i].Lz = operate_scalar_re;

		clear_kets();

		LSKetList[i].CoulombEnerg = new long double[NCoulombHamiltonian];

		for (int iCoulomb = 0; iCoulomb < NCoulombHamiltonian; iCoulomb++) {

			clear_kets();

			InputKet = LSKetList[i].ket;
			InputBra = LSKetList[i].ket;
			operate(Coulomb[iCoulomb].CoulombHamiltonian);
			operate_scalar();
			LSKetList[i].CoulombEnerg[iCoulomb] = operate_scalar_re;

			clear_kets();

		}
		///////////////////////////////////////////////////
		if (EigenvectorOutputMode == 0) {
			writeoutput << i << "\t||\tL=\t" << fixed << setprecision(1) << LSKetList[i].L << "\t||\tLz=\t" << LSKetList[i].Lz << "\t||\tS=\t" << LSKetList[i].S << "\t||\tSz=\t" << LSKetList[i].Sz << "\t||\tU=\t";
			writeoutput.flush();
		}
		else {
			writeoutput << i << "\t||\tL=\t" << fixed << setprecision(1) << LSKetList[i].L << "\t||\tLz=\t" << LSKetList[i].Lz << "\t||\tS=\t" << LSKetList[i].S << "\t||\tSz=\t" << LSKetList[i].Sz << "\n";
			writeoutput.flush();
			fwrite(&i, sizeof(signed long long int), 1, _writeoutput);
			fwrite(&LSKetList[i].L, sizeof(long double), 1, _writeoutput);
			fwrite(&LSKetList[i].Lz, sizeof(long double), 1, _writeoutput);
			fwrite(&LSKetList[i].S, sizeof(long double), 1, _writeoutput);
			fwrite(&LSKetList[i].Sz, sizeof(long double), 1, _writeoutput);
		}
		
		long double enrg = 0.0;
		for (int iCoulomb = 0; iCoulomb < NCoulombHamiltonian; iCoulomb++) {
			if (EigenvectorOutputMode == 0) {
				writeoutput << fixed << setprecision(20) << LSKetList[i].CoulombEnerg[iCoulomb] << "\t" << Coulomb[iCoulomb].CoulombName << "\t(=\t" << Coulomb[iCoulomb].CoulombCoeff << "\t)";
				writeoutput.flush();
			}
			enrg += LSKetList[i].CoulombEnerg[iCoulomb] * Coulomb[iCoulomb].CoulombCoeff;
			if (EigenvectorOutputMode == 0) {
				if (iCoulomb != NCoulombHamiltonian - 1) {
					writeoutput << "\t+\t";
					writeoutput.flush();
				}
				else {
					writeoutput << "\t=\t";
					writeoutput.flush();
				}
			}
		}
		LSKetList[i].CoulombEnergySum = enrg;
		if (EigenvectorOutputMode == 0) {
			writeoutput << fixed << setprecision(20) << LSKetList[i].CoulombEnergySum << "\t||\tNbases=\t" << (signed)LSKetList[i].ket.KetElement.size() << "\t||\t";
			writeoutput.flush();
		}
		else {
			signed long long int size_buf = LSKetList[i].ket.KetElement.size();
			fwrite(&LSKetList[i].CoulombEnergySum, sizeof(long double), 1, _writeoutput);
			fwrite(&size_buf, sizeof(signed long long int), 1, _writeoutput);
		}

		for (int j = 0; j < (signed)LSKetList[i].ket.KetElement.size(); j++) {

			if (EigenvectorOutputMode == 0) {
				writeoutput << fixed << setprecision(20) << LSKetList[i].ket.KetElement[j].coeffcientRe << "\t" << LSKetList[i].ket.KetElement[j].coeffcientIm << "\t";
				writeoutput.flush();
			}
			else {
				fwrite(&LSKetList[i].ket.KetElement[j].coeffcientRe, sizeof(long double), 1, _writeoutput);
				fwrite(&LSKetList[i].ket.KetElement[j].coeffcientIm, sizeof(long double), 1, _writeoutput);
			}

			Secondization(LSKetList[i].ket.KetElement[j].base);
			if (EigenvectorOutputMode == 0) {
				for (int k = 0; k < NOrbitalRoom; k++) {
					writeoutput << SecondQuantizedKet[k];
					writeoutput.flush();
				}
				if (j != (signed)LSKetList[i].ket.KetElement.size() - 1) {
					writeoutput << "\t|\t";
					writeoutput.flush();
				}
				else {
					if (i != LSKetList_end) {
						writeoutput << "\t||\n";
						writeoutput.flush();
					}
				}
			}
			else {
				fwrite(&SecondQuantizedKet[0], sizeof(int), NOrbitalRoom, _writeoutput);
			}

		}

		if (EigenvectorOutputMode == 0 && (signed)LSKetList[i].ket.KetElement.size() == 0 && i != LSKetList_end) {
			writeoutput << "\t||\n";
			writeoutput.flush();
		}

	}

}

void roll_the_dice() {

	writeoutput << "\n=====================================roll_the_dice()=====================================";
	writeoutput.flush();

	if (EigenvectorOutputMode == 1) {
		writeoutput << "\n";
		writeoutput.flush();
	}

	bool GramSchmidt_initialize = true;

	LSKET_size_buf = 0;
	for (int i = 0; i < (signed)_SLlist.SLListLine.size(); i++) {
		for (int j = 0; j < (signed)_SLlist.SLListLine[i].L.size(); j++) {
			if (j == 0) {
				_LmSm(_SLlist.SLListLine[i].S, _SLlist.SLListLine[i].L[j]);
			}
			else {

				// provided that there is the same _SLlist.SLListLine[i].L[j] (j!=0) with _SLlist.SLListLine[i].L[0]

				if (GramSchmidt_initialize == true) {
					iGramSchmidt = 0;
					GramSchmidtKet.clear();
					same_L_series.clear();
					for (int k = j; k < (signed)_SLlist.SLListLine[i].L.size(); k++) {
						if (ChopBool(_SLlist.SLListLine[i].L[k] - _SLlist.SLListLine[i].L[j])) {
							same_L_series.push_back(k);
						}
					}
					GramSchmidt_initialize = false;
				}

				LmSm(_SLlist.SLListLine[i].S, _SLlist.SLListLine[i].L[j]);

				if (iGramSchmidt == (signed)same_L_series.size()) {
					GramSchmidt_initialize = true;
				}

			}
			write_LSKetList(LSKET_size_buf, (signed long long int)LSKetList.size() - 1, _SLlist.SLListLine[i].S, _SLlist.SLListLine[i].L[j]);
			LSKET_size_buf = (signed long long int)LSKetList.size();
		}
	}

	if (EigenvectorOutputMode == 0) {
		writeoutput << "\n";
		writeoutput.flush();
	}
	writeoutput << "***************************************roll_the_dice()***************************************\n";
	writeoutput.flush();

	for (signed long long int i = 0; i < (signed long long int)LSKetList.size(); i++) {
		
		if (EigenvectorOutputMode == 0) {
			writeoutput << i << "\t||\tL=\t" << fixed << setprecision(1) << LSKetList[i].L << "\t||\tLz=\t" << LSKetList[i].Lz << "\t||\tS=\t" << LSKetList[i].S << "\t||\tSz=\t" << LSKetList[i].Sz << "\t||\tU=\t";
			writeoutput.flush();
		}
		else {
			writeoutput << i << "\t||\tL=\t" << fixed << setprecision(1) << LSKetList[i].L << "\t||\tLz=\t" << LSKetList[i].Lz << "\t||\tS=\t" << LSKetList[i].S << "\t||\tSz=\t" << LSKetList[i].Sz << "\n";
			writeoutput.flush();
			fwrite(&i, sizeof(signed long long int), 1, _writeoutput);
			fwrite(&LSKetList[i].L, sizeof(long double), 1, _writeoutput);
			fwrite(&LSKetList[i].Lz, sizeof(long double), 1, _writeoutput);
			fwrite(&LSKetList[i].S, sizeof(long double), 1, _writeoutput);
			fwrite(&LSKetList[i].Sz, sizeof(long double), 1, _writeoutput);
		}

		if (EigenvectorOutputMode == 0) {
			for (int iCoulomb = 0; iCoulomb < NCoulombHamiltonian; iCoulomb++) {
				writeoutput << fixed << setprecision(20) << LSKetList[i].CoulombEnerg[iCoulomb] << "\t" << Coulomb[iCoulomb].CoulombName << "\t(=\t" << Coulomb[iCoulomb].CoulombCoeff << "\t)";
				writeoutput.flush();
				if (iCoulomb != NCoulombHamiltonian - 1) {
					writeoutput << "\t+\t";
					writeoutput.flush();
				}
				else {
					writeoutput << "\t=\t";
					writeoutput.flush();
				}
			}
			writeoutput << fixed << setprecision(20) << LSKetList[i].CoulombEnergySum << "\t||\t";
			writeoutput.flush();
		}
		else {
			fwrite(&LSKetList[i].CoulombEnergySum, sizeof(long double), 1, _writeoutput);
		}
		
			
		vector<long double> CubicHarmonicRe;
		vector<long double> CubicHarmonicIm;
		vector<signed long long int> CubicHarmonicBases;
		
		for (signed long long int j = 0; j < (signed long long int)listRepetitiveNCCubicHarmonic.size(); j++) {

			///////////////////////////////////////////////////

			clear_kets();

			InputKet = LSKetList[i].ket;
			InputBra = listRepetitiveNCCubicHarmonic[j].ket;
			_operate_scalar();

			if (!ChopBool(operate_scalar_re) || !ChopBool(operate_scalar_im)) {
				CubicHarmonicRe.push_back(operate_scalar_re);
				CubicHarmonicIm.push_back(operate_scalar_im);
				CubicHarmonicBases.push_back(j);
			}
			
			clear_kets();

			///////////////////////////////////////////////////

		}
				
		if (EigenvectorOutputMode == 0) {
			writeoutput << "*Nbases=\t" << (signed)CubicHarmonicBases.size() << "\t||\t";
			writeoutput.flush();
		}
		else {
			signed long long int size_buf = (signed)CubicHarmonicBases.size();
			fwrite(&size_buf, sizeof(signed long long int), 1, _writeoutput);
		}

		for (signed long long int j = 0; j < (signed)CubicHarmonicBases.size(); j++) {

			if (EigenvectorOutputMode == 0) {
				writeoutput << fixed << setprecision(20) << CubicHarmonicRe[j] << "\t" << CubicHarmonicIm[j] << "\t";
				writeoutput.flush();
				writeoutput << listRepetitiveNCCubicHarmonic[CubicHarmonicBases[j]].starket_rep;
				writeoutput.flush();
				if (j != (signed)CubicHarmonicBases.size() - 1) {
					writeoutput << "\t|\t";
					writeoutput.flush();
				}
				else {
					writeoutput << "\t||\n";
					writeoutput.flush();
				}
			}
			else {
				fwrite(&CubicHarmonicRe[j], sizeof(long double), 1, _writeoutput);
				fwrite(&CubicHarmonicIm[j], sizeof(long double), 1, _writeoutput);
				fwrite(&listRepetitiveNCCubicHarmonic[CubicHarmonicBases[j]].starket, sizeof(signed long long int), 1, _writeoutput);
			}

		}
		
	}
	
	writeoutput << "=====================================roll_the_dice()=====================================\n";
	writeoutput.flush();

	writeoutput << LSKetList.size() << " LS-coupling bases (all in a ket form) are created; ";
	writeoutput.flush();
	if ((signed long long int)LSKetList.size() == nCr((signed long long int)NOrbitalRoom, (signed long long int)NOirbtalOccupied)) {
		writeoutput << "as expected.";
		writeoutput.flush();
	}
	else {
		writeoutput << "something wrong! (signed)LSKetList.size() != nCr(NOrbitalRoom, NOirbtalOccupied).";
		writeoutput.flush();
	}

	LSBLOCK obj_block;
	obj_block.istart = 0;
	for (signed long long int i = 0; i < (signed long long int)LSKetList.size() - 1; i++) {
		if (!ChopBool(LSKetList[i].L - LSKetList[i + 1].L) || !ChopBool(LSKetList[i].S - LSKetList[i + 1].S)) {
			obj_block.iend = i;
			obj_block.multiplicity = (int)((int)(obj_block.iend - obj_block.istart + 1)/((int)(2 * LSKetList[i].L + 1) * (int)(2 * LSKetList[i].S + 1)));
			obj_block.L = LSKetList[i].L;
			obj_block.S = LSKetList[i].S;
			LSBlock.push_back(obj_block);
			obj_block.istart = i + 1;
		}
	}
	obj_block.iend = (signed long long int)LSKetList.size() - 1;
	obj_block.multiplicity = (int)((int)(obj_block.iend - obj_block.istart + 1) / ((int)(2 * LSKetList[(signed long long int)LSKetList.size() - 1].L + 1) * (int)(2 * LSKetList[(signed long long int)LSKetList.size() - 1].S + 1)));
	obj_block.L = LSKetList[(signed long long int)LSKetList.size() - 1].L;
	obj_block.S = LSKetList[(signed long long int)LSKetList.size() - 1].S;
	LSBlock.push_back(obj_block);
	writeoutput << "\nLSBlock\n";
	writeoutput.flush();
	for (int i = 0; i < (signed)LSBlock.size(); i++) {
		writeoutput << i << "\tL=\t" << fixed << setprecision(1) << LSBlock[i].L << "\tS=\t" << LSBlock[i].S << "\tistart=\t" << LSBlock[i].istart << "\tiend=\t" << LSBlock[i].iend << "\tmultiplicity=\t" << LSBlock[i].multiplicity << "\t";
		writeoutput.flush();
		if (LSBlock[i].multiplicity * (int)(2.0 * LSBlock[i].L + 1) * (int)(2.0 * LSBlock[i].S + 1) == LSBlock[i].iend - LSBlock[i].istart + 1) {
			writeoutput << "as expected";
			writeoutput.flush();
		}
		else {
			writeoutput << "wrong count!";
			writeoutput.flush();
		}
		if (i != (signed)LSBlock.size() - 1) {
			writeoutput << "\n";
			writeoutput.flush();
		}
	}

	/*
	writeoutput << "\nchecking <bra_i|H_int|ket_j>_(i!=j)...";
	writeoutput.flush();
	bool braHintKet_stable = true;
	for (int iLSBlock = 0; iLSBlock < (signed)LSBlock.size(); iLSBlock++) {
		if (LSBlock[iLSBlock].multiplicity != 1) {
			for (signed long long int i = LSBlock[iLSBlock].istart; i <= LSBlock[iLSBlock].iend; i++) {
				for (signed long long int j = i + 1; j <= LSBlock[iLSBlock].iend; j++) {

					double ur = 0.0;
					double ui = 0.0;
					for (int iCoulomb = 0; iCoulomb < NCoulombHamiltonian; iCoulomb++) {
						///////////////////////////////////////////////////
						clear_kets();

						InputKet = LSKetList[j].ket;
						InputBra = LSKetList[i].ket;
						operate(Coulomb[iCoulomb].CoulombHamiltonian);
						operate_scalar();
						ur += Coulomb[iCoulomb].CoulombCoeff * operate_scalar_re;
						ui += Coulomb[iCoulomb].CoulombCoeff * operate_scalar_im;

						clear_kets();
						///////////////////////////////////////////////////
					}

					if (!ChopBool(ur) || !ChopBool(ui)) {
						writeoutput << "\n<\t" << i << "\t|H_int|\t" << j << "\t>=\t" << fixed << setprecision(20) << ur << "\t+I*\t" << ui;
						braHintKet_stable = false;
					}

				}
			}
		}
	}
	if (braHintKet_stable == true) {
		writeoutput << "\nOK.";
		writeoutput.flush();
	}
	*/

}

void Jm(signed long long int istart, signed long long int iend) {

	long double L = LSKetList[istart].L;
	long double S = LSKetList[istart].S;
	int J_count = (int)((L + S) - abs(L - S) + 1.0);
	vector<long double> J_series;
	for (int i = 0; i < J_count; i++) {
		J_series.push_back(L + S - (long double)i);
	}
	signed long long int isize = iend - istart + 1;

	signed long long int iJstart;
	for (signed long long int i = istart; i <= iend; i++) {
		if (ChopBool(LSKetList[i].Lz + LSKetList[i].Sz - J_series[0])) {
			iJstart = i;
			break;
		}
	}
	
	for (int i = 0; i < J_count; i++) {

		///////////////////////////////////////////////////
		clear_kets();

		if (i == 0) {
			InputKet = LSKetList[iJstart].ket;
		}
		else {
			
			GramSchmidtKet.clear();
			catched.clear();
			bool catchedlocked = false;

			for (signed long long int j = istart; j <= (signed)LSJKetList.size() - 1; j++) {
				if (ChopBool(LSJKetList[j].Lz + LSJKetList[j].Sz - J_series[i])) {
					catched.push_back(j);
					catchedlocked = true;
				}
			}

			if (catchedlocked == true) {
				global_istart = istart;
				global_iend = iend;
				GramSchmidt(1);
			}
			else {
				invalid_catchedlocked_exit();
			}
						
			InputKet = GramSchmidtKet[0];
		}
				
		for (int j = 0; j < (int)(2.0 * J_series[i] + 1.0); j++) {

			LSJKETLIST obj_LSJKet;

			if (j == 0) {
				obj_LSJKet.ket = InputKet;
			}
			else {
				clear_kets();
				InputKet = LSJKetList[(signed)LSJKetList.size() - 1].ket;
				operate(OperatorJm);
				normalized_simplifiedoutputket();
				obj_LSJKet.ket = NormalizedSimplifiedOutputKet;
			}
			
			///////////////////////////////////////////////////
			clear_kets();

			InputKet = obj_LSJKet.ket;
			InputBra = obj_LSJKet.ket;
			operate(OperatorSx);
			operate_scalar();
			obj_LSJKet.Sx = operate_scalar_re;

			clear_kets();

			InputKet = obj_LSJKet.ket;
			InputBra = obj_LSJKet.ket;
			operate(OperatorSy);
			operate_scalar();
			obj_LSJKet.Sy = operate_scalar_re;

			clear_kets();

			InputKet = obj_LSJKet.ket;
			InputBra = obj_LSJKet.ket;
			operate(OperatorSz);
			operate_scalar();
			obj_LSJKet.Sz = operate_scalar_re;

			clear_kets();

			InputKet = obj_LSJKet.ket;
			InputBra = obj_LSJKet.ket;
			operate(OperatorLx);
			operate_scalar();
			obj_LSJKet.Lx = operate_scalar_re;

			clear_kets();

			InputKet = obj_LSJKet.ket;
			InputBra = obj_LSJKet.ket;
			operate(OperatorLy);
			operate_scalar();
			obj_LSJKet.Ly = operate_scalar_re;

			clear_kets();

			InputKet = obj_LSJKet.ket;
			InputBra = obj_LSJKet.ket;
			operate(OperatorLz);
			operate_scalar();
			obj_LSJKet.Lz = operate_scalar_re;

			clear_kets();

			obj_LSJKet.Jx = obj_LSJKet.Lx + obj_LSJKet.Sx;
			obj_LSJKet.Jy = obj_LSJKet.Ly + obj_LSJKet.Sy;
			obj_LSJKet.Jz = obj_LSJKet.Lz + obj_LSJKet.Sz;

			obj_LSJKet.L = L;
			obj_LSJKet.S = S;
			obj_LSJKet.J = J_series[i];

			obj_LSJKet.CoulombEnerg = new long double[NCoulombHamiltonian];
			for (int ii = 0; ii < NCoulombHamiltonian; ii++) {
				obj_LSJKet.CoulombEnerg[ii] = LSKetList[(signed)LSJKetList.size()].CoulombEnerg[ii];
			}

			obj_LSJKet.CoulombEnergySum = LSKetList[(signed)LSJKetList.size()].CoulombEnergySum;

			clear_kets();

			InputKet = obj_LSJKet.ket;
			InputBra = obj_LSJKet.ket;
			operate(OperatorSOC);
			operate_scalar();
			obj_LSJKet.SOCEnerg = operate_scalar_re;

			clear_kets();
			///////////////////////////////////////////////////

			LSJKetList.push_back(obj_LSJKet);
			
		}
		
		clear_kets();
		///////////////////////////////////////////////////

	}

	for (signed long long int ii = istart; ii <= iend; ii++) {

		if (EigenvectorOutputMode == 0) {
			writeoutput << ii << "\t||\tL=\t" << fixed << setprecision(1) << LSJKetList[ii].L << "\t||\tS=\t" << LSJKetList[ii].S << "\t||\tJ=\t" << LSJKetList[ii].J << "\t||\tJz=\t" << LSJKetList[ii].Jz << "\t||\tEstimatedEnrg=\t";
			writeoutput.flush();
		}
		else {
			writeoutput << ii << "\t||\tL=\t" << fixed << setprecision(1) << LSJKetList[ii].L << "\t||\tS=\t" << LSJKetList[ii].S << "\t||\tJ=\t" << LSJKetList[ii].J << "\t||\tJz=\t" << LSJKetList[ii].Jz << "\n";
			writeoutput.flush();
			fwrite(&ii, sizeof(signed long long int), 1, _writeoutput);
			fwrite(&LSJKetList[ii].L, sizeof(long double), 1, _writeoutput);
			fwrite(&LSJKetList[ii].S, sizeof(long double), 1, _writeoutput);
			fwrite(&LSJKetList[ii].J, sizeof(long double), 1, _writeoutput);
			fwrite(&LSJKetList[ii].Jz, sizeof(long double), 1, _writeoutput);
		}
				
		long double enrg = 0.0;
		for (int iCoulomb = 0; iCoulomb < NCoulombHamiltonian; iCoulomb++) {
			if (EigenvectorOutputMode == 0) {
				writeoutput << fixed << setprecision(20) << LSJKetList[ii].CoulombEnerg[iCoulomb] << "\t" << Coulomb[iCoulomb].CoulombName << "\t(=\t" << Coulomb[iCoulomb].CoulombCoeff << "\t)";
				writeoutput.flush();
			}
			enrg += LSJKetList[ii].CoulombEnerg[iCoulomb] * Coulomb[iCoulomb].CoulombCoeff;
			if (EigenvectorOutputMode == 0) {
				writeoutput << "\t+\t";
				writeoutput.flush();
			}
		}
		LSJKetList[ii].CoulombEnergySum = enrg;

		LSJKetList[ii].CoulombAndSOCEnerg = LSJKetList[ii].CoulombEnergySum + LSJKetList[ii].SOCEnerg * SOC_coefficient;
		if (EigenvectorOutputMode == 0) {
			writeoutput << fixed << setprecision(20) << LSJKetList[ii].SOCEnerg << "\t" << SOC_coefficient_name << "\t(=\t" << SOC_coefficient << "\t)=\t" << LSJKetList[ii].CoulombAndSOCEnerg << "\t||\tNbases=\t" << (signed)LSJKetList[ii].ket.KetElement.size() << "\t||\t";
			writeoutput.flush();
		}
		else {
			signed long long int size_buf = (signed)LSJKetList[ii].ket.KetElement.size();
			fwrite(&LSJKetList[ii].CoulombAndSOCEnerg, sizeof(long double), 1, _writeoutput);
			fwrite(&size_buf, sizeof(signed long long int), 1, _writeoutput);
		}
		
		for (int j = 0; j < (signed)LSJKetList[ii].ket.KetElement.size(); j++) {

			if (EigenvectorOutputMode == 0) {
				writeoutput << fixed << setprecision(20) << LSJKetList[ii].ket.KetElement[j].coeffcientRe << "\t" << LSJKetList[ii].ket.KetElement[j].coeffcientIm << "\t";
				writeoutput.flush();
			}
			else {
				fwrite(&LSJKetList[ii].ket.KetElement[j].coeffcientRe, sizeof(long double), 1, _writeoutput);
				fwrite(&LSJKetList[ii].ket.KetElement[j].coeffcientIm, sizeof(long double), 1, _writeoutput);
			}

			Secondization(LSJKetList[ii].ket.KetElement[j].base);
			if (EigenvectorOutputMode == 0) {
				for (int k = 0; k < NOrbitalRoom; k++) {
					writeoutput << SecondQuantizedKet[k];
					writeoutput.flush();
				}
				if (j != (signed)LSJKetList[ii].ket.KetElement.size() - 1) {
					writeoutput << "\t|\t";
					writeoutput.flush();
				}
				else {
					writeoutput << "\t||\n";
					writeoutput.flush();
				}
			}
			else {
				fwrite(&SecondQuantizedKet[0], sizeof(int), NOrbitalRoom, _writeoutput);
			}

		}
	
	}
	
}

void prtb_with_SOC() {

	writeoutput << "\n==========================prtb_with_SOC()==========================\n";
	writeoutput.flush();

	for (int iLSBlock = 0; iLSBlock < (signed)LSBlock.size(); iLSBlock++) {
		signed long long int BlockSize = LSBlock[iLSBlock].iend - LSBlock[iLSBlock].istart + 1;
		signed long long int BlockSizeElement = BlockSize / (signed long long int)LSBlock[iLSBlock].multiplicity;
		for (signed long long int ii = 0; ii < (signed long long int)LSBlock[iLSBlock].multiplicity; ii++) {
			Jm(LSBlock[iLSBlock].istart + ii * BlockSizeElement, LSBlock[iLSBlock].istart + (ii + 1) * BlockSizeElement - 1);
		}	
	}

	writeoutput << "****************************prtb_with_SOC()****************************\n";
	writeoutput.flush();

	for (signed long long int ii = 0; ii < (signed long long int)LSJKetList.size(); ii++) {

		if (EigenvectorOutputMode == 0) {
			writeoutput << ii << "\t||\tL=\t" << fixed << setprecision(1) << LSJKetList[ii].L << "\t||\tS=\t" << LSJKetList[ii].S << "\t||\tJ=\t" << LSJKetList[ii].J << "\t||\tJz=\t" << LSJKetList[ii].Jz << "\t||\tEstimatedEnrg=\t";
			writeoutput.flush();
		}
		else {
			writeoutput << ii << "\t||\tL=\t" << fixed << setprecision(1) << LSJKetList[ii].L << "\t||\tS=\t" << LSJKetList[ii].S << "\t||\tJ=\t" << LSJKetList[ii].J << "\t||\tJz=\t" << LSJKetList[ii].Jz << "\n";
			writeoutput.flush();
			fwrite(&ii, sizeof(signed long long int), 1, _writeoutput);
			fwrite(&LSJKetList[ii].L, sizeof(long double), 1, _writeoutput);
			fwrite(&LSJKetList[ii].S, sizeof(long double), 1, _writeoutput);
			fwrite(&LSJKetList[ii].J, sizeof(long double), 1, _writeoutput);
			fwrite(&LSJKetList[ii].Jz, sizeof(long double), 1, _writeoutput);
		}

		if (EigenvectorOutputMode == 0) {
			for (int iCoulomb = 0; iCoulomb < NCoulombHamiltonian; iCoulomb++) {
				writeoutput << fixed << setprecision(20) << LSJKetList[ii].CoulombEnerg[iCoulomb] << "\t" << Coulomb[iCoulomb].CoulombName << "\t(=\t" << Coulomb[iCoulomb].CoulombCoeff << "\t)";
				writeoutput.flush();
				writeoutput << "\t+\t";
				writeoutput.flush();
			}

			writeoutput << fixed << setprecision(20) << LSJKetList[ii].SOCEnerg << "\t" << SOC_coefficient_name << "\t(=\t" << SOC_coefficient << "\t)=\t" << LSJKetList[ii].CoulombAndSOCEnerg << "\t||\t";
			writeoutput.flush();
		}
		else {
			fwrite(&LSJKetList[ii].CoulombAndSOCEnerg, sizeof(long double), 1, _writeoutput);
		}

		vector<long double> CubicHarmonicRe;
		vector<long double> CubicHarmonicIm;
		vector<signed long long int> CubicHarmonicBases;

		for (signed long long int jj = 0; jj < (signed long long int)listRepetitiveNCCubicHarmonic.size(); jj++) {

			///////////////////////////////////////////////////

			clear_kets();

			InputKet = LSJKetList[ii].ket;
			InputBra = listRepetitiveNCCubicHarmonic[jj].ket;
			_operate_scalar();

			if (!ChopBool(operate_scalar_re) || !ChopBool(operate_scalar_im)) {
				CubicHarmonicRe.push_back(operate_scalar_re);
				CubicHarmonicIm.push_back(operate_scalar_im);
				CubicHarmonicBases.push_back(jj);
			}

			clear_kets();

			///////////////////////////////////////////////////

		}
		
		if (EigenvectorOutputMode == 0) {
			writeoutput << "*Nbases=\t" << (signed)CubicHarmonicBases.size() << "\t||\t";
			writeoutput.flush();
		}
		else {
			signed long long int size_buf = (signed)CubicHarmonicBases.size();
			fwrite(&size_buf, sizeof(signed long long int), 1, _writeoutput);
		}

		for (signed long long int jj = 0; jj < (signed)CubicHarmonicBases.size(); jj++) {

			if (EigenvectorOutputMode == 0) {
				writeoutput << fixed << setprecision(20) << CubicHarmonicRe[jj] << "\t" << CubicHarmonicIm[jj] << "\t";
				writeoutput.flush();
			}
			else {
				fwrite(&CubicHarmonicRe[jj], sizeof(long double), 1, _writeoutput);
				fwrite(&CubicHarmonicIm[jj], sizeof(long double), 1, _writeoutput);
			}
			if (EigenvectorOutputMode == 0) {
				writeoutput << listRepetitiveNCCubicHarmonic[CubicHarmonicBases[jj]].starket_rep;
				writeoutput.flush();
				if (jj != (signed)CubicHarmonicBases.size() - 1) {
					writeoutput << "\t|\t";
					writeoutput.flush();
				}
				else {
					writeoutput << "\t||\n";
					writeoutput.flush();
				}
			}
			else {
				fwrite(&listRepetitiveNCCubicHarmonic[CubicHarmonicBases[jj]].starket, sizeof(signed long long int), 1, _writeoutput);
			}

		}
		
	}

	writeoutput << "==========================prtb_with_SOC()==========================\n";
	writeoutput.flush();
	writeoutput << LSJKetList.size() << " LSJ-coupling bases (all in a ket form) are created; ";
	writeoutput.flush();
	if ((signed long long int)LSJKetList.size() == nCr((signed long long int)NOrbitalRoom, (signed long long int)NOirbtalOccupied)) {
		writeoutput << "as expected.";
		writeoutput.flush();
	}
	else {
		writeoutput << "something wrong! (signed)LSJKetList.size() != nCr(NOrbitalRoom, NOirbtalOccupied).";
		writeoutput.flush();
	}

	LSJBLOCK obj_block;
	obj_block.istart = 0;
	for (signed long long int i = 0; i < (signed long long int)LSJKetList.size() - 1; i++) {
		
		int LSBlock_Buf[2];
		for (int j = 0; j < (signed)LSBlock.size(); j++) {
			if (i >= LSBlock[j].istart && i <= LSBlock[j].iend) {
				LSBlock_Buf[0] = j;
			}
			if (i + 1 >= LSBlock[j].istart && i + 1 <= LSBlock[j].iend) {
				LSBlock_Buf[1] = j;
			}
		}

		bool LSBlock_cross = false;
		bool LSBlock_subcross = false;
		if (LSBlock_Buf[0] != LSBlock_Buf[1]) {
			LSBlock_cross = true;
		}
		else {
			if (LSBlock[LSBlock_Buf[0]].multiplicity > 1) {
				int jcount = 0;
				int jmultiplictycount = 0;
				int jmultiplicty_buf[2];
				for (int j = LSBlock[LSBlock_Buf[0]].istart; j <= (signed)LSBlock[LSBlock_Buf[0]].iend; j++) {
					jcount++;
					if (i == j) {
						jmultiplicty_buf[0] = jmultiplictycount;
					}
					if (i + 1 == j) {
						jmultiplicty_buf[1] = jmultiplictycount;
					}
					if (jcount == (int)(2.0 * LSBlock[LSBlock_Buf[0]].L + 1.0) * (int)(2.0 * LSBlock[LSBlock_Buf[0]].S + 1.0)) {
						jcount = 0;
						jmultiplictycount++;
					}
				}
				if (jmultiplicty_buf[0] != jmultiplicty_buf[1]) {
					LSBlock_subcross = true;
				}
			}
		}
		
		if (!ChopBool(LSJKetList[i].L - LSJKetList[i + 1].L) || !ChopBool(LSJKetList[i].S - LSJKetList[i + 1].S) || !ChopBool(LSJKetList[i].J - LSJKetList[i + 1].J) || LSBlock_cross == true || LSBlock_subcross == true) {
			obj_block.iend = i;
			obj_block.L = LSJKetList[i].L;
			obj_block.S = LSJKetList[i].S;
			obj_block.J = LSJKetList[i].J;
			LSJBlock.push_back(obj_block);
			obj_block.istart = i + 1;
		}

	}
	obj_block.iend = (signed long long int)LSJKetList.size() - 1;
	obj_block.L = LSJKetList[(signed long long int)LSJKetList.size() - 1].L;
	obj_block.S = LSJKetList[(signed long long int)LSJKetList.size() - 1].S;
	obj_block.J = LSJKetList[(signed long long int)LSJKetList.size() - 1].J;
	LSJBlock.push_back(obj_block);
	
	writeoutput << "\nLSJBlock\n";
	writeoutput.flush();
	for (int i = 0; i < (signed)LSJBlock.size(); i++) {
		writeoutput << i << "\tL=\t" << fixed << setprecision(1) << LSJBlock[i].L << "\tS=\t" << LSJBlock[i].S << "\tJ=\t" << LSJBlock[i].J << "\tistart=\t" << LSJBlock[i].istart << "\tiend=\t" << LSJBlock[i].iend << "\t";
		writeoutput.flush();
		if ((int)(2.0 * LSJBlock[i].J + 1) == LSJBlock[i].iend - LSJBlock[i].istart + 1) {
			writeoutput << "as expected";
			writeoutput.flush();
		}
		else {
			writeoutput << "wrong count!";
			writeoutput.flush();
		}
		for (int j = LSJBlock[i].istart; j <= LSJBlock[i].iend; j++) {
			if (abs(LSJKetList[j].Jz) > LSJBlock[i].J + 2.0 * Chop) {
				writeoutput << "\twrong Jz included!(\t" << j << "\t)";
				writeoutput.flush();
			}
		}
		writeoutput << "\tEstimatedEnrg=\t";
		writeoutput.flush();
		for (int iCoulomb = 0; iCoulomb < NCoulombHamiltonian; iCoulomb++) {
			writeoutput << fixed << setprecision(20) << LSJKetList[LSJBlock[i].istart].CoulombEnerg[iCoulomb] << "\t" << Coulomb[iCoulomb].CoulombName << "\t(=\t" << Coulomb[iCoulomb].CoulombCoeff << "\t)";
			writeoutput.flush();
			writeoutput << "\t+\t";
			writeoutput.flush();
		}
		writeoutput << fixed << setprecision(20) << LSJKetList[LSJBlock[i].istart].SOCEnerg << "\t" << SOC_coefficient_name << "\t(=\t" << SOC_coefficient << "\t)=\t" << LSJKetList[LSJBlock[i].istart].CoulombAndSOCEnerg;
		writeoutput.flush();
		if (i != (signed)LSJBlock.size() - 1) {
			writeoutput << "\n";
			writeoutput.flush();
		}
	}

	writeoutput << "\n==========================prtb_with_SOC()==========================";
	writeoutput.flush();

}

int main(int argc, char *argv[]) {


	ierr = MPI_Init(&argc, &argv);
	if (ierr != 0)
	{
		cout << "  MPI_Init returned nonzero ierr.\n";
		exit(1);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &pr_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &pr_id);

	read_para();
	
	read_SLlist();
	read_operators();

	create_vac();
	
	pr_job_count = new int[pr_size];
	pr_iconditionstart = new int[pr_size];
	pr_iconditionend = new int[pr_size];
	for (int i = 0; i < pr_size; i++) {
		pr_job_count[i] = 0;
		pr_iconditionstart[i] = 0;
		pr_iconditionend[i] = 0;
	}

	int pr_job_count_round = 0;
	for (int i = 0; i < NConditions; i++) {
		pr_job_count[pr_job_count_round]++;
		pr_job_count_round++;
		if (pr_job_count_round == pr_size) {
			pr_job_count_round = 0;
		}
	}

	for (int i = 0; i < pr_size; i++) {
		for (int j = 0; j < i; j++) {
			pr_iconditionstart[i] += pr_job_count[j];
		}
	}

	for (int i = 0; i < pr_size - 1; i++) {
		pr_iconditionend[i] = pr_iconditionstart[i + 1] - 1;
	}
	pr_iconditionend[pr_size - 1] = NConditions - 1;

	for (int pr_job = pr_iconditionstart[pr_id]; pr_job <= pr_iconditionend[pr_id]; pr_job++) {

		writeoutput.open(OutputFilePrefix + "_" + Conditions[pr_job].condition_name + ".txt");

		writeoutput << "pr_distribution\npr_id\tpr_iconditionstart[pr_id]\tpr_iconditionend[pr_id]\n";
		writeoutput.flush();

		for (int i = 0; i < pr_size; i++) {
			writeoutput << i << "\t" << pr_iconditionstart[i] + 1 << "\t" << pr_iconditionend[i] + 1;
			writeoutput.flush();
			if (i == pr_id) {
				writeoutput << "\t***";
				writeoutput.flush();
			}
			if (i != pr_size - 1) {
				writeoutput << "\n";
				writeoutput.flush();
			}
		}

		if (EigenvectorOutputMode == 1) {

			_writeoutput_filename = OutputFilePrefix + "_" + Conditions[pr_job].condition_name + ".dat";

#ifdef _WIN32
			errno_t error_fwrite;
			error_fwrite = fopen_s(&_writeoutput, _writeoutput_filename.c_str(), "wb");
#elif _WIN64
			errno_t error_fwrite;
			error_fwrite = fopen_s(&_writeoutput, _writeoutput_filename.c_str(), "wb");
#else
			_writeoutput = fopen(_writeoutput_filename.c_str(), "wb");
#endif

		}

		writeoutput << "\nConditions...\n";
		writeoutput.flush();
		writeoutput << Conditions[pr_job].condition_name << "\n";
		writeoutput.flush();
		for (int i = 0; i < NCoulombHamiltonian; i++) {
			Coulomb[i].CoulombCoeff = Conditions[pr_job].CoulombCoeff[i];
			writeoutput << Coulomb[i].CoulombName << "\t" << fixed << setprecision(20) << Coulomb[i].CoulombCoeff << "\n";
			writeoutput.flush();
		}
		SOC_coefficient = Conditions[pr_job].SOCCoeff;
		writeoutput << SOC_coefficient_name << "\t" << fixed << setprecision(20) << SOC_coefficient << "\n=========================================\n";
		writeoutput.flush();

		listRepetitiveNC.clear();
		listRepetitiveNCCubicHarmonic.clear();
		LSKetList.clear();
		LSBlock.clear();
		LSJKetList.clear();
		LSJBlock.clear();

		write_subset();
		roll_the_dice();
		prtb_with_SOC();
		writeoutput.close();

		if (EigenvectorOutputMode == 1) {
			fclose(_writeoutput);
		}

	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;

}