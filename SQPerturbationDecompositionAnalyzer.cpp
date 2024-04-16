#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <cmath>
#include <mpi.h>

typedef struct {

	long double L, S, J, absJz;
	std::vector<int> whichindex;
	std::string symbol_long, symbol_short;

}LSJJz;

std::vector<LSJJz> Intm1LSJJz;
std::vector<LSJJz> Intm2LSJJz;

typedef struct {

	int Intm1LSJJzIndex;
	int Intm2LSJJzIndex;
	long double local_v, global_v;
	bool counted;
	std::vector<std::string> labels;

}DOUBLELSJJz;

std::vector<DOUBLELSJJz> DoubleLSJJz;
std::vector<int> sDoubleLSJJz;

int pr_size, ierr, pr_id;
auto startefr = std::chrono::system_clock::now();

std::string CouplingConstFilename;
long double SizeMinimumLimit;
std::string IntmFilename1;
int MaximumBlocks1;
std::string IntmFilename2;
int MaximumBlocks2;
std::string OutputFilename;
int VerboseMode;
int MonitoringIndex;

std::string sbuf;
std::string couplingconst_name;
long double couplingconst;

std::ifstream readpara, readIntmFile, readdecomposition;
std::ofstream writeoutput;

bool ChopBool(long double v, long double choplevel) {

	if (abs(v) < choplevel) {
		return true;
	}
	else {
		return false;
	}

}

std::string to_string_with_precision(long double value) {

	std::ostringstream out;
	out.precision(0);
	if (ChopBool(remainder(value, 1.0), 0.0001) == true) {
		out << std::fixed << value;
	}
	else {
		out << std::fixed << 2.0 * value << "/2";
	}
	return out.str();

}

std::string symbolL(int L) {

	if (L == 0) {
		return "S";
	}
	else if (L == 1) {
		return "P";
	}
	else if (L == 2) {
		return "D";
	}
	else if (L == 3) {
		return "F";
	}
	else if (L == 4) {
		return "G";
	}
	else if (L == 5) {
		return "H";
	}
	else if (L == 6) {
		return "I";
	}
	else if (L == 7) {
		return "K";
	}
	else if (L == 8) {
		return "L";
	}
	else if (L == 9) {
		return "M";
	}
	else if (L == 10) {
		return "N";
	}
	else if (L == 11) {
		return "O";
	}
	else if (L == 12) {
		return "Q";
	}
	else {
		return "X";
	}

}

void read_para() {

	readpara.open("basic_properties_SQPerturbationDecompositionAnalyzer.txt");

	readpara >> sbuf >> CouplingConstFilename >> SizeMinimumLimit;
	readpara >> sbuf >> IntmFilename1 >> MaximumBlocks1;
	readpara >> sbuf >> IntmFilename2 >> MaximumBlocks2;
	readpara >> sbuf >> OutputFilename;
	readpara >> sbuf >> VerboseMode;
	readpara >> sbuf >> MonitoringIndex;

	readpara.close();

}

void read_IntmFile(int which) {

	int MaximumBlocks;
	if (which == 1) {
		readIntmFile.open(IntmFilename1);
		MaximumBlocks = MaximumBlocks1;
	}
	else {
		readIntmFile.open(IntmFilename2);
		MaximumBlocks = MaximumBlocks2;
	}

	while (true) {
		readIntmFile >> sbuf;
		if (sbuf == "LSJBlock") {
			break;
		}
	}
	getline(readIntmFile, sbuf);

	int NBlocks = 0;
	while (true) {
		getline(readIntmFile, sbuf);
		if (sbuf != "==========================prtb_with_SOC()==========================") {
			long double L, S, J;
			int istart, iend;
			std::stringstream lcstm;
			lcstm << sbuf;
			lcstm >> sbuf >> sbuf >> L >> sbuf >> S >> sbuf >> J >> sbuf >> istart >> sbuf >> iend;
			long double Jz = J;
			int iistart = istart;
			int iiend = iend;
			while (true) {
				LSJJz objLSJJz;
				objLSJJz.L = L;
				objLSJJz.S = S;
				objLSJJz.J = J;
				objLSJJz.absJz = Jz;
				objLSJJz.whichindex.push_back(iistart);
				objLSJJz.whichindex.push_back(iiend);
				std::string sS = to_string_with_precision(2 * S + 1);
				std::string sL = symbolL(L);
				std::string sJ = to_string_with_precision(J);
				std::string sJz = to_string_with_precision(Jz);
				objLSJJz.symbol_long = "2S+1=\t" + sS + "\tL=\t" + sL + "\tJ=\t" + sJ + "\t|Jz|=\t" + sJz + "\t";
				objLSJJz.symbol_short = sS + sL + sJ + ":|Jz|=" + sJz;
				if (which == 1) {
					Intm1LSJJz.push_back(objLSJJz);
				}
				else {
					Intm2LSJJz.push_back(objLSJJz);
				}
				iistart += 1;
				iiend -= 1;
				Jz -= 1.0;
				if (ChopBool(Jz, 0.0001) == false && Jz < 0.0) {
					break;
				}
			}
			NBlocks++;
			if (MaximumBlocks >=0 && NBlocks > MaximumBlocks) {
				break;
			}
		}
		else {
			break;
		}
	}

	readIntmFile.close();

}

void get_VACFULL(int VAC0FULL1, int which) {

	LSJJz objLSJJz;
	objLSJJz.L = 0.0;
	objLSJJz.S = 0.0;
	objLSJJz.J = 0.0;
	objLSJJz.absJz = 0.0;
	objLSJJz.whichindex.push_back(0);
	objLSJJz.whichindex.push_back(0);
	if (VAC0FULL1 == 0) {
		objLSJJz.symbol_long = "VAC";
		objLSJJz.symbol_short = "VAC";
	}
	else {
		objLSJJz.symbol_long = "FULL";
		objLSJJz.symbol_short = "FULL";
	}
	if (which == 1) {
		Intm1LSJJz.push_back(objLSJJz);
	}
	else {
		Intm2LSJJz.push_back(objLSJJz);
	}

}

void create_doubleLSJJz() {

	for (int i = 0; i < (signed)Intm1LSJJz.size(); i++) {
		for (int j = 0; j < (signed)Intm2LSJJz.size(); j++) {
			DOUBLELSJJz objDOUBLELSJJz;
			objDOUBLELSJJz.Intm1LSJJzIndex = i;
			objDOUBLELSJJz.Intm2LSJJzIndex = j;
			objDOUBLELSJJz.local_v = 0.0;
			objDOUBLELSJJz.global_v = 0.0;
			objDOUBLELSJJz.counted = false;
			for (int ii = 0; ii < (signed)Intm1LSJJz[i].whichindex.size(); ii++) {
				for (int jj = 0; jj < (signed)Intm2LSJJz[j].whichindex.size(); jj++) {
					std::string label1 = "1-0-0-" + std::to_string(Intm1LSJJz[i].whichindex[ii]) + "^1-0-1-" + std::to_string(Intm2LSJJz[j].whichindex[jj]);
					std::string label2 = "1-1-0-" + std::to_string(Intm2LSJJz[j].whichindex[jj]) + "^1-1-1-" + std::to_string(Intm1LSJJz[i].whichindex[ii]);
					objDOUBLELSJJz.labels.push_back(label1);
					objDOUBLELSJJz.labels.push_back(label2);
				}
			}
			DoubleLSJJz.push_back(objDOUBLELSJJz);
		}
	}

}

void read_decomposition() {

	readdecomposition.open(CouplingConstFilename);

	readdecomposition >> couplingconst_name >> couplingconst;
	int linecount = 0;
	int finitecount = 0;
	while (true) {
		if (readdecomposition.eof()) {
			break;
		}
		long double ldbuf;
		readdecomposition >> sbuf >> ldbuf;
		if (ldbuf != 0.0) {
			finitecount++;
		}
		linecount++;
	}
	if (pr_id == 0) {
		writeoutput << linecount << "\tlines_detected...\t" << finitecount << "\tfinite_lines_detected...\n";
		writeoutput.flush();
	}

	readdecomposition.close();

	int* pr_Njob, * pr_istart, * pr_iend;
	pr_Njob = new int[pr_size];
	pr_istart = new int[pr_size];
	pr_iend = new int[pr_size];
	for (int ipr = 0; ipr < pr_size; ipr++) {
		pr_Njob[ipr] = 0;
		pr_istart[ipr] = 0;
		pr_iend[ipr] = 0;
	}

	int ipr_point = 0;
	for (int ipr = 0; ipr < linecount; ipr++) {
		pr_Njob[ipr_point]++;
		ipr_point++;
		if (ipr_point == pr_size) {
			ipr_point = 0;
		}
	}

	for (int ipr = 0; ipr < pr_size; ipr++) {
		if (pr_Njob[ipr] > 0) {
			if (ipr > 0) {
				pr_istart[ipr] = pr_iend[ipr - 1] + 1;
			}
			pr_iend[ipr] = pr_istart[ipr] + pr_Njob[ipr] - 1;
		}
	}

	if (VerboseMode == 1 && pr_id == 0) {
		writeoutput << "pr_id\tpr_Njob\tpr_istart\tpr_iend\n";
		for (int ipr = 0; ipr < pr_size; ipr++) {
			writeoutput << ipr << "\t" << pr_Njob[ipr] << "\t" << pr_istart[ipr] + 1 << "\t" << pr_iend[ipr] + 1 << "\n";
			writeoutput.flush();
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (pr_Njob[pr_id] > 0) {
	
		readdecomposition.open(CouplingConstFilename);
		readdecomposition >> couplingconst_name >> couplingconst;

		int count = 0;
		while (true) {
			readdecomposition >> sbuf;
			if (!readdecomposition.eof()) {

				long double ldbuf;
				readdecomposition >> ldbuf;

				bool pass = false;
				if (SizeMinimumLimit < 0) {
					pass = true;
				}
				else {
					if (abs(ldbuf) >= SizeMinimumLimit) {
						pass = true;
					}
				}

				if (count >= pr_istart[pr_id] && count <= pr_iend[pr_id]) {
					bool foundlabel = false;
					for (int i = 0; i < (signed)DoubleLSJJz.size(); i++) {
						for (int j = 0; j < (signed)DoubleLSJJz[i].labels.size(); j++) {
							if (sbuf == DoubleLSJJz[i].labels[j]) {
								DoubleLSJJz[i].local_v += ldbuf;
								foundlabel = true;
								break;
							}
						}
						if (foundlabel == true) {
							break;
						}
					}
					if (VerboseMode == 1 && foundlabel == false) {
						std::cout << "[ERROR]\tnot_found:\t" << sbuf << "\t" << ldbuf << "\n";
						std::cout.flush();
					}
					if (pr_id == 0 && VerboseMode == 1 && count + 1 % MonitoringIndex == 0) {
						auto endefr = std::chrono::system_clock::now();
						std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
						writeoutput << count + 1 << "\tdetected...\t" << elapsed_secondsfr.count() << "\n";
						writeoutput.flush();
					}
				}

				count++;
				
			}
			else {
				break;
			}
		}

		readdecomposition.close();
	
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for (int i = 0; i < (signed)DoubleLSJJz.size(); i++) {
		MPI_Allreduce(&DoubleLSJJz[i].local_v, &DoubleLSJJz[i].global_v, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

}

void sort_decomposition() {

	int NCount = 0;
	for (int i = 0; i < (signed)DoubleLSJJz.size(); i++) {
		if (DoubleLSJJz[i].global_v == 0.0) {
			DoubleLSJJz[i].counted = true;
		}
		else {
			NCount++;
		}
	}
	
	int count = 0;
	while (true) {
		long double max = -100000000000000000.0;
		int whichmax;
		for (int i = 0; i < (signed)DoubleLSJJz.size(); i++) {
			if (DoubleLSJJz[i].counted == false) {
				if (DoubleLSJJz[i].global_v > max) {
					whichmax = i;
					max = DoubleLSJJz[i].global_v;
				}
			}
		}
		sDoubleLSJJz.push_back(whichmax);
		DoubleLSJJz[whichmax].counted = true;
		count++;
		if (count == NCount) {
			break;
		}
	}

}

void write_decomposition() {

	if (VerboseMode == 1) {
		writeoutput << "======================================================================\n";
		writeoutput.flush();
	}
	writeoutput << couplingconst_name << "\t" << couplingconst << "\n";
	writeoutput.flush();

	for (int i = 0; i < (signed)sDoubleLSJJz.size(); i++) {
		int ii = sDoubleLSJJz[i];
		int ii1 = DoubleLSJJz[ii].Intm1LSJJzIndex;
		int ii2 = DoubleLSJJz[ii].Intm2LSJJzIndex;
		writeoutput << Intm1LSJJz[ii1].symbol_short << "\t" << Intm2LSJJz[ii2].symbol_short << "\t" << DoubleLSJJz[ii].global_v << "\t" << Intm1LSJJz[ii1].symbol_long << "\t" << Intm2LSJJz[ii2].symbol_long;
		writeoutput.flush();
		if (i != (signed)sDoubleLSJJz.size() - 1) {
			writeoutput << "\n";
			writeoutput.flush();
		}
	}

}

int main(int argc, char* argv[]) {

	ierr = MPI_Init(&argc, &argv);
	if (ierr != 0)
	{
		std::cout << "  MPI_Init returned nonzero ierr.\n";
		exit(1);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &pr_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &pr_id);

	read_para();
	if (pr_id == 0) {
		writeoutput.open(OutputFilename);
		writeoutput << std::fixed << std::setprecision(20);
	}
	
	if (IntmFilename1 != "VAC" && IntmFilename1 != "FULL") {
		read_IntmFile(1);
	}
	else if (IntmFilename1 == "VAC") {
		get_VACFULL(0, 1);
	}
	else {
		get_VACFULL(1, 1);
	}
	if (IntmFilename2 != "VAC" && IntmFilename2 != "FULL") {
		read_IntmFile(2);
	}
	else if (IntmFilename2 == "VAC") {
		get_VACFULL(0, 2);
	}
	else {
		get_VACFULL(1, 2);
	}
	create_doubleLSJJz();
	MPI_Barrier(MPI_COMM_WORLD);
	
	read_decomposition();
	MPI_Barrier(MPI_COMM_WORLD);
	sort_decomposition();
	MPI_Barrier(MPI_COMM_WORLD);
	if (pr_id == 0) {
		write_decomposition();
		writeoutput.close();
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	return 0;

}