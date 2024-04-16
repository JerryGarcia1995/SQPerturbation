#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <mpi.h>
#include <cmath>

typedef struct {

	long double coeffcientRe;
	long double coeffcientIm;
	std::vector<int> CRList;
	std::vector<int> ANList;

}OPERATORLINE;

typedef struct {

	std::vector<OPERATORLINE> OperatorLine;

}OPERATOR;

std::vector<OPERATOR> Ht;

typedef struct {

	long double coeffcientRe;
	long double coeffcientIm;
	signed long long int base;

}KETELEMENT;

typedef struct {

	std::vector<KETELEMENT> KetElement;

}KET;

typedef struct {

	KET ket;
	long double energy;
	unsigned long long int label;

}KETS;

typedef struct {

	std::vector<KETS> kets;

}KETSS;

typedef struct {

	std::vector<KETSS> ketss;
	int* NOs;

}KETSSS;

typedef struct {

	std::vector<KETSSS> ketsss;
	std::vector<std::string> labels;

}KETSERIES;

std::vector<KETSERIES> ketseries;

typedef struct {

	std::vector<KETS>* ketsptrs;

}KETSPTRS;

typedef struct {

	long double* v;

}PTRB;

typedef struct {

	std::vector<PTRB> ptrb;

}PTRBMATRIX;

typedef struct {

	std::vector<PTRBMATRIX> ptrbmatrix;

}PTRBMATRIXSRS;

std::vector<PTRBMATRIXSRS> ptrbmatrixsrs;

std::vector<PTRB> ptrbbuf;

std::vector<signed long long int> bracomb, ketcomb, Fdiagram;

typedef struct {

	long double coefficient;
	int index1, index2, reim;

}COMPOSITE;

typedef struct {

	std::string TermName;
	std::vector<COMPOSITE> composite;

}USERDEFINEDTERM;

std::vector<USERDEFINEDTERM> userdefinedterm;

typedef struct {

	std::vector<std::string> strterm;
	long double reim[2];

}INTERMEDIATECOMP;

typedef struct {

	std::vector<INTERMEDIATECOMP> intermediatecomp;

}INTERMEDIATEANALYSIS;

INTERMEDIATEANALYSIS** intermeidateanalysis;

int pr_size, ierr, pr_id;
auto startefr = std::chrono::system_clock::now();
auto endefr_buf = std::chrono::system_clock::now();

std::string OutputFilename;
int OutputPrecision;
int ChopLevel;
int ReadPrtb;
int WritePrtb;
int NOrbitalRooms;
std::string HoppintHamiltonianFilename;
int ManualMPI;
std::string ManualMPIFilename;
int SaveMemory;
int GetLabel;
int ReadSHorCH;
signed long long int MonitoringIndex;
int FullMonitoring;
int NUserdefinedTerms;

long double Chop;
std::string sbuf;
std::vector<std::string> globalketstring;
std::vector<std::string> globalbrastring;
int globalketstring_max;
int globalbrastring_max;
bool ptrbmatrixaddup;
long double*** final_prtb_terms;
long double*** global_final_prtb_terms;
bool mpi_activated = true;
signed long long int prtb_count = 0;
bool just_thrown = true;

std::ifstream readpara, readmanualmpi, readoutput_sub;
std::ofstream writeoutput, writeoutput_sub, writeoutput_decomp;

FILE* read_prtb = nullptr;
FILE* write_prtb = nullptr;

void ComplexMTwoRe(long double* a1, long double* a2, long double* b1, long double* b2, long double* answer) {
	*answer = *a1 * *b1 - *a2 * *b2;
}

void ComplexMTwoIm(long double* a1, long double* a2, long double* b1, long double* b2, long double* answer) {
	*answer = *a2 * *b1 + *a1 * *b2;
}

void ComplexMThreeRe(long double* a1, long double* a2, long double* b1, long double* b2, long double* c1, long double* c2, long double* answer) {
	*answer = *a1 * *b1 * *c1 - *a2 * *b2 * *c1 - *a2 * *b1 * *c2 - *a1 * *b2 * *c2;
}

void ComplexMThreeIm(long double* a1, long double* a2, long double* b1, long double* b2, long double* c1, long double* c2, long double* answer) {
	*answer = *a2 * *b1 * *c1 + *a1 * *b2 * *c1 + *a1 * *b1 * *c2 - *a2 * *b2 * *c2;
}

bool ChopBool(long double* v) {

	if (fabs(*v) < Chop) {
		return true;
	}
	else {
		return false;
	}

}

void ConvertKetNumber(int* ConvertedQuantizedKet, int* NO, signed long long int* Bin) {

	*Bin = 0;

	for (int i = 0; i < *NO; i++) {
		*Bin += ConvertedQuantizedKet[*NO - (i + 1)] * (int)pow(2, i);
	}

}

void Secondize(signed long long int* ketNumber, int* NO, int* SecondQuantizedKet) {

	signed long long int ketNumberBin = *ketNumber;

	for (int i = 0; i < *NO; i++) {
		SecondQuantizedKet[*NO - (i + 1)] = ketNumberBin % 2;
		ketNumberBin = ketNumberBin / 2;
	}

}

void PrintKet(KET* prKet, int* NO) {

	for (int i = 0; i < (signed)prKet->KetElement.size(); i++) {
		int* SQKet;
		SQKet = new int[*NO];
		Secondize(&prKet->KetElement[i].base, NO, SQKet);
		if (i != 0) {
			std::cout << " + ";
			std::cout.flush();
		}
		std::cout << " ( " << prKet->KetElement[i].coeffcientRe << " , " << prKet->KetElement[i].coeffcientIm << " ) ";
		std::cout.flush();
		for (int j = 0; j < *NO; j++) {
			std::cout << SQKet[j];
			std::cout.flush();
		}
		delete[] SQKet;
	}

}

void simplify(KET* OutputKet, KET* SSimplifiedOutputKet) {

	SSimplifiedOutputKet->KetElement.clear();
	KET SimplifiedOutputKet;

	bool* ElementChecked;
	ElementChecked = new bool[(signed)OutputKet->KetElement.size()];
	for (int i = 0; i < (signed)OutputKet->KetElement.size(); i++) {
		ElementChecked[i] = false;
	}

	for (int i = 0; i < (signed)OutputKet->KetElement.size(); i++) {
		if (ElementChecked[i] == false) {
			KETELEMENT obj_ketelement;
			obj_ketelement.coeffcientRe = 0.0;
			obj_ketelement.coeffcientIm = 0.0;
			obj_ketelement.base = OutputKet->KetElement[i].base;
			for (int j = i; j < (signed)OutputKet->KetElement.size(); j++) {
				if (ElementChecked[j] == false) {
					if (OutputKet->KetElement[j].base == OutputKet->KetElement[i].base) {
						obj_ketelement.coeffcientRe += OutputKet->KetElement[j].coeffcientRe;
						obj_ketelement.coeffcientIm += OutputKet->KetElement[j].coeffcientIm;
						ElementChecked[j] = true;
					}
				}
			}
			SimplifiedOutputKet.KetElement.push_back(obj_ketelement);
		}
	}

	for (int i = 0; i < (signed)SimplifiedOutputKet.KetElement.size(); i++) {
		if ((ChopBool(&SimplifiedOutputKet.KetElement[i].coeffcientRe) && ChopBool(&SimplifiedOutputKet.KetElement[i].coeffcientIm)) == false) {
			SSimplifiedOutputKet->KetElement.push_back(SimplifiedOutputKet.KetElement[i]);
		}
	}

	delete[] ElementChecked;

}

void SumKets(long double* res, long double* ims, std::vector<KETS>* kets, int* HowMany, KET* kett) {

	KET objKET;
	objKET.KetElement.clear();
	
	for (int i = 0; i < *HowMany; i++) {
		for (int j = 0; j < (signed)(*kets)[i].ket.KetElement.size(); j++) {
			KETELEMENT objKETELEMENT;
			ComplexMTwoRe(&res[i], &ims[i], &(*kets)[i].ket.KetElement[j].coeffcientRe, &(*kets)[i].ket.KetElement[j].coeffcientIm, &objKETELEMENT.coeffcientRe);
			ComplexMTwoIm(&res[i], &ims[i], &(*kets)[i].ket.KetElement[j].coeffcientRe, &(*kets)[i].ket.KetElement[j].coeffcientIm, &objKETELEMENT.coeffcientIm);
			objKETELEMENT.base = (*kets)[i].ket.KetElement[j].base;
			objKET.KetElement.push_back(objKETELEMENT);
		}
	}

	simplify(&objKET, kett);

}

void KetTensorProduct(signed long long int* ketNumber1, int* NO1, signed long long int* ketNumber2, int* NO2, signed long long int* Bin) {

	int* SQKet1;
	SQKet1 = new int[*NO1];
	Secondize(ketNumber1, NO1, SQKet1);

	int* SQKet2;
	SQKet2 = new int[*NO2];
	Secondize(ketNumber2, NO2, SQKet2);

	int* SQKetT;
	int NOT = *NO1 + *NO2;
	SQKetT = new int[NOT];
	for (int i = 0; i < NOT; i++) {
		if (i < *NO1) {
			SQKetT[i] = SQKet1[i];
		}
		else {
			SQKetT[i] = SQKet2[i - *NO1];
		}
	}

	ConvertKetNumber(SQKetT, &NOT, Bin);

	delete[] SQKet1;
	delete[] SQKet2;
	delete[] SQKetT;

}

void KetsTensorProduct(KET* KET1, int* NO1, KET* KET2, int* NO2, KET* KETT) {

	KETT->KetElement.clear();

	for (int i = 0; i < KET1->KetElement.size(); i++) {
		for (int j = 0; j < KET2->KetElement.size(); j++) {
			KETELEMENT objKETELEMENT;
			KetTensorProduct(&KET1->KetElement[i].base, NO1, &KET2->KetElement[j].base, NO2, &objKETELEMENT.base);
			ComplexMTwoRe(&KET1->KetElement[i].coeffcientRe, &KET1->KetElement[i].coeffcientIm, &KET2->KetElement[j].coeffcientRe, &KET2->KetElement[j].coeffcientIm, &objKETELEMENT.coeffcientRe);
			ComplexMTwoIm(&KET1->KetElement[i].coeffcientRe, &KET1->KetElement[i].coeffcientIm, &KET2->KetElement[j].coeffcientRe, &KET2->KetElement[j].coeffcientIm, &objKETELEMENT.coeffcientIm);
			KETT->KetElement.push_back(objKETELEMENT);
		}
	}

}

void ConcatenateKetsTensorProduct(KET* KETS, int* NOs, int* HowMany, KET* KETT) {

	KETT->KetElement.clear();
	*KETT = KETS[0];
	int NOCount = NOs[0];

	for (int i = 1; i < *HowMany; i++) {
		KET BufferKET = *KETT;
		KetsTensorProduct(&BufferKET, &NOCount, &KETS[i], &NOs[i], KETT);
		NOCount += NOs[i];
	}

}

void KetsAndKetsTensorProduct(std::vector<KETS>* KETS1, int* NO1, std::vector<KETS>* KETS2, int* NO2, std::vector<KETS>* KETST) {

	KETST->clear();

	for (int i = 0; i < (signed)KETS1->size(); i++) {
		for (int j = 0; j < (signed)KETS2->size(); j++) {
			KET* objKETs;
			objKETs = new KET[2];
			objKETs[0] = (*KETS1)[i].ket;
			objKETs[1] = (*KETS2)[j].ket;
			int* NOs;
			NOs = new int[2];
			NOs[0] = *NO1;
			NOs[1] = *NO2;
			int HowMany = 2;
			KET objKET;
			ConcatenateKetsTensorProduct(objKETs, NOs, &HowMany, &objKET);
			KETS objKETS;
			objKETS.energy = (*KETS1)[i].energy + (*KETS2)[j].energy;
			objKETS.ket = objKET;
			KETST->push_back(objKETS);
			delete[] objKETs;
			delete[] NOs;
		}
	}

}

void ConcatenateKetsAndKetsTensorProduct(std::vector<KETSPTRS>* ksptrs, int* NOs, std::vector<KETS>* KETST) {

	KETST->clear();
	*KETST = *((*ksptrs)[0].ketsptrs);
	int NOCount = NOs[0];
	
	for (int i = 1; i < (signed)ksptrs->size(); i++) {
		std::vector<KETS> BufferKets = *KETST;
		KetsAndKetsTensorProduct(&BufferKets, &NOCount, (*ksptrs)[i].ketsptrs, &NOs[i], KETST);
		NOCount += NOs[i];
	}

}

void create_vac(KET* objKET) {

	objKET->KetElement.clear();
	KETELEMENT objKETELEMENT;
	objKETELEMENT.coeffcientRe = 1.0;
	objKETELEMENT.coeffcientIm = 0.0;
	objKETELEMENT.base = 0;
	objKET->KetElement.push_back(objKETELEMENT);

}

void create_full(int *NO, KET* objKET) {

	objKET->KetElement.clear();
	KETELEMENT objKETELEMENT;
	objKETELEMENT.coeffcientRe = 1.0;
	objKETELEMENT.coeffcientIm = 0.0;
	signed long long int Bin;
	int* SQKet;
	SQKet = new int[*NO];
	for (int i = 0; i < *NO; i++) {
		SQKet[i] = 1;
	}
	ConvertKetNumber(SQKet, NO, &Bin);
	objKETELEMENT.base = Bin;
	objKET->KetElement.push_back(objKETELEMENT);
	delete[] SQKet;

}

void operate(OPERATOR* _operator, KET* InputKet, int* NO, bool* simple, KET* SSimplifiedOutputKet) {

	SSimplifiedOutputKet->KetElement.clear();
	KET OutputKet;
	KET SimplifiedOutputKet;

	for (int i = 0; i < (signed)InputKet->KetElement.size(); i++) {

		int* SecondQuantizedKet;
		SecondQuantizedKet = new int[*NO];
		Secondize(&InputKet->KetElement[i].base, NO, SecondQuantizedKet);
		/// It might be not so significant to discuss, but for those who are confused...
		/// Definition goes like Array(SecondQuantizedKet/ConvertedQuantizedKet){[0], [1], .... , [NSite-1]} = Exists/Nope at {Site_(NSite-1), Site_(NSite-2), ... , Site_0}
		/// Exists = 1; Nope = 0 

		/// OperatorFiles also should be prepared like...
		/// Coeff. CR{Site_(NSite-1), Site_(NSite-2), ... , Site_0} AN{Site_(NSite-1), Site_(NSite-2), ... , Site_0}

		for (int j = 0; j < (signed)_operator->OperatorLine.size(); j++) {

			bool IfAn = true;
			for (int iSite = *NO - 1; iSite >= 0; iSite--) {
				if (_operator->OperatorLine[j].ANList[iSite] != -1 && _operator->OperatorLine[j].ANList[iSite] != SecondQuantizedKet[iSite]) {
					IfAn = false;
				}
				if (IfAn == false) {
					break;
				}
			}

			if (IfAn == true) {

				int MeaningfulAn = 0;
				for (int iSite = *NO - 1; iSite >= 0; iSite--) {
					if (_operator->OperatorLine[j].ANList[iSite] != -1) {
						MeaningfulAn++;
					}
				}

				std::vector<int> omitted;
				bool invalid_operate = false;

				for (int iSite = *NO - 1; iSite >= 0; iSite--) {
					if (_operator->OperatorLine[j].ANList[iSite] == -1 && SecondQuantizedKet[iSite] == 1) {
						if (_operator->OperatorLine[j].CRList[iSite] != SecondQuantizedKet[iSite]) {
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
							if (_operator->OperatorLine[j].ANList[iSite] != -1) {
								count1++;
							}
						}
					}

					for (int iOmitted = 0; iOmitted < (signed)omitted.size(); iOmitted++) {
						for (int iSite = 0; iSite <= omitted[iOmitted] - 1; iSite++) {
							if (_operator->OperatorLine[j].CRList[iSite] != -1) {
								count2++;
							}
						}
					}

					total_count = MeaningfulAn * (MeaningfulAn - 1) / 2 + count1 + count2;

					KETELEMENT obj_element;
					long double buf_re, buf_im;
					ComplexMTwoRe(&InputKet->KetElement[i].coeffcientRe, &InputKet->KetElement[i].coeffcientIm, &_operator->OperatorLine[j].coeffcientRe, &_operator->OperatorLine[j].coeffcientIm, &buf_re);
					ComplexMTwoIm(&InputKet->KetElement[i].coeffcientRe, &InputKet->KetElement[i].coeffcientIm, &_operator->OperatorLine[j].coeffcientRe, &_operator->OperatorLine[j].coeffcientIm, &buf_im);
					obj_element.coeffcientRe = pow(-1, total_count) * buf_re;
					obj_element.coeffcientIm = pow(-1, total_count) * buf_im;
					int* ConvertedQuantizedKet;
					ConvertedQuantizedKet = new int[*NO];
					for (int iSite = *NO - 1; iSite >= 0; iSite--) {
						ConvertedQuantizedKet[iSite] = 0;
						if (_operator->OperatorLine[j].CRList[iSite] != -1) {
							ConvertedQuantizedKet[iSite] = 1;
						}
						for (int iOmitted = 0; iOmitted < (signed)omitted.size(); iOmitted++) {
							if (omitted[iOmitted] == iSite) {
								ConvertedQuantizedKet[iSite] = 1;
								break;
							}
						}
					}
					ConvertKetNumber(ConvertedQuantizedKet, NO, &obj_element.base);
					OutputKet.KetElement.push_back(obj_element);
					delete[] ConvertedQuantizedKet;
				}

			}

		}

		delete[] SecondQuantizedKet;

	}

	if (*simple == true) {
		simplify(&OutputKet, SSimplifiedOutputKet);
	}
	else {
		*SSimplifiedOutputKet = OutputKet;
	}

}

void operate_scalar(KET* InputBra, KET* SSimplifiedOutputKet, long double* scalars) {

	scalars[0] = 0.0;
	scalars[1] = 0.0;

	for (int i = 0; i < (signed)SSimplifiedOutputKet->KetElement.size(); i++) {
		for (int j = 0; j < (signed)InputBra->KetElement.size(); j++) {
			if (SSimplifiedOutputKet->KetElement[i].base == InputBra->KetElement[j].base) {
				long double re_buf, im_buf;
				long double scndargbuf = -InputBra->KetElement[j].coeffcientIm;
				ComplexMTwoRe(&InputBra->KetElement[j].coeffcientRe, &scndargbuf, &SSimplifiedOutputKet->KetElement[i].coeffcientRe, &SSimplifiedOutputKet->KetElement[i].coeffcientIm, &re_buf);
				ComplexMTwoIm(&InputBra->KetElement[j].coeffcientRe, &scndargbuf, &SSimplifiedOutputKet->KetElement[i].coeffcientRe, &SSimplifiedOutputKet->KetElement[i].coeffcientIm, &im_buf);
				scalars[0] += re_buf;
				scalars[1] += im_buf;
			}
		}
	}

}

void operate_scalar_fullet(KET* InputBra, OPERATOR* _operator, KET* InputKet, int* NO, bool* simple, long double* scalars) {

	scalars[0] = 0.0;
	scalars[1] = 0.0;
	KET SSimplifiedOutputKet;
	operate(_operator, InputKet, NO, simple, &SSimplifiedOutputKet);
	operate_scalar(InputBra, &SSimplifiedOutputKet, scalars);

}

void read_operator(std::string* operator_filename, int* NO, OPERATOR* objOPERATOR) {

	std::ifstream readoperator;
	readoperator.open(*operator_filename);

	while (true) {

		OPERATORLINE objOPERATORLINE;
		getline(readoperator, sbuf);
		std::stringstream line_reader;
		line_reader << sbuf;
		line_reader >> objOPERATORLINE.coeffcientRe >> objOPERATORLINE.coeffcientIm;
		for (int i = 0; i < *NO; i++) {
			int ibuf;
			line_reader >> ibuf;
			objOPERATORLINE.CRList.push_back(ibuf);
		}
		for (int i = 0; i < *NO; i++) {
			int ibuf;
			line_reader >> ibuf;
			objOPERATORLINE.ANList.push_back(ibuf);
		}
		objOPERATOR->OperatorLine.push_back(objOPERATORLINE);

		if (readoperator.eof()) {
			break;
		}

	}

	readoperator.close();

}

void nCr(signed long long int* n, signed long long int* r, signed long long int* getnCr) {

	signed long long int nf, nrf, rf;

	nf = 1;
	for (signed long long int inf = 1; inf <= *n; inf++) {
		nf *= inf;
	}

	nrf = 1;
	for (signed long long int inrf = 1; inrf <= *n - *r; inrf++) {
		nrf *= inrf;
	}

	rf = 1;
	for (signed long long int irf = 1; irf <= *r; irf++) {
		rf *= irf;
	}

	*getnCr = nf / (nrf * rf);

}

void read_LSCoupling_SH(FILE *_readinput, signed long long int* Naimed, int* NO) {

	//index L Lz S Sz U Nbases ket(Re,Im,BinaryKet)
	for (signed long long int i = 0; i < *Naimed; i++) {
		signed long long int sllibuf;
		long double ldbuf;
		fread(&sllibuf, sizeof(signed long long int), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&sllibuf, sizeof(signed long long int), 1, _readinput);
		for (signed long long int ibases = 0; ibases < sllibuf; ibases++) {
			int* SecondQuantizedKet;
			SecondQuantizedKet = new int[*NO];
			fread(&ldbuf, sizeof(long double), 1, _readinput);
			fread(&ldbuf, sizeof(long double), 1, _readinput);
			fread(&SecondQuantizedKet[0], sizeof(int), *NO, _readinput);
			delete[] SecondQuantizedKet;
		}
	}

}

void read_LSCoupling_CH(FILE* _readinput, signed long long int* Naimed) {

	//index L Lz S Sz U Nbases ket(Re,Im,BinaryKet)
	for (signed long long int i = 0; i < *Naimed; i++) {
		signed long long int sllibuf;
		long double ldbuf;
		fread(&sllibuf, sizeof(signed long long int), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&sllibuf, sizeof(signed long long int), 1, _readinput);
		for (signed long long int ibases = 0; ibases < sllibuf; ibases++) {
			fread(&ldbuf, sizeof(long double), 1, _readinput);
			fread(&ldbuf, sizeof(long double), 1, _readinput);
			signed long long int testket;
			fread(&testket, sizeof(signed long long int), 1, _readinput);
		}
	}

}

void read_LSJCoupling_SH(FILE* _readinput, signed long long int* Naimed, int* NO) {

	//index L S J Jz E Nbases ket(Re,Im,BinaryKet)
	for (signed long long int i = 0; i < *Naimed; i++) {
		signed long long int sllibuf;
		long double ldbuf;
		fread(&sllibuf, sizeof(signed long long int), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&ldbuf, sizeof(long double), 1, _readinput);
		fread(&sllibuf, sizeof(signed long long int), 1, _readinput);
		for (signed long long int ibases = 0; ibases < sllibuf; ibases++) {
			int* SecondQuantizedKet;
			SecondQuantizedKet = new int[*NO];
			fread(&ldbuf, sizeof(long double), 1, _readinput);
			fread(&ldbuf, sizeof(long double), 1, _readinput);
			fread(&SecondQuantizedKet[0], sizeof(int), *NO, _readinput);
			delete[] SecondQuantizedKet;
		}
	}

}

void read_LSJCoupling_SH_and_save(FILE* _readinput, signed long long int* Naimed, int* NO, bool* condition, long double* L, long double* S, long double* J, long double* Jz, int* istart, int* iend, long double* energy_offset, std::vector<KETS>* kets) {

	//index L S J Jz E Nbases ket(Re,Im,BinaryKet)
	for (signed long long int i = 0; i < *Naimed; i++) {
		KETS objKETS;
		signed long long int sllibuf;
		long double vL, vS, vJ, vJz, vE, ldbuf;
		fread(&sllibuf, sizeof(signed long long int), 1, _readinput);
		fread(&vL, sizeof(long double), 1, _readinput); // L
		fread(&vS, sizeof(long double), 1, _readinput); // S
		fread(&vJ, sizeof(long double), 1, _readinput); // J
		fread(&vJz, sizeof(long double), 1, _readinput); // Jz
		vL -= *L;
		vS -= *S;
		vJ -= *J;
		vJz -= *Jz;
		if ((*condition == false && ((*istart == -1 || *iend == -1) || (i >= *istart && i <= *iend))) || (*condition == true && ChopBool(&vL) && ChopBool(&vS) && ChopBool(&vJ) && ChopBool(&vJz))) {
			fread(&vE, sizeof(long double), 1, _readinput); // E
			objKETS.energy = vE + *energy_offset;
			fread(&sllibuf, sizeof(signed long long int), 1, _readinput);
			for (signed long long int ibases = 0; ibases < sllibuf; ibases++) {
				KETELEMENT objKETELEMENT;
				fread(&ldbuf, sizeof(long double), 1, _readinput);
				objKETELEMENT.coeffcientRe = ldbuf;
				fread(&ldbuf, sizeof(long double), 1, _readinput);
				objKETELEMENT.coeffcientIm = ldbuf;
				int* SecondQuantizedKet;
				SecondQuantizedKet = new int[*NO];
				fread(&SecondQuantizedKet[0], sizeof(int), *NO, _readinput);
				signed long long int testket;
				ConvertKetNumber(SecondQuantizedKet, NO, &testket);
				objKETELEMENT.base = testket;
				objKETS.ket.KetElement.push_back(objKETELEMENT);
				delete[] SecondQuantizedKet;
			}
			kets->push_back(objKETS);
		}
		else {
			fread(&vE, sizeof(long double), 1, _readinput); // E
			fread(&sllibuf, sizeof(signed long long int), 1, _readinput);
			for (signed long long int ibases = 0; ibases < sllibuf; ibases++) {
				int* SecondQuantizedKet;
				SecondQuantizedKet = new int[*NO];
				fread(&ldbuf, sizeof(long double), 1, _readinput);
				fread(&ldbuf, sizeof(long double), 1, _readinput);
				fread(&SecondQuantizedKet[0], sizeof(int), *NO, _readinput);
				delete[] SecondQuantizedKet;
			}
		}

	}

}

// kets are not initialized INTENTIONALLY!
void read_LSJCoupling_CH_and_save(FILE* _readinput, signed long long int* Naimed, bool* condition, long double* L, long double* S, long double* J, long double *Jz, int* istart, int* iend, long double* energy_offset, std::vector<KETS>* kets) {

	//index L S J Jz E Nbases ket(Re,Im,BinaryKet)
	for (signed long long int i = 0; i < *Naimed; i++) {
		KETS objKETS;
		signed long long int sllibuf;
		long double vL, vS, vJ, vJz, vE, ldbuf;
		fread(&sllibuf, sizeof(signed long long int), 1, _readinput);
		fread(&vL, sizeof(long double), 1, _readinput); // L
		fread(&vS, sizeof(long double), 1, _readinput); // S
		fread(&vJ, sizeof(long double), 1, _readinput); // J
		fread(&vJz, sizeof(long double), 1, _readinput); // Jz
		vL -= *L;
		vS -= *S;
		vJ -= *J;
		vJz -= *Jz;
		if ((*condition == false && ((*istart == -1 || *iend == -1) || (i >= *istart && i <= *iend))) || (*condition == true && ChopBool(&vL) && ChopBool(&vS) && ChopBool(&vJ) && ChopBool(&vJz))) {
			fread(&vE, sizeof(long double), 1, _readinput); // E
			objKETS.energy = vE + *energy_offset;
			fread(&sllibuf, sizeof(signed long long int), 1, _readinput);
			for (signed long long int ibases = 0; ibases < sllibuf; ibases++) {
				KETELEMENT objKETELEMENT;
				fread(&ldbuf, sizeof(long double), 1, _readinput);
				objKETELEMENT.coeffcientRe = ldbuf;
				fread(&ldbuf, sizeof(long double), 1, _readinput);
				objKETELEMENT.coeffcientIm = ldbuf;
				signed long long int testket;
				fread(&testket, sizeof(signed long long int), 1, _readinput);
				objKETELEMENT.base = testket;
				objKETS.ket.KetElement.push_back(objKETELEMENT);
			}
			kets->push_back(objKETS);
		}
		else {
			fread(&vE, sizeof(long double), 1, _readinput); // E
			fread(&sllibuf, sizeof(signed long long int), 1, _readinput);
			for (signed long long int ibases = 0; ibases < sllibuf; ibases++) {
				fread(&ldbuf, sizeof(long double), 1, _readinput);
				fread(&ldbuf, sizeof(long double), 1, _readinput);
				signed long long int testket;
				fread(&testket, sizeof(signed long long int), 1, _readinput);
			}
		}
	}

}

// kets are not initialized INTENTIONALLY!
void read_kets(std::string* filename, signed long long int* Naimed, int *NO, bool* condition, long double* L, long double* S, long double* J, long double* Jz, int* istart, int* iend, long double* energy_offset, bool* SHtCHf, std::vector<KETS>* kets) {

	FILE* _readinput = nullptr;
	const char* ccfilename = filename->c_str();
#ifdef _WIN32
	errno_t error;
	error = fopen_s(&_readinput, ccfilename, "rb");
#elif _WIN64
	errno_t error;
	error = fopen_s(&_readinput, ccfilename, "rb");
#else
	_readinput = fopen(ccfilename, "rb");
#endif

	read_LSCoupling_SH(_readinput, Naimed, NO);
	read_LSCoupling_CH(_readinput, Naimed);
	if (*SHtCHf == true) {
		read_LSJCoupling_SH_and_save(_readinput, Naimed, NO, condition, L, S, J, Jz, istart, iend, energy_offset, kets);
	}
	else {
		read_LSJCoupling_SH(_readinput, Naimed, NO);
		read_LSJCoupling_CH_and_save(_readinput, Naimed, condition, L, S, J, Jz, istart, iend, energy_offset, kets);
	}
	
	fclose(_readinput);

}

// kets are not initialized INTENTIONALLY!
void read_all_kets(std::string* filename, signed long long int* Naimed, int* NO, long double* energy_offset, bool* SHtCHf, std::vector<KETS>* kets) {

	bool condition = false;
	long double L = 0.0;
	long double S = 0.0;
	long double J = 0.0;
	long double Jz = 0.0;
	int istart = -1;
	int iend = -1;
	read_kets(filename, Naimed, NO, &condition, &L, &S, &J, &Jz, &istart, &iend, energy_offset, SHtCHf, kets);

}

// kets are not initialized INTENTIONALLY!
void read_some_kets(std::string* filename, signed long long int* Naimed, int* NO, int* istart, int* iend, long double* energy_offset, bool* SHtCHf, std::vector<KETS>* kets) {

	bool condition = false;
	long double L = 0.0;
	long double S = 0.0;
	long double J = 0.0;
	long double Jz = 0.0;
	read_kets(filename, Naimed, NO, &condition, &L, &S, &J, &Jz, istart, iend, energy_offset, SHtCHf, kets);

}

// kets are not initialized INTENTIONALLY!
void read_conditional_kets(std::string* filename, signed long long int* Naimed, int* NO, long double* L, long double* S, long double* J, long double* Jz, long double* energy_offset, bool* SHtCHf, std::vector<KETS>* kets) {

	bool condition = true;
	int istart = -1;
	int iend = -1;
	read_kets(filename, Naimed, NO, &condition, L, S, J, Jz, &istart, &iend, energy_offset, SHtCHf, kets);

}

void read_ketsfile(bool initial_or_final, bool add_up, bool save_memory, bool SHtCHf) {

	readpara >> sbuf;
	int NStateFile;
	readpara >> sbuf;

	if (sbuf != "RE") {

		readpara >> NStateFile;
		if (initial_or_final == true) {
			std::string vHtFilenmae;
			readpara >> sbuf >> vHtFilenmae;
			if (HoppintHamiltonianFilename == "VA") {
				OPERATOR objHt;
				read_operator(&vHtFilenmae, &NOrbitalRooms, &objHt);
				Ht.push_back(objHt);
			}
			readpara >> sbuf;
		}

		KETSSS objKETSSS;
		for (int i = 0; i < NStateFile; i++) {
			KETSS objKETSS;
			objKETSSS.ketss.push_back(objKETSS);
		}
		objKETSSS.NOs = new int[NStateFile];

		for (int i = 0; i < NStateFile; i++) {
			std::string filename;
			int NOrbitalRoom, NOccupied, ReadStart, ReadEnd;
			long double energy_offset;
			readpara >> filename >> NOrbitalRoom >> NOccupied;
			if (initial_or_final == false) {
				readpara >> ReadStart >> ReadEnd >> energy_offset;
			}
			objKETSSS.NOs[i] = NOrbitalRoom;
			signed long long int slliNOrbitalRoom = NOrbitalRoom;
			signed long long int slliNOccupied = NOccupied;
			signed long long int NAimed;
			nCr(&slliNOrbitalRoom, &slliNOccupied, &NAimed);
			if (initial_or_final == true) {
				readpara >> sbuf;
				int NPickupState;
				readpara >> sbuf >> NPickupState;
				std::vector<KETS> objVKETS;
				for (int j = 0; j < NPickupState; j++) {
					int NKets;
					readpara >> sbuf >> NKets;
					long double* Res;
					Res = new long double[NKets];
					long double* Ims;
					Ims = new long double[NKets];
					std::vector<KETS> KetsElement;
					for (int k = 0; k < NKets; k++) {
						long double vL, vS, vJ, vJz;
						readpara >> sbuf >> Res[k] >> sbuf >> Ims[k] >> sbuf >> vL >> sbuf >> vS >> sbuf >> vJ >> sbuf >> vJz;
						long double dummy_energy = 0.0;
						read_conditional_kets(&filename, &NAimed, &NOrbitalRoom, &vL, &vS, &vJ, &vJz, &dummy_energy, &SHtCHf, &KetsElement);
					}
					KET PickedKet;
					SumKets(Res, Ims, &KetsElement, &NKets, &PickedKet);
					KETS objKETS;
					objKETS.ket = PickedKet;
					objKETS.energy = 0.0;
					objVKETS.push_back(objKETS);
				}
				objKETSSS.ketss[i].kets = objVKETS;
				readpara >> sbuf;
			}
			else {
				if (!(filename == "VAC" || filename == "FULL")) {
					if (ReadStart == -1 || ReadEnd == -1) {
						read_all_kets(&filename, &NAimed, &NOrbitalRoom, &energy_offset, &SHtCHf, &objKETSSS.ketss[i].kets);
					}
					else {
						read_some_kets(&filename, &NAimed, &NOrbitalRoom, &ReadStart, &ReadEnd, &energy_offset, &SHtCHf, &objKETSSS.ketss[i].kets);
					}
				}
				else if(filename == "VAC") {
					KET objvacKET;
					create_vac(&objvacKET);
					KETS objvacKETS;
					objvacKETS.energy = energy_offset;
					objvacKETS.ket = objvacKET;
					std::vector<KETS> objvacVKETS;
					objvacVKETS.push_back(objvacKETS);
					objKETSSS.ketss[i].kets = objvacVKETS;
				}
				else {
					KET objfullKET;
					create_full(&NOrbitalRoom, &objfullKET);
					KETS objfullKETS;
					objfullKETS.energy = energy_offset;
					objfullKETS.ket = objfullKET;
					std::vector<KETS> objfullVKETS;
					objfullVKETS.push_back(objfullKETS);
					objKETSSS.ketss[i].kets = objfullVKETS;
				}
			}
		}

		if (save_memory == true) {
			if (add_up == false) {
				KETSERIES objKETSERIES;
				objKETSERIES.ketsss.push_back(objKETSSS);
				ketseries.push_back(objKETSERIES);
			}
			else {
				ketseries[(signed)ketseries.size() - 1].ketsss.push_back(objKETSSS);
			}
		}
		else {
			std::vector<KETSPTRS> vobjKETSPTRS;
			for (int i = 0; i < NStateFile; i++) {
				KETSPTRS objKETSPTRS;
				objKETSPTRS.ketsptrs = &objKETSSS.ketss[i].kets;
				vobjKETSPTRS.push_back(objKETSPTRS);
			}
			std::vector<KETS> vKETS;
			ConcatenateKetsAndKetsTensorProduct(&vobjKETSPTRS, objKETSSS.NOs, &vKETS);
			if (add_up == false) {
				KETSSS objKETSSSf;
				KETSS objKETSSf;
				objKETSSf.kets = vKETS;
				objKETSSSf.ketss.push_back(objKETSSf);
				objKETSSSf.NOs = new int[1];
				objKETSSSf.NOs[0] = 0;
				for (int ifile = 0; ifile < NStateFile; ifile++) {
					objKETSSSf.NOs[0] += objKETSSS.NOs[ifile];
				}
				KETSERIES objKETSERIES;
				objKETSERIES.ketsss.push_back(objKETSSSf);
				ketseries.push_back(objKETSERIES);
			}
			else {
				for (int i = 0; i < (signed)vKETS.size(); i++) {
					ketseries[(signed)ketseries.size() - 1].ketsss[0].ketss[0].kets.push_back(vKETS[i]);
				}
			}
		}

		if (initial_or_final == true) {
			readpara >> sbuf;
		}

	}
	else {
	
		ketseries.push_back(ketseries[0]);

	}
	
}

void start_write_prtbterm(std::string* filename) {

	const char* ccfilename = filename->c_str();
#ifdef _WIN32
	errno_t error;
	error = fopen_s(&write_prtb, ccfilename, "wb");
#elif _WIN64
	errno_t error;
	error = fopen_s(&write_prtb, ccfilename, "wb");
#else
	write_prtb = fopen(ccfilename, "wb");
#endif

}

void throw_prtbterms(std::vector<PTRB>* thrownptrb) {

	for (int i = 0; i < (signed)thrownptrb->size(); i++) {
		fwrite((*thrownptrb)[i].v, sizeof(long double), 2, write_prtb);
	}
	fflush(write_prtb);
	std::vector<PTRB>().swap(*thrownptrb);

}

void write_prtbterms(long double* v) {

	fwrite(v, sizeof(long double), 2, write_prtb);
	fflush(write_prtb);

}

void close_write_prtbterm() {

	if (write_prtb != NULL) {
		fclose(write_prtb);
	}

}

void start_read_prtbterm(std::string* filename) {

	const char* ccfilename = filename->c_str();
#ifdef _WIN32
	errno_t error;
	error = fopen_s(&read_prtb, ccfilename, "rb");
#elif _WIN64
	errno_t error;
	error = fopen_s(&read_prtb, ccfilename, "rb");
#else
	read_prtb = fopen(ccfilename, "rb");
#endif

	if (!read_prtb) {
		ReadPrtb = 0;
	}

}

void read_prtbterm(long double *v) {

	if (fread(v, sizeof(long double), 2, read_prtb) != 2) {
		ReadPrtb = 0;
		fclose(read_prtb);
	}

}

void get_prtbterm_simply(int* braindex, int* ketindex, OPERATOR* Ham, std::vector<PTRBMATRIX>* ptrbmatrix) {

	ptrbmatrix->clear();

	PTRB objPTRB;
	objPTRB.v = new long double[2];
	for (int i = 0; i < (signed)ketseries[*ketindex].ketsss[0].ketss[0].kets.size(); i++) {
		PTRBMATRIX objPTRBMATRIX;
		ptrbmatrix->push_back(objPTRBMATRIX);
		for (int j = 0; j < (signed)ketseries[*braindex].ketsss[0].ketss[0].kets.size(); j++) {
			bool simple = false;
			bool readfromdat = true;
			if (ReadPrtb == 0) {
				operate_scalar_fullet(&ketseries[*braindex].ketsss[0].ketss[0].kets[j].ket, Ham, &ketseries[*ketindex].ketsss[0].ketss[0].kets[i].ket, &NOrbitalRooms, &simple, objPTRB.v);
				readfromdat = false;
			}
			else {
				read_prtbterm(objPTRB.v);
				if (ReadPrtb == 0) {
					operate_scalar_fullet(&ketseries[*braindex].ketsss[0].ketss[0].kets[j].ket, Ham, &ketseries[*ketindex].ketsss[0].ketss[0].kets[i].ket, &NOrbitalRooms, &simple, objPTRB.v);
					readfromdat = false;
				}
			}
			if (WritePrtb == 1) {
				ptrbbuf.push_back(objPTRB);
			}
			if (ReadPrtb == 0 && WritePrtb == 1) {
				if (just_thrown == true) {
					std::string filename_writeprtb = "__" + std::to_string(pr_id) + "_" + OutputFilename + ".dat";
					start_write_prtbterm(&filename_writeprtb);
					throw_prtbterms(&ptrbbuf);
					just_thrown = false;
				}
				else {
					write_prtbterms(objPTRB.v);
				}
			}
			if (*braindex != ketseries.size() - 1) {
				if (!ChopBool(&ketseries[*braindex].ketsss[0].ketss[0].kets[j].energy)) {
					for (int k = 0; k < 2; k++) {
						objPTRB.v[k] /= -ketseries[*braindex].ketsss[0].ketss[0].kets[j].energy;
					}
				}
				else {
					std::cout << "[ERROR] at " << *braindex << " " << j << " ";
					std::cout.flush();
					PrintKet(&ketseries[*braindex].ketsss[0].ketss[0].kets[j].ket, &ketseries[*ketindex].ketsss[0].NOs[0]);
					std::cout << " , ENERGY = 0 !\n";
					std::cout.flush();
				}
			}
			PTRB replicaPTRB;
			replicaPTRB.v = new long double[2];
			replicaPTRB.v[0] = objPTRB.v[0];
			replicaPTRB.v[1] = objPTRB.v[1];
			(*ptrbmatrix)[(signed)ptrbmatrix->size() - 1].ptrb.push_back(replicaPTRB);
			prtb_count++;
			if (pr_id == 0 && (prtb_count == 1 || prtb_count % MonitoringIndex == 0)) {
				auto endefr = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
				writeoutput << "<\t" << *braindex + 1 << "\t|...|\t" << *ketindex + 1 << "\t>\tprtb_count=\t" << prtb_count << "\tT=\t" << elapsed_secondsfr.count();
				writeoutput.flush();
				if (FullMonitoring == 2) {
					writeoutput << "\t" << replicaPTRB.v[0] << "\t" << replicaPTRB.v[1];
					writeoutput.flush();
				}
				if (readfromdat == true) {
					writeoutput << "\tfromdat";
					writeoutput.flush();
				}
				writeoutput << "\n";
				writeoutput.flush();
			}
			if (pr_id != 0 && (FullMonitoring == 1 || FullMonitoring == 2)) {
				auto endefr = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
				writeoutput_sub << "<\t" << *braindex + 1 << "\t|...|\t" << *ketindex + 1 << "\t>\tprtb_count=\t" << prtb_count << "\tT=\t" << elapsed_secondsfr.count();
				writeoutput_sub.flush();
				if (FullMonitoring == 2) {
					writeoutput_sub << "\t" << replicaPTRB.v[0] << "\t" << replicaPTRB.v[1];
					writeoutput_sub.flush();
				}
				if (readfromdat == true) {
					writeoutput_sub << "\tfromdat";
					writeoutput_sub.flush();
				}
				writeoutput_sub << "\n";
				writeoutput_sub.flush();
			}
		}
	}
	
}

// ptrbmatrix is not initialized INTENTIONALLY!
// check ketcomb and bracomb before using this!
void get_prtbterm(int* braindex, int* ketindex, int* braiketsss, int* ketiketsss, OPERATOR* Ham, std::vector<PTRBMATRIX>* ptrbmatrix) {

	KET* ketready;
	ketready = new KET[(signed)ketcomb.size()];
	for (int i = 0; i < (signed)ketcomb.size(); i++) {
		ketready[i] = ketseries[*ketindex].ketsss[*ketiketsss].ketss[i].kets[ketcomb[i]].ket;
	}

	KET ketKETT;
	int HowManyKetReady = (signed)ketcomb.size();
	ConcatenateKetsTensorProduct(ketready, ketseries[*ketindex].ketsss[*ketiketsss].NOs, &HowManyKetReady, &ketKETT);

	std::string ketlabelstring;
	std::stringstream ketstream;
	if (GetLabel == 1) {
		for (int i = 0; i < (signed)ketcomb.size(); i++) {
			ketstream << *ketindex << "-" << *ketiketsss << "-" << i << "-" << ketseries[*ketindex].ketsss[*ketiketsss].ketss[i].kets[ketcomb[i]].label;
			if (i != (signed)ketcomb.size() - 1) {
				ketstream << "^";
			}
		}
		std::getline(ketstream, ketlabelstring);
		if (*ketindex == 0 && globalketstring.size() < globalketstring_max) {
			bool addglobalketlabel = true;
			for (int ilabel = 0; ilabel < (signed)globalketstring.size(); ilabel++) {
				if (globalketstring[ilabel] == ketlabelstring) {
					addglobalketlabel = false;
					break;
				}
			}
			if (addglobalketlabel == true) {
				globalketstring.push_back(ketlabelstring);
			}
		}
		if (*ketindex != 0) {
			bool addintermediateketlabel = true;
			for (int ilabel = 0; ilabel < (signed)ketseries[*ketindex].labels.size(); ilabel++) {
				if (ketseries[*ketindex].labels[ilabel] == ketlabelstring) {
					addintermediateketlabel = false;
					break;
				}
			}
			if (addintermediateketlabel == true) {
				ketseries[*ketindex].labels.push_back(ketlabelstring);
			}
		}
	}

	KET* braready;
	long double energy = 0.0;
	braready = new KET[(signed)bracomb.size()];
	for (int i = 0; i < (signed)bracomb.size(); i++) {
		braready[i] = ketseries[*braindex].ketsss[*braiketsss].ketss[i].kets[bracomb[i]].ket;
		energy += ketseries[*braindex].ketsss[*braiketsss].ketss[i].kets[bracomb[i]].energy;
	}
	
	KET braKETT;
	int HowManyBraReady = (signed)bracomb.size();
	ConcatenateKetsTensorProduct(braready, ketseries[*braindex].ketsss[*braiketsss].NOs, &HowManyBraReady, &braKETT);

	std::string bralabelstring;
	std::stringstream brastream;
	if (GetLabel == 1) {
		for (int i = 0; i < (signed)bracomb.size(); i++) {
			brastream << *braindex << "-" << *braiketsss << "-" << i << "-" << ketseries[*braindex].ketsss[*braiketsss].ketss[i].kets[bracomb[i]].label;
			if (i != (signed)bracomb.size() - 1) {
				brastream << "^";
			}
		}
		std::getline(brastream, bralabelstring);
		if (*braindex == (signed)ketseries.size() - 1 && globalbrastring.size() < globalbrastring_max) {
			bool addglobalbralabel = true;
			for (int ilabel = 0; ilabel < (signed)globalbrastring.size(); ilabel++) {
				if (globalbrastring[ilabel] == bralabelstring) {
					addglobalbralabel = false;
					break;
				}
			}
			if (addglobalbralabel == true) {
				globalbrastring.push_back(bralabelstring);
			}
		}
	}

	PTRB objPTRB;
	objPTRB.v = new long double[2];
	bool simple = false;
	bool readfromdat = true;
	if (ReadPrtb == 0) {
		operate_scalar_fullet(&braKETT, Ham, &ketKETT, &NOrbitalRooms, &simple, objPTRB.v);
		readfromdat = false;
	}
	else {
		read_prtbterm(objPTRB.v);
		if (ReadPrtb == 0) {
			operate_scalar_fullet(&braKETT, Ham, &ketKETT, &NOrbitalRooms, &simple, objPTRB.v);
			readfromdat = false;
		}
	}
	if (WritePrtb == 1) {
		ptrbbuf.push_back(objPTRB);
	}
	if (ReadPrtb == 0 && WritePrtb == 1) {
		if (just_thrown == true) {
			std::string filename_writeprtb = "__" + std::to_string(pr_id) + "_" + OutputFilename + ".dat";
			start_write_prtbterm(&filename_writeprtb);
			throw_prtbterms(&ptrbbuf);
			just_thrown = false;
		}
		else {
			write_prtbterms(objPTRB.v);
		}
	}
	if (*braindex != ketseries.size() - 1) {
		for (int i = 0; i < 2; i++) {
			if (!ChopBool(&energy)) {
				objPTRB.v[i] /= -energy;
			}
			else {
				std::cout << "[ERROR] at " << *braindex << " " << braiketsss << " ";
				std::cout.flush();
				for (int j = 0; j < (signed)bracomb.size(); j++) {
					std::cout << bracomb[j];
					std::cout.flush();
				}
				PrintKet(&braKETT, &NOrbitalRooms);
				std::cout << " , ENERGY = 0 !\n";
				std::cout.flush();
			}
			
		}
	}

	PTRB replicaPTRB;
	replicaPTRB.v = new long double[2];
	replicaPTRB.v[0] = objPTRB.v[0];
	replicaPTRB.v[1] = objPTRB.v[1];
	if (ptrbmatrixaddup == true) {
		PTRBMATRIX objPTRBMATRIX;
		objPTRBMATRIX.ptrb.push_back(replicaPTRB);
		ptrbmatrix->push_back(objPTRBMATRIX);
	}
	else {
		int get_size = (signed)ptrbmatrix->size();
		(*ptrbmatrix)[get_size - 1].ptrb.push_back(replicaPTRB);
	}

	delete[] ketready;
	delete[] braready;

	prtb_count++;
	if (pr_id == 0 && (prtb_count == 1 || prtb_count % MonitoringIndex == 0)) {
		auto endefr = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
		if (GetLabel == 0) {
			writeoutput << "<\t" << *braindex + 1 << "\t|...|\t" << *ketindex + 1 << "\t>\tprtb_count=\t" << prtb_count << "\tT=\t" << elapsed_secondsfr.count();
			writeoutput.flush();
		}
		else {
			writeoutput << "<\t" << bralabelstring << "\t|...|\t" << ketlabelstring << "\t>\tprtb_count=\t" << prtb_count << "\tT=\t" << elapsed_secondsfr.count();
			writeoutput.flush();
		}
		if (FullMonitoring == 2) {
			writeoutput << "\t" << replicaPTRB.v[0] << "\t" << replicaPTRB.v[1];
			writeoutput.flush();
		}
		if (readfromdat == true) {
			writeoutput << "\tfromdat";
			writeoutput.flush();
		}
		writeoutput << "\n";
		writeoutput.flush();
	}
	if (pr_id != 0 && (FullMonitoring == 1 || FullMonitoring == 2)) {
		auto endefr = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
		if (GetLabel == 0) {
			writeoutput_sub << "<\t" << *braindex + 1 << "\t|...|\t" << *ketindex + 1 << "\t>\tprtb_count=\t" << prtb_count << "\tT=\t" << elapsed_secondsfr.count();
			writeoutput_sub.flush();
		}
		else {
			writeoutput_sub << "<\t" << bralabelstring << "\t|...|\t" << ketlabelstring << "\t>\tprtb_count=\t" << prtb_count << "\tT=\t" << elapsed_secondsfr.count();
			writeoutput_sub.flush();
		}
		if (FullMonitoring == 2) {
			writeoutput_sub << "\t" << replicaPTRB.v[0] << "\t" << replicaPTRB.v[1];
			writeoutput_sub.flush();
		}
		if (readfromdat == true) {
			writeoutput_sub << "\tfromdat";
			writeoutput_sub.flush();
		}
		writeoutput_sub << "\n";
		writeoutput_sub.flush();
	}

}

// ptrbmatrix is not initialized INTENTIONALLY!
void rc_get_bracomb(int* braindex, int* ketindex, int* braiketsss, int* ketiketsss, int depth, OPERATOR* Ham, std::vector<PTRBMATRIX>* ptrbmatrix) {

	for (signed long long int i = 0; i < (signed)ketseries[*braindex].ketsss[*braiketsss].ketss[depth].kets.size(); i++) {
		bracomb[depth] = i;
		if (depth != bracomb.size() - 1) {
			rc_get_bracomb(braindex, ketindex, braiketsss, ketiketsss, depth + 1, Ham, ptrbmatrix);
		}
		else {
			get_prtbterm(braindex, ketindex, braiketsss, ketiketsss, Ham, ptrbmatrix);
			ptrbmatrixaddup = false;
		}
	}

}

// ptrbmatrix is not initialized INTENTIONALLY!
void rc_get_ketcomb(int* braindex, int* ketindex, int* braiketsss, int* ketiketsss, int depth, OPERATOR* Ham, std::vector<PTRBMATRIX>* ptrbmatrix) {

	for (signed long long int i = 0; i < (signed)ketseries[*ketindex].ketsss[*ketiketsss].ketss[depth].kets.size(); i++) {
		ketcomb[depth] = i;
		if (depth != ketcomb.size() - 1) {
			rc_get_ketcomb(braindex, ketindex, braiketsss, ketiketsss, depth + 1, Ham, ptrbmatrix);
		}
		else {
			ptrbmatrixaddup = true;
			for (int j = 0; j < (signed)ketseries[*braindex].ketsss.size(); j++) {
				bracomb.clear();
				for (int k = 0; k < (signed)ketseries[*braindex].ketsss[j].ketss.size(); k++) {
					bracomb.push_back(-1);
				}
				rc_get_bracomb(braindex, ketindex, &j, ketiketsss, 0, Ham, ptrbmatrix);
			}
		}
	}

}

void get_perturbation(int* braindex, int* ketindex, OPERATOR* Ham, std::vector<PTRBMATRIX>* ptrbmatrix) {

	ptrbmatrix->clear();

	for (int i = 0; i < (signed)ketseries[*ketindex].ketsss.size(); i++) {
		ketcomb.clear();
		bracomb.clear();
		for (int j = 0; j < (signed)ketseries[*ketindex].ketsss[i].ketss.size(); j++) {
			ketcomb.push_back(-1);
		}
		int braiketsssdummy = -1;
		rc_get_ketcomb(braindex, ketindex, &braiketsssdummy, &i, 0, Ham, ptrbmatrix);
	}

}

void read_para() {

	readpara.open("basic_properties_SQPerturbation.txt");

	readpara >> sbuf >> OutputFilename;
	readpara >> sbuf >> OutputPrecision;
	if (pr_id == 0) {
		writeoutput.open(OutputFilename + ".txt");
		writeoutput << std::fixed << std::setprecision(OutputPrecision);
	}

	readpara >> sbuf >> ChopLevel;
	Chop = pow(10.0, (long double)(-ChopLevel));

	readpara >> sbuf >> ReadPrtb;
	readpara >> sbuf >> WritePrtb;
	readpara >> sbuf >> NOrbitalRooms;
	readpara >> sbuf >> HoppintHamiltonianFilename;
	if (HoppintHamiltonianFilename != "VA") {
		OPERATOR objHt;
		read_operator(&HoppintHamiltonianFilename, &NOrbitalRooms, &objHt);
		Ht.push_back(objHt);
	}

	readpara >> sbuf >> ManualMPI;
	readpara >> sbuf >> ManualMPIFilename;
	readpara >> sbuf >> SaveMemory;
	readpara >> sbuf >> GetLabel;
	readpara >> sbuf >> ReadSHorCH;
	bool SHtCHf;
	if (ReadSHorCH == 1) {
		SHtCHf = true;
	}
	else {
		SHtCHf = false;
	}

	readpara >> sbuf >> MonitoringIndex;
	readpara >> sbuf >> FullMonitoring;
	if (pr_id != 0 && (FullMonitoring == 1 || FullMonitoring == 2)) {
		writeoutput_sub.open("_" + std::to_string(pr_id) + "_" + OutputFilename + ".txt");
		writeoutput_sub << std::fixed << std::setprecision(OutputPrecision);
	}

	if (SaveMemory == 0) {
		read_ketsfile(true, false, false, SHtCHf);
	}
	else {
		read_ketsfile(true, false, true, SHtCHf);
	}
	if (pr_id == 0) {
		auto endefr = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
		writeoutput << "initial file(s) read...\t" << elapsed_secondsfr.count() << "\n";
		writeoutput.flush();
	}
	if (pr_id != 0 && (FullMonitoring == 1 || FullMonitoring == 2)) {
		auto endefr = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
		writeoutput_sub << "initial file(s) read...\t" << elapsed_secondsfr.count() << "\n";
		writeoutput_sub.flush();
	}
	
	readpara >> sbuf;
	int NIntermediates;
	readpara >> sbuf >> NIntermediates;
	for (int i = 0; i < NIntermediates; i++) {

		readpara >> sbuf;
		int NIntermediateStateGroup;
		readpara >> sbuf >> NIntermediateStateGroup;

		std::string vHtFilenmae;
		readpara >> sbuf >> vHtFilenmae;
		if (HoppintHamiltonianFilename != "VA") {
			Ht.push_back(Ht[Ht.size() - 1]);
		}
		else {
			OPERATOR objHt;
			read_operator(&vHtFilenmae, &NOrbitalRooms, &objHt);
			Ht.push_back(objHt);
		}

		for (int j = 0; j < NIntermediateStateGroup; j++) {
			if (j == 0) {
				if (SaveMemory == 0) {
					read_ketsfile(false, false, false, SHtCHf);
				}
				else {
					read_ketsfile(false, false, true, SHtCHf);
				}
			}
			else {
				if (SaveMemory == 0) {
					read_ketsfile(false, true, false, SHtCHf);
				}
				else {
					read_ketsfile(false, true, true, SHtCHf);
				}
			}
		}
		
		if (pr_id == 0) {
			auto endefr = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
			writeoutput << "intermediate file(s) #" << i + 1 << " read...\t" << elapsed_secondsfr.count() << "\n";
			writeoutput.flush();
		}
		if (pr_id != 0 && (FullMonitoring == 1 || FullMonitoring == 2)) {
			auto endefr = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
			writeoutput_sub << "intermediate file(s) #" << i + 1 << " read...\t" << elapsed_secondsfr.count() << "\n";
			writeoutput_sub.flush();
		}

	}
	readpara >> sbuf;

	if (SaveMemory == 0) {
		read_ketsfile(true, false, false, SHtCHf);
	}
	else {
		read_ketsfile(true, false, true, SHtCHf);
	}
	if (pr_id == 0) {
		auto endefr = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
		writeoutput << "final file(s) read...\t" << elapsed_secondsfr.count() << "\n";
		writeoutput.flush();
	}
	if (pr_id != 0 && (FullMonitoring == 1 || FullMonitoring == 2)) {
		auto endefr = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
		writeoutput_sub << "final file(s) read...\t" << elapsed_secondsfr.count() << "\n";
		writeoutput_sub.flush();
	}

	readpara >> sbuf;

	readpara >> sbuf >> NUserdefinedTerms;
	for (int i = 0; i < NUserdefinedTerms; i++) {
		USERDEFINEDTERM objUSERDEFINEDTERM;
		readpara >> sbuf >> objUSERDEFINEDTERM.TermName;
		int NTermComposites;
		readpara >> sbuf >> NTermComposites;
		readpara >> sbuf;
		for (int j = 0; j < NTermComposites; j++) {
			COMPOSITE objCOMPOSITE;
			readpara >> sbuf >> objCOMPOSITE.coefficient >> objCOMPOSITE.index1 >> objCOMPOSITE.index2;
			readpara >> sbuf;
			if (sbuf == "re") {
				objCOMPOSITE.reim = 0;
			}
			else {
				objCOMPOSITE.reim = 1;
			}
			objUSERDEFINEDTERM.composite.push_back(objCOMPOSITE);
		}
		userdefinedterm.push_back(objUSERDEFINEDTERM);
	}

	readpara.close();

}

void get_label() {

	for (int i = 0; i < (signed)ketseries.size(); i++) {
		for (int j = 0; j < (signed)ketseries[i].ketsss.size(); j++) {
			for (int k = 0; k < (signed)ketseries[i].ketsss[j].ketss.size(); k++) {
				for (unsigned long long int l = 0; l < (signed)ketseries[i].ketsss[j].ketss[k].kets.size(); l++) {
					ketseries[i].ketsss[j].ketss[k].kets[l].label = l;
				}
			}
		}
	}

	globalketstring_max = 1;
	for (int j = 0; j < (signed)ketseries[0].ketsss.size(); j++) {
		for (int k = 0; k < (signed)ketseries[0].ketsss[j].ketss.size(); k++) {
			globalketstring_max *= (signed)ketseries[0].ketsss[j].ketss[k].kets.size();
		}
	}

	globalbrastring_max = 1;
	for (int j = 0; j < (signed)ketseries[ketseries.size() - 1].ketsss.size(); j++) {
		for (int k = 0; k < (signed)ketseries[ketseries.size() - 1].ketsss[j].ketss.size(); k++) {
			globalbrastring_max *= (signed)ketseries[ketseries.size() - 1].ketsss[j].ketss[k].kets.size();
		}
	}

}

void set_MPI() {

	if (ManualMPI == 0) {

		signed long long int max_checker = -10000;
		int target_i;
		for (int i = 1; i < (signed)ketseries.size() - 1; i++) {
			for (int j = 0; j < (signed)ketseries[i].ketsss.size(); j++) {
				for (int k = 0; k < (signed)ketseries[i].ketsss[j].ketss.size(); k++) {
					if ((signed)ketseries[i].ketsss[j].ketss[k].kets.size() > max_checker) {
						max_checker = (signed)ketseries[i].ketsss[j].ketss[k].kets.size();
						target_i = i;
					}
				}
			}
		}

		int target_k;
		for (int j = 0; j < (signed)ketseries[target_i].ketsss.size(); j++) {
			signed long long int max_checker2 = -10000;
			for (int k = 0; k < (signed)ketseries[target_i].ketsss[j].ketss.size(); k++) {
				if ((signed)ketseries[target_i].ketsss[j].ketss[k].kets.size() > max_checker2) {
					max_checker2 = (signed)ketseries[target_i].ketsss[j].ketss[k].kets.size();
					target_k = k;
				}
			}
			std::vector<KETS> ketsbuf;
			for (signed long long int idist = pr_id; idist < (signed)ketseries[target_i].ketsss[j].ketss[target_k].kets.size(); idist += pr_size) {
				ketsbuf.push_back(ketseries[target_i].ketsss[j].ketss[target_k].kets[idist]);
			}
			ketseries[target_i].ketsss[j].ketss[target_k].kets.clear();
			for (signed long long int idistagain = 0; idistagain < (signed)ketsbuf.size(); idistagain++) {
				ketseries[target_i].ketsss[j].ketss[target_k].kets.push_back(ketsbuf[idistagain]);
			}
			if (ketsbuf.size() == 0) {
				mpi_activated = false;
				std::cout << "[ERROR] pr_id at " << pr_id << " ! Decrease #MPI!\n";
				std::cout.flush();
			}
		}

	}
	else {

		readmanualmpi.open(ManualMPIFilename);

		while (true) {
			int ibuf;
			readmanualmpi >> sbuf >> ibuf;
			if (sbuf == "pr_id") {
				if (ibuf == pr_id) {
					break;
				}
			}
		}

		for (int i = 1; i < (signed)ketseries.size() - 1; i++) {
			for (int j = 0; j < (signed)ketseries[i].ketsss.size(); j++) {
				for (int k = 0; k < (signed)ketseries[i].ketsss[j].ketss.size(); k++) {
					std::vector<KETS> ketsbuf;
					int iset0, iset1;
					readmanualmpi >> iset0 >> iset1;
					if (!(iset0 == -1 || iset1 == -1)) {
						for (signed long long int idist = iset0; idist <= iset1; idist++) {
							ketsbuf.push_back(ketseries[i].ketsss[j].ketss[k].kets[idist]);
						}
						ketseries[i].ketsss[j].ketss[k].kets.clear();
						for (signed long long int idistagain = 0; idistagain < (signed)ketsbuf.size(); idistagain++) {
							ketseries[i].ketsss[j].ketss[k].kets.push_back(ketsbuf[idistagain]);
						}
					}
					else {
						ketseries[i].ketsss[j].ketss[k].kets.clear();
					}
				}
			}
		}

		readmanualmpi.close();

		bool zerokets = true;
		for (int i = 1; i < (signed)ketseries.size() - 1; i++) {
			for (int j = 0; j < (signed)ketseries[i].ketsss.size(); j++) {
				for (int k = 0; k < (signed)ketseries[i].ketsss[j].ketss.size(); k++) {
					if (ketseries[i].ketsss[j].ketss[k].kets.size() != 0) {
						zerokets = false;
						break;
					}
				}
				if (zerokets == false) {
					break;
				}
			}
			if (zerokets == false) {
				break;
			}
		}

		if (zerokets == true) {
			mpi_activated = false;
			std::cout << "[ERROR] pr_id at " << pr_id << " ! Set Manual_#MPI again!\n";
			std::cout.flush();
		}
	
	}

}

void get_prtbmatrices() {

	for (int i = 0; i < (signed)ketseries.size() - 1; i++) {
		PTRBMATRIXSRS objPTRBMATRIXSRS;
		int iket = i;
		int ibra = i + 1;
		if (SaveMemory == 0) {
			get_prtbterm_simply(&ibra, &iket, &Ht[i], &objPTRBMATRIXSRS.ptrbmatrix);
		}
		else {
			get_perturbation(&ibra, &iket, &Ht[i], &objPTRBMATRIXSRS.ptrbmatrix);
		}
		ptrbmatrixsrs.push_back(objPTRBMATRIXSRS);
		if (pr_id == 0) {
			auto endefr = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
			writeoutput << "ptrbmatrixsrs[" << i << "] created...\t" << elapsed_secondsfr.count() << "\n";
			writeoutput.flush();
		}
		if (pr_id != 0 && (FullMonitoring == 1 || FullMonitoring == 2)) {
			auto endefr = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
			writeoutput_sub << "ptrbmatrixsrs[" << i << "] created...\t" << elapsed_secondsfr.count() << "\n";
			writeoutput_sub.flush();
		}
	}

}

void get_final_prtb_term() {

	long double prtb_re = 1.0;
	long double prtb_im = 0.0;

	for (int i = 0; i < (signed)ptrbmatrixsrs.size(); i++) {
		long double buf_re = prtb_re;
		long double buf_im = prtb_im;
		signed long long int get_ii = Fdiagram[i];
		signed long long int get_jj = Fdiagram[i + 1];
		long double term_re = ptrbmatrixsrs[i].ptrbmatrix[get_ii].ptrb[get_jj].v[0];
		long double term_im = ptrbmatrixsrs[i].ptrbmatrix[get_ii].ptrb[get_jj].v[1];
		ComplexMTwoRe(&buf_re, &buf_im, &term_re, &term_im, &prtb_re);
		ComplexMTwoIm(&buf_re, &buf_im, &term_re, &term_im, &prtb_im);
	}

	signed long long int get_i = Fdiagram[0];
	signed long long int get_j = Fdiagram[(signed)Fdiagram.size() - 1];

	final_prtb_terms[get_i][get_j][0] += prtb_re;
	final_prtb_terms[get_i][get_j][1] += prtb_im;

	if (FullMonitoring == 2) {
		for (int i = (signed)Fdiagram.size() - 1; i >= 0; i--) {
			if (SaveMemory == 1 && GetLabel == 1) {
				if (i == 0) {
					if (pr_id == 0) {
						writeoutput << globalketstring[Fdiagram[i]] << "\t";
						writeoutput.flush();
					}
					else {
						writeoutput_sub << globalketstring[Fdiagram[i]] << "\t";
						writeoutput_sub.flush();
					}
				}
				else if (i == (signed)Fdiagram.size() - 1) {
					if (pr_id == 0) {
						writeoutput << globalbrastring[Fdiagram[i]] << "\t";
						writeoutput.flush();
					}
					else {
						writeoutput_sub << globalbrastring[Fdiagram[i]] << "\t";
						writeoutput_sub.flush();
					}
				}
				else {
					if (pr_id == 0) {
						writeoutput << ketseries[i].labels[Fdiagram[i]] << "\t";
						writeoutput.flush();
					}
					else {
						writeoutput_sub << ketseries[i].labels[Fdiagram[i]] << "\t";
						writeoutput_sub.flush();
					}
				}
			}
			else {
				if (pr_id == 0) {
					writeoutput << Fdiagram[i] << "\t";
					writeoutput.flush();
				}
				else {
					writeoutput_sub << Fdiagram[i] << "\t";
					writeoutput_sub.flush();
				}
			}
		}
		if (pr_id == 0) {
			writeoutput << ":\t" << prtb_re << "\t+I*\t" << prtb_im << "\n";
			writeoutput.flush();
		}
		else {
			writeoutput_sub << ":\t" << prtb_re << "\t+I*\t" << prtb_im << "\n";
			writeoutput_sub.flush();
		}
	}

}

void get_Fdiagram_rc(int depth) {

	signed long long int get_vector_size;
	if (depth < (signed)ptrbmatrixsrs.size()) {
		get_vector_size = (signed)ptrbmatrixsrs[depth].ptrbmatrix.size();
	}
	else {
		get_vector_size = (signed)ptrbmatrixsrs[depth - 1].ptrbmatrix[0].ptrb.size();
	}

	for (signed long long int i = 0; i < get_vector_size; i++) {
		Fdiagram[depth] = i;
		if (depth != (signed)ptrbmatrixsrs.size()) {
			get_Fdiagram_rc(depth + 1);
		}
		else {
			get_final_prtb_term();
		}
	}

}

void calculate_paths() {

	if (FullMonitoring == 2) {
		if (pr_id == 0) {
			writeoutput << "======================calculate_paths()======================\n";
			writeoutput.flush();
		}
		else {
			writeoutput_sub << "======================calculate_paths()======================\n";
			writeoutput_sub.flush();
		}
	}

	final_prtb_terms = new long double** [(signed)ptrbmatrixsrs[0].ptrbmatrix.size()];
	global_final_prtb_terms = new long double** [(signed)ptrbmatrixsrs[0].ptrbmatrix.size()];
	for (signed long long int i = 0; i < (signed)ptrbmatrixsrs[0].ptrbmatrix.size(); i++) {
		final_prtb_terms[i] = new long double* [(signed)ptrbmatrixsrs[(signed)ptrbmatrixsrs.size() - 1].ptrbmatrix[0].ptrb.size()];
		global_final_prtb_terms[i] = new long double* [(signed)ptrbmatrixsrs[(signed)ptrbmatrixsrs.size() - 1].ptrbmatrix[0].ptrb.size()];
		for (signed long long int j = 0; j < (signed)ptrbmatrixsrs[(signed)ptrbmatrixsrs.size() - 1].ptrbmatrix[0].ptrb.size(); j++) {
			final_prtb_terms[i][j] = new long double[2];
			global_final_prtb_terms[i][j] = new long double[2];
			for (int k = 0; k < 2; k++) {
				final_prtb_terms[i][j][k] = 0.0;
				global_final_prtb_terms[i][j][k] = 0.0;
			}
		}
	}

	Fdiagram.clear();
	for (int i = 0; i < (signed)ptrbmatrixsrs.size() + 1; i++) {
		Fdiagram.push_back(-1);
	}

	get_Fdiagram_rc(0);

	MPI_Barrier(MPI_COMM_WORLD);

	for (signed long long int i = 0; i < (signed)ptrbmatrixsrs[0].ptrbmatrix.size(); i++) {
		for (signed long long int j = 0; j < (signed)ptrbmatrixsrs[(signed)ptrbmatrixsrs.size() - 1].ptrbmatrix[0].ptrb.size(); j++) {
			for (int k = 0; k < 2; k++) {
				MPI_Allreduce(&final_prtb_terms[i][j][k], &global_final_prtb_terms[i][j][k], 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for (signed long long int i = 0; i < (signed)ptrbmatrixsrs[0].ptrbmatrix.size(); i++) {
		for (signed long long int j = 0; j < (signed)ptrbmatrixsrs[(signed)ptrbmatrixsrs.size() - 1].ptrbmatrix[0].ptrb.size(); j++) {
			for (int k = 0; k < 2; k++) {
				if (ChopBool(&global_final_prtb_terms[i][j][k])) {
					global_final_prtb_terms[i][j][k] = 0.0;
				}
			}
		}
	}

	if (pr_id == 0) {
		auto endefr = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
		writeoutput << "finished calculating paths...\t " << elapsed_secondsfr.count() << "\n============================================\n";
		writeoutput.flush();
	}
	if (pr_id != 0 && (FullMonitoring == 1 || FullMonitoring == 2)) {
		auto endefr = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
		writeoutput_sub << "finished calculating paths...\t " << elapsed_secondsfr.count() << "\n============================================\n";
		writeoutput_sub.flush();
		writeoutput_sub.close();
	}

}

void write_prtb_results() {

	if (FullMonitoring == 2) {

		intermeidateanalysis = new INTERMEDIATEANALYSIS *[ptrbmatrixsrs[0].ptrbmatrix.size()];
		for (signed long long int i = 0; i < (signed)ptrbmatrixsrs[0].ptrbmatrix.size(); i++) {
			intermeidateanalysis[i] = new INTERMEDIATEANALYSIS [ptrbmatrixsrs[(signed)ptrbmatrixsrs.size() - 1].ptrbmatrix[0].ptrb.size()];
		}

		for (int ipr = 0; ipr < pr_size; ipr++) {

			if (ipr == 0) {
				readoutput_sub.open(OutputFilename + ".txt");
			}
			else {
				readoutput_sub.open("_" + std::to_string(ipr) + "_" + OutputFilename + ".txt");
			}
			
			while (true) {

				bool catched_finish = false;
				readoutput_sub >> sbuf;
				if (sbuf == "======================calculate_paths()======================") {

					while (true) {

						readoutput_sub >> sbuf;

						if (sbuf == "finished" || readoutput_sub.eof()) {

							catched_finish = true;
							break;

						}
						else {

							int target_j, target_i;
							std::vector<std::string> strbuf;
							INTERMEDIATECOMP objINTERMEDIATECOMP;

							if (SaveMemory == 1 && GetLabel == 1) {
								for (signed long long int j = 0; j < (signed)ptrbmatrixsrs[(signed)ptrbmatrixsrs.size() - 1].ptrbmatrix[0].ptrb.size(); j++) {
									if (sbuf == globalbrastring[j]) {
										target_j = j;
										break;
									}
								}
							}
							else {
								std::stringstream lcstream;
								lcstream << sbuf;
								lcstream >> target_j;
							}

							for (signed long long int j = 0; j < (signed)ptrbmatrixsrs.size() - 1; j++) {
								readoutput_sub >> sbuf;
								objINTERMEDIATECOMP.strterm.push_back(sbuf);
							}

							readoutput_sub >> sbuf;
							if (SaveMemory == 1 && GetLabel == 1) {
								for (signed long long int j = 0; j < (signed)ptrbmatrixsrs[0].ptrbmatrix.size(); j++) {
									if (sbuf == globalketstring[j]) {
										target_i = j;
										break;
									}
								}
							}
							else {
								std::stringstream lcstream;
								lcstream << sbuf;
								lcstream >> target_i;
							}

							readoutput_sub >> sbuf >> objINTERMEDIATECOMP.reim[0] >> sbuf >> objINTERMEDIATECOMP.reim[1];

							for (signed long long int j = 0; j < (signed)ptrbmatrixsrs.size() - 1; j++) {
								intermeidateanalysis[target_i][target_j].intermediatecomp.push_back(objINTERMEDIATECOMP);
							}

						}
						
					}

					if (catched_finish == true || readoutput_sub.eof()) {
						break;
					}

				}

				if (catched_finish == true || readoutput_sub.eof()) {
					break;
				}

			}
			readoutput_sub.close();
		}

	}

	for (signed long long int i = 0; i < (signed)ptrbmatrixsrs[0].ptrbmatrix.size(); i++) {
		for (signed long long int j = 0; j < (signed)ptrbmatrixsrs[(signed)ptrbmatrixsrs.size() - 1].ptrbmatrix[0].ptrb.size(); j++) {

			if (FullMonitoring == 2) {
				writeoutput_decomp.open(OutputFilename + ".ptrbmatrixsrs." + std::to_string(i + 1) + "." + std::to_string(j + 1) + ".decomposition.txt");
				writeoutput_decomp << std::fixed << std::setprecision(OutputPrecision);
			}
			
			if (!(SaveMemory == 1 && GetLabel == 1)) {
				writeoutput << "<\t" << j + 1 << "\t|...(<-)...|\t" << i + 1 << "\t>=\t" << global_final_prtb_terms[i][j][0] << "\t+I*\t" << global_final_prtb_terms[i][j][1];
				writeoutput.flush();
				if (FullMonitoring == 2) {
					writeoutput_decomp << "<\t" << j + 1 << "\t|...(<-)...|\t" << i + 1 << "\t>=\t" << global_final_prtb_terms[i][j][0] << "\t+I*\t" << global_final_prtb_terms[i][j][1];
					writeoutput_decomp.flush();
				}
			}
			else {
				writeoutput << "<\t" << globalbrastring[j] << "\t|...(<-)...|\t" << globalketstring[i] << "\t>=\t" << global_final_prtb_terms[i][j][0] << "\t+I*\t" << global_final_prtb_terms[i][j][1];
				writeoutput.flush();
				if (FullMonitoring == 2) {
					writeoutput_decomp << "<\t" << globalbrastring[j] << "\t|...(<-)...|\t" << globalketstring[i] << "\t>=\t" << global_final_prtb_terms[i][j][0] << "\t+I*\t" << global_final_prtb_terms[i][j][1];
					writeoutput_decomp.flush();
				}
			}

			if (FullMonitoring == 2) {
				writeoutput_decomp << "\n";
				writeoutput_decomp.flush();
				for (signed long long int k = 0; k < (signed)intermeidateanalysis[i][j].intermediatecomp.size(); k++) {
					for (signed long long int l = 0; l < (signed)intermeidateanalysis[i][j].intermediatecomp[k].strterm.size(); l++) {
						writeoutput_decomp << intermeidateanalysis[i][j].intermediatecomp[k].strterm[l] << "\t";
						writeoutput_decomp.flush();
					}
					writeoutput_decomp << intermeidateanalysis[i][j].intermediatecomp[k].reim[0] << "\t+I*\t" << intermeidateanalysis[i][j].intermediatecomp[k].reim[1] << "\n";
					writeoutput_decomp.flush();
				}
			}

			if (!(i == (signed)ptrbmatrixsrs[0].ptrbmatrix.size() - 1 && j == (signed)ptrbmatrixsrs[(signed)ptrbmatrixsrs.size() - 1].ptrbmatrix[0].ptrb.size() - 1)) {
				writeoutput << "\n";
				writeoutput.flush();
				if (FullMonitoring == 2) {
					writeoutput_decomp << "\n";
					writeoutput_decomp.flush();
				}
			}

			if (FullMonitoring == 2) {
				writeoutput_decomp.close();
			}
		
		}
	}

}

void write_userdefinedterms() {

	std::vector<long double> terms;

	writeoutput << "\n============================================\n";
	writeoutput.flush();
	for (int i = 0; i < (signed)userdefinedterm.size(); i++) {
		writeoutput << userdefinedterm[i].TermName << "\t";
		writeoutput.flush();
		long double term = 0.0;
		for (int j = 0; j < (signed)userdefinedterm[i].composite.size(); j++) {
			long double addterm = userdefinedterm[i].composite[j].coefficient;
			addterm *= global_final_prtb_terms[userdefinedterm[i].composite[j].index2 - 1][userdefinedterm[i].composite[j].index1 - 1][userdefinedterm[i].composite[j].reim];
			term += addterm;
		}
		if (ChopBool(&term) == true) {
			term = 0.0;
		}
		terms.push_back(term);
		writeoutput << term;
		writeoutput.flush();
		if (i != (signed)userdefinedterm.size() - 1) {
			writeoutput << "\n";
			writeoutput.flush();
		}
	}

	if (FullMonitoring == 2) {
		for (int i = 0; i < (signed)userdefinedterm.size(); i++) {
			std::ofstream writeoutput_decomp_sub;
			writeoutput_decomp_sub.open(OutputFilename + "." + userdefinedterm[i].TermName + ".decomposition.txt");
			writeoutput_decomp_sub << std::fixed << std::setprecision(OutputPrecision);
			writeoutput_decomp_sub << userdefinedterm[i].TermName << "\t" << terms[i] << "\n";
			writeoutput_decomp_sub.flush();
			for (int ii = 0; ii < (signed)intermeidateanalysis[0][0].intermediatecomp.size(); ii++) {
				long double term = 0.0;
				for (int j = 0; j < (signed)userdefinedterm[i].composite.size(); j++) {
					long double addterm = userdefinedterm[i].composite[j].coefficient;
					addterm *= intermeidateanalysis[userdefinedterm[i].composite[j].index2 - 1][userdefinedterm[i].composite[j].index1 - 1].intermediatecomp[ii].reim[userdefinedterm[i].composite[j].reim];
					term += addterm;
				}
				if (ChopBool(&term) == true) {
					term = 0.0;
				}
				for (int iii = 0; iii < (signed)intermeidateanalysis[0][0].intermediatecomp[ii].strterm.size(); iii++) {
					writeoutput_decomp_sub << intermeidateanalysis[0][0].intermediatecomp[ii].strterm[iii] << "\t";
				}
				writeoutput_decomp_sub << term;
				writeoutput_decomp_sub.flush();
				if (ii != (signed)intermeidateanalysis[0][0].intermediatecomp.size() - 1) {
					writeoutput_decomp_sub << "\n";
					writeoutput_decomp_sub.flush();
				}
			}
			writeoutput_decomp_sub.close();
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

	srand((unsigned)time(0));
	startefr = std::chrono::system_clock::now();

#ifdef _WIN32
	if (pr_id == 0) {
		std::cout << "[LOG] WINDOWS detected...\n";
		std::cout.flush();
	}
#elif _WIN64
	if (pr_id == 0) {
		std::cout << "[LOG] WINDOWS detected...\n";
		std::cout.flush();
	}
#else
	if (pr_id == 0) {
		std::cout << "[LOG] WINDOWS not detected...\n";
		std::cout.flush();
	}
#endif

	read_para();
	MPI_Barrier(MPI_COMM_WORLD);

	if (SaveMemory == 1 && GetLabel == 1) {
		get_label();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	set_MPI();
	MPI_Barrier(MPI_COMM_WORLD);

	if (mpi_activated == false) {

		return -1;
		MPI_Finalize();

	}
	else {

		if (ReadPrtb == 1) {
			std::string filename_readprtb = "__" + std::to_string(pr_id) + "_" + OutputFilename + ".dat";
			start_read_prtbterm(&filename_readprtb);
		}

		get_prtbmatrices();
		writeoutput << "waiting on the next MPI_Barrier...\n";
		writeoutput.flush();
		if (pr_id != 0 && (FullMonitoring == 1 || FullMonitoring == 2)) {
			writeoutput_sub << "waiting on the next MPI_Barrier...\n";
			writeoutput_sub.flush();
		}
		MPI_Barrier(MPI_COMM_WORLD);
		auto endefr = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsfr = endefr - startefr;
		writeoutput << "passed the MPI_Barrier...\t" << elapsed_secondsfr.count() << "\n";
		writeoutput.flush();
		if (pr_id != 0 && (FullMonitoring == 1 || FullMonitoring == 2)) {
			writeoutput_sub << "passed the MPI_Barrier...\t" << elapsed_secondsfr.count() << "\n";
			writeoutput_sub.flush();
		}
		calculate_paths();
		MPI_Barrier(MPI_COMM_WORLD);

		if (pr_id == 0) {
			write_prtb_results();
			if (userdefinedterm.size() > 0) {
				write_userdefinedterms();
			}
			writeoutput.close();
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (WritePrtb == 1) {
			close_write_prtbterm();
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();

	}

	return 0;

}