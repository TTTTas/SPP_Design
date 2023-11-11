#pragma once
#include<Eigen/Dense>
#include"transform.h"
#include<string>
#include<vector>
using namespace Eigen;
using namespace std;

#define GPS_SAT_QUAN 32
#define BDS_SAT_QUAN 63

/*?????????*/
#define UN_Solve 0
#define Success 1
#define Epoch_Loss -1
#define OBS_DATA_Loss -2
#define Set_UP_B_fail -3

class Result_DATA
{
public:
	GPSTIME* OBSTIME;
	MatrixXd* Pos;
	MatrixXd* Vel;
	MatrixXd* Q_Pos;
	MatrixXd* Q_Vel;
	double* thegma_Pos;
	double* thegma_Vel;
	double* PDOP;
	double* VDOP;
	int GPS_num;
	int BDS_num;
	string* SATES;
	int solve_result;
	XYZ* Real_Pos;

	Result_DATA();
	int OUTPUT();
	int WRITEOUTPUT(FILE* fpr);
};

class Configure
{
public:
	static const char* NetIP;
	static const unsigned short NetPort;
	const char* ObsDatFile;
	const char* ResDatFile;
};

int createDirectory(string path);
