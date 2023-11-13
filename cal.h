#pragma once
#include"read.h"
#include"transform.h"
#include"data.h"
using namespace std;

/*大地常数*/
#define WGS84_GM 3.986005E+14
#define CGCS2000_GM 3.986004418E+14
#define omiga_earth 7.2921151467E-05
/*光速*/
#define velocity_c 2.99792458E8

/*code IDs*/
#define UNKOWN          0
/*GPS*/
#define CODE_L1C		1
#define CODE_L2P		2
#define CODE_L2W		3
#define CODE_L5Q		4
#define CODE_L1L		5
#define CODE_L2S		6
/*BDS*/
#define CODE_L2I		7
#define CODE_L7I		8
#define CODE_L6I		9
#define CODE_L1P		10
#define CODE_L5P		11

/*FREQUNCY*/          //MHz
/*GPS*/
#define L1   1575.42
#define L2   1227.60
#define L5   1176.45
/*BDS*/
#define B1   1561.098
#define B1_C 1575.42
#define B2   1207.14
#define B2_a 1176.45
#define B3   1268.52

/*Hopefiled*/
#define H0 0
#define T0 288.16
#define P0 1013.25
#define RH0 0.5

/*粗茶探测阈值*/
#define GF_THRESH 0.05
#define MW_THRESH 10

/*初始化相关变量*/
unsigned int initial();
/*平方数*/
double SQR(double x);
/*模长*/
double Len(XYZ* pos);

/*角度单位转换*/
double degree2rad(double degree);

double rad2degree(double rad);
/*计算DOP值*/
double Cal_PDOP(MatrixXd Qxx);
/*最小二乘计算*/
unsigned int Cal_LEAST_SQR(MatrixXd B, MatrixXd l, MatrixXd P, MatrixXd* Qxx, MatrixXd& x, double* thegma, double* DOP);
/*解码频点*/
unsigned int decode_SYN(int sys, int signal);
/*解码频率*/
double CODE2FREQ(int code);

//钟差改正
double CORRECT_CLK(double t, EPHEMERIS* eph);

//星历位置
unsigned int SAT_POS_CAL(double t, EPHEMERIS* eph, XYZ* xyz, double& clk, double dt, int SYS);

//卫星高度角计算
double Ele_Angle(XYZ* SatPos, XYZ* RcvPos, int sys);

//Hopefiled对流层改正(m)
double Hopefield(double E, double H);

//Hopefiled对流层改正(m)
double Hopefield(XYZ* SatPos, XYZ* RcvPos, int sys);

/*Klobuchar模型改正电离层*/
double Klobuchar(XYZ* RcvPos, double E, double A, double alpha[4], double beta[4], double UT, double code, int sys);

/*一致性检验*/
int CheckOBSConsist(Satellate* sate, int sys, double t, int index, bool& PSE_flag, bool& PHA_flag);

/*粗差探测*/
int DetectOutlier(Satellate* sate, int sys, double t, int index1, int index2);

/*文件流下SPP解算*/
//搭建位置解算构造矩阵
unsigned int setup_Pos(GPSTIME* OBS_TIME, MatrixXd Pos, vector<Satellate*> Sates, EPOCH* eph, bool first_flag, double f1, double f2, MatrixXd* B_Pos, MatrixXd* l_Pos, MatrixXd* P_Pos, string* sate_used);

//搭建速度解算构造矩阵
unsigned int setup_Vel(GPSTIME* OBS_TIME, MatrixXd Pos, vector<Satellate*> Sates, EPOCH* eph, MatrixXd* B_Vel, MatrixXd* l_Vel, MatrixXd* P_Vel);

//GPS、BDS双星解算
unsigned int Cal_2(Result_DATA* data, OBS_DATA* obs, EPOCH* gpseph, EPOCH* bdseph, bool first_flag);

//SPP单点定位
unsigned int Cal_SPP(Result_DATA* data, OBS_DATA* obs, EPOCH* gpseph, EPOCH* bdseph, double dt_e, bool first_flag);

/*网口下位置解算*/
//搭建位置解算构造矩阵
unsigned int setup_Pos(GPSTIME* OBS_TIME, MatrixXd Pos, vector<Satellate*> Sates, EPHEMERIS** eph, bool first_flag, double f1, double f2, MatrixXd* B_Pos, MatrixXd* l_Pos, MatrixXd* P_Pos, string* sate_used);

//搭建速度解算构造矩阵
unsigned int setup_Vel(GPSTIME* OBS_TIME, MatrixXd Pos, vector<Satellate*> Sates, EPHEMERIS** eph, MatrixXd* B_Vel, MatrixXd* l_Vel, MatrixXd* P_Vel);

//GPS、BDS双星解算
unsigned int Cal_2(Result_DATA* data, OBS_DATA* obs, EPHEMERIS** gpseph, EPHEMERIS** bdseph, bool first_flag);

//SPP单点定位
unsigned int Cal_SPP(Result_DATA* data, OBS_DATA* obs, EPHEMERIS** gpseph, EPHEMERIS** bdseph, double dt_e, bool first_flag);

/*网口下解算*/
int decodestream(Result_DATA* result, unsigned char Buff[], int& d);