#pragma once
#include<Eigen/Dense>
#include"transform.h"
#include<string>
#include<vector>
using namespace Eigen;
using namespace std;

#define GPS_SAT_QUAN 32
#define BDS_SAT_QUAN 63

/*������*/
#define UN_Solve 0
#define Success 1
#define Epoch_Loss -1
#define OBS_DATA_Loss -2
#define Set_UP_B_fail -3
/*�������洢��*/
class Result_DATA
{
public:
	GPSTIME* OBSTIME;	//�۲�ʱ��
	MatrixXd* Pos;		//λ��+�Ӳ�
	MatrixXd* Vel;		//�ٶ�+����
	MatrixXd* Q_Pos;	//Pos��Q����
	MatrixXd* Q_Vel;	//Vel��Q����
	double* thegma_Pos;	//Pos�ĵ�λȨ�����
	double* thegma_Vel;	//Vel�ĵ�λȨ�����
	double* PDOP;		//λ�õ�DOPֵ
	double* VDOP;		//�ٶȵ�DOPֵ
	int GPS_num;		//GPS������
	int BDS_num;		//BDS������
	string* SATES;		//ʹ�õ����Ǽ���
	int solve_result;	//������
	XYZ* Real_Pos;		//�ο���ֵ

	Result_DATA();
	int OUTPUT();					//����̨���
	int WRITEOUTPUT(FILE* fpr);		//�ļ����
};
/*���ںʹ洢�ṹ�������*/
class Configure
{
public:
	static const char* NetIP;				//IP
	static const unsigned short NetPort;	//�˿�
	const char* ObsDatFile;					//log�ļ�·��
	const char* ResDatFile;					//pos�ļ�·��
};

/*�����ļ���*/
int createDirectory(string path);