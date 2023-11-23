#include <iostream>
#include <Eigen/Dense>
#include <vector>
using namespace std;
using namespace Eigen;

/*����������ز���*/
#define WGS84_e2 0.0066943799013
#define WGS84_a 6378137.0
#define CGCS2000_e2 0.00669438002290
#define CGCS2000_a 6378137.0
#define Pi 3.1415926
/*system IDs*/
#define SYS_GPS         0
#define SYS_BDS         4

#pragma region ʱ����
/*����ʱ*/
struct UTC
{
	unsigned short Year;		
	unsigned short Month;
	unsigned short Day;
	unsigned short Hour;
	unsigned short Min;
	unsigned short Sec;

	UTC(unsigned short y, unsigned short m, unsigned short d, unsigned short h, unsigned short min, unsigned short s)
	{
		Year = y;
		Month = m;
		Day = d;
		Hour = h;
		Min = min;
		Sec = s;
	}
	UTC()
	{
		Year = Month = Day = Hour = Min = Sec = 0;
	}
};
/*��������*/
struct MJD
{
	int Days;
	double FracDay;

	MJD()
	{
		Days = 0;
		FracDay = 0.0;
	}
};
/*GPSʱ*/
struct GPSTIME
{
	unsigned short Week;
	double SecOfWeek;

	GPSTIME()
	{
		Week = 0;
		SecOfWeek = 0.0;
	}
	GPSTIME(unsigned short w, double s)
	{
		Week = w;
		SecOfWeek = s;
	}
};
/*����ʱ*/
struct BDSTIME
{
	unsigned short Week;
	double SecOfWeek;

	BDSTIME()
	{
		Week = 0;
		SecOfWeek = 0.0;
	}
};
/*ʱ���໥ת������*/
void TIME2MJD(UTC* utc, MJD* mjd);

void MJD2TIME(UTC* utc, MJD* mjd);

void MJD2GPSTIME(MJD* mjd, GPSTIME* gpst);

void GPSTIME2MJD(MJD* mjd, GPSTIME* gpst);

void GPSTIME2BDSTIME(GPSTIME* gpst, GPSTIME* bdst);

void GPSTIME2BDSTIME(GPSTIME* gpst, BDSTIME* bdst);

void BDSTIME2GPSTIME(GPSTIME* gpst, BDSTIME* bdst);

void GPSTIME2TIME(GPSTIME* gpst, UTC* utc);


#pragma endregion

#pragma region ������
/*XYZ����ϵ*/
struct XYZ
{
	double X;		//m
	double Y;		//m
	double Z;		//m
	XYZ()
	{
		X = Y = Z = 0;
	}
	XYZ(double x, double y, double z)
	{
		X = x;
		Y = y;
		Z = z;
	}
};
/*�������ϵ*/
struct BLH
{
	double Lon;				//����,deg
	double Lat;				//γ��,deg
	double Height;			//�߳�,m

	BLH()
	{
		Lon = Lat = Height = 0;
	}
	BLH(double lon, double lat, double h)
	{
		Lon = lon;
		Lat = lat;
		Height = h;
	}
};
/*����ϵת��*/
void BLH2XYZ(BLH* blh, XYZ* xyz, double e2, double a);

void XYZ2BLH(BLH* blh, XYZ* xyz, double e2, double a);

//ת��Ϊ��������ϵ
/*****************************/
//	xyz1: վ����xyz����ϵ������
//	xyz2: Ŀ����xyz����ϵ������
double XYZ2ENU(XYZ* xyz1, XYZ* xyz2, XYZ* enu, int sys);
#pragma endregion


#pragma once
