#include <iostream>
#include <Eigen/Dense>
#include <vector>
using namespace std;
using namespace Eigen;

/*地球椭球相关参数*/
#define WGS84_e2 0.0066943799013
#define WGS84_a 6378137.0
#define CGCS2000_e2 0.00669438002290
#define CGCS2000_a 6378137.0
#define Pi 3.1415926
/*system IDs*/
#define SYS_GPS         0
#define SYS_BDS         4



#pragma region 时间域
/*世界时*/
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
/*简化儒略日*/
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
/*GPS时*/
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
/*北斗时*/
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
/*时间相互转化函数*/
void TIME2MJD(UTC* utc, MJD* mjd);

void MJD2TIME(UTC* utc, MJD* mjd);

void MJD2GPSTIME(MJD* mjd, GPSTIME* gpst);

void GPSTIME2MJD(MJD* mjd, GPSTIME* gpst);

void GPSTIME2BDSTIME(GPSTIME* gpst, GPSTIME* bdst);

void GPSTIME2BDSTIME(GPSTIME* gpst, BDSTIME* bdst);

void BDSTIME2GPSTIME(GPSTIME* gpst, BDSTIME* bdst);

void GPSTIME2TIME(GPSTIME* gpst, UTC* utc);


#pragma endregion

#pragma region 坐标域
/*XYZ坐标系*/
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
/*大地坐标系*/
struct BLH
{
	double Lon;				//经度,deg
	double Lat;				//纬度,deg
	double Height;			//高程,m

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
/*坐标系转换*/
void BLH2XYZ(BLH* blh, XYZ* xyz, double e2, double a);

void XYZ2BLH(BLH* blh, XYZ* xyz, double e2, double a);

//转化为当地坐标系
/*****************************/
//	xyz1: 站心在xyz坐标系下坐标
//	xyz2: 目标在xyz坐标系下坐标
double XYZ2ENU(XYZ* xyz1, XYZ* xyz2, XYZ* enu, int sys);
#pragma endregion


#pragma once
