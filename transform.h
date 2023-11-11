#include <iostream>
#include <Eigen/Dense>
#include <vector>
using namespace std;
using namespace Eigen;

#define WGS84_e2 0.0066943799013
#define WGS84_a 6378137.0
#define CGCS2000_e2 0.00669438002290
#define CGCS2000_a 6378137.0
#define Pi 3.1415926
/*system IDs*/
#define SYS_GPS         0
#define SYS_BDS         4



#pragma region 时间域
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
struct XYZ
{
	double X;
	double Y;
	double Z;
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

struct BLH
{
	double Lon;//经度
	double Lat;//纬度
	double Height;

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

void BLH2XYZ(BLH* blh, XYZ* xyz, double e2, double a);

void XYZ2BLH(BLH* blh, XYZ* xyz, double e2, double a);

//转化为当地坐标系
double XYZ2ENU(XYZ* xyz1, XYZ* xyz2, XYZ* enu, int sys);
#pragma endregion


#pragma once
