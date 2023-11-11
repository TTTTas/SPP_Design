#include"transform.h"
#include"cal.h"

void TIME2MJD(UTC* utc, MJD* mjd)//常规时间转化至简化儒略日
{
	double y, m;
	if (utc->Month > 2)
	{
		y = utc->Year;
		m = utc->Month;
	}
	else
	{
		y = utc->Year - 1;
		m = utc->Month + 12;
	}
	double UT = (double)utc->Hour + (double)utc->Min / 60 + (double)utc->Sec / 3600;
	double JD = (int)(365.25 * y) + (int)(30.6001 * (m + 1)) + utc->Day + UT / 24 + 1720981.5;
	mjd->FracDay = JD - 2400000.5;
	mjd->Days = (int)mjd->FracDay;
	mjd->FracDay = mjd->FracDay - mjd->Days;
}

void MJD2TIME(UTC* utc, MJD* mjd)//简化儒略日转通用时
{
	double JD = (double)mjd->Days + mjd->FracDay + 2400000.5;
	int a = (int)(JD + 0.5);
	int b = a + 1537;
	int c = (int)((b - 122.1) / 365.25);
	int d = (int)(365.25 * c);
	int e = (int)((b - d) / 30.6001);
	double D = b - d - (int)(30.6001 * e) + JD + 0.5 - a;
	utc->Day = (int)D;
	double h = 24 * (D - utc->Day);
	utc->Hour = (int)(h);
	double m = 60 * (h - utc->Hour);
	utc->Min = (int)(m);
	double s = 60 * (m - utc->Min);
	utc->Sec = (int)(s);
	utc->Month = e - 1 - 12 * (int)(e / 14);
	utc->Year = c - 4715 - (int)((7 + utc->Month) / 10);
}

void MJD2GPSTIME(MJD* mjd, GPSTIME* gpst)
{
	int gpsw = int((mjd->Days + mjd->FracDay - 44244) / 7);
	double gpss = double((mjd->Days + mjd->FracDay - 44244 - gpsw * 7)) * 86400.0;
	gpst->Week = gpsw;
	gpst->SecOfWeek = gpss;
}

void GPSTIME2MJD(MJD* mjd, GPSTIME* gpst)
{
	double Mjd = 44244 + gpst->Week * 7 + gpst->SecOfWeek / 86400.0;
	mjd->Days = int(Mjd);
	mjd->FracDay = double(Mjd - mjd->Days);
}

void GPSTIME2BDSTIME(GPSTIME* gpst, GPSTIME* bdst)
{
	bdst->Week = gpst->Week - 1356;
	bdst->SecOfWeek = gpst->SecOfWeek - 14;
}

void GPSTIME2BDSTIME(GPSTIME* gpst, BDSTIME* bdst)
{
	bdst->Week = gpst->Week - 1356;
	bdst->SecOfWeek = gpst->SecOfWeek - 14;
}

void BDSTIME2GPSTIME(GPSTIME* gpst, BDSTIME* bdst)
{
	gpst->Week = bdst->Week + 1356;
	gpst->SecOfWeek = bdst->SecOfWeek + 14;
}

void GPSTIME2TIME(GPSTIME* gpst, UTC* utc)
{
	UTC* u = new UTC();
	MJD* mjd = new MJD();
	GPSTIME2MJD(mjd, gpst);
	MJD2TIME(u, mjd);
	utc->Year = u->Year;
	utc->Month = u->Month;
	int day = gpst->SecOfWeek / 86400;
	utc->Day = u->Day + day;
	int hour = (gpst->SecOfWeek - 86400 * day) / 3600;
	utc->Hour = hour;
	int min = (gpst->SecOfWeek - 86400 * day - 3600 * hour) / 60;
	utc->Min = min;
	utc->Sec = gpst->SecOfWeek - 86400 * day - 3600 * hour - 60 * min;
}

void BLH2XYZ(BLH* blh, XYZ* xyz, double e2, double a)
{
	double B = blh->Lat * Pi / 180;
	double L = blh->Lon * Pi / 180;
	double N = a / sqrt(1 - e2 * sin(B) * sin(B));
	xyz->X = (blh->Height + N) * cos(B) * cos(L);
	xyz->Y = (blh->Height + N) * cos(B) * sin(L);
	xyz->Z = (N * (1 - e2) + blh->Height) * sin(B);
}

void XYZ2BLH(BLH* blh, XYZ* xyz, double e2, double a)
{
	double L = 180 * atan2(xyz->Y, xyz->X) / Pi;
	while (L < 0 || L > 360)
	{
		if (L < 0)
		{
			L += 360;
		}
		if (L > 360)
		{
			L -= 360;
		}
	}
	blh->Lon = L;

	double B = atan(xyz->Z / sqrt(xyz->X * xyz->X + xyz->Y * xyz->Y));
	double B0 = 0;
	double n;
	int count = 0;
	do
	{
		B0 = B;
		n = a / sqrt(1 - e2 * sin(B0) * sin(B0));
		B = atan((xyz->Z + n * e2 * sin(B0)) / sqrt(xyz->X * xyz->X + xyz->Y * xyz->Y));
		count++;
	} while (abs(B - B0) > 1e-10 && count < 100);
	n = a / sqrt(1 - e2 * sin(B) * sin(B));

	blh->Lat = 180 * B / Pi;

	blh->Height = xyz->Z / sin(B) - n * (1 - e2);
}

//转化为当地坐标系
double XYZ2ENU(XYZ* xyz1, XYZ* xyz2, XYZ* enu, int sys)
{
	BLH blh;
	switch (sys)
	{
	case SYS_GPS:
		XYZ2BLH(&blh, xyz1, WGS84_e2, WGS84_a);
		break;
	case SYS_BDS:
		XYZ2BLH(&blh, xyz1, CGCS2000_e2, CGCS2000_a);
		break;
	default:
		break;
	}
	MatrixXd R = MatrixXd::Zero(3, 3);
	double B = degree2rad(blh.Lat);
	double L = degree2rad(blh.Lon);
	R(0, 0) = -sin(L);
	R(0, 1) = cos(L);
	R(0, 2) = 0;
	R(1, 0) = -sin(B) * cos(L);
	R(1, 1) = -sin(B) * sin(L);
	R(1, 2) = cos(B);
	R(2, 0) = cos(B) * cos(L);
	R(2, 1) = cos(B) * sin(L);
	R(2, 2) = sin(B);
	MatrixXd x0(3, 1);
	MatrixXd x_local(3, 1);
	x0(0, 0) = xyz2->X - xyz1->X;
	x0(1, 0) = xyz2->Y - xyz1->Y;
	x0(2, 0) = xyz2->Z - xyz1->Z;
	x_local = R * x0;
	enu->X = x_local(0, 0);
	enu->Y = x_local(1, 0);
	enu->Z = x_local(2, 0);
	return x_local(2, 0);
}