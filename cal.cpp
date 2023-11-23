#include "cal.h"
#include "data.h"

/*卫星粗差探测存储变量*/
double GPS_GF[GPS_SAT_QUAN];
double GPS_MW[GPS_SAT_QUAN];
double GPS_PSE[6][GPS_SAT_QUAN];
double GPS_PHA[6][GPS_SAT_QUAN];
double GPS_DOP[6][GPS_SAT_QUAN];
int GPS_COUNT[GPS_SAT_QUAN];

double BDS_GF[BDS_SAT_QUAN];
double BDS_MW[BDS_SAT_QUAN];
double BDS_PSE[5][BDS_SAT_QUAN];
double BDS_PHA[5][BDS_SAT_QUAN];
double BDS_DOP[5][BDS_SAT_QUAN];
int BDS_COUNT[BDS_SAT_QUAN];

/*单历元星历数据存储*/
EPHEMERIS *GPS_eph[GPS_SAT_QUAN];
EPHEMERIS *BDS_eph[BDS_SAT_QUAN];

/*频点数、卫星种数配置*/
int phase_num; // 单频or双频
int SYS_num;   // 单星or双星
int User_SYS;  // 单星下指定系统，默认GPS
int Hop_used;  // 是否启用对流层改正

unsigned int initial()
{
	for (int i = 0; i < GPS_SAT_QUAN; i++)
	{
		GPS_GF[i] = 0;
		GPS_MW[i] = 0;
		GPS_eph[i] = new EPHEMERIS();
	}
	for (int i = 0; i < BDS_SAT_QUAN; i++)
	{
		BDS_GF[i] = 0;
		BDS_MW[i] = 0;
		BDS_eph[i] = new EPHEMERIS();
	}
	phase_num = 0;
	SYS_num = 0;
	User_SYS = SYS_GPS;
	printf("请选择单频或双频\n1. 单频\t2. 双频\n");
	cin >> phase_num;
	printf("请选择单系统或双系统\n1. 单系统\t2. 双系统\n");
	cin >> SYS_num;
	if (SYS_num == 1)
	{
		printf("请选择解算系统\n1. GPS\t2. BDS\n");
		cin >> User_SYS;
		switch (User_SYS)
		{
		case 1:
			User_SYS = SYS_GPS;
			break;
		case 2:
			User_SYS = SYS_BDS;
		default:
			break;
		}
	}
	printf("是否改正对流层\n1. 是\t2. 否\n");
	cin >> Hop_used;
	return 1;
}

double SQR(double x)
{
	return x * x;
}

double Len(XYZ *pos)
{
	return sqrt(SQR(pos->X) + SQR(pos->Y) + SQR(pos->Z));
}

double degree2rad(double degree)
{
	while (degree > 360 || degree < 0)
	{
		if (degree > 360)
			degree -= 360;
		if (degree < 0)
			degree += 360;
	}
	return degree * Pi / 180;
}

double rad2degree(double rad)
{
	return rad * 180 / Pi;
}

double Cal_PDOP(MatrixXd Qxx)
{
	if (Qxx.rows() < 3 || Qxx.cols() < 3)
		return 0;
	return sqrt(Qxx(0, 0) + Qxx(1, 1) + Qxx(2, 2));
}

unsigned int Cal_LEAST_SQR(MatrixXd B, MatrixXd l, MatrixXd P, MatrixXd *Qxx, MatrixXd &x, double *thegma, double *DOP)
{
	if (B.rows() != P.rows() || B.rows() != l.rows())
		return 0;
	if (B.rows() < B.cols())
		return 0;
	*Qxx = (B.transpose() * P * B).inverse();
	x = *Qxx * B.transpose() * P * l;
	MatrixXd v = B * x - l;
	int m = B.rows() - B.cols();
	*thegma = sqrt(((v.transpose() * P * v) / m)(0, 0));
	*DOP = Cal_PDOP(*Qxx);
	return 1;
}

unsigned int decode_SYN(int sys, int signal)
{
	switch (sys)
	{
	case SYS_GPS:
		switch (signal)
		{
		case 0:
			return CODE_L1C; /* L1C/A */
		case 5:
			return CODE_L2P; /* L2P    (OEM7) */
		case 9:
			return CODE_L2W; /* L2P(Y),semi-codeless */
		case 14:
			return CODE_L5Q; /* L5Q    (OEM6) */
		case 16:
			return CODE_L1L; /* L1C(P) (OEM7) */
		case 17:
			return CODE_L2S; /* L2C(M) (OEM7) */
		default:
			return UNKOWN;
			break;
		}
	case SYS_BDS:
		switch (signal)
		{
		case 0:
			return CODE_L2I; /* B1I with D1 (OEM6) */
		case 1:
			return CODE_L7I; /* B2I with D1 (OEM6) */
		case 2:
			return CODE_L6I; /* B3I with D1 (OEM7) */
		case 4:
			return CODE_L2I; /* B1I with D2 (OEM6) */
		case 5:
			return CODE_L7I; /* B2I with D2 (OEM6) */
		case 6:
			return CODE_L6I; /* B3I with D2 (OEM7) */
		case 7:
			return CODE_L1P; /* B1C(P) (OEM7) */
		case 9:
			return CODE_L5P; /* B2a(P) (OEM7) */
		default:
			return UNKOWN;
			break;
		}
	default:
		return UNKOWN;
		break;
	}
}

double CODE2FREQ(int code)
{
	switch (code)
	{
	case 0:
		return 0;
	case 1:
		return L1;
	case 2:
		return L2;
	case 3:
		return L2;
	case 4:
		return L5;
	case 5:
		return L1;
	case 6:
		return L2;
	case 7:
		return B1;
	case 8:
		return B2;
	case 9:
		return B3;
	case 10:
		return B1_C;
	case 11:
		return B2_a;
	default:
		return 0;
		break;
	}
}

// 钟差改正
double CORRECT_CLK(double t, EPHEMERIS *eph)
{
	double correct_clk = t;
	for (int i = 0; i < 10; i++)
	{
		correct_clk = t - (eph->a_f0 + eph->a_f1 * (correct_clk - eph->toc) + eph->a_f2 * (correct_clk - eph->toc) * (correct_clk - eph->toc));
	}
	return eph->a_f0 + eph->a_f1 * (correct_clk - eph->toc) + eph->a_f2 * (correct_clk - eph->toc) * (correct_clk - eph->toc);
}

double TGD(EPHEMERIS *e, double f, int sys)
{
	switch (sys)
	{
	case SYS_GPS:
		return SQR(f / L1) * e->T_GD1;
		break;
	case SYS_BDS:
		if (f == B1)
			return e->T_GD1;
		if (f == B2)
			return e->T_GD2;
		if (f == B3)
			return 0.0;
	default:
		return -1;
		break;
	}
}

// 星历位置
unsigned int SAT_POS_CAL(double t, EPHEMERIS *eph, XYZ *Sat_Pos, double &clk, double dt, int SYS)
{
	double n, delt_t, M, E, E0, V, u_, u, r, i, dt0, F;
	switch (SYS)
	{
	case SYS_GPS:
		n = sqrt(WGS84_GM) / (eph->sqrt_A * eph->sqrt_A * eph->sqrt_A) + eph->delt_n;
		dt0 = delt_t = t - eph->toe_tow;
		F = -2 * sqrt(WGS84_GM) / (velocity_c * velocity_c);
		break;
	case SYS_BDS:
		n = sqrt(CGCS2000_GM) / (eph->sqrt_A * eph->sqrt_A * eph->sqrt_A) + eph->delt_n;
		dt0 = delt_t = t - eph->toe_tow;
		F = -2 * sqrt(CGCS2000_GM) / (velocity_c * velocity_c);
		break;
	default:
		break;
	}

	double dtr = 0;
	while (abs(delt_t) > 302400)
	{
		if (delt_t > 302400)
		{
			delt_t -= 604800;
		}
		else if (delt_t < -302400)
		{
			delt_t += 604800;
		}
	}

	for (int i = 0; i < 10; i++)
	{
		M = eph->M0 + n * delt_t;
		E = M;
		E0 = M;
		int count = 0;
		do
		{
			E0 = E;
			E = M + eph->e * sin(E0);
			count++;
		} while (abs(E - E0) > 0.000000000001 && count < 10);
		dtr = F * eph->e * eph->sqrt_A * sin(E);
		delt_t = dt0 - dtr;
	}
	clk += dtr;
	dt += (dtr + clk);
	V = atan2((sqrt(1 - eph->e * eph->e) * sin(E)), (cos(E) - eph->e));
	u_ = eph->omiga + V;
	u = u_ + eph->Cuc * cos(2 * u_) + eph->Cus * sin(2 * u_);
	r = eph->sqrt_A * eph->sqrt_A * (1 - eph->e * cos(E)) + eph->Crc * cos(2 * u_) + eph->Crs * sin(2 * u_);
	i = eph->i0 + eph->Cic * cos(2 * u_) + eph->Cis * sin(2 * u_) + eph->dot_i * delt_t;
	double x = r * cos(u);
	double y = r * sin(u);
	double z = 0;
	double L;
	double x0, y0, z0;
	switch (SYS)
	{
	case SYS_GPS:
		L = eph->Omiga0 + (eph->dot_Omiga - omiga_earth) * t - eph->dot_Omiga * eph->toe_tow;
		Sat_Pos->X = (x * cos(L) - y * cos(i) * sin(L)) + omiga_earth * dt * (x * sin(L) + y * cos(i) * cos(L));
		Sat_Pos->Y = x * sin(L) + y * cos(i) * cos(L) - omiga_earth * dt * (x * cos(L) - y * cos(i) * sin(L));
		Sat_Pos->Z = y * sin(i);
		break;
	case SYS_BDS:
		if (fabs(eph->i0 - 0.0873) < 0.1 && fabs(eph->sqrt_A - 6493) < 1) // 通过轨道倾角和轨道根数判断是否为GEO卫星  i: 5/deg sqrt_A: 6493/sqrt_meter
		{
			L = eph->Omiga0 + eph->dot_Omiga * delt_t - omiga_earth * eph->toe_tow;
			x0 = x * cos(L) - y * cos(i) * sin(L);
			y0 = x * sin(L) + y * cos(i) * cos(L);
			z0 = y * sin(i);
			MatrixXd P_GK(3, 1);
			MatrixXd R_Z(3, 3);
			MatrixXd R_X(3, 3);
			MatrixXd P(3, 1);
			P_GK << x0,
				y0,
				z0;
			R_X << 1, 0, 0,
				0, cos(-5 * Pi / 180), sin(-5 * Pi / 180),
				0, -sin(-5 * Pi / 180), cos(-5 * Pi / 180);
			R_Z << cos(omiga_earth * delt_t), sin(omiga_earth * delt_t), 0,
				-sin(omiga_earth * delt_t), cos(omiga_earth * delt_t), 0,
				0, 0, 1;
			P = R_Z * R_X * P_GK;
			Sat_Pos->X = cos(omiga_earth * dt) * P(0, 0) + sin(omiga_earth * dt) * P(1, 0);
			Sat_Pos->Y = cos(omiga_earth * dt) * P(1, 0) - sin(omiga_earth * dt) * P(0, 0);
			Sat_Pos->Z = P(2, 0);
		}
		else
		{
			L = eph->Omiga0 + (eph->dot_Omiga - omiga_earth) * t - eph->dot_Omiga * eph->toe_tow;
			Sat_Pos->X = cos(omiga_earth * dt) * (x * cos(L) - y * cos(i) * sin(L)) + sin(omiga_earth * dt) * (x * sin(L) + y * cos(i) * cos(L));
			Sat_Pos->Y = cos(omiga_earth * dt) * (x * sin(L) + y * cos(i) * cos(L)) - sin(omiga_earth * dt) * (x * cos(L) - y * cos(i) * sin(L));
			Sat_Pos->Z = y * sin(i);
		}
		break;
	default:
		break;
	}

	return 1;
}

// 卫星高度角计算
double Ele_Angle(XYZ *SatPos, XYZ *RcvPos, int sys)
{
	XYZ Satenu;
	XYZ2ENU(RcvPos, SatPos, &Satenu, sys);

	return asin(Satenu.Z / Len(&Satenu));
}

// Hopefiled对流层改正(m)
double Hopefield(double E, double H)
{
	if (Hop_used == 2)
		return 0;
	double Ts = T0 - 0.0065 * (H - H0);
	double hd = 40136 + 148.72 * (Ts - 273.16);
	double Ps = P0 * pow(1 - 0.000026 * (H - H0), 5.225);
	double RH = RH0 * exp(-0.0006396 * (H - H0));
	double es = RH * exp(-37.2465 + 0.213166 * Ts - 0.000256908 * Ts * Ts);
	double md = sin(degree2rad(sqrt(E * E + 6.25)));
	double mw = sin(degree2rad(sqrt(E * E + 2.25)));
	double hw = 11000;
	double ZHD = 155.2 * 1e-7 * Ps * (hd - H) / Ts;
	double ZWD = 155.2 * 1e-7 * 4810 * es * (hw - H) / (Ts * Ts);
	return ZHD / md + ZWD / mw;
}

double Hopefield(XYZ *SatPos, XYZ *RcvPos, int sys)
{
	if (Hop_used == 2)
		return 0;
	if (SatPos->X == 0 && SatPos->Y == 0 && SatPos->Z == 0)
	{
		return 0;
	}
	double E = rad2degree(Ele_Angle(SatPos, RcvPos, sys));
	BLH *Rcvblh = new BLH();
	double H = 0;
	switch (sys)
	{
	case SYS_GPS:
		XYZ2BLH(Rcvblh, RcvPos, WGS84_e2, WGS84_a);
		H = Rcvblh->Height;
		break;
	case SYS_BDS:
		XYZ2BLH(Rcvblh, RcvPos, CGCS2000_e2, CGCS2000_a);
		H = Rcvblh->Height;
		break;
	default:
		break;
	}
	if (H < 20e3 && H > -100)
	{
		double Ts = T0 - 0.0065 * (H - H0);
		double hd = 40136 + 148.72 * (T0 - 273.16);
		double Ps = P0 * pow(1 - 0.000026 * (H - H0), 5.225);
		double RH = RH0 * exp(-0.0006396 * (H - H0));
		double es = RH * exp(-37.2465 + 0.213166 * Ts - 0.000256908 * Ts * Ts);
		double md = sin(degree2rad(sqrt(E * E + 6.25)));
		double mw = sin(degree2rad(sqrt(E * E + 2.25)));
		double hw = 11000;
		double ZHD = 155.2 * 1e-7 * Ps * (hd - H) / Ts;
		double ZWD = 155.2 * 1e-7 * 4810 * es * (hw - H) / (Ts * Ts);
		delete Rcvblh;
		return ZHD / md + ZWD / mw;
	}
	else
	{
		delete Rcvblh;
		return 0;
	}
}

double Klobuchar(XYZ *RcvPos, double E, double A, double alpha[4], double beta[4], double UT, double code, int sys)
{
	if (!(alpha[0] * alpha[1] * alpha[2] * alpha[3] * beta[0] * beta[1] * beta[2] * beta[3]))
		return -1;
	BLH RcvBLH;
	double T_g = 0;
	double EA = 0;
	double B_IPP = 0;
	double L_IPP = 0;
	double B_m = 0;
	double t = 0;
	double A_I = 0;
	double P_I = 0;
	double Phase_I = 0;
	double F = 0;
	switch (sys)
	{
	case SYS_GPS:
		XYZ2BLH(&RcvBLH, RcvPos, WGS84_e2, WGS84_a);
		EA = 0.0137 / (E + 0.11) - 0.022;
		B_IPP = RcvBLH.Lat + EA * cos(A);
		if (B_IPP < -0.416)
			B_IPP = -0.416;
		if (B_IPP > 0.416)
			B_IPP = 0.416;
		L_IPP = RcvBLH.Lon + EA * sin(A) / cos(B_IPP);
		B_m = B_IPP + 0.064 * cos(L_IPP - 1.617);
		t = 43200 * L_IPP + UT;
		while (t > 86400 || t < 0)
		{
			if (t > 86400)
				t -= 86400;
			if (t < 0)
				t += 86400;
		}
		A_I = 0;
		P_I = 0;
		for (int i = 0; i < 4; i++)
		{
			A_I += alpha[i] * pow(B_m, i);
			P_I += beta[i] * pow(B_m, i);
		}
		if (A_I < 0)
			A_I = 0;
		if (P_I < 72000)
			P_I = 72000;
		Phase_I = 2 * Pi * (t - 50400) / P_I;
		F = 1 + 16 * pow((0.53 - E), 3);
		if (abs(Phase_I) < 1.57)
		{
			T_g = F * (5e-9 + A_I * (1 - pow(Phase_I, 2) / 2 + pow(Phase_I, 4) / 24));
		}
		else
		{
			T_g = F * 5e-9;
		}
		return pow(L1 / code, 2) * T_g;
		break;
	case SYS_BDS:
		XYZ2BLH(&RcvBLH, RcvPos, CGCS2000_e2, CGCS2000_a);
		EA = Pi / 2 - E - asin(cos(E) * 6378 / (6378 + 375));
		B_IPP = asin(sin(RcvBLH.Lat) * cos(EA) + cos(RcvBLH.Lat) * sin(EA) * cos(A));
		L_IPP = RcvBLH.Lon + asin(sin(EA) * sin(A) / cos(B_IPP));
		t = UT + L_IPP * 43200 / Pi;
		A_I = 0;
		P_I = 0;
		for (int i = 0; i < 4; i++)
		{
			A_I += alpha[i] * pow(B_IPP / Pi, i);
			P_I += beta[i] * pow(B_IPP / Pi, i);
		}
		if (A_I < 0)
			A_I = 0;
		if (P_I < 72000)
			P_I = 72000;
		if (P_I > 172800)
			P_I = 172800;
		if (abs(t - 50400) < P_I / 4)
		{
			T_g = 5e-9 + A_I * cos(2 * Pi * (t - 50400) / P_I);
		}
		else
		{
			5e-9;
		}
		return T_g / sqrt(1 - pow(cos(E) * 6378 / (6378 + 375), 2));
		break;
	default:
		return 0;
		break;
	}
}

int CheckOBSConsist(Satellate *sate, int sys, double t, int index, bool &PSE_flag, bool &PHA_flag)
{
	int prn = sate->PRN;
	int Index = decode_SYN(sys, sate->SYG_TYPE[index]);
	double f = CODE2FREQ(Index);
	double lamda = 1e-6 * velocity_c / f;
	double C_D = 0;
	double mean_D = 0;
	double C_P = 0;
	double C_L = 0;
	if (sate->DOPPLER[0] == 0)
	{
		cout << "DOPPLER DATA LOSS!" << endl;
		return 0;
	}
	switch (sys)
	{
	case SYS_GPS:
		Index -= 1;
		if (GPS_PHA[Index][prn - 1] == 0 && GPS_PSE[Index][prn - 1] == 0 && GPS_DOP[Index][prn - 1] == 0)
		{
			GPS_PHA[Index][prn - 1] = sate->PHASE[index];
			GPS_PSE[Index][prn - 1] = sate->PSERA[index];
			GPS_DOP[Index][prn - 1] = sate->DOPPLER[index];
			return 1;
		}
		C_D = lamda * abs(GPS_DOP[Index][prn - 1] - sate->DOPPLER[index]);
		if (C_D > 20)
			return -1;
		mean_D = (GPS_DOP[Index][prn - 1] + sate->DOPPLER[index]) / 2;
		C_P = abs((sate->PSERA[index] - GPS_PSE[Index][prn - 1]) + lamda * mean_D * t);
		C_L = abs(lamda * (sate->PHASE[index] - GPS_PHA[Index][prn - 1]) + lamda * mean_D * t);
		if (C_P > 8)
			PSE_flag = false;
		if (C_L > 0.5)
			PHA_flag = false;
		GPS_PSE[Index][prn - 1] = sate->PSERA[index];
		GPS_PHA[Index][prn - 1] = sate->PHASE[index];
		GPS_DOP[Index][prn - 1] = sate->DOPPLER[index];

		return 1;
	case SYS_BDS:
		Index -= 7;
		if (BDS_PHA[Index][prn - 1] == 0 && BDS_PSE[Index][prn - 1] == 0 && BDS_DOP[Index][prn - 1] == 0)
		{
			BDS_PHA[Index][prn - 1] = sate->PHASE[index];
			BDS_PSE[Index][prn - 1] = sate->PSERA[index];
			BDS_DOP[Index][prn - 1] = sate->DOPPLER[index];
			return 1;
		}
		C_D = lamda * abs(BDS_DOP[Index][prn - 1] - sate->DOPPLER[index]);
		if (C_D > 20)
			return -1;
		mean_D = (BDS_DOP[Index][prn - 1] + sate->DOPPLER[index]) / 2;
		C_P = abs((sate->PSERA[index] - BDS_PSE[Index][prn - 1]) + lamda * mean_D * t);
		C_L = abs(lamda * (sate->PHASE[index] - BDS_PHA[Index][prn - 1]) + lamda * mean_D * t);
		if (C_P > 8)
			PSE_flag = false;
		if (C_L > 0.5)
			PHA_flag = false;
		BDS_PSE[Index][prn - 1] = sate->PSERA[index];
		BDS_PHA[Index][prn - 1] = sate->PHASE[index];
		BDS_DOP[Index][prn - 1] = sate->DOPPLER[index];

		return 1;
	default:
		return 0;
		break;
	}
}

int DetectOutlier(Satellate *sate, int sys, double t, int index1, int index2)
{
	int prn = sate->PRN;
	double f1 = CODE2FREQ(decode_SYN(sys, sate->SYG_TYPE[index1]));
	double f2 = CODE2FREQ(decode_SYN(sys, sate->SYG_TYPE[index2]));
	double lamda1 = 1e-6 * velocity_c / f1;
	double lamda2 = 1e-6 * velocity_c / f2;
	bool PSE_flag1 = true;
	bool PHA_flag1 = true;
	bool PSE_flag2 = true;
	bool PHA_flag2 = true;
	double GF = 0;
	double MW = 0;
	double dGF = 0;
	double dMW = 0;
	switch (sys)
	{
	case SYS_GPS:
		if (CheckOBSConsist(sate, sys, t, index1, PSE_flag1, PHA_flag1) && CheckOBSConsist(sate, sys, t, index2, PSE_flag2, PHA_flag2) && PSE_flag1 && PHA_flag1 && PSE_flag2 && PHA_flag2)
		{
			GF = sate->PSERA[index1] - sate->PSERA[index2];
			MW = (f1 - f2) * (sate->PSERA[index1] / lamda1 + sate->PSERA[index2] / lamda2) / (f1 + f2) - (sate->PHASE[index1] - sate->PHASE[index2]);
			if (GPS_GF[prn - 1] == 0 && GPS_MW[prn - 1] == 0)
			{
				GPS_GF[prn - 1] = GF;
				GPS_MW[prn - 1] = MW;
				GPS_COUNT[prn - 1]++;
				return 1;
			}
			dGF = abs(GF - GPS_GF[prn - 1]);
			dMW = abs(MW - GPS_MW[prn - 1]);
			GPS_GF[prn - 1] = GF;
			GPS_MW[prn - 1] = (GPS_MW[prn - 1] * (GPS_COUNT[prn - 1]++) + MW);
			GPS_MW[prn - 1] /= GPS_COUNT[prn - 1];
			if (dGF > GF_THRESH || dMW > MW_THRESH)
				return 0;

			return 1;
		}
		else
		{
			return 0;
		}
	case SYS_BDS:
		if (CheckOBSConsist(sate, sys, t, index1, PSE_flag1, PHA_flag1) && CheckOBSConsist(sate, sys, t, index2, PSE_flag2, PHA_flag2) && PSE_flag1 && PHA_flag1 && PSE_flag2 && PHA_flag2)
		{
			GF = sate->PSERA[index1] - sate->PSERA[index2];
			MW = (f1 - f2) * (sate->PSERA[index1] / lamda1 + sate->PSERA[index2] / lamda2) / (f1 + f2) - (sate->PHASE[index1] - sate->PHASE[index2]);
			if (BDS_GF[prn - 1] == 0 && BDS_MW[prn - 1] == 0)
			{
				BDS_GF[prn - 1] = GF;
				BDS_MW[prn - 1] = MW;
				BDS_COUNT[prn - 1]++;
				return 1;
			}
			dGF = abs(GF - BDS_GF[prn - 1]);
			dMW = abs(MW - BDS_MW[prn - 1]);
			BDS_GF[prn - 1] = GF;
			BDS_MW[prn - 1] = (BDS_MW[prn - 1] * (BDS_COUNT[prn - 1]++) + MW);
			BDS_MW[prn - 1] /= BDS_COUNT[prn - 1];
			if (dGF > GF_THRESH || dMW > MW_THRESH)
				return 0;
			return 1;
		}
		else
		{
			return 0;
		}
	default:
		return 0;
		break;
	}
}
#pragma region vector_EPOCH
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPOCH *eph, bool first_flag, double f1, double f2, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used)
{
	XYZ *sate_pos = new XYZ();
	double clk = 0;
	XYZ *RcvPos = new XYZ(Pos(0, 0), Pos(1, 0), Pos(2, 0));
	double dt_Rcv = Pos(3, 0);
	int ROWS = 0;
	MatrixXd x_Pos = MatrixXd::Zero(4, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	for (int i = 0; i < Sates.size(); i++)
	{
		if (Sates[i]->Phase_NUM < 2)
			continue;
		if (Sates[i]->Outlier)
			continue;
		int prn = Sates[i]->PRN;
		if (eph[prn - 1].num == 0)
			continue;
		double IF = 0;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index1])) != f1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index2])) != f2 && Index2 < MAXNUM)
			Index2++;
		if (Index1 == MAXNUM || Index2 == MAXNUM)
			continue;
		if (!(Sates[i]->LOCK_PSE[Index1] && Sates[i]->LOCK_PSE[Index2] && Sates[i]->LOCK_PHA[Index1] && Sates[i]->LOCK_PHA[Index2]))
			continue;

		// 计算卫星位置、钟差
		double ts = OBS_TIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;
		if (Sates[i]->SYS == SYS_BDS)
			ts -= 14;
		double dt = abs(ts - eph[prn - 1].epoch[0]->toe_tow);
		int index = 0;
		for (int j = 1; j < eph[prn - 1].num; j++)
		{
			if (abs(ts - eph[prn - 1].epoch[j]->toe_tow))
			{
				dt = abs(ts - eph[prn - 1].epoch[j]->toe_tow);
				index = j;
			}
		}

		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, eph[prn - 1].epoch[index]);
			SAT_POS_CAL(ts - clk, eph[prn - 1].epoch[index], sate_pos, clk, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		if (!first_flag)
		{
			if (Ele_Angle(sate_pos, RcvPos, Sates[i]->SYS) < degree2rad(10))
				continue;
		}
		if (Sates[i]->SYS == SYS_GPS)
			IF = SQR(f1) * Sates[i]->PSERA[Index1] / (SQR(f1) - SQR(f2)) - SQR(f2) * Sates[i]->PSERA[Index2] / (SQR(f1) - SQR(f2));
		else if (Sates[i]->SYS == SYS_BDS)
		{
			double k_1_3 = SQR(f1 / f2);
			IF = (Sates[i]->PSERA[Index2] - k_1_3 * Sates[i]->PSERA[Index1]) / (1 - k_1_3) + velocity_c * k_1_3 * eph[prn - 1].epoch[index]->T_GD1 / (1 - k_1_3);
		}

		double len = sqrt(SQR(RcvPos->X - sate_pos->X) + SQR(RcvPos->Y - sate_pos->Y) + SQR(RcvPos->Z - sate_pos->Z));
		double w_pos = IF - (len + dt_Rcv - velocity_c * clk + Hopefield(sate_pos, RcvPos, Sates[i]->SYS));
		if (!first_flag)
		{
			if (abs(w_pos) > 10)
				continue;
		}
		double l = (RcvPos->X - sate_pos->X) / len;
		double m = (RcvPos->Y - sate_pos->Y) / len;
		double n = (RcvPos->Z - sate_pos->Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		l_pos_new(0, 0) = w_pos;
		B_Pos->conservativeResize(B_Pos->rows() + 1, B_Pos->cols());
		B_Pos->bottomRows(1) = B_pos_new;
		l_Pos->conservativeResize(l_Pos->rows() + 1, l_Pos->cols());
		l_Pos->bottomRows(1) = l_pos_new;
		switch (Sates[i]->SYS)
		{
		case SYS_GPS:
			*sate_used += "G" + to_string(Sates[i]->PRN);
			break;
		case SYS_BDS:
			*sate_used += "C" + to_string(Sates[i]->PRN);
			break;
		default:
			break;
		}
		ROWS++;
	}
	*P_Pos = MatrixXd::Identity(ROWS, ROWS);
	delete sate_pos;
	delete RcvPos;
	return ROWS;
}

unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPOCH *eph, bool first_flag, double f, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used)
{
	XYZ *sate_pos = new XYZ();
	double clk = 0;
	XYZ *RcvPos = new XYZ(Pos(0, 0), Pos(1, 0), Pos(2, 0));
	double dt_Rcv = Pos(3, 0);
	int ROWS = 0;
	MatrixXd x_Pos = MatrixXd::Zero(4, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	for (int i = 0; i < Sates.size(); i++)
	{
		if (Sates[i]->Phase_NUM < 2)
			continue;
		if (Sates[i]->Outlier)
			continue;
		int prn = Sates[i]->PRN;
		if (eph[prn - 1].num == 0)
			continue;
		int Index = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index])) != f && Index < MAXNUM)
			Index++;
		if (Index == MAXNUM)
			continue;
		if (!(Sates[i]->LOCK_PSE[Index] && Sates[i]->LOCK_PHA[Index]))
			continue;

		// 计算卫星位置、钟差
		double ts = OBS_TIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;
		if (Sates[i]->SYS == SYS_BDS)
			ts -= 14;
		double dt = abs(ts - eph[prn - 1].epoch[0]->toe_tow);
		int index = 0;
		for (int j = 1; j < eph[prn - 1].num; j++)
		{
			if (abs(ts - eph[prn - 1].epoch[j]->toe_tow))
			{
				dt = abs(ts - eph[prn - 1].epoch[j]->toe_tow);
				index = j;
			}
		}
		double tgd = 0;
		tgd = TGD(eph[prn - 1].epoch[index], f, Sates[i]->SYS);

		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, eph[prn - 1].epoch[index]);
			SAT_POS_CAL(ts - clk + tgd, eph[prn - 1].epoch[index], sate_pos, clk, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		if (!first_flag)
		{
			if (Ele_Angle(sate_pos, RcvPos, Sates[i]->SYS) < degree2rad(10))
				continue;
		}

		double len = sqrt(SQR(RcvPos->X - sate_pos->X) + SQR(RcvPos->Y - sate_pos->Y) + SQR(RcvPos->Z - sate_pos->Z));
		double w_pos = Sates[i]->PSERA[Index] - (len + dt_Rcv - velocity_c * clk + Hopefield(sate_pos, RcvPos, Sates[i]->SYS));
		if (!first_flag)
		{
			if (abs(w_pos) > 10)
				continue;
		}
		double l = (RcvPos->X - sate_pos->X) / len;
		double m = (RcvPos->Y - sate_pos->Y) / len;
		double n = (RcvPos->Z - sate_pos->Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		l_pos_new(0, 0) = w_pos;
		B_Pos->conservativeResize(B_Pos->rows() + 1, B_Pos->cols());
		B_Pos->bottomRows(1) = B_pos_new;
		l_Pos->conservativeResize(l_Pos->rows() + 1, l_Pos->cols());
		l_Pos->bottomRows(1) = l_pos_new;
		switch (Sates[i]->SYS)
		{
		case SYS_GPS:
			*sate_used += "G" + to_string(Sates[i]->PRN);
			break;
		case SYS_BDS:
			*sate_used += "C" + to_string(Sates[i]->PRN);
			break;
		default:
			break;
		}
		ROWS++;
	}
	*P_Pos = MatrixXd::Identity(ROWS, ROWS);
	delete sate_pos;
	delete RcvPos;
	return ROWS;
}

unsigned int setup_Vel(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPOCH *eph, MatrixXd *B_Vel, MatrixXd *l_Vel, MatrixXd *P_Vel)
{
	XYZ *sate_pos0 = new XYZ();
	XYZ *sate_pos1 = new XYZ();
	double clk0 = 0;
	double clk1 = 0;
	double velocity[4] = {0, 0, 0, 0};
	XYZ *RcvPos = new XYZ(Pos(0, 0), Pos(1, 0), Pos(2, 0));
	MatrixXd x_vel = MatrixXd::Zero(4, 1);
	MatrixXd B_vel_new = MatrixXd::Zero(1, 4);
	MatrixXd l_vel_new = MatrixXd::Zero(1, 1);
	int ROWS = 0;
	for (int i = 0; i < Sates.size(); i++)
	{
		if (Sates[i]->Phase_NUM < 2)
			continue;
		if (Sates[i]->Outlier)
			continue;
		int prn = Sates[i]->PRN;
		if (eph[prn - 1].num == 0)
			continue;
		double ts = OBS_TIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;
		if (Sates[i]->SYS == SYS_BDS)
			ts -= 14;
		double dt = abs(ts - eph[prn - 1].epoch[0]->toe_tow);
		int index = 0;
		for (int j = 1; j < eph[prn - 1].num; j++)
		{
			if (abs(ts - eph[prn - 1].epoch[j]->toe_tow) < dt)
			{
				dt = abs(ts - eph[prn - 1].epoch[j]->toe_tow);
				index = j;
			}
		}

		for (int j = 0; j < 3; j++)
		{
			clk0 = CORRECT_CLK(ts - clk0, eph[prn - 1].epoch[index]);
			clk1 = CORRECT_CLK(ts + 1e-3 - clk1, eph[prn - 1].epoch[index]);
			SAT_POS_CAL(ts - clk0, eph[prn - 1].epoch[index], sate_pos0, clk0, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
			SAT_POS_CAL(ts + 1e-3 - clk1, eph[prn - 1].epoch[index], sate_pos1, clk1, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		velocity[0] = (sate_pos1->X - sate_pos0->X) / 1e-3;
		velocity[1] = (sate_pos1->Y - sate_pos0->Y) / 1e-3;
		velocity[2] = (sate_pos1->Z - sate_pos0->Z) / 1e-3;
		velocity[3] = (clk1 - clk0) * velocity_c / 1e-3;
		double f1 = CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[0]));
		double lamda = (1e-6 * velocity_c / f1);
		double len = sqrt(SQR(RcvPos->X - sate_pos0->X) + SQR(RcvPos->Y - sate_pos0->Y) + SQR(RcvPos->Z - sate_pos0->Z));
		double l = (RcvPos->X - sate_pos0->X) / len;
		double m = (RcvPos->Y - sate_pos0->Y) / len;
		double n = (RcvPos->Z - sate_pos0->Z) / len;
		double v0 = ((sate_pos0->X - RcvPos->X) * velocity[0] + (sate_pos0->Y - RcvPos->Y) * velocity[1] + (sate_pos0->Z - RcvPos->Z) * velocity[2]) / len;
		double w_Vel = -lamda * Sates[i]->DOPPLER[0] - (v0 - velocity[3]);
		if (abs(w_Vel) > 10)
			continue;
		B_vel_new(0, 0) = l;
		B_vel_new(0, 1) = m;
		B_vel_new(0, 2) = n;
		B_vel_new(0, 3) = 1;
		l_vel_new(0, 0) = w_Vel;
		B_Vel->conservativeResize(B_Vel->rows() + 1, B_Vel->cols());
		B_Vel->bottomRows(1) = B_vel_new;
		l_Vel->conservativeResize(l_Vel->rows() + 1, l_Vel->cols());
		l_Vel->bottomRows(1) = l_vel_new;
		ROWS++;
	}
	*P_Vel = MatrixXd::Identity(ROWS, ROWS);
	delete sate_pos0;
	delete sate_pos1;
	delete RcvPos;
	return ROWS;
}

unsigned int Cal_1(Result_DATA *data, OBS_DATA *obs, EPOCH *eph, bool first_flag)
{
	MatrixXd B_Pos;
	MatrixXd l_Pos;
	MatrixXd P_Pos;
	MatrixXd x_Pos;
	MatrixXd B_Vel;
	MatrixXd l_Vel;
	MatrixXd P_Vel;
	MatrixXd x_Vel;
	double *thegma_Pos = new double;
	double *PDOP = new double;
	double *thegma_Vel = new double;
	double *VDOP = new double;
	int ROWS = 0;
	for (int i = 0; i < 5; i++)
	{
		MatrixXd *temp_B = new MatrixXd(0, 4);
		MatrixXd *temp_l = new MatrixXd(0, 1);
		MatrixXd *temp_P = new MatrixXd();
		MatrixXd Rcvpos(4, 1);
		Rcvpos.block(0, 0, 3, 1) = data->Pos->block<3, 1>(0, 0);
		*data->SATES = "";
		switch (User_SYS)
		{
		case SYS_GPS:
			Rcvpos(3, 0) = (*data->Pos)(3, 0);
			if (phase_num == 1)
				ROWS = setup_Pos(obs->OBS_TIME, Rcvpos, obs->GPS_SATE, eph, first_flag, L1, temp_B, temp_l, temp_P, data->SATES);
			else if (phase_num == 2)
				ROWS = setup_Pos(obs->OBS_TIME, Rcvpos, obs->GPS_SATE, eph, first_flag, L1, L2, temp_B, temp_l, temp_P, data->SATES);
			break;
		case SYS_BDS:
			Rcvpos(3, 0) = (*data->Pos)(4, 0);
			if (phase_num == 1)
				ROWS = setup_Pos(obs->OBS_TIME, Rcvpos, obs->BDS_SATE, eph, first_flag, B3, temp_B, temp_l, temp_P, data->SATES);
			else if (phase_num == 2)
				ROWS = setup_Pos(obs->OBS_TIME, Rcvpos, obs->BDS_SATE, eph, first_flag, B1, B3, temp_B, temp_l, temp_P, data->SATES);
		default:
			break;
		}
		B_Pos = *temp_B;
		l_Pos = *temp_l;
		P_Pos = *temp_P;
		if (ROWS == 0)
		{
			delete temp_B;
			delete temp_l;
			delete temp_P;
			return 0;
		}
		if (Cal_LEAST_SQR(B_Pos, l_Pos, P_Pos, data->Q_Pos, x_Pos, thegma_Pos, PDOP))
		{
			data->thegma_Pos = thegma_Pos;
			data->PDOP = PDOP;
			switch (User_SYS)
			{
			case SYS_GPS:
				(*data->Pos).block(0, 0, 4, 1) += x_Pos;
				data->GPS_num = ROWS;
				data->BDS_num = 0;
				break;
			case SYS_BDS:
				(*data->Pos).block(0, 0, 3, 1) += x_Pos.block(0, 0, 3, 1);
				(*data->Pos)(4, 0) += x_Pos(3, 0);
				data->BDS_num = ROWS;
				data->GPS_num = 0;
				break;
			default:
				break;
			}
		}
		else
		{
			data->solve_result = OBS_DATA_Loss;
			delete temp_B;
			delete temp_l;
			delete temp_P;
			return 0;
		}
		if (*thegma_Pos < 2)
		{
			first_flag = false;
		}
		delete temp_B;
		delete temp_l;
		delete temp_P;
	}

	MatrixXd *temp_B = new MatrixXd(0, 4);
	MatrixXd *temp_l = new MatrixXd(0, 1);
	MatrixXd *temp_P = new MatrixXd();
	ROWS = setup_Vel(obs->OBS_TIME, *data->Pos, obs->GPS_SATE, eph, temp_B, temp_l, temp_P);
	B_Vel = *temp_B;
	l_Vel = *temp_l;
	P_Vel = *temp_P;
	temp_B->conservativeResize(0, 4);
	temp_l->conservativeResize(0, 1);
	if (Cal_LEAST_SQR(B_Vel, l_Vel, P_Vel, data->Q_Vel, x_Vel, thegma_Vel, VDOP))
	{
		(*data->Vel) = x_Vel;
		data->thegma_Vel = thegma_Vel;
	}
	delete temp_B;
	delete temp_l;
	delete temp_P;
	return 1;
}

unsigned int Cal_2(Result_DATA *data, OBS_DATA *obs, EPOCH *gpseph, EPOCH *bdseph, bool first_flag)
{
	MatrixXd B_Pos;
	MatrixXd l_Pos;
	MatrixXd P_Pos;
	MatrixXd x_Pos;
	MatrixXd B_Vel;
	MatrixXd l_Vel;
	MatrixXd P_Vel;
	MatrixXd x_Vel;
	double *thegma_Pos = new double;
	double *PDOP = new double;
	double *thegma_Vel = new double;
	double *VDOP = new double;
	int GPS_ROWS = 0;
	int BDS_ROWS = 0;
	bool GPS_set, BDS_set;
	for (int i = 0; i < 5; i++)
	{
		GPS_set = BDS_set = true;
		MatrixXd *temp_B = new MatrixXd(0, 4);
		MatrixXd *temp_l = new MatrixXd(0, 1);
		MatrixXd *temp_P = new MatrixXd();
		MatrixXd Rcvpos(4, 1);
		Rcvpos = data->Pos->block<4, 1>(0, 0);
		*data->SATES = "";
		if (phase_num == 1)
			GPS_ROWS = setup_Pos(obs->OBS_TIME, Rcvpos, obs->GPS_SATE, gpseph, first_flag, L1, temp_B, temp_l, temp_P, data->SATES);
		else if (phase_num == 2)
			GPS_ROWS = setup_Pos(obs->OBS_TIME, Rcvpos, obs->GPS_SATE, gpseph, first_flag, L1, L2, temp_B, temp_l, temp_P, data->SATES);
		B_Pos = *temp_B;
		l_Pos = *temp_l;
		P_Pos = *temp_P;
		Rcvpos(3, 0) = (*data->Pos)(4, 0);
		temp_B->conservativeResize(0, 4);
		temp_l->conservativeResize(0, 1);
		if (phase_num == 1)
			BDS_ROWS = setup_Pos(obs->OBS_TIME, Rcvpos, obs->BDS_SATE, bdseph, first_flag, B3, temp_B, temp_l, temp_P, data->SATES);
		else if (phase_num == 2)
			BDS_ROWS = setup_Pos(obs->OBS_TIME, Rcvpos, obs->BDS_SATE, bdseph, first_flag, B1, B3, temp_B, temp_l, temp_P, data->SATES);
		if (GPS_ROWS == 0 || BDS_ROWS == 0)
		{
			delete temp_B;
			delete temp_l;
			delete temp_P;
			return 0;
		}
		if (GPS_ROWS != 0 && BDS_ROWS == 0)
		{
			BDS_set = false;
		}
		else if (GPS_ROWS != 0 && BDS_ROWS == 0)
		{
			GPS_set = false;
			B_Pos = *temp_B;
			l_Pos = *temp_l;
			P_Pos = *temp_P;
		}
		else
		{
			B_Pos.conservativeResize(GPS_ROWS + BDS_ROWS, 5);
			l_Pos.conservativeResize(GPS_ROWS + BDS_ROWS, 1);
			P_Pos.conservativeResize(GPS_ROWS + BDS_ROWS, GPS_ROWS + BDS_ROWS);
			B_Pos.block(0, 4, GPS_ROWS, 1) = MatrixXd::Zero(GPS_ROWS, 1);
			B_Pos.block(GPS_ROWS, 0, BDS_ROWS, 3) = temp_B->block(0, 0, BDS_ROWS, 3);
			B_Pos.block(GPS_ROWS, 3, BDS_ROWS, 1) = MatrixXd::Zero(BDS_ROWS, 1);
			B_Pos.block(GPS_ROWS, 4, BDS_ROWS, 1) = temp_B->block(0, 3, BDS_ROWS, 1);
			l_Pos.bottomRows(BDS_ROWS) = *temp_l;
			P_Pos.block(GPS_ROWS, GPS_ROWS, BDS_ROWS, BDS_ROWS) = *temp_P;
			P_Pos.block(0, GPS_ROWS, GPS_ROWS, BDS_ROWS) = MatrixXd::Zero(GPS_ROWS, BDS_ROWS);
			P_Pos.block(GPS_ROWS, 0, BDS_ROWS, GPS_ROWS) = MatrixXd::Zero(BDS_ROWS, GPS_ROWS);
		}
		if (Cal_LEAST_SQR(B_Pos, l_Pos, P_Pos, data->Q_Pos, x_Pos, thegma_Pos, PDOP))
		{
			if (GPS_set && BDS_set)
				(*data->Pos) += x_Pos;
			else if (GPS_set && !BDS_set)
				(*data->Pos).block(0, 0, 4, 1) += x_Pos;
			else if (!GPS_set && BDS_set)
			{
				(*data->Pos).block(0, 0, 3, 1) += x_Pos.block(0, 0, 3, 1);
				(*data->Pos)(4, 0) += x_Pos(3, 0);
			}
			data->thegma_Pos = thegma_Pos;
			data->PDOP = PDOP;
			data->GPS_num = GPS_ROWS;
			data->BDS_num = BDS_ROWS;
		}
		else
		{
			data->solve_result = OBS_DATA_Loss;
			delete temp_B;
			delete temp_l;
			delete temp_P;
			return 0;
		}
		if (*thegma_Pos < 2)
		{
			first_flag = false;
		}
		delete temp_B;
		delete temp_l;
		delete temp_P;
	}

	MatrixXd *temp_B = new MatrixXd(0, 4);
	MatrixXd *temp_l = new MatrixXd(0, 1);
	MatrixXd *temp_P = new MatrixXd();
	GPS_ROWS = setup_Vel(obs->OBS_TIME, *data->Pos, obs->GPS_SATE, gpseph, temp_B, temp_l, temp_P);
	B_Vel = *temp_B;
	l_Vel = *temp_l;
	P_Vel = *temp_P;
	temp_B->conservativeResize(0, 4);
	temp_l->conservativeResize(0, 1);
	BDS_ROWS = setup_Vel(obs->OBS_TIME, *data->Pos, obs->BDS_SATE, bdseph, temp_B, temp_l, temp_P);
	B_Vel.conservativeResize(GPS_ROWS + BDS_ROWS, 4);
	l_Vel.conservativeResize(GPS_ROWS + BDS_ROWS, 1);
	P_Vel.conservativeResize(GPS_ROWS + BDS_ROWS, GPS_ROWS + BDS_ROWS);
	B_Vel.bottomRows(BDS_ROWS) = *temp_B;
	l_Vel.bottomRows(BDS_ROWS) = *temp_l;
	P_Vel.block(GPS_ROWS, GPS_ROWS, BDS_ROWS, BDS_ROWS) = *temp_P;
	P_Vel.block(0, GPS_ROWS, GPS_ROWS, BDS_ROWS) = MatrixXd::Zero(GPS_ROWS, BDS_ROWS);
	P_Vel.block(GPS_ROWS, 0, BDS_ROWS, GPS_ROWS) = MatrixXd::Zero(BDS_ROWS, GPS_ROWS);
	if (Cal_LEAST_SQR(B_Vel, l_Vel, P_Vel, data->Q_Vel, x_Vel, thegma_Vel, VDOP))
	{
		(*data->Vel) = x_Vel;
		data->thegma_Vel = thegma_Vel;
	}
	delete temp_B;
	delete temp_l;
	delete temp_P;
	return 1;
}

unsigned int Cal_SPP(Result_DATA *data, OBS_DATA *obs, EPOCH *gpseph, EPOCH *bdseph, double dt_e, bool first_flag)
{
	data->OBSTIME = obs->OBS_TIME;
	if (obs->GPS_SATE.size() + obs->BDS_SATE.size() < 5)
	{
		cout << "DATA_LOSS!" << endl;
		return -1;
	}
	for (int i = 0; i < obs->GPS_SATE.size(); i++)
	{
		if (obs->GPS_SATE[i]->Phase_NUM < 2)
			continue;
		int prn = obs->GPS_SATE[i]->PRN;
		double IF = 0;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(obs->GPS_SATE[i]->SYS, obs->GPS_SATE[i]->SYG_TYPE[Index1])) != L1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(obs->GPS_SATE[i]->SYS, obs->GPS_SATE[i]->SYG_TYPE[Index2])) != L2 && Index2 < MAXNUM)
			Index2++;

		if (!DetectOutlier(obs->GPS_SATE[i], obs->GPS_SATE[i]->SYS, dt_e, Index1, Index2) || Index1 == MAXNUM || Index2 == MAXNUM)
		{
			obs->GPS_SATE[i]->Outlier = true;
			continue;
		}
	}
	for (int i = 0; i < obs->BDS_SATE.size(); i++)
	{
		if (obs->BDS_SATE[i]->Phase_NUM < 2)
			continue;
		int prn = obs->BDS_SATE[i]->PRN;
		double IF = 0;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(obs->BDS_SATE[i]->SYS, obs->BDS_SATE[i]->SYG_TYPE[Index1])) != B1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(obs->BDS_SATE[i]->SYS, obs->BDS_SATE[i]->SYG_TYPE[Index2])) != B3 && Index2 < MAXNUM)
			Index2++;

		if (!DetectOutlier(obs->BDS_SATE[i], obs->BDS_SATE[i]->SYS, dt_e, Index1, Index2) || Index1 == MAXNUM || Index2 == MAXNUM)
		{
			obs->BDS_SATE[i]->Outlier = true;
			continue;
		}
	}
	if (SYS_num == 1)
	{
		switch (User_SYS)
		{
		case SYS_GPS:
			if (Cal_1(data, obs, gpseph, first_flag))
				data->solve_result = Success;
			else
				return 0;
			break;
		case SYS_BDS:
			if (Cal_1(data, obs, bdseph, first_flag))
				data->solve_result = Success;
			else
				return 0;
		default:
			break;
		}
	}
	else if (SYS_num == 2)
	{
		if (Cal_2(data, obs, gpseph, bdseph, first_flag))
			data->solve_result = Success;
		else
			return 0;
	}

	return 1;
}
#pragma endregion

unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, bool first_flag, double f1, double f2, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used)
{
	XYZ *sate_pos = new XYZ();
	double clk = 0;
	XYZ *RcvPos = new XYZ(Pos(0, 0), Pos(1, 0), Pos(2, 0));
	double dt_Rcv = Pos(3, 0);
	int ROWS = 0;
	MatrixXd x_Pos = MatrixXd::Zero(4, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	for (int i = 0; i < Sates.size(); i++)
	{
		if (Sates[i]->Phase_NUM < 2)
			continue;
		if (Sates[i]->Outlier)
			continue;
		int prn = Sates[i]->PRN;
		// if (prn == 26)
		//	continue;
		if (eph[prn - 1]->PRN != prn)
			continue;
		double IF = 0;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index1])) != f1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index2])) != f2 && Index2 < MAXNUM)
			Index2++;
		if (Index1 == MAXNUM || Index2 == MAXNUM)
			continue;
		if (!(Sates[i]->LOCK_PSE[Index1] && Sates[i]->LOCK_PSE[Index2] && Sates[i]->LOCK_PHA[Index1] && Sates[i]->LOCK_PHA[Index2]))
			continue;

		// 计算卫星位置、钟差
		double ts = OBS_TIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;
		double dt = abs(ts - eph[prn - 1]->toe_tow + (OBS_TIME->Week - eph[prn - 1]->toe_wn) * 604800);
		if (dt > 14400)
			continue;

		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, eph[prn - 1]);
			SAT_POS_CAL(ts - clk, eph[prn - 1], sate_pos, clk, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		if (!first_flag)
		{
			if (Ele_Angle(sate_pos, RcvPos, Sates[i]->SYS) < degree2rad(10))
				continue;
		}
		if (Sates[i]->SYS == SYS_GPS)
			IF = SQR(f1) * Sates[i]->PSERA[Index1] / (SQR(f1) - SQR(f2)) - SQR(f2) * Sates[i]->PSERA[Index2] / (SQR(f1) - SQR(f2));
		else if (Sates[i]->SYS == SYS_BDS)
		{
			double k_1_3 = SQR(f1 / f2);
			IF = (Sates[i]->PSERA[Index2] - k_1_3 * Sates[i]->PSERA[Index1]) / (1 - k_1_3) + velocity_c * k_1_3 * eph[prn - 1]->T_GD1 / (1 - k_1_3);
		}

		double len = sqrt(SQR(RcvPos->X - sate_pos->X) + SQR(RcvPos->Y - sate_pos->Y) + SQR(RcvPos->Z - sate_pos->Z));
		double w_pos = IF - (len + dt_Rcv - velocity_c * clk + Hopefield(sate_pos, RcvPos, Sates[i]->SYS));
		if (!first_flag)
		{
			if (abs(w_pos) > 10)
				continue;
		}
		double l = (RcvPos->X - sate_pos->X) / len;
		double m = (RcvPos->Y - sate_pos->Y) / len;
		double n = (RcvPos->Z - sate_pos->Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		l_pos_new(0, 0) = w_pos;
		B_Pos->conservativeResize(B_Pos->rows() + 1, B_Pos->cols());
		B_Pos->bottomRows(1) = B_pos_new;
		l_Pos->conservativeResize(l_Pos->rows() + 1, l_Pos->cols());
		l_Pos->bottomRows(1) = l_pos_new;
		switch (Sates[i]->SYS)
		{
		case SYS_GPS:
			*sate_used += "G" + to_string(Sates[i]->PRN);
			break;
		case SYS_BDS:
			*sate_used += "C" + to_string(Sates[i]->PRN);
			break;
		default:
			break;
		}
		ROWS++;
	}
	*P_Pos = MatrixXd::Identity(ROWS, ROWS);
	delete sate_pos;
	delete RcvPos;
	return ROWS;
}

unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, bool first_flag, double f, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used)
{
	XYZ *sate_pos = new XYZ();
	double clk = 0;
	XYZ *RcvPos = new XYZ(Pos(0, 0), Pos(1, 0), Pos(2, 0));
	double dt_Rcv = Pos(3, 0);
	int ROWS = 0;
	MatrixXd x_Pos = MatrixXd::Zero(4, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	for (int i = 0; i < Sates.size(); i++)
	{
		if (Sates[i]->Phase_NUM < 2)
			continue;
		if (Sates[i]->Outlier)
			continue;
		int prn = Sates[i]->PRN;
		if (eph[prn - 1]->PRN != prn)
			continue;
		int Index = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index])) != f && Index < MAXNUM)
			Index++;
		if (Index == MAXNUM)
			continue;
		if (!(Sates[i]->LOCK_PSE[Index] && Sates[i]->LOCK_PHA[Index]))
			continue;

		// 计算卫星位置、钟差
		double ts = OBS_TIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;
		double dt = abs(ts - eph[prn - 1]->toe_tow + (OBS_TIME->Week - eph[prn - 1]->toe_wn) * 604800);
		if (dt > 14400)
			continue;

		double tgd = 0;
		tgd = TGD(eph[prn - 1], f, Sates[i]->SYS);

		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, eph[prn - 1]);
			SAT_POS_CAL(ts - clk + tgd, eph[prn - 1], sate_pos, clk, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		if (!first_flag)
		{
			if (Ele_Angle(sate_pos, RcvPos, Sates[i]->SYS) < degree2rad(10))
				continue;
		}

		double len = sqrt(SQR(RcvPos->X - sate_pos->X) + SQR(RcvPos->Y - sate_pos->Y) + SQR(RcvPos->Z - sate_pos->Z));
		double w_pos = Sates[i]->PSERA[Index] - (len + dt_Rcv - velocity_c * clk + Hopefield(sate_pos, RcvPos, Sates[i]->SYS));
		if (!first_flag)
		{
			if (abs(w_pos) > 10)
				continue;
		}
		double l = (RcvPos->X - sate_pos->X) / len;
		double m = (RcvPos->Y - sate_pos->Y) / len;
		double n = (RcvPos->Z - sate_pos->Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		l_pos_new(0, 0) = w_pos;
		B_Pos->conservativeResize(B_Pos->rows() + 1, B_Pos->cols());
		B_Pos->bottomRows(1) = B_pos_new;
		l_Pos->conservativeResize(l_Pos->rows() + 1, l_Pos->cols());
		l_Pos->bottomRows(1) = l_pos_new;
		switch (Sates[i]->SYS)
		{
		case SYS_GPS:
			*sate_used += "G" + to_string(Sates[i]->PRN);
			break;
		case SYS_BDS:
			*sate_used += "C" + to_string(Sates[i]->PRN);
			break;
		default:
			break;
		}
		ROWS++;
	}
	*P_Pos = MatrixXd::Identity(ROWS, ROWS);
	delete sate_pos;
	delete RcvPos;
	return ROWS;
}

unsigned int setup_Vel(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, MatrixXd *B_Vel, MatrixXd *l_Vel, MatrixXd *P_Vel)
{
	XYZ *sate_pos0 = new XYZ();
	XYZ *sate_pos1 = new XYZ();
	double clk0 = 0;
	double clk1 = 0;
	double velocity[4] = {0, 0, 0, 0};
	XYZ *RcvPos = new XYZ(Pos(0, 0), Pos(1, 0), Pos(2, 0));
	MatrixXd x_vel = MatrixXd::Zero(4, 1);
	MatrixXd B_vel_new = MatrixXd::Zero(1, 4);
	MatrixXd l_vel_new = MatrixXd::Zero(1, 1);
	int ROWS = 0;
	for (int i = 0; i < Sates.size(); i++)
	{
		if (Sates[i]->Phase_NUM < 2)
			continue;
		if (Sates[i]->Outlier)
			continue;
		int prn = Sates[i]->PRN;
		if (eph[prn - 1]->PRN != prn)
			continue;
		double ts = OBS_TIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;

		double dt = abs(ts - eph[prn - 1]->toe_tow + (OBS_TIME->Week - eph[prn - 1]->toe_wn) * 604800);
		if (dt > 14400)
			continue;

		for (int j = 0; j < 3; j++)
		{
			clk0 = CORRECT_CLK(ts - clk0, eph[prn - 1]);
			clk1 = CORRECT_CLK(ts + 1e-3 - clk1, eph[prn - 1]);
			SAT_POS_CAL(ts - clk0, eph[prn - 1], sate_pos0, clk0, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
			SAT_POS_CAL(ts + 1e-3 - clk1, eph[prn - 1], sate_pos1, clk1, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		velocity[0] = (sate_pos1->X - sate_pos0->X) / 1e-3;
		velocity[1] = (sate_pos1->Y - sate_pos0->Y) / 1e-3;
		velocity[2] = (sate_pos1->Z - sate_pos0->Z) / 1e-3;
		velocity[3] = (clk1 - clk0) * velocity_c / 1e-3;
		double f1 = CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[0]));
		double lamda = (1e-6 * velocity_c / f1);
		double len = sqrt(SQR(RcvPos->X - sate_pos0->X) + SQR(RcvPos->Y - sate_pos0->Y) + SQR(RcvPos->Z - sate_pos0->Z));
		double l = (RcvPos->X - sate_pos0->X) / len;
		double m = (RcvPos->Y - sate_pos0->Y) / len;
		double n = (RcvPos->Z - sate_pos0->Z) / len;
		double v0 = ((sate_pos0->X - RcvPos->X) * velocity[0] + (sate_pos0->Y - RcvPos->Y) * velocity[1] + (sate_pos0->Z - RcvPos->Z) * velocity[2]) / len;
		double w_Vel = -lamda * Sates[i]->DOPPLER[0] - (v0 - velocity[3]);
		if (abs(w_Vel) > 10)
			continue;
		B_vel_new(0, 0) = l;
		B_vel_new(0, 1) = m;
		B_vel_new(0, 2) = n;
		B_vel_new(0, 3) = 1;
		l_vel_new(0, 0) = w_Vel;
		B_Vel->conservativeResize(B_Vel->rows() + 1, B_Vel->cols());
		B_Vel->bottomRows(1) = B_vel_new;
		l_Vel->conservativeResize(l_Vel->rows() + 1, l_Vel->cols());
		l_Vel->bottomRows(1) = l_vel_new;
		ROWS++;
	}
	*P_Vel = MatrixXd::Identity(ROWS, ROWS);
	delete sate_pos0;
	delete sate_pos1;
	delete RcvPos;
	return ROWS;
}

unsigned int Cal_1(Result_DATA *data, OBS_DATA *obs, EPHEMERIS **eph, bool first_flag)
{
	MatrixXd B_Pos;
	MatrixXd l_Pos;
	MatrixXd P_Pos;
	MatrixXd x_Pos;
	MatrixXd B_Vel;
	MatrixXd l_Vel;
	MatrixXd P_Vel;
	MatrixXd x_Vel;
	int ROWS = 0;
	GPSTIME *OBS_TIME_bdst = new GPSTIME();
	GPSTIME2BDSTIME(obs->OBS_TIME, OBS_TIME_bdst);
	for (int i = 0; i < 5; i++)
	{
		MatrixXd *temp_B = new MatrixXd(0, 4);
		MatrixXd *temp_l = new MatrixXd(0, 1);
		MatrixXd *temp_P = new MatrixXd();
		MatrixXd Rcvpos(4, 1);
		Rcvpos = data->Pos->block<4, 1>(0, 0);
		*data->SATES = "";
		switch (User_SYS)
		{
		case SYS_GPS:
			if (phase_num == 1)
				ROWS = setup_Pos(obs->OBS_TIME, Rcvpos, obs->GPS_SATE, eph, first_flag, L1, temp_B, temp_l, temp_P, data->SATES);
			else if (phase_num == 2)
				ROWS = setup_Pos(obs->OBS_TIME, Rcvpos, obs->GPS_SATE, eph, first_flag, L1, L2, temp_B, temp_l, temp_P, data->SATES);
			break;
		case SYS_BDS:
			if (phase_num == 1)
				ROWS = setup_Pos(OBS_TIME_bdst, Rcvpos, obs->BDS_SATE, eph, first_flag, B3, temp_B, temp_l, temp_P, data->SATES);
			else if (phase_num == 2)
				ROWS = setup_Pos(OBS_TIME_bdst, Rcvpos, obs->BDS_SATE, eph, first_flag, B1, B3, temp_B, temp_l, temp_P, data->SATES);
		default:
			break;
		}
		B_Pos = *temp_B;
		l_Pos = *temp_l;
		P_Pos = *temp_P;
		if (ROWS == 0)
		{
			delete temp_B;
			delete temp_l;
			delete temp_P;
			return 0;
		}
		if (Cal_LEAST_SQR(B_Pos, l_Pos, P_Pos, data->Q_Pos, x_Pos, data->thegma_Pos, data->PDOP))
		{
			switch (User_SYS)
			{
			case SYS_GPS:
				(*data->Pos).block(0, 0, 4, 1) += x_Pos;
				data->GPS_num = ROWS;
				data->BDS_num = 0;
				break;
			case SYS_BDS:
				(*data->Pos).block(0, 0, 3, 1) += x_Pos.block(0, 0, 3, 1);
				(*data->Pos)(4, 0) += x_Pos(3, 0);
				data->BDS_num = ROWS;
				data->GPS_num = 0;
				break;
			default:
				break;
			}
		}
		else
		{
			data->solve_result = OBS_DATA_Loss;
			delete temp_B;
			delete temp_l;
			delete temp_P;
			return 0;
		}
		if (*data->thegma_Pos < 2)
		{
			first_flag = false;
		}
		delete temp_B;
		delete temp_l;
		delete temp_P;
	}

	MatrixXd *temp_B = new MatrixXd(0, 4);
	MatrixXd *temp_l = new MatrixXd(0, 1);
	MatrixXd *temp_P = new MatrixXd();
	ROWS = setup_Vel(obs->OBS_TIME, *data->Pos, obs->GPS_SATE, eph, temp_B, temp_l, temp_P);
	B_Vel = *temp_B;
	l_Vel = *temp_l;
	P_Vel = *temp_P;
	temp_B->conservativeResize(0, 4);
	temp_l->conservativeResize(0, 1);
	if (Cal_LEAST_SQR(B_Vel, l_Vel, P_Vel, data->Q_Vel, x_Vel, data->thegma_Vel, data->VDOP))
	{
		(*data->Vel) = x_Vel;
	}
	delete temp_B;
	delete temp_l;
	delete temp_P;
	delete OBS_TIME_bdst;
	return 1;
}

unsigned int Cal_2(Result_DATA *data, OBS_DATA *obs, EPHEMERIS **gpseph, EPHEMERIS **bdseph, bool first_flag)
{
	MatrixXd B_Pos;
	MatrixXd l_Pos;
	MatrixXd P_Pos;
	MatrixXd x_Pos;
	MatrixXd B_Vel;
	MatrixXd l_Vel;
	MatrixXd P_Vel;
	MatrixXd x_Vel;
	int GPS_ROWS = 0;
	int BDS_ROWS = 0;
	GPSTIME *OBS_TIME_bdst = new GPSTIME();
	GPSTIME2BDSTIME(obs->OBS_TIME, OBS_TIME_bdst);
	bool GPS_set, BDS_set;
	for (int i = 0; i < 5; i++)
	{
		GPS_set = BDS_set = true;
		MatrixXd *temp_B = new MatrixXd(0, 4);
		MatrixXd *temp_l = new MatrixXd(0, 1);
		MatrixXd *temp_P = new MatrixXd();
		MatrixXd Rcvpos(4, 1);
		Rcvpos = data->Pos->block<4, 1>(0, 0);
		*data->SATES = "";
		if (phase_num == 1)
			GPS_ROWS = setup_Pos(obs->OBS_TIME, Rcvpos, obs->GPS_SATE, gpseph, first_flag, L1, temp_B, temp_l, temp_P, data->SATES);
		else if (phase_num == 2)
			GPS_ROWS = setup_Pos(obs->OBS_TIME, Rcvpos, obs->GPS_SATE, gpseph, first_flag, L1, L2, temp_B, temp_l, temp_P, data->SATES);
		B_Pos = *temp_B;
		l_Pos = *temp_l;
		P_Pos = *temp_P;
		Rcvpos(3, 0) = (*data->Pos)(4, 0);
		temp_B->conservativeResize(0, 4);
		temp_l->conservativeResize(0, 1);
		if (phase_num == 1)
			BDS_ROWS = setup_Pos(OBS_TIME_bdst, Rcvpos, obs->BDS_SATE, bdseph, first_flag, B3, temp_B, temp_l, temp_P, data->SATES);
		else if (phase_num == 2)
			BDS_ROWS = setup_Pos(OBS_TIME_bdst, Rcvpos, obs->BDS_SATE, bdseph, first_flag, B1, B3, temp_B, temp_l, temp_P, data->SATES);
		if (GPS_ROWS == 0 || BDS_ROWS == 0)
		{
			delete temp_B;
			delete temp_l;
			delete temp_P;
			delete OBS_TIME_bdst;
			data->solve_result = Set_UP_B_fail;
			return 0;
		}
		if (GPS_ROWS != 0 && BDS_ROWS == 0)
		{
			BDS_set = false;
		}
		else if (GPS_ROWS != 0 && BDS_ROWS == 0)
		{
			GPS_set = false;
			B_Pos = *temp_B;
			l_Pos = *temp_l;
			P_Pos = *temp_P;
		}
		else
		{
			B_Pos.conservativeResize(GPS_ROWS + BDS_ROWS, 5);
			l_Pos.conservativeResize(GPS_ROWS + BDS_ROWS, 1);
			P_Pos.conservativeResize(GPS_ROWS + BDS_ROWS, GPS_ROWS + BDS_ROWS);
			B_Pos.block(0, 4, GPS_ROWS, 1) = MatrixXd::Zero(GPS_ROWS, 1);
			B_Pos.block(GPS_ROWS, 0, BDS_ROWS, 3) = temp_B->block(0, 0, BDS_ROWS, 3);
			B_Pos.block(GPS_ROWS, 3, BDS_ROWS, 1) = MatrixXd::Zero(BDS_ROWS, 1);
			B_Pos.block(GPS_ROWS, 4, BDS_ROWS, 1) = temp_B->block(0, 3, BDS_ROWS, 1);
			l_Pos.bottomRows(BDS_ROWS) = *temp_l;
			P_Pos.block(GPS_ROWS, GPS_ROWS, BDS_ROWS, BDS_ROWS) = *temp_P;
			P_Pos.block(0, GPS_ROWS, GPS_ROWS, BDS_ROWS) = MatrixXd::Zero(GPS_ROWS, BDS_ROWS);
			P_Pos.block(GPS_ROWS, 0, BDS_ROWS, GPS_ROWS) = MatrixXd::Zero(BDS_ROWS, GPS_ROWS);
		}
		if (Cal_LEAST_SQR(B_Pos, l_Pos, P_Pos, data->Q_Pos, x_Pos, data->thegma_Pos, data->PDOP))
		{
			if (GPS_set && BDS_set)
				(*data->Pos) += x_Pos;
			else if (GPS_set && !BDS_set)
				(*data->Pos).block(0, 0, 4, 1) += x_Pos;
			else if (!GPS_set && BDS_set)
			{
				(*data->Pos).block(0, 0, 3, 1) += x_Pos.block(0, 0, 3, 1);
				(*data->Pos)(4, 0) += x_Pos(3, 0);
			}
			data->GPS_num = GPS_ROWS;
			data->BDS_num = BDS_ROWS;
		}
		else
		{
			data->solve_result = OBS_DATA_Loss;
			delete temp_B;
			delete temp_l;
			delete temp_P;
			return 0;
		}
		if (*data->thegma_Pos < 2)
		{
			first_flag = false;
		}
		delete temp_B;
		delete temp_l;
		delete temp_P;
	}

	MatrixXd *temp_B = new MatrixXd(0, 4);
	MatrixXd *temp_l = new MatrixXd(0, 1);
	MatrixXd *temp_P = new MatrixXd();
	GPS_ROWS = setup_Vel(obs->OBS_TIME, *data->Pos, obs->GPS_SATE, gpseph, temp_B, temp_l, temp_P);
	B_Vel = *temp_B;
	l_Vel = *temp_l;
	P_Vel = *temp_P;
	temp_B->conservativeResize(0, 4);
	temp_l->conservativeResize(0, 1);
	BDS_ROWS = setup_Vel(obs->OBS_TIME, *data->Pos, obs->BDS_SATE, bdseph, temp_B, temp_l, temp_P);
	B_Vel.conservativeResize(GPS_ROWS + BDS_ROWS, 4);
	l_Vel.conservativeResize(GPS_ROWS + BDS_ROWS, 1);
	P_Vel.conservativeResize(GPS_ROWS + BDS_ROWS, GPS_ROWS + BDS_ROWS);
	B_Vel.bottomRows(BDS_ROWS) = *temp_B;
	l_Vel.bottomRows(BDS_ROWS) = *temp_l;
	P_Vel.block(GPS_ROWS, GPS_ROWS, BDS_ROWS, BDS_ROWS) = *temp_P;
	P_Vel.block(0, GPS_ROWS, GPS_ROWS, BDS_ROWS) = MatrixXd::Zero(GPS_ROWS, BDS_ROWS);
	P_Vel.block(GPS_ROWS, 0, BDS_ROWS, GPS_ROWS) = MatrixXd::Zero(BDS_ROWS, GPS_ROWS);
	if (Cal_LEAST_SQR(B_Vel, l_Vel, P_Vel, data->Q_Vel, x_Vel, data->thegma_Vel, data->VDOP))
	{
		(*data->Vel) = x_Vel;
	}
	delete temp_B;
	delete temp_l;
	delete temp_P;
	delete OBS_TIME_bdst;
	return 1;
}

// 位置解算
unsigned int Cal_SPP(Result_DATA *data, OBS_DATA *obs, EPHEMERIS **gpseph, EPHEMERIS **bdseph, double dt_e, bool first_flag)
{
	data->OBSTIME = obs->OBS_TIME;
	if (obs->GPS_SATE.size() + obs->BDS_SATE.size() < 5)
	{
		data->solve_result = OBS_DATA_Loss;
		return -1;
	}
	for (int i = 0; i < obs->GPS_SATE.size(); i++)
	{
		if (obs->GPS_SATE[i]->Phase_NUM < 2)
			continue;
		int prn = obs->GPS_SATE[i]->PRN;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(obs->GPS_SATE[i]->SYS, obs->GPS_SATE[i]->SYG_TYPE[Index1])) != L1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(obs->GPS_SATE[i]->SYS, obs->GPS_SATE[i]->SYG_TYPE[Index2])) != L2 && Index2 < MAXNUM)
			Index2++;

		if (!DetectOutlier(obs->GPS_SATE[i], obs->GPS_SATE[i]->SYS, dt_e, Index1, Index2) || Index1 == MAXNUM || Index2 == MAXNUM)
		{
			obs->GPS_SATE[i]->Outlier = true;
			continue;
		}
	}
	for (int i = 0; i < obs->BDS_SATE.size(); i++)
	{
		if (obs->BDS_SATE[i]->Phase_NUM < 2)
			continue;
		int prn = obs->BDS_SATE[i]->PRN;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(obs->BDS_SATE[i]->SYS, obs->BDS_SATE[i]->SYG_TYPE[Index1])) != B1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(obs->BDS_SATE[i]->SYS, obs->BDS_SATE[i]->SYG_TYPE[Index2])) != B3 && Index2 < MAXNUM)
			Index2++;

		if (!DetectOutlier(obs->BDS_SATE[i], obs->BDS_SATE[i]->SYS, dt_e, Index1, Index2) || Index1 == MAXNUM || Index2 == MAXNUM)
		{
			obs->BDS_SATE[i]->Outlier = true;
			continue;
		}
	}
	if (SYS_num == 1)
	{
		switch (User_SYS)
		{
		case SYS_GPS:
			if (Cal_1(data, obs, gpseph, first_flag))
				data->solve_result = Success;
			else
				return 0;
			break;
		case SYS_BDS:
			if (Cal_1(data, obs, bdseph, first_flag))
				data->solve_result = Success;
			else
				return 0;
		default:
			break;
		}
	}
	else if (SYS_num == 2)
	{
		if (Cal_2(data, obs, gpseph, bdseph, first_flag))
			data->solve_result = Success;
		else
			return 0;
	}

	return 1;
}

int decodestream(Result_DATA *result, unsigned char Buff[], int &d)
{
	unsigned char TempBuff[MAXRAWLEN];
	int len;
	int msgID, msgTYPE;
	GPSTIME *gpstime;
	int key = 0;
	OBS_DATA *range = new OBS_DATA();
	double dt_e = 0;
	double temp_t = 0;
	bool first = true;
	int i, j;
	int val;
	i = 0;
	val = 0;
	while (1)
	{
		/*文件预处理*/
		for (; i < d - 2; i++) // 同步
		{
			if (Buff[i] == OEM4SYNC1 && Buff[i + 1] == OEM4SYNC2 && Buff[i + 2] == OEM4SYNC3)
			{
				break;
			}
		}
		key++;
		if (i + OEM4HLEN >= d) // 消息不完整，跳出
			break;

		for (j = 0; j < OEM4HLEN; j++)
			TempBuff[j] = Buff[i + j]; // 拷贝消息头到待解码缓存中

		len = U2(TempBuff + 8) + OEM4HLEN; // 消息头+消息体长度
		// if (OEM4HLEN != U1(TempBuff + 3))
		//{
		//	i += len + 4;
		//	continue;
		// }

		if ((len + 4 + i) > d || len > MAXRAWLEN) // 消息不完整，跳出
			break;

		for (j = OEM4HLEN; j < len + 4; j++) // 拷贝消息体到待解码缓存中
			TempBuff[j] = Buff[i + j];

		msgID = U2(TempBuff + 4);

		if (CRC32(TempBuff, len) != UCRC32(TempBuff + len, 4)) // 检验CRC32
		{
			i += len + 4;
			continue;
		}

		/*录入文件头信息*/
		msgTYPE = (U1(TempBuff + 6) >> 4) & 0X3;
		gpstime = new GPSTIME(U2(TempBuff + 14), (double)(U4(TempBuff + 16) * 1e-3));
		if (msgTYPE != 0)
			continue;
		int prn = 0;
		/*处理不同数据信息*/
		switch (msgID)
		{
		case ID_RANGE:
			range->OBS_TIME->Week = gpstime->Week;
			range->OBS_TIME->SecOfWeek = gpstime->SecOfWeek;
			range->Sate_Num = U4(TempBuff + OEM4HLEN);
			if (first)
				dt_e = 1;
			else
				dt_e = gpstime->SecOfWeek - temp_t;
			if (dt_e = 0)
				break;
			temp_t = gpstime->SecOfWeek;
			val = decode_RANGE(TempBuff + OEM4HLEN + 4, range->Sate_Num, range);
			if (Cal_SPP(result, range, GPS_eph, BDS_eph, dt_e, first))
				first = false;
			for (vector<Satellate *>::iterator it = range->GPS_SATE.begin(); it != range->GPS_SATE.end(); it++)
			{
				if (NULL != *it)
				{
					delete *it;
					*it = NULL;
				}
			}
			range->GPS_SATE.clear();
			for (vector<Satellate *>::iterator it = range->BDS_SATE.begin(); it != range->BDS_SATE.end(); it++)
			{
				if (NULL != *it)
				{
					delete *it;
					*it = NULL;
				}
			}
			range->BDS_SATE.clear();
			break;
		case ID_GPSEPHEMERIS:
			prn = U4(TempBuff + OEM4HLEN);
			decode_GPSEPH_STAT(TempBuff + OEM4HLEN, GPS_eph[prn - 1]);
			break;
		case ID_BDSEPHEMERIS:
			prn = U4(TempBuff + OEM4HLEN);
			decode_BDSEPH_STAT(TempBuff + OEM4HLEN, BDS_eph[prn - 1]);
			break;
		default:
			break;
		}
		delete gpstime;
		i += len + 4;
	}
	//---------------解码后，缓存的处理-------------------//
	for (j = 0; j < d - i; j++)
		Buff[j] = Buff[i + j];

	d = j; // 解码后，缓存中剩余的尚未解码的字节数
	//---------------解码后，缓存的处理-------------------//
	return val;
}