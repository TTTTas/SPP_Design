#pragma once
#pragma region header
#include<iostream>
#include<fstream>
#include<vector>
#include<bitset>
#include<string>
#include<iterator>
#include"transform.h"

using namespace std;
#define MAXRAWLEN 40960			//����ȡ���ݳ���
#define MAXNUM 8				//������
#define POLYCRC32 0xEDB88320	//CRC32У�������
#define OEM4SYNC1       0xAA    /* oem7/6/4 message start sync code 1 */
#define OEM4SYNC2       0x44    /* oem7/6/4 message start sync code 2 */
#define OEM4SYNC3       0x12    /* oem7/6/4 message start sync code 3 */
#define OEM4HLEN        28      /* oem7/6/4 message header length (bytes) */

/* message IDs */
#define ID_RANGE        43      /* oem7/6/4 range measurement */
#define ID_GPSEPHEMERIS 7       /* oem7/6 decoded gps ephemeris */
#define ID_BDSEPHEMERIS 1696    /* oem7/6 decoded bds ephemeris */

/*system IDs*/
#define SYS_GPS         0
#define SYS_BDS         4

#pragma endregion

struct Satellate
{
	unsigned short PRN = 0;//��ʶ��
	unsigned short SYS = -1;//ϵͳ
	int Phase_NUM = 0;//������

	/*����*/
	int SYG_TYPE[MAXNUM] = { -1,-1,-1,-1,-1,-1,-1,-1 };
	/*α��(m)*/
	double PSERA[MAXNUM] = { 0,0,0,0,0,0,0,0 };
	/*��λ(cycle)*/
	double PHASE[MAXNUM] = { 0,0,0,0,0,0,0,0 };
	/*������(HZ)*/
	double DOPPLER[MAXNUM] = { 0,0,0,0,0,0,0,0 };
	/*�����*/
	double SNR[MAXNUM] = { 0,0,0,0,0,0,0,0 };
	/*α�ྫ��*/
	double PSE_PREC[MAXNUM] = { 0,0,0,0,0,0,0,0 };
	/*�ز���λ����*/
	double PHA_PREC[MAXNUM] = { 0,0,0,0,0,0,0,0 };
	/*ͨ��״̬*/
	int LOCK_PSE[MAXNUM] = { 0,0,0,0,0,0,0,0 };
	int LOCK_PHA[MAXNUM] = { 0,0,0,0,0,0,0,0 };
	int PARITY[MAXNUM] = { 0,0,0,0,0,0,0,0 };

	/*�Ƿ��дֲ�*/
	bool Outlier = false;

};

struct EPHEMERIS
{
	/* data */
	//���Ǳ��
	int PRN = 0;
	//�����������
	double toe_tow;
	double toe_wn;
	double toc;
	double sqrt_A;
	double e;
	double M0;
	double omiga;
	double i0;
	double Omiga0;
	//�㶯������
	double delt_n;
	double dot_i;
	double dot_Omiga;
	double Cus;
	double Cuc;
	double Crs;
	double Crc;
	double Cis;
	double Cic;
	//ʱ��������
	/*��Ʈ*/
	double a_f0;
	double a_f1;
	double a_f2;
	/*Ⱥ�Ӳ�*/
	double T_GD1;
	double T_GD2;
	/*��������*/
	double AODC;
	int AODE;
	int IODC;
	//����ָ��
	unsigned int health;
	double URA;
	GPSTIME* Z_T;
	int IODE1;
	int IODE2;
};

struct EPOCH
{
	int PRN;
	int num = 0;
	vector<EPHEMERIS*> epoch;
};

struct OBS_DATA
{
	GPSTIME* OBS_TIME = new GPSTIME();
	int Sate_Num = 0;
	vector<Satellate*> GPS_SATE;
	vector<Satellate*> BDS_SATE;
};


#define U1(p) (*((uint8_t *)(p)))
#define I1(p) (*((int8_t  *)(p)))
static uint16_t U2(uint8_t* p) { uint16_t u; memcpy(&u, p, 2); return u; }//2�ֽ�����
static uint32_t U4(uint8_t* p) { uint32_t u; memcpy(&u, p, 4); return u; }//4�ֽ�����
static int32_t  I4(uint8_t* p) { int32_t  i; memcpy(&i, p, 4); return i; }//4�ֽ�����
static float    R4(uint8_t* p) { float    r; memcpy(&r, p, 4); return r; }//4�ֽ�float
static double   R8(uint8_t* p) { double   r; memcpy(&r, p, 8); return r; }//8�ֽ�double
/*���ұ�ʶ��AA4412*/
int check_syn(unsigned char* buff, unsigned char data);
/*����У����*/
unsigned int CRC32(const unsigned char* buff, int len);
/*����У����*/
unsigned int UCRC32(const unsigned char* buff, int len);
/*����ϵͳ*/
unsigned short decode_SYS(unsigned char* buff);

/*
* α��
* α�ྫ��
* �ز�
* �ز�����
* ������
* �����
* ����
* ����
* ͨ��״̬
*/
/*�������ǵ�Ƶ��۲�����*/
unsigned int decode_RANGE_STAT(unsigned char* buff, Satellate* sate);
/*����������������Ƶ��۲�����*/
unsigned int decode_RANGE(unsigned char* buff, int num, OBS_DATA* obs);
/*���뵥��GPS��������*/
unsigned int decode_GPSEPH_STAT(unsigned char* buff, EPHEMERIS* epoch);
/*��������GPS��������*/
unsigned int decode_GPSEPH(unsigned char* buff, EPOCH* gpse);
/*���뵥��BDS��������*/
unsigned int decode_BDSEPH_STAT(unsigned char* buff, EPHEMERIS* bdse);
/*��������BDS��������*/
unsigned int decode_BDSEPH(unsigned char* buff, EPOCH* bdse);
/*��ȡ�ļ�����*/
int readfile(const char* filepath, vector<OBS_DATA*>& range, EPOCH* gps, EPOCH* bds);

