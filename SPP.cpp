#include<iostream>
#include<fstream>
#include<vector>
#include<Eigen/Dense>
#include<iomanip>
#include<io.h>
#include<direct.h>

#include"transform.h"
#include"read.h"
#include"cal.h"
#include"sockets.h"

using namespace std;
using namespace Eigen;

#define _CRT_SECURE_NO_WARNINGS

int main()
{
    vector<OBS_DATA *> RANGE;
    EPOCH GPS_EPH[GPS_SAT_QUAN];
    EPOCH BDS_EPH[BDS_SAT_QUAN];
    Result_DATA *result = new Result_DATA();

    initial();

    double t = 1;
    bool first = true;
    FILE *DATA_Fobs;
    FILE *Pos_Fobs;
    int lenR, lenD;
    unsigned char curbuff[MAXRAWLEN];
    lenD = 0;
    unsigned char decBuff[2 * MAXRAWLEN];
    SOCKET NetGps;
    Configure CfgInfo;
    int ReadFlag;

    int choice = 0;
    printf("��ѡ�����뷽ʽ\n1. �ļ�\t2. ����\n");
    cin >> choice;

    time_t nowtime;
    time(&nowtime); // ��ȡ1970��1��1��0��0��0�뵽���ھ���������
    tm p;
    localtime_s(&p, &nowtime); // ������ת��Ϊ����ʱ��,���1900����,��Ҫ+1900,��Ϊ0-11,����Ҫ+1
    string filetime = to_string(p.tm_year + 1900) + "_" + to_string(p.tm_mon + 1) + "_" + to_string(p.tm_mday) + "_" + to_string(p.tm_hour) + "_" + to_string(p.tm_min) + "_" + to_string(p.tm_sec);
    string s1 = "C:\\Users\\Surface\\Desktop\\data\\logs\\";
    createDirectory(s1);
    s1 += filetime + string(".log");
    string s2 = "C:\\Users\\Surface\\Desktop\\data\\Pos\\";
    createDirectory(s2);
    s2 += filetime + string(".pos");
    CfgInfo.ObsDatFile = s1.c_str();
    CfgInfo.ResDatFile = s2.c_str();

    switch (choice)
    {
    case 1:
        readfile("NovatelOEM20211114-01.log", RANGE, GPS_EPH, BDS_EPH);
        for (int i = 0; i < RANGE.size(); i++)
        {
            if (i > 0)
                t = RANGE[i]->OBS_TIME->SecOfWeek - RANGE[i - 1]->OBS_TIME->SecOfWeek;
            if (Cal_SPP(result, RANGE[i], GPS_EPH, BDS_EPH, t, first))
                first = false;
            result->OUTPUT();
        }
        break;
    case 2:
        if (OpenSocket(NetGps, CfgInfo.NetIP, CfgInfo.NetPort) == false)
        {
            printf("The ip %s was not opened\n", CfgInfo.NetIP);
            return 0;
        }
        if ((DATA_Fobs = fopen(CfgInfo.ObsDatFile, "wb")) == NULL)
        {
            printf("The obs file %s was not opened\n", CfgInfo.ObsDatFile);
            exit(0);
        }
        if ((Pos_Fobs = fopen(CfgInfo.ResDatFile, "w")) == NULL)
        {
            printf("The obs file %s was not opened\n", CfgInfo.ResDatFile);
            exit(0);
        }
        while (1)
        {
            Sleep(980);
            if ((lenR = recv(NetGps, (char *)curbuff, MAXRAWLEN, 0)) > 0) // ��ȡ����
            {
                printf("%5d ", lenR);
                fwrite(curbuff, sizeof(unsigned char), lenR, DATA_Fobs); // ��¼���������������ļ���

                if ((lenD + lenR) > 2 * MAXRAWLEN)
                    lenD = 0;
                memcpy(decBuff + lenD, curbuff, lenR); // ����ƴ��
                lenD += lenR;

                ReadFlag = decodestream(result, decBuff, lenD); // ����
                if (ReadFlag != 1)
                {
                    printf("Data acquisition and decode failed \n");
                    continue;
                }
                else
                {
                    result->OUTPUT();
                    result->WRITEOUTPUT(Pos_Fobs);
                }
            }
            else
            {
                printf("NO MESSAGES IN!\n");
            }
        }

        break;
    default:
        break;
    }

    std::system("pause");
    return 0;
}