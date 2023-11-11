#include"read.h"
#include"cal.h"
#include"data.h"

int check_syn(unsigned char* buff, unsigned char data)
{
	buff[0] = buff[1]; buff[1] = buff[2]; buff[2] = data;
	return (buff[0] == OEM4SYNC1 && buff[1] == OEM4SYNC2 && buff[2] == OEM4SYNC3);
}

unsigned int CRC32(const unsigned char* buff, int len)
{
	int i, j;
	unsigned int CRC = 0;
	for (i = 0; i < len; i++)
	{
		CRC ^= buff[i];
		for (j = 0; j < 8; j++)
		{
			if (CRC & 1)
				CRC = (CRC >> 1) ^ POLYCRC32;
			else
				CRC >>= 1;
		}
	}
	return CRC;
}

unsigned int UCRC32(const unsigned char* buff, int len)
{
	unsigned int UCRC32 = 0;
	for (int i = len; i > 0; i--)
	{
		UCRC32 += buff[i - 1] << 8 * (i - 1);
	}
	return UCRC32;
}

unsigned short decode_SYS(unsigned char* buff)
{
	unsigned short SYS = (U4(buff) >> 16) & 7;
	return SYS;
}

/*
* 伪距
* 伪距精度
* 载波
* 载波精度
* 多普勒
* 载噪比
* 波段
* 波长
* 通道状态
*/
unsigned int decode_RANGE_STAT(unsigned char* buff, Satellate* sate)
{
	int index = sate->Phase_NUM - 1;
	sate->PSERA[index] = R8(buff + 4);
	sate->PSE_PREC[index] = R4(buff + 12);
	sate->PHASE[index] = abs(R8(buff + 16));
	sate->PHA_PREC[index] = R4(buff + 24);
	sate->DOPPLER[index] = R4(buff + 28);
	sate->SNR[index] = R4(buff + 32);

	uint32_t stat = U4(buff + 40);
	sate->LOCK_PHA[index] = (stat >> 10) & 0X1;
	sate->PARITY[index] = (stat >> 11) & 0X1;
	sate->LOCK_PSE[index] = (stat >> 12) & 0X1;
	sate->SYG_TYPE[index] = (stat >> 21) & 0x1F;
	return 1;
}

unsigned int decode_RANGE(unsigned char* buff, int num, OBS_DATA* obs)
{
	vector<unsigned short>gps_prn;
	vector<unsigned short>bds_prn;
	vector<unsigned short>::iterator it;
	for (int i = 0; i < num; i++)
	{
		unsigned short sys = decode_SYS(buff + 44 * i + 40);
		unsigned short prn = U2(buff + 44 * i);
		switch (sys)
		{
		case SYS_GPS:
			it = find(gps_prn.begin(), gps_prn.end(), prn);
			if (it == gps_prn.end())
			{
				obs->GPS_SATE.push_back(new Satellate());
				obs->GPS_SATE.back()->PRN = prn;
				obs->GPS_SATE.back()->SYS = SYS_GPS;
				obs->GPS_SATE.back()->Phase_NUM = 1;
				decode_RANGE_STAT(buff + 44 * i, obs->GPS_SATE.back());
				gps_prn.push_back(prn);
			}
			else
			{
				int pos = it - gps_prn.begin();
				obs->GPS_SATE[pos]->Phase_NUM++;
				decode_RANGE_STAT(buff + 44 * i, obs->GPS_SATE[pos]);
			}
			break;
		case SYS_BDS:
			it = find(bds_prn.begin(), bds_prn.end(), prn);
			if (it == bds_prn.end())
			{
				obs->BDS_SATE.push_back(new Satellate());
				obs->BDS_SATE.back()->PRN = prn;
				obs->BDS_SATE.back()->SYS = SYS_BDS;
				obs->BDS_SATE.back()->Phase_NUM = 1;
				decode_RANGE_STAT(buff + 44 * i, obs->BDS_SATE.back());
				bds_prn.push_back(prn);
			}
			else
			{
				int pos = it - bds_prn.begin();
				obs->BDS_SATE[pos]->Phase_NUM++;
				decode_RANGE_STAT(buff + 44 * i, obs->BDS_SATE[pos]);
			}
			break;
		default:
			break;
		}
	}
	return 1;
}

unsigned int decode_GPSEPH_STAT(unsigned char* buff, EPHEMERIS* epoch)
{
	epoch->PRN = U4(buff);
	epoch->health = U4(buff + 12);
	epoch->IODE1 = U4(buff + 16);
	epoch->IODE2 = U4(buff + 20);
	epoch->toe_wn = U4(buff + 24);
	GPSTIME* temp_Z = new GPSTIME(U4(buff + 28), R8(buff + 4));
	epoch->Z_T = temp_Z;
	delete temp_Z;
	epoch->toe_tow = R8(buff + 32);
	epoch->sqrt_A = sqrt(R8(buff + 40));
	epoch->delt_n = R8(buff + 48); //rad/s
	epoch->M0 = R8(buff + 56);  //rad
	epoch->e = R8(buff + 64);
	epoch->omiga = R8(buff + 72); //rad
	epoch->Cuc = R8(buff + 80);
	epoch->Cus = R8(buff + 88);
	epoch->Crc = R8(buff + 96);
	epoch->Crs = R8(buff + 104);
	epoch->Cic = R8(buff + 112);
	epoch->Cis = R8(buff + 120);
	epoch->i0 = R8(buff + 128);
	epoch->dot_i = R8(buff + 136);
	epoch->Omiga0 = R8(buff + 144);
	epoch->dot_Omiga = R8(buff + 152);
	epoch->IODC = U4(buff + 160);
	epoch->toc = R8(buff + 164);
	epoch->T_GD1 = R8(buff + 172);
	epoch->a_f0 = R8(buff + 180);
	epoch->a_f1 = R8(buff + 188);
	epoch->a_f2 = R8(buff + 196);
	epoch->URA = R8(buff + 216);

	return 1;
}

unsigned int decode_GPSEPH(unsigned char* buff, EPOCH* gpse)
{
	int prn = U4(buff);
	gpse[prn - 1].PRN = prn;
	int wn = U4(buff + 24);
	double toe = R8(buff + 32);
	int check = 1;//默认不存在重复数据
	for (int i = 0; i < gpse[prn - 1].num; i++)
	{
		if (wn == gpse[prn - 1].epoch[i]->toe_wn && toe == gpse[prn - 1].epoch[i]->toe_tow)
			check = 0;
	}
	if (check)
	{
		EPHEMERIS* e = new EPHEMERIS();
		decode_GPSEPH_STAT(buff, e);
		gpse[prn - 1].epoch.push_back(e);
		gpse[prn - 1].num++;
		return 1;
	}
	else
	{
		return 0;
	}

}

unsigned int decode_BDSEPH_STAT(unsigned char* buff, EPHEMERIS* bdse)
{
	bdse->PRN = U4(buff);
	bdse->toe_wn = U4(buff + 4);
	bdse->URA = R8(buff + 8);
	bdse->health = U4(buff + 16);
	bdse->T_GD1 = R8(buff + 20);
	bdse->T_GD2 = R8(buff + 28);
	bdse->AODC = U4(buff + 36);
	bdse->toc = U4(buff + 40);
	bdse->a_f0 = R8(buff + 44);
	bdse->a_f1 = R8(buff + 52);
	bdse->a_f2 = R8(buff + 60);
	bdse->AODE = U4(buff + 68);
	bdse->toe_tow = U4(buff + 72);
	bdse->sqrt_A = R8(buff + 76);
	bdse->e = R8(buff + 84);
	bdse->omiga = R8(buff + 92);
	bdse->delt_n = R8(buff + 100);
	bdse->M0 = R8(buff + 108);
	bdse->Omiga0 = R8(buff + 116);
	bdse->dot_Omiga = R8(buff + 124);
	bdse->i0 = R8(buff + 132);
	bdse->dot_i = R8(buff + 140);
	bdse->Cuc = R8(buff + 148);
	bdse->Cus = R8(buff + 156);
	bdse->Crc = R8(buff + 164);
	bdse->Crs = R8(buff + 172);
	bdse->Cic = R8(buff + 180);
	bdse->Cis = R8(buff + 188);

	return 1;
}

unsigned int decode_BDSEPH(unsigned char* buff, EPOCH* bdse)
{
	int prn = U4(buff);
	bdse[prn - 1].PRN = prn;
	int wn = U4(buff + 4);
	int toe = U4(buff + 72);
	int check = 1;//默认不存在重复数据
	for (int i = 0; i < bdse[prn - 1].num; i++)
	{
		if (wn == bdse[prn - 1].epoch[i]->toe_wn && toe == bdse[prn - 1].epoch[i]->toe_tow)
			check = 0;
	}
	if (check)
	{
		EPHEMERIS* e = new EPHEMERIS();
		decode_BDSEPH_STAT(buff, e);
		bdse[prn - 1].epoch.push_back(e);
		bdse[prn - 1].num++;
		return 1;
	}
	else
	{
		return 0;
	}
}

int readfile(const char* filepath, vector<OBS_DATA*>& range, EPOCH* gps, EPOCH* bds)
{
	unsigned char data, buff[MAXRAWLEN];
	FILE* FObs;
	int len, headlen;
	int msgID, msgTYPE;
	GPSTIME* gpstime;
	int key = 0;

	if ((FObs = fopen(filepath, "rb")) == NULL)
	{
		printf("未能打开文件\n");
		return -1;
	}
	int i;
	for (i = 0;; i++)
	{
		/*文件预处理*/
		if (fread(&data, 1, 1, FObs) < 1)
		{
			cout << "END OF FILE!" << endl;
			fclose(FObs);
			return -2;
		}
		if (check_syn(buff, data))
		{
			key++;
			if (fread(buff + 3, 1, 7, FObs) < 7)
			{
				cout << "END OF FILE!" << endl;
				fclose(FObs);
				return -2;
			}
			if (headlen = U1(buff + 3) != OEM4HLEN)
			{
				cout << "HEADLENGTH ERROR!" << endl;
				fclose(FObs);
				return -1;
			}
			if ((len = U2(buff + 8) + OEM4HLEN) > MAXRAWLEN - 4)//检查数据长度，为数据头+数据，不包括CRC校验码
			{
				cout << "END OF FILE!" << endl;
				fclose(FObs);
				return -2;
			}

			if (fread(buff + 10, len - 6, 1, FObs) < 1)//读取包括CRC在内的所有内容
			{
				cout << "END OF FILE!" << endl;
				fclose(FObs);
				return -2;
			}



			if (CRC32(buff, len) != UCRC32(buff + len, 4))//检验CRC32
			{
				cout << "DATA ERROR! " << i << endl;
				continue;
			}

			/*录入文件头信息*/
			msgID = U2(buff + 4);
			msgTYPE = (U1(buff + 6) >> 4) & 0X3;
			gpstime = new GPSTIME(U2(buff + 14), (double)(U4(buff + 16) * 1e-3));

			OBS_DATA* obs = new OBS_DATA();
			/*处理不同数据信息*/
			switch (msgID)
			{
			case ID_RANGE:
				obs->OBS_TIME = gpstime;
				obs->Sate_Num = U4(buff + OEM4HLEN);
				decode_RANGE(buff + OEM4HLEN + 4, obs->Sate_Num, obs);
				range.push_back(obs);
				break;
			case ID_GPSEPHEMERIS:
				decode_GPSEPH(buff + OEM4HLEN, gps);
				break;
			case ID_BDSEPHEMERIS:
				decode_BDSEPH(buff + OEM4HLEN, bds);
				break;
			default:
				break;
			}
		}
	}

	return 1;
}

