#pragma once
//constants.h - ���������� ��������� � ������������ � �������� ��������� � ���������
const double MASS = 66.335; //����� ����� ������� (���� ������)
const double K_B = 1.380648528; //���������� ���������

#pragma region ��������� ���������� �.-��. ��� ������
const double EPS = 1.712; //��������� ������� (������� ������������� ���)
const double SIGMA = 0.3418; //�������� ����� ��������������
const double RCUT =  2.5 * SIGMA; //������ ��������� ����������
const double RCUT2 = RCUT * RCUT; 
const double ACRIST = 1.2; //����� ����� ������������ ������
#pragma endregion

#pragma region ������� ������� �� ���� ���������
const double LX = NUMCRIST_X * ACRIST;
const double LY = NUMCRIST_Y * ACRIST;
const double LZ = NUMCRIST_Z * ACRIST;
#pragma endregion

//LX = LY = LZ = ;
const double VOLUME = LX * LY * LZ; //����� �������
;