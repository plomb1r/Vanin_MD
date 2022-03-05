#include <iostream>
#include <string>
#include <fstream>
#include "params.h"
#include "constants.h"
#include "global_var.h"
#include "start_cond.h"
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <vector>
using std::vector;
using std::string;


void GetMemoryForArrays() { //выделяем память
	coordx = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	coordy = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	coordz = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	vx = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	vy = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	vz = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	Fx = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	Fy = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	Fz = (double*)malloc(NUMBERPARTICLES * sizeof(double));
}

void ClearMemory() { //освобождаем память
	free(coordx);
	free(coordy);
	free(coordz);
	free(vx);
	free(vy);
	free(vz);
	free(Fx);
	free(Fy);
	free(Fz);
}

vector<double> CalculationVectors() { //вычисление векторов и модулей rx12, ry12, rz12, rx21, ry21, rz21, r12_abs, r21_abs, потценциала и сил взаимодействия
	vector<double> resultVector;
#pragma region Декларирование переменных, участвующих в расчетах потцениала и сил взаимодействия
	double rx12 = 0.0;
	double ry12 = 0.0;
	double rz12 = 0.0;
	double rx21 = 0.0;
	double ry21 = 0.0;
	double rz21 = 0.0;
	double r12_abs = 0.0;
	double r21_abs = 0.0;
	double m, m2, m4, m6, m12 = 0.0;
	double n, n2, n4, n6, n12 = 0.0;
	double U12 = 0.0;
	double U21 = 0.0;
	double F12 = 0.0;	
	double F21 = 0.0;
#pragma endregion
	
#pragma region Нахождение координат векторов rx12,rx21
	rx12 = coordx[0] - coordx[1];
	ry12 = coordy[0] - coordy[1];
	rz12 = coordz[0] - coordz[1];
	rx21 = coordx[1] - coordx[0];
	ry21 = coordy[1] - coordy[0];
	rz21 = coordz[1] - coordz[0];
	r12_abs = sqrt(rx12 * rx12 + ry12 * ry12 + rz12 * rz12);  //модуль вектора r12
	r21_abs = sqrt(rx21 * rx21 + ry21 * ry21 + rz21 * rz21);  //модуль вектора r21
#pragma endregion

#pragma region Рассчет сил взаимодействия и потенциалов
	//(можно воспользоваться pow(), но было решено менять на * при кратности 2)
	m = (SIGMA / r12_abs);
	m2 = m * m;
	m4 = m2 * m2;
	m6 = m4 * m2;
	m12 = m6 * m6; //возведение (SIGMA/r12_abs) в 12 степень 
	n = (SIGMA / r21_abs);
	n2 = n * n;
	n4 = n2 * n2;
	n6 = n4 * n2;
	n12 = n6 * n6;
	U12 = 4 * EPS * (m12 - m6); //потенциал Леннарда-Джонса
	U21 = 4 * EPS * (n12 - n6); //потенциал Леннарда-Джонса
	F12 = (24 * EPS / r12_abs) * (2 * m12 - m6); //сила взаимодействия 2 частицы на 1	
	F21 = (24 * EPS / r21_abs) * (2 * n12 - n6);
	Fx[0] = (F12 * rx12) / r12_abs;
	Fx[1] = (F21 * rx21) / r21_abs;
	Fy[0] = (F12 * ry12) / r12_abs;
	Fy[1] = (F21 * ry21) / r21_abs;
	Fz[0] = (F12 * rz12) / r12_abs;
	Fz[1] = (F21 * rz21) / r21_abs;
	
	
#pragma endregion
	resultVector = { rx12, ry12, rz12, rx21, ry21, rz21, r12_abs, r21_abs, U12, U21, F12, Fx[0], Fy[0], Fz[0], Fx[1], Fy[1], Fz[1] }; //результирующие данные 
	return resultVector;
}

void PrintAndCalculateVectorsAndValuesData(vector<double> v) { 
	FILE* filew;
	errno_t error;
	error = fopen_s(&filew, "G:\\Vanin_MD_5.txt", "w");
	if (error != 0) 
	{
		std::cout << "Error:" + error << std::endl;
	}
	else 
	{
		fprintf(filew, "Step = 0 \n");
		fprintf(filew, "r1 = (rx1; ry1; rz1) = (%1.8f; %1.8f; %1.8f) \n", coordx[0], coordy[0], coordz[0]);
		fprintf(filew, "r2 = (rx2; ry2; rz2) = (%1.8f; %1.8f; %1.8f) \n", coordx[1], coordy[1], coordz[1]);
		fprintf(filew, "r12 = (rx12; ry12; rz12) = (%1.8f; %1.8f; %1.8f) \n", v[0], v[1], v[2]);
		fprintf(filew, "r21 = (rx21; ry21; rz21) = (%1.8f; %1.8f; %1.8f) \n", v[3], v[4], v[5]);
		fprintf(filew, "r12_abs = %1.8f \n", v[6]);
		fprintf(filew, "r21_abs = %1.8f \n", v[7]);
		fprintf(filew, "(rx12; ry12; rz12)/r12_abs = (%1.8f; %1.8f; %1.8f) \n", v[0] / v[6], v[1] / v[6], v[2] / v[6]);
		fprintf(filew, "(rx21; ry21; rz21)/r21_abs =  (%1.8f; %1.8f; %1.8f) \n", v[3] / v[6], v[4] / v[6], v[5] / v[6]);		
		fprintf(filew, "U12 = %1.8f \n", v[8]); //потенциал Леннарда-Джонса
		fprintf(filew, "U21 = %1.8f \n", v[9]); //потенциал Леннарда-Джонса
		fprintf(filew, "F12 = %1.8f \n", v[10]);//сила взаимодействия 
		fprintf(filew, "F1 = (F1x1; F1y1; F1z1) = (%1.8f; %1.8f; %1.8f) \n", v[11], v[12] + 0.0, v[13] + 0.0); 
		fprintf(filew, "F2 = (F2x2; F2y2; F2z2) = (%1.8f; %1.8f; %1.8f) \n", v[14], v[15] + 0.0, v[16] + 0.0); 
	}
	fclose(filew);
}

void freeMemory() { //освобождаем память
	free(coordx);
	free(coordy);
	free(coordz);
	free(vx);
	free(vy);
	free(vz);
	free(Fx);
	free(Fy);
	free(Fz);
}

void MD() { //создание функции Molecular Dynamics, в которой создается и декларируется экземпляр класса vector, вызывается функция н.у. для 2 частиц, а также проивзодятся манипуляции с памятью	
	GetMemoryForArrays();
	vector<double> v;
	start_cond_two_particles();
	v = CalculationVectors(); //вызов метода подсчета данных для вектора
	PrintAndCalculateVectorsAndValuesData(v);
	ClearMemory();
}


int main()
{
	MD();
	return 0;
}



