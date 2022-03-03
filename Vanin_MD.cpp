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
using std::setprecision;

void getStrings() {
	string empty;
	std::ofstream StreamOut;
	StreamOut.open("G:\\Vanin_MD_2.txt");
	if (StreamOut.is_open())
	{
		StreamOut << empty << "\n"
			<< "Number of Particles = " << NUMBERPARTICLES << "\n"
			<< "Number of Steps =  " << NSTEPS << "\n"
			<< "Size of System on X =  " << LX << "\n"
			<< "Boltzmann Constant  =  " << setprecision(10) << K_B << "\n"
			<< "Volume  =  " << VOLUME;
	}
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
	double m = 0.0;
	double n = 0.0;
	double m12 = 0.0;
	double m6 = 0.0;
	double n12 = 0.0;
	double n6 = 0.0;
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
	r12_abs = sqrt(rx12 * rx12 + ry12 * ry12 + rz12 * rz12);  //длина вектора r12
	r21_abs = sqrt(rx21 * rx21 + ry21 * ry21 + rz21 * rz21);  //длина вектора r21
#pragma endregion

#pragma region Рассчет сил взаимодействия и потенциалов
	//(можно воспользоваться pow(), но было решено менять на *)
	m = (SIGMA / r12_abs);
	m12 = m * m * m * m * m * m * m * m * m * m * m * m; //возведение (SIGMA/r12_abs) в 12 степень 
	m6 = m * m * m * m * m * m; //возведение (SIGMA/r12_abs) в 6 степень
	n = (SIGMA / r21_abs);
	n12 = n * n * n * n * n * n * n * n * n * n * n *n; //возведение (SIGMA/r21_abs) в 12 степень
	n6 = n * n * n * n * n * n; //возведение (SIGMA/r21_abs) в 6 степень
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
	std::ofstream StreamOut;
	StreamOut.open("G:\\Vanin_MD_5.txt");
	if (StreamOut.is_open())
	{
		StreamOut << "Step = 0" << "\n"
			<< "r1 = (rx1; ry1; rz1) = " << "(" << coordx[0] << "; " << coordy[0] << "; " << coordz[0] << ")" << "\n"
			<< "r2 = (rx2; ry2; rz2) = " << "(" << coordx[1] << "; " << coordy[1] << "; " << coordz[1] << ")" << "\n"
			<< "r12 = (rx12; ry12; rz12) = " << "(" << v[0] << "; " << v[1] << "; " << v[2] << ")" << "\n"
			<< "r21 = (rx21; ry21; rz21) = " << "(" << v[3] << "; " << v[4] << "; " << v[5] << ")" << "\n"
			<< "r12_abs = " << v[6] << "\n"
			<< "r21_abs = " << v[7] << "\n"
			<< "(rx21; ry21; rz21) / r12_abs = " << "(" << v[3] / v[6] << "; " << v[4] / v[6] << "; " << v[4] / v[6] << ")" << "\n"
			<< "(rx12; ry12; rz12) / r12_abs = " << "(" << v[0] / v[6] << "; " << v[1] / v[6] << "; " << v[2] / v[6] << ")" << "\n"
			<< "U12 = " << v[8] << "\n" //потенциал Леннарда-Джонса
			<< "U21 = " << v[9] << "\n" //потенциал Леннарда-Джонса
			<< "F12 = " << v[10] << "\n" //сила взаимодействия
			<< "F1 = (Fx1; Fy1; Fz1) = " << "(" << setprecision(6) << v[11] << "; " << v[12]+0.0 << "; " << v[13]+0.0 << ")" << "\n"
			<< "F2 = (Fx2; Fy2; Fz2) = " << "(" << setprecision(6) << v[14] << "; " << v[15]+0.0 << "; " << v[16] + 0.0 << ")";
	}
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

void MD() { //создание функции Molecular Dynamics, в которой создается и декларируется экземпляр класса vector, вызывается функция н.у. для 2 частиц
	vector<double> v;
	start_cond_two_particles();
	v = CalculationVectors(); //вызов метода подсчета данных для вектора
	PrintAndCalculateVectorsAndValuesData(v);
	//getStrings();
	freeMemory();
}



int main()
{
	setlocale(LC_ALL, "Russian");
	MD();
	return 0;
}



