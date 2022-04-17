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

void fileWrite(int);

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
	Ekin = (double*)malloc(NUMBERPARTICLES * sizeof(double));
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
#pragma region Декларирование переменных, участвующих в расчетах потцениала и сил взаимодействия и макропараметров
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
double vm_x, vm_y, vm_z;
#pragma endregion
double recursive(double x, double y, double z) {
	return x * x + y * y + z * z;
}
double Pressure()
// computes pressure of the system
{
	// computes v_mean
	double v_mean[] = { 0, 0, 0 };
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		v_mean[0] += vx[i];
		v_mean[1] += vy[i];
		v_mean[2] += vz[i];
	}

	v_mean[0] /= NUMBERPARTICLES;
	v_mean[1] /= NUMBERPARTICLES;
	v_mean[2] /= NUMBERPARTICLES;

	// computes P
	double P = 0;
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		P += MASS * recursive(vx[i] - v_mean[0], vy[i] - v_mean[1], vz[i] - v_mean[2]);
		for (int j = i + 1; j < NUMBERPARTICLES; ++j)
		{
			double r_ij[] = { coordx[i] - coordx[j], coordy[i] - coordy[j], coordz[i] - coordz[j] };
			double r_ij_abs = sqrt(recursive(r_ij[0], r_ij[1], r_ij[2]));
			double F_ij = F12;
			for (int alpha = 0; alpha < 3; ++alpha)
			{
				double F_ij_alpha = F_ij * r_ij[alpha] / r_ij_abs;
				P += r_ij[alpha] * F_ij_alpha;
			}
		}
	}
	P = P / (3 * VOLUME);
	return P;
}
double E_thermal()
// computes thermal energy of the system
{
	double v_mean[] = { 0, 0, 0 };
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		v_mean[0] += vx[i];
		v_mean[1] += vy[i];
		v_mean[2] += vz[i];
	}

	v_mean[0] /= NUMBERPARTICLES;
	v_mean[1] /= NUMBERPARTICLES;
	v_mean[2] /= NUMBERPARTICLES;

	double E_term = 0;
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		E_term += recursive(vx[i] - v_mean[0], vy[i] - v_mean[1], vz[i] - v_mean[2]);
	}

	E_term = E_term * MASS / 2;
	return E_term;
}
void CalculationVectors() { //вычисление векторов и модулей rx12, ry12, rz12, rx21, ry21, rz21, r12_abs, r21_abs, потценциала и сил взаимодействия	
#pragma region Создание массивов под предыдущие значения для вычисления разностной схемы Верле
	double Fx_pre[NUMBERPARTICLES];
	double Fy_pre[NUMBERPARTICLES];
	double Fz_pre[NUMBERPARTICLES];
#pragma endregion
	double Pxx, Pyy, Pzz, mv2_x, mv2_y, mv2_z, rF_x, rF_y, rF_z;
	double v_T;
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
	//на 0 шаге фиксируем значения
	m = (SIGMA / r12_abs);
	m2 = m * m;
	m4 = m2 * m2;
	m6 = m4 * m2;
	m12 = m6 * m6;
	n = (SIGMA / r21_abs);
	n2 = n * n;
	n4 = n2 * n2;
	n6 = n4 * n2;
	n12 = n6 * n6;
	U12 = 4 * EPS * (m12 - m6);
	U21 = 4 * EPS * (n12 - n6); //потенциал Леннарда-Джонса
	F12 = (24 * EPS / r12_abs) * (2 * m12 - m6); //сила взаимодействия 2 частицы на 1	
	F21 = (24 * EPS / r21_abs) * (2 * n12 - n6);
	Fx[0] = (F12 * rx12) / r12_abs;
	Fx[1] = (F21 * rx21) / r21_abs;
	Fy[0] = (F12 * ry12) / r12_abs;
	Fy[1] = (F21 * ry21) / r21_abs;
	Fz[0] = (F12 * rz12) / r12_abs;
	Fz[1] = (F21 * rz21) / r21_abs;

	Ekin_sys = 0;
	for (int i = 0; i < NUMBERPARTICLES; i++) {
		Ekin[i] = MASS * (pow(vx[i], 2) + pow(vy[i], 2) + pow(vz[i], 2)) / 2;
		Ekin_sys += Ekin[i];
	}
	double v_mean[] = { 0, 0, 0 };
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		v_mean[0] += vx[i];
		v_mean[1] += vy[i];
		v_mean[2] += vz[i];
	}

	v_mean[0] /= NUMBERPARTICLES;
	v_mean[1] /= NUMBERPARTICLES;
	v_mean[2] /= NUMBERPARTICLES;

	double Eterm = 0.0;
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		Eterm += recursive(vx[i] - v_mean[0], vy[i] - v_mean[1], vz[i] - v_mean[2]);
	}

	Eterm = Eterm * MASS / 2;
	Epot = (U12 + U21) / 2;
	Eint = Epot + Eterm;
	E = Ekin_sys + Epot;
	for (int i = 0; i < NUMBERPARTICLES; i++) {
		vm_x = 0;
		vm_y = 0;
		vm_z = 0;
	}
	for (int i = 0; i < NUMBERPARTICLES; i++) {
		vm_x += vx[i];
		vm_y += vy[i];
		vm_z += vz[i];
	}
	vm_x /= NUMBERPARTICLES;
	vm_y /= NUMBERPARTICLES;
	vm_z /= NUMBERPARTICLES;

	v_T = 0;
	for (int i = 0; i < NUMBERPARTICLES; i++) {
		v_T += pow(vx[i] - vm_x, 2) + pow(vy[i] - vm_y, 2) + pow(vz[i] - vm_z, 2);
		T = MASS * v_T / (3 * NUMBERPARTICLES * K_B);
	}

	/*mv2_x = 0;
	mv2_y = 0;
	mv2_z = 0;
	for (int i = 0; i < NUMBERPARTICLES; i++) {	
		mv2_x += MASS * pow(vx[i] - vm_x, 2);
		mv2_y += MASS * pow(vy[i] - vm_y, 2);
		mv2_z += MASS * pow(vz[i] - vm_z, 2);
	}
	rF_x = rx12 * Fx[0];	
	rF_y = ry12 * Fy[0];
	rF_z = rz12 * Fz[0];

	Pxx = (1 / VOLUME) * (mv2_x + rF_x);
	Pyy = (1 / VOLUME) * (mv2_y + rF_y);
	Pzz = (1 / VOLUME) * (mv2_z + rF_z);

	P = (Pxx + Pyy + Pzz) / 3;*/
	fileWrite(0);

	//считаем схему Верле
	for (int i = 1; i < NSTEPS; i++)
	{
		Fx_pre[0] = Fx[0];
		Fx_pre[1] = Fx[1];
		Fy_pre[0] = Fy[0];
		Fy_pre[1] = Fy[1];
		Fz_pre[0] = Fz[0];
		Fz_pre[1] = Fz[1];
		for (int i = 0; i < NUMBERPARTICLES; ++i) {
			coordx[i] += vx[i] * STEP + Fx[i] / (2 * MASS) * (STEP * STEP);
			coordy[i] += vy[i] * STEP + Fy[i] / (2 * MASS) * (STEP * STEP);
			coordz[i] += vz[i] * STEP + Fz[i] / (2 * MASS) * (STEP * STEP);
		}
		//учет ПГУ
		for (int i = 0; i < NUMBERPARTICLES; ++i) {
			if (coordx[i] >= LX) {
				coordx[i] -= LX;
			}
			if (coordy[i] >= LY) {
				coordy[i] -= LY;
			}
			if (coordz[i] >= LZ) {
				coordz[i] -= LZ;
			}

			if (coordx[i] < 0) {
				coordx[i] += LX;
			}
			if (coordy[i] < 0) {
				coordy[i] += LY;
			}

			if (coordz[i] < 0) {
				coordz[i] += LZ;
			}
		}
		rx12 = coordx[0] - coordx[1];
		ry12 = coordy[0] - coordy[1];
		rz12 = coordz[0] - coordz[1];
		rx21 = coordx[1] - coordx[0];
		ry21 = coordy[1] - coordy[0];
		rz21 = coordz[1] - coordz[0];
		r12_abs = sqrt(rx12 * rx12 + ry12 * ry12 + rz12 * rz12);
		r21_abs = sqrt(rx21 * rx21 + ry21 * ry21 + rz21 * rz21);
		m = (SIGMA / r12_abs);
		m2 = m * m;
		m4 = m2 * m2;
		m6 = m4 * m2;
		m12 = m6 * m6;
		n = (SIGMA / r21_abs);
		n2 = n * n;
		n4 = n2 * n2;
		n6 = n4 * n2;
		n12 = n6 * n6;
		U12 = 4 * EPS * (m12 - m6);
		U21 = 4 * EPS * (n12 - n6);
		F12 = ((24 * EPS) / r12_abs) * ((2 * m12) - (m6));
		F21 = ((24 * EPS) / r21_abs) * ((2 * n12) - (n6));
		Fx[0] = F12 * rx12 / r12_abs;
		Fy[0] = F12 * ry12 / r12_abs;
		Fz[0] = F12 * rz12 / r12_abs;
		Fx[1] = F21 * rx21 / r21_abs;
		Fy[1] = F21 * ry21 / r21_abs;
		Fz[1] = F21 * rz21 / r21_abs;
		vx[0] = vx[0] + (Fx[0] + Fx_pre[0]) / (2 * MASS) * STEP;
		vy[0] = vy[0] + (Fy[0] + Fy_pre[0]) / (2 * MASS) * STEP;
		vz[0] = vz[0] + (Fz[0] + Fz_pre[0]) / (2 * MASS) * STEP;
		vx[1] = vx[1] + (Fx[1] + Fx_pre[1]) / (2 * MASS) * STEP;
		vy[1] = vy[1] + (Fy[1] + Fy_pre[1]) / (2 * MASS) * STEP;
		vz[1] = vz[1] + (Fz[1] + Fz_pre[1]) / (2 * MASS) * STEP;
		Ekin_sys = 0;
		for (int i = 0; i < NUMBERPARTICLES; i++) {
			Ekin[i] = MASS * (pow(vx[i], 2) + pow(vy[i], 2) + pow(vz[i], 2)) / 2;
			Ekin_sys += Ekin[i];
		}
		Epot = (U12 + U21) / 2;
		E = Ekin_sys + Epot;
		vm_x = 0;
		vm_y = 0;
		vm_z = 0;
		for (int i = 0; i < NUMBERPARTICLES; i++) {
			vm_x += vx[i];
			vm_y += vy[i];
			vm_z += vz[i];
		}
		vm_x /= NUMBERPARTICLES;
		vm_y /= NUMBERPARTICLES;
		vm_z /= NUMBERPARTICLES;
		v_T = 0;
		for (int i = 0; i < NUMBERPARTICLES; i++) {
			v_T += pow(vx[i] - vm_x, 2) + pow(vy[i] - vm_y, 2) + pow(vz[i] - vm_z, 2);
			T = MASS * v_T / (3 * NUMBERPARTICLES * K_B);
		}	
		/*Fx[0] = F12 * rx12 / r12_abs;
		Fy[0] = F12 * ry12 / r12_abs;
		Fz[0] = F12 * rz12 / r12_abs;
		Fx[1] = F12 * rx21 / r21_abs;
		Fy[1] = F12 * ry21 / r21_abs;
		Fz[1] = F12 * rz21 / r21_abs;
		Pxx = (MASS / VOLUME) * ((vx[0] - ((vx[0] + vx[1]) / NUMBERPARTICLES)) * (vx[0] - ((vx[0] + vx[1]) / NUMBERPARTICLES)) + (vx[1] - ((vx[0] + vx[1]) / NUMBERPARTICLES)) * (vx[1] - ((vx[0] + vx[1]) / NUMBERPARTICLES))) + (1 / VOLUME) * (rx12 * (Fx[0] - Fx[1]) + rx21 * (Fx[1] - Fx[0]));
		Pyy = (MASS / VOLUME) * ((vy[0] - ((vy[0] + vy[1]) / NUMBERPARTICLES)) * (vy[0] - ((vy[0] + vy[1]) / NUMBERPARTICLES)) + (vy[1] - ((vy[0] + vy[1]) / NUMBERPARTICLES)) * (vy[1] - ((vy[0] + vy[1]) / NUMBERPARTICLES))) + (1 / VOLUME) * (ry12 * (Fy[0] - Fy[1]) + ry21 * (Fy[1] - Fy[0]));
		Pzz = (MASS / VOLUME) * ((vz[0] - ((vz[0] + vz[1]) / NUMBERPARTICLES)) * (vz[0] - ((vz[0] + vz[1]) / NUMBERPARTICLES)) + (vz[1] - ((vz[0] + vz[1]) / NUMBERPARTICLES)) * (vz[1] - ((vz[0] + vz[1]) / NUMBERPARTICLES))) + (1 / VOLUME) * (rz12 * (Fz[0] - Fz[1]) + rz21 * (Fz[1] - Fz[0]));
		P = (Pxx + Pyy + Pzz) / 3;*/
		fileWrite(i);
	}
#pragma endregion	
}

void fileWrite(int iter = 0) {
	FILE* filew;
	errno_t error;
	error = fopen_s(&filew, "G:\\Vanin_MD_11.txt", "a");
	if (error != 0)
	{
		std::cout << "Error:" + error << std::endl;
	}
	else
	{
		double Eterm = E_thermal();
		double P = Pressure();
		fprintf(filew, "Step = %d\n", iter);
		//fprintf(filew, "%1.8f %1.8f %1.8f \n", coordx[0] + 0.0, coordy[0] + 0.0, coordz[0] + 0.0);
		//fprintf(filew, "r2 = (rx2; ry2; rz2) = (%1.8f; %1.8f; %1.8f) \n", coordx[1] + 0.0, coordy[1] + 0.0, coordz[1] + 0.0);
		//fprintf(filew, "r12_abs = %1.8f \n", r12_abs);
		fprintf(filew, "Ekin = %1.8f \n", Ekin_sys);
		fprintf(filew, "Eterm = %1.8f \n", Eterm);
		fprintf(filew, "Epot = % 1.8f \n", Epot);
		fprintf(filew, "Eint = %1.8f \n", Eint = Epot + Eterm);
		fprintf(filew, "E = %1.8f \n", E);
		fprintf(filew, "T = %1.8f \n", T);
		fprintf(filew, "P = %1.8f \n", P);
		fprintf(filew, "\n");
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
	start_cond_two_particles();
	CalculationVectors();
	ClearMemory();
}


int main()
{
	MD();
	return 0;
}