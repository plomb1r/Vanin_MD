﻿#include <stdio.h> 
#define _CRT_SECURE_NO_DEPRECATE 
#define PLOTTING_DATA
#include <math.h>  
#include <stdlib.h> 
#include "params.h"
#include "constants.h"
#include "global_var.h"
#include "start_cond.h"



void getMemoryForArrays() { 
	coordx = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	coordy = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	coordz = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	vx = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	vy = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	vz = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	Fx = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	Fy = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	Fz = (double*)malloc(NUMBERPARTICLES * sizeof(double));
	for (int i = 0; i < NUMBERPARTICLES; ++i) {	
		coordx[i] = 0;
		coordy[i] = 0;
		coordz[i] = 0;
		vx[i] = 0;
		vy[i] = 0;
		vz[i] = 0;
		Fx[i] = 0;
		Fy[i] = 0;
		Fz[i] = 0;
	}
}

void clearMemory() { //освобождаем память
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


double U(double r) { 
	return 4 * EPS * ((pow(SIGMA/r, 10) - pow(SIGMA/r, 6)));
}

double F(double r) {
	return -4 * EPS * (-10 * pow((SIGMA / r), 11) - 6 * pow((SIGMA / r), 7));
}

double recursive(double x, double y, double z) {
	return x * x + y * y + z * z;
}

double sqrtRecursive(double x, double y, double z) {
	return sqrt(recursive(x, y, z));
}


void LJpotentional() { 
	for (int i = 0; i < NUMBERPARTICLES; ++i) {
		Fx[i] = 0;
		Fy[i] = 0;
		Fz[i] = 0; 
#pragma region Виртуальные частицы Lx/2
		for (int j = 0; j < NUMBERPARTICLES; ++j) {
			if (i == j) continue;
			double r_ij[] = { coordx[i] - coordx[j], coordy[i] - coordy[j], coordz[i] - coordz[j] };
			if (r_ij[0] < -LX / 2) {
				r_ij[0] += LX;
			}
			if (r_ij[0] >= LX / 2) {
				r_ij[0] -= LX;
			}
			if (r_ij[1] < -LY / 2) {
				r_ij[1] += LY;
			}
			if (r_ij[1] >= LY / 2) {
				r_ij[1] -= LY;
			}
			if (r_ij[2] < -LZ / 2) {
				r_ij[2] += LZ;
			}
			if (r_ij[2] >= LZ / 2) {
				r_ij[2] -= LZ;
			}
			double r_ij_square = recursive(r_ij[0], r_ij[1], r_ij[2]);
			double F_ij = F(r_ij_square);
			Fx[i] += F_ij * r_ij[0] / r_ij_square;
			Fy[i] += F_ij * r_ij[1] / r_ij_square;
			Fz[i] += F_ij * r_ij[2] / r_ij_square;
		}
	}
#pragma endregion	
}

	void start_cond_two_particles();

	double Ek() { 
		double Ek = 0;
		for (int i = 0; i < NUMBERPARTICLES; ++i) {
			Ek += recursive(vx[i], vy[i], vz[i]);
		}
		Ek = Ek * MASS / 2;
		return Ek;
	}

	double Et() { 
		double vCentral[] = { 0, 0, 0 }; 
		for (int i = 0; i < NUMBERPARTICLES; ++i) {
			vCentral[0] += vx[i];
			vCentral[1] += vy[i];
			vCentral[2] += vz[i];
		}
		vCentral[0] /= NUMBERPARTICLES;
		vCentral[1] /= NUMBERPARTICLES;
		vCentral[2] /= NUMBERPARTICLES;
		double Et = 0;
		for (int i = 0; i < NUMBERPARTICLES; ++i) {
			Et += recursive(vx[i] - vCentral[0], vy[i] - vCentral[1], vz[i] - vCentral[2]);
		}
		Et = Et * MASS / 2; 
		return Et;
	}

	double Ep() { 
		double Ep = 0;
		for (int i = 0; i < NUMBERPARTICLES; ++i) {
			for (int j = 0; j < NUMBERPARTICLES; ++j) {
				if (i == j) continue;
				double r_ij[] = { coordx[i] - coordx[j], coordy[i] - coordy[j], coordz[i] - coordz[j] };
				double r_ij_abs = sqrtRecursive(r_ij[0], r_ij[1], r_ij[2]);
				Ep += U(r_ij_abs);
			}
		}
		Ep /= 2; 
		return Ep;
	}

	double temperatureOfSystem(double Et) { 
		return (2 * Et) / (3 * NUMBERPARTICLES * K_B);
	}

	double pressuareOfSystem() { 
		double v_cetral[] = { 0, 0, 0 }; 
		for (int i = 0; i < NUMBERPARTICLES; ++i) {
			v_cetral[0] += vx[i];
			v_cetral[1] += vy[i];
			v_cetral[2] += vz[i];
		}
		v_cetral[0] /= NUMBERPARTICLES;
		v_cetral[1] /= NUMBERPARTICLES;
		v_cetral[2] /= NUMBERPARTICLES;
		double pressuareOfSystem = 0;
		for (int i = 0; i < NUMBERPARTICLES; ++i) {
			pressuareOfSystem += MASS * recursive(vx[i] - v_cetral[0], vy[i] - v_cetral[1], vz[i] - v_cetral[2]);
			for (int j = i + 1; j < NUMBERPARTICLES; ++j) {
				double r_ij[] = { coordx[i] - coordx[j], coordy[i] - coordy[j], coordz[i] - coordz[j] };
				double r_ij_square = recursive(r_ij[0], r_ij[1], r_ij[2]);
				double F_ij = F(r_ij_square);
				for (int alpha = 0; alpha < 3; ++alpha)
				{
					double F_ij_alpha = F_ij * r_ij[alpha] / r_ij_square;
					pressuareOfSystem += r_ij[alpha] * F_ij_alpha;
				}
			}
		}
		pressuareOfSystem = pressuareOfSystem / (3 * VOLUME);
		return pressuareOfSystem;
	}


void writeToFile(FILE* outStream, int iter) {
#if defined(PLOTTING_DATA)
	fprintf(outStream, "%d, %.8f, %.8f, %.8f, %.8f \n", iter, coordx[0], coordy[0], coordx[1], coordy[1]);
#else
	fprintf(outStream, "Step = %d\n", iter);

	for (int i = 0; i < NUMBERPARTICLES; ++i) { 
		fprintf(outStream, "r%d = (rx%d; ry%d; rz%d) = (%.8f; %.8f; %.8f)\n", i + 1, i + 1, i + 1, i + 1, coordx[i], coordy[i], coordz[i]);
	}

	double r12[] = { coordx[0] - coordx[1], coordy[0] - coordy[1], coordz[0] - coordz[1] };

	double r12_abs = sqrtRecursive(r12[0], r12[1], r12[2]);
	double r12_square = recursive(r12[0], r12[1], r12[2]);
	fprintf(outStream, "r12_abs = %.8f\n", r12_abs);

	double U12 = U(r12_square);
	fprintf(outStream, "U12 = %.8f\n", U12);

	double F12 = F(r12_square) / r12_abs;
	fprintf(outStream, "F12 = %.8f\n", F12);

	for (int i = 0; i < 1; ++i) { 
		fprintf(outStream, "F%d = (Fx%d; Fy%d; Fz%d) = (%.8f; %.8f; %.8f)\n", i + 1, i + 1, i + 1, i + 1, Fx[i], Fy[i], Fz[i]);
	}

	for (int i = 0; i < NUMBERPARTICLES; ++i) { 
		fprintf(outStream, "v%d = (vx%d; vy%d; vz%d) = (%.8f; %.8f; %.8f)\n", i + 1, i + 1, i + 1, i + 1, vx[i], vy[i], vz[i]);
	}
	double Ekin = Ek();
	fprintf(outStream, "Ekin = %.8f\n", Ekin);

	double Eterm = Et();
	fprintf(outStream, "Eterm = %.8f\n", Eterm);

	double Epot = Ep();
	fprintf(outStream, "Epot = %.8f\n", Epot);

	double Ein = Eterm + Epot;
	fprintf(outStream, "Eint = %.8f\n", Ein);

	double E = Ekin + Epot;
	fprintf(outStream, "E = %.8f\n", E);

	double T = temperatureOfSystem(Eterm);
	fprintf(outStream, "T = %.8f\n", T);

	double P = pressuareOfSystem();
	fprintf(outStream, "P = %.8f\n", P);
	
	fprintf(outStream, "\n"); 
#endif
}

void PBC() {
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
}
void verletAlgorithm() {
	for (int i = 0; i < NUMBERPARTICLES; ++i) {
		coordx[i] += vx[i] * STEP + Fx[i] / (2 * MASS) * (STEP * STEP);
		coordy[i] += vy[i] * STEP + Fy[i] / (2 * MASS) * (STEP * STEP);
		coordz[i] += vz[i] * STEP + Fz[i] / (2 * MASS) * (STEP * STEP);
	}
	PBC();
	for (int i = 0; i < NUMBERPARTICLES; ++i) {
		vx[i] += Fx[i] / (2 * MASS) * STEP;
		vy[i] += Fy[i] / (2 * MASS) * STEP;
		vz[i] += Fz[i] / (2 * MASS) * STEP;
	}
	LJpotentional();
	for (int i = 0; i < NUMBERPARTICLES; ++i) {
		vx[i] += Fx[i] / (2 * MASS) * STEP;
		vy[i] += Fy[i] / (2 * MASS) * STEP;
		vz[i] += Fz[i] / (2 * MASS) * STEP;
	}
}



void MD() {
	getMemoryForArrays();
	FILE* outStream;
#if defined(PLOTTING_DATA)
	char output_file_name[] = "plot_data.csv";
#else
	char output_file_name[] = "Vanin_MD_14.txt";
#endif
	outStream = fopen(output_file_name, "w");
		start_cond_two_particles(); 
		LJpotentional(); 
		writeToFile(outStream, 0); 
		for (int step = 1; step <= LASTSTEP; ++step)
		{
			verletAlgorithm(); 
			writeToFile(outStream, step); 
		}
		fclose(outStream); 
		clearMemory(); 

}

int main() {
	MD();
	return 0;
}
