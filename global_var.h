#pragma once
#pragma region Кооординаты (массивы от числа частиц, где i - номер частицы)
double* coordx = (double*)malloc(NUMBERPARTICLES * sizeof(double)); 
double* coordy = (double*)malloc(NUMBERPARTICLES * sizeof(double));
double* coordz = (double*)malloc(NUMBERPARTICLES * sizeof(double));
#pragma endregion

#pragma region Скорости (массивы от числа частиц, где i - номер частицы)
double* vx = (double*)malloc(NUMBERPARTICLES * sizeof(double));	 
double* vy = (double*)malloc(NUMBERPARTICLES * sizeof(double));
double* vz = (double*)malloc(NUMBERPARTICLES * sizeof(double));
#pragma endregion

#pragma region Потенциал Леннарда-Джонса
double U;
#pragma endregion

#pragma region Силы взаимодействия (массивы от числа частиц и вспомогательный скаляр)
double* Fx = (double*)malloc(NUMBERPARTICLES * sizeof(double)); 
double* Fy = (double*)malloc(NUMBERPARTICLES * sizeof(double));
double* Fz = (double*)malloc(NUMBERPARTICLES * sizeof(double));
 
#pragma endregion
double F;
 ;
