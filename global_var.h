#pragma once
#pragma region Кооординаты (массивы от числа частиц, где i - номер частицы)
double* coordx = (double*)malloc(NUMBERPARTICLES * sizeof(double)); // координаты
double* coordy = (double*)malloc(NUMBERPARTICLES * sizeof(double));
double* coordz = (double*)malloc(NUMBERPARTICLES * sizeof(double));
#pragma endregion

#pragma region Скорости (массивы от числа частиц, где i - номер частицы)
double* vx = (double*)malloc(NUMBERPARTICLES * sizeof(double));	 // скорости 
double* vy = (double*)malloc(NUMBERPARTICLES * sizeof(double));
double* vz = (double*)malloc(NUMBERPARTICLES * sizeof(double));

#pragma endregion

// double U потенциальная энергия(скаляр)


#pragma region Силы взаимодействия (массивы от числа частиц и вспомогательный скаляр)
 /*
double* Fx = (double*)malloc(NUMBERPARTICLES * sizeof(double)); // силы
double* Fy = (double*)malloc(NUMBERPARTICLES * sizeof(double));
double* Fz = (double*)malloc(NUMBERPARTICLES * sizeof(double));
 */
 
#pragma endregion
 ;
