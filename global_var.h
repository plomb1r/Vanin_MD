#pragma once
#pragma region Кооординаты (массивы от числа частиц, где i - номер частицы)
double* coordx; 
double* coordy;
double* coordz;
#pragma endregion

#pragma region Скорости (массивы от числа частиц, где i - номер частицы)
double* vx;	 
double* vy;
double* vz;
#pragma endregion

#pragma region Потенциал Леннарда-Джонса
double U;
#pragma endregion

#pragma region Силы взаимодействия (массивы от числа частиц и вспомогательный скаляр)
double* Fx; 
double* Fy;
double* Fz;
 
#pragma endregion
double F;
 ;
