#pragma once
#pragma region Кооординаты (массивы от числа частиц, где i - номер частицы)
double coordx[NUMBERPARTICLES];
double coordy[NUMBERPARTICLES];
double coordz[NUMBERPARTICLES];
#pragma endregion

#pragma region Скорости (массивы от числа частиц, где i - номер частицы)

 double vx[NUMBERPARTICLES];
 double vy[NUMBERPARTICLES];
 double vz[NUMBERPARTICLES];

#pragma endregion

//const double U потенциальная энергия(скаляр)


#pragma region Силы взаимодействия (массивы от числа частиц и вспомогательный скаляр)

 double Fx[NUMBERPARTICLES];
 double Fy[NUMBERPARTICLES];
 double Fz[NUMBERPARTICLES];
 double F;

#pragma endregion
 ;
