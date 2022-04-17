#pragma once
//constants.h - физические константы и используемые в расчетах константы и параметры
const double MASS = 66.335; //масса одной частицы (атом аргона)
const double K_B = 1.380648528; //постоянная Больцмана

#pragma region Параметры потенциала Л.-Дж. для аргона
const double EPS = 1.712; //параметры энергии (глубина потенциальной ямы)
const double SIGMA = 0.3418; //параметр длины взаимодействия
#pragma endregion

const double RCUT = 2.5 * SIGMA; //радиус обрезания потенциала
const double RCUT2 = RCUT * RCUT;
const double ACRIST = 2.0; //длина ребра элементарной ячейки

#pragma region Размеры системы по осям координат
const double LX = NUMCRIST_X * ACRIST;
const double LY = NUMCRIST_Y * ACRIST;
const double LZ = NUMCRIST_Z * ACRIST;
#pragma endregion

//LX = LY = LZ = ;
const double VOLUME = LX * LY * LZ; //объем системы
;