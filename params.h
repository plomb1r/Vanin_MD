#pragma once
//params.h - параметры (константы), значения которых можем менять от запуска к запуску
const int NSTEPS = 4; //число шагов
const int LASTSTEP = NSTEPS - 1; //последний шаг

#pragma region Число элементарных ячеек (кристаллов) по осям координат
const int NUMCRIST_X = 4; 
const int NUMCRIST_Y = 4;
const int NUMCRIST_Z = 4;
#pragma endregion

//число частиц для примитивной решетки с учетом ПГУ
//const int NUMBERPARTICLES = NUMCRIST_X * NUMCRIST_Y * NUMCRIST_Z;
//NUMBERPARTICLES = ;
const int NUMBERPARTICLES = 2;
const double STEP = 0.002 // шаг интегрирования разностной схемы

#pragma region Пояснения
/*ACRIST - длина ребера элементарной ячейки, зависит от термодинамического состояни и модели
* ПГУ - периодиечские граничные условия
*/
#pragma endregion
;