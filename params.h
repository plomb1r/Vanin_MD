#pragma once
//params.h - параметры (константы), значения которых можем менять от запуска к запуску
const int NSTEPS = 1000; //число шагов
const int LASTSTEP = NSTEPS - 1; //последний шаг

#pragma region Число элементарных ячеек (кристаллов) по осям координат
const int NUMCRIST_X = 1; 
const int NUMCRIST_Y = 1;
const int NUMCRIST_Z = 1;
#pragma endregion

//число частиц для примитивной решетки с учетом ПГУ
//const int NUMBERPARTICLES = NUMCRIST_X * NUMCRIST_Y * NUMCRIST_Z;
//NUMBERPARTICLES = ;
const int NUMBERPARTICLES = 2;
const double STEP = 0.002 // шаг интегрирования разностной схемы

#pragma region Пояснения
/*ACRIST - длина ребра элементарной ячейки, зависит от термодинамического состояни и модели
* ПГУ - периодиечские граничные условия
*/
#pragma endregion
;