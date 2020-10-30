/*************************
   * File: utils.h
   * Description: Содержит методы для проведения математических операций с точками на кривой Монтрогмери
   * Created: 29 oct 2020
   * Author: Роман Бесчастных

*************************/

#ifndef UTILS_H
#define UTILS_H

#include "point.h"
#include <gcrypt.h>

/**
 * @param mec - точка на кривой в форме Монтгомери
 * @return void
 * @note удваивает точку, полученную на входе. Результат записывается в неё же
**/
void doubleCurrentPoint(struct montgomeryEllipticCurve* currentPoint);



/**
 * @param firstPoint - первое слагаемое
 * @param secondPoint - второе слагаемое
 * @param initialPoint - начальная точка
 * @param p - модуль, по которому складываем
 * @return void
 * @note Складывает две точки и записывает результат в первый аргумент
**/
void sumPoints(struct point* firstPoint, struct point* secondPoint, struct point* initialPoint, gcry_mpi_t* p);



/**
 * @param mec - точка на кривой в форме Монтгомери
 * @param point - точка (тоже, что и mec, но без констант из ГОСТа)
 * @param k - целое число
 * @return void
 * @note Реализует алгоритм "лесенка Монтгомери" и приводит результат к афинным координатам
**/
void montgomeryLadder(struct montgomeryEllipticCurve* mec, struct point* point, gcry_mpi_t* k);



#endif
