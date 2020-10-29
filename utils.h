#ifndef UTILS_H
#define UTILS_H

#include "point.h"
#include <gcrypt.h>

//Проверка точки на принадлежность кривой в форме Монтгомери
//Реализовано простой подстановкой значений коэф-ов и переменных в уравнение кривой
int isMontCurvePoint(struct montgomeryEllipticCurve* mec);




//
void doubleCurrentPoint(struct montgomeryEllipticCurve* currentPoint);


#endif
