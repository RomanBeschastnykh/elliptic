#ifndef UTILS_H
#define UTILS_H

#include "point.h"
#include <gcrypt.h>

// Метод для удвоения текущей точки
void doubleCurrentPoint(struct montgomeryEllipticCurve* currentPoint);



// Метод для сложения двух точек
void sumPoints(struct point* firstPoint, struct point* secondPoint, struct point* initialPoint, gcry_mpi_t* p);



// 
void montgomeryLadder(struct montgomeryEllipticCurve* mec, struct point* point, gcry_mpi_t* k);


#endif
