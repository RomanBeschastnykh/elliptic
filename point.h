#ifndef POINT_H
#define POINT_H

#include <gcrypt.h>

// точка
struct point {
    gcry_mpi_t x;
    gcry_mpi_t y;    
    gcry_mpi_t z;
};

// структура содержащая параметры для точки на кривой Монтгомери
struct montgomeryEllipticCurve {
    gcry_mpi_t A;
    gcry_mpi_t B;
    gcry_mpi_t p;
    struct point currPoint;
};

//вывод координат точки в стандратный вывод
void pprint(struct point p);



//Создание точки на эллиптической кривой вида Монтгомери
//Значения взяты из https://tc26.ru/standard/rs/Р 50.1.114-2016.pdf
//Использован набор параметров id-tc26-gost-3410-2012-256-paramSetA
struct montgomeryEllipticCurve* createGostCurve256();




//Использован набор параметров id-tc26-gost-3410-12-512-paramSetC
struct montgomeryEllipticCurve* createGostCurve512();



//Общий метод для создания точки на кривой
void createAnyCurveByParameters(gcry_mpi_t p, gcry_mpi_t e, gcry_mpi_t d, gcry_mpi_t q, gcry_mpi_t u, gcry_mpi_t v, struct montgomeryEllipticCurve* mec);

#endif

