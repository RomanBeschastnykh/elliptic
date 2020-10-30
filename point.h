/*************************
   * File: point.h
   * Description: Содержит структуры для простой точки, точки на кривой и базовые методы создания и ввода-вывода, метод проверки принадлежности точки кривой
   * Created: 28 oct 2020
   * Author: Роман Бесчастных

*************************/

#ifndef POINT_H
#define POINT_H

#include <gcrypt.h>

/**
 * @note структура, представляющая точку
**/
struct point {
    gcry_mpi_t x;
    gcry_mpi_t y;    
    gcry_mpi_t z;
};

/**
 * @note структура содержащая параметры для точки на кривой Монтгомери
**/
struct montgomeryEllipticCurve {
    gcry_mpi_t A;
    gcry_mpi_t B;
    gcry_mpi_t p;
    struct point currPoint;
};

/**
 * @param p - точка
 * @return void
 * @note выводит координаты точки в базовый ввыод
**/
void pprint(struct point p);


//Значения взяты из https://tc26.ru/standard/rs/Р 50.1.114-2016.pdf
//Использован набор параметров id-tc26-gost-3410-2012-256-paramSetA
/**
 * @param mec - точка на кривой в форме Монтгомери
 * @return void
 * @note заполняет параметр mec константами из ГОСТ для 256 бит
**/
void createGostCurve256(struct montgomeryEllipticCurve* mec);



/**
 * @param mec - точка на кривой в форме Монтгомери
 * @return void
 * @note заполняет параметр mec константами из ГОСТ для 512 бит
**/
//Использован набор параметров id-tc26-gost-3410-12-512-paramSetC
void createGostCurve512(struct montgomeryEllipticCurve* mec);



/**
 * @param p - точка на кривой в форме Монтгомери
 * @param e - точка на кривой в форме Монтгомери
 * @param d - точка на кривой в форме Монтгомери
 * @param u - точка на кривой в форме Монтгомери
 * @param v - точка на кривой в форме Монтгомери
 * @return void
 * @note заполняет параметр mec константами p, e, d, u, v
**/
void createAnyCurveByParameters(gcry_mpi_t* p, gcry_mpi_t* e, gcry_mpi_t* d, gcry_mpi_t* u, gcry_mpi_t* v, struct montgomeryEllipticCurve* mec);




/**
 * @param mec - точка на кривой в форме Монтгомери
 * @return 1, если точка лежит на кривой и 0 в противном случае 
 * @note проверят, лежит ли точка на кривой в форме Монтгомери
**/
int isMontCurvePoint(struct montgomeryEllipticCurve* mec);

#endif

