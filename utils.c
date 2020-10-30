#include <gcrypt.h>
#include "utils.h"

void doubleCurrentPoint(struct montgomeryEllipticCurve* currentPoint) {
   
   gcry_mpi_t summ = gcry_mpi_new(0);
   gcry_mpi_t diff = gcry_mpi_new(0);
   gcry_mpi_t two = gcry_mpi_new(0);
   gcry_mpi_t four = gcry_mpi_new(0);
   gcry_mpi_t coeff = gcry_mpi_new(0);

   gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, 0);
   gcry_mpi_scan(&four, GCRYMPI_FMT_HEX, "4", 0, 0);
 	

   // Вычисляем значение координаты X1 = (X+Z)^2 * (X-Z)^2
   gcry_mpi_addm(summ, currentPoint->currPoint.x, currentPoint->currPoint.z, currentPoint->p);
   gcry_mpi_mulm(summ, summ, summ, currentPoint->p);
   gcry_mpi_subm(diff, currentPoint->currPoint.x, currentPoint->currPoint.z, currentPoint->p);
   gcry_mpi_mulm(diff, diff, diff, currentPoint->p);
   gcry_mpi_mulm(currentPoint->currPoint.x, summ, diff, currentPoint->p);


   // Вычисляем значение координаты Y1 = ((X+Z)^2 - (X-Z)^2)*((X-Z)^2 + (A+2)/4 * ((X+Z)^2 - (X-Z)^2))
   gcry_mpi_subm(summ, summ, diff, currentPoint->p);                      // (X+Z)^2 - (X-Z)^2)
   gcry_mpi_invm(four, four,  currentPoint->p);                           // 4^(-1)
   gcry_mpi_add(coeff, currentPoint->A, two);                             // A+2 
   gcry_mpi_mulm(coeff, coeff, four, currentPoint->p);                    // (A+2)/4
   gcry_mpi_mulm(coeff, coeff, summ, currentPoint->p);                    // (A+2)/4 * ((X+Z)^2 - (X-Z)^2)
   gcry_mpi_addm(coeff, diff, coeff, currentPoint->p);                    // (X-Z)^2 + (A+2)/4 * ((X+Z)^2 - (X-Z)^2)
   gcry_mpi_mulm(currentPoint->currPoint.z, summ, coeff, currentPoint->p);// итог

   gcry_mpi_release(summ);
   gcry_mpi_release(diff);	
   gcry_mpi_release(coeff);
   gcry_mpi_release(two);
   gcry_mpi_release(four);               
   
};

void sumPoints(struct point* firstPoint, struct point* secondPoint, struct point* initialPoint, gcry_mpi_t* p){ 

    gcry_mpi_t fDiff = gcry_mpi_new(0);
    gcry_mpi_t sSumm = gcry_mpi_new(0);
    gcry_mpi_t fSumm = gcry_mpi_new(0);
    gcry_mpi_t sDiff = gcry_mpi_new(0);
 
    gcry_mpi_subm(fDiff, firstPoint->x, firstPoint->z, *p);
    gcry_mpi_addm(sSumm, secondPoint->x, secondPoint->z, *p);
    gcry_mpi_addm(fSumm, firstPoint->x, firstPoint->z, *p);
    gcry_mpi_subm(sDiff, secondPoint->x, secondPoint->z, *p);
    gcry_mpi_mulm(fDiff, fDiff, sSumm, *p);
    gcry_mpi_mulm(fSumm, fSumm, sDiff, *p);

    gcry_mpi_addm(sSumm, fDiff, fSumm, *p);                   // (X2-Z2)*(X3+Z3) + (X3-Z3)*(X2+Z2)
    gcry_mpi_mulm(sSumm, sSumm, sSumm, *p);                   // ((X2-Z2)*(X3+Z3) + (X3-Z3)*(X2+Z2))^2
    gcry_mpi_mulm(firstPoint->x, initialPoint->z, sSumm, *p); // Z1 * ((X2-Z2)*(X3+Z3) + (X3-Z3)*(X2+Z2))^2
		
    gcry_mpi_subm(sSumm, fDiff, fSumm, *p);                   // (X2-Z2)*(X3+Z3) - (X3-Z3)*(X2+Z2)
    gcry_mpi_mulm(sSumm, sSumm, sSumm, *p);                   // ((X2-Z2)*(X3+Z3) - (X3-Z3)*(X2+Z2))^2
    gcry_mpi_mulm(firstPoint->z, initialPoint->x, sSumm, *p); // X1 * ((X2-Z2)*(X3+Z3) - (X3-Z3)*(X2+Z2))^2

    gcry_mpi_release(fSumm);
    gcry_mpi_release(sSumm);
    gcry_mpi_release(fDiff);
    gcry_mpi_release(sDiff);

}

void montgomeryLadder(struct montgomeryEllipticCurve* mec, struct point* point, gcry_mpi_t* k){

    gcry_mpi_t inverted = gcry_mpi_new(0);
    gcry_mpi_t zero = gcry_mpi_new(0);
    gcry_mpi_scan(&zero, GCRYMPI_FMT_HEX, "0", 0, 0);

    unsigned int bits = gcry_mpi_get_nbits(*k);

    //копия входной точки на кривой
    struct montgomeryEllipticCurve r;

    r.p = gcry_mpi_copy(mec->p);
    r.A = gcry_mpi_copy(mec->A);
    r.B = gcry_mpi_copy(mec->B);
    r.currPoint.x = gcry_mpi_copy(point->x);
    r.currPoint.y = gcry_mpi_new(0);
    r.currPoint.z = gcry_mpi_copy(point->z);

    //Нейтральный элемент
    struct montgomeryEllipticCurve q;

    q.p = gcry_mpi_copy(mec->p);
    q.A = gcry_mpi_copy(mec->A);
    q.B = gcry_mpi_copy(mec->B);
    gcry_mpi_scan(&q.currPoint.x, GCRYMPI_FMT_HEX, "1", 0, 0);
    q.currPoint.y = gcry_mpi_new(0);
    gcry_mpi_scan(&q.currPoint.z, GCRYMPI_FMT_HEX, "0", 0, 0);

    for(int i = bits - 1; i >= 0; i--){
        if(gcry_mpi_test_bit(*k, i) == 1){
	    sumPoints(&q.currPoint, &r.currPoint, point, &mec->p);
	    doubleCurrentPoint(&r);
	}else{
	    sumPoints(&r.currPoint, &q.currPoint, point, &mec->p);
	    doubleCurrentPoint(&q);
	}
    }

    point->x = gcry_mpi_copy(q.currPoint.x);
    point->y = gcry_mpi_copy(q.currPoint.y);
    point->z = gcry_mpi_copy(q.currPoint.z);

    //Проверка на ноль, чтобы на него не поделить случайно
    if(gcry_mpi_cmp(point->z, zero) != 0){
        gcry_mpi_invm(inverted, point->z, mec->p);
	gcry_mpi_mulm(point->x, point->x, inverted, mec->p);
	gcry_mpi_mulm(point->y, point->y, inverted, mec->p);
	gcry_mpi_mulm(point->z, point->z, inverted, mec->p);
    }else{
	gcry_mpi_invm(inverted, point->x, mec->p);
	gcry_mpi_mulm(point->x, point->x, inverted, mec->p);
    }

    gcry_mpi_release(inverted);  
    gcry_mpi_release(zero);
    gcry_mpi_release(r.currPoint.x);
    gcry_mpi_release(r.currPoint.y);
    gcry_mpi_release(r.currPoint.z);
    gcry_mpi_release(r.A);
    gcry_mpi_release(r.B);
    gcry_mpi_release(r.p);
    gcry_mpi_release(q.currPoint.x);
    gcry_mpi_release(q.currPoint.y);
    gcry_mpi_release(q.currPoint.z);
    gcry_mpi_release(q.A);
    gcry_mpi_release(q.B);
    gcry_mpi_release(q.p);

}

