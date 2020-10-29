#include <gcrypt.h>
#include "utils.h"


int isMontCurvePoint(struct montgomeryEllipticCurve* mec){
   gcry_mpi_t leftPart = gcry_mpi_new(0);
   gcry_mpi_t rightPart = gcry_mpi_new(0);
   gcry_mpi_t tmp = gcry_mpi_new(0);

   // Вычисляем левую часть выражения в уравнении кривой Монтгомери
   gcry_mpi_mulm(leftPart, mec->currPoint.y, mec->currPoint.y, mec->p); // y*y
   gcry_mpi_mulm(leftPart, leftPart, mec->B, mec->p); // B*y^2
   
   // Вычисляем правую часть выражения в уравнении кривой Монтгомери
   gcry_mpi_mulm(rightPart, mec->currPoint.x, mec->currPoint.x, mec->p); // x^2
   gcry_mpi_mulm(rightPart, rightPart, mec->currPoint.x, mec->p); // x^3
   gcry_mpi_mulm(tmp, mec->currPoint.x, mec->currPoint.x, mec->p); // x^2
   gcry_mpi_mulm(tmp, tmp, mec->A, mec->p); // A*x^2
   gcry_mpi_addm(rightPart, rightPart, tmp, mec->p); // x^3+A*x^2
   gcry_mpi_addm(rightPart, rightPart, mec->currPoint.x, mec->p); // x^3+A*x^2+x
   	
   if(gcry_mpi_cmp(leftPart, rightPart) == 0) {
	printf("On curve\n");
	return 1;
   }
   printf("Not on curve\n");
   return 0;
};

//удвоение основано на алгоритме, предложенном в статье "Speeding the Pollard and elliptic curve methods of factorization", Montgomery, 1987, стр. 261 
// https://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866113-7/S0025-5718-1987-0866113-7.pdf
void doubleCurrentPoint(struct montgomeryEllipticCurve* currentPoint) {
   
   gcry_mpi_t squareX = gcry_mpi_new(0);
   gcry_mpi_t squareZ = gcry_mpi_new(0);
   gcry_mpi_t diff = gcry_mpi_new(0);
   gcry_mpi_t tmp = gcry_mpi_new(0);
   gcry_mpi_t coeff = gcry_mpi_new(0);   
   gcry_mpi_t four = gcry_mpi_new(0);
   gcry_mpi_t multXZ = gcry_mpi_new(0);

   gcry_mpi_scan(&four, GCRYMPI_FMT_HEX, "4", 0, 0);

   //
   gcry_mpi_mulm(squareX, currentPoint->currPoint.x, currentPoint->currPoint.x, currentPoint->p); // X1 ^ 2
   gcry_mpi_mulm(squareZ, currentPoint->currPoint.z, currentPoint->currPoint.z, currentPoint->p);  // Z1 ^ 2
   gcry_mpi_subm(diff, squareX, squareZ, currentPoint->p); // X1^2 - Z1^2
   gcry_mpi_mulm(currentPoint->currPoint.x, diff, diff, currentPoint->p);  // (X1^2 - Z1^2)^2
   
   //
   gcry_mpi_mulm(multXZ, currentPoint->currPoint.x, currentPoint->currPoint.z, currentPoint->p); // X1 * Z1
   gcry_mpi_mulm(currentPoint->currPoint.z, four, multXZ, currentPoint->p); // 4 * X1 * Z1
   gcry_mpi_addm(tmp, squareX, squareZ, currentPoint->p); // X1^2 + Z1^2
   gcry_mpi_mulm(coeff, multXZ, currentPoint->A, currentPoint->p);// A * X1 * Z1
   gcry_mpi_addm(tmp, tmp, coeff, currentPoint->p); // X1^2 + A * X1 * Z1 + Z1^2 
   gcry_mpi_mulm(currentPoint->currPoint.z, currentPoint->currPoint.z, coeff, currentPoint->p); // 4 * X1 * Z1 * (X1^2 + A * X1 * Z1 + Z1^2)
   
   gcry_mpi_release(multXZ);
   gcry_mpi_release(diff);
   gcry_mpi_release(four);
   gcry_mpi_release(squareZ);
   gcry_mpi_release(squareX);
   gcry_mpi_release(coeff);
   gcry_mpi_release(tmp);               
   
}; 

//transform point from plain point to affinian

//add point to point

//Montgomery ladder


