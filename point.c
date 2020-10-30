#include "point.h"
#include <gcrypt.h>
#include <stdio.h>

void pprint(struct point p){
   gcry_mpi_dump(p.x);
   printf("\n X-axis point     \n\n");
   gcry_mpi_dump(p.y);
   printf("\n Y-axis point     \n\n");
   gcry_mpi_dump(p.z);
   printf("\n Z-axis point     \n\n");

};


void createGostCurve256(struct montgomeryEllipticCurve* mec){

   // модуль эллиптической кривой
   gcry_mpi_t p = gcry_mpi_new(0);
   gcry_mpi_scan(&p, GCRYMPI_FMT_HEX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97", 0, 0);
   
   // Коэффициенты в кривой формы Эдвардса - e, d
   gcry_mpi_t e = gcry_mpi_new(0);
   gcry_mpi_scan(&e, GCRYMPI_FMT_HEX, "1", 0, 0);
   
   gcry_mpi_t d = gcry_mpi_new(0);
   gcry_mpi_scan(&d, GCRYMPI_FMT_HEX, "605F6B7C183FA81578BC39CFAD518132B9DF62897009AF7E522C32D6DC7BFFB", 0, 0);
   
   //coordinaty
   gcry_mpi_t u = gcry_mpi_new(0);
   gcry_mpi_scan(&u, GCRYMPI_FMT_HEX, "D", 0, 0);
      
   gcry_mpi_t v = gcry_mpi_new(0);
   gcry_mpi_scan(&v, GCRYMPI_FMT_HEX, "60CA1E32AA475B348488C38FAB07649CE7EF8DBE87F22E81F92B2592DBA300E7", 0, 0);
   
   createAnyCurveByParameters(p, e, d, u, v, mec); 
};



void createGostCurve512(struct montgomeryEllipticCurve* mec){
   
   // модуль эллиптической кривой
   gcry_mpi_t p = gcry_mpi_new(0);
   gcry_mpi_scan(&p, GCRYMPI_FMT_HEX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC7", 0, 0);
   
   // Коэффициенты в кривой формы Эдвардса - e, d
   gcry_mpi_t e = gcry_mpi_new(0);
   gcry_mpi_scan(&e, GCRYMPI_FMT_HEX, "1", 0, 0);
   
   gcry_mpi_t d = gcry_mpi_new(0);
   gcry_mpi_scan(&d, GCRYMPI_FMT_HEX, "9E4F5D8C017D8D9F13A5CF3CDF5BFE4DAB402D54198E31EBDE28A0621050439CA6B39E0A515C06B304E2CE43E79E369E91A0CFC2BC2A22B4CA302DBB33EE7550", 0, 0);
   
   //coordinaty
   gcry_mpi_t u = gcry_mpi_new(0);
   gcry_mpi_scan(&u, GCRYMPI_FMT_HEX, "12", 0, 0);
      
   gcry_mpi_t v = gcry_mpi_new(0);
   gcry_mpi_scan(&v, GCRYMPI_FMT_HEX, "469AF79D1FB1F5E16B99592B77A01E2A0FDFB0D01794368D9A56117F7B38669522DD4B650CF789EEBF068C5D139732F0905622C04B2BAAE7600303EE73001A3D", 0, 0);
   
   createAnyCurveByParameters(p, e, d, u, v, mec);
};


//https://core.ac.uk/download/pdf/146445895.pdf - переход из формы Эдвардса в форму Монтгомери (используются утверждения 2.9 и 2.10, стр. 48)
void createAnyCurveByParameters(gcry_mpi_t p, gcry_mpi_t e, gcry_mpi_t d, gcry_mpi_t u, gcry_mpi_t v, struct montgomeryEllipticCurve* mec){

   gcry_mpi_t one = gcry_mpi_new(0);
   gcry_mpi_t two = gcry_mpi_new(0);   
   gcry_mpi_t four = gcry_mpi_new(0);
   gcry_mpi_t tmp = gcry_mpi_new(0);
   gcry_mpi_scan(&one, GCRYMPI_FMT_HEX, "1", 0, 0);
   gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, 0);   
   gcry_mpi_scan(&four, GCRYMPI_FMT_HEX, "4", 0, 0);
   //   
   mec->currPoint.x = gcry_mpi_new(0);
   mec->currPoint.y = gcry_mpi_new(0);
   mec->currPoint.z = gcry_mpi_new(0);
   mec->A = gcry_mpi_new(0);
   mec->B = gcry_mpi_new(0);    
   
   //вычислим точку X на кривой в форме Монтгомери: X = (1+v)/(1-v)
   gcry_mpi_addm(mec->currPoint.x, one, v, p); //1 + v
   gcry_mpi_subm(tmp, one, v, p); //1 - v
   gcry_mpi_invm(tmp, tmp, p); // 1/(1 - v)
   gcry_mpi_mulm(mec->currPoint.x, mec->currPoint.x, tmp, p); // итог X
   
      
   //вычислим точку Y на кривой в форме Монтгомери: Y = (1 + v)/u*(1 - v) 
   gcry_mpi_addm(mec->currPoint.y, one, v, p); // 1 + v 
   gcry_mpi_subm(tmp, one, v, p); // 1 - v 
   gcry_mpi_mulm(tmp, tmp, u, p); //u*(1 -v)
   gcry_mpi_invm(tmp, tmp, p); //(u*(1 -v))^(-1)
   gcry_mpi_mulm(mec->currPoint.y, mec->currPoint.y, tmp, p); //итог Y
   
   //вычислим точку Z на кривой в форме Монтгомери
   mec->p = p;
   mec->currPoint.z = one;  
   
   
   gcry_mpi_addm(mec->A, e, d, p); // e + d 
   gcry_mpi_mulm(mec->A, two, mec->A, p); // 2*(e + d)
   gcry_mpi_subm(tmp, e, d, p); // e - d 
   gcry_mpi_invm(tmp, tmp, p); // (e - d)^(-1)
   gcry_mpi_mulm(mec->A, mec->A, tmp, p); // 2*(e + d)/(e - d) 
 
   gcry_mpi_mulm(mec->B, four, tmp, p); // 4/(e - d)
   
   gcry_mpi_release(two);
   gcry_mpi_release(e);
   gcry_mpi_release(d);
   gcry_mpi_release(u);
   gcry_mpi_release(v);   
   gcry_mpi_release(four);
   gcry_mpi_release(tmp);       

};



int isMontCurvePoint(struct montgomeryEllipticCurve* mec){
   gcry_mpi_t leftPart = gcry_mpi_new(0);
   gcry_mpi_t rightPart = gcry_mpi_new(0);
   gcry_mpi_t tmp = gcry_mpi_new(0);

   // Вычисляем левую часть выражения в уравнении кривой Монтгомери
   gcry_mpi_mulm(leftPart, mec->currPoint.y, mec->currPoint.y, mec->p); // y*y
   gcry_mpi_mulm(leftPart, leftPart, mec->B, mec->p);                   // B*y^2
   
   // Вычисляем правую часть выражения в уравнении кривой Монтгомери
   gcry_mpi_mulm(rightPart, mec->currPoint.x, mec->currPoint.x, mec->p);// x^2
   gcry_mpi_mulm(rightPart, rightPart, mec->currPoint.x, mec->p);       // x^3
   gcry_mpi_mulm(tmp, mec->currPoint.x, mec->currPoint.x, mec->p);      // x^2
   gcry_mpi_mulm(tmp, tmp, mec->A, mec->p);                             // A*x^2
   gcry_mpi_addm(rightPart, rightPart, tmp, mec->p);                    // x^3+A*x^2
   gcry_mpi_addm(rightPart, rightPart, mec->currPoint.x, mec->p);       // x^3+A*x^2+x
   	
   if(gcry_mpi_cmp(leftPart, rightPart) == 0) {
	printf("Точка лежит на кривой\n\n");
	return 1;
   }
   printf("Точка не лежит на кривой\n\n");
   return 0;

   gcry_mpi_release(leftPart);
   gcry_mpi_release(rightPart);
   gcry_mpi_release(tmp);
};
