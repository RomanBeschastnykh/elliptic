#include <gcrypt.h>


// точка
struct point {
    gcry_mpi_t x;
    gcry_mpi_t y;    
    gcry_mpi_t z;
};

// структура содержащая параметры для порождения кривой Монтгомери
struct montgomeryEllipticCurve {
    gcry_mpi_t A;
    gcry_mpi_t B;
    gcry_mpi_t p;
    struct point currPoint;
};

struct montgomeryEllipticCurve createMontgomeryEllipticCurve() {
//использовался набор параметров id-tc26-gost-3410-2012-512-paramSetC
//источник: https://tc26.ru/standard/rs/Р 50.1.114-2016.pdf
//


   struct montgomeryEllipticCurve mec;

   // модуль эллиптической кривой
   gcry_mpi_t p = gcry_mpi_new(0);
   gcry_mpi_scan(&p, GCRYMPI_FMT_HEX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC7", 0, 0);
   
   // Коэффициенты в кривой формы Эдвардса - e, d
   gcry_mpi_t e = gcry_mpi_new(0);
   gcry_mpi_scan(&e, GCRYMPI_FMT_HEX, "1", 0, 0);
   
   gcry_mpi_t d = gcry_mpi_new(0);
   gcry_mpi_scan(&d, GCRYMPI_FMT_HEX, "9E4F5D8C017D8D9F13A5CF3CDF5BFE4DAB402D54198E31EBDE28A0621050439CA6B39E0A515C06B304E2CE43E79E369E91A0CFC2BC2A22B4CA302DBB33EE7550", 0, 0);
     
   //   
   gcry_mpi_t q = gcry_mpi_new(0);
   gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC98CDBA46506AB004C33A9FF5147502CC8EDA9E7A769A12694623CEF47F023ED", 0, 0);
   
   gcry_mpi_t u = gcry_mpi_new(0);
   gcry_mpi_scan(&u, GCRYMPI_FMT_HEX, "12", 0, 0);
      
   gcry_mpi_t v = gcry_mpi_new(0);
   gcry_mpi_scan(&v, GCRYMPI_FMT_HEX, "469AF79D1FB1F5E16B99592B77A01E2A0FDFB0D01794368D9A56117F7B38669522DD4B650CF789EEBF068C5D139732F0905622C04B2BAAE7600303EE73001A3D", 0, 0);
   
   //https://core.ac.uk/download/pdf/146445895.pdf - переход из формы Эдвардса в форму Монтгомери (используются утверждения 2.9 и 2.10, стр. 48)
   gcry_mpi_t one = gcry_mpi_new(0);
   gcry_mpi_t tmp = gcry_mpi_new(0);
   gcry_mpi_scan(&one, GCRYMPI_FMT_HEX, "1", 0, 0);
   
   //вычислим точку X на кривой в форме Монтгомери: X = (1+v)/(1-v)
   gcry_mpi_addm(mec.currPoint.x, one, v, p); //1 + v
   gcry_mpi_subm(tmp, one, v, p); //1 - v
   gcry_mpi_invm(tmp, tmp, p); // 1/(1 - v)
   gcry_mpi_mulm(mec.currPoint.x, mec.currPoint.x, tmp, p); // итог X
   
   //вычислим точку Y на кривой в форме Монтгомери: Y = (1 + v)/u*(1 - v) 
   gcry_mpi_addm(mec.currPoint.y, one, v, p); // 1 + v 
   gcry_mpi_subm(tmp, one, v, p); // 1 - v 
   gcry_mpi_mulm(tmp, tmp, u, p); //u*(1 -v)
   gcry_mpi_invm(tmp, tmp, p); //(u*(1 -v))^(-1)
   gcry_mpi_mulm(mec.currPoint.y, mec.currPoint.y, tmp, p); //итог Y
   
   //вычислим точку Z на кривой в форме Монтгомери
   mec.p = p;
   mec.currPoint.z = one;
   
   
   //вычислим параметры A и B
   gcry_mpi_t two = gcry_mpi_new(0);
   gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, 0);
   gcry_mpi_addm(mec.A, e, d, p); // e + d 
   gcry_mpi_mulm(mec.A, two, mec.A, p); // 2*(e + d)
   gcry_mpi_subm(tmp, e, d, p); // e - d 
   gcry_mpi_invm(tmp, tmp, p); // (e-d)^(-1)
   gcry_mpi_mulm(mec.A, mec.A, tmp, p); // итог для A = 2*(e + d)/(e - d)
   
   gcry_mpi_t four = gcry_mpi_new(0);
   gcry_mpi_scan(&four, GCRYMPI_FMT_HEX, "4", 0, 0);
   gcry_mpi_mulm(mec.B, four, tmp, p); // Итог для B = 4/(e - d)
   
   gcry_mpi_release(p);
   gcry_mpi_release(e);
   gcry_mpi_release(d);         
   gcry_mpi_release(q);
   gcry_mpi_release(u);
   gcry_mpi_release(v);
   gcry_mpi_release(one);   
   gcry_mpi_release(two);
   gcry_mpi_release(four);      
                  
   return(mec);         
}

//удвоение основано на алгоритме, предложенном в книге "Speeding the Pollard and elliptic curve methods of factorization", Montgomery, 1987, стр. 261 
// https://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866113-7/S0025-5718-1987-0866113-7.pdf
struct montgomeryEllipticCurve doubleCurrentPoint(struct montgomeryEllipticCurve currentPoint) {
   
   struct montgomeryEllipticCurve doubledPoint;
   
   gcry_mpi_t sum = gcry_mpi_new(0);
   gcry_mpi_t diff = gcry_mpi_new(0);
   gcry_mpi_t two = gcry_mpi_new(0);
   gcry_mpi_t four = gcry_mpi_new(0);
   gcry_mpi_t coeff = gcry_mpi_new(0);

   gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, 0);
   gcry_mpi_scan(&four, GCRYMPI_FMT_HEX, "4", 0, 0);


   gcry_mpi_addm(sum, currentPoint.currPoint.x, currentPoint.currPoint.z, currentPoint.p); // X1 + Z1
   gcry_mpi_mulm(sum, sum, sum, currentPoint.p);  // (X1 + Z1)^2
   gcry_mpi_subm(diff, currentPoint.currPoint.x, currentPoint.currPoint.z, currentPoint.p); // X1 - Z1
   gcry_mpi_mulm(diff, diff, diff, currentPoint.p);  // (X1 - Z1)^2
   gcry_mpi_mulm(doubledPoint.currPoint.x, sum, diff, currentPoint.p);// X3 = (X1 + Z1)^2 * (X1 - Z1)^2
   gcry_mpi_subm(sum, sum, diff, currentPoint.p); // (X1 + Z1)^2 - (X1 - Z1)^2
   gcry_mpi_invm(four, four, currentPoint.p); // 4^(-1)
   gcry_mpi_addm(coeff, currentPoint.A, two, currentPoint.p); // a + 2
   gcry_mpi_mulm(coeff, coeff, four, currentPoint.p); // (a+2)/4
   gcry_mpi_mulm(coeff, coeff, sum, currentPoint.p);// (a+2)/4 * ((X1 + Z1)^2 - (X1 - Z1)^2)
   gcry_mpi_addm(coeff, diff, coeff, currentPoint.p); // (X1 - Z1)^2 + (a+2)/4 * ((X1 + Z1)^2 - (X1 - Z1)^2)
   gcry_mpi_mulm(doubledPoint.currPoint.z, sum, coeff, currentPoint.p); // ((X1 + Z1)^2 - (X1 - Z1)^2) * ((X1 - Z1)^2 + (a+2)/4 * ((X1 + Z1)^2 - (X1 - Z1)^2))
   
   gcry_mpi_release(sum);
   gcry_mpi_release(diff);
   gcry_mpi_release(one);
   gcry_mpi_release(two);
   gcry_mpi_release(four);
   gcry_mpi_release(coeff);               
   
   return(doubledPoint);
   
} 

//transform point from plain point to affinian

//check if point is on curve

//add point to point

//Montgomery ladder


