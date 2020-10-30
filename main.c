#include "utils.h"
#include <gcrypt.h>
#include <stdio.h>

int main() {  

    struct montgomeryEllipticCurve* c;
    createGostCurve512(c);

    gcry_mpi_t q = gcry_mpi_new(0);
    gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC98CDBA46506AB004C33A9FF5147502CC8EDA9E7A769A12694623CEF47F023ED", 0, 0);

    gcry_mpi_t zero = gcry_mpi_new(0);
    gcry_mpi_scan(&zero, GCRYMPI_FMT_HEX, "0", 0, 0);

    gcry_mpi_t one = gcry_mpi_new(0);
    gcry_mpi_scan(&one, GCRYMPI_FMT_HEX, "1", 0, 0);
    
   
    //gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67", 0, 0);
    //pprint(c->currPoint);
    //int a = isMontCurvePoint(c);
    
    //doubleCurrentPoint(c);
    
    // ТЕСТ 1: проверить, что q[P] = 0
    printf("ТЕСТ 1: проверить, что q[P] = 0 \n\n");
    montgomeryLadder(c, &c->currPoint, &q);
    if(gcry_mpi_cmp(c->currPoint.z, zero) == 0 && gcry_mpi_cmp(c->currPoint.x, one) == 0){
        printf("------ ТЕСТ 1 завершился успешно\n\n");
    }
    //pprint(c->currPoint);

    // ТЕСТ 2: проверить, что k[P] принадлежит кривой
    printf("ТЕСТ 2: проверить, что k[P] принадлежит кривой\n\n");
    gcry_mpi_t k = gcry_mpi_new(0);
    gcry_mpi_randomize(k, 47, GCRY_STRONG_RANDOM);
    createGostCurve512(c);
    //pprint(c->currPoint);
    montgomeryLadder(c, &c->currPoint, &k);
    printf("------ ");
    isMontCurvePoint(c);
    gcry_mpi_release(k);
    

    //ТЕСТ 3: проверить, что [q + 1][P] = p
    printf("ТЕСТ 3: проверить, что [q + 1][P] = P\n\n");
    //struct point* testPoint1;
    //testPoint1->x = gcry_mpi_new(0);
    //testPoint1->x = gcry_mpi_copy(one);
    createGostCurve512(c);
    //pprint(c->currPoint);
    gcry_mpi_addm(q, q, one, c->p);
    montgomeryLadder(c, &c->currPoint, &q);
    
    
    
    // ТЕСТ 4: прочерить, что [q - 1][P] = -P
    // ТЕСТ 5: проверить аддитивность операции сложения [k1]P + [k2]P = [k1 + k2]P

    gcry_mpi_release(q);
    gcry_mpi_release(zero);
    //gcry_mpi_release(one);
    //gcry_mpi_release(testPoint1->z);
    //gcry_mpi_release(testPoint1->x);
    //gcry_mpi_release(testPoint1->y);
    gcry_mpi_release(c->currPoint.z);
    gcry_mpi_release(c->currPoint.x);
    gcry_mpi_release(c->currPoint.y);
    
      

     
}
