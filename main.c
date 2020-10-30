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

    gcry_mpi_t two = gcry_mpi_new(0);
    gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, 0);
    
   
    //gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67", 0, 0);
    
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
    gcry_mpi_randomize(k, 23, GCRY_STRONG_RANDOM);
    createGostCurve512(c);
    //pprint(c->currPoint);
    montgomeryLadder(c, &c->currPoint, &k);
    printf("------ ");
    isMontCurvePoint(c);
    gcry_mpi_release(k);
    

    //ТЕСТ 3: проверить, что [q + 1][P] = p
    printf("ТЕСТ 3: проверить, что [q + 1][P] = P\n\n");
    struct point* testPoint1;
    testPoint1->x = gcry_mpi_new(0);
    testPoint1->x = gcry_mpi_copy(c->currPoint.x);
    testPoint1->y = gcry_mpi_new(0);
    testPoint1->y = gcry_mpi_copy(c->currPoint.y);
    testPoint1->z = gcry_mpi_new(0);
    testPoint1->z = gcry_mpi_copy(c->currPoint.z);
    createGostCurve512(c);
    //pprint(c->currPoint);
    gcry_mpi_addm(q, q, one, c->p);
    montgomeryLadder(c, testPoint1, &q);
    if(gcry_mpi_cmp(testPoint->x, c->currPoint.x) == 0){
        printf("----- ТЕСТ 3 завершился успешно");
    }
    
    // ТЕСТ 4: прочерить, что [q - 1][P] = -P
    createGostCurve512(c);
    testPoint1->x = gcry_mpi_new(0);
    testPoint1->x = gcry_mpi_copy(c->currPoint.x);
    testPoint1->y = gcry_mpi_new(0);
    testPoint1->y = gcry_mpi_copy(c->currPoint.y);
    testPoint1->z = gcry_mpi_new(0);
    testPoint1->z = gcry_mpi_copy(c->currPoint.z);
    //pprint(c->currPoint);
    gcry_mpi_subm(q, q, two, c->p);
    montgomeryLadder(c, testPoint1, &q);
    if(gcry_mpi_cmp(testPoint->x, c->currPoint.x) == 0){
        printf("----- ТЕСТ 4 завершился успешно");
    }

    // ТЕСТ 5: проверить аддитивность операции сложения [k1]P + [k2]P = [k1 + k2]P

    gcry_mpi_release(q);
    gcry_mpi_release(zero);
    gcry_mpi_release(one);
    gcry_mpi_release(testPoint1->z);
    gcry_mpi_release(testPoint1->x);
    gcry_mpi_release(testPoint1->y);
    gcry_mpi_release(c->currPoint.z);
    gcry_mpi_release(c->currPoint.x);
    gcry_mpi_release(c->currPoint.y);

    gcry_control(GCRYCTL_FINALIZE);
    
      

     
}
