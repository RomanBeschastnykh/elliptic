#include "utils.h"
#include <gcrypt.h>
#include <stdio.h>

int main() {  

    gcry_mpi_t q = gcry_mpi_new(0);
    gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC98CDBA46506AB004C33A9FF5147502CC8EDA9E7A769A12694623CEF47F023ED", 0, 0);

    //gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67", 0, 0);


    struct montgomeryEllipticCurve* c;
    createGostCurve512(c);
    //pprint(c->currPoint);
    //int a = isMontCurvePoint(c);
    
    //doubleCurrentPoint(c);
    
    // ТЕСТ 1: проверить, что q[P] = 0
    //montgomeryLadder(c, &c->currPoint, &q);
    //pprint(c->currPoint);

    // ТЕСТ 2: проверить, что [q + 1][P] = p
    
    // ТЕСТ 3: прочерить, что [q - 1][P] = -P
    // ТЕСТ 4: проверить аддитивность операции сложения [k1]P + [k2]P = [k1 + k2]P
    // ТЕСТ 5: проверить, что k[P] принадлежит кривой
      

     
}
