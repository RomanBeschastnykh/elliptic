#include "utils.h"
#include <gcrypt.h>
#include <stdio.h>

int main() {  

    struct montgomeryEllipticCurve c;
    createGostCurve256(&c);

    gcry_mpi_t q = gcry_mpi_new(0);
    //gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC98CDBA46506AB004C33A9FF5147502CC8EDA9E7A769A12694623CEF47F023ED", 0, 0);
    gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67", 0, 0);

    gcry_mpi_t zero = gcry_mpi_new(0);
    gcry_mpi_scan(&zero, GCRYMPI_FMT_HEX, "0", 0, 0);

    gcry_mpi_t one = gcry_mpi_new(0);
    gcry_mpi_scan(&one, GCRYMPI_FMT_HEX, "1", 0, 0);

    gcry_mpi_t two = gcry_mpi_new(0);
    gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, 0);
    
   



    
    // ТЕСТ 1: проверить, что q[P] = 0
    printf("ТЕСТ 1: проверить, что q[P] = 0 \n\n");
    struct point testPoint1;
    testPoint1.x = gcry_mpi_new(0);
    testPoint1.x = gcry_mpi_copy(c.currPoint.x);
    testPoint1.y = gcry_mpi_new(0);
    testPoint1.y = gcry_mpi_copy(c.currPoint.y);
    testPoint1.z = gcry_mpi_new(0);
    testPoint1.z = gcry_mpi_copy(c.currPoint.z);
    montgomeryLadder(&c, &testPoint1, &q);
    if(gcry_mpi_cmp(testPoint1.z, zero) == 0 && gcry_mpi_cmp(testPoint1.x, one) == 0){
        printf("------ ТЕСТ 1 завершился успешно\n\n");
    }




    // ТЕСТ 2: проверить, что k[P] принадлежит кривой
    printf("ТЕСТ 2: проверить, что k[P] принадлежит кривой\n\n");
    gcry_mpi_t k = gcry_mpi_new(0);
    gcry_mpi_randomize(k, 23, GCRY_STRONG_RANDOM);
    createGostCurve256(&c);
    testPoint1.x = gcry_mpi_new(0);
    testPoint1.x = gcry_mpi_copy(c.currPoint.x);
    testPoint1.y = gcry_mpi_new(0);
    testPoint1.y = gcry_mpi_copy(c.currPoint.y);
    testPoint1.z = gcry_mpi_new(0);
    testPoint1.z = gcry_mpi_copy(c.currPoint.z);
    montgomeryLadder(&c, &testPoint1, &k);
    printf("------ ");
    isMontCurvePoint(&c);



    

    //ТЕСТ 3: проверить, что [q + 1][P] = p
    printf("ТЕСТ 3: проверить, что [q + 1][P] = P\n\n");
    testPoint1.x = gcry_mpi_new(0);
    testPoint1.x = gcry_mpi_copy(c.currPoint.x);
    testPoint1.y = gcry_mpi_new(0);
    testPoint1.y = gcry_mpi_copy(c.currPoint.y);
    testPoint1.z = gcry_mpi_new(0);
    testPoint1.z = gcry_mpi_copy(c.currPoint.z);
    createGostCurve256(&c);
    //pprint(testPoint1);
    gcry_mpi_addm(q, q, one, c.p);
    montgomeryLadder(&c, &testPoint1, &q);
    if(gcry_mpi_cmp(testPoint1.x, c.currPoint.x) == 0){
        printf("------ ТЕСТ 3 завершился успешно\n\n");
    }







    
    //ТЕСТ 4: прочерить, что [q - 1][P] = -P
    printf("ТЕСТ 4: проверить, что [q - 1][P] = -P\n\n");
    createGostCurve256(&c);
    testPoint1.x = gcry_mpi_new(0);
    testPoint1.x = gcry_mpi_copy(c.currPoint.x);
    testPoint1.y = gcry_mpi_new(0);
    testPoint1.y = gcry_mpi_copy(c.currPoint.y);
    testPoint1.z = gcry_mpi_new(0);
    testPoint1.z = gcry_mpi_copy(c.currPoint.z);
    gcry_mpi_subm(q, q, two, c.p); // После предыдущего теста значение стало q + 1, поэтому надо вычесть двойку
    montgomeryLadder(&c, &testPoint1, &q);
    gcry_mpi_dump(testPoint1.x);
    printf("qqqqq");
    if(gcry_mpi_cmp(testPoint1.x, c.currPoint.x) == 0){
         printf("------ ТЕСТ 4 завершился успешно\n\n");
    }









    // ТЕСТ 5: проверить аддитивность операции сложения [k1]P + [k2]P = [k1 + k2]P
    printf("ТЕСТ 5: проверить аддитивность операции сложения [k1]P + [k2]P = [k1 + k2]P\n\n");
    
    createGostCurve256(&c);
    
    //k1
    gcry_mpi_t k1 = gcry_mpi_new(0);
    gcry_mpi_randomize(k1, 32, GCRY_STRONG_RANDOM);

    testPoint1.x = gcry_mpi_new(0);
    testPoint1.x = gcry_mpi_copy(c.currPoint.x);
    testPoint1.y = gcry_mpi_new(0);
    testPoint1.y = gcry_mpi_copy(c.currPoint.y);
    testPoint1.z = gcry_mpi_new(0);
    testPoint1.z = gcry_mpi_copy(c.currPoint.z);

    montgomeryLadder(&c, &testPoint1, &k1);


    //k2
    struct point testPoint2;

    testPoint2.x = gcry_mpi_new(0);
    testPoint2.x = gcry_mpi_copy(c.currPoint.x);
    testPoint2.y = gcry_mpi_new(0);
    testPoint2.y = gcry_mpi_copy(c.currPoint.y);
    testPoint2.z = gcry_mpi_new(0);
    testPoint2.z = gcry_mpi_copy(c.currPoint.z);

    gcry_mpi_t k2 = gcry_mpi_new(0);
    gcry_mpi_randomize(k2, 32, GCRY_STRONG_RANDOM);

    montgomeryLadder(&c, &testPoint2, &k2);
    
    if(gcry_mpi_cmp(k1, k2) < 0){
	gcry_mpi_swap(k1,k2);
    }   


    //k1+k2
    struct point testPoint3;

    testPoint3.x = gcry_mpi_new(0);
    testPoint3.x = gcry_mpi_copy(c.currPoint.x);
    testPoint3.y = gcry_mpi_new(0);
    testPoint3.y = gcry_mpi_copy(c.currPoint.y);
    testPoint3.z = gcry_mpi_new(0);
    testPoint3.z = gcry_mpi_copy(c.currPoint.z);

    gcry_mpi_t k1k2 = gcry_mpi_new(0);
    gcry_mpi_addm(k1k2, k1, k2, c.p);

    montgomeryLadder(&c, &testPoint3, &k1k2);

    gcry_mpi_t inverted = gcry_mpi_new(0);


    //k1-k2
    struct point testPoint4;

    testPoint4.x = gcry_mpi_new(0);
    testPoint4.x = gcry_mpi_copy(c.currPoint.x);
    testPoint4.y = gcry_mpi_new(0);
    testPoint4.y = gcry_mpi_copy(c.currPoint.y);
    testPoint4.z = gcry_mpi_new(0);
    testPoint4.z = gcry_mpi_copy(c.currPoint.z);

    gcry_mpi_t mink1k2 = gcry_mpi_new(0);
    gcry_mpi_subm(mink1k2, k1, k2, c.p);

    montgomeryLadder(&c, &testPoint4, &mink1k2);

    sumPoints(&testPoint1, &testPoint2, &testPoint4, &c.p);  

    if(gcry_mpi_cmp(testPoint1.z, zero) != 0){
        gcry_mpi_invm(inverted, testPoint1.z, c.p);
	gcry_mpi_mulm(testPoint1.x, testPoint1.x, inverted, c.p);
	gcry_mpi_mulm(testPoint1.y, testPoint1.y, inverted, c.p);
	gcry_mpi_mulm(testPoint1.z, testPoint1.z, inverted, c.p);
    }else{
	gcry_mpi_invm(inverted, testPoint1.x, c.p);
	gcry_mpi_mulm(testPoint1.x, testPoint1.x, inverted, c.p);
    }

    if( (gcry_mpi_cmp(testPoint1.x, testPoint3.x) == 0)  && (gcry_mpi_cmp(testPoint1.z, testPoint3.z) == 0)){
	printf("------ ТЕСТ 5 завершился успешно\n\n");
    }















    gcry_mpi_release(q);
    gcry_mpi_release(zero);
    gcry_mpi_release(one);
    gcry_mpi_release(testPoint1.z);
    gcry_mpi_release(testPoint1.x);
    gcry_mpi_release(testPoint1.y);
    gcry_mpi_release(testPoint2.z);
    gcry_mpi_release(testPoint2.x);
    gcry_mpi_release(testPoint2.y);
    gcry_mpi_release(testPoint3.z);
    gcry_mpi_release(testPoint3.x);
    gcry_mpi_release(testPoint3.y);
    gcry_mpi_release(c.currPoint.z);
    gcry_mpi_release(c.currPoint.x);
    gcry_mpi_release(c.currPoint.y);
    gcry_mpi_release(k1);
    gcry_mpi_release(k2);
    gcry_mpi_release(k1k2);
    gcry_mpi_release(mink1k2);
    gcry_mpi_release(k);
     
}
