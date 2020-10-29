#include "utils.h"
#include <gcrypt.h>
#include <stdio.h>

int main() {
    struct montgomeryEllipticCurve* c = createGostCurve512();
    pprint(c->currPoint);
    int a = isMontCurvePoint(c);
    
    doubleCurrentPoint(c);
    pprint(c->currPoint);    
}
