
#ifndef UNIT_TEST_H
#define UNIT_TEST_H

#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

class UnitTest
{
    public:
        UnitTest();

        void run();

        bool report(bool   result,
                    string message);
};
#endif


//~ Formatted by Jindent --- http://www.jindent.com
