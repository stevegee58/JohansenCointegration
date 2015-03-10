/*
 * common_types.h
 *
 *  Created on: Jan 2, 2015
 *      Author: user
 */

#ifndef COMMON_TYPES_H_
#define COMMON_TYPES_H_

#include <string>
#include <vector>
using namespace std;

namespace CommonTypes
{
    typedef struct
    {
        double TestStatistic;
        double CriticalValue90;
        double CriticalValue95;
        double CriticalValue99;
    } MaxEigenData;

    typedef vector<double> DoubleVector;
    typedef vector<DoubleVector> DoubleMatrix;
}

#endif /* COMMON_TYPES_H_ */
