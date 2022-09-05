#ifndef JOHANSEN_HELPER_DEFINED
#define JOHANSEN_HELPER_DEFINED

#include "CommonTypes.h"
#include <gsl/gsl_matrix.h>

using namespace std;
using namespace CommonTypes;

class JohansenHelper
{
public:
    JohansenHelper(DoubleMatrix xMat);
    ~JohansenHelper();

    void DoMaxEigenValueTest(int nlags);
    int CointegrationCount();
    vector<MaxEigenData> GetOutStats();
    DoubleVector GetEigenValues();
    DoubleMatrix GetEigenVecMatrix();

private:
    // Data members
    gsl_matrix* xMat_gsl;
    vector<MaxEigenData> outStats;
    DoubleVector eigenValuesVec;
    DoubleMatrix eigenVecMatrix;

    // Methods
    gsl_matrix* GetSubMatrix(gsl_matrix* xMat,
        int nBeginRow, int nEndRow,
        int nBeginCol, int nEndCol);
        gsl_matrix* GetMatrixDifference(gsl_matrix* xMat);
    gsl_matrix* DeMean(gsl_matrix* xMat);
    gsl_matrix* GetMatrixLagged(gsl_matrix* xMat, int nlags);
    double GetAverage(gsl_vector* x);
    void MatrixDivideByElem(gsl_matrix* xMat, double val);
    gsl_matrix* MatrixTransposeImpl(gsl_matrix* m);
    gsl_matrix* MatrixMultiply(gsl_matrix* A, gsl_matrix* B);
    gsl_matrix* GetMatrixInverse(gsl_matrix* inMat);
    gsl_matrix* MatrixDivide(gsl_matrix* xMat, gsl_matrix* yMat);
    gsl_matrix* DoubleMatrixToGSLMatrix(DoubleMatrix doubleMat);
    gsl_vector* DoubleVectorToGSLVector(DoubleVector doubleVec);
    DoubleMatrix GSLMatrixToDoubleMatrix(gsl_matrix* gslMat);
    DoubleVector GSLVectorToDoubleVector(gsl_vector* gslMat);
    DoubleVector GSLComplexVecToAbsDoubleVector(gsl_vector_complex* gslVec);
    DoubleMatrix GSLComplexMatToAbsDoubleMatrix(gsl_matrix_complex* gslMat);
};
#endif
