#include "JohansenHelper.h"
#include "gsl/gsl_eigen.h"
#include <math.h>
#include "gsl/gsl_math.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_linalg.h"

namespace
{
    double maxEigenCriticalvalues[][3] =
        {
        {2.7055,    3.8415,     6.6349},
        {12.2971,   14.2639,    18.52},
        {18.8928,   21.1314,    25.865},
        {25.1236,   27.5858,    32.7172},
        {31.2379,   33.8777,    39.3693},
        {37.2786,   40.0763,    45.8662},
        {43.2947,   46.2299,    52.3069},
        {49.2855,   52.3622,    58.6634},
        {55.2412,   58.4332,    64.996},
        {61.2041,   64.504,     71.2525},
        {67.1307,   70.5392,    77.4877},
        {73.0563,   76.5734,    83.7105}
        };
}

JohansenHelper::JohansenHelper(DoubleMatrix xMat)
{
    // Convert input matrix to gsl_matrix
    xMat_gsl = DoubleMatrixToGSLMatrix(xMat);
}

JohansenHelper::~JohansenHelper()
{
    // No processing required.
}

vector<MaxEigenData> JohansenHelper::GetOutStats()
{
    return outStats;
}

DoubleVector JohansenHelper::GetEigenValues()
{
    return eigenValuesVec;
}

DoubleMatrix JohansenHelper::GetEigenVecMatrix()
{
    return eigenVecMatrix;
}

int JohansenHelper::CointegrationCount()
{
    int retVal = 0;
    int foundCointegrations = 0;
    bool conclusive = false;
    for (int i = 0; i < (int)outStats.size(); i++)
    {
        if (conclusive == false)
        {
            if ((outStats[i].TestStatistic < outStats[i].CriticalValue90) &&
                (outStats[i].TestStatistic < outStats[i].CriticalValue95) &&
                (outStats[i].TestStatistic < outStats[i].CriticalValue99))
            {
                foundCointegrations = i;
                conclusive = true;
            }
        }
    }

    if (conclusive == false)
    {
        retVal = -1;
    }
    else
    {
        retVal = foundCointegrations;
    }

    return retVal;
}

void JohansenHelper::DoMaxEigenValueTest(int nlags)
{
    // Demean input data in place
    gsl_matrix* xMat_temp = NULL;
    xMat_temp = DeMean(xMat_gsl);
    gsl_matrix_memcpy(xMat_gsl, xMat_temp);
    gsl_matrix_free(xMat_temp);

    // Get inter-sample differences
    gsl_matrix* dxMat_gsl = GetMatrixDifference(xMat_gsl);
        
    // Apply lag to the differenced data
    gsl_matrix* dxLaggedMatrix_gsl = GetMatrixLagged(dxMat_gsl, nlags);

    // Demean the lagged differenced data
    gsl_matrix* dxLaggedDemeanedMatrix_gsl = DeMean(dxLaggedMatrix_gsl);
    gsl_matrix_free(dxLaggedMatrix_gsl);

    // Pull out the difference data excluding the lagged samples, then demean them
    gsl_matrix* dxDemeanedMatrix_gsl = NULL;
    gsl_matrix* dxSubMat_gsl = GetSubMatrix(dxMat_gsl, nlags, -1, -1, -1);
    dxDemeanedMatrix_gsl = DeMean(dxSubMat_gsl);
    gsl_matrix_free(dxSubMat_gsl);
    gsl_matrix_free(dxMat_gsl);

    int nrx = dxLaggedDemeanedMatrix_gsl->size1;
    int ncx = dxLaggedDemeanedMatrix_gsl->size2;
    int nry = dxDemeanedMatrix_gsl->size1;
    int ncy = dxDemeanedMatrix_gsl->size2;

    // Mat divide the lagged diff data by the differenced data
    gsl_matrix* dmLagDemean = gsl_matrix_alloc(nrx, ncx);
    gsl_matrix_memcpy(dmLagDemean, dxLaggedDemeanedMatrix_gsl);
    gsl_matrix* dmDemean = gsl_matrix_alloc(nry, ncy);
    gsl_matrix_memcpy(dmDemean, dxDemeanedMatrix_gsl);
    gsl_matrix* tmp1 = MatrixDivide(dmLagDemean, dmDemean);
    gsl_matrix_free(dmLagDemean);
    gsl_matrix_free(dmDemean);

    // Mat multiply the lagged diff data by the differenced data to yield the fitted regression
    dmLagDemean = gsl_matrix_alloc(dxLaggedDemeanedMatrix_gsl->size1, dxLaggedDemeanedMatrix_gsl->size2);
    gsl_matrix_memcpy(dmLagDemean, dxLaggedDemeanedMatrix_gsl);
    gsl_matrix* fittedRegressionDX = MatrixMultiply(dmLagDemean, tmp1);
    gsl_matrix_free(dmLagDemean);
    gsl_matrix_free(tmp1);

    int nr = dxDemeanedMatrix_gsl->size1;
    int nc = dxDemeanedMatrix_gsl->size2;

    // Subtract the fitted regression from the residu
    gsl_matrix* ResidualsRegressionDX = gsl_matrix_alloc(nr, nc);
    gsl_matrix_memcpy(ResidualsRegressionDX, dxDemeanedMatrix_gsl);
    gsl_matrix_sub(ResidualsRegressionDX, fittedRegressionDX);

    nrx = xMat_gsl->size1;
    gsl_matrix* tmp4_gsl = NULL;
    tmp4_gsl = GetSubMatrix(xMat_gsl, 1, nrx - nlags - 1, -1, -1);
    gsl_matrix* xDemeanedMatrix_gsl = DeMean(tmp4_gsl);
    gsl_matrix_free(tmp4_gsl);

    dmLagDemean = gsl_matrix_alloc(dxLaggedDemeanedMatrix_gsl->size1, dxLaggedDemeanedMatrix_gsl->size2);
    gsl_matrix_memcpy(dmLagDemean, dxLaggedDemeanedMatrix_gsl);
    dmDemean = gsl_matrix_alloc(xDemeanedMatrix_gsl->size1, xDemeanedMatrix_gsl->size2);
    gsl_matrix_memcpy(dmDemean, xDemeanedMatrix_gsl);
    gsl_matrix* tmp6 = MatrixDivide(dmLagDemean, dmDemean);
    gsl_matrix_free(dmLagDemean);
    gsl_matrix_free(dmDemean);

    dmLagDemean = gsl_matrix_alloc(dxLaggedDemeanedMatrix_gsl->size1, dxLaggedDemeanedMatrix_gsl->size2);
    gsl_matrix_memcpy(dmLagDemean, dxLaggedDemeanedMatrix_gsl);
    gsl_matrix* fittedRegressionX = MatrixMultiply(dmLagDemean, tmp6);
    gsl_matrix_free(dmLagDemean);

    gsl_matrix_free(tmp6);

    gsl_matrix* ResidualsRegressionX = gsl_matrix_alloc(xDemeanedMatrix_gsl->size1, xDemeanedMatrix_gsl->size2);
    gsl_matrix_memcpy(ResidualsRegressionX, xDemeanedMatrix_gsl);
    gsl_matrix_sub(ResidualsRegressionX, fittedRegressionX);

    gsl_matrix* tposeResid = MatrixTransposeImpl(ResidualsRegressionX);
    gsl_matrix* tmp9 = MatrixMultiply(tposeResid,
        ResidualsRegressionX);
    gsl_matrix_free(tposeResid);
    MatrixDivideByElem(tmp9, (double)ResidualsRegressionX->size1);
    gsl_matrix* Skk = tmp9;

    tposeResid = MatrixTransposeImpl(ResidualsRegressionX);
    tmp9 = MatrixMultiply(tposeResid,
        ResidualsRegressionDX);
    gsl_matrix_free(tposeResid);
    MatrixDivideByElem(tmp9, (double)ResidualsRegressionX->size1);
    gsl_matrix* Sk0 = tmp9;

    tposeResid = MatrixTransposeImpl(ResidualsRegressionDX);
    tmp9 = MatrixMultiply(tposeResid,
        ResidualsRegressionDX);
    gsl_matrix_free(tposeResid);
    MatrixDivideByElem(tmp9, (double)ResidualsRegressionDX->size1);
    gsl_matrix* S00 = tmp9;

    tposeResid = MatrixTransposeImpl(Sk0);
    gsl_matrix* skkInverse = GetMatrixInverse(Skk);
    gsl_matrix* s00Inverse = GetMatrixInverse(S00);
    gsl_matrix* matMultTemp1 = MatrixMultiply(Sk0, s00Inverse);
    gsl_matrix* matMultTemp2 = MatrixMultiply(matMultTemp1, tposeResid);
    gsl_matrix* eigenInputMat = MatrixMultiply(skkInverse, matMultTemp2);
    gsl_matrix_free(tposeResid);
    gsl_matrix_free(skkInverse);
    gsl_matrix_free(s00Inverse);
    gsl_matrix_free(matMultTemp1);
    gsl_matrix_free(matMultTemp2);

    gsl_matrix_free(Skk);
    gsl_matrix_free(Sk0);
    gsl_matrix_free(S00);
    int n = eigenInputMat->size1;

    gsl_vector_complex* evalPtr = gsl_vector_complex_alloc(n);
    gsl_matrix_complex* ematPtr = gsl_matrix_complex_alloc(n, n);
    gsl_eigen_nonsymmv_workspace* worspacePtr = gsl_eigen_nonsymmv_alloc(n);
    gsl_eigen_nonsymmv(eigenInputMat, evalPtr, ematPtr, worspacePtr);
    gsl_matrix_free(eigenInputMat);
    gsl_eigen_nonsymmv_free(worspacePtr);

    gsl_eigen_nonsymmv_sort(evalPtr, ematPtr,
        GSL_EIGEN_SORT_ABS_DESC);

    eigenValuesVec = GSLComplexVecToAbsDoubleVector(evalPtr);
    eigenVecMatrix = GSLComplexMatToAbsDoubleMatrix(ematPtr);

    gsl_vector_complex_free(evalPtr);
    gsl_matrix_complex_free(ematPtr);

    int nSamples = ResidualsRegressionX->size1;
    int nVariables = ResidualsRegressionX->size2;

    int counter = 0;
    outStats.clear();
    for (int i = 0; i < (int)eigenValuesVec.size(); i++)
    {
        MaxEigenData eigData;

        double LR_maxeigenvalue = -nSamples * log(1 - eigenValuesVec[i]);
        eigData.TestStatistic = LR_maxeigenvalue;
        eigData.CriticalValue90 = maxEigenCriticalvalues[nVariables - counter - 1][0];
        eigData.CriticalValue95 = maxEigenCriticalvalues[nVariables - counter - 1][1];
        eigData.CriticalValue99 = maxEigenCriticalvalues[nVariables - counter - 1][2];
        counter++;
        outStats.push_back(eigData);
    }

    gsl_matrix_free(fittedRegressionDX);
    gsl_matrix_free(fittedRegressionX);
    gsl_matrix_free(ResidualsRegressionDX);
    gsl_matrix_free(ResidualsRegressionX);
    gsl_matrix_free(dxDemeanedMatrix_gsl);
    gsl_matrix_free(dxLaggedDemeanedMatrix_gsl);
    gsl_matrix_free(xDemeanedMatrix_gsl);
}

gsl_matrix* JohansenHelper::GetSubMatrix(gsl_matrix* xMat,
        int nBeginRow, int nEndRow,
        int nBeginCol, int nEndCol)
{
    int nr = xMat->size1;
    int nc = xMat->size2;

    if (nBeginRow == -1) nBeginRow = 0;
    if (nEndRow == -1) nEndRow = nr - 1;
    if (nBeginCol == -1) nBeginCol = 0;
    if (nEndCol == -1) nEndCol = nc - 1;

    int newRows = nEndRow - nBeginRow + 1;
    int newColumns = nEndCol - nBeginCol + 1;

    gsl_matrix_view subMatrix = gsl_matrix_submatrix(xMat, nBeginRow, nBeginCol, newRows, newColumns);
    gsl_matrix* retVal = gsl_matrix_alloc(subMatrix.matrix.size1, subMatrix.matrix.size2);
    gsl_matrix_memcpy(retVal, &subMatrix.matrix);

    return retVal;
}

gsl_matrix* JohansenHelper::GetMatrixDifference(gsl_matrix* xMat)
{
    int nr = xMat->size1;
    int nc = xMat->size2;

    gsl_matrix* diffMatrix = gsl_matrix_alloc(nr-1,nc);

    for (int i = 0; i < nc; i++)
    {
        for (int j = 0; j < nr - 1; j++)
        {
            double diff = gsl_matrix_get(xMat, j + 1, i) - gsl_matrix_get(xMat, j, i);
            gsl_matrix_set(diffMatrix, j, i, diff);
        }
    }

    return diffMatrix;
}

gsl_matrix* JohansenHelper::DeMean(gsl_matrix* xMat)
{
    int nr = xMat->size1;
    int nc = xMat->size2;
    gsl_matrix* retVal = gsl_matrix_alloc(nr,nc);

    for (int i = 0; i < nc; i++)
    {
        gsl_vector* currentCol = gsl_vector_alloc(nr);
        gsl_matrix_get_col(currentCol, xMat, i);

        double average = GetAverage(currentCol);
        for (int j = 0; j < nr; j++)
        {
            double demeanedValue = gsl_vector_get(currentCol, j) -average;
            gsl_vector_set(currentCol, j, demeanedValue);
        }
        gsl_matrix_set_col(retVal, i, currentCol);
        gsl_vector_free(currentCol);
    }

    return retVal;
}

gsl_matrix* JohansenHelper::GetMatrixLagged(gsl_matrix* xMat, int nlags)
{
    int nr = xMat->size1;
    int nc = xMat->size2;

    int nRows = nr - nlags;
    int nCols = nc * nlags;
    gsl_matrix* retVal = gsl_matrix_alloc(nRows, nCols);

    int counter = 0;
    int counter2 = 0;
    for (int i = 0; i < nCols; i++)
    {
        for (int j = 0; j < nRows; j++)
        {
            int laggedRow = j + nlags - counter - 1;
            double laggedValue = gsl_matrix_get(xMat, laggedRow, counter2);
            gsl_matrix_set(retVal, j, i, laggedValue);
        }
        counter++;
        if (counter >= nlags)
        {
            counter = 0;
            counter2++;
        }
    }

    return retVal;
}

double JohansenHelper::GetAverage(gsl_vector* x)
{
    double retVal = 0;
    int n = x->size;

    for (int i = 0; i < n; i++)
    {
        double element = gsl_vector_get(x, i);
        retVal += element;
    }

    retVal = retVal / n;

    return retVal;
}

void JohansenHelper::MatrixDivideByElem(gsl_matrix* xMat, double val)
{
    size_t nr = xMat->size1;
    size_t nc = xMat->size2;
    gsl_matrix* tmpMat = gsl_matrix_alloc(nr, nc);
    gsl_matrix_set_all(tmpMat, val);
    gsl_matrix_div_elements(xMat, tmpMat);
    gsl_matrix_free(tmpMat);
}

gsl_matrix* JohansenHelper::MatrixTransposeImpl(gsl_matrix* m)
{
    int nr = m->size1;
    int nc = m->size2;

    gsl_matrix* retmat = gsl_matrix_alloc(nc, nr);
    gsl_matrix_transpose_memcpy(retmat, m);

    return retmat;
}

gsl_matrix* JohansenHelper::MatrixMultiply(gsl_matrix* A, gsl_matrix* B)
{
    int nrA = A->size1;
    int ncB = B->size2;

    gsl_matrix* resMat = gsl_matrix_alloc(nrA, ncB);
    gsl_blas_dgemm(CblasNoTrans,
        CblasNoTrans,
        1.0, A, B,
        0.0, resMat);

    return resMat;
}

gsl_matrix* JohansenHelper::GetMatrixInverse(gsl_matrix* inMat)
{
    int size1 = inMat->size1;
    gsl_matrix* outMat = gsl_matrix_alloc(size1, size1);
    gsl_matrix* invert_me = gsl_matrix_alloc(size1, size1);
    gsl_permutation* perm = gsl_permutation_alloc(size1);
    int signum;
    gsl_matrix_memcpy(invert_me, inMat);
    gsl_linalg_LU_decomp(invert_me, perm, &signum);
    gsl_linalg_LU_invert(invert_me, perm, outMat);
    gsl_matrix_free(invert_me);
    gsl_permutation_free(perm);

    return outMat;
}

gsl_matrix* JohansenHelper::MatrixDivide(gsl_matrix* xMat, gsl_matrix* yMat)
{
    //Dim xTranspose As Variant
    gsl_matrix* xTranspose = MatrixTransposeImpl(xMat);
    gsl_matrix* tmp1 = MatrixMultiply(xTranspose, xMat);
    gsl_matrix* tmp2 = GetMatrixInverse(tmp1);
    gsl_matrix* tmp3 = MatrixMultiply(tmp2, xTranspose);
    gsl_matrix* tmp4 = MatrixMultiply(tmp3, yMat);

    gsl_matrix_free(tmp1);
    gsl_matrix_free(tmp2);
    gsl_matrix_free(tmp3);
    gsl_matrix_free(xTranspose);
    return tmp4;
}

gsl_matrix* JohansenHelper::DoubleMatrixToGSLMatrix(CommonTypes::DoubleMatrix doubleMat)
{
    int nc = doubleMat.size();
    int nr = doubleMat[0].size();
    gsl_matrix* x = gsl_matrix_alloc(nr, nc);

    for (int i = 0; i < nc; i++)
    {
        for (int j = 0; j < nr; j++)
        {
            gsl_matrix_set(x, j, i, doubleMat[i][j]);
        }
    }

    return x;
}

gsl_vector* JohansenHelper::DoubleVectorToGSLVector(CommonTypes::DoubleVector doubleVec)
{
    int nr = doubleVec.size();
    gsl_vector* x = gsl_vector_alloc(nr);

    for (int i = 0; i < nr; i++)
    {
        gsl_vector_set(x, i, doubleVec[i]);
    }

    return x;
}

CommonTypes::DoubleMatrix JohansenHelper::GSLMatrixToDoubleMatrix(gsl_matrix* gslMat)
{
    int nr = gslMat->size1;
    int nc = gslMat->size2;
    CommonTypes::DoubleMatrix retmat(nc);
    for (int i = 0; i < nc; i++)
    {
        for (int j = 0; j < nr; j++)
        {
            double element = gsl_matrix_get(gslMat, j, i);
            retmat[i].push_back(element);
        }
    }

    return retmat;
}

CommonTypes::DoubleVector JohansenHelper::GSLComplexVecToAbsDoubleVector(gsl_vector_complex* gslVec)
{
    int n = gslVec->size;
    CommonTypes::DoubleVector retvec(n);
    for (int i = 0; i < n; i++)
    {
        gsl_complex matComp = gsl_vector_complex_get(gslVec, i);
        double realval = matComp.dat[0];
        double imagval = matComp.dat[1];
        retvec[i] = realval;
    }

    return retvec;
}

CommonTypes::DoubleMatrix JohansenHelper::GSLComplexMatToAbsDoubleMatrix(gsl_matrix_complex* gslMat)
{
    int nr = gslMat->size1;
    int nc = gslMat->size2;
    CommonTypes::DoubleMatrix retmat(nr);
    for (int i = 0; i < nr; i++)
    {
        retmat[i].resize(nc);
        for (int j = 0; j < nc; j++)
        {
            gsl_complex matComp = gsl_matrix_complex_get(gslMat, i, j);
            double realval = matComp.dat[0];
            double imagval = matComp.dat[1];

            retmat[i][j] = realval;
        }
    }

    return retmat;
}

CommonTypes::DoubleVector JohansenHelper::GSLVectorToDoubleVector(gsl_vector* gslMat)
{
    CommonTypes::DoubleVector retVal;

    return retVal;
}
