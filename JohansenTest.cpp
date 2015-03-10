// JohansenTest.cpp : Defines the entry point for the console application.
//

#include "JohansenHelper.h"
#include <time.h>
#include <iostream>
#include <fstream>
using namespace std;

void ClearDoubleMatrix(DoubleMatrix& xMat)
{
    for (int i = 0; i < (int)xMat.size();i++)
    {
        xMat[i].erase(xMat[i].begin(),xMat[i].end());
    }
    xMat.erase(xMat.begin(), xMat.end());
}

void FileData(DoubleMatrix& xMat)
{
    ClearDoubleMatrix(xMat);
    xMat.resize(2);
    ifstream inFile;
    inFile.open("gld.csv");
    if (inFile.fail())
    {
        printf("input file not found\n");
        exit(0);
    }
    for (int i = 0; i < 5000; i++)
    {
        string line;
        if (getline(inFile, line))
        {
            double value = atof(line.c_str());
            xMat[0].push_back(value);
        }
        else
        {
            break;
        }
    }
    inFile.close();

    inFile.open("slv.csv");
    if (inFile.fail())
    {
        printf("input file not found\n");
        exit(0);
    }
    for (;;)
    {
        string line;
        if (getline(inFile, line))
        {
            double value = atof(line.c_str());
            xMat[1].push_back(value);
        }
        else
        {
            break;
        }
    }
    inFile.close();
}

void StationarySeries(int seriesCount, int sampleCount, DoubleMatrix& xMat)
{
    ClearDoubleMatrix(xMat);
    xMat.resize(seriesCount);

    srand((unsigned int)time(NULL));
    for (int i = 0; i < seriesCount; i++)
    {
        for (int j = 0; j < sampleCount; j++)
        {
            double randomValue = rand();
            xMat[i].push_back(randomValue);
        }
    }
}

/*
 * Instead of generating a number from 0 to RAND_MAX, shift the range down so that half are negative.
 * Range -RAND_MAX/2 to + RAND_MAX/2
*/
int MyRandom()
{
    int retVal = rand();

    retVal = retVal - ((unsigned int)RAND_MAX + 1) / 2;

    return retVal;
}

typedef enum
{
    SCENARIO_NO_COINTEGRATION,
    SCENARIO_1_COINTEGRATION,
    SCENARIO_2_COINTEGRATIONS
} RandomWalkScenarioType;

void RandomWalk(int seriesCount, int sampleCount, DoubleMatrix& xMat, int scenario)
{
    ClearDoubleMatrix(xMat);
    DoubleMatrix randomValues(seriesCount);
    xMat.resize(seriesCount);

    // Fill the random values and initialize output to 0
    srand((unsigned int)time(NULL));
    for (int i = 0; i < seriesCount; i++)
    {
        for (int j = 0; j < sampleCount; j++)
        {
            randomValues[i].push_back(MyRandom());
            xMat[i].push_back(0);
        }
    }

    for (int i = 0; i < seriesCount; i++)
    {
        xMat[i][0] = 0; // Every series starts walking from 0
        for (int j = 1; j < sampleCount; j++)
        {
            xMat[i][j] = xMat[i][j - 1] + randomValues[i][j - 1];
        }
    }

    int cointegratedSeriesStart = 0;
    int seriesLength = 0;
    switch (scenario)
    {
    case SCENARIO_NO_COINTEGRATION:
        cointegratedSeriesStart = seriesCount;
        break;

    case SCENARIO_1_COINTEGRATION:
        cointegratedSeriesStart = 3;
        seriesLength = 1;
        break;

    case SCENARIO_2_COINTEGRATIONS:
        cointegratedSeriesStart = 3;
        seriesLength = 2;
        break;

    default:
        cointegratedSeriesStart = seriesCount;
        seriesLength = 0;
        break;
    }

    // For scenario 1, there is 1 coint series at index 3
    // For scenario 2, there are 2 coint series at indexes 2 and 3
    for (int i = cointegratedSeriesStart; i < (cointegratedSeriesStart + seriesLength); i++)
    {
        for (int j = 1; j < sampleCount; j++)
        {
            if (i == 3)
            {
                xMat[i][j] = 0.5 * xMat[0][j] + randomValues[i][j - 1];
            }
            else
            {
                xMat[i][j] = 1 * xMat[0][j] + 2 * xMat[1][j] + 3 * xMat[2][j] + randomValues[i][j - 1];
            }
        }
    }
}

int main(int argc, char* argv[])
{
    DoubleMatrix xMat;
    int nlags = 1;
    int cointegrationCount = 0;
#if 1
    for (int i = 0; i < 4; i++)
    {
        cout << endl;
        switch (i)
        {
        default:
        case 0:
            cout << "Stationary Series" << endl;
            StationarySeries(5, 1000, xMat);
            break;

        case 1:
            cout << "Random Walk, No Cointegration" << endl;
            RandomWalk(5, 1000, xMat, SCENARIO_NO_COINTEGRATION);
            break;

        case 2:
            cout << "Random Walk, 1 Cointegration" << endl;
            RandomWalk(5, 1000, xMat, SCENARIO_1_COINTEGRATION);
            break;

        case 3:
            cout << "Random Walk, 2 Cointegrations" << endl;
            RandomWalk(5, 1000, xMat, SCENARIO_2_COINTEGRATIONS);
            break;
        }

        vector<MaxEigenData> outStats;
        DoubleVector eigenValuesVec;
        DoubleMatrix eigenVecMatrix;
        JohansenHelper johansenHelper(xMat);
        johansenHelper.DoMaxEigenValueTest(nlags);
        int cointCount = johansenHelper.CointegrationCount();
        outStats = johansenHelper.GetOutStats();
        for (int i = 0; i < (int)outStats.size(); i++)
        {
            printf("%d Test Stat: %f 90: %f 95: %f 99: %f\n",
                i, outStats[i].TestStatistic, outStats[i].CriticalValue90, outStats[i].CriticalValue95, outStats[i].CriticalValue99);
        }

        cout << "coint count: " << cointCount << endl;
    }
#else
    FileData(xMat);
    JohansenHelper johansenHelper(xMat);

    for (int i = 1; i < 365; i++)
    {
        johansenHelper.DoMaxEigenValueTest(i);
        int cointCount = johansenHelper.CointegrationCount();
        if (cointCount > 0)
        {
            printf("lags: %d, count: %d\n", i, cointCount);
        }
        cout << "lags: " << i << " coint count: " << cointCount << endl;
    }
#endif
    return 0;
}
