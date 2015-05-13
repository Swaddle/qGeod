/*
    * Helper class so we can send and recieve armadillo matrices over MPI
    *
*/

#include <mpi.h>

#include <cstdlib>
#include <armadillo>
#include <vector>
#include <iostream>

using arma::cx_mat;
using std::vector;

class ArmadilloMPI
{
public:
    ArmadilloMPI(int nRows, int nCols)
    {
        this->nRows = nRows;
        this->nCols = nCols;
        double *data = new double[nCols * nRows ];
        realArray = new double*[( nRows )];
        for(int i = 0; i < nRows; ++i)
        {
            realArray[i] = &data[i*nCols];
        }
        double *datai = new double[(nCols * nRows )];
        imArray = new double*[( nRows )];
        for(int i = 0; i < nRows; ++i)
        {
            imArray[i] = &datai[i*nCols];
        }

    }

    ~ArmadilloMPI()
    {
        delete[] realArray[0];
        delete[] realArray;
        delete[] imArray[0];
        delete[] imArray;
    }

    cx_mat matConstructRecv( int src, int tag)
    {
        mat realMat(nRows, nCols);
        mat imMat(nRows, nCols);

        MPI_Recv(&(realArray[0][0]), nRows * nCols, MPI_DOUBLE, src, tag, MPI_COMM_WORLD,0);
        MPI_Recv(&(imArray[0][0]),  nRows * nCols, MPI_DOUBLE, src, tag + 1, MPI_COMM_WORLD,0);

        for(int  i = 0; i < nRows; ++i )
        {
            for(int j = 0; j < nCols; ++j)
            {
                realMat(i,j) = realArray[i][j];
                imMat(i,j) = imArray[i][j];
            }
        }

        cx_mat A = cx_mat(realMat, imMat);

        return A;
    }

    void matDestroySend(cx_mat &A, int dest, int tag)
    {
        for(int  i = 0; i < nRows; ++i )
        {
            for(int j = 0; j < nCols; ++j)
            {
                realArray[i][j] = real((A(i,j)));
                imArray[i][j] = imag((A(i,j)));
            }
        }
        MPI_Send(&(realArray[0][0]), nRows * nCols, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&(imArray[0][0]), nRows * nCols, MPI_DOUBLE, dest, tag + 1, MPI_COMM_WORLD);
    }

private:

    int nCols;
    int nRows;
    double **realArray;
    double **imArray;

};
