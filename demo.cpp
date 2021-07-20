#include "matrix.hpp"
#include "matrix_exception.hpp"
#include <vector>
#include <iostream>

using namespace std;

int main()
{
    std::vector<double> elem_a = {2,5,2,3,6,1,4,5,8};
    std::vector<double> elem_b = {2,7,2,4,3,10,5,2,2};
    Matrix matrixA("A", 3, 3, elem_a);
    Matrix matrixB("B", 3, 3, elem_b);
    Matrix identity = matrixA.identity();

    std::cout << (matrixA * matrixB).toString() << std::endl;
    std::cout << (matrixA.augmented(identity)).toString() << std::endl;

    std::vector<double> elem_c = {2,-1,0,-1,2,-1,0,-1,2};
    Matrix matrixC("C", 3, 3, elem_c);
    std::cout << (matrixC.inverse()).toString() <<  std::endl;


    std::vector<double> elem_d = {2,1,5,3};
    Matrix matrixD("D", 2, 2, elem_d);
    std::cout << (matrixD.inverse()).toString() <<  std::endl;

    std::vector<double> elem_e = {1,2,3,0,1,4,0,0,1};
    Matrix matrixE("E", 3, 3, elem_e);
    std::cout << (matrixE.inverse()).toString() <<  std::endl;


    std::vector<double> aux = {0,0,0};
    matrixA.setRow(1,aux);
    std::cout << (matrixA).toString() << std::endl;
    matrixA.setCol(1,aux);
    std::cout << (matrixA).toString() << std::endl;
    aux = matrixA.getRow(2);
    for(unsigned i =0; i < 3; i++)
        std::cout << aux[i];
    aux = matrixA.getCol(2);
    for(unsigned i =0; i < 3; i++)
        std::cout << aux[i];


    return 0;
}