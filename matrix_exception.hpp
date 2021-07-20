/*
    matrix.h
    That defines a class in order to contains exceptions for matrix operations

    Created by Vitor Manoel on 2021/04/14
    
    Copyright (c) 2021 Vitor-M. All rights reserved.

*/
#ifndef _MATRIX_EXCEPTION_HPP_
#define _MATRIX_EXCEPTION_HPP_

#include <string>
#include <iostream>

class MatrixException
{
    public:
        MatrixException(std::string message){this->message = message;};
        void getMessage(){std::cout << "Matrix Error: " << message << std::endl;};

    private:
        std::string message;
};

#endif