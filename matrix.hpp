/*
    matrix.h
    That defines a class in order to represent a matrix and its elementary operations

    Created by Vitor Manoel on 2021/04/14
    
    Copyright (c) 2021 Vitor-M. All rights reserved.

*/
#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

#include <ostream>
#include <vector>
#include <string>

using namespace std;
class Matrix
{
    private:
        string nametag; 
        unsigned int n_rows;
        unsigned int n_cols;
        vector<vector<double>> matrix;

    public:
        //Build Methods
        Matrix(string, unsigned, unsigned, vector<double>);
        Matrix(const string);
        Matrix(const Matrix &);
        Matrix();
        ~Matrix();

        //Matrix Operations
        Matrix operator+(Matrix );
        Matrix operator-(Matrix );
        Matrix operator*(Matrix );
        
        //Other Operations
        Matrix inverse();
        Matrix augmented(Matrix &);
        Matrix tranpose();
        Matrix identity();
        Matrix opposite();
        Matrix diagonal();

        //Scallar Operations
        Matrix operator+(double);
        Matrix operator-(double);
        Matrix operator*(double);
        Matrix operator/(double);

        //Comparison Operations
        bool operator!=(Matrix );
        bool operator==(Matrix );

        //Access Methods
        double& operator()(const unsigned&, const unsigned&);
        string toString() const;
        unsigned getRowsSize() const;
        unsigned getColsSize() const;
        string getNametag() const;
        void setNametag(string);

        vector<double> getRow(const unsigned&);
        vector<double> getCol(const unsigned&);
        void setRow(const unsigned&, vector<double>);
        void setCol(const unsigned&, vector<double>);

        //Error Matrix
        Matrix getError() const;
        
        //Firend class
        friend ostream& operator<< (ostream&, Matrix); 
};
#endif