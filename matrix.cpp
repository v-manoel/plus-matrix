#include "matrix.hpp"
#include "matrix_exception.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

Matrix::Matrix(string _nametag, unsigned _n_rows, unsigned _n_cols, vector<double> cels_list)
{
    nametag = _nametag;
    n_rows = _n_rows;
    n_cols = _n_cols;
    matrix.resize(n_rows);

    if(cels_list.size() < (n_rows * n_cols))
    {
        cels_list.resize((n_cols * n_rows),0);
    }

    unsigned entries_rows = cels_list.size() / n_rows; 
    for(unsigned i=0; i < n_rows; i++)
    {
        matrix[i].resize(n_cols);
        for(unsigned j=0; j < n_cols; j++)
            matrix[i][j] = cels_list[j + (i * entries_rows)];
    }   
    
}

Matrix::Matrix(const string filename)
{

}

Matrix::Matrix(const Matrix &matrix_copy)
{
    nametag = matrix_copy.getNametag();
    n_rows = matrix_copy.getRowsSize();
    n_cols = matrix_copy.getColsSize();
    matrix = matrix_copy.matrix;
}

Matrix::Matrix()
{

}

Matrix::~Matrix()
{

}

//Matrix Operations
Matrix Matrix::operator+(Matrix matrix_b)
{
    try
    {
        if(n_cols != matrix_b.getColsSize() || n_rows != matrix_b.getRowsSize())
        {
            throw(MatrixException("Matrizes de tipos diferentes"));
        
        }else{
            Matrix sum((nametag + '+' + matrix_b.getNametag()), n_rows, n_cols, vector<double>());
            
            for(unsigned i = 0; i < n_rows; i++)
                for(unsigned j = 0; j < n_cols; j++)
                    sum(i,j) = matrix[i][j] + matrix_b(i,j);
            
            return sum;
        }  
    }
    catch(MatrixException& e)
    {
        e.getMessage();
        return getError();
    }
    

}

Matrix Matrix::operator-(Matrix matrix_b)
{
    try
    {
        if(n_cols != matrix_b.getColsSize() || n_rows != matrix_b.getRowsSize())
        {
            throw(MatrixException("Matrizes de tipos diferentes"));
        
        }else{
            Matrix dif((nametag + '-' + matrix_b.getNametag()), n_rows, n_cols, vector<double>());
            
            for(unsigned i = 0; i < n_rows; i++)
                for(unsigned j = 0; j < n_cols; j++)
                    dif(i,j) = matrix[i][j] - matrix_b(i,j);
            
            return dif;
        }  
    }
    catch(MatrixException& e)
    {
        e.getMessage();
        return getError();
    }
}

Matrix Matrix::operator*(Matrix matrix_b)
{
    //multiplication of the rows of the matrix by the columns of the matrix b
    try
    {
        if(n_cols !=  matrix_b.getRowsSize())
        {
            throw(MatrixException("O número de colunas da matriz " + nametag + " é diferente do numero das linhas de " + matrix_b.getNametag()));
        
        }else{
            Matrix mul("", n_rows, matrix_b.getColsSize(), vector<double>());
            Matrix aux("", n_rows, matrix_b.getColsSize(), vector<double>());

            unsigned n = 0;
            for(unsigned n=0; n < n_cols; n++)
            {
                for(unsigned i = 0; i < n_rows; i++)
                {
                    for(unsigned j = 0; j < matrix_b.getColsSize(); j++)
                            aux(i,j) = matrix[i][n] * matrix_b(n,j);
                }
                mul = mul + aux;
            }
            mul.setNametag(nametag + '*' + matrix_b.getNametag());
            return mul;
        }  
    }
    catch(MatrixException& e)
    {
        e.getMessage();
        return getError();
    }
}

//Other Operations
Matrix Matrix::inverse()
{
    //Uses Gauss-Jordan method to find the inverse matrix
    try
    {
        if(n_cols !=  n_rows)
        {
            throw(MatrixException("Matriz não quadrada"));
        
        }else{
            //Builds the augmented matrix using the identity matrix
            Matrix ide = identity();
            Matrix aug = augmented(ide);

            //to check pivot status
            bool pivo_is_zero = false;

            //swap lines to organize pivos
            for(unsigned row = 0; row < n_rows; row++)
            {
                /*
                Checks if there are any lines where the element below the pivot is different from zero. 
                Then replace the pivot line with this new line, as long as the previous condition is 
                satisfied in both directions. 
                */
                if(!aug(row,row))
                {
                    pivo_is_zero = true;
                    for(unsigned swp_row = 0; swp_row  < n_rows; swp_row ++)
                    {
                        if(aug(swp_row ,row) && aug(row,swp_row ))
                        {
                            vector<double> aux_row = aug.getRow(row);
                            aug.setRow(row, aug.getRow(swp_row));
                            aug.setRow(swp_row , aux_row);
                            pivo_is_zero = false;
                        }
                    }
                }

                if(pivo_is_zero)
                {
                    throw(MatrixException("The matrix has pivot with zero value"));
                }

                //Fill all cells bellow the pivto with zero
                for(unsigned next_row = row +1; next_row < n_rows; next_row++)
                {
                    double k = aug(next_row,row) / aug(row,row);
                    //Sum lines
                    for(unsigned col = 0; col < aug.getColsSize(); col++)
                        aug(next_row,col) -= k * aug(row,col);
                
                }


            }

            for(unsigned row = n_rows -1; row > 0; row--)
                //Fill all cells above the pivto with zero
                for(int previous_row = row -1; previous_row >= 0; previous_row--)
                {
                    double k = aug(previous_row,row) / aug(row,row);
                    //Sum lines
                    for(unsigned col = 0; col < aug.getColsSize(); col++)
                    {
                        aug(previous_row,col) -= k * aug(row,col);
                    }
                }

           
             Matrix inv((nametag + " - inverse"), n_rows, n_cols, vector<double>());

            //Divides each line by its pivot
            for(unsigned row = 0; row < n_rows; row++)
            {
                double pivot = aug(row,row);
                for(unsigned col = 0; col < aug.getColsSize(); col++)
                    aug(row,col) = aug(row,col) / pivot;
                    

                //buildes the inverse matrix
                for(unsigned col = 0; col < n_cols; col++)
                {
                    unsigned aux_col = aug.getColsSize()/2;
                    inv(row,col) = aug(row, (aux_col + col));
                }
                        
            }
        
            return inv;
        }  
    }
    catch(MatrixException& e)
    {
        e.getMessage();
        return getError();
    }
}

Matrix Matrix::augmented(Matrix &matrix_b)
{
    try
    {
        if(n_rows !=  matrix_b.getRowsSize())
        {
            throw(MatrixException("Matrizes nao possuem mesmo numero de linhas"));
        
        }else{
            Matrix aug((nametag + " - augmented"), n_rows, (n_cols + matrix_b.getColsSize()), vector<double>());
            
            for(unsigned i = 0; i < n_rows; i++)
            {
                for(unsigned j = 0; j < n_cols; j++)    //This matrix elements
                    aug(i,j) = matrix[i][j];
                for(unsigned j = 0; j < matrix_b.getColsSize(); j++)    //Aux matrix elements
                    aug(i,(n_cols + j)) = matrix_b(i,j);
            }
            
            return aug;
        }  
    }
    catch(MatrixException& e)
    {
        e.getMessage();
        return getError();
    }
}

Matrix Matrix::tranpose()
{
    Matrix tra((nametag + " - transpose"), n_rows, n_cols, vector<double>());

    for(unsigned i = 0; i < n_cols; i++)
        for(unsigned j = 0; j < n_rows; j++)
            tra(i,j) = matrix[j][i];

    return tra;
}

Matrix Matrix::identity()
{
    try
    {
        if(n_cols !=  n_rows)
        {
            throw(MatrixException("Matriz não quadrada"));
        
        }else{
            Matrix ide((nametag + " - identity"), n_rows, n_cols, vector<double>());
            
            for(unsigned i = 0; i < n_rows; i++)
                for(unsigned j = 0; j < n_cols; j++)
                    if(i == j)
                        ide(i,j) = 1;
            
            return ide;
        }  
    }
    catch(MatrixException& e)
    {
        e.getMessage();
        return getError();
    }
}

Matrix Matrix::opposite()
{
    Matrix opp((nametag + " - opposite"), n_rows, n_cols, vector<double>());
    
    for(unsigned i = 0; i < n_rows; i++)
        for(unsigned j = 0; j < n_cols; j++)
            opp(i,j) = -1* matrix[i][j];
    
    return opp;
}

Matrix Matrix::diagonal()
{
    unsigned matrix_size;
    if(n_cols < n_rows)
        matrix_size = n_cols;
    else
        matrix_size = n_rows;

    Matrix dia((nametag + " - diagonal"), matrix_size, matrix_size, vector<double>());

    for(unsigned i =0; i < n_rows; i++)
    {
        if(i < n_cols)
            dia(i,i) = matrix[i][i];
    }
    return dia;
}

//Scalar Operations
Matrix Matrix::operator+(double num)
{
    Matrix sum((nametag), n_rows, n_cols, vector<double>());

    for(unsigned i = 0; i < n_rows; i++)
        for(unsigned j = 0; j < n_cols; j++)
            sum(i,j) = matrix[i][j] + num;

    return sum;

}

Matrix Matrix::operator-(double num)
{
    Matrix dif((nametag), n_rows, n_cols, vector<double>());

    for(unsigned i = 0; i < n_rows; i++)
    for(unsigned j = 0; j < n_cols; j++)
        dif(i,j) = matrix[i][j] - num;

    return dif;
}

Matrix Matrix::operator*(double num)
{
    Matrix mul((nametag), n_rows, n_cols, vector<double>());

    for(unsigned i = 0; i < n_rows; i++)
    for(unsigned j = 0; j < n_cols; j++)
        mul(i,j) = matrix[i][j] * num;

    return mul;
}

Matrix Matrix::operator/(double num)
{
    try
    {
        if(!num)
        {
            throw(MatrixException("Divisão por zero"));
        
        }else{
            Matrix dif((nametag), n_rows, n_cols, vector<double>());
            
            for(unsigned i = 0; i < n_rows; i++)
                for(unsigned j = 0; j < n_cols; j++)
                    dif(i,j) = matrix[i][j] / num;
            
            return dif;
        }  
    }
    catch(MatrixException& e)
    {
        e.getMessage();
        return getError();
    }
}

bool Matrix::operator!=(Matrix matrix_b)
{
    if(n_rows != matrix_b.getRowsSize() || n_cols != matrix_b.getColsSize())
        return true;
    for(unsigned i = 0; i < n_rows; i++)
        for(unsigned j = 0; j < n_cols; j++)
            if(matrix[i][j] != matrix_b(i,j))
                return true;
    
    return false;
}

bool Matrix::operator==(Matrix matrix_b)
{
    if(n_rows != matrix_b.getRowsSize() || n_cols != matrix_b.getColsSize())
        return false;
    for(unsigned i = 0; i < n_rows; i++)
        for(unsigned j = 0; j < n_cols; j++)
            if(matrix[i][j] != matrix_b(i,j))
                return false;
    
    return true;
}

//Access Methods
double&  Matrix::operator()(const unsigned& row, const unsigned& col)
{
    return matrix[row][col];
}

string Matrix::toString() const
{
    stringstream out;
    out << nametag << endl;
    for(unsigned i = 0; i < n_rows; i++)
    {
        for(unsigned j = 0; j < n_cols; j++)
            out << matrix[i][j] << " ";
        out << endl;
    }

    return out.str();
}

void Matrix::setRow(const unsigned& row_id, vector<double> row)
{
    try
    {
        if(row_id >= n_rows)
        {
            throw(MatrixException("A matrix nao possui uma linha correspondente a informada"));
        }else if(row.size() != n_cols)
        {
            throw(MatrixException("Numero de elementos não coicidem"));
        }else{
            matrix[row_id] = row; 
        }  
    }
    catch(MatrixException& e)
    {
        e.getMessage();
    }
}

vector<double> Matrix::getRow(const unsigned& row_id)
{
    try
    {
        if(row_id >= n_rows)
        {
            throw(MatrixException("A matrix nao possui uma linha correspondente a informada"));

        }else{
            return matrix[row_id];
        }  
    }
    catch(MatrixException& e)
    {
        e.getMessage();
        vector<double> row;
        row.resize(0);
        return row;
    }
}

void Matrix::setCol(const unsigned& col_id, vector<double> col)
{
    try
    {
        if(col_id >= n_cols)
        {
            throw(MatrixException("A matrix nao possui uma coluna correspondente a informada"));
        }else if(col.size() != n_rows)
        {
            throw(MatrixException("Numero de elementos não coicidem"));
        }else{
            for(unsigned i = 0; i < n_rows; i++)
            {
                matrix[i][col_id] = col[i];
            } 
        }  
    }
    catch(MatrixException& e)
    {
        e.getMessage();
    }
}

vector<double> Matrix::getCol(const unsigned& col_id)
{
    vector<double> col;
    try
    {
        if(col_id >= n_cols)
        {
            throw(MatrixException("A matrix nao possui uma linha correspondente a informada"));

        }else{
            col.resize(n_rows);
            for(unsigned i = 0; i < n_rows; i++)
            {
                col[i] = matrix[i][col_id];
            } 
        }  
    }
    catch(MatrixException& e)
    {
        e.getMessage();
        col.resize(0);
    }
    
    return col;
}

unsigned Matrix::getRowsSize() const
{
    return n_rows;
}

unsigned Matrix::getColsSize() const
{
    return n_cols;
}

string Matrix::getNametag() const
{
    return nametag;
}

void Matrix::setNametag(string new_name)
{
    nametag = new_name;
}

 Matrix Matrix::getError() const
 {

     return Matrix("Error",0,0,vector<double>());
 }

 ostream& operator<< (ostream& os, Matrix matrix_)
 {
   os << matrix_.toString();
   return os;
 } 