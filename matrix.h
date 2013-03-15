//
//  matrix.h
//  SubmatrixQueries
//
//  Created by Raphael Bost on 02/01/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#ifndef __SubmatrixQueries__matrix__
#define __SubmatrixQueries__matrix__

#include <valarray>
#include <cassert>
#include <iostream>

namespace matrix {

    /*
     * class Matrix
     *
     * The matrix virtual class provides the blueprint for the matrices operation we need.
     * It also implements monotone and Monge property checkers.
     *
     */
    
    template <typename T>
    class Matrix {
    public:        
        virtual ~Matrix() {}
        
        virtual size_t rows() const = 0;
        virtual size_t cols() const = 0;
        
        virtual T& operator()(size_t i, size_t j) = 0;
        virtual T operator()(size_t i, size_t j) const = 0;
        
        bool isMonotone() const;
        bool isMonge() const;
        bool isInverseMonge() const;
        
        void print() const;
    };
    
    // Monotone & Monge propery
    
    template <typename T> bool Matrix<T>::isMonotone() const
    {
        for (size_t i = 0; i < rows(); i++) {
            for (size_t k = 0; k < cols(); k++) {
                for (size_t j = i+1; j < rows(); j++) {
                    for (size_t l = k+1; l < cols(); l++) {
                        if (operator()(i,k) > operator()(i,l) || operator()(j,k) <= operator()(j,l)) {
                            return false;
                        }
                    }
                }
            }
        }
        
        return true;
    }
    
    template <typename T> bool Matrix<T>::isMonge() const
    {
        for (size_t i = 0; i < rows()-1; i++) {
            for (size_t j = 0; j < cols()-1; j++) {
                // equivalent to the regular definition but uses substractions to avoid overflows
                if (operator()(i,j) - operator()(i+1,j)  > operator()(i,j+1) - operator()(i+1,j+1)) {
                    return false;
                }
            }
        }
        return true;
    }
    
    template <typename T> bool Matrix<T>::isInverseMonge() const
    {
        for (size_t i = 0; i < rows()-1; i++) {
            for (size_t j = 0; j < cols()-1; j++) {
                // equivalent to the regular definition but uses substractions to avoid overflows
                if (operator()(i,j) - operator()(i+1,j) < operator()(i,j+1)  - operator()(i+1,j+1) ) {
                    return false;
                }
            }
        }
        return true;
    }

    template <typename T> void Matrix<T>::print() const
    {
        for (size_t i = 0; i < rows(); i++) {
            for (size_t j = 0; j < cols(); j++) {
                std::cout << (*this)(i,j) << "  ;  ";
            }
            std::cout<<"\n";
        }        
    }
    
    /* class SimpleMatrix
     * 
     * This is a very simple matrix implemenation. Its purpose is to compare its performance
     * with the ComplexMatrix class based on valarray.
     */
    
    /*
     * Apparently, for the creation of a big matrix, this version is faster but
     * it uses a little bit more memory than the ComplexMatrix. Furthermore, for the
     * naive algorithm, the access performances seem to be worse than the ComplexMatrix.
     */
    
    template <typename T>
    class SimpleMatrix : public Matrix<T>
    {
    private:
        size_t _rows, _cols;
        T **_data;
                
        void initializeData()
        {
            _data = new T*[this->rows()];
            
            for (size_t i = 0; i < this->rows(); i++) {
                _data[i] = new T[this->cols()];
            }
        }
        
    public:
        SimpleMatrix(size_t rows, size_t cols) : _rows(rows), _cols(cols)
        {
            initializeData();
        }

        SimpleMatrix(size_t rows, size_t cols, T **data)
        {
            initializeData();
            for (size_t i = 0; i < this->rows(); i++) {
                for (size_t j = 0; j < this->rows(); j++) {
                    _data[i][j] = data[i][j];
                }
            }
        }
        SimpleMatrix(Matrix<T> *m) : _rows(m->rows()), _cols(m->cols())
        {
            initializeData();
            for (size_t i = 0; i < this->rows(); i++) {
                for (size_t j = 0; j < this->rows(); j++) {
                    _data[i][j] = (*m)(i,j);
                }
            }    
        }
        
        ~SimpleMatrix()
        {
            for (size_t i = 0; i < this->rows(); i++) {
                delete [] _data[i];
            }
            
            delete [] _data;
        }
        
        size_t rows() const { return _rows; }
        size_t cols() const { return _cols; }
        
        T& operator()(size_t i, size_t j)
        {
            return _data[i][j];
        }
        T operator()(size_t i, size_t j) const
        {
            return _data[i][j];
        }
    };
    
    /* class ComplexMatrix
     * This is a basic implementation of matrices based on the STL's valarray.
     * In the implementation of the SubmatrixQueries article, we are supposed to have a matrix implementation that 
     * can answer data queries on one entry in O(1) time. This is the case here.
     *
     */

    template <typename T>
    class ComplexMatrix : public Matrix<T> {
    public:
        ComplexMatrix(size_t rows, size_t cols);
        ComplexMatrix(size_t rows, size_t cols, std::valarray<T> data);
        ComplexMatrix(Matrix<T> *m);
        ComplexMatrix(ComplexMatrix<T> *m);
        
        size_t rows() const;
        size_t cols() const;
        
        std::valarray<T> row(size_t r) const;
        std::slice_array<T> row(size_t r);
        std::valarray<T> row(size_t r, size_t start, size_t end) const;
        std::slice_array<T> row(size_t r, size_t start, size_t end);
        
        std::valarray<T> col(size_t c) const;
        std::slice_array<T> col(size_t c);
        std::valarray<T> col(size_t c, size_t start, size_t end) const;
        std::slice_array<T> col(size_t c, size_t start, size_t end);
        
        std::valarray<T> submatrix(size_t minRow, size_t maxRow, size_t minCol, size_t maxCol) const;
        std::gslice_array<T> submatrix(size_t minRow, size_t maxRow, size_t minCol, size_t maxCol);
        
        T& operator()(size_t i, size_t j);
        T operator()(size_t i, size_t j) const;
        
        ComplexMatrix<T> transpose() const;
        
    private:
        size_t _rows;
        size_t _cols;
        std::valarray<T> _data;
    };

    template <typename T> ComplexMatrix<T>::ComplexMatrix(size_t rows, size_t cols) : _rows(rows),
    _cols(cols),
    _data(rows * cols)
    {
        assert(cols > 0);
    }

    template <typename T> ComplexMatrix<T>::ComplexMatrix(size_t rows, size_t cols, std::valarray<T> data) : _rows(rows),
    _cols(cols),
    _data(data)
    {
        assert(rows > 0);
        assert(cols > 0);
    }

    template <typename T> ComplexMatrix<T>::ComplexMatrix(ComplexMatrix<T> *m) : _rows(m->rows()), _cols(m->cols()), _data(m->_data)
    {
        
    }
    
    template <typename T> ComplexMatrix<T>::ComplexMatrix(Matrix<T> *m) : _rows(m->rows()), _cols(m->cols()), _data(m->rows() * m->cols())
    {
        for (size_t i = 0; i < this->rows(); i++) {
            for (size_t j = 0; j < this->cols(); j++) {
                (*this)(i,j) = (*m)(i,j);
            }
        }
    }
    // Accessors
    
    template <typename T> size_t ComplexMatrix<T>::rows() const{
        return _rows;
    }
    
    template <typename T> size_t ComplexMatrix<T>::cols() const{
        return _cols;
    }
    
    template<typename T>
    std::valarray<T> ComplexMatrix<T>::row(size_t r) const {
        assert(r >= 0 && r < rows());
        return _data[std::slice(r * cols(), cols(), 1)];
    }

    template<typename T>
    std::slice_array<T> ComplexMatrix<T>::row(size_t r) {
        assert(r >= 0 && r < rows());
        return _data[std::slice(r * cols(), cols(), 1)];
    }
    
    template<typename T>
    std::valarray<T> ComplexMatrix<T>::row(size_t r, size_t start, size_t end) const{
        assert(r >= 0 && r < rows());
        return _data[std::slice(r*cols()+start,end-start+1,1)];
    }

    template<typename T>
    std::slice_array<T> ComplexMatrix<T>::row(size_t r, size_t start, size_t end){
        assert(r >= 0 && r < rows());
        return _data[std::slice(r*cols()+start,end-start+1,1)];
    }
    
    template<typename T>
    std::valarray<T> ComplexMatrix<T>::col(size_t c) const {
        assert(c >= 0 && c < cols());
        return _data[std::slice(c, rows(), cols())];
    }
    
    template<typename T>
    std::slice_array<T> ComplexMatrix<T>::col(size_t c) {
        assert(c >= 0 && c < cols());
        return _data[std::slice(c, rows(), cols())];
    }
    
    template<typename T>
    std::valarray<T> ComplexMatrix<T>::col(size_t c, size_t start, size_t end) const{
        assert(c >= 0 && c < cols());
        return _data[std::slice(c+start*cols(),end-start+1,cols())];
    }
    
    template<typename T>
    std::slice_array<T> ComplexMatrix<T>::col(size_t c, size_t start, size_t end){
        assert(c >= 0 && c < cols());
        return _data[std::slice(c+start*cols(),end-start+1,cols())];
    }
    
    
    template<typename T>
    std::valarray<T> ComplexMatrix<T>::submatrix(size_t minRow, size_t maxRow, size_t minCol, size_t maxCol) const{
        assert(minRow <= maxRow);
        assert(minRow >= 0);
        assert(maxRow < rows());
        
        assert(minCol <= maxCol);
        assert(minCol >= 0);
        assert(maxCol < cols());
        
        size_t start = minRow*cols() + minCol;
        size_t lengths[]= {maxCol-minCol+1,maxRow-minRow+1};
        size_t strides[]= {1,cols()};
        
        
        return _data[std::gslice(start,std::valarray<size_t>(lengths,2),std::valarray<size_t>(strides,2))];
    }
    
    template<typename T>
    std::gslice_array<T> ComplexMatrix<T>::submatrix(size_t minRow, size_t maxRow, size_t minCol, size_t maxCol){
        assert(minRow <= maxRow);
        assert(minRow >= 0);
        assert(maxRow < rows());
        
        assert(minCol <= maxCol);
        assert(minCol >= 0);
        assert(maxCol < cols());

        size_t start = minRow*cols() + minCol;
        size_t lengths[]= {maxCol-minCol+1,maxRow-minRow+1};
        size_t strides[]= {1,cols()};
        
        
        return _data[std::gslice(start,std::valarray<size_t>(lengths,2),std::valarray<size_t>(strides,2))];
    }
    
    // Operators
    
    template <typename T> T& ComplexMatrix<T>::operator()(size_t i, size_t j)
    {
        assert(i >= 0 && j >= 0);
        assert(i<rows() && j < cols());
        return _data[i * _cols + j];
    }

    template <typename T> T ComplexMatrix<T>::operator()(size_t i, size_t j) const
    {
        assert(i >= 0 && j >= 0);
        assert(i<rows() && j < cols());
        return _data[i * _cols + j];
    }
    
    template <typename T> ComplexMatrix<T> ComplexMatrix<T>::transpose() const
    {
        ComplexMatrix<T> m = ComplexMatrix<T>(cols(), rows());
        
        for (size_t i = 0; i < rows(); i++) {
            std::valarray<T> r = row(i);
            m.col(i) = r;
        }
        return m;
    }
}
#endif /* defined(__SubmatrixQueries__matrix__) */
