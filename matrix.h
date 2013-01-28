//
//  matrix.h
//  KMNS
//
//  Created by Raphael Bost on 02/01/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#ifndef __KMNS__matrix__
#define __KMNS__matrix__

#include <valarray>
#include <cassert>

namespace matrix {

    /* class Matrix
     * This is a basic implementation of matrices based on the STL's valarray.
     * In the implementation of the KMNS article, we are supposed to have a matrix implementation that 
     * can answer data queries on one entry in O(1) time. This is the case here.
     *
     * I also added methods to test the Monge and monotonicity properties of the matrix.
     */

    template <typename T>
    class Matrix {
    public:
        Matrix(size_t rows, size_t cols);
        Matrix(size_t rows, size_t cols, std::valarray<T> data);
        
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
        
        Matrix<T> transpose() const;
        
        bool isMonotone() const;
        bool isMonge() const;
        bool isInverseMonge() const;
    private:
        size_t _rows;
        size_t _cols;
        std::valarray<T> _data;
    };

    template <typename T> Matrix<T>::Matrix(size_t rows, size_t cols) : _rows(rows),
    _cols(cols),
    _data(rows * cols)
    {
        assert(cols > 0);
    }

    template <typename T> Matrix<T>::Matrix(size_t rows, size_t cols, std::valarray<T> data) : _rows(rows),
    _cols(cols),
    _data(data)
    {
        assert(rows > 0);
        assert(cols > 0);
    }

    // Accessors
    
    template <typename T> size_t Matrix<T>::rows() const{
        return _rows;
    }
    
    template <typename T> size_t Matrix<T>::cols() const{
        return _cols;
    }
    
    template<typename T>
    std::valarray<T> Matrix<T>::row(size_t r) const {
        assert(r >= 0 && r < rows());
        return _data[std::slice(r * cols(), cols(), 1)];
    }

    template<typename T>
    std::slice_array<T> Matrix<T>::row(size_t r) {
        assert(r >= 0 && r < rows());
        return _data[std::slice(r * cols(), cols(), 1)];
    }
    
    template<typename T>
    std::valarray<T> Matrix<T>::row(size_t r, size_t start, size_t end) const{
        assert(r >= 0 && r < rows());
        return _data[std::slice(r*cols()+start,end-start+1,1)];
    }

    template<typename T>
    std::slice_array<T> Matrix<T>::row(size_t r, size_t start, size_t end){
        assert(r >= 0 && r < rows());
        return _data[std::slice(r*cols()+start,end-start+1,1)];
    }
    
    template<typename T>
    std::valarray<T> Matrix<T>::col(size_t c) const {
        assert(c >= 0 && c < cols());
        return _data[std::slice(c, rows(), cols())];
    }
    
    template<typename T>
    std::slice_array<T> Matrix<T>::col(size_t c) {
        assert(c >= 0 && c < cols());
        return _data[std::slice(c, rows(), cols())];
    }
    
    template<typename T>
    std::valarray<T> Matrix<T>::col(size_t c, size_t start, size_t end) const{
        assert(c >= 0 && c < cols());
        return _data[std::slice(c+start*cols(),end-start+1,cols())];
    }
    
    template<typename T>
    std::slice_array<T> Matrix<T>::col(size_t c, size_t start, size_t end){
        assert(c >= 0 && c < cols());
        return _data[std::slice(c+start*cols(),end-start+1,cols())];
    }
    
    
    template<typename T>
    std::valarray<T> Matrix<T>::submatrix(size_t minRow, size_t maxRow, size_t minCol, size_t maxCol) const{
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
    std::gslice_array<T> Matrix<T>::submatrix(size_t minRow, size_t maxRow, size_t minCol, size_t maxCol){
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
    
    template <typename T> T& Matrix<T>::operator()(size_t i, size_t j)
    {
        assert(i >= 0 && j >= 0);
        assert(i<rows() && j < cols());
        return _data[i * _cols + j];
    }

    template <typename T> T Matrix<T>::operator()(size_t i, size_t j) const
    {
        assert(i >= 0 && j >= 0);
        assert(i<rows() && j < cols());
        return _data[i * _cols + j];
    }
    
    template <typename T> Matrix<T> Matrix<T>::transpose() const
    {
        Matrix<T> m = Matrix<T>(cols(), rows());
        
        for (size_t i = 0; i < rows(); i++) {
            std::valarray<T> r = row(i);
            m.col(i) = r;
        }
        return m;
    }
    
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
                if (operator()(i,j) + operator()(i+1,j+1) > operator()(i,j+1) + operator()(i+1,j)) {
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
                if (operator()(i,j) + operator()(i+1,j+1) < operator()(i,j+1) + operator()(i+1,j)) {
                    return false;
                }
            }
        }
        return true;
    }
}
#endif /* defined(__KMNS__matrix__) */
