#pragma once
#include <string>
#include <vector>
#include <iostream>

namespace LinA
{
	static constexpr float accuracy = 0.0001f;
	struct Matrix 
	{
		Matrix(std::string str);
		Matrix(int m, int n, float val);
		Matrix() = default;
		~Matrix() = default;
	
		std::vector<std::vector<float>> A;
		size_t m = 0; // number of rows
		size_t n = 0; // number of columns
	};

	static Matrix Zeros(int m, int n);
	static Matrix Eye(int n);
	
	std::ostream& operator<<(std::ostream& stream, LinA::Matrix A);
	Matrix operator*(LinA::Matrix lhs, LinA::Matrix rhs);
	Matrix operator*(LinA::Matrix lhs, float rhs);
	Matrix operator-(LinA::Matrix lhs, LinA::Matrix rhs);
	Matrix operator+(LinA::Matrix lhs, LinA::Matrix rhs);

	
	Matrix Transpose(const Matrix& A);
	void RowOp(Matrix& A, float multiplier, int subtract, int into);
	void RowExchange(Matrix& A, int row1, int row2);
	std::pair<Matrix, Matrix> RowReduce(Matrix& A, int col);
	Matrix InverseEliminationMatrix(const Matrix& E);
	std::pair<Matrix, int> FixLowTriangular(const Matrix& E);
	float DetOfTriangular(const Matrix& E);
}




