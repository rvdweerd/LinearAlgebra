#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <tuple>

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

	// Basic matrix factories
	static Matrix Zeros(int m, int n);
	Matrix Eye(int n);
	
	// Operator overloads
	std::ostream& operator<<(std::ostream& stream, LinA::Matrix A);
	Matrix operator*(LinA::Matrix lhs, LinA::Matrix rhs);
	Matrix operator*(LinA::Matrix lhs, float rhs);
	Matrix operator-(LinA::Matrix lhs, LinA::Matrix rhs);
	Matrix operator+(LinA::Matrix lhs, LinA::Matrix rhs);

	// Basic matrix operations
	Matrix Transpose(const Matrix& A);
	void RowOp(Matrix& A, float multiplier, int subtract, int into);
	void RowExchange(Matrix& A, int row1, int row2);
	std::pair<Matrix, Matrix> RowReduce(Matrix& A, int col);
	std::pair<Matrix, Matrix> RowReduceUp(Matrix& A, int col);
	Matrix InverseEliminationMatrix(const Matrix& E);
	std::pair<Matrix, int> FixLowTriangular(const Matrix& E);
	std::pair<Matrix, Matrix> ProjectVec(Matrix b, Matrix a);
	Matrix GetColumn(const Matrix& A, int n);
	void ReplaceColumn(Matrix& A, const Matrix& sourceCol, int n);
	
	// Norms
	float VecNorm_L2(Matrix vec);

	// Matrix properties	
	Matrix SingularValues(Matrix A, float precision);
	float Det(const Matrix& A);
	float DetOfTriangular(const Matrix& E);
	
	// Decomposition functions
	std::pair<Matrix, Matrix> QR(Matrix A);
	Matrix Eig(Matrix A, float precision);
	std::tuple<Matrix, Matrix,Matrix> LU(Matrix A, bool printLU = true); // returns {L,U,E} (E=Elimination matrix)
	std::tuple<Matrix, Matrix> RREF(Matrix A); // returns {rref(A),E}
	Matrix Inverse(Matrix A);
}




