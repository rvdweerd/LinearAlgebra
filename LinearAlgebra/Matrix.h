#pragma once
#include <string>
#include <vector>
#include <iostream>

namespace LinA
{
	class Matrix 
	{
	public:
		Matrix(std::string str);
		Matrix()
		{
			std::cout << "default ctor\n";
		}// = default;
		Matrix Transpose();
		void RowOp(int multiplier, int subtract, int into);
		void RowExchange(int row1, int row2);
	private:
		std::vector<std::vector<float>> A;
		size_t m = 0; // number of rows
		size_t n = 0; // number of columns

		
		friend std::ostream& operator<<(std::ostream& stream, LinA::Matrix A);
		friend Matrix operator*(LinA::Matrix lhs, LinA::Matrix rhs);
		friend void RowOp(Matrix& A, int multiplier, int subtract, int into);
		friend void RowExchange(Matrix& A, int row1, int row2);
	};

	void RowOp(Matrix& A, int multiplier, int subtract, int into);
	void RowExchange(Matrix& A, int row1, int row2);
}




