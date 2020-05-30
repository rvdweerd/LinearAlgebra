#include "Matrix.h"
#include <strstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <numeric>
#include <iterator>
#include <cassert>

LinA::Matrix::Matrix(std::string str)
{
	std::istrstream stream_in(str.c_str());
	std::string s;
	size_t x = 0;
	size_t y = 0;
	std::vector<float> row;
	while (stream_in >> s)
	{
		if (s == ";")
		{
			//n = x;
			x = 0;
			A.push_back(row);
			row.clear();
			y++;
		}
		else if (s == ",")
		{ }
		else
		{
			std::istrstream sub(s.c_str());
			float f;
			sub >> f;
			row.push_back(f);
			x++;
		}
	}
	A.push_back(row);
	m = y+1;
	n = A[0].size();
}

LinA::Matrix::Matrix(int m, int n, float val)
	:
	m(m),
	n(n)
{
	std::vector<float> row(n, val);
	for (size_t y = 0; y < (size_t)m; y++)
	{
		A.push_back(row);
	}
}

LinA::Matrix LinA::Zeros(int m, int n)
{
	return LinA::Matrix(m,n,0);
}

LinA::Matrix LinA::Eye(int n)
{
	LinA::Matrix X = LinA::Matrix(n, n, 0);
	for (size_t i = 0; i < (size_t)n; i++)
	{
		X.A[i][i] = 1;
	}
	return X;
}

std::ostream& LinA::operator<<(std::ostream& stream, LinA::Matrix A)
{
	stream.precision(3);
	stream << '\n';
	for (size_t row = 0; row < A.m; row++)
	{
		for (size_t col = 0; col < A.n; col++)
		{
			stream << " " << std::setw(10) << A.A[row][col] << ", ";
		}
		stream << "\n";
	}
	stream << '\n';
	return stream;
}

LinA::Matrix LinA::operator*(LinA::Matrix lhs, LinA::Matrix rhs)
{
	if (lhs.n == rhs.m)
	{
		LinA::Matrix Mret;
		Mret.m = lhs.m;
		Mret.n = rhs.n;
		for (size_t row = 0; row < lhs.m; row++)
		{
			std::vector<float> v;
			for (size_t col = 0; col < rhs.n; col++)
			{
				float prodsum = 0;
				for (size_t i = 0; i < lhs.n; i++)
				{
					prodsum += lhs.A[row][i] * rhs.A[i][col];
				}
				v.push_back(prodsum);
			}
			Mret.A.push_back(v);
			v.clear();
		}
		return Mret;
	}
	return LinA::Matrix();
}

LinA::Matrix LinA::operator*(LinA::Matrix lhs, float rhs)
{
	for (auto& vec : lhs.A)
	{
		for (auto& v : vec)
		{
			v *= rhs;
		}
	}
	return lhs;
}

LinA::Matrix LinA::operator-(LinA::Matrix lhs, LinA::Matrix rhs)
{
	assert(lhs.A.size() == rhs.A.size());
	assert(lhs.A[0].size() == rhs.A[0].size());
	for (size_t row = 0; row<lhs.A.size();row++)
	{
		for (size_t col = 0; col < lhs.A[0].size(); col++)
		{
			lhs.A[row][col] -= rhs.A[row][col];
		}
	}
	return lhs;
	return Matrix();
}

LinA::Matrix LinA::operator+(LinA::Matrix lhs, LinA::Matrix rhs)
{
	assert(lhs.A.size() == rhs.A.size());
	assert(lhs.A[0].size() == rhs.A[0].size());
	for (size_t row = 0; row < lhs.A.size(); row++)
	{
		for (size_t col = 0; col < lhs.A[0].size(); col++)
		{
			lhs.A[row][col] += rhs.A[row][col];
		}
	}
	return lhs;
	return Matrix();
}

LinA::Matrix LinA::Transpose(const Matrix& A)
{
	Matrix A_T;
	A_T.m = A.n;
	A_T.n = A.m;
	std::vector<float> v;
	for (size_t col = 0; col < A.n; col++)
	{
		for (size_t row = 0; row < A.m; row++)
		{
			v.push_back(A.A[row][col]);
		}
		A_T.A.push_back(v);
		v.clear();
	}
	return A_T;
}

void LinA::RowOp(LinA::Matrix& A,float multiplier, int subtract, int into)
{
	std::transform(A.A[into - 1].begin(), A.A[into - 1].end(),			// transform this row (elimination)
		A.A[subtract - 1].begin(),										// by subtracting elements from this row
		A.A[into - 1].begin(), [&multiplier](float& val1, float& val2)	// using the multiplier
		{
			return val1 - val2 * multiplier;
		});

}

void LinA::RowExchange(LinA::Matrix& A, int row1, int row2)
{
	std::swap(A.A[row1 - 1], A.A[row2 - 1]);
}

std::pair<LinA::Matrix, LinA::Matrix> LinA::RowReduce(Matrix& A, int col)
{
	LinA::Matrix E_rowexch = LinA::Eye(A.m);
	LinA::Matrix E_rowexch_inv = LinA::Eye(A.m);
	if (col == A.n) return { E_rowexch,E_rowexch_inv };
	int rowrunner = col+1;
	while ( abs(A.A[col - 1][col - 1]) < LinA::accuracy && rowrunner <= (int)A.n)
	{
		LinA::RowExchange(A,col, rowrunner);
		LinA::RowExchange(E_rowexch, col, rowrunner);
		LinA::RowExchange(E_rowexch_inv, col, rowrunner);
		rowrunner++;
	}

	LinA::Matrix E_rowmult = LinA::Eye(A.m);
	LinA::Matrix E_rowmult_inv = LinA::Eye(A.m);
	if (abs(A.A[col - 1][col - 1]) > LinA::accuracy )//&& abs(A.A[col][col-1]) > LinA::accuracy) // pivot value acceptable
	{
		for (size_t i = col; i < A.m; i++)
		{
			if ( abs(A.A[i][col - 1]) > LinA::accuracy)
			{
				float m = A.A[i][col - 1] / A.A[col - 1][col - 1];
				std::cout << "m=" << m << ", ";
				RowOp(A,m, col, (int)i + 1);
				E_rowmult.A[i][(int)col - 1] = -m;
			}
		}
		E_rowmult_inv = InverseEliminationMatrix(E_rowmult);
	}
	return { E_rowmult*E_rowexch,E_rowexch_inv*E_rowmult_inv };
}

LinA::Matrix LinA::InverseEliminationMatrix(const Matrix& E)
{
	LinA::Matrix E_inv(E);
	for (size_t i = 1; i <= E.m; i++)
	{
		for (size_t j = 1; j < i; j++)
		{
			float entry = E_inv.A[i - 1][j - 1];
			if (std::abs(entry) > 1e-4)
			{
				E_inv.A[i - 1][j - 1] = -E_inv.A[i - 1][j - 1];
			}
		}
	}
	return E_inv;
}

std::pair<LinA::Matrix, int> LinA::FixLowTriangular(const Matrix& E)
{
	int n = 0;
	LinA::Matrix E_cpy = E;
	for (size_t rows = 1; rows < E.m; rows++)
	{
		if (abs(E_cpy.A[rows - 1][E.n - 1]) > 0.004) // check if last column entries are 0 except for A[m][n]
		{
			std::swap(E_cpy.A[rows - 1], E_cpy.A[rows]);
			n++;
		}
	}
	return std::pair<Matrix, int>(E_cpy,n);
}

float LinA::DetOfTriangular(const Matrix& E)
{
	if (E.n == E.m)
	{
		//check if triangular...
		float det = 1.f;
		for (size_t i = 1; i <= E.n; i++)
		{
			det *= E.A[i - 1][i - 1];
		}
		return det;
	}
	assert(false);
	return -9999;
}

