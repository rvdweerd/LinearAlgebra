#include "Matrix.h"
#include <strstream>
#include <iostream>
#include <algorithm>
#include <ostream>
#include <numeric>
#include <iterator>

LinA::Matrix LinA::Zeros(int m, int n)
{
	return LinA::Matrix(m,n,0);
}

LinA::Matrix LinA::Eye(int n)
{
	LinA::Matrix X = LinA::Matrix(n, n, 0);
	for (size_t i = 0; i < n; i++)
	{
		X.A[i][i] = 1;
	}
	return X;
}

void LinA::RowOp(LinA::Matrix& A,float multiplier, int subtract, int into)
{
	std::transform(A.A[into - 1].begin(), A.A[into - 1].end(),							// transform this row (elimination)
		A.A[subtract - 1].begin(),										// by subtracting elements from this row
		A.A[into - 1].begin(), [&multiplier](float& val1, float& val2)		// using the multiplier
		{
			return val1 - val2 * multiplier;
		});

}

void LinA::RowExchange(LinA::Matrix& A, int row1, int row2)
{
	std::swap(A.A[row1 - 1], A.A[row2 - 1]);
}

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
			n = x;
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
}

LinA::Matrix::Matrix(int m, int n, float val)
	:
	m(m),
	n(n)
{
	std::vector<float> row(n, val);
	for (size_t y = 0; y < m; y++)
	{
		A.push_back(row);
	}
}

LinA::Matrix LinA::Matrix::Transpose()
{
	Matrix A_T;
	A_T.m = n;
	A_T.n = m;
	std::vector<float> v;
	for (size_t col = 0; col < n; col++)
	{
		for (size_t row = 0; row < m; row++)
		{
			v.push_back(A[row][col]);
		}
		A_T.A.push_back(v);
		v.clear();
	}
	return A_T;
}

void LinA::Matrix::RowOp(float multiplier, int subtract, int into)
{
	std::transform( A[into - 1].begin(), A[into - 1].end(),							// transform this row (elimination)
					A[subtract - 1].begin(),										// by subtracting elements from this row
					A[into-1].begin(), [&multiplier](float& val1, float& val2)		// using the multiplier
					{
						return val1 - val2 * multiplier; 
					});
	std::cout << '\n';
}

void LinA::Matrix::RowExchange(int row1, int row2)
{
	std::swap(A[row1 - 1], A[row2 - 1]);
}

LinA::Matrix LinA::Matrix::RowReduce(int col)
{
	LinA::Matrix E = LinA::Eye(m);

	if (A[col - 1][col - 1] > 1e-4) // pivot value acceptable
	{
		for (size_t i = col; i < m; i++)
		{
			//if (A[i - 1][col - 1] > 1e-4)
			{
				float m = A[i][col - 1] / A[col - 1][col - 1];
				std::cout << "m=" << m << ", ";
				RowOp(m, col , (int)i+1);
				E.A[i][(int)col - 1] = -m;
			}
		}
	}
	return E;
}


std::ostream& LinA::operator<<(std::ostream& stream, LinA::Matrix A)
{
	stream << '\n';
	for (size_t row = 0; row < A.m; row++)
	{
		for (size_t col = 0; col < A.n; col++)
		{
			stream << " " << A.A[row][col] << ", ";
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