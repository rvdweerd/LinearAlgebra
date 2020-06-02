#include "Matrix.h"
#include <strstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <numeric>
#include <iterator>
#include <cassert>
#include <tuple>

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
	stream.precision(2);
	stream << '\n';
	for (size_t row = 0; row < A.m; row++)
	{
		for (size_t col = 0; col < A.n; col++)
		{
			if (abs(A.A[row][col]) > LinA::accuracy)
			{
				stream << " " << std::setw(10) << A.A[row][col] << ", ";
			}
			else
			{
				stream << " " << std::setw(10) << 0.0f << ", ";
			}
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
				//std::cout << "m=" << m << ", ";
				RowOp(A,m, col, (int)i + 1);
				E_rowmult.A[i][(int)col - 1] = -m;
			}
		}
		E_rowmult_inv = InverseEliminationMatrix(E_rowmult);
	}
	return { E_rowmult*E_rowexch , E_rowexch_inv*E_rowmult_inv };
}
std::pair<LinA::Matrix, LinA::Matrix> LinA::RowReduceUp(Matrix& A, int col)
{
	assert(A.m == A.n);
	LinA::Matrix E_rowexch = LinA::Eye(A.m);
	LinA::Matrix E_rowexch_inv = LinA::Eye(A.m);
	//if (col == A.n) return { E_rowexch,E_rowexch_inv };
	int rowrunner = col - 2;
	while (abs(A.A[col - 1][col - 1]) < LinA::accuracy && rowrunner >= 0 )
	{
		assert(false);
		//LinA::RowExchange(A, col, rowrunner);
		//LinA::RowExchange(E_rowexch, col, rowrunner);
		//LinA::RowExchange(E_rowexch_inv, col, rowrunner);
		//rowrunner++;
	}

	LinA::Matrix E_rowmult = LinA::Eye(A.m);
	LinA::Matrix E_rowmult_inv = LinA::Eye(A.m);
	if (abs(A.A[col - 1][col - 1]) > LinA::accuracy)//&& abs(A.A[col][col-1]) > LinA::accuracy) // pivot value acceptable
	{
		for (int i = col-2; i >=0; i--)
		{
			if (abs(A.A[i][col - 1]) > LinA::accuracy)
			{
				float m = A.A[i][col - 1] / A.A[col - 1][col - 1];
				//std::cout << "m=" << m << ", ";
				RowOp(A, m, col, (int)i+1 );
				E_rowmult.A[i][(int)col - 1] = -m;
			}
		}
		E_rowmult_inv = InverseEliminationMatrix(E_rowmult);
	}
	//std::cout << E_rowmult * E_rowexch;
	//std::cout << E_rowexch_inv * E_rowmult_inv;
	return { E_rowmult * E_rowexch,E_rowexch_inv * E_rowmult_inv };
}
LinA::Matrix LinA::InverseEliminationMatrix(const Matrix& E)
{
	assert(E.m == E.n);
	LinA::Matrix E_inv(E);
	for (size_t i = 1; i <= E.m; i++)
	{
		for (size_t j = 1; j <= E.n; j++)
		{
			//float entry = E_inv.A[i - 1][j - 1];
			if (std::abs(E_inv.A[i - 1][j - 1]) > 1e-4 && i!=j)
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
LinA::Matrix LinA::GetColumn(const LinA::Matrix& A, int n)
{
	assert(A.A.size() > 0);
	assert(A.A[0].size() > 0);
	assert(A.n >= n);
	LinA::Matrix col(A.m, 1, 0);
	for (size_t row = 0; row < A.m; row++)
	{
		col.A[row][0] = A.A[row][n - 1];
	}
	return col;
}
void LinA::ReplaceColumn(LinA::Matrix& A, const LinA::Matrix& sourceCol, int n)
{
	assert(A.m == sourceCol.m);
	assert(n <= A.n);
	for (size_t col = 0; col < A.m; col++)
	{
		A.A[col][n - 1] = sourceCol.A[col][0];
	}
}


std::pair<LinA::Matrix, LinA::Matrix> LinA::ProjectVec(LinA::Matrix b, LinA::Matrix a) // project vector b onto vector a, returns { p , e }  
{
	assert(a.n == 1 && b.n == 1);
	assert(a.m == b.m);
	LinA::Matrix aT = Transpose(a);
	float x_hat = (aT * b).A[0][0] / (aT * a).A[0][0];
	LinA::Matrix p = a * x_hat;
	LinA::Matrix e = b - p;
	return { p,e };
}

float LinA::VecNorm_L2(LinA::Matrix vec)
{
	assert(vec.m == 1 || vec.n == 1);
	if (vec.n == 1)
	{
		return sqrt((Transpose(vec) * vec).A[0][0]);
	}
	if (vec.m == 1)
	{
		return sqrt((vec * Transpose(vec)).A[0][0]);
	}
	return -1;
}

std::pair<LinA::Matrix, LinA::Matrix> LinA::QR(LinA::Matrix A)
{
	assert(A.m == A.n);
	assert(A.m > 1);
	LinA::Matrix Q(A.m, A.n, 0);
	LinA::Matrix col0 = GetColumn(A, 1);
	col0 = col0 * (1 / VecNorm_L2(col0));
	ReplaceColumn(Q, col0, 1);
	for (int i = 2; i <= A.n; i++)
	{
		LinA::Matrix Ai = GetColumn(A, i);
		for (int j = 1; j < i; j++)
		{
			LinA::Matrix Qj = GetColumn(Q, j);
			Ai = Ai - Qj * (Transpose(Ai) * Qj);
		}
		Ai = Ai * (1 / VecNorm_L2(Ai));
		ReplaceColumn(Q, Ai, i);
	}
	LinA::Matrix R = Transpose(Q) * A;
	return { Q,R };

}
LinA::Matrix LinA::Eig(LinA::Matrix A, float precision)
{
	assert(A.m == A.n);
	int count = 0;
	while (true)
	{
		auto Fi = QR(A);
		LinA::Matrix Qi = Fi.first;
		LinA::Matrix Ri = Fi.second;
		LinA::Matrix Aip = Ri * Qi;
		float eps = 0.f;
		for (size_t i = 0; i < A.n; i++)
		{
			eps = std::max(eps,abs(A.A[i][i] - Aip.A[i][i]));
		}
		
		A = Aip;
		count++;
		if (eps < precision) break;
	}
	LinA::Matrix eig(A.m, 1, 0);
	for (size_t i = 0; i < A.m; i++)
	{
		eig.A[i][0] = A.A[i][i];
	}
	//std::cout << "Eigenvalues: \n" << eig;
	std::cout << "Eigenvalue calculation using "<< count << " QR iterations to reach precision " << precision << "\n";
	return eig;
}
std::tuple<LinA::Matrix, LinA::Matrix, LinA::Matrix> LinA::LU(LinA::Matrix A, bool printLU)
{
	LinA::Matrix A_cpy = A;
	std::vector<LinA::Matrix> vecE; // elimination matrices container
	std::vector<LinA::Matrix> vecE_inverse;

	for (int i = 1; i <= A.m; i++)
	{
		std::pair<LinA::Matrix, LinA::Matrix> E_i = RowReduce(A, i);
		vecE.push_back(E_i.first);
		vecE_inverse.push_back(E_i.second);
		//std::cout << E_i.first;
	}

	LinA::Matrix L = LinA::Eye(A.n);
	//LinA::Matrix U = A;// _cpy;
	LinA::Matrix Elim = LinA::Eye(A.n);
	for (auto it = vecE_inverse.rbegin(); it != vecE_inverse.rend(); ++it)
	{
		//std::cout << "Multipltying " << *it;
		L = *it * L;
	}
	for (auto it = vecE.begin(); it != vecE.end(); ++it)
	{
		//std::cout << "Multipltying " << *it;
		Elim = *it * Elim;
	}

	if (printLU)
	{
		std::cout << "L\n" << L;
		std::cout << "U\n" << A; // note: A has been transformed to U
	}
	return std::make_tuple(L, A, Elim);
}
std::tuple<LinA::Matrix, LinA::Matrix> LinA::RREF(Matrix A)
{
	LinA::Matrix L, U, Elim;
	std::tie(L, U, Elim) = LinA::LU(A, false);

	auto pair = LinA::RowReduceUp(U, 3);
	Elim = pair.first * Elim;
	
	pair = LinA::RowReduceUp(U, 2);
	Elim = pair.first * Elim;
	

	// normalize diagonals
	LinA::Matrix S = LinA::Eye(3);
	S.A[0][0] = 1 / U.A[0][0];
	S.A[1][1] = 1 / U.A[1][1];
	S.A[2][2] = 1 / U.A[2][2];
	Elim = S * Elim;

	return std::make_tuple( Elim*A,Elim );
}
LinA::Matrix LinA::Inverse(LinA::Matrix A)
{
	LinA::Matrix Rref, Elim;
	std::tie(Rref, Elim) = LinA::RREF(A);
	return Elim;
}


LinA::Matrix LinA::SingularValues(LinA::Matrix A, float precision)
{
	// Estimates singular values based on eigenvalues of (At * A)
	assert(A.m == A.n);
	auto v = Eig(LinA::Transpose(A) * A, 1e-3);
	LinA::Matrix Sigma = LinA::Eye(A.n);
	for (size_t i = 0; i < v.m; i++)
	{
		Sigma.A[i][i] = sqrt(v.A[i][0]);
	}
	return Sigma;
}
float LinA::Det(const Matrix& A)
{
	LinA::Matrix L, U, Elim;
	std::tie(L, U, Elim) = LinA::LU(A, false);
	float DetU = LinA::DetOfTriangular(U);
	std::pair<LinA::Matrix, int> L_fixedPair = LinA::FixLowTriangular(L);
	float DetL = LinA::DetOfTriangular(L_fixedPair.first) * (1 + L_fixedPair.second % 2 * -2);
	std::cout << "Determinant of U = " << DetU << std::endl;
	std::cout << "Determinant of L = " << DetL << std::endl;
	std::cout << "Determinant of A = |LU| = |L||U| = " << DetL * DetU << std::endl;
	return DetL* DetU;
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
