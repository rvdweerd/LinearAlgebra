#pragma once
#include <iostream>
#include "Matrix.h"
#include <cassert>
#include <tuple>

namespace Tests
{
	void TestEigenvalues(const LinA::Matrix& A)
	{
		std::cout << "======================================================\n";
		std::cout << "TEST: Eigenvalues (using QR iteration)\n";
		std::cout << "======================================================\n";
		std::cout << "A: \n" << A;
		std::cout << "\nEigenvalues of A:\n";
		auto evals = Eig(A, 1e-3);
		std::cout << evals;
		std::cin.get();
	}
	void TestSingularValues(const LinA::Matrix& A)
	{
		std::cout << "======================================================\n";
		std::cout << "TEST: Singular values (naive based on At * A)\n";
		std::cout << "======================================================\n";
		std::cout << "A: \n" << A;
		std::cout << "\nSingular values of A (sqrt of eigenvalues of At * A):\n";
		auto svals = LinA::SingularValues(A, 1e-3);
		std::cout << "Sigma: \n" << svals;
		std::cin.get();
	}
	void TestVectorProjections()
	{
		std::cout << "======================================================\n";
		std::cout << "TEST: Vector projection\n";
		std::cout << "======================================================\n";
		LinA::Matrix a(" 1 ; 4 ; -3 ");
		LinA::Matrix b(" 2 ; 1 ; 2 ");
		std::cout << "aT :\n" << LinA::Transpose(a);
		std::cout << "bT :\n" << LinA::Transpose(b);
		auto pair = ProjectVec(b, a);
		std::cout << "pT (projection of b onto a):\n" << LinA::Transpose(pair.first);
		std::cout << "eT (error vector e = b - p):\n" << LinA::Transpose(pair.second);
		std::cin.get();
	}
	void TestQR(LinA::Matrix A)
	{
		std::cout << "======================================================\n";
		std::cout << "TEST: QR Decomposition\n";
		std::cout << "======================================================\n";
		std::cout << "A: \n" << A;
		auto qr = QR(A);
		std::cout << "QR decomposition of A.\n";
		std::cout << "Q: \n" << qr.first;
		std::cout << "R: \n" << qr.second;
		std::cout << "Q*R: \n" << qr.first * qr.second;
		std::cin.get();
	}
	void TestLU(LinA::Matrix A)
	{
		std::cout << "======================================================\n";
		std::cout << "TEST: LU Decomposition\n";
		std::cout << "======================================================\n";
		std::cout << "A: \n" << A;
		LinA::Matrix L, U, Elim;
		std::tie(L, U, Elim) = LU(A, false);
		std::cout << "LU decomposition of A.\n";
		std::cout << "L: \n" << L;
		std::cout << "U: \n" << U;
		std::cout << "L*U: \n" << L * U;
		std::cin.get();
	}

	void TestInverse(const LinA::Matrix& A)
	{
		std::cout << "======================================================\n";
		std::cout << "TEST: INVERSE\n";
		std::cout << "======================================================\n";
		std::cout << "A:\n" << A;
		LinA::Matrix Ainv = LinA::Inverse(A);
		std::cout << "inv(A):\n" << Ainv;
		std::cout << "A*inv(A):" << A * Ainv;
		std::cout << "inv(A)*A:" << Ainv * A;
		std::cin.get();
	}
	void TestRREF(const LinA::Matrix& A)
	{
		std::cout << "======================================================\n";
		std::cout << "TEST: RREF\n";
		std::cout << "======================================================\n";
		std::cout << "A:\n" << A;
		LinA::Matrix Rref, Elim;
		std::tie(Rref, Elim) = LinA::RREF(A);
		std::cout << "rref(A):\n";
		std::cout << Rref;
		std::cout << "Elimination matrix E ( = inv(A) when A is invertible):\n";
		std::cout << Elim;
		std::cin.get();
	}
}