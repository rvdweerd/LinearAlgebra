#include <iostream>
#include "Matrix.h"
#include <cassert>
#include <tuple>
//#include <Eigen/Dense>
//#include <Eigen/src/LU/InverseImpl.h>

void test1()
{
	{
		//LinA::Matrix A(" 1  1 1 ; 1 1 2 ; -1 1 3 "); // row exchanges test
		//LinA::Matrix A(" 1 4 -1 ; -2 1 1 ; -12 4 5 "); // no row exchanges
		LinA::Matrix A(" 2 -1 0 0 ; -1 2 -1 0 ; 0 -1 2 -1 ; 0 0 -1 2 "); // no row exchanges
		//LinA::Matrix A(" 1  1 1 ; 1 1 2 ; 2 2 4 "); // dependent columns
		std::cout << "A: \n" << A;

		auto lu = LU(A, false);
		std::cout << "LU decomposition of A.\n";
		std::cout << "L: \n" << lu.first;
		std::cout << "U: \n" << lu.second;
		std::cout << "L*U: \n" << lu.first * lu.second;

		auto qr = QR(A);
		std::cout << "QR decomposition of A.\n";
		std::cout << "Q: \n" << qr.first;
		std::cout << "R: \n" << qr.second;

		std::cout << "\nEigenvalues of A:\n";
		auto evals = Eig(A, 1e-3);
		std::cout << evals;

		std::cout << "\nSingular values of A (sqrt of eigenvalues of At * A):\n";
		auto svals = LinA::SingularValues(A, 1e-3);
		std::cout << "Sigma: \n" << svals;
	}
	std::cin.get();
}
void test2()
{
	LinA::Matrix a(" 1 ; 4  ");
	LinA::Matrix b(" 2 ; 1  ");
	auto f = ProjectVec(b, a);
}

//std::tuple<int, int, int> GetTuple()
//{
//	return std::make_tuple(3, 3, 3);
//}
int main()
{
	//int p, q, r;
	//std::tie(p,q,r) = GetTuple();

	test1();
	//LinA::Matrix A(" 1 4 -1 ; -2 1 1 ; -12 4 5 "); // no row exchanges
	//LinA::Matrix A(" 2 1 ; 2 -2 "); // no row exchanges
	//LinA::Matrix A(" 1 4 -1 ; -2 1 1 ; -1 5 0 "); // no row exchanges
	//std::cout <<"A= \n"<< A;
	//auto pair = LU(A);
	
	std::cin.get();


	//test1();
	//test2();
	
}