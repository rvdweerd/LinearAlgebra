#include <iostream>
#include "Matrix.h"
#include <cassert>
//#include <Eigen/Dense>
//#include <Eigen/src/LU/InverseImpl.h>

void test()
{
	{
		//LinA::Matrix A(" 1  1 1 ; 1 1 2 ; -1 1 3 ");
		LinA::Matrix A(" 1 4 -1 ; -2 1 1 ; -12 4 5 ");
		LinA::Matrix AT;
		AT = LinA::Transpose(A);
		//A = AT * A;
		LinA::Matrix A_cpy = A;

		std::cout << "Matrix A: " << A;

		std::cout << "Test3a: RowReduce() for column 1 of A: \n";
		std::pair<LinA::Matrix,LinA::Matrix> E1 = RowReduce(A, 1);
		std::cout << "Matrix A_new: " << A;
		std::cout << "Matrix E1: " << E1.first;
		std::cout << "Matrix E1*A: " << E1.first * A_cpy;
		std::cout << "Matrix E1_inv: " << E1.second;


		std::cout << "Test3b: RowReduce() for column 2 of A: \n";
		std::pair<LinA::Matrix, LinA::Matrix> E2 = RowReduce(A, 2);
		std::cout << "Matrix A_new: " << A;
		std::cout << "Matrix E2: " << E2.first;
		std::cout << "Matrix E2*E1*A: " << E2.first * (E1.first * A_cpy);
		std::cout << "Matrix E2_inv: " << E2.second;

		std::cout << "Test3c: RowReduce() for column 3 of A: \n";
		std::pair<LinA::Matrix, LinA::Matrix> E3 = LinA::RowReduce(A, 3);
		std::cout << "Matrix A_new: " << A;
		std::cout << "Matrix E3: " << E3.first;
		std::cout << "Matrix E3*E2*E1*A: " << E3.first * (E2.first * (E1.first * A_cpy));
		std::cout << "Matrix E3_inv: " << E3.second;

		LinA::Matrix E = E3.first * (E2.first * (E1.first));
		std::cout << "Matrix E: " << E;
		LinA::Matrix E_inv = E1.second * (E2.second * E3.second);
		std::cout << "Matrix E_inv = L: " << E_inv;
		LinA::Matrix U = E3.first * (E2.first * (E1.first * A_cpy));
		std::cout << "Matrix U: " << U;
		std::cout << "Matrix LU: " << E_inv*U;

		std::pair<LinA::Matrix, int> L_fixedPair = LinA::FixLowTriangular(E_inv);
		std::cout << "Fixed L with " << L_fixedPair.second << " row exchanges\n";
		std::cout << "L_fixed = " << L_fixedPair.first;

		float DetU = LinA::DetOfTriangular(U);
		float DetL = LinA::DetOfTriangular(L_fixedPair.first) * (1+ L_fixedPair.second%2 * -2);
		std::cout << "Determinant of U = " << DetU << std::endl;
		std::cout << "Determinant of L = " << DetL << std::endl;
		std::cout << "Determinant of A = |LU| = |L||U| = " << DetL*DetU << std::endl;



		//std::cout << "Test1a: RowOp(A,3,1,2) as free function on A: \n";
		//LinA::RowOp(A, 3, 1, 2);
		//std::cout << "Matrix A: " << A;

		//std::cout << "Test1b: invoke A.RowOp(2,2,3) as member function of A: \n";
		//A.RowOp(2, 2, 3);
		//std::cout << "Matrix A: " << A;
		//
		//std::cout << "Test2a: invoke A.RowExchange(1,2) as member function of A: \n";
		//A.RowExchange(1, 2);
		//std::cout << "Matrix A: " << A;

		//std::cout << "Test2b: RowExchange(A,2,3) as free function on A: \n";
		//LinA::RowExchange(A, 2, 3);
		//std::cout << "Matrix A: " << A;
	}
	std::cin.get();
}

std::pair<LinA::Matrix,LinA::Matrix> ProjectVec(LinA::Matrix b, LinA::Matrix a) // project vector b onto vector a, returns { p , e }  
{
	assert(a.n == 1 && b.n == 1);
	assert(a.m == b.m);
	LinA::Matrix aT = Transpose(a);
	float x_hat = (aT * b).A[0][0] / (aT * a).A[0][0];
	LinA::Matrix p = a * x_hat;
	LinA::Matrix e = b-p;
	return { p,e };
}


int main()
{
	LinA::Matrix a(" 1 ; 1  ");
	LinA::Matrix b(" 0 ; 1  ");
	auto f = ProjectVec(b,a);
	auto g = a * b;
}