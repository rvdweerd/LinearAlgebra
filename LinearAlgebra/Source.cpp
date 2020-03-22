#include <iostream>
#include "Matrix.h"
//#include <Eigen/Dense>
//#include <Eigen/src/LU/InverseImpl.h>


int main()
{
	{
		LinA::Matrix A(" 1 2 1 ; 3 8 1 ; 1 4 1 ; 1 1 1");
		//A=LinA::Transpose(A);
		LinA::Matrix A_cpy = A;

		std::cout << "Matrix A: " << A;

		std::cout << "Test3a: RowReduce() for column 1 of A: \n";
		LinA::Matrix E1 = RowReduce(A, 1);
		LinA::Matrix E1_inv = LinA::InverseEliminationMatrix(E1);
		std::cout << "Matrix A_new: " << A;
		std::cout << "Matrix E1: " << E1;
		std::cout << "Matrix E1*A: " << E1 * A_cpy;
		std::cout << "Matrix E1_inv: " << E1_inv;


		std::cout << "Test3b: RowReduce() for column 2 of A: \n";
		LinA::Matrix E2 = RowReduce(A, 2);
		LinA::Matrix E2_inv = LinA::InverseEliminationMatrix(E2);
		std::cout << "Matrix A_new: " << A;
		std::cout << "Matrix E2: " << E2;
		std::cout << "Matrix E2*E1*A: " << E2 * (E1 * A_cpy);
		std::cout << "Matrix E2_inv: " << E2_inv;

		std::cout << "Test3c: RowReduce() for column 3 of A: \n";
		LinA::Matrix E3 = RowReduce(A, 3);
		LinA::Matrix E3_inv = LinA::InverseEliminationMatrix(E3);
		std::cout << "Matrix A_new: " << A;
		std::cout << "Matrix E3: " << E3;
		std::cout << "Matrix E3*E2*E1*A: " << E3 * (E2 * (E1 * A_cpy));
		std::cout << "Matrix E2_inv: " << E3_inv;

		LinA::Matrix E = E3 * (E2 * (E1));
		std::cout << "Matrix E: " << E;
		LinA::Matrix E_inv = E1_inv * (E2_inv * E3_inv);
		std::cout << "Matrix E_inv = L: " << E_inv;
		LinA::Matrix U = E3 * (E2 * (E1 * A_cpy));
		std::cout << "Matrix U: " << U;
		std::cout << "Matrix LU: " << E_inv*U;


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