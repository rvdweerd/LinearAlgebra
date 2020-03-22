#include <iostream>
#include "Matrix.h"
//#include <Eigen/Dense>
//#include <Eigen/src/LU/InverseImpl.h>


int main()
{
	LinA::Matrix A(" 1 2 1 ; 3 8 1 ; 0 4 1");
	std::cout << "Matrix A: " << A;
	
	std::cout << "Test1a: RowOp(A,3,1,2) as free function on A: \n";
	LinA::RowOp(A, 3, 1, 2);
	std::cout << "Matrix A: " << A;

	std::cout << "Test1b: invoke A.RowOp(2,2,3) as member function of A: \n";
	A.RowOp(2, 2, 3);
	std::cout << "Matrix A: " << A;
	
	std::cout << "Test2a: invoke A.RowExchange(1,2) as member function of A: \n";
	A.RowExchange(1, 2);
	std::cout << "Matrix A: " << A;

	std::cout << "Test2b: RowExchange(A,2,3) as free function on A: \n";
	LinA::RowExchange(A, 2, 3);
	std::cout << "Matrix A: " << A;

	std::cin.get();
}