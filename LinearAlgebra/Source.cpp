#include <iostream>
#include "Matrix.h"
#include <Eigen/Dense>
#include <Eigen/src/LU/InverseImpl.h>


int main()
{
	LinA::Matrix A(" 1 2 1 ; 3 8 1 ; 0 4 1");
	
	
	std::cout << "Matrix A: " << A;
	A.RowOp(3, 1, 2);
	std::cout << "Matrix A: " << A;
	A.RowOp(2, 2, 3);
	std::cout << "Matrix A: " << A;

	std::cin.get();
}