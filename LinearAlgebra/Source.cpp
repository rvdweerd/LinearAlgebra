#include <iostream>
#include "Matrix.h"
#include <cassert>
#include <tuple>
#include "Tests.h"

int main()
{
	// TEST MATRICES
	LinA::Matrix A(" 1 4 -1 ; -2 1 1 ; -12 4 5 "); // no row exchanges
	//LinA::Matrix A(" 2 1 ; 2 -2 "); // no row exchanges
	//LinA::Matrix A(" 1 4 -1 ; -2 1 1 ; -1 5 0 "); // no row exchanges
	//LinA::Matrix A(" 1  1 1 ; 1 1 2 ; -1 1 3 "); // row exchanges test
	//LinA::Matrix A(" 1  1 1 ; 1 1 2 ; 2 2 4 "); // dependent columns

	Tests::TestLU(A);
	Tests::TestRREF(A);
	Tests::TestInverse(A);

	Tests::TestVectorProjections();
	Tests::TestQR(A);
	Tests::TestEigenvalues(A);
	Tests::TestSingularValues(A);


	
}