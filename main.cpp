// Importing Eigen library and trying transpose and null space of a matrix

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

int main() {
    MatrixXf A{6, 3};
    A << 0, 0, 1, 0, 1, 0, 1, 0, 0, -1, 0, 0, 0, 0, -1, 0, -1, 0;
    MatrixXf B = A;
    std::cout << "Here is the initial matrix m:\n" << A << std::endl;
    MatrixXf b{6, 1};
    b << 1, 1, 1, 1, 1, 1;
    std::cout << "Here is the initial matrix d:\n" << b << std::endl;

    std::cout << "The matrix m is of size "
              << A.rows() << "x" << A.cols() << std::endl;
    std::cout << "It has " << A.size() << " coefficients" << std::endl;

    A.transposeInPlace();
    std::cout << "and after being transposed:\n" << A << std::endl;


    FullPivLU<MatrixXf> lu(A);
    MatrixXf A_null_space = lu.kernel();
    std::cout << "Null space:\n" << A_null_space << std::endl;
    A_null_space.transposeInPlace();
    std::cout << "Null space Transposed_A:\n" << A_null_space << std::endl;

    /*
    FullPivLU<MatrixXf> luB(A);
    MatrixXf B_null_space = luB.kernel();
    B_null_space.transposeInPlace();
    std::cout << "Null space Transposed_B:\n" << B_null_space << std::endl;
     */

    return 0;
}

/*
    --  OUTPUT  --
    
    /Users/soner/CLionProjects/proje/cmake-build-debug/proje
Here is the initial matrix m:
 0  0  1
 0  1  0
 1  0  0
-1  0  0
 0  0 -1
 0 -1  0
Here is the initial matrix d:
1
1
1
1
1
1
The matrix m is of size 6x3
It has 18 coefficients
and after being transposed:
 0  0  1 -1  0  0
 0  1  0  0  0 -1
 1  0  0  0 -1  0
Null space:
-0  1 -0
-0 -0  1
 1 -0 -0
 1  0  0
 0  1  0
 0  0  1
Null space Transposed_A:
-0 -0  1  1  0  0
 1 -0 -0  0  1  0
-0  1 -0  0  0  1

Process finished with exit code 0

*/
