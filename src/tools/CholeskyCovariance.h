#ifndef CHOLESKYCOVARIANCE_H_
#define CHOLESKYCOVARIANCE_H_

#include <algorithm>
#include <numeric>
#include <cassert>

namespace SLOM {



namespace CholeskyMode{
/**
 * How to initialize a CholeskyCovariance.
 * COPY_* means source already is a Cholesky factor,
 * CHOLESKY_* does the decomposition itself.
 * *_UPPER means only the upper half of the matrix is stored,
 * *_FULL means zeros or symmetric values from lower half are stored.
 */
enum CM {
	COPY_UPPER,
	COPY_FULL,
	CHOLESKY_UPPER,
	CHOLESKY_FULL
};
};

template<int dim>
struct CholeskyCovariance
{
	enum {DIM = dim, SIZE = (dim*(dim+1))/2};
	// chol is a lower triangular matrix, saved in row major order
	double chol[SIZE];


	CholeskyCovariance() {};

	CholeskyCovariance(const double *A, CholeskyMode::CM mode){
		bool skip = false;
		switch(mode){
		case CholeskyMode::COPY_UPPER:
			skip = true;
		case CholeskyMode::COPY_FULL:
			copyCholesky(A, skip);
			break;
		case CholeskyMode::CHOLESKY_UPPER:
			skip = true;
		case CholeskyMode::CHOLESKY_FULL:
			calculateCholesky(A, skip);
			break;
		default:
			assert(false); //Unknown CholeskyMode
		}
	}
	/**
	 * Initialize by copying values from data.
	 */
	void copyCholesky(const double *A, bool skip){
		if(skip){
			std::copy(A, A+SIZE, chol);
			return;
		}
		const double *a = A;
		double *coli = chol;
		for(int i=1; i<=DIM; i++){
			for(int j=0; j<i; j++){
				*coli++ = *a++;
			}
			a += DIM-i;
		}
		assert(a==A+DIM*DIM);
		assert(coli == chol+SIZE);
	}

	/**
	 * Initialize by Cholesky decomposition of matrix A.
	 * A is stored row major and only the upper half of A is considered.
	 * If parameter skip == true, it is assumed, that only the upper half of A is stored.
	 * I.e: A = [a00 a01 .. a0n a11 a12 .. a1n a22 ...]
	 * (As A is symmetric, you can substitute "row major" by "column major"
	 * and every "upper" by "lower")
	 */
	void calculateCholesky(const double *A, bool skip){
		const double *a = A;
		double *coli = chol;
		for(int i=0; i<DIM; i++){
			if(!skip){ // skip the first elements of this row in A:
				a += i;
			}
			coli += i; // start of column i
			// C(i,i) = sqrt(A(i,i) - C(0:i,i)'*C(0:i,i))
			double sum = *a++ - std::inner_product(coli, coli+i, coli, 0.0 );
			assert(sum>0);
			double sqrtSum = coli[i] = sqrt(sum);
			sqrtSum = 1/sqrtSum;
			double *colj = coli;
			for(int j=i+1; j<DIM; j++){
				colj += j;
				// C(i,j) = (A(i,j) - C(0:i,i)'*C(0:i,j))/C(i,i)
				colj[i] = (*a++ - std::inner_product(coli, coli+i, colj, 0.0)) * sqrtSum;
			}
		}
		assert(skip ? (a==A+SIZE) : (a==A+DIM*DIM));
	}

	/**
	 * multiplies arr[0,dim) by inverse of this Cholesky factor.
	 * I.e. computes chol\arr;
	 */
	void invApply(double* arr) const {
		const double *c = chol;
		for(double *a = arr; a < arr + DIM; a++){
			double ai = *a;
			for(double *b=arr; b<a; b++){
				ai -= *c++ * *b;
			}
			*a = ai / *c++;
		}
	}

	/**
	 * multiplies arr[0,dim) by this Cholesky factor.
	 * I.e. computes chol*arr;
	 */
	void apply(double *arr) const {
		const double *c = chol + SIZE -1;
		for(double *a=arr + DIM - 1; a >= arr; --a){
			*a *= *c--;
			for(double *b = a-1; b>= arr; b--){
				*a += *c-- * *b;
			}
		}

	}
};

}  // namespace SLOM

#endif /*CHOLESKYCOVARIANCE_H_*/
