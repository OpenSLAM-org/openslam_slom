#include "Estimator.h"

#include <algorithm>
#include <cmath>
#include <cassert>

#include <iostream>


//#include "tools/cs_extension.h"

namespace SLOM {



void Estimator::freeWorkspace(){
	jacobian=cs_spfree(jacobian);
	JtJ = cs_spfree(JtJ);
	symbolic=cs_sfree(symbolic);
	numeric=cs_nfree(numeric);

	delete[] res;            res = 0;
	delete[] workspace;      workspace = 0;
	delete[] cholCovariance; cholCovariance = 0;
}


Estimator::~Estimator(){
	freeWorkspace();
}


double Estimator::evaluate(double * result) const{
	IdxVector<IMeasurement>::const_iterator it = measurements.begin();
	double sum = 0;
	double *ptr = result;
	for(; it != measurements.end(); it++){
		assert(result + (*it)->idx == ptr);
		(*it)->eval(ptr);
		double *end = ptr+(*it)->getDim();
		for( ; ptr < end; ptr++){
			sum += std::pow(*ptr, 2);
		}
	}
	return sum;
}


void Estimator::createSparse(){
	freeWorkspace();
	int M = measurements.getDim();
	int N = variables.getDim();
	int size = nnz;
	int skip = 0;
	switch(usedAlgorithm){
	case GaussNewton: break;
	// For Levenberg(-Marquardt) put a diagonal matrix below the actual Jacobian:
	case Levenberg:
	case LevenbergMarquardt:
		size += N;
		M += N;
		skip = 1;
	}
	//allocate a compressed column matrix:
	jacobian = cs_spalloc(M, N,   // size
	                    size,    // number of non-zeros
	                    true,   // allocate space for values
	                    false); // compressed-column

	//create the structure of the matrix:
	int *cIdx=jacobian->p; // column pointer
	int *rIdx=jacobian->i; // row indices
	*cIdx = 0; //first column starts at 0.
	int n=measurements.getDim();
	for(IdxVector<IRVWrapper>::const_iterator v = variables.begin(); v!= variables.end(); v++){
		const IRVWrapper* var = *v;
		int vDOF = var->getDOF();
		assert(vDOF>0);
		if(var->begin()==var->end()){
			std::cerr << v - variables.begin() << std::endl;
			assert(var->begin() != var->end()); // No measurements

		}
		for(IRVWrapper::const_iterator meas= var->begin(); meas!= var->end(); meas++){
			int row = (*meas)->idx;
			for(int mDim = (*meas)->getDim(); mDim>0; mDim--){
				*rIdx++ = row++;
			}
		}
		if(usedAlgorithm != GaussNewton){
			*rIdx++ = n++;
		}
		*++cIdx = rIdx - jacobian->i; //start of next column
		while(--vDOF > 0){ // copy the current column for each DOF of the variable
			rIdx = std::copy(jacobian->i + cIdx[-1], // start of previous column
					rIdx-skip, rIdx);
			if(usedAlgorithm != GaussNewton){
				*rIdx++ = n++;
			}
			*++cIdx = rIdx - jacobian->i;
		}
	}
	//*cIdx = rIdx - matrix->i; //end of last column;
	assert(*cIdx == size);

	switch(usedSolver){
	case QR:
		symbolic = cs_sqr(3, jacobian, true);
		assert(symbolic);
		break;
	case Cholesky:
		break;
	}
	res=new double[M];
	workspace = new double[M];
}


void Estimator::initCovariance(){
	delete[](cholCovariance);
	//cs_spfree(cholCovariance);
	int n=variables.getDim();
	cholCovariance = new double[n];
	std::fill_n(cholCovariance, n, 0);
}

void Estimator::calculateJacobian(){
	assert(jacobian);   // matrix is allocated
	assert(res);
	assert(cholCovariance);
	int m = jacobian->m, n = jacobian->n;
	double *x = jacobian->x;

	int skip=0;
	switch(usedAlgorithm){
	case GaussNewton:
		assert(m == measurements.getDim());
		assert(n == variables.getDim());
		skip = 0;
		break;
	// For Levenberg(-Marquardt) put a diagonal matrix below the Jacobian:
	case Levenberg:
	case LevenbergMarquardt:
		assert(m-n == measurements.getDim());
		assert(n == variables.getDim());
		skip = 1;
		break;
	}

	double *chol = cholCovariance;

	for(IdxVector<IRVWrapper>::iterator v = variables.begin(); v!= variables.end(); v++){
		IRVWrapper* var = *v;
		int vDOF = var->getDOF();
		assert(vDOF>0);
		double add[vDOF]; // temp-array for adding
		std::fill_n(add, vDOF, 0);
		for(int k=0; k<vDOF; k++){
			double d = 1e6;
			add[k] = (d ? 1/d : 0);

			// store $f(\mu \mplus 1/d)$ in res:
			var->add(add);
			double *temp=workspace;
			for(IRVWrapper::const_iterator meas= var->begin(); meas!= var->end(); meas++){
				temp=(*meas)->eval(temp);
			}
			var->restore();

			// store $f(\mu \mplus -1/d)$ directly in the matrix:
			var->add(add, -1);
			temp = x;
			for(IRVWrapper::const_iterator meas= var->begin(); meas!= var->end(); meas++){
				temp=(*meas)->eval(temp);
			}
			var->restore();

			// calculate difference and multiply by $0.5d$
			// also accumulate results for new inverse covariance
			double c = 0;
			for(double *xP=workspace ;x<temp; x++){
				*x = 0.5*d*(*xP++ - *x);
				assert(std::isfinite(*x));
				c += std::pow(*x,2);
			}
			*chol++ = std::sqrt(c);
			add[k] = 0; // reset delta-vector
			x += skip; // skip one entry for Levenberg(-Marquardt)
		}

	}



}

void Estimator::updateDiagonal() {
	if(usedAlgorithm == GaussNewton) return;
	// for LMA set the last entry of each column to lamda or lamda*cholCovariance;
	int n=jacobian->n;
	int *p = jacobian->p;
	double *x = jacobian->x;

	for(int k=0; k<n; k++){
		x[p[k+1]-1] = (usedAlgorithm == Levenberg ?
				lamda : lamda*cholCovariance[k]);
	}

}

void Estimator::updateSparse(){
	calculateJacobian();
	updateDiagonal();
}

void Estimator::qrSolve(double* delta){
	assert(symbolic);
	cs_nfree(numeric);
	numeric = cs_qr(jacobian, symbolic);
	assert(numeric);
	int m=jacobian->m, n=jacobian->n;

	// The following code essentially does a cs_qrsol(3, matrix, workspace);
	// but it doesn't recalculate the symbolic decomposition

	cs_ipvec(symbolic->pinv, res, workspace, m) ; /* x(0:m-1) = b(p(0:m-1) */
	for (int k = 0; k < n; k++) /* apply Householder refl. to x */
	{
		cs_happly(numeric->L, k, numeric->B [k], workspace) ;
	}
	cs_usolve(numeric->U, workspace) ; /* x = R\x */

	cs_ipvec(symbolic->q, workspace, delta, n) ; /* b(q(0:n-1)) = x(0:n-1) */
	// end of qrsol
}

void Estimator::choleskySolve(double *delta){
	// TODO reduce number of copies and mallocs.
	// TODO reuse symbolic Cholesky decomp.
	// TODO reuse structure of JtJ, only store upper half
	cs* Jt = cs_transpose(jacobian, true);
	cs_spfree(JtJ);
	JtJ = cs_multiply(Jt, jacobian);
	assert(Jt && JtJ);
	std::fill(workspace, workspace+jacobian->m, 0);
	cs_gaxpy(Jt, res, workspace);
	Jt = cs_spfree(Jt);
	int ok=cs_cholsol(1, JtJ, workspace);
	std::copy(workspace, workspace + jacobian->n, delta);
	JtJ = cs_spfree(JtJ);
	assert(ok && !Jt && !JtJ);
}


void Estimator::initialize(){
	freeWorkspace();
	createSparse();
	initCovariance();
}

double Estimator::optimizeStep(){
	// TODO better parameter control for LMA
	assert(jacobian);

	int m=jacobian->m, n=jacobian->n;

	updateSparse();


	//std::fill(workspace, workspace+m, 0);
	double sum=evaluate(res);
	std::cout << sum;
	if(usedAlgorithm != GaussNewton){
		std::fill(res + m-n, res+m, 0);
	}

	double* delta = usedAlgorithm == GaussNewton ? res : res + m-n;

	switch(usedSolver){
	case QR:
		qrSolve(delta);
		break;
	case Cholesky:
		choleskySolve(delta);
		break;
	}

	const double *temp=delta;
	for(IdxVector<IRVWrapper>::iterator v = variables.begin(); v!= variables.end(); v++){
		IRVWrapper* var = *v;
		if(var->optimize){
			temp = var->add(temp, -1);
		} else {
			temp += var->getDOF();
		}
		//var->store();
	}
	double sum2 = evaluate(workspace);
	double gain = (sum - sum2)/sum2;
	std::cout << ", RSS: " << sum2 << ", RMS: " << std::sqrt(sum2/m) << ", Gain: " << gain;
	if(gain > 0 || usedAlgorithm == GaussNewton){
		for(IdxVector<IRVWrapper>::iterator v = variables.begin(); v!= variables.end(); v++){
			(*v)->store();
		}
		std::swap(workspace, res);
		lamda *= sqrt(0.1);
	} else {
		for(IdxVector<IRVWrapper>::iterator v = variables.begin(); v!= variables.end(); v++){
			(*v)->restore();
		}
		lamda *= sqrt(10.0);
	}
	if(usedAlgorithm != GaussNewton){
		std::cout << ", lamda = " << lamda;
	}
	std::cout << std::endl;
	return gain;
}

}  // namespace SLOM
