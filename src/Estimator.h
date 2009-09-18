#ifndef ESTIMATOR_H_
#define ESTIMATOR_H_


#include "types/RandomVariable.h"
#include "types/Measurement.h"

#include "types/IdxVector.h"

#include <cs.h>

namespace SLOM {

class Estimator
{
public:
	enum Algorithm{ // use the following "damping term" $N$ in $(J^T J + N)\delta = J^T [y - f(\beta)]}$
		GaussNewton,       // $N=0$
		Levenberg,         // $N= \lambda*I$
		LevenbergMarquardt // $N= \lambda*diag(J^T J)$
	};
	
	enum Solver{
		QR,
		Cholesky
	};
private:
	
	enum Algorithm usedAlgorithm;
	enum Solver usedSolver;
	
	// input data:
	IdxVector<IRVWrapper> variables;
	IdxVector<IMeasurement> measurements;
	
	
	// the big matrix:
	int nnz; // number of non-zeroes in Jacobian
	cs* jacobian;
	cs* JtJ;  // $J^T J$ for CholeskySolve, this is also the information matrix
	
	css* symbolic; // symbolic decomposition of jacobian or JtJ 
	csn* numeric;  // numeric decomposition of jacobian or JtJ
	
	// the current residuum:
	double *res;
	
	// workspace for solving:
	double *workspace;
	
	
	// lamda parameter for LMA:
	double lamda; 
	
	
	/** the cholesky factor of the current covariance.
	 */
	double* cholCovariance;
	
	/** the last Residual Sum of Squares
	 */
	double lastRSS;

	
	/**
	 * Creates "the big matrix", 
	 */
	void createSparse();
	void updateSparse();
	void updateDiagonal();
	void freeWorkspace();
	void qrSolve(double* delta);
	void choleskySolve(double* delta);

	/**
	 * Numerically calculates the Jacobian of the function with the 
	 * and updates cholCovariance. The result is stored in matrix.
	 */
	void calculateJacobian();
	
	void initCovariance();
	
public:
	
	Estimator(Algorithm alg=GaussNewton, double lamda0=1e-3) : 
		usedAlgorithm(alg), usedSolver(Cholesky),
		nnz(0), jacobian(0), JtJ(0), symbolic(0), numeric(0), res(0), workspace(0), lamda(lamda0), cholCovariance(0) {};

	Estimator(Solver solver, Algorithm alg=GaussNewton, double lamda0=1e-3) : 
		usedAlgorithm(alg), usedSolver(solver),
		nnz(0), jacobian(0), JtJ(0), symbolic(0), numeric(0), res(0), workspace(0), lamda(lamda0), cholCovariance(0) {};
		
	
	
	virtual ~Estimator();
	
	typedef const IRVWrapper* RVId;
	
	/**
	 * Adds a copy of var into the Estimator. Returns an id.
	 * DEPRECATED
	 */ 
	//RVId insertRV(const IRandomVar &var); //TODO
	
	/**
	 * Adds the RandomVar *var itself to the Estimator. 
	 * The user is responsible for data holding. 
	 */
	RVId insertRV(IRVWrapper *var){
		if(var->optimize){
			var->registered=true;
			variables.push_back(var);
		}
		return var;
	}
	
	/**
	 * Removes the RandomVar with id from the estimator.
	 * Returns true on success. If the RV is still used by measurements, 
	 * it returns false.
	 */
	bool removeRV(const RVId id); //TODO
	
	/**
	 * Inserts a new measurement. 
	 * The Measurement itself registers the variables it depends on. 
	 */
	void insertMeasurement(IMeasurement* meas){
		//MeasId id=
		measurements.push_back(meas);
		int depend = meas->registerVariables();
		nnz += depend * meas->getDim();
		//return id;
	}
	
	void printJacobian(bool brief=false) const {
		cs_print(jacobian, brief);
	}
	
	double evaluate(double *result) const;
	
	void initialize();
	
	/**
	 * Optimizes 
	 */
	double optimizeStep(); //TODO parameters
	
	int getM() const {
		return measurements.getDim();
	}
	
	int getN() const {
		return variables.getDim();
	}
	
	double getLastRSS() const {
		return lastRSS;
	}
	
	const double * getCholCovariance() const {
		return cholCovariance;
	}
	
	void changeAlgorithm(Algorithm algo, double lamdaNew=-1){
		usedAlgorithm = algo;
		if(lamdaNew > 0) lamda = lamdaNew;
		freeWorkspace();
		initialize();
	}
	
	
};

}  // namespace SLOM

#endif /*ESTIMATOR_H_*/
