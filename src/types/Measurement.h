#ifndef MEASUREMENT_H_
#define MEASUREMENT_H_

//#include <vector>

//#include "RandomVariable.h"


namespace SLOM {
template<typename T>
class IdxVector;



struct IMeasurement {
	/**
	 * eval(res) evaluates the function storing the result in res.
	 * The result is expected to be normalized i.e. having mean 0, and unit-covariance.
	 * Also it must write exactly getDim() values to the block starting at res.
	 */
	virtual double* eval(double* res) const = 0;

private:
	friend class Estimator;
	friend class IdxVector<IMeasurement>;
	/**
	 * getDim() returns the dimension of the result
	 */
	virtual int getDim() const = 0;
	/**
	 * getDepend() returns the sum of DOFs of the variables the measurement depends on
	 */
	virtual int getDepend() const = 0;
	
	/**
	 * registerVariables() registers the current measurement to the Variables it depends on.
	 * It returns the sum of DOFs of the variables, the measurement depends on.
	 */
	virtual int registerVariables() const = 0;
	
	/** 
	 * idx is the starting index of the Measurement in "the big matrix".
	 */
	int idx;
};

//typedef const IMeasurement* MeasId;



}  // namespace SLOM



#endif /*MEASUREMENT_H_*/
