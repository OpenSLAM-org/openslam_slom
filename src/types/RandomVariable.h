#ifndef RANDOMVARIABLE_H_
#define RANDOMVARIABLE_H_

#include <deque>

namespace SLOM {

template<typename T>
class IdxVector;

class IMeasurement;


struct IRVWrapper : private std::deque<const IMeasurement*>{
	IRVWrapper(bool optimize=true) : idx(-1), optimize(optimize), registered(false) {}
	virtual ~IRVWrapper() {}
	/**
	 * Gets the DOF of the enclosed RandomVariable.
	 */
	virtual int getDOF() const = 0;
	/**
	 * Adds vec to backup and stores result in var. backup is unaltered.
	 */
	virtual const double* add(const double* vec, double scale=1) = 0;
	/**
	 * Stores var in backup.
	 */
	virtual void store() = 0;
	/**
	 * Restores var from backup.
	 */
	virtual void restore() = 0;
	/**
	 * Registers the passed IMeasurement.
	 */
	int registerMeasurement(const IMeasurement* m){
		push_back(m);
		return registered ? getDOF() : 0;
	}
	
private:
	friend class IdxVector<IRVWrapper>;
	/// Requirements for IdxVector:
	//friend class IdxVector;
	/**
	 * The start index of this variable in the complete variable vector.
	 */
	int idx;
	/**
	 * Just a wrapper function for IdxVector.
	 */
	inline int getDim() const { return getDOF(); }
	/// \IdxVector
	
	friend class Estimator;
	/**
	 * Decide whether this variable should be optimized.
	 */
	bool optimize;
	bool registered;
};



/**
 * The templated class RVWrapper wraps RV to implement the interface IRVWrapper.
 * Requirements to RV are:
 * - An enum DOF which gives its degrees of freedom,
 * - A method add, which adds a scaled vector to the RV.
 * - Being CopyConstructable and Assignable. 
 */
template<typename RV>
class RVWrapper : public IRVWrapper{
	RV var;
	RV backup;
public:
	enum {DOF = RV::DOF};
	RVWrapper(const RV& v=RV(), bool optimize=true) : IRVWrapper(optimize), var(v), backup(v) {}
	int getDOF() const {return DOF;}
	const double* add(const double* vec, double scale=1) {
		var = backup;
		return var.add(vec, scale);
	}
	void store() {backup = var;}
	void restore() {var = backup;}

	// Getters and setters:
	const RV& operator*() const { return var; }
	const RV* operator->() const { return &var; } 
//	RV& operator*() { return var; }
//	RV* operator->() { return &var; } 
	const RV& operator=(const RV& v){
		return var = backup = v;
	}
};


}  // namespace SLSQ

#endif /*RANDOMVARIABLE_H_*/
