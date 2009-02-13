#ifndef SON_H_
#define SON_H_

#include "../types/Manifold.h"
#include "../tools/quaternion.h"
#include "Vect.h"

namespace SLOM {


// Special Orthogonal Groups

/**
 * Overloading operators for rotations.
 * Operators are overloaded as if rotations were represented as matrix.
 * A/B calculates A*B^-1, A%B calculates A^-1*B (like e.g. Matlabs '\')
 */
template<typename T, int N>
struct RotationGroup {
	typedef Vect<N> VectN;
	
	T& operator*=(const T& oth){
		T& ths = static_cast<T&>(*this);
		ths.mult(oth, false, false);
		return ths;
	}
	T& operator/=(const T& oth){
		T& ths = static_cast<T&>(*this);
		ths.mult(oth, false, true);
		return ths;
	}
	T& operator%=(const T& oth){
		T& ths = static_cast<T&>(*this);
		ths.mult(oth, true, false);
		return ths;
	}
	friend T operator*(const T& ths, const T& oth) {
		T res(ths);
		return res*=oth;
	}
	friend T operator/(const T& ths, const T& oth) {
		T res(ths);
		return res/=oth;
	}
	friend T operator%(const T& ths, const T& oth) {
		T res(ths);
		return res%=oth;
	}
	friend VectN operator*(const T& r, const VectN &vec) {
		VectN res;
		r.rotate(res, vec, false);
		return res;
	}
	friend VectN operator%(const T& r, const VectN &vec) {
		VectN res;
		r.rotate(res, vec, true);
		return res;
	}
};




/**
 * Rotations in 3D.
 * In this implementation rotations are represented as quaternions.
 */
struct SO3 : public Manifold<SO3, 3>, public RotationGroup<SO3, 3>{
	Quaternion quat;
	// Requirements for Manifold:
	void add_(const double vec[3], double scale=1){
		Quaternion q(vec, scale);
		quat*=q;
	}
	void sub_(double res[3], const SO3& oth) const{
		Quaternion tmp = quat % oth.quat; // this^-1 * oth
		tmp.toScaledAxis(res);
	}
	///
	
	SO3(const Quaternion &q=Quaternion()) : quat(q) {}
	
	
	/**
	 * multiplies (*this) by oth, inverting this or oth priorly if requested.
	 */
	void mult(const SO3& oth, bool invThis, bool invOth){
		quat.multiply(oth.quat, invThis, invOth);
	}

	/**
	 * Rotates a 3D Vector (backwards if requested).
	 */
	void rotate(Vect<3> &res, const Vect<3> &vec, bool back=false) const {
		quat.rotate(res.data, vec.data, back);
	}
};


	

/**
 * Rotations in 2D (simple Rotations).
 * 
 */
struct SO2 : public Manifold<SO2, 1>, public RotationGroup<SO2,2>{
	double angle;
	SO2(double angle = 0) : angle(angle) {}
	SO2(const Vect<2> &dir) : angle(atan2(dir[1], dir[0])) {}
	
	void add_(const double vec[1], double scale=1){
		angle += scale*vec[0];
	}
	void sub_(double res[1], const SO2& oth) const{
		res[0] = normalize(angle-oth.angle);
	}
	
	void mult(const SO2& oth, bool invThis, bool invOth){
		if(invThis) angle = -angle;
		angle += invOth ? -oth.angle : oth.angle;
	}
	
	void rotate(double res[2], const double vec[2], bool back=false) const{
		double c=cos(angle), s= back ? -sin(angle) : sin(angle);
		double x=vec[0], y = vec[1];
		res[0] = c*x - s*y;
		res[1] = s*x + c*y;
	}

	void rotate(Vect<2> &res, const Vect<2> &vec, bool back=false) const{
		rotate(res.data, vec.data, back);
	}
	
	operator double() const {return angle;}
private:
	static inline double normalize(double x){
		if(fabs(x) <= M_PI) return x;
		int r = (int)(x*M_1_PI);
		return x - ((r + (r>>31) + 1) & ~1)*M_PI; 
	}
};




}  // namespace SLOM

#endif /*SON_H_*/

