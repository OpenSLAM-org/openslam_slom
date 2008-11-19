#ifndef MANIFOLD_H_
#define MANIFOLD_H_


namespace SLOM {


template<class Derived, int D>
struct Manifold
{
	enum{DOF=D};
	int getDOF() const {return DOF;}
	const double * add(const double *vec, double scale=1){
		static_cast<Derived*>(this)->add_(vec, scale);
		return vec+D;
	}
	double * sub(double *res, const Derived& oth) const{
		static_cast<const Derived*>(this)->sub_(res, oth);
		return res+D;
	}
};
	


}  // namespace SLSQ

#endif /*MANIFOLD_H_*/
