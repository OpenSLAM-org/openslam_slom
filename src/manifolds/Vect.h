#ifndef VECT_H_
#define VECT_H_

namespace SLOM {


template<int D>
struct Vect : public Manifold<Vect<D>,D>{
	double data[D];
	
	Vect(){
		for(int i=0; i<D; i++) data[i]=0;
	}
	
	Vect(const double* src){
		for(int i=0; i<D; i++)
			data[i] = *src++;
	}
	
	void add_(const double vec[D], double scale=1){
		for(int i=0; i<D; i++){
			data[i]+= scale * *vec++;
		}
	}
	void sub_(double res[D], const Vect<D>& oth) const {
		for(int i=0; i<D; i++){
			*res++ = data[i]-oth.data[i];
		}
	}
	
	double& operator[](int idx) {return data[idx]; }
	const double& operator[](int idx) const {return data[idx]; }
	
};


}  // namespace SLOM


#endif /*VECT_H_*/
