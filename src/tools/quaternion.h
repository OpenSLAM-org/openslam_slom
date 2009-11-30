#ifndef QUATERNION_H_
#define QUATERNION_H_

#include <cmath>


struct Quaternion{
	double w,x,y,z;
	
	
	//! Default constructor
	Quaternion(){w=1; x=y=z=0;}
	
	//! Init by scaled axis: scale*vec(1:3)
	Quaternion(const double* vec, double scale=1){
		x=*vec++, y=*vec++, z=*vec;
		double norm = sqrt(x*x + y*y + z*z);
		
		double alpha = scale*norm;
		if(norm < 1e-12) norm = 1e-12;
		
		w = cos(0.5*alpha);
		
		double mult = sin(0.5*alpha)/norm;
		
		x*=mult; y*=mult; z*=mult;
	}
	
	Quaternion& conj(){
		x = -x; y = -y; z = -z;
		return *this;
	}
	
	/**
	 * Convert Euler-Angle to Quaternion.
	 * Code based on: 
	 * http://www.euclideanspace.com/maths/geometry/rotations/conversions/eulerToQuaternion/
	 */
	void rollPitchYaw(double r, double p, double j){
		// TODO check order of rotations!!!
		double c1=cos(r/2), c2=cos(p/2), c3=cos(j/2);
		double s1=sin(r/2), s2=sin(p/2), s3=sin(j/2);
		
		w = c1 * c2 * c3 - s1 * s2 * s3;
		x = s1 * s2 * c3 + c1 * c2 * s3;
		y = s1 * c2 * c3 + c1 * s2 * s3;
		z = c1 * s2 * c3 - s1 * c2 * s3;
	}
	
	/**
	 * Multiplication routine.
	 * q is the second factor, 
	 * c1 and c2 determine if either factor is to be conjugated.
	 */
	Quaternion& multiply(const Quaternion& q, bool c1=false, bool c2=false){
		double a=w,  b=x,  c=y,  d=z;
		double e=q.w,f=q.x,g=q.y,h=q.z;
		if(c1){
			b=-b; c=-c; d=-d;
		}
		if(c2){
			f=-f; g=-g; h=-h;
		}
		
		w = a*e - b*f - c*g - d*h;
		x = a*f + b*e + c*h - d*g;
		y = a*g - b*h + c*e + d*f;
		z = a*h + b*g - c*f + d*e;
		
		return *this;
	}

	Quaternion& operator*=(const Quaternion& q){
		return multiply(q, false, false);
	}

	Quaternion& operator/=(const Quaternion& q){
		return multiply(q, false, true);
	}
	
	Quaternion operator*(const Quaternion& q) const {
		Quaternion p(*this);
		return p.multiply(q, false, false);
	}

	Quaternion operator%(const Quaternion& q) const {
		Quaternion p(*this);
		return p.multiply(q, true, false);
	}
	
	Quaternion operator/(const Quaternion& q) const {
		Quaternion p(*this);
		return p.multiply(q, false, true);
	}
	
	Quaternion& operator*=(double s){
		w*=s; x*=s; y*=s; z*=s;
		return *this;
	}
	
	double norm(){
		return sqrt(w*w + x*x + y*y + z*z);
	}
	
	Quaternion& normalize(){
		return (*this)*=1/norm();
	}
	
	void toScaledAxis(double* res){
		double nv = sqrt(x*x + y*y + z*z);
		if( nv < 1e-12) nv = 1e-12;
		// BUGFIX 2009-11-30: q and -q are represent the same rotation 
		//                    and must lead to the same result!
		// Note that singularity for w==0 is not dramatic, as 
		// atan(+/-inf) = +/- pi/2 (in this case both solutions are acceptable)
		double s = 2*atan(nv/w)/nv; // != 2*atan2(nv, w)/nv;
		*res++ = x*s; *res++ = y*s; *res++ = z*s;
	}
	
	/** 
	 * rotates the vector vec and stores the result in res
	 * vec and res can point to the same or overlapping addresses
	 * If back==true, it rotates backwards.
	 */
	void rotate(double res[3], const double vec[3], bool back=false) const {
		double 
			//ww=w*w, 
			wx=w*x, xx=x*x,
			wy=w*y, xy=x*y, yy=y*y,
			wz=w*z, xz=x*z, yz=y*z, zz=z*z;
		
		if(back){
			wx=-wx; wy=-wy; wz=-wz;
		}
		double X= *vec++, Y=*vec++, Z=*vec++;
		*res++ = (1-2*zz-2*yy)*X + 2*(xy-wz)*Y + 2*(wy+xz)*Z;
		*res++ = 2*(xy+wz)*X + (1-2*zz-2*xx)*Y + 2*(yz-wx)*Z;
		*res++ = 2*(xz-wy)*X + 2*(yz+wx)*Y + (1-2*yy-2*xx)*Z;
	}
};

#endif /*QUATERNION_H_*/
