#ifndef DLR_TYPES_H_
#define DLR_TYPES_H_

#include <deque>
#include <boost/ptr_container/ptr_map.hpp>


#include <types/Measurement.h>
#include <tools/MakePose.h>
#include <tools/CholeskyCovariance.h>

#include "../tools.h"

// Landmark is a very simply random var, so BUILD_RANDOMVAR isn't necessary
//BUILD_RANDOMVAR(LandMark, ((Vect<2>, pos)))
typedef SLOM::RVWrapper<SLOM::Vect<2> > LandMark;


MAKE_POSE2D(Pose, pos , orientation, )


/// Measurements:

BUILD_MEASUREMENT(Odo, 3, ((Pose, t0)) ((Pose, t1)),
                  ((Pose_T, odo)) ((SLOM::CholeskyCovariance<3>, cov)) )
double* Odo::eval(double ret[3]) const
{
	Pose_T diff = t0->world2Local(*t1);
	diff.sub(ret, odo);
	cov.invApply(ret);
	return ret+3;
}


BUILD_MEASUREMENT(LM_observation, 2, ((Pose, pose)) ((LandMark, lm)),
		((SLOM::Vect<2>, rel_coord )) ((SLOM::CholeskyCovariance<2>, cov)) )
double* LM_observation::eval(double ret[2]) const
{
	SLOM::Vect<2> landmark = pose->world2Local(*lm);
	rel_coord.sub(ret, landmark);
	cov.invApply(ret);
	return ret+2;
}

// Optional: A calibration matrix, for the Landmark measurements
BUILD_RANDOMVAR(Calibration, ((SLOM::Vect<4>, mat)) ((SLOM::Vect<2>, off)));

BUILD_MEASUREMENT(LM_observation_Calib, 2,
		((Pose, pose)) ((LandMark, lm)) ((Calibration, cal)),
		((SLOM::Vect<2>, rel_coord )) ((SLOM::CholeskyCovariance<2>, cov)) )
double* LM_observation_Calib::eval(double ret[2]) const
{
	SLOM::Vect<2> landmark = pose->world2Local(*lm);
	SLOM::Vect<2> calibrated;
	matrixSubMulDiff(calibrated.data, cal->mat, cal->off, rel_coord);
	calibrated.sub(ret, landmark);
	cov.invApply(ret);
	return ret+2;
}



typedef boost::ptr_map<int, LandMark> LM_storage;
typedef std::deque<Pose> Poses;



#endif /* DLR_TYPES_H_ */
