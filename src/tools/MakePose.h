#ifndef MAKEPOSE_H_
#define MAKEPOSE_H_

#include "AutoConstruct.h"
#include "../manifolds/SOn.h"
#include "../manifolds/Vect.h"
#include "../types/RandomVariable.h"

/**
 * The MAKE_POSE macro gets a baseclass and derivates a class, which 
 * provides local2World and world2Local methods, that convert Vect<>s and 
 * Poses between corresponding coordinate systems.
 * 
 */


#define MAKE_POSE(name, posT, posN, orientT, orientN, otherVars) \
BUILD_RANDOMVAR(name ## _Base, ((posT, posN)) ((orientT, orientN)) otherVars) \
struct name ## _T : public name ## _Base_T { \
	name ## _T(const posT &pos=posT(), const orientT& orientation=orientT()) : name ## _Base_T(pos, orientation) {} \
	name ## _T local2World(const name ## _T& oth) const { \
		name ## _T pose; \
		pose.orientN = orientN * oth.orientN; \
		pose.posN = orientN * oth.posN; \
		pose.posN.add(posN.data); \
		return pose; \
	} \
	posT local2World(const posT& oth) const { \
		posT position; \
		position = orientN * oth; \
		position.add(posN.data); \
		return position; \
	} \
	name ## _T world2Local(const name ## _T& oth) const { \
		name ## _T pose(oth); \
		pose.posN.add(posN.data, -1); \
		pose.orientN = orientN % oth.orientN; \
		pose.posN = orientN % pose.posN; \
		return pose; \
	} \
	posT world2Local(const posT& oth) const { \
		posT position(oth); \
		position.add(posN.data, -1); \
		position = orientN % position; \
		return position; \
	} \
}; \
typedef SLOM::RVWrapper<name ## _T> name;

#define MAKE_POSE2D(name, posN, orientN, otherVars) \
MAKE_POSE(name, SLOM::Vect<2>, posN, SLOM::SO2, orientN, otherVars)

#define MAKE_POSE3D(name, posN, orientN, otherVars) \
MAKE_POSE(name, SLOM::Vect<3>, posN, SLOM::SO3, orientN, otherVars)

#endif /*MAKEPOSE_H_*/
