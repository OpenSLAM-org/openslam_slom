#ifndef AUTOCONSTRUCT_H_
#define AUTOCONSTRUCT_H_

#include <boost/preprocessor/seq.hpp>


/**
 * Macros to automatically construct RandomVariables and Measurements.
 */





/**
 * BUILD_RANDOMVAR(name, entries)
 * makro to generate RandomVariables (derived from RandomVar, IRandomVar)
 * name is the class-name of the RV, entries is the list of entries like this:
 * BUILD_RANDOMVAR(Pose3D, (( Vect<3>, pos)) (( SO3, orientation)) )
 * Whitespace is optional, but the double parentheses are necessary.
 * Construction is done entirely in preprocessor.
 */
#define BUILD_RANDOMVAR(name, entries) \
struct name ## _T{ \
	SEQ_FOR_EACH(RANDOMVAR_GENERATE_VARLIST, entries) \
	enum {DOF = 0 \
	SEQ_FOR_EACH(RANDOMVAR_GENERATE_GETDOF, entries) \
	}; \
	name ## _T( \
		SEQ_TRANSFORM(RANDOMVAR_CONSTRUCTOR_ARG, entries) \
		) : \
		SEQ_TRANSFORM(RANDOMVAR_GENERATE_CONSTRUCTOR, entries) {}\
	int getDOF() const { return DOF; } \
	const double* add(const double* vec, double scale=1) { \
		SEQ_FOR_EACH(RANDOMVAR_GENERATE_ADD, entries) \
		return vec; \
	} \
	double* sub(double *res, const name ## _T& oth) const { \
		SEQ_FOR_EACH(RANDOMVAR_GENERATE_SUB, entries) \
		return res; \
	} \
}; \
typedef SLOM::RVWrapper<name ## _T> name;


/**
 * BUILD_MEASUREMENT(name, dim, variables, data)
 * makro to generate Measurements (derived from IMeasurement)
 * name is the class-name of the Measurement, dim is its dimension.
 * variables are the RVs the measurement depends on and data can be
 * some arbitrary data variables.
 * Both variables and data must be given in a list like this:
 * BUILD_MEASUREMENT(OdoMeas, 3,
 *    ((Pose, start)) ((Pose, end)),
 *    ((double, linear)) ((double, theta))
 * )
 * Whitespace is optional, but the double parentheses are necessary.
 * Construction is done entirely in preprocessor.
 * After declaration the function double* name::eval(double *ret) const
 * has to be defined
 */
#define BUILD_MEASUREMENT(name, dim, variables, data) \
struct name : public SLOM::IMeasurement { \
	SEQ_FOR_EACH(MEASUREMENT_GENERATE_VARLIST, variables) \
	SEQ_FOR_EACH(MEASUREMENT_GENERATE_DATALIST, data)     \
	name( \
		SEQ_TRANSFORM(MEASUREMENT_CONSTRUCTOR_ARG, variables) \
		SEQ_FOR_EACH(MEASUREMENT_CONSTRUCTOR_ARG_D, data) \
		) : \
		SEQ_TRANSFORM(MEASUREMENT_GENERATE_CONSTRUCTOR, variables) \
		SEQ_FOR_EACH(MEASUREMENT_GENERATE_CONSTRUCTOR_D, data) {}\
	int getDepend() const { return 0 \
		SEQ_FOR_EACH(MEASUREMENT_GENERATE_DEPEND, variables);} \
	int registerVariables() const{ return 0\
		SEQ_FOR_EACH(MEASUREMENT_GENERATE_REGISTER, variables);} \
	name& operator=(const name& f){ return *(new(this)name(f)); }\
	int getDim() const { return dim; } \
	double* eval(double *ret) const; \
};


#define SEQ_TRANSFORM(op, seq) BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM_S(1, op, , seq))

#define SEQ_FOR_EACH(macro, entries) BOOST_PP_SEQ_FOR_EACH_R(1, macro, , entries)



// Internal Macros:
#define TYPE_ID_OF_TYPEID(type, id) type id
#define TYPEREF_ID_OF_TYPEID(type, id) type& id
#define TYPE_OF_TYPEID(type, id) type
#define ID_OF_TYPEID(type, id) id


#define RANDOMVAR_GENERATE_VARLIST(r, data, type_id) \
	TYPE_ID_OF_TYPEID type_id;

#define RANDOMVAR_GENERATE_GETDOF(r, data, type_id) \
	+ TYPE_OF_TYPEID type_id ::DOF

#define RANDOMVAR_CONSTRUCTOR_ARG(r, data, type_id) \
	const TYPEREF_ID_OF_TYPEID type_id = TYPE_OF_TYPEID type_id()

#define RANDOMVAR_GENERATE_CONSTRUCTOR(r, data, type_id) \
	ID_OF_TYPEID type_id(ID_OF_TYPEID type_id)


#define RANDOMVAR_GENERATE_ADD(r, data, type_id) \
	vec = ID_OF_TYPEID type_id .add(vec, scale);

#define RANDOMVAR_GENERATE_SUB(r, data, type_id) \
	res = ID_OF_TYPEID type_id . sub(res, oth. ID_OF_TYPEID type_id);


#define MEASUREMENT_GENERATE_VARLIST(r, data, type_id) \
	TYPEREF_ID_OF_TYPEID type_id;

#define MEASUREMENT_GENERATE_DATALIST(r, data, type_id) \
	TYPE_ID_OF_TYPEID type_id;

#define MEASUREMENT_CONSTRUCTOR_ARG(r, data, type_id) \
	TYPEREF_ID_OF_TYPEID type_id

#define MEASUREMENT_CONSTRUCTOR_ARG_D(r, data, type_id) \
	, TYPE_ID_OF_TYPEID type_id

#define MEASUREMENT_GENERATE_CONSTRUCTOR(r, data, type_id) \
	ID_OF_TYPEID type_id(ID_OF_TYPEID type_id)

#define MEASUREMENT_GENERATE_CONSTRUCTOR_D(r, data, type_id) \
	, ID_OF_TYPEID type_id(ID_OF_TYPEID type_id)

#define MEASUREMENT_GENERATE_DEPEND(r, data, type_id) \
	+ TYPE_OF_TYPEID type_id ::DOF

#define MEASUREMENT_GENERATE_REGISTER(r, data, type_id) \
	+ ID_OF_TYPEID type_id .registerMeasurement(this)



#endif /*AUTOCONSTRUCT_H_*/
