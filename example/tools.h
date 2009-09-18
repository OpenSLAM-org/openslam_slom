#ifndef TOOLS_H_
#define TOOLS_H_

#include <string>
#include <sstream>
#include <iomanip>

#include <manifolds/Vect.h>

// calculate res-=matrix*(vec-off)*factor
template<int n>
void matrixSubMulDiff(double res[n],
		const SLOM::Vect<n*n> &matrix, const SLOM::Vect<n> &vec,
		const SLOM::Vect<n> &off = SLOM::Vect<n>(), double factor=1){
	const double *m=matrix.data;
	for(int i=0; i<n; i++){
		double x=factor*(vec[i] - off[i]);
		for(int j=0; j<n; j++){
			res[j] -= x * (*m++);
		}
	}
}




inline std::string make_filename( const std::string& basename, int index,
		const std::string& ext ) {
	std::ostringstream result;
	result << basename << std::setfill('0') << std::setw(3) << index << ext;
	return result.str();
}


struct nullstream : std::ostream {
	struct nullbuf: std::streambuf {
		int overflow(int c) { return traits_type::not_eof(c); }
	} m_sbuf;
	nullstream(): std::ios(&m_sbuf), std::ostream(&m_sbuf) {}
};




#endif /* TOOLS_H_ */
