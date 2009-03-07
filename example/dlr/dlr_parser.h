#ifndef DLR_PARSER_H_
#define DLR_PARSER_H_
#include <vector>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>

class DLR_Data_Parser {
public:
	typedef std::vector<double> dvector;

	DLR_Data_Parser(const char* filename, std::ostream &out = std::cout): fname(filename), out(out) {
		open();
	}

	void open() {
		f.open(fname.c_str());
		out << "opening " <<  fname << "...";
		if ( f.good() ) {
			out << "okay\n";
			out << "Position: " << f.tellg() << std::endl;
		} else {
			out << "failed" << std::endl;
		}
	}

	int n_step() { return count("STEP"); }

    int n_lm() { return count("LAND"); }

	int max_lm_id()
	{
		size_t pos = f.tellg();
		f.seekg(std::ios_base::beg);
		std::string line;
		int count = 0;
		while (getline(f, line) ) {
			if ( line.compare(0, 4, "LAND") == 0) {
				std::stringstream str(line);
				std::string token;
				double x,y, q, cov_x, cov_y, cov_xy;
				str >> token >> x >> y >> q >> cov_x >> cov_y >> cov_xy;
				while( str.good() ) {
					int id;
					str >> id;
					count = std::max(count, id);
				}
			}
		}
		f.seekg(pos);
		return count;
	}

	void reset() {
		f.close();
		f.open(fname.c_str());
	}

	bool eof() {
		return f.eof();
	}

	bool get_next_step(dvector& step) {
		step.clear();
		if ( f.tellg() == std::istream::pos_type(-1)) return false;
		do {
			if ( line.compare(0, 4, "LAND") == 0) return false;
			if ( line.compare(0, 4, "STEP") == 0) {
				std::stringstream str(line);
				std::string token, file;
				double dx, dy, dphi, cov_xx, cov_xy, cov_yy, cov_xphi, cov_yphi, cov_phiphi;
				str >> token >> file >> dx >> dy >> dphi >> cov_xx >> cov_xy >> cov_yy >> cov_xphi >> cov_yphi >> cov_phiphi;
				step.push_back(dx);
				step.push_back(dy);
				step.push_back(dphi);
				step.push_back(cov_xx);
				step.push_back(cov_xy);
				step.push_back(cov_xphi);
				step.push_back(cov_yy);
				step.push_back(cov_yphi);
				step.push_back(cov_phiphi);
				getline(f, line);
				return true;
			}
		} while ( getline(f, line) );
		return false;
	}

	bool get_next_landmark(dvector& lm) {
		lm.clear();
		do {
			if ( line[0] == '#' ) continue;
			if ( line.compare(0, 4, "STEP") == 0) return false;
			if ( line.compare(0, 4, "LAND") == 0) {
				std::stringstream str(line);
				std::string token;
				double x,y,q, cov_x, cov_y, cov_xy;
				str >> token >> x >> y >> q >> cov_x >> cov_xy >> cov_y;
				lm.push_back(x);
				lm.push_back(y);
				lm.push_back(q);
				lm.push_back(cov_x);
				lm.push_back(cov_xy);
				lm.push_back(cov_y);
				int id;
				while (str >> id, str.good() ) {
					lm.push_back(id);
				}
				getline(f, line);
				return true;
			}
		} while ( getline(f, line) );
		return false;
	}

protected:
	std::string fname;
	std::ifstream f;
	std::string line;
	std::ostream &out;

	int count(const char token[4])    {
		size_t pos = f.tellg();
		f.seekg(std::ios_base::beg);
		std::string line;
		int count = 0;
		while(getline(f, line)){
			if(line.compare(0, 4, token) == 0)
				++count;
		}
		f.seekg(pos);
		return count;
	}



};


#endif /*DLR_PARSER_H_*/
