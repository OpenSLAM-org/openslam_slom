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
#include <exception>
#include <map>
#include <functional>

class DLR_Data_Parser
{
public:
	DLR_Data_Parser(const char* filename): fname(filename)
	{
		open(); 
	}
	
	void open()
	{
		f.open(fname.c_str());
		cout << "opening " <<  fname << "...";
		if ( f.good() )
		{
			cout << "okay" << endl;
			cout << "Position: " << f.tellg() << endl;
		}
		else
		{
			cout << "failed" << endl;
		}
	}
	
	int n_step()
	{
		size_t pos = f.tellg();
		f.seekg(ios_base::beg);
		string line;
		int count = 0;
		while ( getline(f, line) )
		{
			if ( string( &line[0], &line[4])  == string("STEP")) ++count;
		}
		f.seekg(pos);
		return count;
	}
	
	int n_lm()
	{
		size_t pos = f.tellg();
		f.seekg(ios_base::beg);
		string line;
		int count = 0;
		while ( getline(f, line) )
		{
			if ( string( &line[0], &line[4])  == string("LAND")) ++count;
		}
		f.seekg(pos);
		return count;
	}
	
	int max_lm_id()
	{
		size_t pos = f.tellg();
		f.seekg(ios_base::beg);
		string line;
		int count = 0;
		while (getline(f, line) )
		{
			if ( string( &line[0], &line[4])  == string("LAND"))
			{
				stringstream str(line);
				string token;
				double x,y, q, cov_x, cov_y, cov_xy;
				str >> token >> x >> y >> q >> cov_x >> cov_y >> cov_xy;
				while( str.good() )
				{
					int id;
					str >> id;
					count = max(count, id);
				}
			}
		}
		f.seekg(pos);
		return count;
	}
	
	void reset()
	{
		f.close();
		f.open(fname.c_str());
	}
	
	bool eof()
	{
		return f.eof();
	}
	
	bool get_next_step(vector<double>& step)
	{
		step.clear();
		//cout << "get_next_step()" << endl;
		//cout << "position: " << f.tellg() << endl;
		if ( f.tellg() == -1) return false; // uiuiui alle gehackt. mal aufräumen hier....
		do
		{
			//cout << "Reading line: (Should be step) " << line << "\n" << endl;
			if ( string(&line[0], &line[4]) == string("LAND") ) return false;
			if ( string(&line[0], &line[4]) == string("STEP") )
			{
				stringstream str(line);
				string token, file;
				double dx,dy, dphi, cov_xx, cov_xy, cov_yy, cov_xphi, cov_yphi, cov_phiphi;
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
				//copy(step.begin(), step.end(), ostream_iterator<double>(cout, " "));
				//cout << endl;
				getline(f, line);
				return true;
			}
		} while ( getline(f, line) );
		return false;
	}
	
	bool get_next_landmark(vector<double>& lm)
	{
		lm.clear();
		do
		{
			if ( line[0] == '#' ) continue;
			if ( string(&line[0], &line[4]) == string("STEP") ) return false;
			if ( string(&line[0], &line[4]) == string("LAND") )
			{
				stringstream str(line);
				string token;
				double x,y,q, cov_x, cov_y, cov_xy;
				str >> token >> x >> y >> q >> cov_x >> cov_xy >> cov_y;
				lm.push_back(x);
				lm.push_back(y);
				lm.push_back(q);
				lm.push_back(cov_x);
				lm.push_back(cov_xy);
				lm.push_back(cov_y);
				str.exceptions( ios::failbit| ios::badbit );
				try
				{
					for (;;)
					{
						int id;
						str >> id;
						lm.push_back(id);
					}
				}
				catch ( std::exception& e)
				{
					cout << "no more ids" << endl;
					getline(f, line);
					return true;
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
};


#endif /*DLR_PARSER_H_*/
