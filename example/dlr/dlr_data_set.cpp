#include <Estimator.h>
#include <manifolds/SOn.h>
#include <manifolds/Vect.h>
#include <types/RandomVariable.h>
#include <types/Measurement.h>
#include <tools/AutoConstruct.h>
#include <tools/MakePose.h>
#include <tools/CholeskyCovariance.h>

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
#include <sys/time.h>


using namespace std;
using namespace SLOM;



/**
 * 2D SLAM on DLR data set
 * 
 */

// Landmark is a very simply random var, so BUILD_RANDOMVAR isn't necessary
//BUILD_RANDOMVAR(LandMark, ((Vect<2>, pos)))
typedef RVWrapper<Vect<2> > LandMark;


MAKE_POSE2D(Pose, pos , orientation, )


/// Measurements:

BUILD_MEASUREMENT(Odo,3, ((Pose, t0)) ((Pose, t1)), ((Pose_T, odo)) ((CholeskyCovariance<3>, cov)) )
double* Odo::eval(double ret[3]) const
{
	Pose_T diff = t0->world2Local(*t1);
	diff.sub(ret, odo); 
	cov.invApply(ret);
	return ret+3; 
}


BUILD_MEASUREMENT(LM_observation, 2, ((Pose, pose)) ((LandMark, lm)), ((Vect<2>, rel_coord )) ((CholeskyCovariance<2>, cov)) )
double* LM_observation::eval(double ret[2]) const
{
	Vect<2> landmark = pose->world2Local(*lm);
	rel_coord.sub(ret, landmark);
	cov.invApply(ret);
	return ret+2;
}






#include "dlr_parser.h"



typedef map<int, LandMark*> LM_storage;
typedef LM_storage::value_type LM_storage_item;

void write_lm(ofstream& f, LM_storage_item item)
{
	f << (**item.second)[0] << " " << (**item.second)[1] << " " << item.first << endl;
}

void write_pos(ofstream& f, LM_storage_item item)
{
}

void dalloc(LM_storage_item item)
{
	delete item.second;
}

string make_filename( const string& basename, int index, const string& ext )
{
	ostringstream result;
	result << basename << setfill ('0') << setw (3) << index << ext;
	return result.str();
}

void parseLandmarks(LM_storage& landmarks, const char* filename){
	ifstream logfile(filename);
	string line;
	
	while(getline(logfile, line)){
		stringstream str(line);
		int id;
		double pos[2];
		str >> pos[0] >> pos[1] >> id;
		if(id >= 0){
			delete landmarks[id];
			landmarks[id] = new LandMark(pos, false);
		}
	}
}


int main(int argc, char* argv[])
{
	DLR_Data_Parser parser(argv[1]);
	int n_step = parser.n_step();
	parser.reset();
	int n_lm = parser.n_lm();
	parser.reset();
	
	Estimator e(Estimator::Cholesky, Estimator::GaussNewton);
	
	// State
	vector<Pose> poses; poses.reserve(n_step+1);
	LM_storage landmarks;
	
	// Observation
	vector<Odo> odometry; odometry.reserve( n_step );
	vector<LM_observation> lm_obs;  lm_obs.reserve( n_lm );
	
	
	if(argc>=3){ // read Landmarks from file
		parseLandmarks(landmarks, argv[2]);
		for(LM_storage::iterator it = landmarks.begin(); it != landmarks.end(); it++){
			e.insertRV(it->second);
		}
		poses.push_back(Pose(Pose_T(), true));
	} else {
		poses.push_back(Pose(Pose_T(), false));
	}
	// Start pose
	e.insertRV( &poses.back());
	
	vector<double> step;
	vector<double> lmob;
	// read data from data set initialize poses
	parser.reset();
	while ( parser.get_next_step(step) )
	{
		//step.clear();
		lmob.clear();
		Pose_T tmp_pose((Vect<2>(&step[0])), SO2(step[2]));
		Pose_T next_pos = poses.back()->local2World( tmp_pose );
		poses.push_back( next_pos );
		odometry.push_back( Odo( *(poses.end()-2), *(poses.end()-1),
				tmp_pose, 
				CholeskyCovariance<3>(&step[3], CholeskyMode::CHOLESKY_UPPER) )
		);
		e.insertRV( &poses.back());
		e.insertMeasurement( &odometry.back() );
		cout << "added STEP" << endl;
		while ( parser.get_next_landmark(lmob) )
		{
			Vect<2> t_vec2( &lmob[0] );
			// insert landmark for every ID
			vector<double>::iterator it= lmob.begin();
			//copy ( lmob.begin(), lmob.end(), ostream_iterator<double>(cout, " ") );
			it += 6;
			Vect<2> g_pos = next_pos.local2World( t_vec2 );
			while ( it != lmob.end() )
			{
				int id = static_cast<int>(*it);
				cout << "Landmark ID: " << id << endl;
				if ( id != -1 )
				{
					if (landmarks.find(id) == landmarks.end() )
					{
						LandMark* n = new LandMark(g_pos);
						landmarks.insert(make_pair( id, n ) );
						e.insertRV( landmarks[id]);
					}
					lm_obs.push_back( LM_observation(poses.back(), *landmarks[id], t_vec2, 
							CholeskyCovariance<2>( &lmob[3], CholeskyMode::CHOLESKY_UPPER ) ));
					e.insertMeasurement( &lm_obs.back() );
				}
				++it;
			}
			cout << "Added Landmark" << endl;
		}
	}
	cout << "read data " << endl;
	
	cout << "Poses: " << poses.size() << "\nLandmarks: " << landmarks.size() << "\nLM-Observations: " << lm_obs.size() << endl;
	
	
	
	
	ofstream outp("output.initpos");
	for (vector<Pose>::iterator it = poses.begin(); it != poses.end(); ++it)
	{
		outp << (*it)->pos[0] << " " << (*it)->pos[1] << endl;
	}
	outp.close();
	
	ofstream outi("output.init");
	for_each( landmarks.begin(), landmarks.end(), bind1st(ptr_fun(write_lm), outi) );
	outi.close();
	
	struct timeval ts, te;
	gettimeofday(&ts,0);
	e.initialize();
	std::cout << "initialisation done" << std::endl;
	int steps_arg = 20;
	for(int k=0; k< steps_arg; k++){
		std::cout << "\nStep " << k << std::endl;
		
		
		// here the actual optimization happens:
		double gain = e.optimizeStep();
		
		// rest is for storing intermediate results:
		ofstream out(make_filename("output",k,".lm").c_str());
		for_each( landmarks.begin(), landmarks.end(), bind1st(ptr_fun(write_lm), out) );
		out.close();
		ofstream outpos(make_filename("output",k,".pos").c_str());
		for (vector<Pose>::iterator it = poses.begin(); it != poses.end(); ++it)
		{
			outpos << (*it)->pos[0] << " " << (*it)->pos[1] << endl;
		}
		outpos.close();
		if( 0 <= gain && gain < 1e-9) break;
	}
	
	gettimeofday(&te,0);
	cout << "**** Optimization Done ****" << endl;
	
	double dts=(te.tv_sec-ts.tv_sec)+1e-6*(te.tv_usec-ts.tv_usec);
	cout << "TOTAL TIME= " << dts << " s." << endl;
	for_each( landmarks.begin(), landmarks.end(), dalloc);
	std::cout << "\n\nDone\n";
}
