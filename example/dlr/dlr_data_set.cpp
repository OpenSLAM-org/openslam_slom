/**
 * 2D SLAM on DLR data set
 *
 */


#include <Estimator.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

#include <sys/time.h>


#include "dlr_parser.h"
#include "dlr_types.h"

#include "../tools.h"

using namespace SLOM;
using namespace std;


//#define DLR_CALIBRATE


#ifndef DLR_CALIBRATE
typedef LM_observation LM_obs;
#else
typedef LM_observation_Calib LM_obs;
#endif

void outputLandmarks(const LM_storage &lms, int k){
	std::ofstream out(make_filename("output",k,".lm").c_str());
	for(LM_storage::const_iterator it=lms.begin(); it!=lms.end(); ++it){
		out << (**it)[0] << " " << (**it)[1] << " " << it.key() << endl;
	}
	out.close();
}

void outputPoses(const Poses &poses, int k){
	std::ofstream out(make_filename("output",k,".pos").c_str());
	for (Poses::const_iterator it = poses.begin(); it != poses.end(); ++it) {
		out << (*it)->pos[0] << " " << (*it)->pos[1] << endl;
	}
	out.close();
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
			landmarks[id] = LandMark(pos, false);
		}
	}
}

void parse(Poses &poses, LM_storage &landmarks,
		deque<Odo> &odometry, deque<LM_obs> &lm_obs,
		Estimator &e, DLR_Data_Parser &p,
#ifdef DLR_CALIBRATE
		Calibration& cal,
#endif
		std::ostream &out = std::cout){

	vector<double> step;
	// read data from data set initialize poses
	while ( p.get_next_step(step) ) {
		Pose_T tmp_pose((Vect<2>(&step[0])), SO2(step[2]));
		Pose_T next_pos = poses.back()->local2World( tmp_pose );
		poses.push_back( next_pos );
		odometry.push_back(
				Odo( *(poses.end()-2), *(poses.end()-1),
						tmp_pose,
						CholeskyCovariance<3>(&step[3], CholeskyMode::CHOLESKY_UPPER)
				) );
		e.insertRV( &poses.back());
		e.insertMeasurement( &odometry.back() );
		out << "added STEP" << endl;
		vector<double> lmob;
		while ( p.get_next_landmark(lmob) ) {
			Vect<2> t_vec2( &lmob[0] );
			// insert landmark for every ID
			vector<double>::iterator it= lmob.begin();
			it += 6;
			for( ; it != lmob.end(); ++it ) {
				int id = static_cast<int>(*it);
				out << "Landmark ID: " << id << endl;
				if ( id != -1 ) {
					if (landmarks.find(id) == landmarks.end() ) {
						landmarks[id] = next_pos.local2World( t_vec2 );
						e.insertRV( &landmarks[id]);
					}
					lm_obs.push_back(LM_obs(poses.back(), landmarks[id],
#ifdef DLR_CALIBRATE
							cal,
#endif
							t_vec2, CholeskyCovariance<2>( &lmob[3], CholeskyMode::CHOLESKY_UPPER ) )
					);
					e.insertMeasurement( &lm_obs.back() );
				}
			}
			out << "Added Landmark" << endl;
		}
	}
}



int main(int argc, char* argv[])
{

	Estimator e(Estimator::Cholesky, Estimator::GaussNewton);

	// State
	Poses poses;
	LM_storage landmarks;

	// Observation
	deque<Odo> odometry;
	deque<LM_obs> lm_obs;

#ifdef DLR_CALIBRATE
	double cal_init[4] = {1,0,0,1};
	Calibration calib = Calibration_T(cal_init);
	e.insertRV(&calib);
#endif

	nullstream nullstrm;

	DLR_Data_Parser parser(argv[1],  nullstrm);


	if(argc>=3){ // read Landmarks from file
		parseLandmarks(landmarks, argv[2]);
		for(LM_storage::iterator it = landmarks.begin(); it != landmarks.end(); it++){
			e.insertRV(&(*it));
		}
		poses.push_back(Pose(Pose_T(), true));
	} else {
		poses.push_back(Pose(Pose_T(), false));
	}
	// Start pose
	e.insertRV( &poses.back());

	parse(poses, landmarks, odometry, lm_obs, e, parser,
#ifdef DLR_CALIBRATE
			calib,
#endif
			nullstrm);

	cout << "read data " << endl;

	cout << "Poses: " << poses.size()
	<< "\nLandmarks: " << landmarks.size()
	<< "\nLM-Observations: " << lm_obs.size() << endl;


	outputPoses(poses, 0);
	outputLandmarks(landmarks, 0);

	struct timeval ts, te;
	gettimeofday(&ts,0);
	e.initialize();
	std::cout << "Initialization done" << std::endl;
	int steps_arg = 20;
	for(int k=1; k <= steps_arg; k++){
		std::cout << "\nStep " << k << std::endl;

		// here the actual optimization happens:
		double gain = e.optimizeStep();

		// rest is for storing intermediate results:
		outputLandmarks(landmarks, k);
		outputPoses(poses, k);

		if( 0 <= gain && gain < 1e-5) break;
	}

	gettimeofday(&te,0);
	cout << "**** Optimization Done ****" << endl;

	double dts=(te.tv_sec-ts.tv_sec)+1e-6*(te.tv_usec-ts.tv_usec);
	cout << "TOTAL TIME= " << dts << " s." << endl;
#ifdef DLR_CALIBRATE
	cout << "Calibration result:\n";
	for(int i=0; i<4; i++){
		cout << calib->mat[i] << "\t";
		if(i % 2 == 1) cout << "\t" << calib->off[i/2] << "\n";
	}
#endif

	std::cout << "\n\nDone\n";
}
