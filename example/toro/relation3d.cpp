#include <Estimator.h>
#include <tools/MakePose.h>
#include <tools/CholeskyCovariance.h>

#include <algorithm>
#include <deque>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <boost/ptr_container/ptr_map.hpp>
#include <sys/time.h>

#include "../tools.h"

using namespace SLOM;


using namespace std;


/**
 * Reading and optimizing TORO2D .graph files.
 */

MAKE_POSE3D(Pose, pos , orientation, )



// This method is based on
// Quaternion<T>::Quaternion(const T phi, const T theta, const T psi)
// from transformation3.hxx from the toro-Framework.
void euler2quaternion(
		Quaternion &q,
		const double phi, const double theta, const double psi){
	double sphi = sin(phi), stheta = sin(theta), spsi = sin(psi);
	double cphi = cos(phi), ctheta = cos(theta), cpsi = cos(psi);
	
	double _r[3][3] = { //create rotational Matrix
		{cpsi*ctheta, cpsi*stheta*sphi - spsi*cphi, cpsi*stheta*cphi + spsi*sphi},
		{spsi*ctheta, spsi*stheta*sphi + cpsi*cphi, spsi*stheta*cphi - cpsi*sphi},
		{    -stheta,                  ctheta*sphi,                  ctheta*cphi}
	};
	
	q.w = sqrt(std::max(0.0, 1 + _r[0][0] + _r[1][1] + _r[2][2]))/2.0;
	q.x = sqrt(std::max(0.0, 1 + _r[0][0] - _r[1][1] - _r[2][2]))/2.0;
	q.y = sqrt(std::max(0.0, 1 - _r[0][0] + _r[1][1] - _r[2][2]))/2.0;
	q.z = sqrt(std::max(0.0, 1 - _r[0][0] - _r[1][1] + _r[2][2]))/2.0;
	
	if(_r[2][1] - _r[1][2]<0) q.x = -q.x;
	if(_r[0][2] - _r[2][0]<0) q.y = -q.y;
	if(_r[1][0] - _r[0][1]<0) q.z = -q.z;
}


void readPose(Pose_T & pose, istream & inp){
	for(int i=0; i<3; i++) inp >> pose.pos[i];
	double phi, theta, psi;
	inp >> phi >> theta >> psi;
	
	euler2quaternion(pose.orientation.quat, phi, theta, psi);
}



typedef boost::ptr_map<int, Pose> Poses;

/// Measurement model:

BUILD_MEASUREMENT(Odo, 6, ((Pose, t0)) ((Pose, t1)), 
		((Pose_T, odo)) /* ((CholeskyCovariance<3>, cov)) */)
double* Odo::eval(double ret[6]) const
{
	Pose_T diff = t0->world2Local(*t1);
	odo.sub(ret, diff); 
	//cov.apply(ret);
	return ret+6;
}



void outputPoses(const Poses &poses, int k){
	std::ofstream out(make_filename("output",k,".pos").c_str());
	for (Poses::const_iterator it = poses.begin(); it != poses.end(); ++it) {
		out << (*it)->pos[0] << " " << (*it)->pos[1] << " " << (*it)->pos[2] << " ";
		const Quaternion &q = (*it)->orientation.quat;
		out << q.w << " " << q.x << " " << q.y << " " << q.z << endl;
	}
	out.close();
}

int main(int argc, char** argv){
	if(argc < 2){
		cerr << "need input file\n";
		return -1;
	}
	ifstream logfile(argv[1]);
	
	
	Estimator e(Estimator::Cholesky, Estimator::Levenberg, 10);
	Poses poses;
	deque<Odo> odo;
		
	std::string line;
	
	
	bool firstPose = true;
	while(getline(logfile, line)){
		enum Indexes { RGX, RGY, RGPHI};

		istringstream inp(line);
		string tag;
		inp >> tag;
		if(tag=="VERTEX" || tag=="VERTEX3") {
			//if(!firstPose) continue;
			int id;
			inp >> id;
			//check if pose already included:
			if(poses.find(id)!=poses.end()) continue;
			Pose_T pose;
			readPose(pose, inp);
			poses.insert(id, new Pose(pose, !firstPose));  // TODO poses[id] = Pose(pose);
			e.insertRV(&poses[id]);
			firstPose = false;

		} else if (tag == "EDGE" || tag == "EDGE3") {
			int frameA, frameB;
			//EDGE3 observed_vertex_id observing_vertex_id forward sideward rotate inf_ff inf_fs inf_ss inf_rr inf_fr inf_sr 
			inp >> frameA >> frameB;
			
			Pose_T delta;
			readPose(delta, inp);

			if(poses.find(frameA)==poses.end()){
				if(poses.find(frameB)==poses.end()){
					cerr << "PROBLEM: Edge between unknown Vertices " 
					     << frameA  << " and " << frameB << endl;
					poses[frameB] = Pose();
					e.insertRV(&poses[frameB]);
				}
				cout << frameA << " is new";
				poses[frameA]= Pose(poses[frameB]->local2World(delta));
				//assert(poses.size()==(unsigned)frameA+1);
				e.insertRV(&poses[frameA]);
			}
			odo.push_back(Odo(poses[frameA], poses[frameB], delta));
			e.insertMeasurement(&odo.back());
		}
	}
	logfile.close();
	outputPoses(poses, 0);
	
	cout << "\nLogfile read\nInitializing" << endl;
	
	struct timeval ts, te;
	gettimeofday(&ts,0);
	
	e.initialize();

	int kMax = 50; //TODO read from commandline
	for(int k=1; k<=kMax; k++){
		cout << "Step " << k << ": ";
		double gain = e.optimizeStep();
		
		outputPoses(poses, k);
		if( 0 <= gain && gain < 1e-9) break;
		if(k==15){
			cout << "\n Switching to Gauss-Newton\n";
			e.changeAlgorithm(SLOM::Estimator::GaussNewton);
			
		}
	}
	gettimeofday(&te,0);
	cout << "**** Optimization Done ****" << endl;
	
	double dts=(te.tv_sec-ts.tv_sec)+1e-6*(te.tv_usec-ts.tv_usec);
	cout << "TOTAL TIME= " << dts << " s." << endl;
	
	cout << "Done";
}
