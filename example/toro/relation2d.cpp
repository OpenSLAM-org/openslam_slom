#include <Estimator.h>
#include <tools/MakePose.h>
#include <tools/CholeskyCovariance.h>

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

MAKE_POSE2D(Pose, pos , orientation, )

typedef boost::ptr_map<int, Pose> Poses;

/// Measurement model:

BUILD_MEASUREMENT(Odo, 3, ((Pose, t0)) ((Pose, t1)), 
		((Pose_T, odo)) ((CholeskyCovariance<3>, cov)) )
double* Odo::eval(double ret[3]) const
{
	Pose_T diff = t0->world2Local(*t1);
	diff.sub(ret, odo); 
	cov.apply(ret);
	return ret+3;
}



void outputPoses(const Poses &poses, int k){
	std::ofstream out(make_filename("output",k,".pos").c_str());
	for (Poses::const_iterator it = poses.begin(); it != poses.end(); ++it) {
		out << (*it)->pos[0] << " " << (*it)->pos[1] << " " 
		    << (*it)->orientation << endl;
	}
	out.close();
}

int main(int argc, char** argv){
	if(argc < 2){
		cerr << "need input file\n";
		return -1;
	}
	ifstream logfile(argv[1]);
	
	
	Estimator e(Estimator::Cholesky, Estimator::GaussNewton);
	Poses poses;
	deque<Odo> odo;
	
	poses[0] = Pose(Pose_T(), false);
	
	e.insertRV(&poses[0]);
	
	std::string line;
	
	
	while(getline(logfile, line)){
		enum Indexes { RGX, RGY, RGPHI};

		istringstream inp(line);
		string tag;
		inp >> tag;
		if(tag=="VERTEX" || tag=="VERTEX2") {
			int id;
			double xEst[3];
			inp >> id >> xEst[RGX] >> xEst[RGY] >> xEst[RGPHI];
			if(poses.find(id)==poses.end()){
//				cout << id << " is new ";
				poses.insert(id, new Pose(Pose_T(xEst, xEst[RGPHI])));
				e.insertRV(&poses[id]);
			}
		} else if (tag == "EDGE" || tag == "EDGE2") {
			int frameA, frameB;
			double xEst[3];
			double xCov[3][3];
			//EDGE2 observed_vertex_id observing_vertex_id forward sideward rotate inf_ff inf_fs inf_ss inf_rr inf_fr inf_sr 
			inp >> frameA >> frameB
			    >> xEst[RGX] >> xEst[RGY] >> xEst[RGPHI]
			    >> xCov[RGX][RGX] >> xCov[RGX][RGY] >> xCov[RGY][RGY]
			    >> xCov[RGPHI][RGPHI] >> xCov[RGX][RGPHI] >> xCov[RGY][RGPHI];

			Pose_T delta(xEst, xEst[RGPHI]);
//			cout << "EDGE " << frameB << " " << frameA << ": ";
//			cout << frameB << "->" << frameA << "   ";
			if(poses.find(frameA)==poses.end()){
				if(poses.find(frameB)==poses.end()){
					poses[frameB] = Pose();
					e.insertRV(&poses[frameB]);
				}
				cout << frameA << " is new";
				poses[frameA]= Pose(poses[frameB]->local2World(delta));
				//assert(poses.size()==(unsigned)frameA+1);
				e.insertRV(&poses[frameA]);
			}
//			cout << endl;
			odo.push_back(Odo(poses[frameA], poses[frameB], delta, 
					CholeskyCovariance<3>(&xCov[0][0],CholeskyMode::CHOLESKY_FULL)));
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
		
		double gain = e.optimizeStep();
		
		//outputPoses(poses, k);
		if( 0 <= gain && gain < 1e-9) break;
	}
	gettimeofday(&te,0);
	cout << "**** Optimization Done ****" << endl;
	
	double dts=(te.tv_sec-ts.tv_sec)+1e-6*(te.tv_usec-ts.tv_usec);
	cout << "TOTAL TIME= " << dts << " s." << endl;
	
	
	cout << "Done";
}
