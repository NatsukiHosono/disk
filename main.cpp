#include "header.h"
/////////////////////
//Main
/////////////////////
int main(int argc, char* argv[]){
	#ifdef _OPENMP
	omp_set_num_threads(omp_get_max_threads()-2);
	#endif
	typedef ieSPH<double> ptcl_t;
	//alloc
	vector<ptcl_t> ptcl;
	Tree<ptcl_t> tree;
	system_t<double> sys;
	std::map<int, EoS::EoS_t<double>*> EoS;
	boundary_t<double> boundary;
	const kernel_t<double> kernel(PARAM::Ndim);
	FileIO<ptcl_t> IO(PARAM::OUTPUT_TIME_INTERVAL);

	//is NOT Restart
	if(argc == 1){
		//setup.
		SetupIC<ptcl_t>(ptcl, boundary);
		//Initialize
		#pragma omp parallel for //schedule(static)
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].id = i;
			ptcl[i].initialize(PARAM::Ndim);
		}
		//Prepare for the initial kick of time integration loop.
		tree.Plant(ptcl, kernel);
		tree.setCellProperties(kernel);
		tree.setNeighbourList(ptcl, kernel, boundary);
		for(std::size_t loop = 0 ; loop < PARAM::LOOP ; ++ loop){
			#pragma omp parallel for //schedule(static)
			for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
				ptcl[i].setDensity(kernel, boundary);
			}
			#pragma omp parallel for //schedule(static)
			for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
				ptcl[i].setSmoothingLength(PARAM::Ndim);
			}
			tree.setCellProperties(kernel);
			tree.setNeighbourList(ptcl, kernel, boundary);
		}
		#pragma omp parallel for //schedule(static)
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].setPressure();
			ptcl[i].setInverseOfMomentumMatrix(kernel, boundary, PARAM::Ndim);
			ptcl[i].setDivRotVelocity(kernel, boundary, PARAM::Ndim);
		}
		#pragma omp parallel for //schedule(static)
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].setTimeDerivative(kernel, boundary);
			ptcl[i].setSource();
			ptcl[i].setTimestep();
		}
		sys.setTimeStep(ptcl);
		//Output IC
		IO.OutputFile(ptcl, sys, boundary);
		IO.OutputBinary(ptcl, sys, boundary);

		sys.conservativeCheck(ptcl);
		printf("###############\n");
	}else{//is restart
		std::cout << "restart" << std::endl;
		IO.InputBinary(ptcl, &sys, &boundary);
		std::cout << "read done." << std::endl;
	}
	//start time integration loooooooooop;
	while(sys.time < PARAM::T_end){
		CPUtimer_t CPUtimer;
		sys.setTimeStep(ptcl);

		#pragma omp parallel for //schedule(static)
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].initialKick(sys.dt);
			ptcl[i].fullDrift(sys.dt, boundary);
			ptcl[i].predict(sys.dt);
		}
		++ sys.step;
		sys.time += sys.dt;

		removeParticle(ptcl);
		tree.Plant(ptcl, kernel);
		tree.setCellProperties(kernel);
		tree.setNeighbourList(ptcl, kernel, boundary);
		for(std::size_t loop = 0 ; loop < PARAM::LOOP ; ++ loop){
			#pragma omp parallel for //schedule(static)
			for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
				ptcl[i].setDensity(kernel, boundary);
			}
			#pragma omp parallel for //schedule(static)
			for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
				ptcl[i].setSmoothingLength(PARAM::Ndim);
			}
			tree.setCellProperties(kernel);
			tree.setNeighbourList(ptcl, kernel, boundary);
		}
		#pragma omp parallel for //schedule(static)
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].setPressure();
			ptcl[i].setInverseOfMomentumMatrix(kernel, boundary, PARAM::Ndim);
			ptcl[i].setDivRotVelocity(kernel, boundary, PARAM::Ndim);
		}
		#pragma omp parallel for //schedule(static)
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].setTimeDerivative(kernel, boundary);
			ptcl[i].setSource();
			ptcl[i].setTimestep();
		}
		IO.OutputFile(ptcl, sys, boundary);

		#pragma omp parallel for //schedule(static)
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].finalKick(sys.dt);
		}

		printf("step # %lu (dt = %.16e) progress rate: %lf\n", sys.step, sys.dt, sys.time / PARAM::T_end);
		sys.conservativeCheck(ptcl);
		printf("###############\n");
	}
	IO.OutputBinary(ptcl, sys, boundary);
	return 0;
}

