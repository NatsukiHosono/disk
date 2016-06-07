#pragma once

namespace PARAM{
	const double T_end = 10.0;
	const short int Ndim = 2;
	const double OUTPUT_TIME_INTERVAL = T_end / 100.0;
}
template <typename Ptcl> void SetupIC(vector<Ptcl>& ptcl, boundary_t<double>& boundary){
	boundary.setXRange(0.0, 1.0);
	boundary.setYRange(0.0, 1.0);

	const int N = 64;
	const double dx = 1./(double)(N);
	int cnt = 0;
	for(double x = 0 ; x < 1.0 ; x += dx){
		for(double y = 0 ; y < 1.0 ; y += dx){
			++ cnt;
		}
	}
	for(double x = 0 ; x < 1.0 ; x += dx){
		for(double y = 0 ; y < 1.0 ; y += dx){
			const vec3<double> pos = vec3<double>(x, y, 0.0);
			const double dens = 1.0;
			const double mass = 1.0 / (double)(cnt);
			vec3<double> vel;
			if(y < 0.8){
				vel.x = 0.0;
			}else{
				vel.x = 0.01;
			}
			vel.y = 0;
			const double eng = 1.0 /(dens);
			Ptcl ith;
			ith.createParticle(pos, vel, mass, dens, eng);
			if(y < 0.2 || y > 0.8){
				ith.type = FREEZE;
			}
			ith.setEoS(&Diatomic);
			ptcl.push_back(ith);
		}
	}
	return;
}

template <typename real> void ieSPH<real>::setSource(){

}

