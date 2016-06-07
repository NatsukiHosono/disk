#pragma once

namespace PARAM{
	const double T_end = 0.2;
	const double OUTPUT_TIME_INTERVAL = 0.01;
	const short int Ndim = 1;
}

template <typename Ptcl> void SetupIC(vector<Ptcl>& ptcl, boundary_t<double>& boundary){
	boundary.setXRange(-1.0, 1.0);
	int N = 1024;
	double dx = 1./(double)(N);
	double dens_L = 1.0;
	double dens_R = 0.5;
	double vel_L = 0;
	double vel_R = 0;
	double pres_L = 1.0;
	double pres_R = 0.2;
	std::size_t cnt = 0;
	for(double x = 0 ; x < 0.5 ; x += dx){
		++cnt;
	}
	for(double x = 0.5 ; x < 1.0 ; x += dx * 2.0){
		++cnt;
	}
	std::cout << cnt << std::endl;
	for(double x = -1.0 ; x < 0.0 ; x += dx){
		vec3<double> pos = vec3<double>(x, 0.0, 0.0);
		double specific_int_eng = pres_L / (1.4 - 1.0) / dens_L;
		double mass = (dens_L + dens_R) * 0.5 / double(cnt);
		//SET
		//if(ith.r.x < 0.1 || ith.r.x > 0.9) ith.type = FREEZE;
		Ptcl ith;
		ith.createParticle(pos, vel_L, mass, dens_L, specific_int_eng);
		ith.setEoS(&Diatomic);
		ptcl.push_back(ith);
	}
	for(double x = 0.0 ; x < 1.0 ; x += dx * 2.0){
		vec3<double> pos = vec3<double>(x, 0.0, 0.0);
		double specific_int_eng = pres_R / (1.4 - 1.0) / dens_R;
		double mass = (dens_L + dens_R) * 0.5 / double(cnt);
		//SET
		//if(ith.r.x < 0.1 || ith.r.x > 0.9) ith.type = FREEZE;
		Ptcl ith;
		ith.createParticle(pos, vel_R, mass, dens_R, specific_int_eng);
		ith.setEoS(&Diatomic);
		ptcl.push_back(ith);
	}
	return;
}

template <typename real> void ieSPH<real>::setSource(){
}


