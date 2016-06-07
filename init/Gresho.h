#pragma once

namespace PARAM{
	const double T_end = 3.0;
	const short int Ndim = 2;
	const double OUTPUT_TIME_INTERVAL = T_end / 50.0;
}

template <typename Ptcl> void SetupIC(vector<Ptcl>& ptcl, boundary_t<double>& boundary){
	boundary.setXRange(-1.0, 1.0);
	boundary.setYRange(-1.0, 1.0);
	
	const int N = 64;
	const double dx = 1./(double)(N);
	const double hcr = 1.4;
	const double dens = 1.0;
	for(double x = -1.0 ; x < 1.0 ; x += dx){
		for(double y = -1.0 ; y < 1.0 ; y += dx){
			double pres;
			vec3<double> vel;
			const vec3<double> pos = vec3<double>(x, y, 0.0);
			const double R = sqrt(x * x + y * y);
			const double angle = atan2(y, x);
			if(R < 0.2){
				pres = 5.0 + 12.5 * R * R;
				vel = 5.0 * R * R * vec3<double>(- sin(angle), cos(angle), 0.0);
			}else if(R < 0.4){
				pres = 5.0 + 12.5 * R * R + 4.0 - 20.0 * R + 4.0 * log(5.0 * R);
				vel = (2.0 - 5.0 * R) * R * vec3<double>(- sin(angle), cos(angle), 0.0);
			}else{
				pres = 5.0 + 2.0 * (2.0 * log(2.0) - 1.0);
				vel = 0.0;
			}
			const double eng = pres / ((hcr - 1.0) * dens);
			const double mass = 1.0 / (N * N);

			Ptcl ith;
			ith.createParticle(pos, vel, mass, dens, eng);
			ith.setEoS(&Diatomic);
			ptcl.push_back(ith);
		}
	}
	return;
}

template <typename real> void ieSPH<real>::setSource(){

}

