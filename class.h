#pragma once

template <typename real> struct boundary_t{
	private:
	real x_max, x_min, y_max, y_min, z_max, z_min;
	real x_range, y_range, z_range;
	public:
	boundary_t() : x_max(1.0e+30), x_min(-1.0e+30), y_max(1.0e+30), y_min(-1.0e+30), z_max(1.0e+30), z_min(-1.0e+30){
	}
	void setXRange(const real _x_min, const real _x_max){
		x_min = _x_min;
		x_max = _x_max;
		x_range = x_max - x_min;
	}
	void setYRange(const real _y_min, const real _y_max){
		y_min = _y_min;
		y_max = _y_max;
		y_range = y_max - y_min;
	}
	void setZRange(const real _z_min, const real _z_max){
		z_min = _z_min;
		z_max = _z_max;
		z_range = z_max - z_min;
	}
	inline vec3<real> Periodic(vec3<real> r) const{
		r.x = (r.x > 0.5 * x_range) ? r.x - x_range : ((r.x < - 0.5 * x_range) ? r.x + x_range : r.x);
		r.y = (r.y > 0.5 * y_range) ? r.y - y_range : ((r.y < - 0.5 * y_range) ? r.y + y_range : r.y);
		r.z = (r.z > 0.5 * z_range) ? r.z - z_range : ((r.z < - 0.5 * z_range) ? r.z + z_range : r.z);
		return r;
	}
	inline vec3<real> PeriodicWrapping(vec3<real> r) const{
		if (r.x >= x_max) r.x -= x_range;
		if (r.x <  x_min) r.x += x_range;

		if (r.y >= y_max) r.y -= y_range;
		if (r.y <  y_min) r.y += y_range;

		if (r.z >= z_max) r.z -= z_range;
		if (r.z <  z_min) r.z += z_range;
		return r;
	}
};

enum PTCLType{
	HYDRO, //Hydro particle.
	FREEZE,//calculate NOTHING, keep the velocity and specific internal energy.
};

template <typename real> struct ieSPH{
	////////////
	protected:
	EoS::EoS_t<real>* EoS;//Equation of State;
	public:
	vec3<real> r, v, a, v_half;//position vectors;
	real eng, eng_half;
	real eng_dot;
	
	vec3<real> AVa, ext_a;
	real AVeng_dot;
	
	real one;        //this should be one;
	real h;          //smoothing length;
	real gradh;
	real pres;       //PRESsure;
	real snds;       //SouND Speed;
	real dens;       //DENSity
	real mass;       //mass
	real psmth;      //smoothed pres. for DISPH.
	real div_v;       //nabla  . v
	vec3<real> rot_v; //nabla  x v = vorticity
	mat33<real> dya_v;//nabla (x) v
	vec3<real> div_div_v; //nabla nabla . v
	vec3<real> lap_v; //laplacian v
	mat33<real> dya_a;//nabla (x) a
	real qAV;        //AV
	real aAV;        //strength of AV
	real sAV;        //Shear switch
	real isInShock;  //CD10
	mat33<real> Minv;//Moment matrix inverse
	real dt;         //timestep
	std::size_t id;  //index;
	PTCLType type;   //Type of Ptcl;
	std::vector<const ieSPH*> ngb;//neighbour list
	short tag;
	////////////
	ieSPH(): type(HYDRO){
	}
	////////////
	inline void initialize(const int Ndim){
		aAV = PARAM::AV_alpha;
		gradh = 1.0;
		one = 1.0;
		setPressure();
		setSmoothingLength(Ndim);
	}
	inline void setDensity(const kernel_t<real>& kernel, const boundary_t<real>& boundary){
		psmth = 0;
		dens = 0;
		#if __AVX__
		AVX dens_tmp = 0;
		AVX psmth_tmp = 0;
		const std::size_t volatile loopsize = (ngb.size() / 4) * 4;
		const AVX xi = r.x;
		const AVX yi = r.y;
		const AVX zi = r.z;
		for(std::size_t j = 0 ; j < loopsize ; j += 4){
			const AVX xj = (v4df){ngb[j+0]->r.x, ngb[j+1]->r.x, ngb[j+2]->r.x, ngb[j+3]->r.x};
			const AVX yj = (v4df){ngb[j+0]->r.y, ngb[j+1]->r.y, ngb[j+2]->r.y, ngb[j+3]->r.y};
			const AVX zj = (v4df){ngb[j+0]->r.z, ngb[j+1]->r.z, ngb[j+2]->r.z, ngb[j+3]->r.z};
			const AVX drx = xj - xi;
			const AVX dry = yj - yi;
			const AVX drz = zj - zi;
			const AVX dr = (drx * drx + dry * dry + drz * drz).sqrt();
			const AVX hi = h;
			const AVX mj   =      (v4df){ngb[j+0]->mass, ngb[j+1]->mass, ngb[j+2]->mass, ngb[j+3]->mass};
			const AVX Engj = mj * (v4df){ngb[j+0]->eng , ngb[j+1]->eng , ngb[j+2]->eng , ngb[j+3]->eng };
			dens_tmp += mj * kernel(dr, hi);
			psmth_tmp += Engj * kernel(dr, hi);
		}
		dens  = dens_tmp[0]  + dens_tmp[1]  + dens_tmp[2]  + dens_tmp[3];
		psmth = psmth_tmp[0] + psmth_tmp[1] + psmth_tmp[2] + psmth_tmp[3];
		for(std::size_t j = loopsize ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(ngb[j]->r - r);
			dens  += ngb[j]->mass * kernel(dr, h);
			psmth += ngb[j]->mass * ngb[j]->eng * kernel(dr, h);
		}
		#else
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(ngb[j]->r - r);
			dens  += ngb[j]->mass * kernel(dr, h);
			psmth += ngb[j]->mass * ngb[j]->eng * kernel(dr, h);
		}
		#endif
	}
	inline void testUnity(const kernel_t<real>& kernel, const boundary_t<real>& boundary){
		one = 0.0;
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(ngb[j]->r - r);
			one += ngb[j]->mass / ngb[j]->dens * kernel(dr, h);
		}
	}
	inline void setSmoothingLength(const int Ndim){
		h = PARAM::SMTH * math::rootdim(mass / dens, Ndim);
	}
	inline void setInverseOfMomentumMatrix(const kernel_t<real>& kernel, const boundary_t<real>& boundary, const short int Ndim){
		Minv = 0;
		gradh = 0;
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(ngb[j]->r - r);
			const vec3<real> gradW = kernel.gradW(dr, h);
			#if DERIVATIVE == 1
				Minv += tensor(dr, dr) * kernel(dr, h) * ngb[j]->mass / ngb[j]->dens;
			#elif DERIVATIVE == 2
				Minv += tensor(dr, gradW) * ngb[j]->mass / ngb[j]->dens;
			#endif
			#if 0 //ADD gradh
				gradh += - ngb[j]->mass * (Ndim * kernel(dr, h) + inner(dr, gradW)) / h;
				//gradh += - ngb[j]->mass * ngb[j]->eng * (Ndim * kernel(dr, h) + inner(dr, gradW)) / h;
			#endif
		}
		#warning TEMPORARY!!
		if(Ndim == 3){
			Minv = Minv.inv();
		}else if(Ndim == 2){
			Minv.zz = 1.0;
			Minv = Minv.inv();
		}else{
			Minv.yy = Minv.zz = 1.0;
			Minv = Minv.inv();
		}
		gradh = 1.0 + h / (Ndim * dens) * gradh;
		//gradh = 1.0 + h / (Ndim * psmth) * gradh;
	}
	inline void setDivRotVelocity(const kernel_t<real>& kernel, const boundary_t<real>& boundary, const int Ndim){
		div_v = 0;
		rot_v = 0;
		dya_v = 0;
		dya_a = 0;
		div_div_v = 0;
		isInShock = 0.0;
		#if DERIVATIVE == 0
		#warning SPH type Derivative
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(r - ngb[j]->r);
			const vec3<real> dv = ngb[j]->v - v;
			const vec3<real> da = ngb[j]->getSPHAcc() - getSPHAcc();
			const vec3<real> gradW = kernel.gradW(dr, h);			
			div_v += ngb[j]->mass / dens * inner (dv, gradW);
			rot_v += ngb[j]->mass / dens * outer (dv, gradW);
			dya_v += ngb[j]->mass / dens * tensor(dv, gradW);
			dya_a += ngb[j]->mass / dens * inner (da, gradW);
			
		}
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(r - ngb[j]->r);
			div_div_v += ngb[j]->div_v * kernel.gradW(dr, h) * ngb[j]->mass / ngb[j]->dens;
			isInShock += ngb[j]->mass / dens * math::sign(ngb[j]->div_v) * kernel(dr, h);
		}
		#elif DERIVATIVE == 1
		#warning Moment Based Gradient Derivative
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(ngb[j]->r - r);
			const vec3<real> dv = ngb[j]->v - v;
			const vec3<real> da = ngb[j]->getSPHAcc() - getSPHAcc();
			const vec3<real> psi_ji = Minv * (dr * kernel(dr, h)) * ngb[j]->mass / ngb[j]->dens;
			div_v += inner (dv, psi_ji);
			rot_v += outer (dv, psi_ji);
			dya_v += tensor(dv, psi_ji);
			dya_a += tensor(da, psi_ji);
		}
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(ngb[j]->r - r);
			const vec3<real> psi_ji = Minv * (dr * kernel(dr, h)) * ngb[j]->mass / ngb[j]->dens;
			div_div_v += (ngb[j]->div_v - div_v) * psi_ji;
			isInShock += ngb[j]->mass / dens * math::sign(ngb[j]->div_v) * kernel(dr, h);
		}
		#elif DERIVATIVE == 2
		#warning Linear Exact type Derivative
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(ngb[j]->r - r);
			const vec3<real> dv = ngb[j]->v - v;
			const vec3<real> da = ngb[j]->getSPHAcc() - getSPHAcc();
			const vec3<real> psi_ji = Minv * kernel.gradW(dr, h) * ngb[j]->mass / ngb[j]->dens;
			div_v += inner (dv, psi_ji);
			rot_v += outer (dv, psi_ji);
			dya_v += tensor(dv, psi_ji);
			dya_a += tensor(da, psi_ji);
		}
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(ngb[j]->r - r);
			const vec3<real> psi_ji = Minv * kernel.gradW(dr, h) * ngb[j]->mass / ngb[j]->dens;
			div_div_v += (ngb[j]->div_v - div_v) * psi_ji;
			isInShock += math::sign(ngb[j]->div_v) * kernel(dr, h) * ngb[j]->mass / dens;
		}
		#elif DERIVATIVE == 3
		#warning MPS type Derivative
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(ngb[j]->r - r);
			const vec3<real> dv = ngb[j]->v - v;
			const vec3<real> da = ngb[j]->getSPHAcc() - getSPHAcc();
			const real dr2 = abs2(dr) + 1.0e-16;
			div_v += inner (dv, dr) * ngb[j]->mass * / (dr2 * ngb[j]->dens) * kernel(dr, h);
			rot_v += outer (dv, dr) * ngb[j]->mass * / (dr2 * ngb[j]->dens) * kernel(dr, h);
			dya_v += tensor(dv, dr) * ngb[j]->mass * / (dr2 * ngb[j]->dens) * kernel(dr, h);
			dya_a += tensor(da, dr) * ngb[j]->mass * / (dr2 * ngb[j]->dens) * kernel(dr, h);
		}
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(ngb[j]->r - r);
			const real dr2 = abs2(dr) + 1.0e-16;
			div_div_v += (ngb[j]->div_v - div_v) * dr * ngb[j]->mass / (dr2 * ngb[j]->dens) * kernel(dr, h);
			isInShock += math::sign(ngb[j]->div_v) * kernel(dr, h) * ngb[j]->mass / dens;
		}
		#else
		#error DERIVATIVE ERROR
		#endif

		#if AV_LIMITTER == 0
		#warning NO AV. LIMITTER
		sAV = 1.0;
		#elif AV_LIMITTER == 1
		sAV = math::abs(div_v) / (math::abs(div_v) + abs(rot_v) + snds / h * 1.0e-4);
		#elif AV_LIMITTER == 2
		const mat33<real> shear = 0.5 * (dya_v + transpose(dya_v)) - div_v * mat33<real>::I() / (real)(Ndim);
		sAV = math::square(std::abs(2.0 * pow(1.0 - isInShock, 4) * div_v));
		sAV = sAV / (sAV + (shear * transpose(shear)).trace() + snds / h * 1.0e-4);
		#else
		#error AV_LIMITTER ERROR
		#endif
		qAV = (div_v < 0) ? (- aAV * snds * h * div_v + 2.0 * aAV * h * h * div_v * div_v) * dens : 0.0;
		qAV *= sAV;
	}
	inline void setTimestep(){
		dt = std::min(PARAM::C_CFL * sqrt(h / abs(a)), dt);
		//dt = std::min(PARAM::C_CFL * eng / std::abs(eng_dot), dt);
	}
	inline void setTimeDerivative(const kernel_t<real>& kernel, const boundary_t<real>& boundary){
		//clear
		a = 0;
		eng_dot = 0;
		dt = 1.0e+30;

		AVa = 0;
		AVeng_dot = 0;
		if(__builtin_expect(type != HYDRO, 0)) return;
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(r - ngb[j]->r);
			const vec3<real> dv = v - ngb[j]->v;
			real AV = 0;
			int AVflag = 0;
			//mat33<real> AV = 0;
			const real dvdr = inner(dv, dr);
			if(dvdr < 0){
				AVflag = 1;
				#if AV_TYPE == 0
				#warning Monaghan & Gingold 1983
				const real mu = 0.5 * (h + ngb[j]->h) * dvdr / (abs2(dr) + 1.0e-4 * h * ngb[j]->h);
				AV = (- 0.5 * (aAV + ngb[j]->aAV) * mu * 0.5 * (snds + ngb[j]->snds) + 2.0 * 0.5 * (aAV + ngb[j]->aAV) * mu * mu) / (0.5 * (dens + ngb[j]->dens));
				AV *= 0.5 * (sAV + ngb[j]->sAV);
				dt = PARAM::C_CFL * h / (h * math::abs(div_v) + snds + 1.2 * (aAV * snds + 2.0 * aAV * h * math::abs(math::min(div_v, 0.0))));
				#elif AV_TYPE == 1
				#warning Monaghan 1997
				const real w_ij = dvdr / abs(dr);
				const real v_sig = snds + ngb[j]->snds - 3.0 * w_ij;
				AV = -0.5 * 0.5 * (aAV + ngb[j]->aAV) * v_sig * w_ij / (0.5 * (dens + ngb[j]->dens));
				AV *= 0.5 * (sAV + ngb[j]->sAV);
				dt = std::min(2.0 * PARAM::C_CFL * h / v_sig, dt);
				#elif AV_TYPE == 3
				#warning Lattanzio & Monaghan 1985
				const real mu = 0.5 * (h + ngb[j]->h) * dvdr / (abs2(dr) + 1.0e-4 * h * ngb[j]->h) * 0.5 * (sAV + ngb[j]->sAV) / (0.5 * (snds + ngb[j]->snds));
				AV = (pres / (dens * dens) +  ngb[j]->pres / (ngb[j]->dens * ngb[j]->dens)) * (0.5 * (aAV + ngb[j]->aAV)) * (- mu + 2.0 * mu * mu);
				dt = PARAM::C_CFL * h / (h * math::abs(div_v) + snds + 1.2 * (aAV * snds + 2.0 * aAV * h * math::abs(math::min(div_v, 0.0))));
				#else
				//local type AV
				#endif
				//Owen 2004?
				//mat33<real> sigma = 0.5 * (h + ngb[j]->h) * tensor(dr, dv) / (abs2(dr) + 1.0e-4 * h * ngb[j]->h);
				//AV = (- 0.5 * (aAV + ngb[j]->aAV) * 0.5 * (snds + ngb[j]->snds) * sigma + 2.0 * 0.5 * (aAV + ngb[j]->aAV) * sigma * sigma) / (0.5 * (dens + ngb[j]->dens));
			}
			const vec3<real> igradW = kernel.gradW(dr, h) / gradh;
			const vec3<real> jgradW = kernel.gradW(dr, ngb[j]->h) / ngb[j]->gradh;
			const vec3<real> gradW = 0.5 * (igradW + jgradW);
			
			#if (AV_TYPE == 0 || AV_TYPE == 1 || AV_TYPE == 3)
			AVa -= ngb[j]->mass * AV * gradW;
			a   -= (EoS->HeatCapacityRatio() - 1.0) * ngb[j]->mass * eng * ngb[j]->eng * (1.0 / ngb[j]->psmth * jgradW + 1.0 / psmth * igradW) + ngb[j]->mass * AV * gradW;

			AVeng_dot += ngb[j]->mass * 0.5 * inner(AV * dv /*+ TC*/, gradW);
			eng_dot   += (EoS->HeatCapacityRatio() - 1.0) * ngb[j]->mass * eng * ngb[j]->eng * (1.0 / psmth * inner(dv, igradW)) + 0.5 * ngb[j]->mass * inner(AV * dv, gradW);
			#elif (AV_TYPE == 2)
			#warning Local AV
			dt = PARAM::C_CFL * h / (h * math::abs(div_v) + snds + 1.2 * (aAV * snds + 2.0 * aAV * h * math::abs(math::min(div_v, 0.0))));
			AVa -= AVflag * ngb[j]->mass * (ngb[j]->qAV / (ngb[j]->dens * ngb[j]->dens) * jgradW + qAV / (dens * dens) * igradW);
			a   -= (EoS->HeatCapacityRatio() - 1.0) * ngb[j]->mass * eng * ngb[j]->eng * (1.0 / ngb[j]->psmth * jgradW + 1.0 / psmth * igradW) + AVflag * ngb[j]->mass * ngb[j]->qAV / (ngb[j]->dens * ngb[j]->dens) * jgradW + AVflag * mass * qAV / (dens * dens) * igradW;

			AVeng_dot += AVflag * ngb[j]->mass * qAV / (dens * dens) * inner(dv, igradW);
			eng_dot += (EoS->HeatCapacityRatio() - 1.0) * ngb[j]->mass * eng * ngb[j]->eng * (1.0 / psmth * inner(dv, igradW)) + AVflag * ngb[j]->mass * qAV / (dens * dens) * inner(dv, igradW);
			#elif (AV_TYPE == 4)
			#warning NO AV.
			dt = PARAM::C_CFL * h / (h * math::abs(div_v) + snds + 1.2 * (aAV * snds + 2.0 * aAV * h * math::abs(math::min(div_v, 0.0))));
			AV = 0;
			AVa -= ngb[j]->mass * AV * gradW;
			a   -= (EoS->HeatCapacityRatio() - 1.0) * ngb[j]->mass * eng * ngb[j]->eng * (1.0 / ngb[j]->psmth * jgradW + 1.0 / psmth * igradW) + ngb[j]->mass * AV * gradW;

			AVeng_dot += ngb[j]->mass * 0.5 * inner(AV * dv, gradW);
			eng_dot += (EoS->HeatCapacityRatio() - 1.0) * ngb[j]->mass * eng * ngb[j]->eng * (1.0 / psmth * inner(dv, igradW)) + 0.5 * ngb[j]->mass * inner(AV * dv, gradW);
			#else
			#error AV_TYPE ERROR(2).
			#endif
			if(ngb[j]->type == FREEZE){
				#if 1
				const real dr = abs(dr) + 1.0e-4 * h;
				a += 1.0 * dr * pow(sin(2.0 * math::pi * dr / h) / (2.0 * math::pi * dr / h) , 6.0);
				#else
				const real dr2 = abs2(dr);
				const real h_over_dr = h / sqrt(dr2);
				a += 0.5 * (pow(h_over_dr, 4.0) - pow(h_over_dr, 2.0)) * dr / dr2;
				#endif
			}
		}
		return ;
	}
	inline void setPressure(){
		pres = EoS->Pressure(dens, eng);
		snds = EoS->SoundSpeed(dens, eng);
		return ;
	}
	inline void setEoS(EoS::EoS_t<real> * const _EoS){
		EoS = _EoS;
		return ;
	}
	inline void createParticle(const vec3<real> _r, const vec3<real> _v, const real _mass, const real _dens, const real _eng, short _tag = 0){
		r = _r;
		v = _v;
		mass = _mass;
		dens = _dens;
		eng  = _eng;
		tag  = _tag;
	}
	inline real getSpecificInternalEnergy() const{
		return eng;
	}
	inline vec3<real> getHydroAcc() const{
		return a - AVa - ext_a;
	}
	inline vec3<real> getSPHAcc() const{
		return a - ext_a;
	}
	inline const ieSPH* const getPointer() const{
		return this;
	}
	//defined in init/header.
	inline void setSource();
	//integral
	inline void initialKick(const real dt_glb){
		v_half   = v + a * dt_glb * 0.5;
		eng_half = eng + eng_dot * dt_glb * 0.5;
		eng_half = std::max(eng_half, 0.0);
	}
	inline void fullDrift(const real dt_glb, const boundary_t<real>& boundary){
		r += v_half * dt_glb;
		r = boundary.PeriodicWrapping(r);
		if(type != HYDRO) return;

		const real aAVmin = 0.1;
		const real aAVmax = 2.0;
		const real decay  = 1.0;
		#if TIME_VAR_AV == 1
		#warning TIME VAR. AV by MM97
		#if 1
		aAV += (- (aAV - aAVmin) / (2.0 * decay * h / snds) + (aAVmax - aAV) * std::max(- div_v, 0.0)) * dt_glb;
		#else
		real v_sig_max = 0.0;
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(r - ngb[j]->r);
			const vec3<real> dv = v - ngb[j]->v;
			const real dvdr_hat = inner(dv, dr) / abs(dr);
			v_sig_max = std::max(v_sig_max, 0.5 * (snds + ngb[j]->snds) - std::min(dvdr_hat, 0.0));
		}
		aAV += (- (aAV - aAVmin) / (2.0 * decay * h / v_sig_max) + sAV * (aAVmax - aAV) * std::max(- div_v, 0.0)) * dt_glb;
		#endif
		#elif TIME_VAR_AV == 2
		#warning TIME VAR. AV by CD10
		const real A = sAV * std::max(- (dya_a - dya_v * dya_v).trace(), 0.0);//eq. 13 in CD10
		real v_sig_max = 0.0;
		for(std::size_t j = 0 ; j < ngb.size() ; ++ j){
			const vec3<real> dr = boundary.Periodic(r - ngb[j]->r);
			const vec3<real> dv = v - ngb[j]->v;
			const real dvdr_hat = inner(dv, dr) / abs(dr);
			v_sig_max = std::max(v_sig_max, 0.5 * (snds + ngb[j]->snds) - std::min(dvdr_hat, 0.0));
		}
		const real aAVloc = aAVmax * h * h * A / (h * h * A + v_sig_max * v_sig_max);
		if(aAV < aAVloc){
			aAV = aAVloc;
		}else{
			aAV += - (aAV - std::max(aAVmin, aAVloc)) / (2.0 * decay * h / v_sig_max) * dt_glb;
		}
		#elif TIME_VAR_AV == 3
		#warning TIME VAR. AV by RH12
		const real aAVloc = (div_v < 0) ? sAV * aAVmax * (h * h * abs(div_div_v) / (h * h * abs(div_div_v) + h * math::abs(div_v) + 0.05 * snds)) : 0.0;
		if(aAV < aAVloc){
			aAV = aAVloc;
		}else{
			aAV += - (aAV - std::max(aAVmin, aAVloc)) / (2.0 * decay * h / snds) * dt_glb;
		}
		#elif TIME_VAR_AV == 4
		#warning TIME VAR. AV by R15
		const real A = sAV * std::max(h * (dya_a - dya_v * dya_v).trace() * abs(div_div_v), 0.0);
		const real aAVloc = aAVmax * A / (A + 0.1 * pow(snds / h, 3));
		if(aAV < aAVloc){
			aAV = aAVloc;
		}else{
			aAV += - (aAV - std::max(aAVmin, aAVloc)) / (2.0 * decay * h / snds) * dt_glb;
		}
		#else
		#warning TIME CONST AV.
		#endif
	}
	inline void predict(const real dt_glb){
		v   = v   + a * dt_glb;
		eng = eng + eng_dot * dt_glb;
		eng = std::max(eng, 0.0);
	}
	inline void finalKick(const real dt_glb){
		v   = v_half   + a * dt_glb * 0.5;
		eng = eng_half + eng_dot * dt_glb * 0.5;
		eng = std::max(eng, 0.0);
	}
};

template <typename real> struct system_t{
	real time, dt;
	std::size_t step;
	vec3<real> lMom, aMom;//linear and angular momentum;
	real mass, tEng;
	system_t() : time(0.0), dt(math::NaN), step(0){
	}
	template <class Ptcl> void setTimeStep(const std::vector<Ptcl>& ptcl){
		dt = 1.0e+30;
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			dt = std::min(dt, ptcl[i].dt);
		}
	}
	template <class Ptcl> void conservativeCheck(const std::vector<Ptcl>& ptcl){
		mass = 0;
		lMom = 0;
		aMom = 0;
		tEng = 0;
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			mass += ptcl[i].mass;
			lMom += ptcl[i].v * ptcl[i].mass;
			aMom += outer(ptcl[i].r, ptcl[i].v * ptcl[i].mass);
			tEng += (ptcl[i].getSpecificInternalEnergy() + 0.5 * abs2(ptcl[i].v)) * ptcl[i].mass;
		}
		printf("mass = %.16e\n", mass);
		printf("lMom = %.16e, %.16e, %.16e\n", lMom.x, lMom.y, lMom.z);
		printf("aMom = %.16e, %.16e, %.16e\n", aMom.x, aMom.y, aMom.z);
		printf("tEng = %.16e\n", tEng);
		return;
	}
};

class CPUtimer_t{
	double start, end;
	public:
	CPUtimer_t(){
		timeval tv;
		gettimeofday(&tv, NULL);
		start = (double)tv.tv_sec + (double)tv.tv_usec * 1.0e-6;
	}
	~CPUtimer_t(){
		timeval tv;
		gettimeofday(&tv, NULL);
		end = (double)tv.tv_sec + (double)tv.tv_usec * 1.0e-6;
		printf("%e [sec]\n", end - start);
	}
};

template <typename Ptcl> class FileIO{
	std::size_t step;
	double time;
	const double time_interval;
	public:
	FileIO(const double _time_interval) : step(0), time(0), time_interval(_time_interval){
	}
	void OutputFile(const vector<Ptcl>& ptcl, const system_t<double>& sys, const boundary_t<double>& boundary){
		if(sys.time < time && time_interval > 0.0) return;
		std::cout << "OUTPUT... " << step << std::endl;
		{
			char filename[256];
			sprintf(filename, "%s/%lu.txt", PARAM::OUTPUT_DIR, step);
			FILE* out = fopen(filename, "w");
			fprintf(out, "%e\n", sys.time);
			fprintf(out, "%lu\n", ptcl.size());
			for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
				fprintf(out, "%lu\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", ptcl[i].id, ptcl[i].r.x, ptcl[i].r.y, ptcl[i].r.z, ptcl[i].v.x, ptcl[i].v.y, ptcl[i].v.z, ptcl[i].dens, ptcl[i].getSpecificInternalEnergy(), ptcl[i].pres, ptcl[i].div_v, ptcl[i].rot_v.z, ptcl[i].aAV, ptcl[i].sAV, ptcl[i].psmth);
			}
			fclose(out);
		}
		{
			char filename[256];
			sprintf(filename, "%s/test%lu.txt", PARAM::OUTPUT_DIR, step);
			FILE* out = fopen(filename, "w");
			fprintf(out, "%e\n", sys.time);
			fprintf(out, "%lu\n", ptcl.size());
			for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
				const vec3<double> hydro_a = ptcl[i].getHydroAcc();
				fprintf(out, "%d\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", ptcl[i].tag, ptcl[i].r.x, ptcl[i].r.y, ptcl[i].r.z, ptcl[i].v.x, ptcl[i].v.y, ptcl[i].v.z, hydro_a.x, hydro_a.y, hydro_a.z, ptcl[i].AVa.x, ptcl[i].AVa.y, ptcl[i].AVa.z, ptcl[i].ext_a.x, ptcl[i].ext_a.y, ptcl[i].ext_a.z);
			}
			fclose(out);
		}
		if(0){
			char filename[256];
			sprintf(filename, "%s/a%lu.txt", PARAM::OUTPUT_DIR, step);
			FILE* out = fopen(filename, "w");
			fprintf(out, "%e\n", sys.time);
			fprintf(out, "%lu\n", ptcl.size());
			for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
				const vec3<double> hydro_a = ptcl[i].getHydroAcc();
				const double vr     = inner(ptcl[i].r, ptcl[i].v)   / abs(ptcl[i].r);
				const double vtheta = outer(ptcl[i].r, ptcl[i].v).z / abs2(ptcl[i].r);
				//const double ar_hy  = inner(ptcl[i].r, hydro_a) / abs(ptcl[i].r);
				//const double ar_hy  = inner(ptcl[i].r, ptcl[i].a) / abs(ptcl[i].r);
				const double ar_hy  = inner(ptcl[i].r, ptcl[i].ext_a) / abs(ptcl[i].r) + abs(ptcl[i].r) * vtheta * vtheta;
				const double ar_AV  = inner(ptcl[i].r, ptcl[i].AVa) / abs(ptcl[i].r);
				//const double atheta_hy  = outer(ptcl[i].r, hydro_a).z / abs2(ptcl[i].r);
				const double atheta_hy  = outer(ptcl[i].r, ptcl[i].a).z / abs2(ptcl[i].r) - 2.0 * vr * vtheta / abs(ptcl[i].r);
				const double atheta_AV  = outer(ptcl[i].r, ptcl[i].AVa).z / abs2(ptcl[i].r) - 2.0 * vr * vtheta / abs(ptcl[i].r);
				fprintf(out, "%lu\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", ptcl[i].id, ptcl[i].r.x, ptcl[i].r.y, vr, vtheta, ar_hy, ar_AV, atheta_hy, atheta_AV);
			}
			fclose(out);
		}
		++ step;
		time += time_interval;
	}
	void OutputBinary(const vector<Ptcl>& ptcl, const system_t<double>& sys, const boundary_t<double>& boundary){
		std::ofstream fout;
		fout.open("log.bin", std::ios::out | std::ios::binary | std::ios::trunc);
		if(!fout){
			cout << "cannot open restart file." << endl;
			exit(1);
		}
		fout.write(reinterpret_cast<const char * const>(&sys), sizeof(system_t<double>));
		fout.write(reinterpret_cast<const char * const>(&boundary), sizeof(boundary_t<double>));
		for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
			Ptcl ith = ptcl[i];
			//destruct ith.ngb;
			vector<const Ptcl* > dummy;
			ith.ngb.swap(dummy);
			fout.write((char*)&ith, sizeof(Ptcl));
		}
		fout.close();
	}
	void InputBinary(vector<Ptcl>& ptcl, system_t<double>* sys, boundary_t<double>* boundary){
		ptcl.clear();
		std::ifstream fin("log.bin", std::ios::in | std::ios::binary);
		if(!fin){
			cout << "cannot open restart file." << endl;
			exit(1);
		}
		fin.read((char*)sys, sizeof(system_t<double>));
		fin.read((char*)boundary, sizeof(boundary_t<double>));
		while(1){
			Ptcl ith;
			fin.read((char*)&ith, sizeof(Ptcl));
			ith.ngb.clear();
			if(fin.eof() == true) break;
			ptcl.push_back(ith);
		}
		fin.close();
		if(time_interval <= 0.0) return;
		while(this->time < sys->time){
			this->time += time_interval;
			++ this->step;
		}
	}
};

template <class Ptcl> bool removeParticle(std::vector<Ptcl>& ptcl){
	bool isRemoved = false;
	for(std::size_t i = 0 ; i < ptcl.size() ; ++ i){
		if(ptcl[i].eng < 0){
			std::cout << i << ": Negative eng" << std::endl;
			ptcl[i] = ptcl[ptcl.size()];
			ptcl.pop_back();
			isRemoved = true;
		}
	}
	return isRemoved;
}

