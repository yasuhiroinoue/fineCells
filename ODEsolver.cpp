/*
// For contact information
// Yasuhiro Inoue
// inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
// 一階微分方程式のソルバ。1,2,4次精度及び4-5次の埋め込み型ルンゲクッタ。いずれも陽解法。

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "ODEsolver.h"
#include "_class.h"
#include "_parameters.h"
#include "vec.h"

#define OPENMP_ODE

namespace ode{

//Euler Method
void firstMotionSolver(){
	//#ifdef OPENMP_ODE
	//#pragma omp parallel num_threads(num_thread)
	//#endif
	{
		//現在位置を格納
		#ifdef OPENMP_ODE
		//#pragma omp for
		#pragma omp parallel for num_threads(num_thread)	//設計ミス。力の計算をGlobalクラスに突っ込んでいるので分けざるを得ない。
		#endif
		for(int i = 0; i < p_g->n_v; i++){
			p_g->getVertex(i)->copyLocationToLocationT();
		}

		//1st translocation
		p_g->calcForces();
		#ifdef OPENMP_ODE
		#pragma omp parallel for num_threads(num_thread)
		#endif
		for(int i = 0; i < p_g->n_v; i++){
			Vertex* p_v = p_g->getVertex(i);
			//p_v->copyForceRk(0);
			_vec<double> tmp_loc;
			tmp_loc = p_v->getLocationT() + (p_v->getForce()*dt)/coeff_friction;
			//cout << tmp_loc.x << " " << tmp_loc.y << " " << tmp_loc.z << endl;
			p_v->setLocation(tmp_loc);
		}
	}
}

//Modified Euler Method / 2nd-order Runge-Kutta
void secondMotionSolver(){
	//#ifdef OPENMP_ODE
	//#pragma omp parallel num_threads(num_thread)
	//#endif
	{
		//現在位置を格納
		#ifdef OPENMP_ODE
		#pragma omp parallel for num_threads(num_thread)
		#endif
		for(int i = 0; i < p_g->n_v; i++){
			p_g->getVertex(i)->copyLocationToLocationT();
		}
		//1st translocation
		p_g->calcForces();
		#ifdef OPENMP_ODE
		#pragma omp parallel for num_threads(num_thread)
		#endif
		for(int i = 0; i < p_g->n_v; i++){
			Vertex* p_v = p_g->getVertex(i);
			//p_v->copyForceRk(0);
			_vec<double> tmp_loc;
			tmp_loc = p_v->getLocationT() + (p_v->getForce()*0.5*dt)/coeff_friction;
			p_v->setLocation(tmp_loc);
		}
		//2nd translocation
		p_g->updateGeometry();	//追加　2017/07/21
		p_g->calcForces();
		#ifdef OPENMP_ODE
		#pragma omp parallel for num_threads(num_thread)
		#endif
		for(int i = 0; i < p_g->n_v; i++){
			Vertex* p_v = p_g->getVertex(i);
			//p_v->copyForceRk(1);
			_vec<double> tmp_loc;
			tmp_loc = p_v->getLocationT() + (p_v->getForce()*dt)/coeff_friction;
			p_v->setLocation(tmp_loc);
		}
	}
}

//4th-order Runge-Kutta
void fourthMotionSolver(){
	//現在位置を格納
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		p_g->getVertex(i)->copyLocationToLocationT();
	}
	//1st translocation
	p_g->calcForces();
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* p_v = p_g->getVertex(i);
		p_v->copyForceRk(0);
		_vec<double> tmp_loc;
		tmp_loc = p_v->getLocationT() + (p_v->getForce()*0.5*dt)/coeff_friction;
		p_v->setLocation(tmp_loc);
	}
	//2nd translocation
	p_g->updateGeometry();	//追加　2017/07/21
	p_g->calcForces();
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* p_v = p_g->getVertex(i);
		p_v->copyForceRk(1);
		_vec<double> tmp_loc;
		tmp_loc = p_v->getLocationT() + (p_v->getForce()*0.5*dt)/coeff_friction;
		p_v->setLocation(tmp_loc);
	}
	//3rd translocation
	p_g->updateGeometry();	//追加　2017/07/21
	p_g->calcForces();
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* p_v = p_g->getVertex(i);
		p_v->copyForceRk(2);
		_vec<double> tmp_loc;
		tmp_loc = p_v->getLocationT() + (p_v->getForce()*dt)/coeff_friction;
		p_v->setLocation(tmp_loc);
	}
	//4th translocation
	p_g->updateGeometry();	//追加　2017/07/21
	p_g->calcForces();
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* p_v = p_g->getVertex(i);
		p_v->copyForceRk(3);
		_vec<double> tmp_loc;
		tmp_loc = p_v->getForceRk(0) + 2.0*p_v->getForceRk(1) + 2.0*p_v->getForceRk(2) + p_v->getForceRk(3);
		tmp_loc *= dt;
		tmp_loc /= (6.0*coeff_friction);
		tmp_loc += p_v->getLocationT();
		p_v->setLocation(tmp_loc);
	}
	physical_time += dt;
	count_physical_time += dt;
}

//埋め込み型Runge-Kutta(Dormand-Prince Method)。dtが可変の場合。
//バーテックスで作った折り畳みは初期条件がシビアすぎて、刻み幅を小さくしないと発散するし、1e-12とか一定にしてたら時間がかかりすぎるので、これでいい感じの刻み幅にしようという試み。
bool coreDormandPrinceSolver(){
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		p_g->getVertex(i)->copyLocationToLocationT();
	}

	//1st translocation
	p_g->calcForces();
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* p_v = p_g->getVertex(i);
		p_v->copyForceRk(0);
		_vec<double> tmp_loc;
		tmp_loc = p_v->getLocationT() + (p_v->getForce()*0.2*dt)/coeff_friction;
		p_v->setLocation(tmp_loc);
	}
	//2nd translocation
	p_g->updateGeometry();	//追加　2017/07/21
	p_g->calcForces();
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* p_v = p_g->getVertex(i);
		p_v->copyForceRk(1);
		_vec<double> tmp_loc;
		tmp_loc = p_v->getForce()*9.0/40.0 + p_v->getForceRk(0)*3.0/40.0;
		tmp_loc *= dt/coeff_friction;
		tmp_loc += p_v->getLocationT();
		p_v->setLocation(tmp_loc);
	}
	//3rd translocation
	p_g->updateGeometry();	//追加　2017/07/21
	p_g->calcForces();
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* p_v = p_g->getVertex(i);
		p_v->copyForceRk(2);
		_vec<double> tmp_loc;
		tmp_loc = p_v->getForce()*32.0/9.0 - p_v->getForceRk(1)*56.0/15.0 + p_v->getForceRk(0)*44.0/45.0;
		tmp_loc *= dt/coeff_friction;
		tmp_loc += p_v->getLocationT();
		p_v->setLocation(tmp_loc);
	}
	//4th translocation
	p_g->updateGeometry();	//追加　2017/07/21
	p_g->calcForces();
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* p_v = p_g->getVertex(i);
		p_v->copyForceRk(3);
		_vec<double> tmp_loc;
		tmp_loc = p_v->getForce()*(-212.0)/729.0 + p_v->getForceRk(2)*64448.0/6561.0 - p_v->getForceRk(1)*25360.0/2187.0 + p_v->getForceRk(0)*19372.0/6561.0;
		tmp_loc *=dt/coeff_friction;
		tmp_loc += p_v->getLocationT();
		p_v->setLocation(tmp_loc);
	}
	//5th translocation
	p_g->updateGeometry();	//追加　2017/07/21
	p_g->calcForces();
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* p_v = p_g->getVertex(i);
		p_v->copyForceRk(4);
		_vec<double> tmp_loc;
		tmp_loc = p_v->getForce()*(-5103.0)/18656.0 + p_v->getForceRk(3)*49.0/176.0 + p_v->getForceRk(2)*46732.0/5247.0 - p_v->getForceRk(1)*355.0/33.0 + p_v->getForceRk(0)*9017.0/3168.0;
		tmp_loc *= dt/coeff_friction;
		tmp_loc += p_v->getLocationT();
		p_v->setLocation(tmp_loc);
	}
	//6th translocation / 4th-order Runge-Kutta
	p_g->updateGeometry();	//追加　2017/07/21
	p_g->calcForces();
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* p_v = p_g->getVertex(i);
		p_v->copyForceRk(5);
		_vec<double> tmp_loc;
		tmp_loc = p_v->getForce()*11.0/84.0 - p_v->getForceRk(4)*2187.0/6784.0 + p_v->getForceRk(3)*125.0/192.0 + p_v->getForceRk(2)*500.0/1113.0 + p_v->getForceRk(0)*35.0/384.0;
		tmp_loc *= dt/coeff_friction;
		tmp_loc += p_v->getLocationT();
		p_v->setLocation(tmp_loc);
		p_v->setLocation4th(tmp_loc);
	}
	//7th translocation / 5th-order Runge-Kutta
	p_g->updateGeometry();	//追加　2017/07/21
	p_g->calcForces();
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* p_v = p_g->getVertex(i);
		p_v->copyForceRk(6);
		_vec<double> tmp_loc;
		tmp_loc = p_v->getForce()*1.0/40.0 + p_v->getForceRk(5)*187.0/2100.0 - p_v->getForceRk(4)*92097.0/339200.0 + p_v->getForceRk(3)*393.0/640.0 + p_v->getForceRk(2)*7571.0/16695.0 + p_v->getForceRk(0)*5179.0/57600.0;
		tmp_loc *= dt/coeff_friction;
		tmp_loc += p_v->getLocationT();
		p_v->setLocation(tmp_loc);
	}
	//return true;	//for debug.

	//evaluate error of ODE solution & estimate adaptive time step
	double dt_opt;
	double err = 0.0;
	#ifdef OPENMP_ODE
	#pragma omp parallel for num_threads(num_thread) reduction(+:err)
	#endif
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* p_v  = p_g->getVertex(i);
		_vec<double> scale_vec;
		_vec<double> loc_5th = p_v->getLocation();
		_vec<double> loc_4th = p_v->getLocation4th();
		//scale_vec.x = a_tol;
		//scale_vec.y = a_tol;
		//scale_vec.z = a_tol;
		//scale_vec = a_tol + r_tol*(p_v->getLocation().norm());
		scale_vec.x = a_tol + r_tol*std::max(fabs(loc_5th.x), fabs(loc_4th.x));
		scale_vec.y = a_tol + r_tol*std::max(fabs(loc_5th.y), fabs(loc_4th.y));
		scale_vec.z = a_tol + r_tol*std::max(fabs(loc_5th.z), fabs(loc_4th.z));
		_vec<double> err_vec = loc_5th - loc_4th;
		err_vec.x /= scale_vec.x;
		err_vec.y /= scale_vec.y;
		err_vec.z /= scale_vec.z;
		
		err += err_vec.sqr();
	}
	
	err /= (3.0*p_g->n_v);
	err = sqrt(err);
//	std::cout << "error = " << err << ":\t" << "dt = " << dt << std::endl;
	double dt_new = dt;	//次のステップで使う刻み幅
	dt_opt = fac_save*pow(err, -0.2);	//計算される最適な刻み幅はerr^(-1/5)*dtだけど、数値的な安定のために安全率を掛ける。
	dt_opt = std::max(fac_min, std::min(fac_max, dt_opt));	//急激な刻み幅の更新を避ける
	dt_new *= dt_opt;
	if(dt_new > 0.1*dt_vtk){
		dt_new = 0.1*dt_vtk;	// dtが発散しないようにする
	}
	
	if(err < 1.0){
		physical_time += dt;
		count_physical_time += dt;
//		cout << "error = " << err << ":\t" << "time = " << physical_time << ":\t" << "dt = " << dt << endl;

		dt_step = dt;	// v17	// 積分誤差の判定でOKが出たため、このdtは実際のステップ幅として認められる。その値はmakeRandomChoice()で参照される必要があるので、更新されてしまう前に保存しておく。

		dt = dt_new;

		return true;
	} else{
		#ifdef OPENMP_ODE
		#pragma omp parallel for num_threads(num_thread)
		#endif
		for(int i = 0; i < p_g->n_v; i++){
			Vertex* p_v = p_g->getVertex(i);
			p_v->setLocation(p_v->getLocationT());
		}
		p_g->updateGeometry();	//追加　2017/07/21
		dt = dt_new;

		return false;
	}
	
}

void solveMotionVertices(){
	//firstMotionSolver();	// Euler Method.
	//secondMotionSolver();	// 2nd-order Runge-Kutta(Modified Euler Method).
	//fourthMotionSolver();	// 4th-order Runge-Kutta.
	bool flag = false;
	while(!flag){
		flag = coreDormandPrinceSolver();	// 埋め込み型Runge-Kutta(Dormand-Prince Method)．dtが可変の場合．
	}
}

}