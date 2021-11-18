/*
// For contact information
// Yasuhiro Inoue
// inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "_parameters.h"

// _________________________
// __Cells__________________
double k_volume;
double k_area;
double k_angle_pre;
double k_length;
double ene_repulsive;
double ene_adhesion_ii;
double ene_adhesion_io;
double ene_adhesion_oo;
double frc_t;	// Tangential force		// 接線方向の力
double frc_c;	// Centripetal force	// 中心方向の力
double frc_z;	// Dorsalward force	 	// 背側方向の力
double swelling_ratio;
double p_stay_to_move;
double p_move_to_stay;
// _________________________
// __Bounds_________________
double owall_a;
double owall_b;
double iwall_a;
double iwall_b;
// _________________________
// __Settings_______________
unsigned int seed;
int rechoice_step;
// _________________________
// __Other_parameters_______
int end_step;
double dt_vtk;	// 小さくするとvtkが多く出てくる
double a_tol;
double r_tol;
double kappa_d;


//int vtk_step;
double dt = 1e-7;	// 可変刻み時間
double dt_step = 1e-7;	// ステップ間の時間幅	// 本当の初期値は知りようがない
double physical_time = 0.0;
double count_physical_time = 0.0;
long int total_n_revcell = 0;	// 移動細胞数(のべ)
long double total_work_extf = 0.0;	// 外力が系にした仕事の積算(累計)値


void readParametersFile(){
	std::cout << "Reading Parameter File...";

	std::ifstream fin("parameters.dat");
	if(!fin){
		std::cout << "Error. cannot open the paramters file." << std::endl;
		exit(1);
	}
	char dummy[100];


	fin >> dummy;	// _________________________
	fin >> dummy;	// __Cells__________________
	fin >> dummy >> k_volume >> k_area >> k_angle_pre >> k_length;
	fin >> dummy >> ene_repulsive;
	fin >> dummy >> ene_adhesion_ii >> ene_adhesion_io >> ene_adhesion_oo;
	fin >> dummy >> frc_t >> frc_c >> frc_z;
	fin >> dummy >> swelling_ratio;
	fin >> dummy >> p_stay_to_move;
	fin >> dummy >> p_move_to_stay;	

	fin >> dummy;	// _________________________
	fin >> dummy;	// __Bounds_________________
	fin >> dummy >> owall_a >> owall_b;
	fin >> dummy >> iwall_a >> iwall_b;

	fin >> dummy;	// _________________________
	fin >> dummy;	// __Settings_______________
	fin >> dummy >> seed;
	fin >> dummy >> rechoice_step;

	fin >> dummy;	// _________________________
	fin >> dummy;	// __Other_parameters_______	
	fin >> dummy >> end_step;
	fin >> dummy >> dt_vtk;
	fin >> dummy >> a_tol >> r_tol;
	fin >> dummy >> kappa_d;


	std::cout << " done." << std::endl;
}