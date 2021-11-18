/*
// For contact information
// Yasuhiro Inoue
// inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#ifndef _PARAMETERS_H
#define _PARAMETERS_H

// 外部ファイルから読み込む値。
// _________________________
// __Cells__________________
extern double k_volume;
extern double k_area;
extern double k_angle_pre;
extern double k_length;
extern double ene_repulsive;	// 粒子間の反発エネルギー
extern double ene_adhesion_ii;	// 接着のエネルギー (粒子間に定義するが、平面間のエネルギーを粒子に離散化している)
extern double ene_adhesion_io;	// 接着のエネルギー (粒子間に定義するが、平面間のエネルギーを粒子に離散化している)
extern double ene_adhesion_oo;	// 接着のエネルギー (粒子間に定義するが、平面間のエネルギーを粒子に離散化している)
extern double frc_t;	// Tangential force		// 接線方向の力
extern double frc_c;	// Centripetal force	// 中心方向の力
extern double frc_z;	// Dorsalward force	 	// 背側方向の力
extern double swelling_ratio;	// 細胞の膨潤率
extern double p_stay_to_move;
extern double p_move_to_stay;
// _________________________
// __Bounds_________________
extern double owall_a;
extern double owall_b;
extern double iwall_a;
extern double iwall_b;
// _________________________
// __Settings_______________
extern unsigned int seed;
extern int rechoice_step;
// _________________________
// __Other_parameters_______
extern int end_step;
extern double dt_vtk;	// 小さくするとvtkが多く出てくる
extern double a_tol;
extern double r_tol;
extern double kappa_d;


extern double dt;	// 可変刻み時間
extern double dt_step;	// ステップ間の時間幅
extern double physical_time;
extern double count_physical_time;
extern long int total_n_revcell;	// 移動細胞数(のべ)
extern long double total_work_extf;	// 外力が系にした仕事の積算(累計)値


// 一度決めればあまり変えないだろうと思うので直接値を書き込む。
constexpr double sigma = 1.0;	// 細胞の直径
constexpr double dist_rep_cutoff = 0.03;
constexpr double dist_adh_cutoff = 0.1;

constexpr int num_thread = 8;

constexpr double coeff_friction = 1.0;
//constexpr double dt = 1e-12;	// パラメータによっては小さくしないと発散する。長さが極端に短い辺が見られ、その辺がなす面がおかしくなることから、頂点が入れ替わり幾何学的におかしくなる？
//constexpr double dt_vtk = 5e-2;	// 小さくするとvtkが多く出てくる
constexpr double eps = 1e-20;	// 0除算を回避するための無限小。とりあえずdoubleの有効数字より小さくしておく。

// Gridに区切る空間全体の大きさ(大きめにとる)
constexpr double sys_size_x = 7.0;
constexpr double sys_size_y = 4.0;
constexpr double sys_size_z = 4.0;
// Grid(立方体)の一辺の長さ
constexpr double grid_length = 0.105;	// dist_cutoff, dist_adh_cutoff より少し大きい値とする

// for Dormand-Prince Method
// 性能を出すにはここのチューニングが大事。力学パラメータのオーダーとの兼ね合いで変わってくると思われる。
// バーテックスで作った折り畳みの場合、k_area=2.0 k_angle_pre=2.0程度だと許容誤差を1e-7程度にしとくのが効率よさそう。1e-8だと効率があまり良くない。1e-6とか1e-5にしたら発散した。
//constexpr double a_tol = 1e-5;		// absolute tolerance of error
//constexpr double r_tol = 1e-5;		// relative tolerance of error
constexpr double fac_max = 4.0;		// 元のdtの4倍を超える刻み幅の更新を禁止する
constexpr double fac_min = 0.25;	// 元のdtの1/4を下回る刻み幅の更新を禁止する
constexpr double fac_save = 0.80;	// 数値的な安定を考えて安全率をとる。

// パラメータ読み込み関数
void readParametersFile();

#endif