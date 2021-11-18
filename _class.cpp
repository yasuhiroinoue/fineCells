/*
// For contact information
// Yasuhiro Inoue
// inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
//----------------------------------------------------------------------------
//	_class.cpp
//----------------------------------------------------------------------------
//		クラスの名前とかは3Dバーテックスと合わせている。手癖の問題でその方が慣れている。
//		push関数とset関数の名前の違いは、push_backしているか否か。情報取得はget関数
//----------------------------------------------------------------------------

#include <vector>
#include <iostream>
#include <omp.h>
#include <time.h>
#include <algorithm>	// std::max, std::min
#include "_class.h"
#include "_parameters.h"
#include "vec.h"
#include "random.h"

#define OPENMP_CLASS

using namespace std;

Global::~Global(){
	// newで確保している以上は始末をつける。vectorはdelete[]で解放できないのでループで。
	for(int i = 0; i < n_v; i++){
		delete(p_v[i]);
	}
	for(int i = 0; i < n_l; i++){
		delete(p_l[i]);
	}
	for(int i = 0; i < n_s; i++){
		delete(p_s[i]);
	}
	for(int i = 0; i < n_c; i++){
		delete(p_c[i]);
	}
	p_v.clear();
	p_l.clear();
	p_s.clear();
	p_c.clear();
}
void Global::pushVertex(Vertex* p_v_tmp){
	p_v.push_back(p_v_tmp);
}
Vertex* Global::getVertex(int a){
	return p_v[a];
}
void Global::pushLine(Line* p_l_tmp){
	p_l.push_back(p_l_tmp);
}
Line* Global::getLine(int a){
	return p_l[a];
}
void Global::pushSurface(Surface* p_s_tmp){
	p_s.push_back(p_s_tmp);
}
Surface* Global::getSurface(int a){
	return p_s[a];
}
void Global::pushCell(Cell* p_c_tmp){
	p_c.push_back(p_c_tmp);
}
Cell* Global::getCell(int a){
	return p_c[a];
}
void Global::updateGeometry(){
	#ifdef OPENMP_CLASS
	#pragma omp parallel num_threads(num_thread)
	#endif
	{
		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < n_l; i++){
			p_l[i]->calcLength();	// 各辺の長さの計算
		}
		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < n_s; i++){
			p_s[i]->calcNormalAndAreaTriangle();	// 各面の法線方向および面積の計算
			p_s[i]->calcCentroidTriangle();			// 各面の重心の計算
		}
		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < n_c; i++){
			p_c[i]->calcVolumeCell();		// 各細胞の体積の計算
			p_c[i]->calcCentroidCell();		// 各細胞の重心の計算	// 前のstepからの重心移動量もここで求めている
			p_c[i]->judgeCentralOrNot();	// この細胞が中心部にあるかどうか判定
		}

		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < n_v; i++){
			p_v[i]->calcNormalAndAreaVertex();	// 各頂点の法線方向および占有面積の計算
		}
	}
}
void Global::calcForces(){
	#ifdef OPENMP_CLASS
	#pragma omp parallel num_threads(num_thread)
	#endif
	{
		// 力の変数を初期化
		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < n_v; i++){
			p_v[i]->setForceZero();
			p_v[i]->setForceZero(num_thread);
		}

		// 体積保存力
		#ifdef OPENMP_CLASS
		#pragma omp for nowait
		#endif
		for(int i = 0; i < n_c; i++){
			p_c[i]->calcVolumeForce();
		}

		// 面積保存力
		#ifdef OPENMP_CLASS
		#pragma omp for nowait
		#endif
		for(int i = 0; i < n_s; i++){
			p_s[i]->calcAreaForce();
		}

		// 二面角保存力
		#ifdef OPENMP_CLASS
		#pragma omp for nowait
		#endif
		for(int i = 0; i < n_l; i++){
			p_l[i]->calcAngleForce();
		}

		// 辺長保存力
		#ifdef OPENMP_CLASS
		#pragma omp for nowait
		#endif
		for(int i = 0; i < n_l; i++){
			p_l[i]->calcLengthForce();
		}

		// 斥力(排除体積効果)・引力(接着効果)
		#ifdef OPENMP_CLASS
		#pragma omp for nowait
		#endif
		for(int i = 0; i < n_v; i++){
			Vertex* vp = getVertex(i);
			// 同じgridに属するvertex(IDが自分のIDより小さいものだけ)をIDの大きいものから順にみていく
			int next_vid_sm = next_vid[vp->getId()];	// そのgridの中でvpの次に大きいvertexID
			while(next_vid_sm >= 0){
				Vertex* vq = getVertex(next_vid_sm);

				bool check1 = false;	// vpとvqが同じ細胞か？
				bool check2 = false;	// 同じ辺を共有しているか？
				if(vp->getCell() == vq->getCell()){
					check1 = true;
					for(int j = 0; j < vp->getLineSize(); j++){
						Line* lp = vp->getLine(j);
						if(lp->getVertex(0) == vq || lp->getVertex(1) == vq){
							check2 = true;	break;
						}
					}
				}
				if(!check2) calcRepulsiveForce_ij(vp, vq);		// vpとvqが同じ辺を共有していなければ 斥力を計算
				//if(!check1) calcAdhesiveForce_ij(vp, vq);		// vpとvqが同じ細胞でなければ 接着力を計算

				next_vid_sm = next_vid[next_vid_sm];
			}
			// gridIDのより小さい隣接gridに属するvertexをIDの大きいものから順にみていく
			int current_grid_id = vp->getGridId();	// vpが属するgridのID
			int tmp = 13*current_grid_id;

			for(int j = 0; j < 13; j++){
				int neighbor_grid_id = nb_list[tmp + j];
				if(neighbor_grid_id >= 0){
					int next_vid_nb = getHeadId(neighbor_grid_id);
					while(next_vid_nb >= 0){
						Vertex* vq = getVertex(next_vid_nb);

						bool check1 = false;	// vpとvqが同じ細胞か？
						bool check2 = false;	// 同じ辺を共有しているか？
						if(vp->getCell() == vq->getCell()){
							check1 = true;
							for(int j = 0; j < vp->getLineSize(); j++){
								Line* lp = vp->getLine(j);
								if(lp->getVertex(0) == vq || lp->getVertex(1) == vq){
									check2 = true;	break;
								}
							}
						}
						if(!check2) calcRepulsiveForce_ij(vp, vq);		// vpとvqが同じ辺を共有していなければ 斥力を計算
						//if(!check1) calcAdhesiveForce_ij(vp, vq);		// vpとvqが同じ細胞でなければ 接着力を計算

						next_vid_nb = next_vid[next_vid_nb];
					}
				}
			}
		}

		// 外力
		#ifdef OPENMP_CLASS
		#pragma omp for nowait
		#endif
		for(int i = 0; i < n_v; i++){
			p_v[i]->calcExternalForce();
		}
/*
		// 内力
		#ifdef OPENMP_CLASS
		#pragma omp for nowait
		#endif
		for(int i = 0; i < n_c; i++){
			p_c[i]->calcInternalForce();
		}


		// 細胞移動力
		#ifdef OPENMP_CLASS
		#pragma omp for nowait
		#endif
		for(int i = 0; i < n_c; i++){
			p_c[i]->setExternalForceZero();
			if(p_c[i]->getBoolRevcell()){
				p_c[i]->calcCellMigrationForce();
			}
		}

		// 境界条件
		#ifdef OPENMP_CLASS
		#pragma omp for nowait
		#endif
		for(int i = 0; i < n_v; i++){
			p_v[i]->calcBoundaryForce();
		}
*/
		//OPENMPで分けて計算していた力の統合
		#ifdef OPENMP_CLASS
		#pragma omp barrier
		#pragma omp for
		#endif
		for(int i = 0; i < n_v; i++){
			p_v[i]->sumForceOmp();
		}
	}
}
void Global::setEnergy(double a){
	u_total = a;
}
void Global::calcEnergy(){
	double tmp_u_vol = 0.0;
	double tmp_u_are = 0.0;
	double tmp_u_ang = 0.0;
	double tmp_u_len = 0.0;
	double tmp_u_rep = 0.0;
	double tmp_u_adh = 0.0;
	double tmp_u_bnd = 0.0;
	// 直接メンバ変数でreductionをとろうとしてもコンパイルが通らない
	#ifdef OPENMP_CLASS
	#pragma omp parallel num_threads(num_thread) reduction(+:tmp_u_vol, tmp_u_are, tmp_u_ang, tmp_u_len, tmp_u_rep, tmp_u_adh, tmp_u_bnd)
	#endif
	{

		// 体積ポテンシャル
		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < n_c; i++){
			tmp_u_vol += p_c[i]->calcVolumeEnergy();
		}
		// 面積ポテンシャル
		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < n_s; i++){
			tmp_u_are += p_s[i]->calcAreaEnergy();
		}
		// 二面角ポテンシャル
		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < n_l; i++){
			tmp_u_ang += p_l[i]->calcAngleEnergy();
		}
		// 辺長ポテンシャル
		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < n_l; i++){
			tmp_u_len += p_l[i]->calcLengthEnergy();
		}

		// 斥力(排除体積効果)・引力(接着効果)のポテンシャル
		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < n_v; i++){
			Vertex* vp = getVertex(i);
			// 同じgridに属するvertex(IDが自分のIDより小さいものだけ)をIDの大きいものから順にみていく
			int next_vid_sm = next_vid[vp->getId()];	// そのgridの中でvpの次に大きいvertexID
			while(next_vid_sm >= 0){
				Vertex* vq = getVertex(next_vid_sm);

				bool check1 = false;	// vpとvqが同じ細胞か？
				bool check2 = false;	// 同じ辺を共有しているか？
				if(vp->getCell() == vq->getCell()){
					check1 = true;
					for(int j = 0; j < vp->getLineSize(); j++){
						Line* lp = vp->getLine(j);
						if(lp->getVertex(0) == vq || lp->getVertex(1) == vq){
							check2 = true;	break;
						}
					}
				}
				if(!check2){
					tmp_u_rep += calcRepulsiveEnergy_ij(vp, vq);
				}
				if(!check1){
					//tmp_u_adh += calcAdhesiveEnergy_ij(vp, vq);
				}

				next_vid_sm = next_vid[next_vid_sm];
			}
			// gridIDのより小さい隣接gridに属するvertexをIDの大きいものから順にみていく
			int current_grid_id = vp->getGridId();	// vpが属するgridのID
			int tmp = 13*current_grid_id;

			for(int j = 0; j < 13; j++){
				int neighbor_grid_id = nb_list[tmp + j];
				if(neighbor_grid_id >= 0){
					int next_vid_nb = getHeadId(neighbor_grid_id);
					while(next_vid_nb >= 0){
						Vertex* vq = getVertex(next_vid_nb);

						bool check1 = false;	// vpとvqが同じ細胞か？
						bool check2 = false;	// 同じ辺を共有しているか？
						if(vp->getCell() == vq->getCell()){
							check1 = true;
							for(int j = 0; j < vp->getLineSize(); j++){
								Line* lp = vp->getLine(j);
								if(lp->getVertex(0) == vq || lp->getVertex(1) == vq){
									check2 = true;	break;
								}
							}
						}
						if(!check2){
							tmp_u_rep += calcRepulsiveEnergy_ij(vp, vq);
						}
						if(!check1){
							//tmp_u_adh += calcAdhesiveEnergy_ij(vp, vq);
						}

						next_vid_nb = next_vid[next_vid_nb];
					}
				}
			}
		}
/*
		// 境界条件
		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < n_v; i++){
			tmp_u_bnd += p_v[i]->calcBoundaryEnergy();
		}
*/
	}
	u_total = tmp_u_vol + tmp_u_are + tmp_u_ang + tmp_u_len + tmp_u_rep + tmp_u_adh + tmp_u_bnd;
	u_vol = tmp_u_vol;
	u_are = tmp_u_are;
	u_ang = tmp_u_ang;
	u_len = tmp_u_len;
	u_rep = tmp_u_rep;
	u_adh = tmp_u_adh;
	u_bnd = tmp_u_bnd;
}
double Global::getEnergy(){
	return u_total;
}
void Global::calcExtfWork(){
	#ifdef OPENMP_CLASS
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < n_c; i++){
		p_c[i]->calcMovementVector();
	}

	double tmp_work_extf = 0.0;
	#ifdef OPENMP_CLASS
	#pragma omp parallel for num_threads(num_thread) reduction(+:tmp_work_extf)
	#endif
	for(int i = 0; i < n_c; i++){
		tmp_work_extf += (p_c[i]->getExternalForce())*(p_c[i]->getMovementVector());
	}
	work_extf = tmp_work_extf;

	total_work_extf += work_extf;
}
void Global::calcRepulsiveForce_ij(Vertex* vi, Vertex* vj){
	// gridの1辺の長さはdist_rep(adh)_cutoffよりも長くする
	_vec<double> dist_vec = vj->getLocation() - vi->getLocation();	// 相対位置ベクトル
	_vec<double> mean_vec = vi->getNormal() - vj->getNormal();	// 相対法線ベクトル
	double dist = dist_vec.norm();	// viとvjの距離

	// checkがfalseのままなら斥力を計算
	bool check = false;
	if(dist >= dist_rep_cutoff){
		check = true;
	}else if(vi->getCell() == vj->getCell() && dist_vec*mean_vec < 0.){
		// 同じ細胞内では、viとvjの法線ベクトルが向かい合っていなかったら(それらの内積が正だったら)斥力を計算しない
		check = true;
	}

	if(!check){
		double ai = vi->getArea();
		double aj = vj->getArea();
		double ave_a = 0.00325754;	// 点の平均占有面積
		mean_vec *= 0.50;
		mean_vec /= sigma*0.03;	// 規格化

		_vec<double> frc = (-1.0)*ene_repulsive*(1.0 - (dist_vec*mean_vec))*mean_vec*ai*aj/(ave_a*ave_a);
		// 上式は次式をdist_vecで微分したもの: U = 0.50*ene_repulsive*(1 - dist_vec*mean_vec)^2
		vi->addForce(frc, omp_get_thread_num());
		vj->addForce((-1.0)*frc, omp_get_thread_num());

		int vi_cellid = vi->getCell()->getId();
		int vj_cellid = vj->getCell()->getId();
		if(vi_cellid != vj_cellid){
			int room = (n_c-2)*min(vi_cellid, vj_cellid)-min(vi_cellid, vj_cellid)*(min(vi_cellid, vj_cellid)-1)/2+max(vi_cellid, vj_cellid)-1;
			if(cellcell_rep_adh[room] == 0 || cellcell_rep_adh[room] == 2){
				cellcell_rep_adh[room] += 1;
			}
		}

		double vi_surfsize = vi->getSurfaceSize();
		double vj_surfsize = vj->getSurfaceSize();
		double vi_vj_total_energy = 0.50*ene_repulsive*(1 - dist_vec*mean_vec)*(1 - dist_vec*mean_vec);
		double vi_vj_energy_per_surface = vi_vj_total_energy/(vi_surfsize + vj_surfsize);
		for(int i = 0; i < vi_surfsize; i++){
			vi->getSurface(i)->addSurfaceEnergy(vi_vj_energy_per_surface);
		}
		for(int j = 0; j < vj_surfsize; j++){
			vj->getSurface(j)->addSurfaceEnergy(vi_vj_energy_per_surface);
		}
	}
}
double Global::calcRepulsiveEnergy_ij(Vertex* vi, Vertex* vj){
	// gridの1辺の長さはdist_rep(adh)_cutoffよりも長くする
	_vec<double> dist_vec = vj->getLocation() - vi->getLocation();	// 相対位置ベクトル
	_vec<double> mean_vec = vi->getNormal() - vj->getNormal();	// 相対法線ベクトル
	double dist = dist_vec.norm();	// viとvjの距離

	// checkがfalseのままなら斥力を計算
	bool check = false;
	if(dist >= dist_rep_cutoff){
		check = true;
	}else if(vi->getCell() == vj->getCell() && dist_vec*mean_vec < 0.){
		// 同じ細胞内では、viとvjの法線ベクトルが向かい合っていなかったら(それらの内積が正だったら)斥力を計算しない
		check = true;
	}

	double rep_energy = 0.0;
	if(!check){
		double ai = vi->getArea();
		double aj = vj->getArea();
		double ave_a = 0.00325754;	// 点の平均占有面積
		mean_vec *= 0.50;
		mean_vec /= sigma*0.03;	// 規格化

		rep_energy = 0.50*ene_repulsive*(dist_vec*mean_vec - 1.0)*(dist_vec*mean_vec - 1.0)*ai*aj/(ave_a*ave_a);
	}

	return rep_energy;
}
void Global::calcAdhesiveForce_ij(Vertex* vi, Vertex* vj){
	_vec<double> dist_vec = vj->getLocation() - vi->getLocation();	// 相対位置ベクトル
	_vec<double> mean_vec = vi->getNormal() - vj->getNormal();	// 相対法線ベクトル
	double dist = dist_vec.norm();	// viとvjの距離

	// checkがfalseのままなら接着力を計算
	bool check = false;
	if(dist >= dist_adh_cutoff){
		check = true;
	}else if(dist_vec*mean_vec < 0.){
		// viとvjの法線ベクトルが向かい合っていなかったら(それらの内積が正だったら)接着力を計算しない
		check = true;
	}

	if(!check){
		_vec<double> frc;
		Cell* ci = vi->getCell();
		Cell* cj = vj->getCell();
		double ai = vi->getArea();
		double aj = vj->getArea();
		double ave_a = 0.00325754;	// 点の平均占有面積

		double ene_adh;
		if(!(ci->getFlagCentral()) && !(cj->getFlagCentral())){
			ene_adh = ene_adhesion_oo;
		}else if(ci->getFlagCentral() && cj->getFlagCentral()){
			ene_adh = ene_adhesion_ii;
		}else{
			ene_adh = ene_adhesion_io;
		}
		frc = ene_adh*sigma*ai*aj*dist_vec/(dist*ave_a*ave_a);

		vi->addForce(frc, omp_get_thread_num());
		vj->addForce((-1.0)*frc, omp_get_thread_num());

		int vi_cellid = ci->getId();
		int vj_cellid = cj->getId();
		int room = (n_c-2)*min(vi_cellid, vj_cellid)-min(vi_cellid, vj_cellid)*(min(vi_cellid, vj_cellid)-1)/2+max(vi_cellid, vj_cellid)-1;
		if(cellcell_rep_adh[room] == 0 || cellcell_rep_adh[room] == 1){
			cellcell_rep_adh[room] += 2;
		}
	}
}
double Global::calcAdhesiveEnergy_ij(Vertex* vi, Vertex* vj){
	_vec<double> dist_vec = vj->getLocation() - vi->getLocation();	// 相対位置ベクトル
	_vec<double> mean_vec = vi->getNormal() - vj->getNormal();	// 相対法線ベクトル
	double dist = dist_vec.norm();	// viとvjの距離

	// checkがfalseのままなら接着力を計算
	bool check = false;
	if(dist >= dist_adh_cutoff){
		check = true;
	}else if(dist_vec*mean_vec < 0.){
		// viとvjの法線ベクトルが向かい合っていなかったら(それらの内積が正だったら)接着力を計算しない
		check = true;
	}

	double adh_energy = 0.0;
	if(!check){
		Cell* ci = vi->getCell();
		Cell* cj = vj->getCell();
		double ai = vi->getArea();
		double aj = vj->getArea();
		double ave_a = 0.00325754;	// 点の平均占有面積

		double ene_adh;
		if(!(ci->getFlagCentral()) && !(cj->getFlagCentral())){
			ene_adh = ene_adhesion_oo;
		}else if(ci->getFlagCentral() && cj->getFlagCentral()){
			ene_adh = ene_adhesion_ii;
		}else{
			ene_adh = ene_adhesion_io;
		}
		adh_energy = (-1.0)*ene_adh*sigma*(dist - dist_adh_cutoff)*ai*aj/(ave_a*ave_a);
	}

	return adh_energy;
}
void Global::defineGrid(){
	// この関数は1回しか呼び出されない
	// x,y,z方向それぞれのGridの個数を計算する
	grid_size_x = (int)(sys_size_x / grid_length) + 1;
	grid_size_y = (int)(sys_size_y / grid_length) + 1;
	grid_size_z = (int)(sys_size_z / grid_length) + 1;

	cout << grid_size_x*grid_size_y*grid_size_z << " grids (";
	cout << grid_size_x << " * " << grid_size_y << " * " << grid_size_z << ")" << endl;

	int neighbors_size = 13*grid_size_x*grid_size_y*grid_size_z;
	nb_list.assign(neighbors_size, 0);	// 配列nb_listの要素数を変更しそれらの値をすべて0で初期化

	for(int gz = 0; gz < grid_size_z; gz++){
		for(int gy = 0; gy < grid_size_y; gy++){
			for(int gx = 0; gx < grid_size_x; gx++){
				int imap = 13*getGridId(gx, gy, gz);
				// 隣接gridのうち、自分のgridIDよりもIDが小さいgridについて、そのID((27-1)/2 = 13個)を格納
				nb_list[imap] = getGridId(gx-1,	gy-1,	gz-1);
				nb_list[imap+1] = getGridId(gx	,	gy-1,	gz-1);
				nb_list[imap+2] = getGridId(gx+1,	gy-1,	gz-1);

				nb_list[imap+3] = getGridId(gx-1,	gy,		gz-1);
				nb_list[imap+4] = getGridId(gx	,	gy,		gz-1);
				nb_list[imap+5] = getGridId(gx+1,	gy,		gz-1);

				nb_list[imap+6] = getGridId(gx-1,	gy+1,	gz-1);
				nb_list[imap+7] = getGridId(gx	,	gy+1,	gz-1);
				nb_list[imap+8] = getGridId(gx+1,	gy+1,	gz-1);

				nb_list[imap+9] = getGridId(gx-1,	gy-1,	gz);
				nb_list[imap+10] = getGridId(gx	,	gy-1,	gz);
				nb_list[imap+11] = getGridId(gx+1,	gy-1,	gz);

				nb_list[imap+12] = getGridId(gx-1,	gy,		gz);
			}
		}
	}
}
int Global::getGridId(int tmp_i, int tmp_j, int tmp_k){
	if(tmp_i == -1 || tmp_j == -1 || tmp_k == -1 || tmp_i == grid_size_x || tmp_j == grid_size_y || tmp_k == grid_size_z){
		return -1;	// グリッドが存在しないときは-1を返す
	}else{
		int tmp_1 = (tmp_i+grid_size_x)%grid_size_x;
		int tmp_2 = ((tmp_j+grid_size_y)%grid_size_y)*grid_size_x;
		int tmp_3 = ((tmp_k+grid_size_z)%grid_size_z)*grid_size_x*grid_size_y;
		return tmp_1 + tmp_2 + tmp_3;
	}
}
void Global::setupGrid(){
	// 系全体の重心の計算(各細胞の重心の重心)(旧calc/getCentroidSystem())
	sys_centroid = _vec<double>(0., 0., 0.);
	for(int ci = 0; ci < n_c; ci++){
		Cell* cp = getCell(ci);
		sys_centroid += cp->getCentroidCell();
	}
	sys_centroid /= n_c;

	// 全座標を非負にするための移動量の計算(旧calcShift())
	shift = _vec<double>(0., 0., 0.);
	shift.x = 0.5*grid_length*(grid_size_x);
	shift.y = 0.5*grid_length*(grid_size_y);
	shift.z = 0.5*grid_length*(grid_size_z);
	shift -= sys_centroid;

	// vertexとgridの関係づけ(旧generateGrid())
	head_vid.resize(grid_size_x*grid_size_y*grid_size_z);
	#ifdef OPENMP_CLASS
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int i = 0; i < head_vid.size(); i++){
		head_vid[i] = -1;	// すべての要素を-1で初期化
	}
	next_vid.resize(0);
	for(int vi = 0; vi < n_v; vi++){
		Vertex* vp = getVertex(vi);
		_vec<double> copy_loc = vp->getLocation();
		copy_loc += shift;
		copy_loc /= grid_length;
		int x_i = (int)copy_loc.x;
		int y_j = (int)copy_loc.y;
		int z_k = (int)copy_loc.z;
		if(x_i >= grid_size_x || y_j >= grid_size_y || z_k >= grid_size_z){
			cout << "Error: out of grid" << endl;	exit(1);
		}
		// このvertexが属するgridのID	// gridID = i + j*GX + k*GX*GY
		int GRID_ID = x_i + y_j*grid_size_x + z_k*grid_size_x*grid_size_y;
		if(GRID_ID < 0){
			cout << "Error: negative grid ID" << endl;	exit(1);
		}
		vp->setGridId(GRID_ID);
		next_vid.push_back(getHeadId(GRID_ID));
		setHeadId(GRID_ID, vp->getId());	// 上書きされていく
		// updateGeometryのときに、各ステップに固有の next_vid および head_vid が完成する
	}
}
void Global::setHeadId(int a, int b){
	head_vid[a] = b;
}
int Global::getHeadId(int a){
	return head_vid[a];
}
void Global::makeRandomChoice(){
	box_revcell.clear();
	vector<int> chosen;	chosen.clear();

	for(int i = 0; i < n_c; i++){
		if(!(p_c[i]->getFlagCentral())){
			box_revcell.push_back(i);	// 回転領域内のCellのIDを配列boxに入れる
		}else{
			p_c[i]->setBoolRevcell(false);	// v12追加
		}
	}

/*	// ***一定stepおきに一定の割合のCellを抽出する場合***
{
	for(int i = 0; i < n_c; i++){
		p_c[i]->setBoolRevcell(false);	// いったん全Cellのrev_cellをfalseに設定
	}

	// 指定された割合に基づいた抽出Cell数
	int n_need = (int)(box_revcell.size()*0.7);

	// n_need の数だけrandomに抽出
	for(int i = 0; i < n_need; i++){
		int index = (int)(RAND()*box_revcell.size());	// box_revcellの何番目を抽出するかrandomに決定
		chosen.push_back(box_revcell[index]);			// 抽出したrandom番目のCellIDを配列chosenに入れる
		box_revcell.erase(box_revcell.begin() + index);	// 抽出したCellIDをbox_revcell内から消去
	}

	// 抽出したCellに対してrev_cellをfalseからtrueに変更
	// コンソール画面に出力
	cout << chosen.size() << " cells (ID:";
	for(int i = 0; i < chosen.size(); i++){
		int rev_cell_id = chosen[i];
		p_c[rev_cell_id]->setBoolRevcell(true);
		cout << " " << rev_cell_id << ",";
	}
	cout << ") are chosen." << endl;
}*/
	// ***1stepおきに確率的にrev_cellを変更する場合***
{
	for(int i = 0; i < box_revcell.size(); i++){
		Cell* cp = p_c[box_revcell[i]];
		double stay_to_move = (p_stay_to_move*dt_step*0.25)/((4.78*2.0)/120.0);	// 止まっているCellが動き出す確率(0.00～1.00)	// v17にて修正
		double move_to_stay = (p_move_to_stay*dt_step*0.25)/((4.78*2.0)/120.0);	// 動いているCellが止まる確率(0.00～1.00)	// v17にて修正
		if(!(cp->getBoolRevcell())){
			if(RAND() < stay_to_move){
				cp->setBoolRevcell(true);	// stay -> move
			}
		}else if(cp->getBoolRevcell()){
			if(RAND() < move_to_stay){
				cp->setBoolRevcell(false);	// move -> stay
			}
		}
	}
}
}
void Global::countRevCell(){
	n_revcell = 0;
	for(int i = 0; i < n_c; i++){
		if(p_c[i]->getBoolRevcell()){
			n_revcell++;
		}
	}

	total_n_revcell += n_revcell;

	// コンソール画面に出力
/*
	cout << n_revcell << " cells (ID:";
	for(int i = 0; i < n_c; i++){
		if(p_c[i]->getBoolRevcell()){
			cout << " " << p_c[i]->getId() << ",";
		}
	}
	cout << ") are moving." << endl;
*/
}
int Global::getCellCellRepAdh(int room){
	return cellcell_rep_adh[room];
}
void Global::Initialize(){
	cellcell_rep_adh.resize(n_c*(n_c-1)/2);

	#ifdef OPENMP_CLASS
	#pragma omp parallel num_threads(num_thread)
	#endif
	{
		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < cellcell_rep_adh.size(); i++){
			cellcell_rep_adh[i] = 0;
		}

		#ifdef OPENMP_CLASS
		#pragma omp for
		#endif
		for(int i = 0; i < n_s; i++){
			p_s[i]->setSurfaceEnergy(0.0);
		}
	}
}


Vertex::Vertex(int a){
	id = a;
}
void Vertex::pushLine(Line* tmp_l){
	p_l.push_back(tmp_l);
}
Line* Vertex::getLine(int a){
	return p_l[a];
}
void Vertex::pushSurface(Surface* tmp_s){
	p_s.push_back(tmp_s);
}
Surface* Vertex::getSurface(int a){
	return p_s[a];
}
void Vertex::setCell(Cell* cp){
	p_c = cp;
}
Cell* Vertex::getCell(){
	return p_c;
}
void Vertex::setLocation(_vec<double> vec){
	loc = vec;
}
_vec<double> Vertex::getLocation(){
	return loc;
}
void Vertex::copyLocationToLocationT(){
	loc_t = loc;
}
_vec<double> Vertex::getLocationT(){
	return  loc_t;
}
int Vertex::getId(){
	return id;
}
void Vertex::setForce(_vec<double> v){
	frc = v;
}
_vec<double> Vertex::getForce(){
	return frc;
}
void Vertex::addForce(_vec<double> v){
	frc += v;
}
void Vertex::addForce(_vec<double> v, int a){
	frc_omp[a] += v;
}
void Vertex::setForceZero(){
	frc.ZEROS();
}
void Vertex::setForceZero(int a){
	for(int i = 0; i < a; i++){
		frc_omp[i].ZEROS();
	}
}
void Vertex::copyForceRk(int a){
	frc_rk[a] = frc;
}
_vec<double> Vertex::getForceRk(int a){
	return frc_rk[a];
}
void Vertex::sumForceOmp(){
	for(int i = 0; i < num_thread; i++){
		frc += frc_omp[i];
	}
}
void Vertex::setLocation4th(_vec<double> v){
	loc_4th = v;
}
_vec<double> Vertex::getLocation4th(){
	return loc_4th;
}
void Vertex::calcNormalAndAreaVertex(){
	_vec<double> sum_prod_a_n;	sum_prod_a_n.ZEROS();	// 足しこむ
	double sum_surf_area = 0.0;							// 足しこむ
	double max_surf_area = -1.0;						// 更新されていく
	_vec<double> normal_at_max_area;					// 更新されていく

	for(int i = 0; i < getSurfaceSize(); i++){
		Surface* sp = getSurface(i);
		_vec<double> surf_normal = sp->getNormal();
		double surf_area = sp->getArea();

		sum_prod_a_n += surf_area*surf_normal;
		sum_surf_area += surf_area;

		if(surf_area > max_surf_area){
			max_surf_area = surf_area;			// 更新
			normal_at_max_area = surf_normal;	// 更新
		}
	}

	normal = sum_prod_a_n / sum_surf_area;	// 重みづけ平均	//n = Σ(S_i*n_i)/Σ(S_i)
	double norm_normal = normal.norm();		// 法線ベクトルの大きさ
	if(norm_normal > eps){					// eps = 1e-20 <- _parameters.h
		normal /= norm_normal;
	}else{
		normal = normal_at_max_area;
	}

	area = sum_surf_area / 3.0;

}
_vec<double> Vertex::getNormal(){
	return normal;
}
double Vertex::getArea(){
	return area;
}
int Vertex::getLineSize(){
	return p_l.size();
}
int Vertex::getSurfaceSize(){
	return p_s.size();
}
void Vertex::calcExternalForce(){
/*
	Cell* cp = getCell();
	for(int i = 19; i < 38; i++){
		if(cp->getId() == i){
			_vec<double> frc = _vec<double>(0.0, 0.0, -0.2);
			addForce(frc, omp_get_thread_num());
		}
	}
	for(int i = 38; i < 57; i++){
		if(cp->getId() == i){
			_vec<double> frc = _vec<double>(0.0, 0.0, -0.4);
			addForce(frc, omp_get_thread_num());
		}
	}
	for(int i = 57; i < 76; i++){
		if(cp->getId() == i){
			_vec<double> frc = _vec<double>(0.0, 0.0, -0.6);
			addForce(frc, omp_get_thread_num());
		}
	}
*/

// ***2体衝突***
{
	Cell* cp = getCell();	// このvertexが所属するCell
	if(cp->getId() == 0){
		_vec<double> frc = _vec<double>(0.5, 0., 0.);
		addForce(frc, omp_get_thread_num());
	}
	/*
	if(cp->getId() == 1){
		_vec<double> frc = _vec<double>(-0.5, 0., 0.);
		addForce(frc, omp_get_thread_num());
	}
	*/
}
// ***楕円ケプラー運動***
/*{
	Cell* cp = getCell();	// このvertexが所属するCell
	_vec<double> cell_cg = cp->getCentroidCell();	// このvertexが所属するCellの重心

	_vec<double> fpoint = _vec<double>(-owall_a, 0.0, 0.0);	// 焦点位置
	_vec<double> far_point = _vec<double>(owall_a, 0., 0.);
	_vec<double> height_vec = _vec<double>(0.0, 0.0, cell_cg.z);
	double far_dist = (far_point - fpoint).norm();
	double motion_beyond_z = 0.0;

	double phi = atan2(cell_cg.y, cell_cg.x);	// 細胞重心の偏角
	double hantei = (cell_cg.x)*(cell_cg.x)/(iwall_a*iwall_a)+(cell_cg.y)*(cell_cg.y)/(iwall_b*iwall_b);

	// 細胞がiwallよりも外側にあるならば(ifの条件)
	if( hantei > 1.0 ){
		// 接線方向ベクトル
		_vec<double> vec_t = _vec<double>(-owall_a*owall_a*sin(phi), owall_b*owall_b*cos(phi), 0.0);
		vec_t /= vec_t.norm();  // 正規化(単位ベクトル化)

		// v11.3.1, v11.4.1
	//	double rot = (-30.0)*M_PI/180.0;	// この角度だけベクトルを反時計回りに回転させる[rad]
	//	_vec<double> vec_t_rot = _vec<double>(cos(rot)*vec_t.x - sin(rot)*vec_t.y, sin(rot)*vec_t.x + cos(rot)*vec_t.y, 0.0);
	//	vec_t_rot /= vec_t_rot.norm();	// 正規化(単位ベクトル化)

		cp->rot_t = vec_t;  // saigo

		// 法線方向ベクトル(内向き)
		_vec<double> vec_n = _vec<double>(-owall_b*owall_b*cos(phi), -owall_a*owall_a*sin(phi), 0.0);
		//_vec<double> vec_n = (-1.0)*(pp->loc);  // saigo
		vec_n /= vec_n.norm();  // saigo
		cp->rot_n = vec_n;  // saigo

		// 細胞が描く楕円の面積
		//double area = ellipse_area*radi*radi/(e_radi*e_radi);  // saigo
		// 細胞の焦点からの距離(速さはこの距離に依存)
		double dist_f = (cell_cg - (fpoint+height_vec)).norm();  // saigo
		// 細胞の速さ
		//double speed = area*rt_area/dist_f;  // saigo
		dist_f /= far_dist;//規格化

		// distance criterion for elliptic motion
		double dist_c = 0.8;
		double speed = 1. - tanh((dist_f - dist_c)*kappa_d);
		speed *= 0.5;
		speed *= frc_t;

		// 接線方向＋中心方向
		_vec<double> frc = speed*(cp->rot_t) + frc_c*(cp->rot_n);

		//Dorsalward force is applied only if the position-z exceeds motion_beyond_z
		if( cell_cg.z > motion_beyond_z ) frc += _vec<double>(0., 0., frc_z);

		addForce(frc, omp_get_thread_num());
	}
}*/
}
void Vertex::calcBoundaryForce(){
	double phi = atan2(getLocation().y, getLocation().x);	// このvertexの偏角[rad]
	_vec<double> vec_n = _vec<double>(-owall_b*owall_b*cos(phi), -owall_a*owall_a*sin(phi), 0.0);	// 上で求めた偏角をもつ楕円上の点における法線方向(内向き)ベクトル
	vec_n /= vec_n.norm();	// 正規化(単位ベクトル化)

	// Repulsion from ECM wall
	double hantei = (getLocation().x)*(getLocation().x)/(owall_a*owall_a)+(getLocation().y)*(getLocation().y)/(owall_b*owall_b);
	// このvertexが楕円壁の外にあるならば，vec_nの向きに，壁との距離に応じた力を加える
	if(hantei > 1.0){
		double rep_frc_wall = exp(10.0*(hantei - 1.0)) - 1.0;	// 力の大きさ
		_vec<double> frc = rep_frc_wall * vec_n;	// 力ベクトル
		addForce(frc, omp_get_thread_num());
	}
/*
	// 正中面(z=0.0)
	if(getLocation().z < 0.0){
		double rep_frc_floor = exp(10.0*(0.0 - getLocation().z)) - 1.0;
		_vec<double> frc = _vec<double>(0.0, 0.0, rep_frc_floor);
		addForce(frc, omp_get_thread_num());
	}

	// フタ	// 4_layers
	if(getLocation().z > 3.6){
		double rep_frc_ceiling = exp(10.0*(getLocation().z - 3.6)) - 1.0;
		_vec<double> frc = _vec<double>(0.0, 0.0, (-1.0)*rep_frc_ceiling);
		addForce(frc, omp_get_thread_num());
	}

	// フタ	// 1_layer
	if(getLocation().z > 1.0){
		double rep_frc_ceiling = exp(10.0*(getLocation().z - 1.0)) - 1.0;
		_vec<double> frc = _vec<double>(0.0, 0.0, (-1.0)*rep_frc_ceiling);
		addForce(frc, omp_get_thread_num());
	}
*/
}
double Vertex::calcBoundaryEnergy(){
	double bnd_energy = 0.0;

	// Repulsion from ECM wall
	double hantei = (getLocation().x)*(getLocation().x)/(owall_a*owall_a)+(getLocation().y)*(getLocation().y)/(owall_b*owall_b);

	if(hantei > 1.0){
		bnd_energy += 0.1*exp(10.0*(hantei - 1.0))- (hantei - 1.0) - 0.1;
	}
/*
	// 正中面(z=0.0)
	if(getLocation().z < 0.0){
		bnd_energy += 0.1*exp(10.0*(0.0 - getLocation().z)) - (0.0 - getLocation().z) - 0.1;
	}

	// フタ	// 4_layers
	if(getLocation().z > 3.6){
		bnd_energy += 0.1*exp(10.0*(getLocation().z - 3.6)) - (getLocation().z - 3.6) - 0.1;
	}

	// フタ	// 1_layer
	if(getLocation().z > 1.0){
		bnd_energy += 0.1*exp(10.0*(getLocation().z - 1.0)) - (getLocation().z - 1.0) - 0.1;
	}
*/
	return bnd_energy;
}
void Vertex::setGridId(int a){
	grid_id = a;
}
int Vertex::getGridId(){
	return grid_id;
}


void Line::setId(int a){
	id = a;
}
int Line::getId(){
	return id;
}
void Line::setVertex(Vertex* p_v_a, Vertex* p_v_b){
	p_v[0] = p_v_a;
	p_v[1] = p_v_b;
}
Vertex* Line::getVertex(int a){
	if(a < 2){
		return p_v[a];
	} else{
		cout << "Error(Line::getVertex): inadmissible value." << endl;
		exit(1);
	}
}
void Line::setSurface(int a, Surface* p_s_a){
	if(a < 2){
		p_s[a] = p_s_a;
	} else{
		cout << "Error(Line::setSurface): inadmissible value." << endl;
		exit(1);
	}
}
void Line::setSurface(Surface* p_s_a, Surface* p_s_b){
	p_s[0] = p_s_a;
	p_s[1] = p_s_b;
}
Surface* Line::getSurface(int a){
	if(a < 2){
		return p_s[a];
	} else{
		cout << "Error(Line::getSurface): inadmissible value." << endl;
		exit(1);
	}
}
// 初期形状による角度ばね係数の補正。Disney Studioによると、メッシュの初期辺長と面積によって重みづけをする必要があるらしい。
// 今回は、三角形メッシュのサイズがばらついているので、この効果を入れる。
void Line::correctKAngle(){
	// この辺が所属する面が1つしかないなら、この辺は端の辺なので、角度ばねは入らない。
	if(p_s[0] == nullptr || p_s[1] == nullptr){
		k_angle_pos = 0.0;
	} else{
		_vec<double> dist = p_v[1]->getLocation() - p_v[0]->getLocation();
		double dist_sqr = dist.sqr();
		double area_tri = p_s[0]->getArea() + p_s[1]->getArea();
		k_angle_pos = 6.0*k_angle_pre*dist_sqr/area_tri;
	}
	// k_angle_pos = k_angle_pre;	// ***
	// ↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
	// 初期形状による角度ばね係数の補正の効果をなくし、すべて一定値にしたいときは、
	// ***の行のコメントアウトをはずし、この行以外をコメントアウトする。
	// 補正の効果をなくすと、同じ値を設定したとき角度ばねの効果が弱まる。
}
void Line::calcAngleForce(){
	// この辺が所属する面が1つしかないなら、この辺は端の辺なので、角度ばねは入らない。
	if(p_s[0] == nullptr || p_s[1] == nullptr){
		return;
	} else{
		// U = D/2 * 3|e_0(0)|^2/A_i(0) * (2tan(th/2)-0)^2
		// ∇U = D * 3|e_0(0)|^2/A_i(0) * 2 / (1 + cos(th)) * ∇th
		Vertex* p_v_a = nullptr;
		Vertex* p_v_b = nullptr;
		for(int i = 0; i < 3; i++){
			Vertex* tmp_p_v = p_s[0]->getVertex(i);
			if(tmp_p_v != p_v[0] && tmp_p_v != p_v[1]){
				p_v_a = tmp_p_v;
				i = 3;
			}
		}

		for(int i = 0; i < 3; i++){
			Vertex* tmp_p_v = p_s[1]->getVertex(i);
			if(tmp_p_v != p_v[0] && tmp_p_v != p_v[1]){
				p_v_b = tmp_p_v;
				i = 3;
			}
		}

		// 辺の構成の時に、p_s[0]に対してp_l（v:0→1）がCounterClockWise、p_s[1]に対してp_l（v:0→1）がClockWiseになるようにしているはず。

		//				1
		//			  //|\\
		//			 /  |  \
		//		e3	/ a3|a4 \	e4
		//		   / 	|	 \
		//		 |/_ 	|	 _\|
		//		 2_	  e0|	  _3
		//		 |\ 	|	  /|
		//		   \	|	 /
		//		e1	\ a1|a2 /	e2
		//			 \  |  /
		//			  \ | /
		//				0

		// 辺ベクトル
		_vec<double> e_0 = p_v[1]->getLocation() - p_v[0]->getLocation();
		_vec<double> e_1 = p_v_a->getLocation() - p_v[0]->getLocation();
		_vec<double> e_2 = p_v_b->getLocation() - p_v[0]->getLocation();
		_vec<double> e_3 = p_v_a->getLocation() - p_v[1]->getLocation();
		_vec<double> e_4 = p_v_b->getLocation() - p_v[1]->getLocation();

		_vec<double> n_0 = e_0%e_3;
		_vec<double> n_1 = e_4%e_0;
		n_0 = n_0/(n_0.norm() + eps);
		n_1 = n_1/(n_1.norm() + eps);

		double cos_th = n_0*n_1;
		//double tan_hth = sqrt( (1 - cos_th)/(1 + cos_th) );	// 数値誤差とかでわずかにでもcosが1を超えたらnanになる
		double tan_hth = (n_0 - n_1).norm()/(n_0 + n_1).norm();

		// thetaの正負
		_vec<double> relative_loc = p_s[1]->getCentroid() - p_s[0]->getCentroid();
		_vec<double> dn = n_1 - n_0;

		double dU__dth;
		if(relative_loc*dn > 0){
			dU__dth = (-1.0)*k_angle_pos*tan_hth*4.0/(1.0 + cos_th);
		} else{
			dU__dth = k_angle_pos*tan_hth*4.0/(1.0 + cos_th);
		}

		/*
		double dU__dth = (-1.0)*k_angle_pos*tan_hth*2.0/(1.0 + cos_th);
		double dU__dth = k_angle_pos*tan_hth*2.0/(1.0 + cos_th);
		if(dU__dth > 0){
			cout << dU__dth << endl;
		}
		*/

		// 頂点から向かい側の辺に下した垂線の長さ
		double h_01 = 2.0*(p_s[0]->getArea()) / (e_0.norm() + eps);
		double h_02 = 2.0*(p_s[1]->getArea()) / (e_0.norm() + eps);
		double h_1 = 2.0*(p_s[0]->getArea()) / (e_1.norm() + eps);
		double h_2 = 2.0*(p_s[1]->getArea()) / (e_2.norm() + eps);
		double h_3 = 2.0*(p_s[0]->getArea()) / (e_3.norm() + eps);
		double h_4 = 2.0*(p_s[1]->getArea()) / (e_4.norm() + eps);

		// 2面で共有する辺付近の角
		double cos_a1 = e_0*e_1/((e_0.norm()*e_1.norm()) + eps);
		double cos_a2 = e_0*e_2/((e_0.norm()*e_2.norm()) + eps);
		double cos_a3 = (-1.0)*e_0*e_3/((e_0.norm()*e_3.norm()) + eps);
		double cos_a4 = (-1.0)*e_0*e_4/((e_0.norm()*e_4.norm()) + eps);

		_vec<double> frc_r0 = cos_a3/(h_3 + eps)*n_0 + cos_a4/(h_4 + eps)*n_1;
		_vec<double> frc_r1 = cos_a1/(h_1 + eps)*n_0 + cos_a2/(h_2 + eps)*n_1;
		_vec<double> frc_r2 = -1.0/(h_01 + eps)*n_0;
		_vec<double> frc_r3 = -1.0/(h_02 + eps)*n_1;

		frc_r0 *= dU__dth;
		frc_r1 *= dU__dth;
		frc_r2 *= dU__dth;
		frc_r3 *= dU__dth;

		p_v[0]->addForce(frc_r0, omp_get_thread_num());
		p_v[1]->addForce(frc_r1, omp_get_thread_num());
		p_v_a->addForce(frc_r2, omp_get_thread_num());
		p_v_b->addForce(frc_r3, omp_get_thread_num());
	}
}
double Line::calcAngleEnergy(){
	if(p_s[0] == nullptr || p_s[1] == nullptr){
		return 0.0;
	} else{
		// U = D/2 * 3|e_0(0)|^2/A_i(0) * (2tan(th/2)-0)^2
		// ∇U = D * 3|e_0(0)|^2/A_i(0) * 2 / (1 + cos(th)) * ∇th
		Vertex* p_v_a;
		Vertex* p_v_b;
		for(int i = 0; i < 3; i++){
			Vertex* tmp_p_v = p_s[0]->getVertex(i);
			if(tmp_p_v != p_v[0] && tmp_p_v != p_v[1]){
				p_v_a = tmp_p_v;
				i = 3;
			}
		}
		for(int i = 0; i < 3; i++){
			Vertex* tmp_p_v = p_s[1]->getVertex(i);
			if(tmp_p_v != p_v[0] && tmp_p_v != p_v[1]){
				p_v_b = tmp_p_v;
				i = 3;
			}
		}

		// 辺の構成の時に、p_s[0]に対してp_l（v:0→1）がCounterClockWise、p_s[1]に対してp_l（v:0→1）がClockWiseになるようにしているはず。

		//				1
		//			  //|\\
		//			 /  |  \
		//		e3	/ a3|a4 \	e4
		//		   / 	|	 \
		//		 |/_ 	|	 _\|
		//		 2_	  e0|	  _3
		//		 |\ 	|	  /|
		//		   \	|	 /
		//		e1	\ a1|a2 /	e2
		//			 \  |  /
		//			  \ | /
		//				0

		// 辺ベクトル
		_vec<double> e_0 = p_v[1]->getLocation() - p_v[0]->getLocation();
		//_vec<double> e_1 = p_v_a->getLocation() - p_v[0]->getLocation();
		//_vec<double> e_2 = p_v_b->getLocation() - p_v[0]->getLocation();
		_vec<double> e_3 = p_v_a->getLocation() - p_v[1]->getLocation();
		_vec<double> e_4 = p_v_b->getLocation() - p_v[1]->getLocation();

		_vec<double> n_0 = e_0%e_3;
		_vec<double> n_1 = e_4%e_0;
		n_0 = n_0/(n_0.norm() + eps);
		n_1 = n_1/(n_1.norm() + eps);

		double cos_th = n_0*n_1;
		//double tan_hth = sqrt( (1 - cos_th)/(1 + cos_th) );	// 数値誤差とかでわずかにでもcosが1を超えたらnanになる
		double tan_hth = (n_0 - n_1).norm()/(n_0 + n_1).norm();

		return 2.0*k_angle_pos*tan_hth*tan_hth;	//0.5→2.0 phi(theta) = 2tan(theta/2)なので。
	}
}
void Line::setLen0(double a){
	len_0 = a;
}
double Line::getLen0(){
	return len_0;
}
void Line::copyLenToLen0(){
	len_0 = len;
}
void Line::calcLength(){
	len = (p_v[1]->getLocation() - p_v[0]->getLocation()).norm();
}
double Line::getLength(){
	return len;
}
void Line::calcLengthForce(){
	_vec<double> vec_len = p_v[1]->getLocation() - p_v[0]->getLocation();
	double du__dl = (-1.0)*k_length*(len/(len_0+eps) - 1.0)/(len_0+eps);
	//double du__dl = (-1.0)*k_length*(len/(len_0+eps) - 1.0)*0.0638452/((len_0+eps)*(len_0+eps));
	_vec<double> frc = vec_len/len;
	frc *= (-1.0)*du__dl;

	p_v[0]->addForce(frc, omp_get_thread_num());
	p_v[1]->addForce((-1.0)*frc, omp_get_thread_num());
}
double Line::calcLengthEnergy(){
	return 0.5*k_length*(len/(len_0+eps) - 1.0)*(len/(len_0+eps) - 1.0);
	//return 0.5*k_length*(len/(len_0+eps) - 1.0)*(len/(len_0+eps) - 1.0)*0.0638452/(len_0+eps);
}


void Surface::setId(int a){
	id = a;
}
int Surface::getId(){
	return id;
}
void Surface::setVertex(Vertex* p_v_a, Vertex* p_v_b, Vertex* p_v_c){
	p_v[0] = p_v_a;
	p_v[1] = p_v_b;
	p_v[2] = p_v_c;
}
Vertex* Surface::getVertex(int a){
	if(a < 3){
		return p_v[a];
	} else{
		cout << "Error(Surface::getVertex): inadmissible value." << endl;
		exit(1);
	}
}
void Surface::setLine(Line* p_l_a, Line* p_l_b, Line* p_l_c){
	p_l[0] = p_l_a;
	p_l[1] = p_l_b;
	p_l[2] = p_l_c;
}
Line* Surface::getLine(int a){
	if(a < 3){
		return p_l[a];
	} else{
		cout << "Error(Surface::getLine): inadmissible value." << endl;
		exit(1);
	}
}
void Surface::calcCentroidTriangle(){
	centroid = (p_v[0]->getLocation() + p_v[1]->getLocation() + p_v[2]->getLocation())/3.0;
}
_vec<double> Surface::getCentroid(){
	return centroid;
}
void Surface::calcNormalVector(){
	_vec<double> r_01 = p_v[1]->getLocation() - p_v[0]->getLocation();
	_vec<double> r_02 = p_v[2]->getLocation() - p_v[0]->getLocation();
	_vec<double> cross_prod = r_01%r_02;
	normal = cross_prod/(cross_prod.norm() + eps);
}
_vec<double> Surface::getNormal(){
	return normal;
}
void Surface::setArea0(double a){
	area_0 = a;
}
double Surface::getArea0(){
	return area_0;
}
void Surface::setArea(double a){
	area = a;
}
void Surface::calcNormalAndAreaTriangle(){
	_vec<double> r_01 = p_v[1]->getLocation() - p_v[0]->getLocation();
	_vec<double> r_02 = p_v[2]->getLocation() - p_v[0]->getLocation();
	_vec<double> cross_prod = r_01%r_02;

	area = 0.5*(cross_prod.norm());
	normal = cross_prod/(cross_prod.norm() + eps);
}
double Surface::getArea(){
	return area;
}
void Surface::copyAreaToArea0(){
	area_0 = area;
}
void Surface::calcAreaForce(){
	_vec<double> r_01 = p_v[1]->getLocation() - p_v[0]->getLocation();
	_vec<double> r_02 = p_v[2]->getLocation() - p_v[0]->getLocation();
	_vec<double> cross_prod = r_01%r_02;
	double norm_cp = cross_prod.norm();

	// 面積の各頂点に対する勾配。手計算で出てくる結果を計算している
	_vec<double> frc_r0 = 0.5*cross_prod%(r_02-r_01)/(norm_cp+eps);
	_vec<double> frc_r1 = 0.5*r_02%cross_prod/(norm_cp+eps);
	_vec<double> frc_r2 = 0.5*cross_prod%r_01/(norm_cp+eps);

	// 勾配に乗ずる係数。dU/dS = k(S/S_0 - 1)/S_0
	double du__ds = k_area*(0.5*norm_cp/(area_0+eps) - 1.0)/(area_0+eps);
	//double du__ds = k_area*(0.5*norm_cp/(area_0+eps) - 1.0)*0.00163216/((area_0+eps)*(area_0+eps));
	frc_r0 *= -1.0*du__ds;
	frc_r1 *= -1.0*du__ds;
	frc_r2 *= -1.0*du__ds;

	p_v[0]->addForce(frc_r0, omp_get_thread_num());
	p_v[1]->addForce(frc_r1, omp_get_thread_num());
	p_v[2]->addForce(frc_r2, omp_get_thread_num());
}
double Surface::calcAreaEnergy(){
	return 0.5*k_area*(area/(area_0+eps) - 1.0)*(area/(area_0+eps) - 1.0);
	//return 0.5*k_area*(area/(area_0+eps) - 1.0)*(area/(area_0+eps) - 1.0)*0.00163216/(area_0+eps);
}
void Surface::setSurfaceEnergy(double tmp_energy){
	energy = tmp_energy;
}
void Surface::addSurfaceEnergy(double tmp_energy){
	energy += tmp_energy;
}
double Surface::getSurfaceEnergy(){
	return energy;
}


void Cell::setId(int a){
	id = a;
}
int Cell::getId(){
	return id;
}
void Cell::setVertex(vector<Vertex*> vec){
	for(int i = 0; i < vec.size(); i++){
		p_v.push_back(vec[i]);
	}
}
Vertex* Cell::getVertex(int a){
	if(a < getVertexSize()){
		return p_v[a];
	} else{
		cout << "Error(Cell::getVertex): inadmissible value." << endl;
		exit(1);
	}
}
void Cell::setLine(vector<Line*> vec){
	for(int i = 0; i < vec.size(); i++){
		p_l.push_back(vec[i]);
	}
}
Line* Cell::getLine(int a){
	if(a < getLineSize()){
		return p_l[a];
	} else{
		cout << "Error(Cell::getLine): inadmissible value." << endl;
		exit(1);
	}
}
void Cell::setSurface(vector<Surface*> vec){
	for(int i = 0; i < vec.size(); i++){
		p_s.push_back(vec[i]);
	}
}
Surface* Cell::getSurface(int a){
	if(a < getSurfaceSize()){
		return p_s[a];
	} else{
		cout << "Error(Cell::getSurface): inadmissible value." << endl;
		exit(1);
	}
}
void Cell::calcCentroidCell(){
	// 細胞を構成する全頂点座標の平均
	_vec<double> tmp_centroid = _vec<double>(0.0, 0.0, 0.0);
	for(int i = 0; i < getVertexSize(); i++){
		tmp_centroid += getVertex(i)->getLocation();
	}
	centroid = _vec<double>(tmp_centroid.x / getVertexSize(), tmp_centroid.y / getVertexSize(), tmp_centroid.z / getVertexSize());
}
_vec<double> Cell::getCentroidCell(){
	return centroid;
}
void Cell::storeCentroidCell(){
	prev_centroid = centroid;
}
_vec<double> Cell::getStoredCentroidCell(){
	return prev_centroid;
}
void Cell::calcMovementVector(){
	centroid_move = centroid - prev_centroid;
}
_vec<double> Cell::getMovementVector(){
	return centroid_move;
}
void Cell::judgeCentralOrNot(){
	double judgement = ((getCentroidCell().x)*(getCentroidCell().x)/(iwall_a*iwall_a)) + ((getCentroidCell().y)*(getCentroidCell().y)/(iwall_b*iwall_b));
	if(judgement < 1.0){
		flag_central = true;	// central
	}else{
		flag_central = false;	// not central
	}
}
bool Cell::getFlagCentral(){
	return flag_central;
}
void Cell::setBoolRevcell(bool a){
	rev_cell = a;
}
bool Cell::getBoolRevcell(){
	return rev_cell;
}
void Cell::setVolume0(double a){
	volume_0 = a;
}
double Cell::getVolume0(){
	return volume_0;
}
void Cell::setVolume(double a){
	volume = a;
}
void Cell::calcVolumeCell(){
	// 原点(0.,0.,0.)と各面とでできる四面体の符号つき体積を足し合わせて求める
	double tmp_volume = 0.0;
	#ifdef OPENMP_CLASS
	#pragma omp parallel for num_threads(num_thread) reduction(+:tmp_volume)
	#endif
	for(int si = 0; si < getSurfaceSize(); si++){
		_vec<double> r_o0 = p_s[si]->getVertex(0)->getLocation();
		_vec<double> r_o1 = p_s[si]->getVertex(1)->getLocation();
		_vec<double> r_o2 = p_s[si]->getVertex(2)->getLocation();
		double tri_prod = r_o0*(r_o1%r_o2);
		double part_volume = tri_prod/6.0;
		tmp_volume += part_volume;
	}
	volume = tmp_volume;
}
double Cell::getVolume(){
	return volume;
}
void Cell::copyVolumeToVolume0(){
	volume_0 = volume;
}
void Cell::calcVolumeForce(){
	#ifdef OPENMP_CLASS
	#pragma omp parallel for num_threads(num_thread)
	#endif
	for(int si = 0; si < getSurfaceSize(); si++){
		// ここでは重心にしているが別にどこでもよい。
		_vec<double> r_g0 = p_s[si]->getVertex(0)->getLocation() - getCentroidCell();
		_vec<double> r_g1 = p_s[si]->getVertex(1)->getLocation() - getCentroidCell();
		_vec<double> r_g2 = p_s[si]->getVertex(2)->getLocation() - getCentroidCell();

		double tri_prod = r_g0*(r_g1%r_g2);	// 重心から面を眺めたときv0v1v2がClockWiseなら正

		//if(tri_prod <= 0.0){
		//	cout << "inside out" << endl;	// 重心から面を眺めたときv0v1v2がCounterClockWise
		//}

		// 体積の各頂点に対する勾配。手計算で得られる式
		_vec<double> frc_cg = (tri_prod*((r_g0-r_g1)%(r_g2-r_g1)))/(6.0*fabs(tri_prod)+eps);
		_vec<double> frc_r0 = (tri_prod*(r_g1%r_g2))/(6.0*fabs(tri_prod)+eps);
		_vec<double> frc_r1 = (tri_prod*(r_g2%r_g0))/(6.0*fabs(tri_prod)+eps);
		_vec<double> frc_r2 = (tri_prod*(r_g0%r_g1))/(6.0*fabs(tri_prod)+eps);

		// 勾配に乗ずる係数。dU/dV = k(V/V_0 - 1)/V_0
		double du__dv = k_volume*(volume/(volume_0+eps) - 1.0)/(volume_0+eps);

		//frc_cg *= (-1.0)*du__dv;
		frc_r0 *= (-1.0)*du__dv;
		frc_r1 *= (-1.0)*du__dv;
		frc_r2 *= (-1.0)*du__dv;

		if(tri_prod <= 0.0){
			//frc_cg *= -1.0;
			frc_r0 *= -1.0;
			frc_r1 *= -1.0;
			frc_r2 *= -1.0;
		}

		p_s[si]->getVertex(0)->addForce(frc_r0, omp_get_thread_num());
		p_s[si]->getVertex(1)->addForce(frc_r1, omp_get_thread_num());
		p_s[si]->getVertex(2)->addForce(frc_r2, omp_get_thread_num());
	}
}
double Cell::calcVolumeEnergy(){
	return 0.5*k_volume*(volume/(volume_0+eps) - 1.0)*(volume/(volume_0+eps) - 1.0);
}
void Cell::calcInternalForce(){
	int cell_p_id = getId();	// このCellのID
	_vec<double> cell_p_cg = getCentroidCell();	// このCellの重心座標
	_vec<double> z_vec = _vec<double>(0., 0., 0.1);	// 外積用z方向ベクトル
/*
	// ========== 2細胞間の内力(互いを蹴って進む) ==========
	// このCell(p)のIDより小さいIDのCell(q)をpairとする
	for(int cell_q_id = 0; cell_q_id < cell_p_id; cell_q_id++){
		Cell* cq = p_g->getCell(cell_q_id);
		_vec<double> cell_q_cg = cq->getCentroidCell();

		_vec<double> cg_vec = cell_q_cg - cell_p_cg;	// Cellの重心同士をむすぶベクトル(CellIDの大→小の向き)
		double cg_dist = cg_vec.norm();	// Cellの重心間の距離
		cg_dist /= owall_a;	// 規格化

		_vec<double> frc_p = cg_vec % z_vec;
		frc_p /= cg_dist * cg_dist;
		_vec<double> frc_q = z_vec % cg_vec;
		frc_q /= cg_dist * cg_dist;

		for(int p = 0; p < getVertexSize(); p++){
			p_v[p]->addForce(frc_p, omp_get_thread_num());	// このCell(p)に属するvertexにかかる力
		}
		for(int q = 0; q < cq->getVertexSize(); q++){
			cq->getVertex(q)->addForce(frc_q, omp_get_thread_num());	// pairのCell(q)に属するvertexにかかる力
		}
	}
*/
	// ==========壁を蹴る力 ==========
	int counter = 0;	// 壁から出ているvertexの個数
	double outvtx_x = 0.0;	// はみ出たvertexのx座標を足しこむ
	double outvtx_y = 0.0;	// はみ出たvertexのy座標を足しこむ

	for(int i = 0; i < getVertexSize(); i++){
		Vertex* vp = getVertex(i);	// vp = このCellのi番目のvertex
		// vpが楕円壁の外部にあるかどうかの判定
		double hantei = ((vp->getLocation().x)*(vp->getLocation().x))/(owall_a*owall_a) + ((vp->getLocation().y)*(vp->getLocation().y))/(owall_b*owall_b);
		if(hantei > 1.0){	// vpが楕円壁の外部にあるならば
			counter += 1;
			outvtx_x += vp->getLocation().x;
			outvtx_y += vp->getLocation().y;
		}
	}
	if(counter > 0){	// このCellに属するvertexの少なくとも1つが楕円壁の外部にあるならば
		// 楕円壁の外部にあるvertex群の平均座標を求める
		outvtx_x /= counter;
		outvtx_y /= counter;

		double phi = atan2(outvtx_y, outvtx_x);	// 偏角
		//phi -= (5.0)*M_PI/180.0;	// 何度前の情報を参照するか

		_vec<double> vec_t = _vec<double>(-owall_a*owall_a*sin(phi), owall_b*owall_b*cos(phi), 0.0);	// 上で求めた偏角をもつ楕円上の点における接線方向(反時計方向)ベクトル
		vec_t /= vec_t.norm();	// 正規化(単位ベクトル化)

		// 接線方向よりやや外側方向にする
		double rot = (-30.0)*M_PI/180.0;	// この角度だけベクトルを反時計回りに回転させる[rad]
		_vec<double> vec_t_rot = _vec<double>(cos(rot)*vec_t.x - sin(rot)*vec_t.y, sin(rot)*vec_t.x + cos(rot)*vec_t.y, 0.0);
		vec_t_rot /= vec_t_rot.norm();	// 正規化(単位ベクトル化)
		_vec<double> frc = 0.5*vec_t_rot;	// 力の大きさの変更

		for(int i = 0; i < getVertexSize(); i++){
			p_v[i]->addForce(frc, omp_get_thread_num());	// このCellに属するvertexにかける力
		}
	}
}
void Cell::calcCellMigrationForce(){
	// ***楕円回転運動(力一定)***
{
	// 細胞重心位置の偏角(-π ～ π)
	double phi = atan2(getCentroidCell().y, getCentroidCell().x);
	// 偏角がphiであるような楕円上の点の，z軸からの距離
	//double o_dist = sqrt(1.0/ ( (cos(phi)*cos(phi))/(owall_a*owall_a) + (sin(phi)*sin(phi))/(owall_b*owall_b) ));
	// 細胞重心の，z軸からの距離
	//double r_dist = sqrt((getCentroidCell().x)*(getCentroidCell().x) + (getCentroidCell().y)*(getCentroidCell().y));

	// 接線方向ベクトル
	_vec<double> vec_t = _vec<double>(-owall_a*owall_a*sin(phi), owall_b*owall_b*cos(phi), 0.0);
	vec_t /= vec_t.norm();  // 正規化(単位ベクトル化)
	rot_t = vec_t;

	// 法線方向ベクトル(内向き)
	//_vec<double> vec_n = _vec<double>(-owall_b*owall_b*cos(phi), -owall_a*owall_a*sin(phi), 0.0);
	//vec_n /= vec_n.norm();  // 正規化(単位ベクトル化)
	//rot_n = vec_n;

	//double speed = frc_t * owall_a / 2.8;	// 楕円の長短の比率を変えていなければ．(double speed = frc_t * owall_b / 2.2; でも同じ)
	double speed = frc_t;	// 以前のシミュレーションと同じ値とする

	// 接線方向のみ
	//_vec<double> frc = ((speed*r_dist)/o_dist)*rot_t;
	_vec<double> frc = speed*rot_t;
	extf = frc*((double)getVertexSize());	// v17

	// 接線方向＋中心方向
	//_vec<double> frc = speed*rot_t + frc_c*rot_n;

	for(int i = 0; i < getVertexSize(); i++){
		p_v[i]->addForce(frc, omp_get_thread_num());	// このCellに属するvertexにかける力
	}
}
}
void Cell::setExternalForceZero(){
	extf = _vec<double>(0., 0., 0.);
}
_vec<double> Cell::getExternalForce(){
	return extf;
}
int Cell::getVertexSize(){
	return p_v.size();
}
int Cell::getLineSize(){
	return p_l.size();
}
int Cell::getSurfaceSize(){
	return p_s.size();
}
