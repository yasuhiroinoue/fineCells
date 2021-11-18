/*
// For contact information
// Yasuhiro Inoue
// inoue.yasuhiro.4n@kyoto-u.ac.jp
*/

#include <iostream>	// 入出力ライブラリ
#include <cmath>
#include <cstdlib>
#include <vector>
#include <fstream>	// ファイル入出力
#include <omp.h>
#include <string.h>		// 追加
#include <algorithm>	// std::sort, std::unique
#include "vec.h"
#include "_class.h"
#include "_parameters.h"
#include "ODEsolver.h"
#include "random.h"		// 追加

#define OPENMP_MAIN

using namespace std;

Global *p_g;

void readSIG(){
	cout << "Reading SIG File...";

	ifstream fin("two_spheres_distant_2.sig");
	if(!fin){
		cout << "Error: cannot open input sig file." << endl;
		exit(1);
	}

	char dummy[100];
	int num_vertices;
	int num_surfaces;
	int num_lines;
	int num_cells;

	fin >> dummy;
	fin >> num_vertices;
	fin >> num_surfaces;
	fin >> num_lines;
	fin >> num_cells;

	p_g->n_v = num_vertices;
	//p_g->n_l = num_lines;	// sigファイルだとどうせゼロなのでコメントアウト
	p_g->n_s = num_surfaces;
	p_g->n_c = num_cells;

	// 頂点(生成) <-- 位置
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* tmp_v = new Vertex(i);	// ポインタで格納するのでnewで確保しないと消える
		_vec<double> tmp_vec;
		fin >> tmp_vec.x >> tmp_vec.y >> tmp_vec.z;
		tmp_v->setLocation(tmp_vec);
		p_g->pushVertex(tmp_v);
	}

	// 面(生成) <-- 3頂点
	// 頂点の順番は0 1 2の循環から法線ベクトルを作った時にそれが外側を向くように
	for(int i = 0; i < p_g->n_s; i++){
		int tmp;
		fin >> tmp;	// 最初の数字は面の構成点の数。ここでは全部三角形であるので無視する

		Surface* tmp_s = new Surface();
		int tmp_int[3];
		fin >> tmp_int[0] >> tmp_int[1] >> tmp_int[2];
		Vertex* p_v_a = p_g->getVertex(tmp_int[0]);
		Vertex* p_v_b = p_g->getVertex(tmp_int[1]);
		Vertex* p_v_c = p_g->getVertex(tmp_int[2]);
		tmp_s->setVertex(p_v_a, p_v_b, p_v_c);
		tmp_s->setId(i);
		p_g->pushSurface(tmp_s);

		// 3頂点 <-- 面
		p_v_a->pushSurface(tmp_s);
		p_v_b->pushSurface(tmp_s);
		p_v_c->pushSurface(tmp_s);
	}

	// 細胞 <-- n面
	for(int i = 0; i < p_g->n_c; i++){
		vector<Surface*> s_list;
		int num;
		fin >> num; //最初の数字は細胞の構成面の数。
		Cell* tmp_c = new Cell();
		int tmp_int[num];
		for(int j = 0; j < num; j++){
			fin >> tmp_int[j];
		}
		for(int j = 0; j < num; j++){
			Surface* p_s_tmp = p_g->getSurface(tmp_int[j]);
			s_list.push_back(p_s_tmp);
		}
		tmp_c->setSurface(s_list);
		tmp_c->setId(i);
		p_g->pushCell(tmp_c);
	}

	#ifdef OPENMP_MAIN
	#pragma omp parallel num_threads(num_thread)
	#endif
	{
		#ifdef OPENMP_MAIN
		#pragma omp for
		#endif
		// 細胞 <-- n頂点
		for(int ci = 0; ci < p_g->n_c; ci++){
			vector<Vertex*> vertex_list;
			Cell* cp = p_g->getCell(ci);
			for(int si = 0; si < cp->getSurfaceSize(); si++){
				Surface* sp = cp->getSurface(si);
				for(int vi = 0; vi < 3; vi++){
					Vertex* vp = sp->getVertex(vi);
					vertex_list.push_back(vp);
				}
			}
			// ここまでで，ci番目の細胞に属するすべての面について頂点を数え上げ，そのリストがvertex_listに入っている
			// 次の2行はvertex_listの中の頂点の重複をなくすアルゴリズム
			sort(begin(vertex_list), end(vertex_list));
			vertex_list.erase(unique(begin(vertex_list), end(vertex_list)), end(vertex_list));
			// c++11では，～.begin()よりもbegin(～)が推奨されている
			// 以上でci番目の細胞に属する頂点の重複なしリストvertex_listが完成
			cp->setVertex(vertex_list);	// このcpがもつ情報：細胞ciにはvertex_list内にある頂点が属している
		}

		#ifdef OPENMP_MAIN
		#pragma omp for
		#endif
		// n頂点 <-- 細胞
		for(int ci = 0; ci < p_g->n_c; ci++){
			Cell* cp = p_g->getCell(ci);
			for(int vi = 0; vi < cp->getVertexSize(); vi++){
				Vertex* vp = cp->getVertex(vi);
				vp->setCell(cp);
			}
		}
	}
	cout << " done." << endl;
}

// 辺を生成し、面と頂点の情報を入れる
void reconstLines(){
	cout << "Reconstructing lines..." << endl;

	for(int i = 0; i < p_g->n_s; i++){
		if(i%10000 == 0) cout << "1... " << i << "/" << p_g->n_s << endl;
		Surface* p_s = p_g->getSurface(i);
		Vertex* p_v_a = p_s->getVertex(0);
		Vertex* p_v_b = p_s->getVertex(1);
		Vertex* p_v_c = p_s->getVertex(2);

		Line* tmp_l_a = new Line();
		Line* tmp_l_b = new Line();
		Line* tmp_l_c = new Line();
		tmp_l_a->setVertex(p_v_a, p_v_b);
		tmp_l_b->setVertex(p_v_b, p_v_c);
		tmp_l_c->setVertex(p_v_c, p_v_a);

		bool flag_a = true;
		bool flag_b = true;
		bool flag_c = true;
		Line* p_l_to_s[3] = {tmp_l_a, tmp_l_b, tmp_l_c};

		#ifdef OPENMP_MAIN
		#pragma omp parallel for num_threads(num_thread)
		#endif
		// 最初n_l=0なのでこのforの中身は1巡目は実行されない。後ろでn_l++していくのでそれ以降は(2巡目はn_l=3として)実行される
		for(int i = 0; i < p_g->n_l; i++){
			Line* pre_p_l = p_g->getLine(i);

			// 隣接する2面を考えているので、
			// それらの境界となる辺に属する頂点を参照する順番は、(i-1)番目の面とi番目の面とで逆である。
			// なので、逆にならないおかしい場合をelse ifではじいて(エラーが出るようにして)いる。
			if((pre_p_l->getVertex(0) == tmp_l_a->getVertex(1))
			&& pre_p_l->getVertex(1) == tmp_l_a->getVertex(0)){	
				// 既に登録されている辺(pre_p_l)が新規のi番目の面のある1つの辺(tmp_l_a)と一致するならば(ifの条件)、
				// この辺が属するもう1つの面(Line::p_s[1])はi番目の面(p_s=p_g->getSurface(i))ということになる。
				pre_p_l->setSurface(1, p_s);
				// この辺は既にglobalのp_lに登録されている辺なので今回は登録しない。
				flag_a = false;	// 破棄
				// 新規のi番目の面に属する辺のリストを一部更新する。
				p_l_to_s[0] = pre_p_l;
				delete tmp_l_a;
			} else if((pre_p_l->getVertex(0) == tmp_l_a->getVertex(0))
					&& pre_p_l->getVertex(1) == tmp_l_a->getVertex(1)){
				// 破棄だけど、面の頂点が反時計回りに格納されていることを考えると、多分入力ファイルがおかしい。
				cout << "Warning: Check input file." << endl;
				flag_a = false;	// 破棄
				delete tmp_l_a;
			}

			if((pre_p_l->getVertex(0) == tmp_l_b->getVertex(1))
			&& pre_p_l->getVertex(1) == tmp_l_b->getVertex(0)){
				// 既に登録されている辺(pre_p_l)が新規のi番目の面のある1つの辺(tmp_l_b)と一致するならば(ifの条件)、
				// この辺が属するもう1つの面(Line::p_s[1])はi番目の面(p_s=p_g->getSurface(i))ということになる。
				pre_p_l->setSurface(1, p_s);
				flag_b = false;	// 破棄
				// tmp_l_b = pre_p_l;
				p_l_to_s[1] = pre_p_l;
				delete tmp_l_b;
			} else if((pre_p_l->getVertex(0) == tmp_l_b->getVertex(0))
					&& pre_p_l->getVertex(1) == tmp_l_b->getVertex(1)){
				// 破棄だけど、面の頂点が反時計回りに格納されていることを考えると、多分入力ファイルがおかしい。
				cout << "Warning: Check input file." << endl;
				flag_b = false;	// 破棄
				delete tmp_l_b;
			}

			if((pre_p_l->getVertex(0) == tmp_l_c->getVertex(1))
			&& pre_p_l->getVertex(1) == tmp_l_c->getVertex(0)){
				// 既に登録されている辺(pre_p_l)が新規のi番目の面のある1つの辺(tmp_l_c)と一致するならば(ifの条件)、
				// この辺が属するもう1つの面(Line::p_s[1])はi番目の面(p_s=p_g->getSurface(i))ということになる。
				pre_p_l->setSurface(1, p_s);
				flag_c = false;	// 破棄
				// tmp_l_c = pre_p_l;
				p_l_to_s[2] = pre_p_l;
				delete tmp_l_c;
			} else if((pre_p_l->getVertex(0) == tmp_l_c->getVertex(0))
					&& pre_p_l->getVertex(1) == tmp_l_c->getVertex(1)){
				// 破棄だけど、面の頂点が反時計回りに格納されていることを考えると、多分入力ファイルがおかしい。
				cout << "Warning: Check input file." << endl;
				flag_c = false;	// 破棄
				delete tmp_l_c;
			}
		}
		// 面に属する辺のリストp_l_to_sが更新されたので、その辺のリストをglobalのp_sに登録する。
		p_s->setLine(p_l_to_s[0], p_l_to_s[1], p_l_to_s[2]);

		// ()の中がtrueならば実行
		// 既にカウントしている辺については上でfalseに切り替えてあるので、2重カウントはされない。
		if(flag_a){
			// もし新しく登録された辺だったら(ifの条件)、
			// その辺が属する0番目の面(Line::p_s[0])はi番目の面(p_s=p_g->getSurface(i))ということになる。
			// その辺が属する1番目の面(Line::p_s[1])については、
			// それが端辺でなければ、何番目の面か知らないがいずれ上で定まる。
			tmp_l_a->setSurface(0, p_s);
			tmp_l_a->setId(p_g->n_l);
			// globalのp_lに登録
			p_g->pushLine(tmp_l_a);
			// 逆に、頂点に辺をset
			p_v_a->pushLine(tmp_l_a);
			p_v_b->pushLine(tmp_l_a);
			// 辺の総数を1増やす
			p_g->n_l++;
		}
		if(flag_b){
			tmp_l_b->setSurface(0, p_s);
			tmp_l_b->setId(p_g->n_l);
			// globalのp_lに登録
			p_g->pushLine(tmp_l_b);
			// 逆に、頂点に辺をset
			p_v_b->pushLine(tmp_l_b);
			p_v_c->pushLine(tmp_l_b);
			// 辺の総数を1増やす
			p_g->n_l++;
		}
		if(flag_c){
			tmp_l_c->setSurface(0, p_s);
			tmp_l_c->setId(p_g->n_l);
			// globalのp_lに登録
			p_g->pushLine(tmp_l_c);
			// 逆に、頂点に辺をset
			p_v_c->pushLine(tmp_l_c);
			p_v_a->pushLine(tmp_l_c);
			// 辺の総数を1増やす
			p_g->n_l++;
		}
	}

	// check duplication of lines.
	for(int i = 0; i < p_g->n_l; i++){
		if(i%10000 == 0) cout << "2... " << i << "/" << p_g->n_l << endl;
		Line* p_l_i = p_g->getLine(i);

		#ifdef OPENMP_MAIN
		#pragma omp parallel for num_threads(num_thread)
		#endif
		for(int j = 0; j < i; j++){
			Line* p_l_j = p_g->getLine(j);
			// 異なる辺を参照したはずなのに実体が同じならそれはおかしい。
			if(p_l_i == p_l_j){
				cout << "Warning: Line Reconstruction -1" << endl;
			}
			if(p_l_i->getVertex(0) == p_l_j->getVertex(0)
			&& p_l_i->getVertex(1) == p_l_j->getVertex(1)){
				cout << "Warning: Line Reconstruction -2" << endl;
			} else if(p_l_i->getVertex(0) == p_l_j->getVertex(1)
				&& p_l_i->getVertex(1) == p_l_j->getVertex(0)){
				cout << "Warning: Line Reconstruction -2" << endl;
			}
		}
	}
	
	#ifdef OPENMP_MAIN
	#pragma omp parallel for num_threads(num_thread)
	#endif
	// 細胞と辺を過不足なく対応づける
	for(int ci = 0; ci < p_g->n_c; ci++){
		vector<Line*> line_list;
		Cell* cp = p_g->getCell(ci);
		for(int si = 0; si < cp->getSurfaceSize(); si++){
			Surface* sp = cp->getSurface(si);
			for(int li = 0; li < 3; li++){
				Line* lp = sp->getLine(li);
				line_list.push_back(lp);
			}
		}
		// ここまでで，ci番目の細胞に属するすべての面について辺を数え上げ，そのリストがline_listに入っている
		// 次の2行はline_listの中の辺の重複をなくすアルゴリズム
		sort(begin(line_list), end(line_list));
		line_list.erase(unique(begin(line_list), end(line_list)), end(line_list));
		// c++11では，～.begin()よりもbegin(～)が推奨されている
		// 以上でci番目の細胞に属する辺の重複なしリストline_listが完成
		cp->setLine(line_list);	// このcpがもつ情報：細胞ciにはline_list内にある辺が属している
	}
	cout << "done." << endl;
}

void readLineDAT(){
	cout << "Reading Line File...";

	ifstream fin("data/_line.dat");
	if(!fin){
		cout << "Error: cannot open line file." << endl;
		exit(1);
	}

	int lsize;
	fin >> lsize;
	p_g->n_l = lsize;	

	// 辺(生成) <-- 2頂点
	for(int i = 0; i < p_g->n_l; i++){
		int l_vsize;
		fin >> l_vsize;	// 最初の数字は辺の構成点の数(=2)。

		Line* lp_new = new Line();
		int vid[2];
		fin >> vid[0] >> vid[1];
		Vertex* vp_0 = p_g->getVertex(vid[0]);
		Vertex* vp_1 = p_g->getVertex(vid[1]);
		lp_new->setVertex(vp_0, vp_1);
		lp_new->setId(i);
		p_g->pushLine(lp_new);

		// 頂点 <-- n辺
		vp_0->pushLine(lp_new);
		vp_1->pushLine(lp_new);
	}
	// 辺 <-- 2面
	for(int i = 0; i < p_g->n_l; i++){
		int l_ssize;
		fin >> l_ssize;	// 最初の数字は辺の構成面の数(=2)。
		Line* lp = p_g->getLine(i);
		int sid[2];
		fin >> sid[0] >> sid[1];
		Surface* sp_0 = p_g->getSurface(sid[0]);
		Surface* sp_1 = p_g->getSurface(sid[1]);
		lp->setSurface(sp_0, sp_1);
	}
	// 面 <-- 3辺
	for(int i = 0; i < p_g->n_s; i++){
		int s_lsize;
		fin >> s_lsize;	// 最初の数字は面の構成辺の数(=3)。
		Surface* sp = p_g->getSurface(i);
		int lid[3];
		fin >> lid[0] >> lid[1] >> lid[2];
		Line* lp_0 = p_g->getLine(lid[0]);
		Line* lp_1 = p_g->getLine(lid[1]);
		Line* lp_2 = p_g->getLine(lid[2]);
		sp->setLine(lp_0, lp_1, lp_2);
	}
	// 細胞 <-- n辺
	for(int i = 0; i < p_g->n_c; i++){
		int c_lsize;
		fin >> c_lsize;
		Cell* cp = p_g->getCell(i);

		vector<Line*> lp_list;
		for(int j = 0; j < c_lsize; j++){
			int lid;
			fin >> lid;
			Line* lp = p_g->getLine(lid);
			lp_list.push_back(lp);
		}
		cp->setLine(lp_list);
	}

	cout << " done." << endl;
}

void outputBoundVTK(){
// ***outer bound***
{
	char fname0[100];
	sprintf(fname0, "data/_outer_bound.vtk");
	ofstream fout0(fname0);
	if(!fout0){
		cout << "Error: cannot open outer bound file." << endl;
		exit(1);
	}

	fout0 << "# vtk DataFile Version 2.0" << endl;
	fout0 << "Outer_Bound" << endl;
	fout0 << "ASCII" << endl;
	fout0 << "DATASET UNSTRUCTURED_GRID" << endl;
	// 点の位置を書き込む
	//fout0 << "POINTS" << " " << (72+72) << " float" << endl;
	fout0 << "POINTS" << " " << 144 << " float" << endl;
	for(int i = 0; i < 72; i++){
		fout0 << owall_a*cos(5*i*M_PI/180) << " " << owall_b*sin(5*i*M_PI/180) << " " << -2.0 << endl;
	}
	for(int i = 0; i < 72; i++){
		fout0 << owall_a*cos(5*i*M_PI/180) << " " << owall_b*sin(5*i*M_PI/180) << " " << 6.0 << endl;
	}
	// 面の要素を書き込む
	//fout0 << "CELLS" << " " << 72 << " " << (72*5) << endl;
	fout0 << "CELLS" << " " << 72 << " " << 360 << endl;
	for(int i = 0; i < 72; i++){
		fout0 << "4" << " ";
		fout0 << i << " " << (i+1)%72 << " " << (72+(i+73)%72) << " " << (i+72) << endl;
	}
	//CELL_TYPES 図形要素 	// 5は三角形、3は線分
	fout0 << "CELL_TYPES" << " " << 72 << endl;
	for(int i = 0; i < 72; i++){
		fout0 << "9" << endl;
	}
}
// ***inner bound***
{
	char fname1[100];
	sprintf(fname1, "data/_inner_bound.vtk");
	ofstream fout1(fname1);
	if(!fout1){
		cout << "Error: cannot open inner bound file." << endl;
		exit(1);
	}

	fout1 << "# vtk DataFile Version 2.0" << endl;
	fout1 << "Inner_Bound" << endl;
	fout1 << "ASCII" << endl;
	fout1 << "DATASET UNSTRUCTURED_GRID" << endl;
	// 点の位置を書き込む
	//fout1 << "POINTS" << " " << (72+72) << " float" << endl;
	fout1 << "POINTS" << " " << 144 << " float" << endl;
	for(int i = 0; i < 72; i++){
		fout1 << iwall_a*cos(5*i*M_PI/180) << " " << iwall_b*sin(5*i*M_PI/180) << " " << -2.0 << endl;
	}
	for(int i = 0; i < 72; i++){
		fout1 << iwall_a*cos(5*i*M_PI/180) << " " << iwall_b*sin(5*i*M_PI/180) << " " << 6.0 << endl;
	}
	// 面の要素を書き込む
	//fout1 << "CELLS" << " " << (72) << " " << (72*5) << endl;
	fout1 << "CELLS" << " " << 72 << " " << 360 << endl;
	for(int i = 0; i < 72; i++){
		fout1 << "4" << " ";
		fout1 << i << " " << (i+1)%72 << " " << (72+(i+73)%72) << " " << (i+72) << endl;
	}
	//CELL_TYPES 図形要素 	// 5は三角形、3は線分
	fout1 << "CELL_TYPES" << " " << 72 << endl;
	for(int i = 0; i < 72; i++){
		fout1 << "9" << endl;
	}
}
}
void outputVTK(unsigned int step){
	char fname[100];
	sprintf(fname, "output/cm%010u.vtk", step);	// cell migration
	ofstream fout(fname);
	if(!fout){
		cout << "Error: cannot open output vtk file." << endl;
		exit(1);
	}

	// vtkファイルのヘッダ
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "Cell Migration" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET UNSTRUCTURED_GRID" << endl;

	// 点の位置を書き込む
	fout << "POINTS" << " " << p_g->n_v << " float" << endl;
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* vp = p_g->getVertex(i);
		double xx, yy, zz;
		xx = (vp->getLocation()).x;
		yy = (vp->getLocation()).y;
		zz = (vp->getLocation()).z;
		// Paraviewでfloat精度以上の値を扱おうとすると落ちる(double->float)
		if(fabs(xx) < 1e-6) xx = 0.0;
		if(fabs(yy) < 1e-6) yy = 0.0;
		if(fabs(zz) < 1e-6) zz = 0.0;
		fout << xx << " ";
		fout << yy << " ";
		fout << zz << endl;
	}
	// 面の要素を書き込む
	fout << "CELLS" << " " << (p_g->n_s) << " " << (p_g->n_s*4) << endl;
	for(int i = 0; i < p_g->n_s; i++){
		Surface* sp = p_g->getSurface(i);
		fout << "3" << " ";
		fout << sp->getVertex(0)->getId() << " ";
		fout << sp->getVertex(1)->getId() << " ";
		fout << sp->getVertex(2)->getId() << endl;
	}
	//CELL_TYPES 図形要素 	// 5は三角形、3は線分
	fout << "CELL_TYPES" << " " << (p_g->n_s) << endl;
	for(int i = 0; i < p_g->n_s; i++){
		fout << "5" << endl;
	}
	//CELL_DATA
	fout << "CELL_DATA" << " " << (p_g->n_s) << endl;
// ここの順番を入れ替えてしまうと，readVTK()にも影響する
// ***移動力を与えられている細胞かどうかを細胞表面の色で表現***
{
	fout << "SCALARS" << " " << "Rev_Cell" << " " << "int" << endl;
	fout << "LOOKUP_TABLE" << " " << "default" << endl;
	for(int i = 0; i < p_g->n_s; i++){
		Cell* cp = p_g->getSurface(i)->getVertex(0)->getCell();
		int tmp = 0;	// 青色(move)
		if(!(cp->getBoolRevcell())){
			tmp = 1;	// 赤色(stay)
		}
		fout << tmp << endl;
	}
}
// ***細胞のIDを細胞表面の色で表現***
{
	fout << "SCALARS" << " " << "Cell_ID" << " " << "int" << endl;
	fout << "LOOKUP_TABLE" << " " << "default" << endl;
	for(int i = 0; i < p_g->n_s; i++){
		Surface* sp = p_g->getSurface(i);
		fout << sp->getVertex(0)->getCell()->getId() << endl;
	}
}
// ***面のもつエネルギーの大きさをそれぞれの面の色で表現***
{
	fout << "SCALARS" << " " << "Surf_Energy" << " " << "float" << endl;
	fout << "LOOKUP_TABLE" << " " << "default" << endl;
	for(int i = 0; i < p_g->n_s; i++){
		Surface* sp = p_g->getSurface(i);
		// Paraviewでfloat精度以上の値を扱おうとすると落ちる(double->float)
		double tmp = sp->getSurfaceEnergy();
		if(fabs(tmp) < 1e-6) tmp = 0.0;
		fout << tmp << endl;
	}
}
// ***回転領域にある細胞かどうかを細胞表面の色で表現***
{
	fout << "SCALARS" << " " << "Central_Cell" << " " << "int" << endl;
	fout << "LOOKUP_TABLE" << " " << "default" << endl;
	for(int i = 0; i < p_g->n_s; i++){
		Cell* cp = p_g->getSurface(i)->getVertex(0)->getCell();
		int tmp = 1;	// 赤色(central)
		if(!(cp->getFlagCentral())){
			tmp = 0;	// 青色(not_central)
		}
		fout << tmp << endl;
	}
}

}

void outputLineDAT(){
	char fname[100];
	sprintf(fname, "data/_line.dat");
	ofstream fout(fname);
	if(!fout){
		cout << "Error: cannot open output line file." << endl;
		exit(1);
	}

	fout << p_g->n_l << endl;
	// 辺 <-- 2頂点
	for(int i = 0; i < p_g->n_l; i++){
		Line* lp = p_g->getLine(i);
		fout << "2" << " ";
		fout << lp->getVertex(0)->getId() << " ";
		fout << lp->getVertex(1)->getId() << endl;
	}
	// 辺 <-- 2面
	for(int i = 0; i < p_g->n_l; i++){
		Line* lp = p_g->getLine(i);
		fout << "2" << " ";
		fout << lp->getSurface(0)->getId() << " ";
		fout << lp->getSurface(1)->getId() << endl;
	}
	// 面 <-- 3辺
	for(int i = 0; i < p_g->n_s; i++){
		Surface* sp = p_g->getSurface(i);
		fout << "3" << " ";
		fout << sp->getLine(0)->getId() << " ";
		fout << sp->getLine(1)->getId() << " ";
		fout << sp->getLine(2)->getId() << endl;
	}
	// 細胞 <-- 辺
	for(int i = 0; i < p_g->n_c; i++){
		Cell* cp = p_g->getCell(i);
		fout << cp->getLineSize();
		for(int j = 0; j < cp->getLineSize(); j++){
			fout << " " << cp->getLine(j)->getId();
		}
		fout << endl;
	}
}

void outputInitialCellCG(){
	char fname[100];
	sprintf(fname, "data/_initial_cellcg.dat");
	ofstream fout(fname);
	if(!fout){
		cout << "Error: cannot open output initial cg file." << endl;
		exit(1);
	}

	for(int i = 0; i < p_g->n_c; i++){
		_vec<double> cell_cg = p_g->getCell(i)->getCentroidCell();
		fout << cell_cg.x << " ";
		fout << cell_cg.y << " ";
		fout << cell_cg.z << endl;
	}
}

void outputCellCG(unsigned int step){
	char fname[100];
	sprintf(fname, "data/cellcg%010u.cg", step);
	ofstream fout(fname);
	if(!fout){
		cout << "Error: cannot open output cg file." << endl;
		exit(1);
	}

	fout << "dt_step" << "\t" << dt_step << endl;
	fout << "dt_new" << "\t" << dt << endl;
	fout << "count_physical_time" << "\t" << count_physical_time << endl;
	fout << "physical_time" << "\t" << physical_time << endl;
	fout << "work_extf" << "\t" << p_g->work_extf << endl;
	fout << "total_work_extf" << "\t" << total_work_extf << endl;
	fout << "n_revcell" << "\t" << p_g->n_revcell << endl;
	fout << "total_n_revcell" << "\t" << total_n_revcell << endl;
	
	for(int i = 0; i < p_g->n_c; i++){
		_vec<double> cell_cg = p_g->getCell(i)->getCentroidCell();
		fout << cell_cg.x << " ";
		fout << cell_cg.y << " ";
		fout << cell_cg.z << endl;
	}

	for(int i = 0; i < (p_g->n_c)*(p_g->n_c-1)/2; i++){
		fout << p_g->getCellCellRepAdh(i) << endl;
	}

}

void readCellCG(unsigned int step){
	cout << "Reading Cell CG File...";
	char filename[50];
	sprintf(filename, "data/cellcg%010u.cg", step);
	ifstream fin(filename);
	if(!fin){
		cout << "Error: cannot open cell cg file." << endl;
		exit(1);
	}

	char dummy[50];
	double dummy_double;
	int dummy_int;

	fin >> dummy >> dt_step;
	fin >> dummy >> dt;
	fin >> dummy >> count_physical_time;
	fin >> dummy >> physical_time;
	fin >> dummy >> dummy_double;
	fin >> dummy >> total_work_extf;
	fin >> dummy >> dummy_int;
	fin >> dummy >> total_n_revcell;

	cout << " done." << endl;
}

void readVTK(unsigned int step_no){
	cout << "Reading VTK File...";

	char filename[50];
	sprintf(filename, "output/cm%010u.vtk", step_no);
	ifstream fin(filename);
	if(!fin){
		cout << "Error: cannot open input vtk file." << endl;
		exit(1);
	}
	// vtkの通りに座標をset
	char dummy[50];
	int dummy_int;
	for(int i = 0; i < 13; i++){
		fin >> dummy;	// Header部分
	}
	for(int i = 0; i < p_g->n_v; i++){
		Vertex* vp = p_g->getVertex(i);
		_vec<double> tmp_loc;
		fin >> tmp_loc.x >> tmp_loc.y >> tmp_loc.z;	// _ _ _
		vp->setLocation(tmp_loc);
	}
	fin >> dummy >> dummy_int >> dummy_int;	// CELLS _ _
	for(int i = 0; i < p_g->n_s; i++){
		fin >> dummy_int >> dummy_int >> dummy_int >> dummy_int;	// 3 _ _ _
	}
	fin >> dummy >> dummy_int;	// CELL_TYPES _
	for(int i = 0; i < p_g->n_s; i++){
		fin >> dummy_int;	// 5
	}
	fin >> dummy >> dummy_int;	// CELL_DATA _
	fin >> dummy >> dummy >> dummy;	// SCALARS Rev_Cell int
	fin >> dummy >> dummy;	// LOOKUP_TABLE default
	for(int i = 0; i < p_g->n_s; i++){
		Cell* cp = p_g->getSurface(i)->getVertex(0)->getCell();
		int tmp;
		fin >> tmp;	// 0 or 1
		if(tmp == 0){
			cp->setBoolRevcell(true);	// 青色(move)
		}else if(tmp == 1){
			cp->setBoolRevcell(false);	// 赤色(stay)
		}else{
			cout << "tmp = " << tmp << endl;
			cout << "Error: readVTK()--inadmissible value" << endl;
		}
	}

	cout << " done." << endl;
}

// 自然体積、自然面積、自然二面角、自然長の設定
void setEquilibriumState(){
	#ifdef OPENMP_MAIN
	#pragma omp parallel num_threads(num_thread)
	#endif
	{
		#ifdef OPENMP_MAIN
		#pragma omp for
		#endif
		for(int i = 0; i < p_g->n_c; i++){
			p_g->getCell(i)->setVolume0(0.520895*swelling_ratio);
		}

		#ifdef OPENMP_MAIN
		#pragma omp for
		#endif
		for(int i = 0; i < p_g->n_s; i++){
			p_g->getSurface(i)->setArea0(0.00163216);
		}
		
		#ifdef OPENMP_MAIN
		#pragma omp for
		#endif
		for(int i = 0; i < p_g->n_l; i++){
			p_g->getLine(i)->correctKAngle();
		}

		#ifdef OPENMP_MAIN
		#pragma omp for
		#endif
		for(int i = 0; i < p_g->n_l; i++){
			p_g->getLine(i)->setLen0(0.0638452);
		}
	}
}

int main(int argc, char* argv[]){
	if(argc != 2){
		cout << "Input a command-line argument." << endl;	exit(1);	// コマンドラインに引数を入力してください
	}

	p_g = new Global();
	unsigned int step = atoi(argv[1]);	// 文字列ポインタ型をint型に変換
	//unsigned int step_vtk = round(dt_vtk/dt);
	unsigned int init_step = step;	// ./a.out 0←この数字

	char fname0[100];
	sprintf(fname0, "data/_energy%010u.csv", init_step);
	ofstream fout0(fname0);

	if(!fout0){
		cout << "Error: cannot open energy file." << endl;	exit(1);
	}

	readSIG();			// sigファイルからの読み取り
	
	for(int i = 0; i < p_g->n_c; i++){
		p_g->getCell(i)->setBoolRevcell(false);
	}
	if(step != 0){
		readVTK(step);	// vtkファイルから読み込んでloc情報とrev_cell情報を更新
	}

	readParametersFile();	// datファイルからの読み取り

	if(step == 0){
		reconstLines();		// 辺の構成
		outputLineDAT();	// 辺の構成情報の出力
	}else{
		readLineDAT();		// 辺の構成情報の読み込み
		readCellCG(step);	// 中断データの読み込み
	}

	(void)Init_Random(seed);	// 必ずseed値を読み込んでから
	//outputBoundVTK();

	// 幾何学情報の更新(loc情報をもとにarea等を計算)
	// 各Cellが回転領域にあるかどうかの判定もここで行われる
	p_g->updateGeometry();

	// len0, area0など平衡状態の設定(不要loc,area情報)
	setEquilibriumState();

	// gridの定義(この関数の呼び出しは最初1回のみ)
	p_g->defineGrid();
	// 系の重心計算・系の重心位置を参照した仮移動量計算・gridとvertexの対応づけ
	// 斥力エネルギー・引力エネルギーを計算しているcalcEnergy()よりも前に記述(おそらく)
	p_g->setupGrid();

	// 初期状態における系のエネルギーの計算
//	p_g->calcEnergy();

	cout << "Step:" << " " << step << "\t" << "Energy:" << " " << p_g->getEnergy() << endl;

	fout0 << "Step" << "," << "Time" << "," << "TOTAL" << ",";
	fout0 << "Volume" << "," << "Area" << "," << "Angle" << "," << "Length" << ",";
	fout0 << "Repulsion" << "," << "Adhesion" << "," << "Boundary" << endl;
	if(step == 0){
		fout0 << step << "," << physical_time << "," << p_g->getEnergy() << ",";
		fout0 << p_g->u_vol << "," << p_g->u_are << "," << p_g->u_ang << "," << p_g->u_len << ",";
		fout0 << p_g->u_rep << "," << p_g->u_adh << "," << p_g->u_bnd << endl;
	}

	p_g->Initialize();

	//p_g->makeRandomChoice();	// 必ずupdateGeometryの後	// 各Cell*->rev_cellが確率的に変更される

	if(step == 0){
		outputVTK(step);
		outputInitialCellCG();
	}

	//p_g->countRevCell();

	step++;

	cout << "--------------------------------------------------------" << endl;

	while(step <= end_step){
//	while(physical_time < 4.78){
//	while(physical_time < 9.56){
//	while(total_n_revcell < 26882099){

		// stepの経過と移動細胞数を一定間隔おきに表示
		//if(step%100 == 0) cout << "Step: " << step << "\t" << "Migrating: " << p_g->n_revcell << " cells" << endl;

		// 頂点移動前の各cellの重心座標を保存
		/*
		#ifdef OPENMP_MAIN
		#pragma omp parallel for num_threads(num_thread)
		#endif
		for(int i = 0; i < p_g->n_c; i++){
			p_g->getCell(i)->storeCentroidCell();
		}*/

		// 運動方程式を計算し頂点を移動
		ode::solveMotionVertices();

		// 幾何学情報の更新(loc情報をもとにarea等を計算)
		// 各Cellが回転領域にあるかどうかの判定もここで行われる
		p_g->updateGeometry();

		// 各cellの重心移動量を求め、外力が系にした仕事を計算する		// v17
		//p_g->calcExtfWork();

		//if(step%500 == 0) cout << "\t" << "Work_extf: " << p_g->work_extf << "\t" << "Total_work_extf: " << total_work_extf << endl;

		// 系の重心計算・系の重心位置を参照した仮移動量計算・gridとvertexの対応づけ
		// 斥力エネルギー・引力エネルギーを計算しているcalcEnergy()よりも前に記述(おそらく)
		p_g->setupGrid();

		// 頂点移動後の系のエネルギーを計算
		p_g->calcEnergy();

		//p_g->makeRandomChoice();	// 必ずupdateGeometryの後
		
		// 条件設定
		double eps_time = dt*1e-2;	// 出力目安の値から少し小さいくらいでも許容とする。あと、刻み幅一定としたときに、たまに切り捨て誤差の累積で小さく評価される可能性があるため。

		// 条件を満たすstepで，VTK・SAVの出力と，頂点移動後の系のエネルギー・体節長さのファイル出力
		if(count_physical_time - dt_vtk + eps_time >= 0.0){
			outputVTK(step);

			cout << "Step: " << step << "\t" << "Time: " << physical_time << "\t" << "Energy: " << p_g->getEnergy() << endl;

			fout0 << step << "," << physical_time << "," << p_g->getEnergy() << ",";
			fout0 << p_g->u_vol << "," << p_g->u_are << "," << p_g->u_ang << "," << p_g->u_len << ",";
			fout0 << p_g->u_rep << "," << p_g->u_adh << "," << p_g->u_bnd << endl;
		
			count_physical_time -= dt_vtk;	// 修正済み

			outputCellCG(step);			
		}

		//p_g->countRevCell();

		p_g->Initialize();

		// 回転移動力を与えるCellのIDをrandomに再抽出後，配列chosen_cellsに再格納
		/*
		if((step-init_step)%rechoice_step == 0){
			p_g->makeRandomChoice();	// 必ずupdateGeometryの後
		}
		*/


		// 強制終了条件
		Cell* cp0 = p_g->getCell(0);
		Cell* cp1 = p_g->getCell(1);
		double dist = (cp1->getCentroidCell()).x - (cp0->getCentroidCell()).x;
		if(dist < 0){
			exit(0);
		}

		step++;
	}

	cout << "Finished." << endl;

	return 0;
}