/*
// For contact information
// Yasuhiro Inoue
// inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
//----------------------------------------------------------------------------
//	_class.h
//----------------------------------------------------------------------------
//		push関数とset関数の名前の違いは、push_backしているか否か。情報取得はget関数
//		インスタンスは値コピー禁止。ポインタで渡す。
//----------------------------------------------------------------------------

#ifndef _CLASS_H
#define _CLASS_H

#include <vector>
#include <random>
#include <algorithm>
#include "vec.h"
#include "_parameters.h"

using namespace std;

class Global;
class Vertex;
class Line;
class Surface;
class Cell;

extern Global* p_g;

class Global{
	private:
		vector<Vertex*> p_v;
		vector<Line*> p_l;
		vector<Surface*> p_s;
		vector<Cell*> p_c;
		double u_total;
		

		// グリッドの個数
		int grid_size_x;
		int grid_size_y;
		int grid_size_z;
		_vec<double> sys_centroid;	// 系の重心(各細胞の重心の重心)
		_vec<double> shift;
		vector<int> nb_list;		// 隣接grid群のリスト(gridがないところには-1が入っている)
		vector<int> head_vid;		// 各グリッドの中で最も大きいvertexID
		vector<int> next_vid;		// 各vertexが所属するgridのIDリスト	// 同じgrid内のvertexをたどるための配列
		vector<int> cellcell_rep_adh;
		vector<int> box_revcell;

		//コピー禁止
		Global(const Global&);
		Global& operator= (const Global&);
	
	public:
		int n_v;
		int n_l;
		int n_s;
		int n_c;
		int n_revcell;	// v16
		double u_vol;
		double u_are;
		double u_ang;
		double u_len;
		double u_rep;
		double u_adh;
		double u_bnd;
		double work_extf;	// v17

		Global():n_v(0),n_l(0),n_s(0),n_c(0){};	//どうせ初期化しかしないので初期化リストを試してみる。
		~Global();	//基本的にnewで各インスタンスを生成しているのでデストラクタを一応用意する
		void pushVertex(Vertex*);
		Vertex* getVertex(int);
		void pushLine(Line*);
		Line* getLine(int);
		void pushSurface(Surface*);
		Surface* getSurface(int);
		void pushCell(Cell*);
		Cell* getCell(int);
		void updateGeometry();		//幾何学情報の更新
		void calcForces();			//力を計算する関数
		void setEnergy(double);
		void calcEnergy();
		double getEnergy();
		void calcExtfWork();		// 外力が系にした仕事の計算
		void calcRepulsiveForce_ij(Vertex*, Vertex*);	// 斥力計算
		double calcRepulsiveEnergy_ij(Vertex*, Vertex*);	// 斥力ポテンシャルエネルギー計算
		void calcAdhesiveForce_ij(Vertex*, Vertex*);	// 接着力計算
		double calcAdhesiveEnergy_ij(Vertex*, Vertex*);	// 接着力ポテンシャルエネルギー計算

		void defineGrid();
		int getGridId(int, int, int);
		void setupGrid();
		void setHeadId(int, int);
		int getHeadId(int);		// 指定した番号のgridの中の最大vertexIDを返す
		void makeRandomChoice();	// rev_zone=trueを満たすcellIDの中から，実際に回転力を与えるcellIDを，指定した割合に基づいて抽出し，リストchosen_cellsにそれらのIDを格納
		void countRevCell();
		int getCellCellRepAdh(int);

		void Initialize();	// 毎step初期化したいものがあれば
};

class Vertex{
	private:
		int id;
		_vec<double> loc;		//位置
		_vec<double> loc_t;		//時刻tにおける位置をここに格納しておく。精度を上げた場合の数値積分に用いる。
		_vec<double> loc_4th;
		_vec<double> frc;		//力
		//_vec<double> frc_rk[4];	//高次精度における力をこの配列に格納しておく。
		_vec<double> frc_rk[7];	//高次精度における力をこの配列に格納しておく。
		_vec<double> frc_omp[num_thread];	//OPENMP用

		vector<Line*> p_l;
		vector<Surface*> p_s;
		Cell* p_c;

		int grid_id;	// 所属するgridのID

		_vec<double> normal;  // 頂点の単位法線ベクトル(頂点が属する面の法線ベクトルから推定)
		double area;			// 頂点の占有面積

		//コピー禁止
		Vertex(const Vertex&);
		Vertex& operator= (const Vertex&);
	
	public:
		_vec<double> rot_t;	// 回転軌道の接線方向単位ベクトル
		_vec<double> rot_n;	// 回転軌道の法線方向単位ベクトル

		Vertex(int);
		void pushLine(Line*);
		Line* getLine(int);
		void pushSurface(Surface*);
		Surface* getSurface(int);
		void setCell(Cell*);
		Cell* getCell();
		void setLocation(_vec<double>);
		_vec<double> getLocation();
		void copyLocationToLocationT();
		_vec<double> getLocationT();
		int getId();
		void setForce(_vec<double>);
		_vec<double> getForce();
		void addForce(_vec<double>);
		void addForce(_vec<double>, int);
		void setForceZero();
		void setForceZero(int);
		void copyForceRk(int);
		_vec<double> getForceRk(int);
		void sumForceOmp();
		void setLocation4th(_vec<double>);
		_vec<double> getLocation4th();
		void calcNormalAndAreaVertex();	// Vertexの単位法線ベクトルと占有面積の計算
		_vec<double> getNormal();	// calcしてから使う
		double getArea();
		int getLineSize();
		int getSurfaceSize();
		void calcExternalForce();		// 外力計算
		void calcBoundaryForce();		// 境界から受ける力の計算
		double calcBoundaryEnergy();	// 境界から受ける力のポテンシャルエネルギー計算

		void setGridId(int);	// 所属するgridのIDを定める
		int getGridId();		// 所属するgridのIDを返す
};

class Line{
	private:
		int id;
		Surface* p_s[2];
		Vertex* p_v[2];
		double k_angle_pos;
		double len;
		double len_0;
	
		//コピー禁止
		Line(const Line&);
		Line& operator= (const Line&);

	public:
		Line():p_s{nullptr, nullptr},p_v{nullptr, nullptr}{};	//せっかくなのでnullptrで初期化したいのでこうしているが、もっとスマートな書き方がある気がする。
		void setId(int);
		int getId();
		void setVertex(Vertex*, Vertex*);
		Vertex* getVertex(int);
		void setSurface(int, Surface*);
		void setSurface(Surface*, Surface*);
		Surface* getSurface(int);
		void correctKAngle();
		void calcAngleForce();
		double calcAngleEnergy();
		void setLen0(double);
		double getLen0();
		void copyLenToLen0();
		void calcLength();
		double getLength();
		void calcLengthForce();
		double calcLengthEnergy();
};

class Surface{
	private:
		int id;
		_vec<double> normal;
		_vec<double> centroid;
		Line* p_l[3];
		Vertex* p_v[3];
		double area;
		double area_0;
		double energy;
		
		//コピー禁止
		Surface(const Surface&);
		Surface& operator= (const Surface&);

	public:
		Surface():p_l{nullptr, nullptr, nullptr},p_v{nullptr, nullptr, nullptr}{};
		void setId(int);
		int getId();
		void setVertex(Vertex*, Vertex*, Vertex*);
		Vertex* getVertex(int a);
		void setLine(Line*, Line*, Line*);
		Line* getLine(int a);
		void calcCentroidTriangle();
		_vec<double> getCentroid();
		void calcNormalVector();		// 面の単位法線ベクトル の計算
		_vec<double> getNormal();		// 面の単位法線ベクトル の取得	// calcしてから使う
		void setArea0(double);
		double getArea0();
		void setArea(double);
		void calcNormalAndAreaTriangle();
		double getArea();
		void copyAreaToArea0();
		void calcAreaForce();
		double calcAreaEnergy();
		void setSurfaceEnergy(double);
		void addSurfaceEnergy(double);
		double getSurfaceEnergy();
};

class Cell{
	private:
		int id;
		vector<Vertex*> p_v;
		vector<Line*> p_l;
		vector<Surface*> p_s;
		_vec<double> centroid;
		_vec<double> prev_centroid;
		double volume;
		double volume_0;
		bool flag_central;	// 中心部にあるならtrue，周辺部にあるならfalse
		bool rev_cell;		// 回転移動力が与えられているならtrue，与えられていないならfalse
		_vec<double> extf;	// 各細胞重心に加わる外力
		_vec<double> centroid_move;	// 各細胞重心の移動ベクトル

		//コピー禁止
		Cell(const Cell&);
		Cell& operator= (const Cell&);

	public:
		_vec<double> rot_t;	// 回転軌道の接線方向単位ベクトル
		_vec<double> rot_n;	// 回転軌道の法線方向単位ベクトル

		Cell(){};
		void setId(int);
		int getId();
		void setVertex(vector<Vertex*>);
		Vertex* getVertex(int a);
		void setLine(vector<Line*>);
		Line* getLine(int a);
		void setSurface(vector<Surface*>);
		Surface* getSurface(int a);
		void calcCentroidCell();
		_vec<double> getCentroidCell();
		void storeCentroidCell();	// v17
		_vec<double> getStoredCentroidCell();	// v17
		void calcMovementVector();	// v17
		_vec<double> getMovementVector();	// v17
		void judgeCentralOrNot();
		bool getFlagCentral();
		void setBoolRevcell(bool);
		bool getBoolRevcell();
		void setVolume0(double);
		double getVolume0();
		void setVolume(double);
		void calcVolumeCell();
		double getVolume();
		void copyVolumeToVolume0();
		void calcVolumeForce();
		double calcVolumeEnergy();
		void calcInternalForce();		// 内力計算
		void calcCellMigrationForce();
		void setExternalForceZero();
		_vec<double> getExternalForce();	// v17
		int getVertexSize();
		int getLineSize();
		int getSurfaceSize();
};

#endif