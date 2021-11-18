/*
// For contact information
// Yasuhiro Inoue
// inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#ifndef ODE_SOLVER
#define ODE_SOLVER

namespace ode{
	void firstMotionSolver();
	void secondMotionSolver();
	void fourthMotionSolver();
	bool coreDormandPrinceSolver();
	void solveMotionVertices();
}
#endif