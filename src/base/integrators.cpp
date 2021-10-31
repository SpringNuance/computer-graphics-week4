
#include "utility.hpp"
#include "particle_systems.hpp"
#include "integrators.hpp"

void eulerStep(ParticleSystem& ps, float step) {

	// YOUR CODE HERE (R1)
	// Implement an Euler integrator.
	const auto& x0 = ps.state(); // current state of PS aka x0 = X
	auto n = x0.size(); // size of current state of PS
	auto f0 = ps.evalF(x0); // current force calculated from PS's evalF aka f(X, t)
	auto x1 = State(n); // empty state of length n => both x0 and x1 has same length

	for (auto i = 0u; i < n; ++i) {
		x1[i] = x0[i] + step * f0[i]; // X(t + h) = X + h * f(X, t)
	}

	ps.set_state(x1); // Set new calculated state for PS
};

void trapezoidStep(ParticleSystem& ps, float step) {
	// YOUR CODE HERE (R3)
// Implement a trapezoid integrator.
	const auto& x0 = ps.state(); // current state of PS aka x0 = X
	auto n = x0.size(); // size of current state of PS
	auto f0 = ps.evalF(x0); // current force calculated from PS's evalF aka f0 = f(X, t)
	auto x1 = State(n); // empty state of length n => both x0 and x1 has same length
	auto xtrapez = State(n); // trapezoidal state which will be calculated from f0 and f1

	for (auto i = 0u; i < n; ++i) {
		x1[i] = x0[i] + step * f0[i]; // x1 = X + h*f0
	}

	auto f1 = ps.evalF(x1); // next force calculated from PS's evalF aka f1 = f(X + h * f0, t + h)
	for (auto i = 0u; i < n; ++i) {
		xtrapez[i] = x0[i] + (0.5f * step) * (f0[i] + f1[i]); // X(t + h) = X + 0.5 * step * (f0 + f1)
	}

	ps.set_state(xtrapez);

}

void midpointStep(ParticleSystem& ps, float step) {
	const auto& x0 = ps.state();
	auto n = x0.size();
	auto f0 = ps.evalF(x0);
	auto xm = State(n), x1 = State(n);
	for (auto i = 0u; i < n; ++i) {
		xm[i] = x0[i] + (0.5f * step) * f0[i];
	}
	auto fm = ps.evalF(xm);
	for (auto i = 0u; i < n; ++i) {
		x1[i] = x0[i] + step * fm[i];
	}
	ps.set_state(x1);
}

void rk4Step(ParticleSystem& ps, float step) {
	// EXTRA: Implement the RK4 Runge-Kutta integrator.
	
	// Calculating k1
	const auto& x1 = ps.state(); // current state of PS or x1 = X 
	auto n = x1.size(); // size of current state of PS
	auto k1 = ps.evalF(x1); // k1 = f(X, t)

	// Calculating k2
	auto x2 = State(n);
	for (auto i = 0u; i < n; ++i) {
		x2[i] = x1[i] + 0.5f * step * k1[i]; // x2 = X + h * k1/2 = x1 + h * k1/2
	}

	auto k2 = ps.evalF(x2); // k2 = f(X + h * k1/2, t + h/2)

	// Calculating k3
	auto x3 = State(n);
	for (auto i = 0u; i < n; ++i) {
		x3[i] = x1[i] + 0.5f * step * k2[i]; // x3 = X + h * k2/2 = x1 + h * k2/2
	}

	auto k3 = ps.evalF(x3); // k3 = f(X + h * k2/2, t + h/2)
	
	// Calculating k4
	auto x4 = State(n); 
	for (auto i = 0u; i < n; ++i) {
		x4[i] = x1[i] + step * k3[i]; // x4 = X + h * k3 = x1 + h * k3
	}

	auto k4 = ps.evalF(x4); // k4 = f(X + h * k3, t + h)

	// Calculating RK4 Runge-kutta
	auto xRK = State(n); 
	for (auto i = 0u; i < n; ++i) {
		xRK[i] = x1[i] + 1.f/6 * step * (k1[i] + 2.0f * k2[i] + 2.0f * k3[i] + k4[i]); 
		// X(t + h) = X  + 1/6 * h * (k1 + 2 * k2 + 2 * k3 + k4)
	    // or xRK   = x1 + 1/6 * h * (k1 + 2 * k2 + 2 * k3 + k4)
	}

	ps.set_state(xRK);
	
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void implicit_euler_step(ParticleSystem& ps, float step, SparseMatrix& J, SparseLU& solver, bool initial) {
	// EXTRA: Implement the implicit Euler integrator. (Note that the related formula on page 134 on the lecture slides is missing a 'h'; the formula should be (I-h*Jf(Yi))DY=-F(Yi))
}

void implicit_midpoint_step(ParticleSystem& ps, float step, SparseMatrix& J, SparseLU& solver, bool initial) {
	// EXTRA: Implement the implicit midpoint integrator.
}

void crank_nicolson_step(ParticleSystem & ps, float step, SparseMatrix & J, SparseLU & solver, bool initial) {
		// EXTRA: Implement the crank-nicolson integrator.
}
#endif
