#include "particle_systems.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>

using namespace std;
using namespace FW;

namespace {

	inline Vec3f fGravity(float mass) {
		return Vec3f(0, -9.8f * mass, 0);
	}

	// force acting on particle at pos1 due to spring attached to pos2 at the other end
	inline Vec3f fSpring(const Vec3f& pos1, const Vec3f& pos2, float k, float rest_length) {
		// YOUR CODE HERE (R2)
		Vec3f dist = pos1 - pos2; // vec(PiPj)
		// F(Pi, Pj) = K(L0 - ||vec(PiPj)||) * vec(PiPj)/||vec(PiPj)||
		Vec3f fSpring = k * (rest_length - dist.length()) * (dist / dist.length());
		return fSpring;
	}

	inline Vec3f fDrag(const Vec3f& v, float k) {
		// YOUR CODE HERE (R2)
		Vec3f fDrag = -k * v;
		return fDrag;
	}
	
	inline int pos(int index) {
		return 2 * index;
	}

	inline int vel(int index) {
		return (2 * index) + 1;
	}

	inline unsigned idx(unsigned num, unsigned x, unsigned y) {
		return x + y * num;
	}

	

} // namespace

/*Simple system*/
void SimpleSystem::reset() {
	current_state_ = State(1, Vec3f(0, radius_, 0));
}

State SimpleSystem::evalF(const State& state) const {
	State f(1, Vec3f(-state[0].y, state[0].x, 0));
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H
// using the implicit Euler method, the simple system should converge towards origin -- as opposed to the explicit Euler, which diverges outwards from the origin.
void SimpleSystem::evalJ(const State&, SparseMatrix& result, bool initial) const {
	if (initial) {
		result.coeffRef(1, 0) = 1.0f;
		result.coeffRef(0, 1) = -1.0f;
	}
}
#endif

Points SimpleSystem::getPoints() {
	return Points(1, current_state_[0]);
}

Lines SimpleSystem::getLines() {
	static const auto n_lines = 50u;
	auto l = Lines(n_lines * 2);
	const auto angle_incr = 2 * FW_PI / n_lines;
	for (auto i = 0u; i < n_lines; ++i) {
		l[2 * i] = l[2 * i + 1] =
			Vec3f(radius_ * FW::sin(angle_incr * i), radius_ * FW::cos(angle_incr * i), 0);
	}
	rotate(l.begin(), l.begin() + 1, l.end());
	return l;
}



/*Spring system*/
void SpringSystem::reset() {
	const auto start_pos = Vec3f(0.1f, -0.5f, 0.0f);
	const auto spring_k = 30.0f;
	const auto rest_length = 0.5f;
	current_state_ = State(4);
	// YOUR CODE HERE (R2)
	// Set the initial state for a particle system with one particle fixed
	// at origin and another particle hanging off the first one with a spring.
	// Place the second particle initially at start_pos.

	Vec3f origin = Vec3f(0.0f, 0.0f, 0.0f);

	spring_.rlen = (start_pos - origin).length();
	spring_.k = spring_k;
	current_state_[0] = origin;
	current_state_[1] = origin; 
	current_state_[2] = start_pos;
	current_state_[3] = origin; 
}

State SpringSystem::evalF(const State& state) const {
	const auto drag_k = 0.5f;
	const auto mass = 1.0f;
	State f(4);
	// YOUR CODE HERE (R2)
	// Return a derivative for the system as if it was in state "state".
	// You can use the fGravity, fDrag and fSpring helper functions for the forces.
	Vec3f origin = Vec3f(0.0f, 0.0f, 0.0f);
	Vec3f combinedForce = fGravity(mass) + fDrag(state[3], drag_k) + fSpring(state[2], state[0], spring_.k, spring_.rlen);

	f[0] = origin;
	f[1] = origin;
	f[2] = state[3];
	f[3] = combinedForce / mass; // acceleration = F/m

	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

// This is a very useful read for the Jacobians of the spring forces. It deals with spring damping as well, we don't do that -- our drag is simply a linear damping of velocity (that results in some constants in the Jacobian).
// http://blog.mmacklin.com/2012/05/04/implicitsprings/

void SpringSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {
	const auto drag_k = 0.5f;
	const auto mass = 1.0f;
	// EXTRA: Evaluate the Jacobian into the 'result' matrix here. Only the free end of the spring should have any nonzero values related to it.
}
#endif

Points SpringSystem::getPoints() {
	auto p = Points(2);
	p[0] = current_state_[0]; p[1] = current_state_[2];
	return p;
}

Lines SpringSystem::getLines() {
	auto l = Lines(2);
	l[0] = current_state_[0]; l[1] = current_state_[2];
	return l;
}

/*Pendulum system*/
void PendulumSystem::reset() {
	const auto spring_k = 1000.0f;
	const auto start_point = Vec3f(0);
	const auto end_point = Vec3f(0.05, -1.5, 0);
	current_state_ = State(2 * n_);
	// YOUR CODE HERE (R4)
	// Set the initial state for a pendulum system with n_ particles
	// connected with springs into a chain from start_point to end_point with uniform intervals.
	// The rest length of each spring is its length in this initial configuration.
	springs_.clear();
	current_state_[pos(0)] = start_point;
	current_state_[vel(0)] = Vec3f(0.0f, 0.0f, 0.0f);
	

	Vec3f vecStep = (end_point - start_point) / (n_ - 1);

	for (unsigned int i = 1; i < n_; ++i) {
		Spring s = Spring(i, i - 1, spring_k, vecStep.length());
		springs_.push_back(s);
		current_state_[pos(i)] = current_state_[pos(i - 1)] + vecStep;
		current_state_[vel(i)] = Vec3f(0.0f, 0.0f, 0.0f);
	}

}

State PendulumSystem::evalF(const State& state) const {
	const auto drag_k = 0.5f;
	const auto mass = 0.5f;
	auto f = State(2 * n_); // f(velo0, acce0, velo1, acce1, ...)
	// YOUR CODE HERE (R4)
	// As in R2, return a derivative of the system state "state".
	f[pos(0)] = state[vel(0)]; // zero velocity for first particle
	f[vel(0)] = Vec3f(0.0f, 0.0f, 0.0f); // zero acceleration for first particle

    // vector of spring forces
	std::vector<Vec3f> springForce(n_, Vec3f(0.0f, 0.0f, 0.0f));

	// There are n partiles and n - 1 springs that connect these particles together
	for (unsigned int i = 0; i < n_ - 1; ++i) {
		auto belowPar = springs_[i].i1;
		auto abovePar = springs_[i].i2;
		auto stiff = springs_[i].k;
		auto rlen = springs_[i].rlen;
		springForce[abovePar] += fSpring(state[pos(abovePar)], state[pos(belowPar)], stiff, rlen);
		springForce[belowPar] += fSpring(state[pos(belowPar)], state[pos(abovePar)], stiff, rlen);
	}

	for (unsigned int i = 1; i < n_; ++i) {
		Vec3f combinedForce = springForce[i] + fDrag(state[vel(i)], drag_k) + fGravity(mass);
		f[pos(i)] = state[vel(i)]; // velocity
		f[vel(i)] = combinedForce / mass; //acceleration 
	}
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void PendulumSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {

	const auto drag_k = 0.5f;
	const auto mass = 0.5f;

	// EXTRA: Evaluate the Jacobian here. Each spring has an effect on four blocks of the matrix -- both of the positions of the endpoints will have an effect on both of the velocities of the endpoints.
}
#endif


Points PendulumSystem::getPoints() {
	auto p = Points(n_);
	for (auto i = 0u; i < n_; ++i) {
		p[i] = current_state_[i * 2];
	}
	return p;
}

Lines PendulumSystem::getLines() {
	auto l = Lines();
	for (const auto& s : springs_) {
		l.push_back(current_state_[2 * s.i1]);
		l.push_back(current_state_[2 * s.i2]);
	}
	return l;
}

/*CLoth system*/
void ClothSystem::reset() {
	const auto spring_k = 300.0f;
	const auto width = 1.5f, height = 1.5f; // width and height of the whole grid
	current_state_ = State(2 * x_*y_);
	// YOUR CODE HERE (R5)
	// Construct a particle system with a x_ * y_ grid of particles,
	// connected with a variety of springs as described in the handout:
	// structural springs, shear springs and flex springs.
	unsigned xNum = x_;
	unsigned yNum = y_;
	float xGridUnit = width / (xNum - 1);
	float yGridUnit = height / (yNum - 1);
	float diagonalUnit = FW::sqrt(xGridUnit * xGridUnit + yGridUnit * yGridUnit);
	springs_.clear();
	for (unsigned int x = 0; x < xNum; x++) {
		for (unsigned int y = 0; y < yNum; y++) {
			current_state_[pos(idx(xNum, x, y))] = Vec3f(width * (float)x / (xNum - 1) - 0.5f * width, 0.0f, -(height) * ((float)y / (yNum - 1)));
		
			// Structural springs of cloth grid
			if (x + 1 < xNum) {
				Spring spring(idx(xNum, x, y), idx(xNum, x + 1, y), spring_k, xGridUnit);
				springs_.push_back(spring);
			}

			if (y + 1 < yNum) {
				Spring spring(idx(xNum, x, y), idx(xNum, x, y + 1), spring_k, yGridUnit);
				springs_.push_back(spring);
			}

			// Shear springs of the cloth grid
			if (x + 1 < xNum && y > 0) {
				Spring spring(idx(xNum, x, y), idx(xNum, x + 1, y - 1), spring_k, diagonalUnit);
				springs_.push_back(spring);
			}

			if (x + 1 < xNum && y + 1 < yNum) {
				Spring spring (idx(xNum, x, y), idx(xNum, x + 1, y + 1), spring_k, diagonalUnit);
				springs_.push_back(spring);
			}

			// Flex springs of the cloth grid
			if (x + 2 < xNum) {
				Spring spring(idx(xNum, x, y), idx(xNum, x + 2, y), spring_k, xGridUnit * 2);
				springs_.push_back(spring);
			}

			if (y + 2 < yNum) {
				Spring spring(idx(xNum, x, y), idx(xNum, x, y + 2), spring_k, yGridUnit * 2);
				springs_.push_back(spring);
			}
		}
	}
}


State ClothSystem::evalF(const State& state) const {
	const auto drag_k = 0.08f;
	const auto n = x_ * y_;
	static const auto mass = 0.025f;
	auto f = State(2 * n);
	// YOUR CODE HERE (R5)
	// This will be much like in R2 and R4.
	// vector of spring forces
	std::vector<Vec3f> springForce(n, Vec3f(0.0f, 0.0f, 0.0f));
	unsigned xNum = x_;
	unsigned yNum = y_;

	for (int i = 0; i < springs_.size(); i++) {
		auto belowPar = springs_[i].i1;
		auto abovePar = springs_[i].i2;
		auto stiff = springs_[i].k;
		auto rlen = springs_[i].rlen;
		springForce[abovePar] += fSpring(state[pos(abovePar)], state[pos(belowPar)], stiff, rlen);
		springForce[belowPar] += fSpring(state[pos(belowPar)], state[pos(abovePar)], stiff, rlen);
	}

	for (unsigned int i = 0; i < n; i++) {
		Vec3f combinedForce = fGravity(mass) + fDrag(state[vel(i)], drag_k) + springForce[i];

		f[pos(i)] = state[vel(i)]; // velocity
		f[vel(i)] = combinedForce / mass; //acceleration
	}

	// fixing two upper corner particles
	f[vel(idx(xNum, 0, 0))] = Vec3f(0.0f, 0.0f, 0.0f); 
	f[vel(idx(xNum, xNum - 1, 0))] = Vec3f(0.0f, 0.0f, 0.0f);

	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void ClothSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {
	const auto drag_k = 0.08f;
	static const auto mass = 0.025f;

	// EXTRA: Evaluate the Jacobian here. The code is more or less the same as for the pendulum.
}

#endif

Points ClothSystem::getPoints() {
	auto n = x_ * y_;
	auto p = Points(n);
	for (auto i = 0u; i < n; ++i) {
		p[i] = current_state_[2 * i];
	}
	return p;
}

Lines ClothSystem::getLines() {
	auto l = Lines();
	for (const auto& s : springs_) {
		l.push_back(current_state_[2 * s.i1]);
		l.push_back(current_state_[2 * s.i2]);
	}
	return l;
}



/*Fluid system*/
State FluidSystem::evalF(const State&) const {
	return State();
}

