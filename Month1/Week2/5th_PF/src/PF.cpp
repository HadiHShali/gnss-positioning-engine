//==========================================================//
// Week 2, Day 5: PArticle Filters - GPS denied Navigation  //
//==========================================================//

// State : 2D position:  [x, y]
// 500 particles, Sequencial Importance Resampling (SIR) algorithm with resampling

// SECTION 0: Header
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace std;

// SECTION 1: Particle Structure
struct Particle
{
	double x, y; // 2D position
	double weight; // importance weight
};

// SECTION 2: Particle Class
class ParticleFilter
{
public:
	int N; //number of particles

	vector<Particle> particles; // It is the standard way to hold a list of "things" whose count changes during runtime.
	// vector: the C++ standard container — like a smart array
    // <Particle> : what type of items it holds ( Particle struct)
    // particles: the variable name we chose
	// particles:
	//	├──[0] →{ x = 2.3,  y = 5.7,  weight = 0.001 }
	//	├──[1] →{ x = 4.1,  y = 3.2,  weight = 0.003 }
	//   ...

	mt19937 rng;  // random number generator
	double process_noise; // sigma for predict step
	double meas_noise; // sigma for measurement


	// Constructor
	ParticleFilter(int n_particles, double x_min, double x_max, double y_min, double y_max,
		double p_noise, double m_noise)
		:N(n_particles), rng(random_device{}()), process_noise(p_noise), meas_noise(m_noise) 
		// random_device : creates a hardware-based random source
		// random_device{}() : calls it to get one truly random seed number
		// rng(...) : seeds the Mersenne Twister with that number
	{
		// Initialize particles uniformly in search area
		// uniform_real_distribution : random number generators that produce uniform values in a range
		uniform_real_distribution<double> dx(x_min, x_max);
		uniform_real_distribution<double> dy(y_min, y_max);

		// This pre-allocates memory for N particles in the vector. It's a performance optimization
		// — without it, the vector would have to keep growing its memory buffer every few push_back calls.
		// Reserving once is much faster.
		particles.reserve(N);

		// create the particles
		for (int i=0; i < N; i++)
		{
			particles.push_back({ dx(rng), dx(rng), 1.0 / N });
		}
	}
	// Prediction
	// Move each particle: xᵢ_new = f(xᵢ) + wᵢ
	void predict(double dx, double dy) 
		// dx, dy = how much the vehicle is expected to move East and North this epoch
		// These could come from a control input (steering, throttle) or a process model. In real code they often come
		// from IMU integration over the time step.
	{
		normal_distribution<double> noise(0.0, process_noise);
		// create a Gaussain Noise Distribution with Mean=0, and std=process_noise

		//for (auto& p : particles) // loops over every particle in the vector
		//						  // auto: Compiler figures out the type (it is Particle here)
		//{
		//	p.x += dx + noise(rng);
		//	p.y += dy + noise(rng);
		//}
		for (int i=0; i < particles.size(); i++) //loops over every particle in the vector
		{
			Particle& p = particles[i];  // & : Reference — modifies the original particle, not a copy
			p.x += dx + noise(rng); // add mostion plus noise
			p.y += dy + noise(rng);
		}
	}

	// Weight (Importance Sampling)
	// wᵢ = exp(-0.5 * (z - h(xᵢ))² / R)
	void update(double meas_x, double meas_y)
	{
		double sum_w = 0.0;
		for (auto& p : particles)
		{
			double err_x = meas_x - p.x;
			double err_y = meas_y - p.y;
			double d2 = err_x * err_x + err_y * err_y;
			p.weight = exp(-0.5 * d2 / (meas_noise * meas_noise));
			sum_w += p.weight;
		}
		// Normalize weight
		for (auto& p : particles) p.weight /= sum_w ;
	}


	// Resmaple (Synthetic Resampling)
	// What resampling does (the goal):
	//				You have N particles with weights — some big, some tiny.You want to create N new particles where :

	//              - Particles with high weight appear many times(duplicated)
	//              - Particles with low weight are deleted
	//              - All new particles have equal weight(1 / N)
	// 
	// Imagine laying all particles in a line, each occupying space proportional to its weight:
	//	weights:   0.05   0.30      0.10    0.45     0.10
	//              ┌────┬────────┬──────┬────────────┬──────┐
	//	particles : │ P0 │   P1   │  P2  │     P3     │  P4  │
	//	            └────┴────────┴──────┴────────────┴──────┘
	//	position : 0   0.05    0.35   0.45         0.90   1.00
	//	↑ total = 1

	// Then drop N equally-spaced points along this line. Wherever a point lands, that's the particle you pick:
	// N=5 marks at positions:  0.1   0.3   0.5   0.7   0.9
	//                           ↓     ↓     ↓     ↓     ↓
	//                    Pick : P1    P1    P3    P3    P4
	// P1 gets picked 2 times (big weight), P0 and P2 are dropped (small weights). This is exactly what we want.

	void resample()
	{
		// Create empty container and reserve a space in advance
		vector<Particle> new_particles;
		new_particles.reserve(N);

		// Random starting offset : for diversity
		// this draws a random r between 0 and 0.002 (= 1/500) 
		uniform_real_distribution<double> dist(0.0, 1.0 / N);
		double r = dist(rng);

		// Initialize Cumulative Distribution Function tracker
		// cumulative weight up through particle i (starts at particle 0's weight)
		double c = particles[0].weight;

		int i = 0;
		for (int m = 0; m < N; m++)
		{
			// (double)m : it converts m from an int to a double
			double u = r + (double)m / N;
			while (u > c && i < N - 1) //  Find which particle this u falls into
			{
				i++;
				c += particles[i].weight;
			}

			Particle pNew = particles[i]; // Copy particle i (its x, y values)
			pNew.weight = 1.0 / N;        // Reset its weight to 1/N (all new particles are equally weighted)
			new_particles.push_back(pNew);// Add it to the new list
		}
		particles = new_particles;  // Replace the old particles
	}

	// Estimate:

	// x̂ = mean of particles
	void estimate(double& x_est, double& y_est) // &x_est, &y_est — passed by reference, so the function fills them with the result
	{
		x_est = 0.0; y_est = 0.0;
		for (const auto& p : particles)  //const: promise NOT to modify p — read-only access
		{
			x_est += p.x * p.weight;
			y_est += p.y * p.weight;
		}
	}



};

// MAin Program

int main() {
	// Initialize particles in a 100x100 area around (0,0)
	ParticleFilter pf(500, -50.0, 50.0, -50.0, 50.0, 0.5, 2.0);

	// Simulated vehicle motion: moves 1m east, 0.5m north per epoch
	double true_x = 0.0, true_y = 0.0;
	double dx = 1.0, dy = 0.5;

	// Simulated noisy GPS measurements
	mt19937 rng(42);
	normal_distribution<double> gps_noise(0.0, 2.0);

	cout << "=== Particle Filter: 2D Navigation ===" << endl;
	cout << "Epoch  True_X  True_Y  GPS_X   GPS_Y   PF_X    PF_Y    Error" << endl;

	for (int k = 0; k < 20; k++) {
		// True vehicle moves
		true_x += dx; true_y += dy;

		// Step 1: PREDICT — move particles
		pf.predict(dx, dy);

		// Step 2: WEIGHT — incorporate noisy GPS
		double gps_x = true_x + gps_noise(rng);
		double gps_y = true_y + gps_noise(rng);
		pf.update(gps_x, gps_y);

		// Step 3: RESAMPLE
		pf.resample();

		// Step 4: ESTIMATE
		double est_x, est_y;
		pf.estimate(est_x, est_y);
		double err = sqrt(pow(est_x - true_x, 2) + pow(est_y - true_y, 2));

		cout << "  " << k + 1 << "     " << true_x << "     " << true_y
			<< "     " << gps_x << "   " << gps_y
			<< "   " << est_x << "   " << est_y
			<< "   " << err << endl;
	}
	return 0;
}
