/**
 * N-Body Gravitational Simulation
 *
 * Simulates gravitational interactions between celestial bodies using
 * the Velocity Verlet integration scheme. This algorithm is symplectic,
 * meaning it conserves energy over long timescales — crucial for orbital mechanics.
 *
 * Physics: Newton's Law of Universal Gravitation
 *   F = G * m1 * m2 / r²
 *
 * Author: Makar Ulesov
 * License: MIT
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>

// Physical constants (SI units)
constexpr double G = 6.67430e-11;  // Gravitational constant [m³/(kg·s²)]

/**
 * 3D Vector class for position, velocity, and acceleration calculations.
 * Operator overloading provides clean, physics-like syntax.
 */
struct Vec3 {
    double x{}, y{}, z{};

    Vec3() = default;
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }

    double magnitude() const { return std::sqrt(x*x + y*y + z*z); }
    double magnitude_squared() const { return x*x + y*y + z*z; }
};

/**
 * Represents a celestial body with mass, position, velocity, and acceleration.
 * The acceleration is stored to enable the Velocity Verlet algorithm.
 */
struct Body {
    std::string name;
    double mass;        // [kg]
    Vec3 position;      // [m]
    Vec3 velocity;      // [m/s]
    Vec3 acceleration;  // [m/s²] - computed from gravitational forces

    Body(std::string name, double mass, Vec3 pos, Vec3 vel)
        : name(std::move(name)), mass(mass), position(pos), velocity(vel) {}
};

/**
 * Computes gravitational acceleration on body i due to all other bodies.
 *
 * From Newton's second law: a = F/m = G * m_j / r²
 * Direction: unit vector pointing from body i toward body j
 *
 * A softening factor prevents numerical instabilities when bodies
 * approach very close distances (avoids division by near-zero).
 */
Vec3 compute_acceleration(const std::vector<Body>& bodies, size_t i) {
    constexpr double softening = 1e4;  // Softening parameter [m]
    Vec3 acc;

    for (size_t j = 0; j < bodies.size(); ++j) {
        if (i == j) continue;

        Vec3 r = bodies[j].position - bodies[i].position;  // Vector from i to j
        double dist_sq = r.magnitude_squared() + softening * softening;
        double dist = std::sqrt(dist_sq);

        // Gravitational acceleration: a = G * m_j * r_hat / r²
        // Combined: a = G * m_j * r / r³
        double factor = G * bodies[j].mass / (dist_sq * dist);
        acc += r * factor;
    }
    return acc;
}

/**
 * Velocity Verlet Integration Algorithm
 *
 * This is a second-order symplectic integrator, superior to Euler's method because:
 * 1. It conserves energy in Hamiltonian systems (crucial for orbits)
 * 2. It's time-reversible (a physical requirement)
 * 3. Second-order accuracy: error scales as O(dt²)
 *
 * Algorithm:
 *   1. x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt²
 *   2. Compute a(t+dt) from new positions
 *   3. v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
 */
void velocity_verlet_step(std::vector<Body>& bodies, double dt) {
    // Step 1: Update positions using current velocities and accelerations
    for (auto& body : bodies) {
        body.position += body.velocity * dt + body.acceleration * (0.5 * dt * dt);
    }

    // Step 2: Compute new accelerations at updated positions
    std::vector<Vec3> new_accelerations(bodies.size());
    for (size_t i = 0; i < bodies.size(); ++i) {
        new_accelerations[i] = compute_acceleration(bodies, i);
    }

    // Step 3: Update velocities using average of old and new accelerations
    for (size_t i = 0; i < bodies.size(); ++i) {
        bodies[i].velocity += (bodies[i].acceleration + new_accelerations[i]) * (0.5 * dt);
        bodies[i].acceleration = new_accelerations[i];
    }
}

/**
 * Initializes the Sun-Earth-Moon system with realistic parameters.
 * Data sourced from NASA JPL Horizons system.
 */
std::vector<Body> create_sun_earth_moon_system() {
    std::vector<Body> bodies;

    // Sun: At origin, stationary (barycenter approximation)
    bodies.emplace_back("Sun", 1.989e30, Vec3(0, 0, 0), Vec3(0, 0, 0));

    // Earth: 1 AU from Sun, circular orbit velocity
    // v = sqrt(G*M_sun/r) for circular orbit
    double earth_distance = 1.496e11;  // 1 AU in meters
    double earth_velocity = std::sqrt(G * 1.989e30 / earth_distance);
    bodies.emplace_back("Earth", 5.972e24,
                        Vec3(earth_distance, 0, 0),
                        Vec3(0, earth_velocity, 0));

    // Moon: 384,400 km from Earth, orbiting Earth
    double moon_distance = 3.844e8;
    double moon_velocity = std::sqrt(G * 5.972e24 / moon_distance);
    bodies.emplace_back("Moon", 7.342e22,
                        Vec3(earth_distance + moon_distance, 0, 0),
                        Vec3(0, earth_velocity + moon_velocity, 0));

    return bodies;
}

/**
 * Writes simulation state to CSV file.
 * Format: time, body1_x, body1_y, body1_z, body2_x, ...
 */
void write_csv_header(std::ofstream& file, const std::vector<Body>& bodies) {
    file << "time";
    for (const auto& body : bodies) {
        file << "," << body.name << "_x," << body.name << "_y," << body.name << "_z";
    }
    file << "\n";
}

void write_csv_row(std::ofstream& file, double time, const std::vector<Body>& bodies) {
    file << std::scientific << std::setprecision(6) << time;
    for (const auto& body : bodies) {
        file << "," << body.position.x << "," << body.position.y << "," << body.position.z;
    }
    file << "\n";
}

int main() {
    // Simulation parameters
    constexpr double dt = 3600.0;           // Time step: 1 hour [s]
    constexpr double total_time = 365.25 * 24 * 3600;  // 1 year [s]
    constexpr int output_interval = 24;     // Output every 24 hours

    std::cout << "N-Body Gravitational Simulation\n";
    std::cout << "================================\n";
    std::cout << "Integration: Velocity Verlet\n";
    std::cout << "Time step: " << dt << " s\n";
    std::cout << "Duration: 1 year\n\n";

    // Initialize system
    auto bodies = create_sun_earth_moon_system();

    // Initialize accelerations for Verlet algorithm
    for (size_t i = 0; i < bodies.size(); ++i) {
        bodies[i].acceleration = compute_acceleration(bodies, i);
    }

    // Open output file
    std::ofstream output("simulation_output.csv");
    if (!output.is_open()) {
        std::cerr << "Error: Could not open output file\n";
        return 1;
    }

    write_csv_header(output, bodies);

    // Main simulation loop
    int steps = static_cast<int>(total_time / dt);
    int output_count = 0;

    std::cout << "Running simulation...\n";

    for (int step = 0; step <= steps; ++step) {
        double current_time = step * dt;

        // Output at regular intervals
        if (step % output_interval == 0) {
            write_csv_row(output, current_time, bodies);
            ++output_count;

            // Progress indicator
            if (step % (steps / 10) == 0) {
                int progress = (step * 100) / steps;
                std::cout << "  Progress: " << progress << "%\n";
            }
        }

        // Advance simulation
        velocity_verlet_step(bodies, dt);
    }

    output.close();

    std::cout << "\nSimulation complete!\n";
    std::cout << "Output: " << output_count << " data points written to simulation_output.csv\n";
    std::cout << "\nRun 'python3 plot.py' to generate visualization.\n";

    return 0;
}
