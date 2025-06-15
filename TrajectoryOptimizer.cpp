#include "TrajectoryOptimizer.hpp"
#include <limits>
#include <cmath>
#include <chrono>
#include <fstream>

AirspaceConstraint AirspaceConstraint::createSphere(
    double center_x, double center_y, double center_z, double radius,
    double t_start, double t_end, const std::string& desc) {
    
    AirspaceConstraint constraint;
    constraint.constraint_type = AirspaceConstraintType::SPHERE;
    constraint.parameters = {center_x, center_y, center_z, radius};
    constraint.start_time = t_start;
    constraint.end_time = t_end;
    constraint.is_active = true;
    constraint.description = desc;
    return constraint;
}

AirspaceConstraint AirspaceConstraint::createCylinder(
    double center_x, double center_y, double min_z, double max_z, double radius,
    double t_start, double t_end, const std::string& desc) {
    
    AirspaceConstraint constraint;
    constraint.constraint_type = AirspaceConstraintType::CYLINDER;
    constraint.parameters = {center_x, center_y, min_z, max_z, radius};
    constraint.start_time = t_start;
    constraint.end_time = t_end;
    constraint.is_active = true;
    constraint.description = desc;
    return constraint;
}

AirspaceConstraint AirspaceConstraint::createCorridor(
    double start_x, double start_y, double start_z,
    double end_x, double end_y, double end_z, double width,
    double t_start, double t_end, const std::string& desc) {
    
    AirspaceConstraint constraint;
    constraint.constraint_type = AirspaceConstraintType::CORRIDOR;
    constraint.parameters = {start_x, start_y, start_z, end_x, end_y, end_z, width};
    constraint.start_time = t_start;
    constraint.end_time = t_end;
    constraint.is_active = true;
    constraint.description = desc;
    return constraint;
}

AirspaceConstraint AirspaceConstraint::createTimeWindow(
    double target_x, double target_y, double target_z,
    double earliest_arrival, double latest_arrival,
    const std::string& desc) {
    
    AirspaceConstraint constraint;
    constraint.constraint_type = AirspaceConstraintType::TIME_WINDOW;
    constraint.parameters = {target_x, target_y, target_z, earliest_arrival, latest_arrival};
    constraint.start_time = earliest_arrival;
    constraint.end_time = latest_arrival;
    constraint.is_active = true;
    constraint.description = desc;
    return constraint;
}

TrajectoryOptimizer::TrajectoryOptimizer()
    : system(), opt_params() {}

TrajectoryOptimizer::TrajectoryOptimizer(const SystemParams& sys, const OptimizationParams& opt)
    : system(sys), opt_params(opt) {}

TrajectoryOptimizer::TrajectoryOptimizer(double mass, double dt, int N) 
    : system(mass), opt_params(N, dt) {}

void TrajectoryOptimizer::setSystemParams(const SystemParams& params) {
    system = params;
}

void TrajectoryOptimizer::setOptimizationParams(const OptimizationParams& params) {
    opt_params = params;
}

Function TrajectoryOptimizer::createSolver4D(const MX& x, const MX& u, const MX& t, const std::string& objective) {
    int N = opt_params.N;
    double mass = system.mass;
    
    if (fabs(mass) < 1e-6) {
        mass = 80000.0;
    }
    
    std::vector<MX> equations;
    
    for (int k = 0; k < N; k++) {
        MX pos_x = x(0, k);
        MX pos_y = x(1, k);
        MX pos_z = x(2, k);
        MX vel_x = x(3, k);
        MX vel_y = x(4, k);
        MX vel_z = x(5, k);
        
        MX force_x = u(0, k);
        MX force_y = u(1, k);
        MX force_z = u(2, k);
        
        MX dt = t(k);
        MX safe_dt = fmax(dt, 1.0);
        
        MX next_pos_x = pos_x + vel_x * safe_dt;
        MX next_pos_y = pos_y + vel_y * safe_dt;
        MX next_pos_z = pos_z + vel_z * safe_dt;
        
        MX next_vel_x = vel_x + (force_x / mass) * safe_dt;
        MX next_vel_y = vel_y + (force_y / mass) * safe_dt;
        MX next_vel_z = vel_z + (force_z / mass) * safe_dt - 9.81 * safe_dt;
        
        equations.push_back(x(0, k+1) - next_pos_x);
        equations.push_back(x(1, k+1) - next_pos_y);
        equations.push_back(x(2, k+1) - next_pos_z);
        equations.push_back(x(3, k+1) - next_vel_x);
        equations.push_back(x(4, k+1) - next_vel_y);
        equations.push_back(x(5, k+1) - next_vel_z);
    }
    
    std::vector<MX> airspace_constraints;
    
    MX current_time = 0;
    for (int k = 0; k <= N; k++) {
        if (k > 0) {
            current_time = current_time + fmax(t(k-1), 1.0);
        }
        
        MX pos_x = x(0, k);
        MX pos_y = x(1, k);
        MX pos_z = x(2, k);
        
        airspace_constraints.push_back(pos_z - system.min_altitude);
        airspace_constraints.push_back(system.max_altitude - pos_z);
        
        MX speed_squared = pow(x(3, k), 2) + pow(x(4, k), 2) + pow(x(5, k), 2);
        MX speed = sqrt(speed_squared + 1e-6);
        airspace_constraints.push_back(system.max_speed - speed);
        
        for (const auto& constraint : constraints) {
            if (!constraint.is_active) continue;
            
            try {
                MX constraint_val = createDistanceConstraint(pos_x, pos_y, pos_z, constraint, current_time);
                airspace_constraints.push_back(constraint_val);
            } catch (...) {
                continue;
            }
        }
    }
    
    MX obj;
    if (objective == "fuel") {
        obj = 0;
        for (int k = 0; k < N; k++) {
            obj += fabs(u(0, k)) + fabs(u(1, k)) + fabs(u(2, k));
        }
    } else if (objective == "time") {
        obj = 0;
        for (int k = 0; k < N; k++) {
            obj += t(k);
        }
    } else if (objective == "balanced") {
        double time_weight = 0.3;
        double fuel_weight = 0.7;
        
        MX time_obj = 0;
        for (int k = 0; k < N; k++) {
            time_obj += t(k);
        }
        
        MX fuel_obj = 0;
        for (int k = 0; k < N; k++) {
            fuel_obj += fabs(u(0, k)) + fabs(u(1, k)) + fabs(u(2, k));
        }
        
        obj = time_weight * time_obj / 3600.0 + fuel_weight * fuel_obj / 1000000.0;
    } else {
        obj = 0;
        for (int k = 0; k < N; k++) {
            obj += pow(u(0, k), 2) + pow(u(1, k), 2) + pow(u(2, k), 2);
        }
    }
    
    std::vector<MX> vars = {reshape(x, -1, 1), reshape(u, -1, 1), t};
    MX opt_vars = vertcat(vars);
    
    std::vector<MX> all_constraints = equations;
    all_constraints.insert(all_constraints.end(), airspace_constraints.begin(), airspace_constraints.end());
    
    MX g = vertcat(all_constraints);
    
    MXDict nlp = {{"x", opt_vars}, {"f", obj}, {"g", g}};
    
    Dict solver_opts;
    solver_opts["ipopt.print_level"] = 0;
    solver_opts["ipopt.tol"] = 1e-4;
    solver_opts["ipopt.max_iter"] = 1000;
    solver_opts["ipopt.acceptable_tol"] = 1e-3;
    solver_opts["ipopt.hessian_approximation"] = "limited-memory";
    solver_opts["ipopt.mu_strategy"] = "adaptive";
    solver_opts["print_time"] = 0;
    solver_opts["ipopt.linear_solver"] = "mumps";
    
    return nlpsol("solver", opt_params.solver, nlp, solver_opts);
}

MX TrajectoryOptimizer::createDistanceConstraint(const MX& pos_x, const MX& pos_y, const MX& pos_z, 
                                            const AirspaceConstraint& constraint, const MX& time) {
    switch (constraint.constraint_type) {
        case AirspaceConstraintType::SPHERE: {
            double cx = constraint.parameters[0];
            double cy = constraint.parameters[1];
            double cz = constraint.parameters[2];
            double radius = constraint.parameters[3];
            
            MX dist_to_center = sqrt(pow(pos_x - cx, 2) + pow(pos_y - cy, 2) + pow(pos_z - cz, 2));
            
            MX constraint_value = dist_to_center - radius;
            
            double t_start = constraint.start_time;
            double t_end = constraint.end_time;
            double t_margin = 60.0;
            
            MX start_factor = 0.5 * (1 + tanh((time - t_start) / t_margin));
            MX end_factor = 0.5 * (1 + tanh((t_end - time) / t_margin));
            MX time_factor = start_factor * end_factor;
            
            return time_factor * constraint_value + (1.0 - time_factor) * MX(1000.0);
        }
        
        case AirspaceConstraintType::CYLINDER: {
            double cx = constraint.parameters[0];
            double cy = constraint.parameters[1];
            double min_z = constraint.parameters[2];
            double max_z = constraint.parameters[3];
            double radius = constraint.parameters[4];
            
            MX dist_to_axis = sqrt(pow(pos_x - cx, 2) + pow(pos_y - cy, 2));
            
            double smooth_factor = 10.0;
            MX in_height_range = (0.5 * (1 + tanh(smooth_factor * (pos_z - min_z)))) * 
                                (0.5 * (1 - tanh(smooth_factor * (pos_z - max_z))));
            
            MX radial_constraint = dist_to_axis - radius;
            
            MX constraint_value = in_height_range * radial_constraint + (1.0 - in_height_range) * MX(1000.0);
            
            double t_start = constraint.start_time;
            double t_end = constraint.end_time;
            double t_margin = 60.0;
            
            MX start_factor = 0.5 * (1 + tanh((time - t_start) / t_margin));
            MX end_factor = 0.5 * (1 + tanh((t_end - time) / t_margin));
            MX time_factor = start_factor * end_factor;
            
            return time_factor * constraint_value + (1.0 - time_factor) * MX(1000.0);
        }
        
        case AirspaceConstraintType::CORRIDOR: {
            double sx = constraint.parameters[0];
            double sy = constraint.parameters[1];
            double sz = constraint.parameters[2];
            double ex = constraint.parameters[3];
            double ey = constraint.parameters[4];
            double ez = constraint.parameters[5];
            double width = constraint.parameters[6];
            
            double dx = ex - sx;
            double dy = ey - sy;
            double dz = ez - sz;
            double len = sqrt(dx*dx + dy*dy + dz*dz);
            dx /= len; dy /= len; dz /= len;
            
            MX proj = (pos_x - sx) * dx + (pos_y - sy) * dy + (pos_z - sz) * dz;
            
            MX closest_x = sx + proj * dx;
            MX closest_y = sy + proj * dy;
            MX closest_z = sz + proj * dz;
            
            MX dist = sqrt(pow(pos_x - closest_x, 2) + pow(pos_y - closest_y, 2) + pow(pos_z - closest_z, 2));

            double smooth_factor = 10.0;
            
            MX in_corridor_length = (0.5 * (1 + tanh(smooth_factor * proj))) * 
                                   (0.5 * (1 - tanh(smooth_factor * (proj - len))));
            
            MX corridor_metric = width - dist;
            
            MX smooth_constraint = in_corridor_length * corridor_metric + 
                                  (1.0 - in_corridor_length) * MX(1.0);
                                  
            double t_start = constraint.start_time;
            double t_end = constraint.end_time;
            double t_margin = 60.0;
            
            MX start_factor = 0.5 * (1 + tanh((time - t_start) / t_margin));
            MX end_factor = 0.5 * (1 + tanh((t_end - time) / t_margin));
            MX time_factor = start_factor * end_factor;
            
            return time_factor * smooth_constraint + (1.0 - time_factor) * MX(1.0);
        }
        
        case AirspaceConstraintType::TIME_WINDOW: {
            double tx = constraint.parameters[0];
            double ty = constraint.parameters[1];
            double tz = constraint.parameters[2];
            double t_early = constraint.parameters[3];
            double t_late = constraint.parameters[4];
            
            MX dist_to_target = sqrt(pow(pos_x - tx, 2) + pow(pos_y - ty, 2) + pow(pos_z - tz, 2));
            
            double prox_factor = 50.0;
            MX proximity = 0.5 * (1 - tanh(prox_factor * (dist_to_target - 100.0)));
            
            MX time_window = (time - t_early) * (t_late - time);
            
            return proximity * time_window + (1.0 - proximity) * MX(1.0);
        }
        
        default:
            return MX(1000.0);
    }
}

void TrajectoryOptimizer::addConstraint(const AirspaceConstraint& constraint) {
    constraints.push_back(constraint);
}

void TrajectoryOptimizer::clearConstraints() {
    constraints.clear();
}

void TrajectoryOptimizer::enableConstraint(size_t index, bool enable) {
    if (index < constraints.size()) {
        constraints[index].is_active = enable;
    }
}

void TrajectoryOptimizer::extractResult4D(const DMDict& result, TrajectoryResult& traj_result) {
    DM x_opt = result.at("x");
    int N = opt_params.N;
    
    traj_result.objective_value = result.at("f").scalar();
    
    traj_result.positions.resize(N + 1);
    traj_result.velocities.resize(N + 1);
    traj_result.controls.resize(N);
    
    for (int k = 0; k <= N; k++) {
        double pos_x = x_opt(6*k).scalar();
        double pos_y = x_opt(6*k+1).scalar();
        double pos_z = x_opt(6*k+2).scalar();
        double vel_x = x_opt(6*k+3).scalar();
        double vel_y = x_opt(6*k+4).scalar();
        double vel_z = x_opt(6*k+5).scalar();
        
        traj_result.positions[k] = {pos_x, pos_y, pos_z};
        traj_result.velocities[k] = {vel_x, vel_y, vel_z};
        
        if (k < N) {
            double force_x = x_opt(6*(N+1) + 3*k).scalar();
            double force_y = x_opt(6*(N+1) + 3*k+1).scalar();
            double force_z = x_opt(6*(N+1) + 3*k+2).scalar();
            traj_result.controls[k] = {force_x, force_y, force_z};
        }
    }
}

TrajectoryResult TrajectoryOptimizer::optimize4D(
    const std::string& objective, 
    const std::map<std::string, double>& kwargs) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    int N = opt_params.N;
    
    MX x = MX::sym("x", 6, N+1);
    MX u = MX::sym("u", 3, N);
    MX t = MX::sym("t", N);
    
    std::vector<double> x0 = {0.0, 0.0, 3000.0, 0.0, 0.0, 0.0};
    std::vector<double> xf = {100000.0, 80000.0, 10000.0, 0.0, 0.0, 0.0};
    
    for (const auto& kv : kwargs) {
        if (kv.first == "x0") x0[0] = kv.second;
        else if (kv.first == "y0") x0[1] = kv.second;
        else if (kv.first == "z0") x0[2] = kv.second;
        else if (kv.first == "vx0") x0[3] = kv.second;
        else if (kv.first == "vy0") x0[4] = kv.second;
        else if (kv.first == "vz0") x0[5] = kv.second;
        else if (kv.first == "xf") xf[0] = kv.second;
        else if (kv.first == "yf") xf[1] = kv.second;
        else if (kv.first == "zf") xf[2] = kv.second;
        else if (kv.first == "vxf") xf[3] = kv.second;
        else if (kv.first == "vyf") xf[4] = kv.second;
        else if (kv.first == "vzf") xf[5] = kv.second;
    }
    
    Function solver = createSolver4D(x, u, t, objective);
    
    int n_states = 6 * (N + 1);
    int n_controls = 3 * N;
    int n_times = N;
    
    std::vector<double> lbx(n_states + n_controls + n_times, -std::numeric_limits<double>::infinity());
    std::vector<double> ubx(n_states + n_controls + n_times, std::numeric_limits<double>::infinity());
    
    for (int i = n_states; i < n_states + n_controls; i += 3) {
        lbx[i] = -system.max_force_xy;
        ubx[i] = system.max_force_xy;
        lbx[i+1] = -system.max_force_xy;
        ubx[i+1] = system.max_force_xy;
        lbx[i+2] = -system.max_force_z;
        ubx[i+2] = system.max_force_z;
    }
    
    for (int i = n_states + n_controls; i < n_states + n_controls + n_times; i++) {
        lbx[i] = 10.0;
        ubx[i] = 300.0;
    }
    
    for (int i = 0; i < 6; i++) {
        lbx[i] = ubx[i] = x0[i];
    }
    
    for (int i = 0; i < 6; i++) {
        lbx[n_states - 6 + i] = ubx[n_states - 6 + i] = xf[i];
    }
    
    int n_dynamics_constraints = 6 * N;
    
    int active_constraints = 0;
    for (const auto& constraint : constraints) {
        if (constraint.is_active) {
            active_constraints++;
        }
    }
    
    int n_airspace_constraints = 3 * (N + 1) + active_constraints * (N + 1);
    int total_constraints = n_dynamics_constraints + n_airspace_constraints;
    
    std::vector<double> lbg(n_dynamics_constraints, 0.0);
    std::vector<double> ubg(n_dynamics_constraints, 0.0);
    
    for (int i = 0; i < n_airspace_constraints; i++) {
        lbg.push_back(0.0);
        ubg.push_back(INFINITY);
    }
    
    std::vector<double> x0_guess(n_states + n_controls + n_times, 0.0);
    
    for (int k = 0; k <= N; k++) {
        double t = static_cast<double>(k) / N;
        for (int i = 0; i < 3; i++) {
            x0_guess[6*k + i] = x0[i] + t * (xf[i] - x0[i]);
        }
        
        if (k < N) {
            double estimated_flight_time = 3600.0;
            for (int i = 0; i < 3; i++) {
                x0_guess[6*k + 3 + i] = (xf[i] - x0[i]) / estimated_flight_time;
            }
        }
    }
    
    for (int k = 0; k < N; k++) {
        x0_guess[n_states + 3*k] = (xf[0] - x0[0]) > 0 ? 10000.0 : -10000.0;
        x0_guess[n_states + 3*k + 1] = (xf[1] - x0[1]) > 0 ? 10000.0 : -10000.0;
        x0_guess[n_states + 3*k + 2] = (xf[2] - x0[2]) > 0 ? 5000.0 : -5000.0;
    }
    
    double estimated_flight_time = std::sqrt(
        std::pow(xf[0] - x0[0], 2) + 
        std::pow(xf[1] - x0[1], 2) + 
        std::pow(xf[2] - x0[2], 2)
    ) / (0.8 * system.max_speed);

    double dt_guess = std::max(10.0, std::min(300.0, estimated_flight_time / N));

    for (int k = 0; k < N; k++) {
        x0_guess[n_states + n_controls + k] = dt_guess;
    }
    
    DMDict result;
    try {
        if (containsNaN(x0_guess)) {
            std::fill(x0_guess.begin(), x0_guess.end(), 0.0);
            for (int k = 0; k < N; k++) {
                x0_guess[n_states + n_controls + k] = dt_guess;
            }
        }
        
        result = solver(DMDict{
            {"x0", x0_guess},
            {"lbx", lbx},
            {"ubx", ubx},
            {"lbg", lbg},
            {"ubg", ubg}
        });
    } catch (std::exception& e) {
        std::cerr << "Оптимизация не удалась: " << e.what() << std::endl;
        std::cerr << "Пробуем с упрощенными настройками..." << std::endl;
        
        TrajectoryResult debug_result = debug4D(objective, kwargs);
        if (debug_result.computation_time > 0) {
            std::cout << "Траектория успешно создана с использованием упрощенного подхода." << std::endl;
            return debug_result;
        }
        
        TrajectoryResult failed_result;
        failed_result.computation_time = -1.0;
        return failed_result;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    TrajectoryResult traj_result;
    extractResult4D(result, traj_result);
    
    traj_result.times.resize(N + 1, 0.0);
    DM t_opt = result.at("x").nz(Slice(n_states + n_controls, n_states + n_controls + n_times));
    
    for (int k = 1; k <= N; k++) {
        traj_result.times[k] = traj_result.times[k-1] + t_opt(k-1).scalar();
    }
    
    traj_result.computation_time = duration.count() / 1000.0;
    traj_result.total_flight_time = traj_result.times.back();
    
    traj_result.fuel_consumption = 0.0;
    for (size_t k = 0; k < traj_result.controls.size(); k++) {
        double force_magnitude = sqrt(pow(traj_result.controls[k][0], 2) + 
                                      pow(traj_result.controls[k][1], 2) + 
                                      pow(traj_result.controls[k][2], 2));
        
        double dt = k < N ? t_opt(k).scalar() : 0.0;
        traj_result.fuel_consumption += force_magnitude * dt * 0.001;
    }
    
    return traj_result;
}

void TrajectoryOptimizer::visualizeTrajectory(const TrajectoryResult& result) const {
    if (result.positions.empty()) {
        std::cout << "Траектория пуста или не была рассчитана" << std::endl;
        return;
    }
    
    std::cout << "=== Результаты оптимизации траектории ===" << std::endl;
    std::cout << "Время вычислений: " << result.computation_time << " сек." << std::endl;
    if (!result.times.empty() && result.times.size() >= result.positions.size()) {
        std::cout << "Общее время полета: " << result.total_flight_time / 60.0 << " мин." << std::endl;
    } else {
        std::cout << "Общее время полета: Данные о времени недоступны" << std::endl;
    }
    
    std::cout << "Расход топлива: " << result.fuel_consumption << " единиц" << std::endl;
    std::cout << "Значение целевой функции: " << result.objective_value << std::endl;
    std::cout << "------------------" << std::endl;
    
    int num_points = result.positions.size();
    int step = num_points <= 10 ? 1 : num_points / 10;
    
    std::cout << "Ключевые точки траектории (всего " << num_points << " точек):" << std::endl;
    for (size_t k = 0; k < result.positions.size(); k += step) {
        if (!result.times.empty() && k < result.times.size()) {
            std::cout << "Точка " << k << ": Время = " << result.times[k] / 60.0 << " мин.";
        } else {
            std::cout << "Точка " << k << ": Время = N/A";
        }
        
        if (result.positions[k].size() >= 3) {
            std::cout << ", Позиция = (" 
                      << result.positions[k][0] / 1000.0 << " км, " 
                      << result.positions[k][1] / 1000.0 << " км, " 
                      << result.positions[k][2] << " м)";
        } else if (result.positions[k].size() >= 2) {
            std::cout << ", Позиция = (" 
                      << result.positions[k][0] << ", " 
                      << result.positions[k][1] << ")";
        }
        
        if (k < result.velocities.size()) {
            if (result.velocities[k].size() >= 3) {
                std::cout << ", Скорость = (" 
                          << result.velocities[k][0] << ", " 
                          << result.velocities[k][1] << ", " 
                          << result.velocities[k][2] << ") м/с";
            } else if (result.velocities[k].size() >= 2) {
                std::cout << ", Скорость = (" 
                          << result.velocities[k][0] << ", " 
                          << result.velocities[k][1] << ")";
            }
        }
                  
        if (k < result.controls.size()) {
            if (result.controls[k].size() >= 3) {
                std::cout << std::endl << "    Управление = (" 
                          << result.controls[k][0] / 1000.0 << ", " 
                          << result.controls[k][1] / 1000.0 << ", " 
                          << result.controls[k][2] / 1000.0 << ") кН";
            } else if (result.controls[k].size() >= 2) {
                std::cout << std::endl << "    Управление = (" 
                          << result.controls[k][0] << ", " 
                          << result.controls[k][1] << ")";
            }
        }
        std::cout << std::endl;
    }
    
    if (num_points > 1 && (num_points - 1) % step != 0) {
        size_t k = num_points - 1;
        
        if (!result.times.empty() && k < result.times.size()) {
            std::cout << "Точка " << k << " (конечная): Время = " << result.times[k] / 60.0 << " мин.";
        } else {
            std::cout << "Точка " << k << " (конечная): Время = N/A";
        }
        
        if (result.positions[k].size() >= 3) {
            std::cout << ", Позиция = (" 
                      << result.positions[k][0] / 1000.0 << " км, "
                      << result.positions[k][1] / 1000.0 << " км, " 
                      << result.positions[k][2] << " м)";
        } else if (result.positions[k].size() >= 2) {
            std::cout << ", Позиция = (" 
                      << result.positions[k][0] << ", "
                      << result.positions[k][1] << ")";
        }
        
        if (k < result.velocities.size()) {
            if (result.velocities[k].size() >= 3) {
                std::cout << ", Скорость = (" 
                          << result.velocities[k][0] << ", " 
                          << result.velocities[k][1] << ", " 
                          << result.velocities[k][2] << ") м/с";
            } else if (result.velocities[k].size() >= 2) {
                std::cout << ", Скорость = (" 
                          << result.velocities[k][0] << ", " 
                          << result.velocities[k][1] << ")";
            }
        }
        std::cout << std::endl;
    }
}

void TrajectoryOptimizer::exportTrajectory(const TrajectoryResult& result, const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка при открытии файла " << filename << " для записи." << std::endl;
        return;
    }
    
    file << "point_id,time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z";
    if (!result.controls.empty()) {
        file << ",force_x,force_y,force_z";
    }
    file << std::endl;
    
    for (size_t k = 0; k < result.positions.size(); k++) {
        file << k << "," << result.times[k] << ","
             << result.positions[k][0] << "," << result.positions[k][1] << "," << result.positions[k][2] << ","
             << result.velocities[k][0] << "," << result.velocities[k][1] << "," << result.velocities[k][2];
             
        if (k < result.controls.size()) {
            file << "," << result.controls[k][0] << "," << result.controls[k][1] << "," << result.controls[k][2];
        } else if (!result.controls.empty()) {
            file << ",,,"; // Пустые значения для последней точки
        }
        
        file << std::endl;
    }
    
    file.close();
    std::cout << "Траектория экспортирована в файл " << filename << std::endl;
}

bool TrajectoryOptimizer::containsNaN(const std::vector<double>& vec) const {
    for (const double& val : vec) {
        if (std::isnan(val) || std::isinf(val)) {
            return true;
        }
    }
    return false;
}

TrajectoryResult TrajectoryOptimizer::debug4D(const std::string& objective,
                        const std::map<std::string, double>& kwargs) {
    
    int N = 10; 
    
    OptimizationParams saved_params = opt_params;
    opt_params.N = N;
    opt_params.dt = 300.0;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    MX x = MX::sym("x", 6, N+1);
    MX u = MX::sym("u", 3, N);
    
    std::vector<double> x0 = {0.0, 0.0, 3000.0, 0.0, 0.0, 0.0};
    std::vector<double> xf = {100000.0, 80000.0, 5000.0, 0.0, 0.0, 0.0};
    
    for (const auto& kv : kwargs) {
        if (kv.first == "x0") x0[0] = kv.second;
        else if (kv.first == "y0") x0[1] = kv.second;
        else if (kv.first == "z0") x0[2] = kv.second;
        else if (kv.first == "xf") xf[0] = kv.second;
        else if (kv.first == "yf") xf[1] = kv.second;
        else if (kv.first == "zf") xf[2] = kv.second;
    }
    
    std::vector<MX> equations;
    for (int k = 0; k < N; k++) {
        MX pos_x = x(0, k);
        MX pos_y = x(1, k);
        MX pos_z = x(2, k);
        MX vel_x = x(3, k);
        MX vel_y = x(4, k);
        MX vel_z = x(5, k);
        
        MX force_x = u(0, k);
        MX force_y = u(1, k);
        MX force_z = u(2, k);
        
        double dt = opt_params.dt;
        double mass = system.mass;
        if (fabs(mass) < 1e-6) mass = 1e-6;
        
        MX next_pos_x = pos_x + vel_x * dt;
        MX next_pos_y = pos_y + vel_y * dt;
        MX next_pos_z = pos_z + vel_z * dt;
        
        MX next_vel_x = vel_x + (force_x / mass) * dt;
        MX next_vel_y = vel_y + (force_y / mass) * dt;
        MX next_vel_z = vel_z + (force_z / mass) * dt - 9.81 * dt;
        
        equations.push_back(x(0, k+1) - next_pos_x);
        equations.push_back(x(1, k+1) - next_pos_y);
        equations.push_back(x(2, k+1) - next_pos_z);
        equations.push_back(x(3, k+1) - next_vel_x);
        equations.push_back(x(4, k+1) - next_vel_y);
        equations.push_back(x(5, k+1) - next_vel_z);
    }
    
    MX obj;
    if (objective == "fuel" || objective == "energy") {
        obj = 0;
        for (int k = 0; k < N; k++) {
            obj += pow(u(0, k), 2) + pow(u(1, k), 2) + pow(u(2, k), 2);
        }
    } else {
        obj = N * opt_params.dt;
    }
    
    std::vector<MX> vars = {reshape(x, -1, 1), reshape(u, -1, 1)};
    MX opt_vars = vertcat(vars);
    MX g = vertcat(equations);
    MXDict nlp = {{"x", opt_vars}, {"f", obj}, {"g", g}};
    
    Dict solver_opts;
    solver_opts["ipopt.print_level"] = 0;
    solver_opts["ipopt.tol"] = 1e-3;       
    solver_opts["ipopt.max_iter"] = 500;   
    solver_opts["ipopt.acceptable_tol"] = 1e-2;
    solver_opts["print_time"] = 0;

    Function solver = nlpsol("solver", opt_params.solver, nlp, solver_opts);
    
    int n_states = 6 * (N + 1);
    int n_controls = 3 * N;
    std::vector<double> lbx(n_states + n_controls, -std::numeric_limits<double>::infinity());
    std::vector<double> ubx(n_states + n_controls, std::numeric_limits<double>::infinity());
    
    for (int i = n_states; i < n_states + n_controls; i += 3) {
        lbx[i] = lbx[i+1] = -system.max_force_xy;
        ubx[i] = ubx[i+1] = system.max_force_xy;
        lbx[i+2] = -system.max_force_z;
        ubx[i+2] = system.max_force_z;
    }
    
    for (int i = 0; i < 6; i++) {
        lbx[i] = ubx[i] = x0[i];
    }
    
    for (int i = 0; i < 6; i++) {
        lbx[n_states - 6 + i] = ubx[n_states - 6 + i] = xf[i];
    }
    
    std::vector<double> x0_guess(n_states + n_controls, 0.0);
    
    for (int k = 0; k <= N; k++) {
        double t = static_cast<double>(k) / N;
        for (int i = 0; i < 3; i++) {
            x0_guess[6*k + i] = x0[i] + t * (xf[i] - x0[i]);
        }
        
        if (k < N) {
            double estimated_flight_time = 3600.0;
            for (int i = 0; i < 3; i++) {
                x0_guess[6*k + 3 + i] = (xf[i] - x0[i]) / estimated_flight_time;
            }
        }
    }
    
    for (int k = 0; k < N; k++) {
        x0_guess[n_states + 3*k] = (xf[0] - x0[0]) > 0 ? 10000.0 : -10000.0;
        x0_guess[n_states + 3*k + 1] = (xf[1] - x0[1]) > 0 ? 10000.0 : -10000.0;
        x0_guess[n_states + 3*k + 2] = (xf[2] - x0[2]) > 0 ? 5000.0 : -5000.0;
    }
    
    double estimated_flight_time = std::sqrt(
        std::pow(xf[0] - x0[0], 2) + 
        std::pow(xf[1] - x0[1], 2) + 
        std::pow(xf[2] - x0[2], 2)
    ) / (0.8 * system.max_speed);

    double dt_guess = std::max(10.0, std::min(300.0, estimated_flight_time / N));

    for (int k = 0; k < N; k++) {
        x0_guess[n_states + n_controls + k] = dt_guess;
    }
    
    DMDict result;
    try {
        if (containsNaN(x0_guess)) {
            std::fill(x0_guess.begin(), x0_guess.end(), 0.0);
            for (int k = 0; k < N; k++) {
                x0_guess[n_states + n_controls + k] = dt_guess;
            }
        }
        
        result = solver(DMDict{
            {"x0", x0_guess},
            {"lbx", lbx},
            {"ubx", ubx},
            {"lbg", lbg},
            {"ubg", ubg}
        });
    } catch (std::exception& e) {
        std::cerr << "Оптимизация не удалась: " << e.what() << std::endl;
        std::cerr << "Пробуем с упрощенными настройками..." << std::endl;
        
        TrajectoryResult debug_result = debug4D(objective, kwargs);
        if (debug_result.computation_time > 0) {
            std::cout << "Траектория успешно создана с использованием упрощенного подхода." << std::endl;
            return debug_result;
        }
        
        TrajectoryResult failed_result;
        failed_result.computation_time = -1.0;
        return failed_result;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    TrajectoryResult traj_result;
    extractResult4D(result, traj_result);
    
    traj_result.times.resize(N + 1, 0.0);
    DM t_opt = result.at("x").nz(Slice(n_states + n_controls, n_states + n_controls + n_times));
    
    for (int k = 1; k <= N; k++) {
        traj_result.times[k] = traj_result.times[k-1] + t_opt(k-1).scalar();
    }
    
    traj_result.computation_time = duration.count() / 1000.0;
    traj_result.total_flight_time = traj_result.times.back();
    
    traj_result.fuel_consumption = 0.0;
    for (size_t k = 0; k < traj_result.controls.size(); k++) {
        double force_magnitude = sqrt(pow(traj_result.controls[k][0], 2) + 
                                      pow(traj_result.controls[k][1], 2) + 
                                      pow(traj_result.controls[k][2], 2));
        
        double dt = k < N ? t_opt(k).scalar() : 0.0;
        traj_result.fuel_consumption += force_magnitude * dt * 0.001;
    }
    
    return traj_result;
}


// Метод для визуализации результатов 4D-траектории
void TrajectoryOptimizer::visualizeTrajectory(const TrajectoryResult& result) const {
    if (result.positions.empty()) {
        std::cout << "Траектория пуста или не была рассчитана" << std::endl;
        return;
    }
    
    std::cout << "=== Результаты оптимизации траектории ===" << std::endl;
    std::cout << "Время вычислений: " << result.computation_time << " сек." << std::endl;
    if (!result.times.empty() && result.times.size() >= result.positions.size()) {
        std::cout << "Общее время полета: " << result.total_flight_time / 60.0 << " мин." << std::endl;
    } else {
        std::cout << "Общее время полета: Данные о времени недоступны" << std::endl;
    }
    
    std::cout << "Расход топлива: " << result.fuel_consumption << " единиц" << std::endl;
    std::cout << "Значение целевой функции: " << result.objective_value << std::endl;
    std::cout << "------------------" << std::endl;
    
    // Выводим только ключевые точки траектории для компактности
    int num_points = result.positions.size();
    int step = num_points <= 10 ? 1 : num_points / 10;
    
    std::cout << "Ключевые точки траектории (всего " << num_points << " точек):" << std::endl;
    for (size_t k = 0; k < result.positions.size(); k += step) {
        if (!result.times.empty() && k < result.times.size()) {
            std::cout << "Точка " << k << ": Время = " << result.times[k] / 60.0 << " мин.";
        } else {
            std::cout << "Точка " << k << ": Время = N/A";
        }
        
        if (result.positions[k].size() >= 3) {
            std::cout << ", Позиция = (" 
                      << result.positions[k][0] / 1000.0 << " км, " 
                      << result.positions[k][1] / 1000.0 << " км, " 
                      << result.positions[k][2] << " м)";
        } else if (result.positions[k].size() >= 2) {
            // 2D trajectory
            std::cout << ", Позиция = (" 
                      << result.positions[k][0] << ", " 
                      << result.positions[k][1] << ")";
        }
        
        if (k < result.velocities.size()) {
            if (result.velocities[k].size() >= 3) {
                std::cout << ", Скорость = (" 
                          << result.velocities[k][0] << ", " 
                          << result.velocities[k][1] << ", " 
                          << result.velocities[k][2] << ") м/с";
            } else if (result.velocities[k].size() >= 2) {
                std::cout << ", Скорость = (" 
                          << result.velocities[k][0] << ", " 
                          << result.velocities[k][1] << ")";
            }
        }
                  
        if (k < result.controls.size()) {
            if (result.controls[k].size() >= 3) {
                std::cout << std::endl << "    Управление = (" 
                          << result.controls[k][0] / 1000.0 << ", " 
                          << result.controls[k][1] / 1000.0 << ", " 
                          << result.controls[k][2] / 1000.0 << ") кН";
            } else if (result.controls[k].size() >= 2) {
                std::cout << std::endl << "    Управление = (" 
                          << result.controls[k][0] << ", " 
                          << result.controls[k][1] << ")";
            }
        }
        std::cout << std::endl;
    }
    
    // Отображаем последнюю точку всегда
    if (num_points > 1 && (num_points - 1) % step != 0) {
        size_t k = num_points - 1;
        
        if (!result.times.empty() && k < result.times.size()) {
            std::cout << "Точка " << k << " (конечная): Время = " << result.times[k] / 60.0 << " мин.";
        } else {
            std::cout << "Точка " << k << " (конечная): Время = N/A";
        }
        
        if (result.positions[k].size() >= 3) {
            std::cout << ", Позиция = (" 
                      << result.positions[k][0] / 1000.0 << " км, "
                      << result.positions[k][1] / 1000.0 << " км, " 
                      << result.positions[k][2] << " м)";
        } else if (result.positions[k].size() >= 2) {
            std::cout << ", Позиция = (" 
                      << result.positions[k][0] << ", "
                      << result.positions[k][1] << ")";
        }
        
        if (k < result.velocities.size()) {
            if (result.velocities[k].size() >= 3) {
                std::cout << ", Скорость = (" 
                          << result.velocities[k][0] << ", " 
                          << result.velocities[k][1] << ", " 
                          << result.velocities[k][2] << ") м/с";
            } else if (result.velocities[k].size() >= 2) {
                std::cout << ", Скорость = (" 
                          << result.velocities[k][0] << ", " 
                          << result.velocities[k][1] << ")";
            }
        }
        std::cout << std::endl;
    }
}

// Экспорт траектории в CSV файл для последующего анализа
void TrajectoryOptimizer::exportTrajectory(const TrajectoryResult& result, const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка при открытии файла " << filename << " для записи." << std::endl;
        return;
    }
    
    // Записываем заголовок
    file << "point_id,time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z";
    if (!result.controls.empty()) {
        file << ",force_x,force_y,force_z";
    }
    file << std::endl;
    
    // Записываем данные
    for (size_t k = 0; k < result.positions.size(); k++) {
        file << k << "," << result.times[k] << ","
             << result.positions[k][0] << "," << result.positions[k][1] << "," << result.positions[k][2] << ","
             << result.velocities[k][0] << "," << result.velocities[k][1] << "," << result.velocities[k][2];
             
        if (k < result.controls.size()) {
            file << "," << result.controls[k][0] << "," << result.controls[k][1] << "," << result.controls[k][2];
        } else if (!result.controls.empty()) {
            file << ",,,"; // Пустые значения для последней точки
        }
        
        file << std::endl;
    }
    
    file.close();
    std::cout << "Траектория экспортирована в файл " << filename << std::endl;
}

bool TrajectoryOptimizer::containsNaN(const std::vector<double>& vec) const {
    for (const double& val : vec) {
        if (std::isnan(val) || std::isinf(val)) {
            return true;
        }
    }
    return false;
}

TrajectoryResult TrajectoryOptimizer::debug4D(const std::string& objective,
                        const std::map<std::string, double>& kwargs) {
    
    int N = 10; 
    
    OptimizationParams saved_params = opt_params;
    opt_params.N = N;
    opt_params.dt = 300.0;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    MX x = MX::sym("x", 6, N+1);
    MX u = MX::sym("u", 3, N);
    
    std::vector<double> x0 = {0.0, 0.0, 3000.0, 0.0, 0.0, 0.0};
    std::vector<double> xf = {100000.0, 80000.0, 5000.0, 0.0, 0.0, 0.0};
    
    for (const auto& kv : kwargs) {
        if (kv.first == "x0") x0[0] = kv.second;
        else if (kv.first == "y0") x0[1] = kv.second;
        else if (kv.first == "z0") x0[2] = kv.second;
        else if (kv.first == "xf") xf[0] = kv.second;
        else if (kv.first == "yf") xf[1] = kv.second;
        else if (kv.first == "zf") xf[2] = kv.second;
    }
    
    std::vector<MX> equations;
    for (int k = 0; k < N; k++) {
        MX pos_x = x(0, k);
        MX pos_y = x(1, k);
        MX pos_z = x(2, k);
        MX vel_x = x(3, k);
        MX vel_y = x(4, k);
        MX vel_z = x(5, k);
        
        MX force_x = u(0, k);
        MX force_y = u(1, k);
        MX force_z = u(2, k);
        
        double dt = opt_params.dt;
        double mass = system.mass;
        if (fabs(mass) < 1e-6) mass = 1e-6;
        
        MX next_pos_x = pos_x + vel_x * dt;
        MX next_pos_y = pos_y + vel_y * dt;
        MX next_pos_z = pos_z + vel_z * dt;
        
        MX next_vel_x = vel_x + (force_x / mass) * dt;
        MX next_vel_y = vel_y + (force_y / mass) * dt;
        MX next_vel_z = vel_z + (force_z / mass) * dt - 9.81 * dt;
        
        equations.push_back(x(0, k+1) - next_pos_x);
        equations.push_back(x(1, k+1) - next_pos_y);
        equations.push_back(x(2, k+1) - next_pos_z);
        equations.push_back(x(3, k+1) - next_vel_x);
        equations.push_back(x(4, k+1) - next_vel_y);
        equations.push_back(x(5, k+1) - next_vel_z);
    }
    
    MX obj;
    if (objective == "fuel" || objective == "energy") {
        obj = 0;
        for (int k = 0; k < N; k++) {
            obj += pow(u(0, k), 2) + pow(u(1, k), 2) + pow(u(2, k), 2);
        }
    } else {
        obj = N * opt_params.dt;
    }
    
    std::vector<MX> vars = {reshape(x, -1, 1), reshape(u, -1, 1)};
    MX opt_vars = vertcat(vars);
    MX g = vertcat(equations);
    MXDict nlp = {{"x", opt_vars}, {"f", obj}, {"g", g}};
    
    Dict solver_opts;
    solver_opts["ipopt.print_level"] = 0;
    solver_opts["ipopt.tol"] = 1e-3;       
    solver_opts["ipopt.max_iter"] = 500;   
    solver_opts["ipopt.acceptable_tol"] = 1e-2;
    solver_opts["print_time"] = 0;

    Function solver = nlpsol("solver", opt_params.solver, nlp, solver_opts);
    
    int n_states = 6 * (N + 1);
    int n_controls = 3 * N;
    std::vector<double> lbx(n_states + n_controls, -std::numeric_limits<double>::infinity());
    std::vector<double> ubx(n_states + n_controls, std::numeric_limits<double>::infinity());
    
    for (int i = n_states; i < n_states + n_controls; i += 3) {
        lbx[i] = lbx[i+1] = -system.max_force_xy;
        ubx[i] = ubx[i+1] = system.max_force_xy;
        lbx[i+2] = -system.max_force_z;
        ubx[i+2] = system.max_force_z;
    }
    
    for (int i = 0; i < 6; i++) {
        lbx[i] = ubx[i] = x0[i];
    }
    
    for (int i = 0; i < 6; i++) {
        lbx[n_states - 6 + i] = ubx[n_states - 6 + i] = xf[i];
    }
    
    std::vector<double> x0_guess(n_states + n_controls, 0.0);
    
    for (int k = 0; k <= N; k++) {
        double t = static_cast<double>(k) / N;
        for (int i = 0; i < 3; i++) {
            x0_guess[6*k + i] = x0[i] + t * (xf[i] - x0[i]);
        }
        
        if (k < N) {
            double estimated_flight_time = 3600.0;
            for (int i = 0; i < 3; i++) {
                x0_guess[6*k + 3 + i] = (xf[i] - x0[i]) / estimated_flight_time;
            }
        }
    }
    
    for (int k = 0; k < N; k++) {
        x0_guess[n_states + 3*k] = (xf[0] - x0[0]) > 0 ? 10000.0 : -10000.0;
        x0_guess[n_states + 3*k + 1] = (xf[1] - x0[1]) > 0 ? 10000.0 : -10000.0;
        x0_guess[n_states + 3*k + 2] = (xf[2] - x0[2]) > 0 ? 5000.0 : -5000.0;
    }
    
    double estimated_flight_time = std::sqrt(
        std::pow(xf[0] - x0[0], 2) + 
        std::pow(xf[1] - x0[1], 2) + 
        std::pow(xf[2] - x0[2], 2)
    ) / (0.8 * system.max_speed);

    double dt_guess = std::max(10.0, std::min(300.0, estimated_flight_time / N));

    for (int k = 0; k < N; k++) {
        x0_guess[n_states + n_controls + k] = dt_guess;
    }
    
    DMDict result;
    try {
        if (containsNaN(x0_guess)) {
            std::fill(x0_guess.begin(), x0_guess.end(), 0.0);
            for (int k = 0; k < N; k++) {
                x0_guess[n_states + n_controls + k] = dt_guess;
            }
        }
        
        result = solver(DMDict{
            {"x0", x0_guess},
            {"lbx", lbx},
            {"ubx", ubx},
            {"lbg", lbg},
            {"ubg", ubg}
        });
    } catch (std::exception& e) {
        std::cerr << "Оптимизация не удалась: " << e.what() << std::endl;
        std::cerr << "Пробуем с упрощенными настройками..." << std::endl;
        
        TrajectoryResult debug_result = debug4D(objective, kwargs);
        if (debug_result.computation_time > 0) {
            std::cout << "Траектория успешно создана с использованием упрощенного подхода." << std::endl;
            return debug_result;
        }
        
        TrajectoryResult failed_result;
        failed_result.computation_time = -1.0;
        return failed_result;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    TrajectoryResult traj_result;
    extractResult4D(result, traj_result);
    
    traj_result.times.resize(N + 1, 0.0);
    DM t_opt = result.at("x").nz(Slice(n_states + n_controls, n_states + n_controls + n_times));
    
    for (int k = 1; k <= N; k++) {
        traj_result.times[k] = traj_result.times[k-1] + t_opt(k-1).scalar();
    }
    
    traj_result.computation_time = duration.count() / 1000.0;
    traj_result.total_flight_time = traj_result.times.back();
    
    traj_result.fuel_consumption = 0.0;
    for (size_t k = 0; k < traj_result.controls.size(); k++) {
        double force_magnitude = sqrt(pow(traj_result.controls[k][0], 2) + 
                                      pow(traj_result.controls[k][1], 2) + 
                                      pow(traj_result.controls[k][2], 2));
        
        double dt = k < N ? t_opt(k).scalar() : 0.0;
        traj_result.fuel_consumption += force_magnitude * dt * 0.001;
    }
    
    return traj_result;
}