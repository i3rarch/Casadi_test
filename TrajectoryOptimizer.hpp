#pragma once

#include <casadi/casadi.hpp>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <memory>
#include <fstream>
#include <chrono>  
#include <limits>  

using namespace casadi;

// Структура для хранения параметров системы
struct SystemParams {
    double mass;        // масса воздушного судна
    double max_force_xy;// максимальное управляющее воздействие в горизонтальной плоскости
    double max_force_z; // максимальное управляющее воздействие по вертикали
    double max_speed;   // максимальная скорость
    double min_altitude;// минимальная высота
    double max_altitude;// максимальная высота
    
    SystemParams(double m = 80000.0, double fxy = 200000.0, double fz = 100000.0, 
                 double v_max = 250.0, double alt_min = 3000.0, double alt_max = 12000.0) 
        : mass(m), max_force_xy(fxy), max_force_z(fz), max_speed(v_max),
          min_altitude(alt_min), max_altitude(alt_max) {}
};

// Структура для хранения параметров оптимизации
struct OptimizationParams {
    int N;              // количество шагов горизонта
    double dt;          // шаг по времени
    std::string solver; // тип решателя ('ipopt', 'sqpmethod', 'scpgen')
    
    OptimizationParams(int steps = 50, double step_size = 0.1, 
                     const std::string& slv = "ipopt") 
        : N(steps), dt(step_size), solver(slv) {}
};

// Структура для хранения результатов оптимизации
struct TrajectoryResult {
    std::vector<std::vector<double>> positions;  // координаты x,y,z в каждый момент времени
    std::vector<std::vector<double>> velocities; // скорости по x,y,z в каждый момент времени
    std::vector<std::vector<double>> controls;   // управления по x,y,z в каждый момент времени
    std::vector<double> times;                   // время в каждой точке маршрута
    double objective_value;                      // значение целевой функции
    double computation_time;                     // время вычислений
    double total_flight_time;                    // общее время полета
    double fuel_consumption;                     // расход топлива
    
    TrajectoryResult() : objective_value(0.0), computation_time(0.0),
                         total_flight_time(0.0), fuel_consumption(0.0) {}
};

// Переименовываем ConstraintType в AirspaceConstraintType чтобы избежать конфликта с casadi::ConstraintType
enum class AirspaceConstraintType {
    SPHERE,         // Сферическая область (гроза, опасная зона)
    CYLINDER,       // Цилиндрическая область (временно закрытая зона)
    CORRIDOR,       // Коридор для полета
    TIME_WINDOW     // Временное окно для прибытия
};

// Структура для ограничений
struct AirspaceConstraint {
    AirspaceConstraintType constraint_type; // переименованное поле
    std::vector<double> parameters; // Параметры ограничения (зависят от типа)
    double start_time;              // Время начала действия ограничения
    double end_time;                // Время окончания действия ограничения
    bool is_active;                 // Активно ли ограничение
    std::string description;        // Описание ограничения
    
    // Конструктор для сферического ограничения (грозовой фронт, запретная зона)
    static AirspaceConstraint createSphere(
        double center_x, double center_y, double center_z, double radius,
        double t_start = 0.0, double t_end = INFINITY, const std::string& desc = "Sphere constraint") {
        
        AirspaceConstraint constraint;
        constraint.constraint_type = AirspaceConstraintType::SPHERE;
        constraint.parameters = {center_x, center_y, center_z, radius};
        constraint.start_time = t_start;
        constraint.end_time = t_end;
        constraint.is_active = true;
        constraint.description = desc;
        return constraint;
    }
    
    // Конструктор для цилиндрического ограничения (закрытая зона на определенной высоте)
    static AirspaceConstraint createCylinder(
        double center_x, double center_y, double min_z, double max_z, double radius,
        double t_start = 0.0, double t_end = INFINITY, const std::string& desc = "Cylinder constraint") {
        
        AirspaceConstraint constraint;
        constraint.constraint_type = AirspaceConstraintType::CYLINDER;
        constraint.parameters = {center_x, center_y, min_z, max_z, radius};
        constraint.start_time = t_start;
        constraint.end_time = t_end;
        constraint.is_active = true;
        constraint.description = desc;
        return constraint;
    }
    
    // Конструктор для коридора (обязательный маршрут следования)
    static AirspaceConstraint createCorridor(
        double start_x, double start_y, double start_z,
        double end_x, double end_y, double end_z, double width,
        double t_start = 0.0, double t_end = INFINITY, const std::string& desc = "Corridor constraint") {
        
        AirspaceConstraint constraint;
        constraint.constraint_type = AirspaceConstraintType::CORRIDOR;
        constraint.parameters = {start_x, start_y, start_z, end_x, end_y, end_z, width};
        constraint.start_time = t_start;
        constraint.end_time = t_end;
        constraint.is_active = true;
        constraint.description = desc;
        return constraint;
    }
    
    // Конструктор для временного окна (ограничение по времени прибытия)
    static AirspaceConstraint createTimeWindow(
        double target_x, double target_y, double target_z,
        double earliest_arrival, double latest_arrival,
        const std::string& desc = "Time window constraint") {
        
        AirspaceConstraint constraint;
        constraint.constraint_type = AirspaceConstraintType::TIME_WINDOW;
        constraint.parameters = {target_x, target_y, target_z, earliest_arrival, latest_arrival};
        constraint.start_time = earliest_arrival;
        constraint.end_time = latest_arrival;
        constraint.is_active = true;
        constraint.description = desc;
        return constraint;
    }
};

// Основной класс для оптимизации траектории
class TrajectoryOptimizer {
private:
    SystemParams system;           // параметры системы
    OptimizationParams opt_params; // параметры оптимизации
    std::vector<AirspaceConstraint> constraints; // ограничения воздушного пространства
    
    // Приватные методы
    Function createSolver4D(const MX& x, const MX& u, const MX& t, const std::string& objective);
    void extractResult4D(const DMDict& result, TrajectoryResult& traj_result);
    
    MX createDistanceConstraint(const MX& pos_x, const MX& pos_y, const MX& pos_z, 
                              const AirspaceConstraint& constraint, const MX& time);
    
    // Add this helper method
    bool containsNaN(const std::vector<double>& vec) const;
    
public:
    // Конструкторы
    TrajectoryOptimizer();
    TrajectoryOptimizer(const SystemParams& sys, const OptimizationParams& opt);
    TrajectoryOptimizer(double mass, double dt, int N);
    
    // Методы управления ограничениями
    void addConstraint(const AirspaceConstraint& constraint);
    void clearConstraints();
    void enableConstraint(size_t index, bool enable = true);
    
    // Устанавливает параметры системы
    void setSystemParams(const SystemParams& params);
    
    // Устанавливает параметры оптимизации
    void setOptimizationParams(const OptimizationParams& params);
    
    // Метод для построения оптимальной 4D-траектории
    TrajectoryResult optimize4D(const std::string& objective = "fuel", 
                             const std::map<std::string, double>& kwargs = {});
    
    // Метод для визуализации результатов
    void visualizeTrajectory(const TrajectoryResult& result) const;
    
    // Экспорт траектории в файл для анализа
    void exportTrajectory(const TrajectoryResult& result, const std::string& filename) const;
    
    TrajectoryResult debug4D(const std::string& objective = "fuel",
                            const std::map<std::string, double>& kwargs = {}) {
        
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
        
        // Построение уравнений без отладочных сообщений
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
        
        // Целевая функция
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
        
        // Тихие настройки решателя
        Dict solver_opts;
        solver_opts["ipopt.print_level"] = 0;  // Отключаем вывод IPOPT
        solver_opts["ipopt.tol"] = 1e-3;       
        solver_opts["ipopt.max_iter"] = 500;   
        solver_opts["ipopt.acceptable_tol"] = 1e-2;
        solver_opts["print_time"] = 0;         // Отключаем вывод времени

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
        
        // Границы ограничений (все уравнения = 0)
        std::vector<double> lbg(equations.size(), 0.0);
        std::vector<double> ubg(equations.size(), 0.0);
        
        // Начальное приближение - линейная интерполяция
        std::vector<double> x0_guess(n_states + n_controls, 0.0);
        for (int k = 0; k <= N; k++) {
            double t = static_cast<double>(k) / N;
            for (int i = 0; i < 3; i++) {
                x0_guess[6*k + i] = x0[i] + t * (xf[i] - x0[i]);
            }
        }
        
        TrajectoryResult traj_result;
        
        try {
            // Call solver
            DMDict result;
            result = solver(DMDict{
                {"x0", x0_guess},
                {"lbx", lbx},
                {"ubx", ubx},
                {"lbg", lbg},
                {"ubg", ubg}
            });
            
            // Измеряем время после решения
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            
            // Get optimization result
            DM x_opt = result.at("x");
            
            // Fill result structure
            traj_result.positions.resize(N + 1);
            traj_result.velocities.resize(N + 1);
            traj_result.controls.resize(N);
            traj_result.times.resize(N + 1, 0.0);
            
            // Calculate times
            for (int k = 1; k <= N; k++) {
                traj_result.times[k] = traj_result.times[k-1] + opt_params.dt;
            }
            
            // Extract states and controls
            for (int k = 0; k <= N; k++) {
                traj_result.positions[k].resize(3);
                traj_result.velocities[k].resize(3);
                
                for (int i = 0; i < 3; i++) {
                    traj_result.positions[k][i] = x_opt(6*k + i).scalar();
                    traj_result.velocities[k][i] = x_opt(6*k + i + 3).scalar();
                }
                
                if (k < N) {
                    traj_result.controls[k].resize(3);
                    for (int i = 0; i < 3; i++) {
                        traj_result.controls[k][i] = x_opt(n_states + 3*k + i).scalar();
                    }
                }
            }
            
            traj_result.objective_value = result.at("f").scalar();
            traj_result.computation_time = duration.count() / 1000.0;
            traj_result.total_flight_time = traj_result.times.back();
            traj_result.fuel_consumption = 0.0; // Простой расчет можно добавить позже
            
        } catch (const std::exception& e) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            
            traj_result.computation_time = -1.0; // Индикатор ошибки
            std::cerr << "Error in debug4D: " << e.what() << std::endl;
        }
        
        opt_params = saved_params;
        return traj_result;
    }
    
    // Accessor methods for system and optimization parameters
    const SystemParams& getSystemParams() const { return system; }
    const OptimizationParams& getOptimizationParams() const { return opt_params; }
    const std::vector<AirspaceConstraint>& getConstraints() const { return constraints; }
};