#pragma once

#include <casadi/casadi.hpp>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <memory>

using namespace casadi;

// Структура для хранения параметров системы
struct SystemParams {
    double mass;        // масса объекта
    double max_force;   // максимальное управляющее воздействие
    
    SystemParams(double m = 1.0, double f = 5.0) : mass(m), max_force(f) {}
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
    std::vector<std::vector<double>> positions;  // координаты x,y в каждый момент времени
    std::vector<std::vector<double>> velocities; // скорости по x,y в каждый момент времени
    std::vector<std::vector<double>> controls;   // управления по x,y в каждый момент времени
    double objective_value;                      // значение целевой функции
    double computation_time;                     // время вычислений
    
    TrajectoryResult() : objective_value(0.0), computation_time(0.0) {}
};

// Основной класс для оптимизации траектории
class TrajectoryOptimizer {
private:
    SystemParams system;           // параметры системы
    OptimizationParams opt_params; // параметры оптимизации
    
    // Приватные методы
    Function createSolver(const MX& x, const MX& u, const std::string& objective);
    void extractResult(const DMDict& result, TrajectoryResult& traj_result);
    
public:
    // Конструкторы
    TrajectoryOptimizer();
    TrajectoryOptimizer(const SystemParams& sys, const OptimizationParams& opt);
    TrajectoryOptimizer(double mass, double dt, int N);
    
    // Устанавливает параметры системы
    void setSystemParams(const SystemParams& params);
    
    // Устанавливает параметры оптимизации
    void setOptimizationParams(const OptimizationParams& params);
    
    // Метод для построения оптимальной траектории
    TrajectoryResult optimize(const std::string& objective = "fuel", 
                              const std::map<std::string, double>& kwargs = {});
    
    // Метод для визуализации результатов
    void visualizeTrajectory(const TrajectoryResult& result) const;
    
    // Метод обратной совместимости с предыдущей реализацией
    DMDict trajectory(const std::string& objective = "fuel", 
                      const std::map<std::string, double>& kwargs = {});
    void visualizeTrajectory(const DMDict& result);
};