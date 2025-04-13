#pragma once

#include <casadi/casadi.hpp>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <memory>
#include <fstream>

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
    Function createSolver2D(const MX& x, const MX& u, const std::string& objective);
    Function createSolver4D(const MX& x, const MX& u, const MX& t, const std::string& objective);
    
    void extractResult2D(const DMDict& result, TrajectoryResult& traj_result);
    void extractResult4D(const DMDict& result, TrajectoryResult& traj_result);
    
    MX createDistanceConstraint(const MX& pos_x, const MX& pos_y, const MX& pos_z, 
                              const AirspaceConstraint& constraint, const MX& time);
    
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
    
    // Основной метод для 2D-траектории
    TrajectoryResult optimize2D(const std::string& objective = "fuel", 
                             const std::map<std::string, double>& kwargs = {});
    
    // Метод для построения оптимальной 4D-траектории
    TrajectoryResult optimize4D(const std::string& objective = "fuel", 
                             const std::map<std::string, double>& kwargs = {});
    
    // Метод для визуализации результатов
    void visualizeTrajectory(const TrajectoryResult& result) const;
    
    // Экспорт траектории в файл для анализа
    void exportTrajectory(const TrajectoryResult& result, const std::string& filename) const;
    
    // Метод обратной совместимости с предыдущей реализацией
    DMDict trajectory(const std::string& objective = "fuel", 
                      const std::map<std::string, double>& kwargs = {});
    void visualizeTrajectory(const DMDict& result);
};