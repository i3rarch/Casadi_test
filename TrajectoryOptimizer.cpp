#include "TrajectoryOptimizer.hpp"
#include <limits>
#include <cmath>
#include <chrono>
#include <fstream>

// Конструкторы
TrajectoryOptimizer::TrajectoryOptimizer()
    : system(), opt_params() {}

TrajectoryOptimizer::TrajectoryOptimizer(const SystemParams& sys, const OptimizationParams& opt)
    : system(sys), opt_params(opt) {}

TrajectoryOptimizer::TrajectoryOptimizer(double mass, double dt, int N) 
    : system(mass), opt_params(N, dt) {}

// Устанавливает параметры системы
void TrajectoryOptimizer::setSystemParams(const SystemParams& params) {
    system = params;
}

// Устанавливает параметры оптимизации
void TrajectoryOptimizer::setOptimizationParams(const OptimizationParams& params) {
    opt_params = params;
}

// Создание решателя для оптимизации 2D траектории
Function TrajectoryOptimizer::createSolver2D(const MX& x, const MX& u, const std::string& objective) {
    int N = opt_params.N;
    double dt = opt_params.dt;
    double mass = system.mass;
    
    // Уравнения динамики
    std::vector<MX> equations;
    for (int k = 0; k < N; k++) {
        // Распаковка состояния
        MX pos_x = x(0, k); 
        MX pos_y = x(1, k);
        MX vel_x = x(2, k);
        MX vel_y = x(3, k);
        
        // Управляющие воздействия
        MX force_x = u(0, k);
        MX force_y = u(1, k);
        
        // Уравнения движения (дискретизация):
        MX next_pos_x = pos_x + vel_x * dt;
        MX next_pos_y = pos_y + vel_y * dt;
        
        MX next_vel_x = vel_x + (force_x / mass) * dt;
        MX next_vel_y = vel_y + (force_y / mass) * dt;
        
        // Накладываем ограничения
        equations.push_back(x(0, k+1) - next_pos_x);
        equations.push_back(x(1, k+1) - next_pos_y);
        equations.push_back(x(2, k+1) - next_vel_x);
        equations.push_back(x(3, k+1) - next_vel_y);
    }
    
    // Целевая функция
    MX obj;
    if (objective == "fuel") {
        // Минимизация расхода топлива (сумма абсолютных значений управлений)
        obj = 0;
        for (int k = 0; k < N; k++) {
            obj += fabs(u(0, k)) + fabs(u(1, k));
        }
    } else if (objective == "time") {
        // Минимизация времени (аппроксимация)
        obj = N * dt;
        // Штраф за медленное движение
        for (int k = 0; k < N; k++) {
            obj += dt * 0.1 * (1.0 / (fabs(u(0, k)) + fabs(u(1, k)) + 0.1));
        }
    } else { // "energy" по умолчанию
        // Минимизация энергии (сумма квадратов управлений)
        obj = 0;
        for (int k = 0; k < N; k++) {
            obj += pow(u(0, k), 2) + pow(u(1, k), 2);
        }
    }
    
    // Параметризуем переменные оптимизации
    std::vector<MX> vars = {reshape(x, -1, 1), reshape(u, -1, 1)};
    MX opt_vars = vertcat(vars);
    
    // Формируем ограничения
    MX g = vertcat(equations);
    
    // Создаем задачу оптимизации
    MXDict nlp = {{"x", opt_vars}, {"f", obj}, {"g", g}};
    
    // Настройки решателя
    Dict solver_opts;
    solver_opts["ipopt.print_level"] = 0;
    solver_opts["ipopt.tol"] = 1e-6;
    solver_opts["ipopt.max_iter"] = 1000;
    solver_opts["print_time"] = 0;
    
    // Создаем решатель
    return nlpsol("solver", opt_params.solver, nlp, solver_opts);
}

// Создание решателя для 4D-траектории
Function TrajectoryOptimizer::createSolver4D(const MX& x, const MX& u, const MX& t, const std::string& objective) {
    int N = opt_params.N;
    double mass = system.mass;
    
    // Уравнения динамики
    std::vector<MX> equations;
    
    // Временные шаги могут быть переменными для оптимизации времени полета
    for (int k = 0; k < N; k++) {
        // Распаковка состояния: позиция (x,y,z) и скорость (vx,vy,vz)
        MX pos_x = x(0, k);
        MX pos_y = x(1, k);
        MX pos_z = x(2, k);
        MX vel_x = x(3, k);
        MX vel_y = x(4, k);
        MX vel_z = x(5, k);
        
        // Управляющие воздействия
        MX force_x = u(0, k);
        MX force_y = u(1, k);
        MX force_z = u(2, k);
        
        // Временной шаг между точками k и k+1
        MX dt = t(k);
        
        // Уравнения движения самолета (упрощенная модель)
        MX next_pos_x = pos_x + vel_x * dt;
        MX next_pos_y = pos_y + vel_y * dt;
        MX next_pos_z = pos_z + vel_z * dt;
        
        MX next_vel_x = vel_x + (force_x / mass) * dt;
        MX next_vel_y = vel_y + (force_y / mass) * dt;
        MX next_vel_z = vel_z + (force_z / mass) * dt - 9.81 * dt; // С учетом гравитации
        
        // Ограничения на следующие значения
        equations.push_back(x(0, k+1) - next_pos_x);
        equations.push_back(x(1, k+1) - next_pos_y);
        equations.push_back(x(2, k+1) - next_pos_z);
        equations.push_back(x(3, k+1) - next_vel_x);
        equations.push_back(x(4, k+1) - next_vel_y);
        equations.push_back(x(5, k+1) - next_vel_z);
    }
    
    // Дополнительные ограничения для воздушного пространства
    std::vector<MX> airspace_constraints;
    
    // Время начала траектории
    MX current_time = 0;
    
    // Проверяем все ограничения воздушного пространства
    for (int k = 0; k < N; k++) {
        // Текущее время в точке k
        if (k > 0) {
            current_time = current_time + t(k-1);
        }
        
        MX pos_x = x(0, k);
        MX pos_y = x(1, k);
        MX pos_z = x(2, k);
        
        // Ограничения на высоту полета
        airspace_constraints.push_back(pos_z - system.min_altitude); // Высота должна быть >= min_altitude
        airspace_constraints.push_back(system.max_altitude - pos_z); // Высота должна быть <= max_altitude
        
        // Ограничения скорости
        MX speed = sqrt(pow(x(3, k), 2) + pow(x(4, k), 2) + pow(x(5, k), 2));
        airspace_constraints.push_back(system.max_speed - speed); // Скорость должна быть <= max_speed
        
        // Применяем все активные ограничения воздушного пространства
        for (const auto& constraint : constraints) {
            if (!constraint.is_active) continue;
            
            // Проверяем, активно ли ограничение в текущий момент времени
            if (double(current_time) >= constraint.start_time && double(current_time) <= constraint.end_time) {
                
                // Добавляем ограничение в зависимости от его типа
                airspace_constraints.push_back(createDistanceConstraint(pos_x, pos_y, pos_z, constraint, current_time));
            }
        }
    }
    
    // Целевая функция
    MX obj;
    if (objective == "fuel") {
        // Минимизация расхода топлива (сумма абсолютных значений управлений)
        obj = 0;
        for (int k = 0; k < N; k++) {
            obj += fabs(u(0, k)) + fabs(u(1, k)) + fabs(u(2, k));
        }
    } else if (objective == "time") {
        // Минимизация времени полета - сумма всех временных шагов
        obj = 0;
        for (int k = 0; k < N; k++) {
            obj += t(k);
        }
    } else if (objective == "balanced") {
        // Сбалансированная оптимизация (время и расход топлива)
        double time_weight = 0.3;
        double fuel_weight = 0.7;
        
        // Компонента времени
        MX time_obj = 0;
        for (int k = 0; k < N; k++) {
            time_obj += t(k);
        }
        
        // Компонента расхода топлива
        MX fuel_obj = 0;
        for (int k = 0; k < N; k++) {
            fuel_obj += fabs(u(0, k)) + fabs(u(1, k)) + fabs(u(2, k));
        }
        
        // Общая целевая функция
        obj = time_weight * time_obj / 3600.0 + fuel_weight * fuel_obj / 1000000.0;
    } else { // "energy" по умолчанию
        // Минимизация энергии (сумма квадратов управлений)
        obj = 0;
        for (int k = 0; k < N; k++) {
            obj += pow(u(0, k), 2) + pow(u(1, k), 2) + pow(u(2, k), 2);
        }
    }
    
    // Параметризуем переменные оптимизации
    std::vector<MX> vars = {reshape(x, -1, 1), reshape(u, -1, 1), t};
    MX opt_vars = vertcat(vars);
    
    // Формируем все ограничения
    MX g = vertcat(equations);
    
    // Добавляем ограничения воздушного пространства, если они есть
    if (!airspace_constraints.empty()) {
        MX airspace_g = vertcat(airspace_constraints);
        g = vertcat(std::vector<MX>{g, airspace_g});
    }
    
    // Создаем задачу оптимизации
    MXDict nlp = {{"x", opt_vars}, {"f", obj}, {"g", g}};
    
    // Настройки решателя
    Dict solver_opts;
    solver_opts["ipopt.print_level"] = 0;
    solver_opts["ipopt.tol"] = 1e-4;
    solver_opts["ipopt.max_iter"] = 3000;
    solver_opts["ipopt.acceptable_tol"] = 1e-2;
    solver_opts["print_time"] = 0;
    
    // Создаем решатель
    return nlpsol("solver", opt_params.solver, nlp, solver_opts);
}

// Создание ограничений расстояния для различных типов ограничений воздушного пространства
MX TrajectoryOptimizer::createDistanceConstraint(const MX& pos_x, const MX& pos_y, const MX& pos_z, 
                                             const AirspaceConstraint& constraint, const MX& time) {
    switch (constraint.constraint_type) {
        case AirspaceConstraintType::SPHERE: {
            // Параметры: [center_x, center_y, center_z, radius]
            double cx = constraint.parameters[0];
            double cy = constraint.parameters[1];
            double cz = constraint.parameters[2];
            double r = constraint.parameters[3];
            
            // Расстояние от точки до центра сферы должно быть больше радиуса (избегание зоны)
            return pow(pos_x - cx, 2) + pow(pos_y - cy, 2) + pow(pos_z - cz, 2) - pow(r, 2);
        }
        
        case AirspaceConstraintType::CYLINDER: {
            // Параметры: [center_x, center_y, min_z, max_z, radius]
            double cx = constraint.parameters[0];
            double cy = constraint.parameters[1];
            double min_z = constraint.parameters[2];
            double max_z = constraint.parameters[3];
            double r = constraint.parameters[4];
            
            // Находимся ли по высоте в пределах цилиндра
            MX in_height_range = MX::if_else(
                pos_z >= min_z && pos_z <= max_z,
                1.0, -1.0
            );
            
            // Горизонтальное расстояние до оси цилиндра
            MX dist_xy = pow(pos_x - cx, 2) + pow(pos_y - cy, 2);
            
            // Если находимся в диапазоне высот цилиндра, должны быть за пределами радиуса
            return in_height_range * (dist_xy - pow(r, 2));
        }
        
        case AirspaceConstraintType::CORRIDOR: {
            // Параметры: [start_x, start_y, start_z, end_x, end_y, end_z, width]
            double sx = constraint.parameters[0];
            double sy = constraint.parameters[1];
            double sz = constraint.parameters[2];
            double ex = constraint.parameters[3];
            double ey = constraint.parameters[4];
            double ez = constraint.parameters[5];
            double width = constraint.parameters[6];
            
            // Вектор направления коридора
            double dx = ex - sx;
            double dy = ey - sy;
            double dz = ez - sz;
            double len = sqrt(dx*dx + dy*dy + dz*dz);
            dx /= len; dy /= len; dz /= len;
            
            // Проекция вектора (pos - start) на направление коридора
            MX proj = (pos_x - sx) * dx + (pos_y - sy) * dy + (pos_z - sz) * dz;
            
            // Ближайшая точка на коридоре
            MX closest_x = sx + proj * dx;
            MX closest_y = sy + proj * dy;
            MX closest_z = sz + proj * dz;
            
            // Расстояние до оси коридора
            MX dist = sqrt(pow(pos_x - closest_x, 2) + pow(pos_y - closest_y, 2) + pow(pos_z - closest_z, 2));
            
            // Должны находиться внутри коридора (расстояние <= width) и между началом и концом коридора
            MX in_corridor_bounds = MX::if_else(
                proj >= 0 && proj <= len,
                width - dist, // Если внутри диапазона длины коридора - должны быть ближе чем width
                -1.0          // Иначе считаем ограничение выполненным
            );
            
            return in_corridor_bounds;
        }
        
        case AirspaceConstraintType::TIME_WINDOW: {
            // Параметры: [target_x, target_y, target_z, earliest_arrival, latest_arrival]
            double tx = constraint.parameters[0];
            double ty = constraint.parameters[1];
            double tz = constraint.parameters[2];
            double t_early = constraint.parameters[3];
            double t_late = constraint.parameters[4];
            
            // Расстояние до целевой точки
            MX dist_to_target = sqrt(pow(pos_x - tx, 2) + pow(pos_y - ty, 2) + pow(pos_z - tz, 2));
            
            // Если мы близко к целевой точке (меньше 100м), проверяем временное окно
            MX at_target = MX::if_else(
                dist_to_target < 100.0,
                (time - t_early) * (t_late - time), // Должно быть положительно, если внутри окна
                1.0 // Иначе считаем ограничение выполненным
            );
            
            return at_target;
        }
        
        default:
            // По умолчанию возвращаем константу - ограничение всегда выполнено
            return MX(1.0);
    }
}

// Добавление ограничения воздушного пространства
void TrajectoryOptimizer::addConstraint(const AirspaceConstraint& constraint) {
    constraints.push_back(constraint);
}

// Очистка всех ограничений
void TrajectoryOptimizer::clearConstraints() {
    constraints.clear();
}

// Включение/выключение ограничения по индексу
void TrajectoryOptimizer::enableConstraint(size_t index, bool enable) {
    if (index < constraints.size()) {
        constraints[index].is_active = enable;
    }
}

// Извлечение результатов из решения оптимизации
void TrajectoryOptimizer::extractResult2D(const DMDict& result, TrajectoryResult& traj_result) {
    DM x_opt = result.at("x");
    int N = opt_params.N;
    
    // Получаем значение целевой функции
    traj_result.objective_value = result.at("f").scalar();
    
    // Резервируем память для данных
    traj_result.positions.resize(N + 1);
    traj_result.velocities.resize(N + 1);
    traj_result.controls.resize(N);
    
    // Извлекаем значения состояний и управления
    for (int k = 0; k <= N; k++) {
        double pos_x = x_opt(4*k).scalar();
        double pos_y = x_opt(4*k+1).scalar();
        double vel_x = x_opt(4*k+2).scalar();
        double vel_y = x_opt(4*k+3).scalar();
        
        traj_result.positions[k] = {pos_x, pos_y};
        traj_result.velocities[k] = {vel_x, vel_y};
        
        if (k < N) {
            double force_x = x_opt(4*(N+1) + 2*k).scalar();
            double force_y = x_opt(4*(N+1) + 2*k+1).scalar();
            traj_result.controls[k] = {force_x, force_y};
        }
    }
}

// Извлечение результатов из решения оптимизации для 4D-траектории
void TrajectoryOptimizer::extractResult4D(const DMDict& result, TrajectoryResult& traj_result) {
    DM x_opt = result.at("x");
    int N = opt_params.N;
    
    // Получаем значение целевой функции
    traj_result.objective_value = result.at("f").scalar();
    
    // Резервируем память для данных
    traj_result.positions.resize(N + 1);
    traj_result.velocities.resize(N + 1);
    traj_result.controls.resize(N);
    
    // Извлекаем значения состояний и управления
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

// Метод для оптимизации траектории
TrajectoryResult TrajectoryOptimizer::optimize2D(
    const std::string& objective, 
    const std::map<std::string, double>& kwargs) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    int N = opt_params.N;
    
    // Создаем переменные состояния и управления
    MX x = MX::sym("x", 4, N+1);  // состояние: [позиция_x, позиция_y, скорость_x, скорость_y]
    MX u = MX::sym("u", 2, N);    // управление: [сила_x, сила_y]
    
    // Начальные и конечные условия
    std::vector<double> x0 = {0.0, 0.0, 0.0, 0.0};  // начальное положение и скорость
    std::vector<double> xf = {10.0, 5.0, 0.0, 0.0}; // целевое положение и скорость
    
    // Переопределение значений из kwargs
    auto it = kwargs.find("x0");
    if (it != kwargs.end()) x0[0] = it->second;
    
    it = kwargs.find("y0");
    if (it != kwargs.end()) x0[1] = it->second;
    
    it = kwargs.find("xf");
    if (it != kwargs.end()) xf[0] = it->second;
    
    it = kwargs.find("yf"); 
    if (it != kwargs.end()) xf[1] = it->second;
    
    // Создаем решатель
    Function solver = createSolver2D(x, u, objective);
    
    // Границы для переменных оптимизации
    int n_states = 4 * (N + 1);
    int n_controls = 2 * N;
    
    std::vector<double> lbx(n_states + n_controls, -std::numeric_limits<double>::infinity());
    std::vector<double> ubx(n_states + n_controls, std::numeric_limits<double>::infinity());
    
    // Ограничения на управление
    double max_force = system.max_force_xy;  // Changed from max_force to max_force_xy
    for (int i = n_states; i < n_states + n_controls; i++) {
        lbx[i] = -max_force;
        ubx[i] = max_force;
    }
    
    // Начальные условия
    lbx[0] = ubx[0] = x0[0]; // начальная позиция x
    lbx[1] = ubx[1] = x0[1]; // начальная позиция y
    lbx[2] = ubx[2] = x0[2]; // начальная скорость x
    lbx[3] = ubx[3] = x0[3]; // начальная скорость y
    
    // Конечные условия
    lbx[n_states-4] = ubx[n_states-4] = xf[0]; // конечная позиция x
    lbx[n_states-3] = ubx[n_states-3] = xf[1]; // конечная позиция y
    lbx[n_states-2] = ubx[n_states-2] = xf[2]; // конечная скорость x
    lbx[n_states-1] = ubx[n_states-1] = xf[3]; // конечная скорость y
    
    // Границы для ограничений (все уравнения должны быть равны нулю)
    int n_equations = 4 * N;  // 4 уравнения для каждого шага времени
    std::vector<double> lbg(n_equations, 0.0);
    std::vector<double> ubg(n_equations, 0.0);
    
    // Начальное предположение для переменных оптимизации
    std::vector<double> x0_guess(n_states + n_controls, 0.0);
    
    // Вызываем решатель
    DMDict result = solver(DMDict{
        {"x0", x0_guess},
        {"lbx", lbx},
        {"ubx", ubx},
        {"lbg", lbg},
        {"ubg", ubg}
    });
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    // Заполняем структуру результата
    TrajectoryResult traj_result;
    extractResult2D(result, traj_result);
    traj_result.computation_time = duration.count() / 1000.0; // в секундах
    
    return traj_result;
}

// Метод оптимизации 4D-траектории
TrajectoryResult TrajectoryOptimizer::optimize4D(
    const std::string& objective, 
    const std::map<std::string, double>& kwargs) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    int N = opt_params.N;
    
    // Создаем переменные состояния, управления и времени
    MX x = MX::sym("x", 6, N+1);  // состояние: [pos_x, pos_y, pos_z, vel_x, vel_y, vel_z]
    MX u = MX::sym("u", 3, N);    // управление: [force_x, force_y, force_z]
    MX t = MX::sym("t", N);       // временные шаги между точками
    
    // Начальные и конечные условия
    std::vector<double> x0 = {0.0, 0.0, 3000.0, 0.0, 0.0, 0.0};  // начальное положение и скорость
    std::vector<double> xf = {100000.0, 80000.0, 10000.0, 0.0, 0.0, 0.0}; // целевое положение и скорость
    
    // Переопределение значений из kwargs
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
    
    // Создаем решатель
    Function solver = createSolver4D(x, u, t, objective);
    
    // Границы для переменных оптимизации
    int n_states = 6 * (N + 1);     // 6 состояний на каждом шаге
    int n_controls = 3 * N;         // 3 управления на каждом шаге
    int n_times = N;                // N временных шагов
    
    std::vector<double> lbx(n_states + n_controls + n_times, -std::numeric_limits<double>::infinity());
    std::vector<double> ubx(n_states + n_controls + n_times, std::numeric_limits<double>::infinity());
    
    // Ограничения на управление
    double max_force_xy = system.max_force_xy;
    double max_force_z = system.max_force_z;
    
    // Ограничения на управления x, y (горизонтальная тяга)
    for (int i = n_states; i < n_states + n_controls; i += 3) {
        lbx[i] = -max_force_xy;
        ubx[i] = max_force_xy;
        
        lbx[i+1] = -max_force_xy;
        ubx[i+1] = max_force_xy;
    }
    
    // Ограничения на управления z (вертикальная тяга)
    for (int i = n_states + 2; i < n_states + n_controls; i += 3) {
        lbx[i] = -max_force_z;
        ubx[i] = max_force_z;
    }
    
    // Ограничения на временные шаги (минимальное и максимальное время между точками)
    double min_dt = 10.0;  // минимум 10 секунд между точками
    double max_dt = 300.0; // максимум 5 минут между точками
    
    for (int i = n_states + n_controls; i < n_states + n_controls + n_times; i++) {
        lbx[i] = min_dt;
        ubx[i] = max_dt;
    }
    
    // Начальные условия
    for (int i = 0; i < 6; i++) {
        lbx[i] = ubx[i] = x0[i];
    }
    
    // Конечные условия
    for (int i = 0; i < 6; i++) {
        lbx[n_states - 6 + i] = ubx[n_states - 6 + i] = xf[i];
    }
    
    // Оцениваем количество ограничений
    int n_dynamics_constraints = 6 * N;  // 6 уравнений динамики для каждого шага
    
    // Считаем количество ограничений воздушного пространства (высота, скорость + активные зоны)
    int n_airspace_constraints = 0;
    for (int k = 0; k < N; k++) {
        n_airspace_constraints += 3; // Мин. высота, макс. высота, макс. скорость
        
        for (const auto& constraint : constraints) {
            if (constraint.is_active) {
                n_airspace_constraints++;
            }
        }
    }
    
    int total_constraints = n_dynamics_constraints + n_airspace_constraints;
    
    // Границы для ограничений
    std::vector<double> lbg(n_dynamics_constraints, 0.0); // Уравнения динамики должны быть = 0
    std::vector<double> ubg(n_dynamics_constraints, 0.0);
    
    // Границы для ограничений воздушного пространства
    for (int i = 0; i < n_airspace_constraints; i++) {
        lbg.push_back(0.0);       // Ограничения должны быть неотрицательны
        ubg.push_back(INFINITY);  // Верхней границы нет
    }
    
    // Начальное предположение для переменных оптимизации
    std::vector<double> x0_guess(n_states + n_controls + n_times, 0.0);
    
    // Начальное предположение для позиций - линейная интерполяция между начальной и конечной точками
    for (int k = 0; k <= N; k++) {
        double t = static_cast<double>(k) / N;
        for (int i = 0; i < 3; i++) {
            x0_guess[6*k + i] = x0[i] + t * (xf[i] - x0[i]);
        }
    }
    
    // Начальное предположение для временных шагов - равные интервалы
    double estimated_flight_time = 3600.0; // 1 час в секундах
    double dt_guess = estimated_flight_time / N;
    
    for (int k = 0; k < N; k++) {
        x0_guess[n_states + n_controls + k] = dt_guess;
    }
    
    // Вызываем решатель
    DMDict result;
    try {
        result = solver(DMDict{
            {"x0", x0_guess},
            {"lbx", lbx},
            {"ubx", ubx},
            {"lbg", lbg},
            {"ubg", ubg}
        });
    } catch (std::exception& e) {
        std::cerr << "Optimization failed: " << e.what() << std::endl;
        TrajectoryResult failed_result;
        failed_result.computation_time = -1.0;
        return failed_result;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    // Заполняем структуру результата
    TrajectoryResult traj_result;
    extractResult4D(result, traj_result);
    
    // Заполняем время в каждой точке
    traj_result.times.resize(N + 1, 0.0);
    DM t_opt = result.at("x").nz(Slice(n_states + n_controls, n_states + n_controls + n_times));
    
    for (int k = 1; k <= N; k++) {
        traj_result.times[k] = traj_result.times[k-1] + t_opt(k-1).scalar();
    }
    
    traj_result.computation_time = duration.count() / 1000.0; // в секундах
    traj_result.total_flight_time = traj_result.times.back(); // общее время полета
    
    // Рассчитываем расход топлива (упрощенно пропорционально силам)
    traj_result.fuel_consumption = 0.0;
    for (size_t k = 0; k < traj_result.controls.size(); k++) {
        double force_magnitude = sqrt(pow(traj_result.controls[k][0], 2) + 
                                      pow(traj_result.controls[k][1], 2) + 
                                      pow(traj_result.controls[k][2], 2));
        
        double dt = k < N ? t_opt(k).scalar() : 0.0;
        traj_result.fuel_consumption += force_magnitude * dt * 0.001; // коэффициент расхода топлива
    }
    
    return traj_result;
}

// // Метод для визуализации результатов
// void TrajectoryOptimizer::visualizeTrajectory(const TrajectoryResult& result) const {
//     std::cout << "Оптимальная траектория:" << std::endl;
//     std::cout << "------------------" << std::endl;
//     std::cout << "Время вычислений: " << result.computation_time << " сек." << std::endl;
//     std::cout << "Значение целевой функции: " << result.objective_value << std::endl;
//     std::cout << "------------------" << std::endl;
    
//     for (size_t k = 0; k < result.positions.size(); ++k) {
//         std::cout << "Шаг " << k << ": Позиция = (" 
//                   << result.positions[k][0] << ", " << result.positions[k][1]
//                   << "), Скорость = (" 
//                   << result.velocities[k][0] << ", " << result.velocities[k][1] << ")";
        
//         if (k < result.controls.size()) {
//             std::cout << std::endl << "           Управление = (" 
//                       << result.controls[k][0] << ", " << result.controls[k][1] << ")";
//         }
//         std::cout << std::endl;
//     }
// }

// Метод для визуализации результатов 4D-траектории
void TrajectoryOptimizer::visualizeTrajectory(const TrajectoryResult& result) const {
    if (result.positions.empty()) {
        std::cout << "Траектория пуста или не была рассчитана" << std::endl;
        return;
    }
    
    std::cout << "=== Результаты оптимизации 4D-траектории ===" << std::endl;
    std::cout << "Время вычислений: " << result.computation_time << " сек." << std::endl;
    std::cout << "Общее время полета: " << result.total_flight_time / 60.0 << " мин." << std::endl;
    std::cout << "Расход топлива: " << result.fuel_consumption << " единиц" << std::endl;
    std::cout << "Значение целевой функции: " << result.objective_value << std::endl;
    std::cout << "------------------" << std::endl;
    
    // Выводим только ключевые точки траектории для компактности
    int num_points = result.positions.size();
    int step = num_points <= 10 ? 1 : num_points / 10;
    
    std::cout << "Ключевые точки траектории (всего " << num_points << " точек):" << std::endl;
    for (size_t k = 0; k < result.positions.size(); k += step) {
        std::cout << "Точка " << k << ": Время = " << result.times[k] / 60.0 << " мин.";
        std::cout << ", Позиция = (" 
                  << result.positions[k][0] / 1000.0 << " км, " // Конвертируем в км для удобства
                  << result.positions[k][1] / 1000.0 << " км, " 
                  << result.positions[k][2] << " м)";
        std::cout << ", Скорость = (" 
                  << result.velocities[k][0] << ", " 
                  << result.velocities[k][1] << ", " 
                  << result.velocities[k][2] << ") м/с";
                  
        if (k < result.controls.size()) {
            std::cout << std::endl << "    Управление = (" 
                      << result.controls[k][0] / 1000.0 << ", " // Конвертируем в кН для удобства
                      << result.controls[k][1] / 1000.0 << ", " 
                      << result.controls[k][2] / 1000.0 << ") кН";
        }
        std::cout << std::endl;
    }
    
    // Отображаем последнюю точку всегда
    if (num_points > 1 && (num_points - 1) % step != 0) {
        size_t k = num_points - 1;
        std::cout << "Точка " << k << " (конечная): Время = " << result.times[k] / 60.0 << " мин.";
        std::cout << ", Позиция = (" 
                  << result.positions[k][0] / 1000.0 << " км, "
                  << result.positions[k][1] / 1000.0 << " км, " 
                  << result.positions[k][2] << " м)";
        std::cout << ", Скорость = (" 
                  << result.velocities[k][0] << ", " 
                  << result.velocities[k][1] << ", " 
                  << result.velocities[k][2] << ") м/с" << std::endl;
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

// Метод обратной совместимости
DMDict TrajectoryOptimizer::trajectory(
    const std::string& objective, 
    const std::map<std::string, double>& kwargs) {
    
    int N = opt_params.N;
    
    // Создаем переменные состояния и управления
    MX x = MX::sym("x", 4, N+1);  // состояние: [позиция_x, позиция_y, скорость_x, скорость_y]
    MX u = MX::sym("u", 2, N);    // управление: [сила_x, сила_y]
    
    // Начальные и конечные условия
    std::vector<double> x0 = {0.0, 0.0, 0.0, 0.0};  // начальное положение и скорость
    std::vector<double> xf = {10.0, 5.0, 0.0, 0.0}; // целевое положение и скорость
    
    // Переопределение значений из kwargs
    auto it = kwargs.find("x0");
    if (it != kwargs.end()) x0[0] = it->second;
    
    it = kwargs.find("y0");
    if (it != kwargs.end()) x0[1] = it->second;
    
    it = kwargs.find("xf");
    if (it != kwargs.end()) xf[0] = it->second;
    
    it = kwargs.find("yf"); 
    if (it != kwargs.end()) xf[1] = it->second;
    
    // Уравнения динамики
    std::vector<MX> equations;
    for (int k = 0; k < N; k++) {
        // Распаковка состояния
        MX pos_x = x(0, k);
        MX pos_y = x(1, k);
        MX vel_x = x(2, k);
        MX vel_y = x(3, k);
        
        // Управляющие воздействия
        MX force_x = u(0, k);
        MX force_y = u(1, k);
        
        // Уравнения движения (дискретизация):
        MX next_pos_x = pos_x + vel_x * opt_params.dt;
        MX next_pos_y = pos_y + vel_y * opt_params.dt;
        
        MX next_vel_x = vel_x + (force_x / system.mass) * opt_params.dt;
        MX next_vel_y = vel_y + (force_y / system.mass) * opt_params.dt;
        
        // Накладываем ограничения на следующий шаг
        equations.push_back(x(0, k+1) - next_pos_x);
        equations.push_back(x(1, k+1) - next_pos_y);
        equations.push_back(x(2, k+1) - next_vel_x);
        equations.push_back(x(3, k+1) - next_vel_y);
    }
    
    // Целевая функция в зависимости от выбранного критерия
    MX obj;
    if (objective == "fuel") {
        // Минимизация расхода топлива (сумма абсолютных значений управлений)
        obj = 0;
        for (int k = 0; k < N; k++) {
            obj += fabs(u(0, k)) + fabs(u(1, k));
        }
    } else if (objective == "time") {
        // Минимизация времени (аппроксимация)
        obj = N * opt_params.dt;
        // Добавляем штраф за медленное движение
        for (int k = 0; k < N; k++) {
            obj += opt_params.dt * 0.1 * (1.0 / (fabs(u(0, k)) + fabs(u(1, k)) + 0.1));
        }
    } else { // "energy" по умолчанию
        // Минимизация энергии (сумма квадратов управлений)
        obj = 0;
        for (int k = 0; k < N; k++) {
            obj += pow(u(0, k), 2) + pow(u(1, k), 2);
        }
    }
    
    // Параметризуем переменные оптимизации
    std::vector<MX> vars = {reshape(x, -1, 1), reshape(u, -1, 1)};
    MX opt_vars = vertcat(vars);
    
    // Формируем ограничения
    MX g = vertcat(equations);
    
    // Создаем задачу оптимизации
    MXDict nlp = {{"x", opt_vars}, {"f", obj}, {"g", g}};
    
    // Настройки решателя
    Dict solver_opts;
    solver_opts["ipopt.print_level"] = 0;
    solver_opts["ipopt.tol"] = 1e-6;
    solver_opts["ipopt.max_iter"] = 1000;
    solver_opts["print_time"] = 0;
    
    // Создаем решатель
    Function solver = nlpsol("solver", opt_params.solver, nlp, solver_opts);
    
    // Границы для переменных оптимизации
    int n_states = 4 * (N + 1);
    int n_controls = 2 * N;
    
    std::vector<double> lbx(n_states + n_controls, -std::numeric_limits<double>::infinity());
    std::vector<double> ubx(n_states + n_controls, std::numeric_limits<double>::infinity());
    
    // Ограничения на управление
    double max_force = system.max_force_xy;
    for (int i = n_states; i < n_states + n_controls; i++) {
        lbx[i] = -max_force;
        ubx[i] = max_force;
    }
    
    // Начальные условия
    lbx[0] = ubx[0] = x0[0]; // начальная позиция x
    lbx[1] = ubx[1] = x0[1]; // начальная позиция y
    lbx[2] = ubx[2] = x0[2]; // начальная скорость x
    lbx[3] = ubx[3] = x0[3]; // начальная скорость y
    
    // Конечные условия
    lbx[n_states-4] = ubx[n_states-4] = xf[0]; // конечная позиция x
    lbx[n_states-3] = ubx[n_states-3] = xf[1]; // конечная позиция y
    lbx[n_states-2] = ubx[n_states-2] = xf[2]; // конечная скорость x
    lbx[n_states-1] = ubx[n_states-1] = xf[3]; // конечная скорость y
    
    // Границы для ограничений (все уравнения должны быть равны нулю)
    std::vector<double> lbg(equations.size(), 0.0);
    std::vector<double> ubg(equations.size(), 0.0);
    
    // Начальное предположение для переменных оптимизации
    std::vector<double> x0_guess(n_states + n_controls, 0.0);
    
    // Вызываем решатель
    DMDict result = solver(DMDict{
        {"x0", x0_guess},
        {"lbx", lbx},
        {"ubx", ubx},
        {"lbg", lbg},
        {"ubg", ubg}
    });
    
    return result;
}

void TrajectoryOptimizer::visualizeTrajectory(const DMDict& result) {
    DM x_opt = result.at("x");
    int N = opt_params.N;
    
    std::cout << "Оптимальная траектория:" << std::endl;
    std::cout << "------------------" << std::endl;
    
    for (int k = 0; k <= N; k++) {
        double pos_x = x_opt(4*k).scalar();
        double pos_y = x_opt(4*k+1).scalar();
        double vel_x = x_opt(4*k+2).scalar();
        double vel_y = x_opt(4*k+3).scalar();
        
        std::cout << "Шаг " << k << ": Позиция = (" << pos_x << ", " << pos_y 
                  << "), Скорость = (" << vel_x << ", " << vel_y << ")" << std::endl;
        
        if (k < N) {
            double force_x = x_opt(4*(N+1) + 2*k).scalar();
            double force_y = x_opt(4*(N+1) + 2*k+1).scalar();
            std::cout << "           Управление = (" << force_x << ", " << force_y << ")" << std::endl;
        }
    }
}