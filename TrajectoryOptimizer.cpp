#include "TrajectoryOptimizer.hpp"
#include <limits>
#include <cmath>

TrajectoryOptimizer::TrajectoryOptimizer(double mass, double dt, int N) 
    : mass(mass), dt(dt), N(N) {}

DMDict TrajectoryOptimizer::trajectory(const std::string& objective, 
                const std::map<std::string, double>& kwargs) {
    
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
        // следующая позиция = текущая позиция + скорость * dt
        MX next_pos_x = pos_x + vel_x * dt;
        MX next_pos_y = pos_y + vel_y * dt;
        
        // следующая скорость = текущая скорость + ускорение * dt
        MX next_vel_x = vel_x + (force_x / mass) * dt;
        MX next_vel_y = vel_y + (force_y / mass) * dt;
        
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
        obj = N * dt;
        // Добавляем штраф за медленное движение
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
    Function solver = nlpsol("solver", "ipopt", nlp, solver_opts);
    
    // Границы для переменных оптимизации
    int n_states = 4 * (N + 1);
    int n_controls = 2 * N;
    
    std::vector<double> lbx(n_states + n_controls, -std::numeric_limits<double>::infinity());
    std::vector<double> ubx(n_states + n_controls, std::numeric_limits<double>::infinity());
    
    // Ограничения на управление
    double max_force = 5.0; // максимальная сила
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