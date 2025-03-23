#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>
#include <casadi/casadi.hpp>
#include <windows.h>  // For SetConsoleOutputCP

using namespace casadi;
using namespace std::chrono;

// Функция для решения задачи оптимизации с заданным решателем
std::pair<DMDict, double> solve_with_solver(const std::string& solver_name, 
                                           const MXDict& nlp, 
                                           const std::vector<double>& x0,
                                           const std::vector<double>& lbx,
                                           const std::vector<double>& ubx,
                                           const std::vector<double>& lbg,
                                           const std::vector<double>& ubg) {
    Dict opts;
    
    // Базовые настройки для всех солверов
    opts["print_time"] = 0;
    opts["print_out"] = 0;
    opts["print_in"] = 0;
    
    // Специфичные настройки для отдельных солверов
    if (solver_name == "ipopt") {
        opts["ipopt.print_level"] = 0;
        opts["ipopt.tol"] = 1e-8;
        opts["ipopt.max_iter"] = 1000;
    } else if (solver_name == "bonmin") {
        opts["bonmin.tol"] = 1e-8;
        opts["bonmin.print_level"] = 0;
    } else if (solver_name == "knitro") {
        opts["knitro.tol"] = 1e-8;
        opts["knitro.outlev"] = 0;
    } else if (solver_name == "snopt") {
        opts["snopt.Major_optimality_tolerance"] = 1e-8;
    } else if (solver_name == "sqpmethod") {
        opts["qpsol"] = "qpoases";
        opts["max_iter"] = 1000;
        // Удаляем ipopt настройки для sqpmethod
    } else if (solver_name == "scpgen") {
        opts["qpsol"] = "qpoases";
        opts["max_iter"] = 1000;
        // Дополнительные настройки для scpgen
        opts["qpsol_options.printLevel"] = "none"; // Отключаем вывод от QP решателя
        opts["qpsol_options.enableRegularisation"] = true; // Включаем регуляризацию
        opts["reg_threshold"] = 1e-4; // Увеличиваем порог регуляризации
        opts["error_on_fail"] = false; // Не выдавать ошибку при неудаче QP
    }
    
    // Создаем решатель
    Function solver = nlpsol("solver", solver_name, nlp, opts);
    
    // Формируем входные аргументы
    DMDict arg = {
        {"x0", DM(x0)},
        {"lbx", DM(lbx)},
        {"ubx", DM(ubx)},
        {"lbg", DM(lbg)},
        {"ubg", DM(ubg)}
    };
    
    // Измеряем время выполнения
    auto start = high_resolution_clock::now();
    DMDict res = solver(arg);
    auto end = high_resolution_clock::now();
    
    duration<double> elapsed = end - start;
    
    return {res, elapsed.count()};
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    
    std::cout << "=== Сравнение решателей для нелинейной оптимизации ===" << std::endl;
    std::cout << "Задача: минимизация функции Розенброка с ограничениями" << std::endl << std::endl;
    
    // Создаем символическую переменную (3-мерный вектор)
    MX x = MX::sym("x", 3);
    
    // Целевая функция: расширенная "чаша Розенброка"
    MX f = 100*pow(x(1)-pow(x(0),2),2) + pow(1-x(0),2) + 
           100*pow(x(2)-pow(x(1),2),2) + pow(1-x(1),2);
    
    // Выводим формулу целевой функции
    std::cout << "Целевая функция:" << std::endl;
    std::cout << "f(x) = 100*(x₁-x₀²)² + (1-x₀)² + 100*(x₂-x₁²)² + (1-x₁)²" << std::endl << std::endl;
    
    // Определение ограничений
    std::cout << "Ограничения задачи:" << std::endl;
    std::cout << "1. x₀ + x₁ + x₂ = 2         (линейное равенство)" << std::endl;
    std::cout << "2. x₀² + x₁² ≤ 1.5         (нелинейное неравенство)" << std::endl;
    std::cout << "3. x₂ ≥ 0                  (линейное неравенство)" << std::endl << std::endl;
    
    // Ограничение равенства: x0 + x1 + x2 = 2
    MX g_eq = x(0) + x(1) + x(2) - 2;
    
    // Ограничения неравенства: x0^2 + x1^2 <= 1.5, x2 >= 0
    MX g_ineq1 = pow(x(0),2) + pow(x(1),2) - 1.5;
    MX g_ineq2 = -x(2);
    
    // Объединяем все ограничения
    MX g = MX::vertcat({g_eq, g_ineq1, g_ineq2});
    
    // Формируем задачу NLP
    MXDict nlp = {{"x", x}, {"f", f}, {"g", g}};
    
    // Начальная точка и границы
    std::vector<double> x0 = {0.5, 0.5, 1.0};
    
    // Границы переменных
    std::vector<double> lbx = {-10.0, -10.0, -10.0};
    std::vector<double> ubx = {10.0, 10.0, 10.0};
    
    // Границы ограничений
    std::vector<double> lbg = {0.0, -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    std::vector<double> ubg = {0.0, 0.0, 0.0};
    
    // Список солверов для тестирования
    std::vector<std::string> solvers = {"ipopt", "sqpmethod", "scpgen"};
    
    // Попробуем дополнительные солверы, если они доступны
    std::cout << "\nПроверка доступности дополнительных солверов..." << std::endl;
    
    // Нет необходимости в отключении предупреждений, так как методы не существуют
    // Просто используем блок try-catch для обработки ошибок
    bool prev_warnings = false;
    
    try {
        Function test_solver = nlpsol("test", "bonmin", nlp);
        solvers.push_back("bonmin");
        std::cout << "✓ Bonmin доступен" << std::endl;
    } catch (const std::exception&) {
        std::cout << "✗ Bonmin недоступен" << std::endl;
    }
    
    try {
        Function test_solver = nlpsol("test", "knitro", nlp);
        solvers.push_back("knitro");
        std::cout << "✓ Knitro доступен" << std::endl;
    } catch (const std::exception&) {
        std::cout << "✗ Knitro недоступен (требуется коммерческая лицензия)" << std::endl;
    }
    
    try {
        Function test_solver = nlpsol("test", "snopt", nlp);
        solvers.push_back("snopt");
        std::cout << "✓ SNOPT доступен" << std::endl;
    } catch (const std::exception&) {
        std::cout << "✗ SNOPT недоступен (требуется коммерческая лицензия)" << std::endl;
    }
    
    // Таблица для результатов
    std::cout << "\n=== Сравнение решателей ===\n" << std::endl;
    std::cout << std::left << std::setw(12) << "Решатель" 
              << std::setw(20) << "Оптимальное значение" 
              << std::setw(15) << "Время (сек)" 
              << "Решение" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    // Отключаем вывод от нижележащих солверов
    auto prev_stdout = std::cout.rdbuf();
    std::stringstream null_stream;
    
    // Проверяем каждый решатель
    for (const auto& solver_name : solvers) {
        try {
            // Перенаправляем стандартный вывод для подавления сообщений от солверов
            if (solver_name != "ipopt") { // Оставляем вывод только для ipopt как пример
                std::cout.rdbuf(null_stream.rdbuf());
            }
            
            // Решаем задачу и измеряем время
            auto [res, time] = solve_with_solver(solver_name, nlp, x0, lbx, ubx, lbg, ubg);
            
            // Восстанавливаем стандартный вывод
            std::cout.rdbuf(prev_stdout);
            
            // Извлекаем результаты
            DM x_opt = res.at("x");
            DM f_opt = res.at("f");
            
            // Формируем строку с решением
            std::stringstream solution_str;
            solution_str << "[" << x_opt(0).scalar() << ", " 
                         << x_opt(1).scalar() << ", " 
                         << x_opt(2).scalar() << "]";
            
            // Выводим результаты в таблицу
            std::cout << std::left << std::setw(12) << solver_name 
                      << std::setw(20) << f_opt.scalar() 
                      << std::setw(15) << std::fixed << std::setprecision(6) << time 
                      << solution_str.str() << std::endl;
        }
        catch (const std::exception& e) {
            // Восстанавливаем стандартный вывод в случае ошибки
            std::cout.rdbuf(prev_stdout);
            std::cout << std::left << std::setw(12) << solver_name 
                      << "Ошибка: " << e.what() << std::endl;
        }
    }
    
    std::cout << "\n=== Детальный анализ результатов IPOPT ===" << std::endl;
    
    // Запускаем IPOPT с подробным выводом для детального анализа
    Dict opts;
    opts["ipopt.print_level"] = 5;
    opts["ipopt.tol"] = 1e-8;
    opts["ipopt.max_iter"] = 1000;
    opts["print_time"] = 1;
    
    Function solver = nlpsol("solver", "ipopt", nlp, opts);
    
    // Формируем входные аргументы
    DMDict arg = {
        {"x0", DM(x0)},
        {"lbx", DM(lbx)},
        {"ubx", DM(ubx)},
        {"lbg", DM(lbg)},
        {"ubg", DM(ubg)}
    };
    
    std::cout << "\nРешение задачи с помощью IPOPT (подробный вывод)..." << std::endl;
    // Решаем задачу оптимизации
    DMDict res = solver(arg);
    
    // Получаем результаты
    DM x_opt = res.at("x");
    DM f_opt = res.at("f");
    DM g_opt = res.at("g");
    DM lam_g = res.at("lam_g");
    
    // Выводим результаты в удобном формате
    std::cout << "\n=== Результаты оптимизации (IPOPT) ===" << std::endl;
    std::cout << "Оптимальное решение x* = [" 
              << x_opt(0).scalar() << ", " 
              << x_opt(1).scalar() << ", " 
              << x_opt(2).scalar() << "]" << std::endl;
    std::cout << "Минимальное значение f(x*) = " << f_opt << std::endl;
    
    std::cout << "\n=== Анализ ограничений ===" << std::endl;
    std::cout << "Значения ограничений:" << std::endl;
    std::cout << "1. x₀ + x₁ + x₂ - 2 = " << g_opt(0) << std::endl;
    std::cout << "2. x₀² + x₁² - 1.5 = " << g_opt(1) << std::endl;
    std::cout << "3. -x₂ = " << g_opt(2) << std::endl;
    
    std::cout << "\nМножители Лагранжа:" << std::endl;
    std::cout << "λ₁ = " << lam_g(0) << " (для ограничения равенства)" << std::endl;
    std::cout << "λ₂ = " << lam_g(1) << " (для нелинейного ограничения)" << std::endl;
    std::cout << "λ₃ = " << lam_g(2) << " (для ограничения x₂ ≥ 0)" << std::endl;
    
    // Проверка активности ограничений
    std::cout << "\n=== Проверка ограничений ===" << std::endl;
    std::cout << "1. x₀ + x₁ + x₂ = 2: " 
              << (fabs(g_opt(0).scalar()) < 1e-6 ? "Выполнено ✓" : "Не выполнено ✗") << std::endl;
    
    // Для неравенств проверяем, активные они или нет
    bool is_active_ineq1 = fabs(g_opt(1).scalar()) < 1e-6;
    std::cout << "2. x₀² + x₁² ≤ 1.5: " 
              << (g_opt(1).scalar() <= 1e-6 ? "Выполнено ✓" : "Не выполнено ✗")
              << (is_active_ineq1 ? " (активное)" : " (неактивное)") << std::endl;
    
    bool is_active_ineq2 = fabs(g_opt(2).scalar()) < 1e-6;
    std::cout << "3. x₂ ≥ 0: " 
              << (g_opt(2).scalar() <= 1e-6 ? "Выполнено ✓" : "Не выполнено ✗") 
              << (is_active_ineq2 ? " (активное)" : " (неактивное)") << std::endl;
    
    
    return 0;
}