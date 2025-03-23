#include <iostream>
#include <casadi/casadi.hpp>
#include <windows.h>  // For SetConsoleOutputCP

using namespace casadi;

int main() {
    SetConsoleOutputCP(CP_UTF8);
    
    std::cout << "=== Оптимизация с нелинейными ограничениями ===" << std::endl;
    std::cout << "Решаем задачу минимизации функции Розенброка с ограничениями" << std::endl << std::endl;
    
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
    
    // Настраиваем параметры решателя
    Dict opts;
    opts["ipopt.print_level"] = 5;      // Подробный вывод решателя
    opts["ipopt.tol"] = 1e-8;           // Точность решения
    opts["ipopt.max_iter"] = 1000;      // Максимальное число итераций
    opts["print_time"] = 1;
    
    // Формируем задачу NLP
    MXDict nlp = {{"x", x}, {"f", f}, {"g", g}};
    
    std::cout << "Создание решателя оптимизации (IPOPT)..." << std::endl;
    Function solver = nlpsol("solver", "ipopt", nlp, opts);
    
    // Начальная точка и границы
    std::vector<double> x0 = {0.5, 0.5, 1.0};
    std::cout << "Начальная точка: [" << x0[0] << ", " << x0[1] << ", " << x0[2] << "]" << std::endl;
    
    // Границы переменных
    std::vector<double> lbx = {-10.0, -10.0, -10.0};
    std::vector<double> ubx = {10.0, 10.0, 10.0};
    
    // Границы ограничений
    std::vector<double> lbg = {0.0, -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    std::vector<double> ubg = {0.0, 0.0, 0.0};
    
    // Формируем входные аргументы
    DMDict arg = {
        {"x0", DM(x0)},
        {"lbx", DM(lbx)},
        {"ubx", DM(ubx)},
        {"lbg", DM(lbg)},
        {"ubg", DM(ubg)}
    };
    
    std::cout << "\nЗапуск оптимизации..." << std::endl;
    // Решаем задачу оптимизации
    DMDict res = solver(arg);
    
    // Получаем результаты
    DM x_opt = res.at("x");
    DM f_opt = res.at("f");
    DM g_opt = res.at("g");
    DM lam_g = res.at("lam_g");
    
    // Выводим результаты в удобном формате
    std::cout << "\n=== Результаты оптимизации ===" << std::endl;
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
    
    // Геометрическая интерпретация решения
    std::cout << "\n=== Геометрическая интерпретация ===" << std::endl;
    std::cout << "• Ограничение x₀ + x₁ + x₂ = 2 определяет плоскость в трехмерном пространстве" << std::endl;
    std::cout << "• Ограничение x₀² + x₁² ≤ 1.5 определяет цилиндр вдоль оси x₂" << std::endl;
    std::cout << "• Ограничение x₂ ≥ 0 отсекает отрицательную часть пространства по оси x₂" << std::endl;
    
    if (is_active_ineq1) {
        std::cout << "• Решение лежит на границе цилиндра (x₀² + x₁² = 1.5)" << std::endl;
    }
    
    if (is_active_ineq2) {
        std::cout << "• Решение лежит в плоскости x₂ = 0" << std::endl;
    }
    
    return 0;
}