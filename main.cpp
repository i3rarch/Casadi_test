#include <iostream>
#include <casadi/casadi.hpp>
#include <windows.h>  // Add this include for SetConsoleOutputCP

using namespace casadi;

int main() {
    SetConsoleOutputCP(CP_UTF8);
    
    // Создаем символическую переменную (3-мерный вектор)
    MX x = MX::sym("x", 3);
    
    // Создаем трехмерную функцию (например, "чаша Розенброка" с дополнительным измерением)
    // f(x,y,z) = 100*(y-x^2)^2 + (1-x)^2 + (z-y^2)^2 + (1-y)^2
    MX f = 100*pow(x(1)-pow(x(0),2),2) + pow(1-x(0),2) + 
           100*pow(x(2)-pow(x(1),2),2) + pow(1-x(1),2);
    
    // Добавляем ограничения
    // Ограничение равенства: x0 + x1 + x2 = 2
    MX g_eq = x(0) + x(1) + x(2) - 2;
    
    // Ограничения неравенства: x0^2 + x1^2 <= 1.5, x2 >= 0
    MX g_ineq = MX::vertcat({pow(x(0),2) + pow(x(1),2) - 1.5, -x(2)});
    
    // Создаем объект для решения задачи оптимизации
    Dict opts;
    opts["ipopt.print_level"] = 3;  // Увеличиваем уровень вывода для информативности
    opts["print_time"] = 1;
    
    // Формируем задачу с ограничениями
    MXDict nlp = {
        {"x", x}, 
        {"f", f},
        {"g", MX::vertcat({g_eq, g_ineq})}  // Объединяем все ограничения в один вектор
    };
    
    Function solver = nlpsol("solver", "ipopt", nlp, opts);
    
    // Задаем начальную точку и границы
    std::vector<double> x0 = {0.5, 0.5, 1.0};  // Начальная точка
    std::vector<double> lbx = {-10.0, -10.0, -10.0};  // Нижние границы переменных
    std::vector<double> ubx = {10.0, 10.0, 10.0};     // Верхние границы переменных
    
    // Границы для ограничений
    // g_eq = 0 (равенство)
    // Первое неравенство g_ineq[0] <= 0
    // Второе неравенство g_ineq[1] <= 0
    std::vector<double> lbg = {0.0, -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    std::vector<double> ubg = {0.0, 0.0, 0.0};
    
    // Формируем входные аргументы для решателя
    DMDict arg = {
        {"x0", DM(x0)},
        {"lbx", DM(lbx)},
        {"ubx", DM(ubx)},
        {"lbg", DM(lbg)},
        {"ubg", DM(ubg)}
    };
    
    // Решаем задачу оптимизации
    DMDict res = solver(arg);
    
    // Получаем оптимальное решение
    DM x_opt = res.at("x");
    DM f_opt = res.at("f");
    DM g_opt = res.at("g");
    DM lam_g = res.at("lam_g");  // Множители Лагранжа для ограничений
    
    // Выводим результат
    std::cout << "Оптимальное решение x* = " << x_opt << std::endl;
    std::cout << "Минимальное значение f(x*) = " << f_opt << std::endl;
    std::cout << "Значения ограничений g(x*) = " << g_opt << std::endl;
    std::cout << "Множители Лагранжа для ограничений = " << lam_g << std::endl;
    
    // Проверяем выполнение ограничений
    std::cout << "\nПроверка ограничений:" << std::endl;
    std::cout << "x0 + x1 + x2 = 2: " << (fabs(g_opt(0).scalar() - 0.0) < 1e-6 ? "Выполнено" : "Не выполнено") << std::endl;
    std::cout << "x0^2 + x1^2 <= 1.5: " << (g_opt(1).scalar() <= 1e-6 ? "Выполнено" : "Не выполнено") << std::endl;
    std::cout << "x2 >= 0: " << (g_opt(2).scalar() <= 1e-6 ? "Выполнено" : "Не выполнено") << std::endl;
    
    return 0;
}