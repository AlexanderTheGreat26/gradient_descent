#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <array>
#include <algorithm>
#include <sstream>
#include <tuple>


const double x_min = -1.5;
const double x_max = 4.0;
const double y_min = -3.0;
const double y_max = 4.0;
const double step = 0.01;
const double eps = std::numeric_limits<double>::epsilon();
const double l = 0.1;

typedef std::vector<double> grid1;
typedef std::tuple<double, double, double> coordinates;
typedef std::tuple<double, double> vector_coordinates;
typedef std::vector<coordinates> grid3;


std::random_device rd;  // Will be used to obtain a seed for the random number engine.
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd().


grid1 mesh (double left_border, const double & right_border, const double & mesh_step);

grid3 function_nodes (grid1 & x_grid, grid1 & y_grid, double f(double & x, double & y));

double kenny (double & x, double & y);

template<typename... Tp>
void file_creation (const std::string & file_name, std::vector<Tp...> & data);

void plot (const std::string & name, const std::string & title);

grid3 gradient_descent (coordinates q_1st, coordinates q_2nd, double f(double & x, double & y));




int main() {
    grid1 x_grid = std::move(mesh(x_min, x_max, step));
    grid1 y_grid = std::move(mesh(y_min, y_max, step));
    grid3 McCormick = std::move(function_nodes(x_grid, y_grid, kenny));
    file_creation ("function", McCormick);
    std::uniform_int_distribution<> dis_node (1, McCormick.size()-2);
    int i = dis_node(gen);
    int j = dis_node(gen);
    grid3 traj = gradient_descent(McCormick[i], McCormick[j], kenny);
    file_creation ("trajectory", traj);
    //plot("kenny_test", "McCormick function");
    return 0;
}


double kenny (double & x, double & y) {
    return std::sin(x + y) + std::pow(x - y, 2.0) - 1.5*x + 2.5*y + 1.0;
}


grid3 function_nodes (grid1 & x_grid, grid1 & y_grid, double f(double & x, double & y)) {
    grid3 z;
    for (double & x : x_grid)
        for (double & y : y_grid)
            z.emplace_back(std::make_tuple(x, y, f(x, y)));
    return z;
}


double diff (double & f_curr, double & f_prev, double dx) {
    if (std::fabs(dx) < std::numeric_limits<double>::epsilon())
        return 0;
    else
        return (f_curr - f_prev) / dx;
}


vector_coordinates grad (vector_coordinates & current_point, vector_coordinates & new_point, double f (double & x, double & y)) {
    double x_curr = std::get<0>(current_point);
    double y_curr = std::get<1>(current_point);
    double x_new = std::get<0>(new_point);
    double y_new = std::get<1>(new_point);
    double dx = x_new - x_curr;
    double dy = y_new - y_curr;
    double f_curr = f(x_curr, y_curr);
    double f_x_new = f(x_new, y_curr);
    double f_y_new = f(x_curr, y_new);
    return std::make_tuple(diff(f_x_new, f_curr, dx), diff(f_y_new, f_curr, dy));
}


template<size_t Is = 0, typename... Tp>
void new_point (std::tuple<Tp...> & point, std::tuple<Tp...> & function_gradient, const double & lambda) {
        std::get<Is>(point) -= lambda * std::get<Is>(function_gradient);
    if constexpr (Is + 1 != sizeof...(Tp))
        new_point<Is + 1>(point, function_gradient, lambda);
}


grid3 gradient_descent (coordinates q_1st, coordinates q_2nd, double f(double & x, double & y)) {
    double old_value, new_value = std::get<2>(q_1st);
    vector_coordinates current_point = std::make_tuple(std::get<0>(q_1st), std::get<1>(q_1st));
    vector_coordinates next_point = std::make_tuple(std::get<0>(q_2nd), std::get<1>(q_2nd));
    vector_coordinates grad_value = grad(current_point, next_point, f);
    grid3 trajectory {q_1st, q_2nd};
    do {
        int n = trajectory.size()-1;
        old_value = new_value;
        current_point = std::make_tuple(std::get<0>(trajectory[n]), std::get<1>(trajectory[n]));
        next_point = current_point;
        new_point(next_point, grad_value, l);
        grad_value = grad(current_point, next_point, f);
        new_value = f(std::get<0>(next_point), std::get<1>(next_point));
        trajectory.emplace_back(std::make_tuple(std::get<0>(next_point), std::get<1>(next_point), new_value));
    } while (std::fabs(old_value - new_value) > eps);
    return trajectory;
}


grid1 mesh (double left_border, const double & right_border, const double & mesh_step) {
    std::vector <double> xx ((right_border-left_border) / mesh_step + 1);
    xx[0] = left_border;
    std::generate(xx.begin()+1, xx.end(), [&] {left_border += step; return left_border;});
    return xx;
}


template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


template<typename T, size_t... Is>
std::string tuple2string_impl (T const& t, std::index_sequence<Is...>) {
    return ((toString(std::get<Is>(t)) + '\t') + ...);
}

template <class Tuple>
std::string tuple2string (const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return tuple2string_impl(t, std::make_index_sequence<size>{});
}


template<typename... Tp>
void file_creation (const std::string & file_name, std::vector<Tp...> & data) {
    std::ofstream fout;
    fout.open(file_name+".txt", std::ios::trunc);
    for (auto & i : data)
        fout << tuple2string(i) << '\n';
    fout.close();
}


void plot (const std::string & name, const std::string & title) {
    FILE *gp = popen("gnuplot -persist", "w");
    if (!gp) throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term jpeg size 1920, 1080 font \"Helvetica,30\"",
                                      "set output \'" + title + ".jpg\'",
                                      "set title \'" + title + R"(' font 'Helvetica,20')",
                                      "set grid xtics ytics ztics",
                                      "set autoscale xfix",
                                      "set key off",
                                      "set ticslevel 0",
                                      "set hidden3d",
                                      "set dgrid3d 50,50 qnorm 4",
                                      "splot \'" + name + ".txt\' w l",
                                      "set tics font \'Helvetica,14\'",
                                      R"(set xlabel 'x' font 'Helvetica,20')",
                                      R"(set ylabel 'y' font 'Helvetica,20')",
                                      R"(set zlabel 'z' font 'Helvetica,20')",
                                      "set terminal pop",
                                      "set output",
                                      "replot", "q"};
    for (const auto & it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}