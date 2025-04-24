#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>

void write_vtk(std::vector<std::vector<double>> mesh, std::string file_path, int time_step, double dx) {
    std::filesystem::create_directory(file_path);
    std::filesystem::path f_name{ "step_" + std::to_string(time_step) + ".vtk" };
    f_name = file_path / f_name;

    std::ofstream ofs{ f_name };
    size_t Nx{ mesh.size() }, Ny{ mesh.at(0).size() };

    ofs << "# vtk DataFile Version 3.0\n";
    ofs << f_name.string() << std::endl;
    ofs << "ASCII\n";
    ofs << "DATASET STRUCTURED_GRID\n";

    ofs << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << "\n";
    ofs << "POINTS " << Nx * Ny * 1 << " float\n";

    for (size_t i = 0; i < Nx; i++) {
        for (size_t j = 0; j < Ny; j++) {
            ofs << (double)i * dx << " " << (double)j * dx << " " << 1 << std::endl;
        }
    }
    ofs << "POINT_DATA " << Nx * Ny * 1 << std::endl;

    ofs << "SCALARS " << "CON " << "float 1\n";
    ofs << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < Nx; i++) {
        for (size_t j = 0; j < Ny; j++) {
            ofs << mesh.at(i).at(j) << std::endl;
        }
    }

    ofs.close();
}




const int Nx = 300;
const double dx =0.03 , dy = 0.03;
const int nstep = 4000, pstep = 50;
const double dt = 1.0e-4, tau = 0.0003,  mu = 1.0, kappa = 1.8;
const double pi = 4.0 * atan(1.0);

double laplacian(double cl, double cr, double cd, double cu, double cc, double dx) {
    return (cl + cr + cd + cu - 4.0 * cc) / (dx * dx);
}

double f(double phi, double T) {
    double gamma = 10.0, T_eq = 1.0, alpha = 0.9;
    double m = alpha / pi * atan(gamma * (T_eq - T));
    return 0.25 * pow(phi, 4) - (0.5 - m / 3.0) * pow(phi, 3) + (0.25 - m / 2.0) * pow(phi, 2);
}

double dphi_dx(double phir, double phicl, double dx) {
    return (phir - phicl) / (2*dx);
}
double dphi_dy(double phiu, double phil, double dy) {
    return (phiu - phil) /(2*dy);
}

double epsilon( double dphi_dx, double dphi_dy) {
    // 计算梯度方向的角度
    double jj = 4,  theta_0 = 0.2,  delta = 0.02,  epsilon_0 = 0.01;
    double theta = atan2(dphi_dy, dphi_dx);  // atan2 可避免除以零的情况

    // 计算 epsilon
    return epsilon_0 * (1 + delta * cos(jj * (theta - theta_0)));
}
double df_dphi(double phi, double T) {
    double gamma = 10.0, T_eq = 1.0, alpha = 0.9;
    double m = alpha / pi * atan(gamma * (T_eq - T));
    return phi * (1 - phi) * (phi - 0.5 + m);
}
double depisilon_dtheta(double dphi_dx, double dphi_dy) {
    double jj = 4,  theta_0 = 0.2,  delta = 0.02,  epsilon_0 = 0.01;
    double theta = atan2(dphi_dy, dphi_dx);
    return epsilon_0 * (jj * delta * sin(jj * (theta_0 - theta)));
}



int main() {
    //设置四个网格
    std::vector<std::vector<double>> phi_mesh(300, std::vector<double>(300, 0));
    std::vector<std::vector<double>> temp_mesh(300, std::vector<double>(300, 0));
    std::vector<std::vector<double>> custom_1(300, std::vector<double>(300, 0));
    std::vector<std::vector<double>> custom_2(300, std::vector<double>(300, 0));
    std::vector<std::vector<double>> dphi_dt_mesh(300, std::vector<double>(300, 0));



    //初始化网格
    for (int i = 0; i < 300; i++) {
        for (int j = 0; j < 300; j++) {
            if (
                (i - 150) * (i - 150) + (j - 150) * (j - 150) < 5 * 5) {
                phi_mesh.at(i).at(j) = 1.0;

            }
            temp_mesh.at(i).at(j) = 0.0;
        }
    }

    for (int istep = 0; istep < nstep; istep++) {

        for (int i = 0; i < 300; i++) {
            for (int j = 0; j < 300; j++) {
                int im = i - 1;
                if (im == -1)
                    im = Nx - 1;
                int ip = i + 1;
                if (ip == Nx)
                    ip = 0;
                int jm = j - 1;
                if (jm == -1)
                    jm = Nx - 1;
                int jp = j + 1;
                if (jp == Nx)
                    jp = 0;
                double phil{ phi_mesh.at(im).at(j) };
                double phir{ phi_mesh.at(ip).at(j) };
                double phid{ phi_mesh.at(i).at(jm) };
                double phiu{ phi_mesh.at(i).at(jp) };
                double phic{ phi_mesh.at(i).at(j) };
                double templ{ temp_mesh.at(im).at(j) };
                double tempr{ temp_mesh.at(ip).at(j) };
                double tempd{ temp_mesh.at(i).at(jm) };
                double tempu{ temp_mesh.at(i).at(jp) };
                double tempc{ temp_mesh.at(i).at(j) };

                // calculate gradient of phi
                double dphi_dx_n = dphi_dx(phir, phil, dx);
                double dphi_dy_n = dphi_dy(phiu, phid, dy);
                // calculate epsilon
                // derivative of epsilon 
                double epsilon_n = epsilon( dphi_dx_n, dphi_dy_n);
                double depisilon_dtheta_n = depisilon_dtheta(dphi_dx_n, dphi_dy_n);
                // take gradient of phi
                double lap_phi = laplacian(phil, phir, phid, phiu, phic, dx);

                // update custom_1 and custom_2
                double& dterm_1{ custom_1.at(i).at(j) };
                double& dterm_2{ custom_2.at(i).at(j) };
                dterm_1 = epsilon_n * depisilon_dtheta_n * dphi_dx_n;
                dterm_2 = epsilon_n * depisilon_dtheta_n * dphi_dy_n;
                // third term calculation
                double term_3 = epsilon_n * epsilon_n * lap_phi;
                // forth term calculation
                double T_n{ temp_mesh.at(i).at(j) };
                double term_4 = df_dphi(phic, T_n);
                // laplacian of temperature -> T_first term
                phi_mesh.at(i).at(j) += ((term_3 + term_4) * dt) / tau;
                double lap_T = laplacian(templ, tempr, tempd, tempu, tempc, dx);
                temp_mesh.at(i).at(j) += lap_T * dt;
                dphi_dt_mesh.at(i).at(j) = (term_3 + term_4) / tau;

            }
        }
        for (int i = 0; i < 300; i++) {
            for (int j = 0; j < 300; j++) {
                int im = i - 1;
                if (im == -1)
                    im = Nx - 1;
                int ip = i + 1;
                if (ip == Nx)
                    ip = 0;
                int jm = j - 1;
                if (jm == -1)
                    jm = Nx - 1;
                int jp = j + 1;
                if (jp == Nx)
                    jp = 0;

                // gradient of custom_1 -> first term
                // gradient of custom_2 -> second term
                double term1_u{ custom_1.at(i).at(jp) };
                double term1_d{ custom_1.at(i).at(jm) };
                double term2_r{ custom_2.at(ip).at(j) };
                double term2_l{ custom_2.at(im).at(j) };
                double term_1 = dphi_dy(term1_u, term1_d, dy);
                double term_2 = dphi_dx(term2_r, term2_l, dx);
                double term_1_2 = ((term_1 - term_2) * dt) / tau;
                phi_mesh.at(i).at(j) += term_1_2;
                double dphi_dt = dphi_dt_mesh.at(i).at(j) + (term_1 - term_2) / tau;
                temp_mesh.at(i).at(j) += kappa * dphi_dt * dt;
                


            }
        }


        if (istep % pstep == 0) {

            for (int i = 0; i < 300; i++) {
                for (int j = 0; j < 300; j++) {
                    if (std::abs(phi_mesh.at(i).at(j)) < 1e-6) {
                        phi_mesh.at(i).at(j) = 0.0;
                    }
                }
            }
            
            write_vtk(phi_mesh, "./result", istep, dx);

        }

    }

}