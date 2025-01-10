
#include <cmath>
#include <Eigen/Dense>
#include <vector>

#include "const_param.h"

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

#pragma warning(disable : 4244)

using Eigen::Vector3d;
using Eigen::Vector4d;

double Date_to_JD(int day, int month, int year, int hour, int minute, int second) {

	int a = static_cast <int>((14 - month) / 12);
	int y = year + 4800 - a;
	int m = month + 12 * a - 3;

	int jdn = day + 365 * y - 32045
		+ static_cast <int>((153 * m + 2) / 5)
		+ static_cast <int>(y / 4)
		- static_cast <int>(y / 100)
		+ static_cast <int>(y / 400);

	double jd = static_cast<double>(jdn)
		+ static_cast<double>((hour - 12) / 24)
		+ static_cast<double>(minute / 1440)
		+ static_cast<double>(second / 86400);

	return jd;
}

double JD_start = Date_to_JD(11, 6, 2014, 15, 46, 16); // начало данных
double JD_0 = Date_to_JD(20, 3, 2014, 20, 57, 0); // день весеннего равноденствия для года выше
double JD_dif_sec = (JD_start - JD_0) * 86400; // временной промежуток в секундах

// функция поворота кватернионом q
Eigen::Vector3d Rotate_by_q(const Eigen::Vector4d& q, const Eigen::Vector3d& r) {
	Eigen::Vector3d x;

	x(0) = (q(0) * q(0) + q(1) * q(1) - q(2) * q(2) - q(3) * q(3)) * r(0) + 2 * (q(1) * q(2) - q(0) * q(3)) * r(1) + 2 * (q(1) * q(3) + q(0) * q(2)) * r(2);
	x(1) = 2 * (q(1) * q(2) + q(0) * q(3)) * r(0) + (q(0) * q(0) + q(2) * q(2) - q(1) * q(1) - q(3) * q(3)) * r(1) + 2 * (q(2) * q(3) - q(0) * q(1)) * r(2);
	x(2) = 2 * (q(1) * q(3) - q(0) * q(2)) * r(0) + 2 * (q(2) * q(3) + q(0) * q(1)) * r(1) + (q(0) * q(0) + q(3) * q(3) - q(1) * q(1) - q(2) * q(2)) * r(2);

	return x;
}

// функция перемножения кватернионов
Eigen::Vector4d Q_by_q(const Eigen::Vector4d& q1, const Eigen::Vector4d& q2) {
	Eigen::Vector4d q3;
	Eigen::Vector3d q1_vect(q1(1), q1(2), q1(3));
	Eigen::Vector3d q2_vect(q2(1), q2(2), q2(3));
	q3(0) = q1(0) * q2(0) - q1_vect.dot(q2_vect);
	Eigen::Vector3d q3_vect = q1(0) * q2_vect + q2(0) * q1_vect + q1_vect.cross(q2_vect);
	q3(1) = q3_vect(0);
	q3(2) = q3_vect(1);
	q3(3) = q3_vect(2);
	return q3;
}


class Unumbra {
public:

	explicit Unumbra(double& t, const double& unumbra_param)
		: t_(t)
		, unumbra_param_(unumbra_param)
	{

	}

	// Unumbra function
	double Unumbra_function() {

		double unumbra_funct = 0;

		double k1 = -sqrt(r * r - Re * Re) / r;
		double k2 = -sqrt((r + unumbra_param_) * (r + unumbra_param_) - Re * Re) / (r + unumbra_param_);

		Vector3d sun_vect_ECI = Sun_vector();

		PQW_in_ECI  pqw_eci = Get_current_PQW_in_ECI();

		Vector3d eZ = pqw_eci.Zpqw_ECI;

		double unumbra_index = sun_vect_ECI.dot(eZ);

		if (unumbra_index > k1) {
			unumbra_funct = 1;
		}
		else if ((unumbra_index <= k1) && (unumbra_index >= k2)) {
			unumbra_funct = (unumbra_index - k2) / (k1 - k2);
		}

		return unumbra_funct;
	}

private:

	double t_;
	const double unumbra_param_;

	/*

PQW:
x - по трансверсали
y - по нормали к плоскости орбиты
z - по радиус-вектору орбиты

ECI:
вторая геоэкваториальная СК
x - до правой тройки
y - ось вращения Земли
z - точка весеннего равнодействия

*/

	// оси pqw в eci
	struct PQW_in_ECI {
		Vector3d Xpqw_ECI;
		Vector3d Ypqw_ECI;
		Vector3d Zpqw_ECI;
	};

	PQW_in_ECI Get_current_PQW_in_ECI() {

		double u = u0 + t_ * wd;
		double Omega = Omega0 + Omega_p * t_;
		// используем активную точку зрения
		Vector4d q_1(cos(Omega / 2), 0, sin(Omega / 2), 0);
		Vector4d q_2(cos(incl / 2), sin(Omega) * sin(incl / 2), 0, cos(Omega) * sin(incl / 2));
		Vector4d q_3(cos(u / 2), -sin(incl) * cos(Omega) * sin(u / 2), cos(incl) * sin(u / 2), sin(incl) * sin(Omega) * sin(u / 2));

		// кватернион задающий положение осей орбитальной СК в неподвижной ECI
		Vector4d q_pqw_in_ECI = Q_by_q(q_3, Q_by_q(q_2, q_1));

		PQW_in_ECI current_pqw_in_eci;

		current_pqw_in_eci.Xpqw_ECI = Rotate_by_q(q_pqw_in_ECI, { 1, 0, 0 });
		current_pqw_in_eci.Ypqw_ECI = Rotate_by_q(q_pqw_in_ECI, { 0, 1, 0 });
		current_pqw_in_eci.Zpqw_ECI = Rotate_by_q(q_pqw_in_ECI, { 0, 0, 1 });

		return current_pqw_in_eci;
	}

	// вычисление солнечного вектора во второй геоэкваториальной СК
	Eigen::Vector3d Sun_vector() {

		Eigen::Vector3d sun_vector(sin(n_earth * (t_ + JD_dif_sec)) * cos(i_e), sin(n_earth * (t_ + JD_dif_sec)) * sin(i_e), cos(n_earth * (t_ + JD_dif_sec)));
		return sun_vector;
	}

};

int main()
{
	double t_end_days = 1;
	double t_end = t_end_days * 86400.0;
	const double dt = 10;

	double t = 0;
	std::vector<double> unumbra_f;
	std::vector<double> time;

	while (t != t_end) {
		Unumbra umbra(t, 1000);
		unumbra_f.push_back(umbra.Unumbra_function());
		time.push_back(t / 60);
		t += dt;
	}

	plt::plot(time, unumbra_f);
	plt::title("Unumbra function");
	plt::show();

}
