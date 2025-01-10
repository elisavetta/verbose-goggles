#pragma once

const double PI = acos(-1.0);
const double deg2rad = PI / 180;
const double years2sec = 365.25 * 24 * 3600;

// sun parameters
double n_earth = 2 * PI / (86400 * 365.25); // среднее движение земли по орбите
double i_e = 23.5 * deg2rad; //угол наклоненния эклиптики

//CONSTANTS:

// environment parameters----------------------------
const double Re = 6378.245e3;

//Earth gravity parameters----------------------------
const double vmu = 3.986e14;     //the Earth gravity constant
const double vJ2 = 1.082626e-3;  // J2 coefficient
const double re = Re;    // mean radius

// initial orbital values
const double h0_km = 1340;
const double i0_deg = 66;
const double Omega0_deg = 0;
const double u0_deg = 0;

// orbit parameters----------------------------
const double r = h0_km * 1000 + re;
const double r32 = sqrt(r * r * r);
const double r72 = r * r * r * sqrt(r);

// evolution of the orbit----------------------------
const double incl = i0_deg * deg2rad;
const double Omega0 = Omega0_deg * deg2rad;
const double u0 = u0_deg * deg2rad;
const double cos_incl = cos(incl);
const double sin_incl = sin(incl);

const double vmu12 = sqrt(vmu);
const double w0 = vmu12 / r32;

const double Omega_p = -1.5 * vJ2 * vmu12 * re * re * cos_incl / r72;
const double wd = w0 * (1 - 1.5 * vJ2 * re * re * (1 - 4 * cos_incl * cos_incl) / (r * r));