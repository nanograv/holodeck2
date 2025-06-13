/**
 * Numerical Constants
 *
 * All constants are used in CGS units, as raw floats.  Most of the holodeck package works in CGS units
 * whenever possible.  Constants and units should only be added when they are frequently used
 * (i.e. in multiple files/submodules).
 *
 * Notes
 * -------------------
 * - [cm] = centimeter
 * - [g] = gram
 * - [s] = second
 * - [erg] = cm^2 * g / s^2
 * - [Jy] = jansky = [erg/s/cm^2/Hz]
 * - [fr] franklin = statcoulomb = electro-static unit [esu]
 * - [K] Kelvin
 *
 */

#pragma once


// ---- Magic Numbers

constexpr double PI         = 3.14159265E+00;
constexpr double TWO_PI     = 6.28318531E+00;
constexpr double LOG10      = 2.30258509E+00;

// ---- Constants

constexpr double EDDT       = 6.32196208E+04;   //
constexpr double NWTG       = 6.67430000E-08;   // [cm^3/g/sec^2] Newton's gravitational constant
constexpr double SCHW       = 1.48523205E-28;   //
constexpr double SPLC       = 2.99792458E+10;   // [cm/sec] speed of light in vacuum

// ---- Units

constexpr double ARCSEC     = 4.84813681E-06;
constexpr double AU         = 1.49597871E+13;
constexpr double DAY        = 8.64000000E+04;
constexpr double KMPERSEC   = 1.00000000E+05;
constexpr double LSOL       = 3.82800000E+33;
constexpr double MPRT       = 1.67262192E-24;
constexpr double MSOL       = 1.98840987E+33;
constexpr double PC         = 3.08567758E+18;
constexpr double RSOL       = 6.95700000E+10;
constexpr double YR         = 3.15576000E+07;

// ---- Derived

constexpr double MYR        = 1.0E6 * YR;
constexpr double GYR        = 1.0E9 * YR;

constexpr double KPC        = 1.0E3 * PC;
constexpr double MPC        = 1.0E6 * PC;
constexpr double GPC        = 1.0E9 * PC;

constexpr double KM_S_MPC   = 1.0E5 / MPC;    // 1 km/s/Mpc [1/sec]


// constexpr double HPLANCK    = 6.62607015E-27;
// constexpr double JY         = 1.00000000E-23;
// constexpr double KBOLTZ     = 1.38064900E-16;
// constexpr double EVOLT      = 1.60217663E-12;    // Electron Volt
// constexpr double MELC       = 9.10938370E-28;
// constexpr double QELC       = 4.80320471E-10;
// constexpr double SIGMA_SB   = 5.67037442E-05;
// constexpr double SIGMA_T    = 6.65245873E-25;

