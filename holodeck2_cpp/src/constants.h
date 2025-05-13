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

constexpr float ARCSEC     = 4.84813681E-06;
constexpr float AU         = 1.49597871E+13;
constexpr float DAY        = 8.64000000E+04;
constexpr float EDDT       = 6.32196208E+04;
constexpr float EVOLT      = 1.60217663E-12;
constexpr float GPC        = 3.08567758E+27;
constexpr float GYR        = 3.15576000E+16;
constexpr float HPLANCK    = 6.62607015E-27;
constexpr float JY         = 1.00000000E-23;
constexpr float KBOLTZ     = 1.38064900E-16;
constexpr float KMPERSEC   = 1.00000000E+05;
constexpr float KPC        = 3.08567758E+21;
constexpr float LSOL       = 3.82800000E+33;
constexpr float MELC       = 9.10938370E-28;
constexpr float MPC        = 3.08567758E+24;
constexpr float MPRT       = 1.67262192E-24;
constexpr float MSOL       = 1.98840987E+33;
constexpr float MYR        = 3.15576000E+13;
constexpr float NWTG       = 6.67430000E-08;
constexpr float PC         = 3.08567758E+18;
constexpr float QELC       = 4.80320471E-10;
constexpr float RSOL       = 6.95700000E+10;
constexpr float SCHW       = 1.48523205E-28;
constexpr float SIGMA_SB   = 5.67037442E-05;
constexpr float SIGMA_T    = 6.65245873E-25;
constexpr float SPLC       = 2.99792458E+10;
constexpr float TWO_PI     = 6.28318531E+00;
constexpr float YR         = 3.15576000E+07;