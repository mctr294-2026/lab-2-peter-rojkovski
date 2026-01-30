#include <iostream>
#include <cstdlib>
#include <cmath> // for getting the abs of floating-point numbers
//#include <iostream> // for personal debugging
#include "roots.hpp"

// A double is similar to a float, but can hold twice the space (64 bit vs 32 bit) therefore more precision and a greater range
// This is how close we need to get to 0
const double tolerance = 1e-6;
// Max number of iterations before giving up
const int iter_max = 1e6;

/* Tries to find a zero crossing in f() in the interval [a,b] with the bisection method
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found. Finding a root
 * is only guarenteed if f is continous within the interval and
 * a & b have opposite signs
*/
bool bisection(std::function<double(double)> f,
    double a, double b,
    double *root) {

    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0){
        return false;
    }

    if (std::abs(fa) <= tolerance){
        *root = a;
        return true;
    }

    if (std::abs(fb) <= tolerance){
        *root = b;
        return true;
    }

// For loop: initiator; condition; expression
    for (int i = 0; i < iter_max; i++){

        fa = f(a);
        fb = f(b);
        double c = (a+b)/2;
        double fc = f(c);

        if (std::abs(fc) <= tolerance){
            *root = c;
            return true;
        }

        if (fc * fa > 0){
            a = c;
        }
        else{
            b = c;
        }

    }

    return false;
}

/* Tries to find a zero crossing in f() in the interval [a,b] with the
 * false positive / regula falsi method
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found. Finding a root
 * is only guarenteed if f is continous within the interval and
 * a & b have opposite signs
 */

bool regula_falsi(std::function<double(double)> f,
    double a, double b,
    double *root){

    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0){
        return false;
    }

    if (std::abs(fa) <= tolerance){
        *root = a;
        return true;
    }

    if (std::abs(fb) <= tolerance){
        *root = b;
        return true;
    }
    for (int i = 0; i < iter_max; i++){

        fa = f(a);
        fb = f(b);
        double c = a - (fa * (b - a) / (fb - fa));
        double fc = f(c);

        if (std::abs(fc) <= tolerance){
            *root = c;
            return true;
        }

        if (fc * fa < 0){
            b = c;
        }
        else{
            a = c;
        }
    }
    return false;
}

/* Tries to find a zero crossing in f() in the interval [a,b] with
 * the netwon-raphson method, given a function that computes the
 * derivative g() and a starting guess c.
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found, which can happen
 * if iteration leaves the interval, or the derivative is zero.
 */
bool newton_raphson(std::function<double(double)> f,
    std::function<double(double)> g,
    double a, double b, double c,
    double *root) {

    double fa = f(a);
    double fb = f(b);

    if (std::abs(fa) <= tolerance){
        *root = a;
        return true;
    }

    if (std::abs(fb) <= tolerance){
        *root = b;
        return true;
    }

    for (int i = 0; i < iter_max; i++){

        double fc = f(c);
        double gc = g(c);

        if (std::abs(fc) <= tolerance){
            *root = c;
            return true;
        }

        if (gc == 0){
            return false;
        }

        c = c - fc/gc;

    }
    return false;
}

/* Tries to find a zero crossing in f() in the interval [a,b] with
 * the secant method, given a starting guess c.
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found, which can happen
 * if iteration leaves the interval, or derivative is zero.
 */
bool secant(std::function<double(double)> f,
    double a, double b, double c,
    double *root) {
    
    double fa = f(a);
    double fb = f(b);
    double c1 = a;

    if (std::abs(fa) <= tolerance){
        *root = a;
        return true;
    }

    if (std::abs(fb) <= tolerance){
        *root = b;
        return true;
    }

    for (int i = 0; i < iter_max; i++){

        double fc = f(c);
        double fc1 = f(c1);

        if (std::abs(fc1) <= tolerance){
            *root = c1;
            return true;
        }

        if (f(c1) - f(c) == 0){
            return false;
        }

        double c_temp = c1;

        c1 = c1 - f(c1) * (c1 - c) / (f(c1) - f(c));

        c = c_temp;

    }
    return false;
}