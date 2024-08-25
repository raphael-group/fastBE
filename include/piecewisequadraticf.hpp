#ifndef PIECEWISEQUADRATICF_HPP
#define PIECEWISEQUADRATICF_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

class PiecewiseQuadraticF {
public:
    // b_0 = -\infty, b_{k+1} = \infty, so are not stored
    double f0, c0, m0; 
    std::vector<double> slopes;      // m_1 <= ... <= m_k
    std::vector<double> breakpoints; // b_1 <= ... <= b_k

    PiecewiseQuadraticF() {
        f0 = 0.0;
        c0 = 0.0;
        m0 = 0.0;
    }

    // leaf constructor
    PiecewiseQuadraticF(double frequency) {
        breakpoints.push_back(2.0 * frequency);
        slopes.push_back(0.0);
        f0 = 0.0;
        c0 = frequency;
        m0 = -0.5;
    }

    // runs in O(k + k') time
    PiecewiseQuadraticF operator+(const PiecewiseQuadraticF& other) const {
        PiecewiseQuadraticF result;
        result.f0 = f0 + other.f0;
        result.c0 = c0 + other.c0;
        result.m0 = m0 + other.m0;

        std::vector<double> merged_breakpoints;
        std::vector<double> merged_slopes;

        merged_breakpoints.reserve(breakpoints.size() + other.breakpoints.size());
        merged_slopes.reserve(slopes.size() + other.slopes.size());

        for (size_t i = 0, j = 0; i < breakpoints.size() || j < other.breakpoints.size();) {
            if (j == other.breakpoints.size() || (i != breakpoints.size() && breakpoints[i] < other.breakpoints[j])) {
                merged_breakpoints.push_back(breakpoints[i]);
                i++;
            } else if (i == breakpoints.size() || (j != other.breakpoints.size() && breakpoints[i] > other.breakpoints[j])) {
                merged_breakpoints.push_back(other.breakpoints[j]);
                j++;
            } else {
                if (i == breakpoints.size()) {
                    merged_breakpoints.push_back(other.breakpoints[j]);
                } else {
                    merged_breakpoints.push_back(breakpoints[i]);
                }

                i++;
                j++;
            }

            double slope = 0;
            if (i == 0) {
                slope += m0;
            } else {
                slope += slopes[i - 1];
            }

            if (j == 0) {
                slope += other.m0;
            } else {
                slope += other.slopes[j - 1];
            }

            merged_slopes.push_back(slope);
        }

        result.breakpoints = std::move(merged_breakpoints);
        result.slopes = std::move(merged_slopes);
        return result;
    }

    std::vector<double> get_derivative_intercepts() const {
        // compute intercepts of the pieces of the derivative, using continuity
        std::vector<double> cs(breakpoints.size());
        cs[0] = c0 + m0 * (breakpoints[0]);
        for (size_t i = 1; i < breakpoints.size(); i++) {
            cs[i] = cs[i - 1] + slopes[i - 1] * (breakpoints[i] - breakpoints[i - 1]);
        }

        return cs; 
    }

    std::vector<double> evaluate_at_breakpoints() const {
        std::vector<double> cs = get_derivative_intercepts(); 

        // evaluate at all breakpoints, find interval that contains 0.0 and start there
        size_t l = std::lower_bound(breakpoints.begin(), breakpoints.end(), 0.0) - breakpoints.begin();

        std::vector<double> values(breakpoints.size() + 1);
        if (l == 0) {
            values[l] = f0 + c0 * breakpoints[l] + 0.5 * m0 * (breakpoints[l] * breakpoints[l]);
        } else {
            values[l] = f0 + (cs[l - 1] - slopes[l - 1] * breakpoints[l - 1]) * (breakpoints[l]) + 0.5 * slopes[l - 1] * (breakpoints[l] * breakpoints[l]);
        }

        for (size_t i = l + 1; i < breakpoints.size(); i++) {
            values[i] = values[i - 1];
            values[i] += (cs[i - 1] - slopes[i - 1] * breakpoints[i - 1]) * (breakpoints[i] - breakpoints[i - 1]);
            values[i] += 0.5 * slopes[i - 1] * (breakpoints[i] * breakpoints[i] - breakpoints[i - 1] * breakpoints[i - 1]);
        }

        for (size_t i = l; i > 0; i--) {
            values[i - 1] = values[i];
            values[i - 1] -= (cs[i - 1] - slopes[i - 1] * breakpoints[i - 1]) * (breakpoints[i] - breakpoints[i - 1]);
            values[i - 1] -= 0.5 * slopes[i - 1] * (breakpoints[i] * breakpoints[i] - breakpoints[i - 1] * breakpoints[i - 1]);
        }

        return values;
    }

    // runs in O(k) time
    double operator()(double x) const {
        // the above can all be done in O(k) time and precomputed
        std::vector<double> cs = get_derivative_intercepts();
        std::vector<double> values = evaluate_at_breakpoints();

        // now we find the piece that x is in
        size_t l = std::lower_bound(breakpoints.begin(), breakpoints.end(), x) - breakpoints.begin();
        if (l == 0) {
            return values[0] + (c0 * (x - breakpoints[0]) + 0.5 * m0 * (x * x - breakpoints[0] * breakpoints[0]));
        } else {
            return values[l - 1] + (cs[l - 1] - slopes[l - 1] * breakpoints[l - 1]) * (x - breakpoints[l - 1]) + 0.5 * slopes[l - 1] * (x * x - breakpoints[l - 1] * breakpoints[l - 1]);
        }
    }

    // when F = \sum{j \in \delta(i)}J_j, this updates F to be 
    // J_i(\gamma) = max_{x \geq 0}(h_i(x - \gamma) + F(x))
    // really, this is the meat of the algorithm
    PiecewiseQuadraticF update_representation(float frequency) const {
        // compute intercepts of the pieces of the derivative, using continuity
        std::vector<double> cs = get_derivative_intercepts();

        // find first breakpoint x
        size_t l = std::lower_bound(breakpoints.begin(), breakpoints.end(), 0.0) - breakpoints.begin();
        float x = 0.0;
        if (l == 0) {
            x = frequency - c0;
        } else {
            x = frequency + (0.5 * breakpoints[l-1]) - cs[l - 1] + (breakpoints[l-1] * (slopes[l-1] - 0.5));
        }

        x *= 2.0;

        std::vector<double> zs(breakpoints.size());
        for (size_t i = 0; i < breakpoints.size(); i++) {
            zs[i] = 2.0 * (frequency + 0.5 * breakpoints[i] - cs[i]);
        }

        std::vector<double> new_breakpoints;
        new_breakpoints.push_back(x);
        for (size_t i = l; i < breakpoints.size(); i++) {
            new_breakpoints.push_back(zs[i]);
        }

        std::vector<double> new_slopes;
        for (size_t i = 0; i < new_breakpoints.size(); i++) {
            double slope;
            if (i + l == 0) {
                slope = m0;
            } else {
                slope = slopes[i + l - 1];
            }

            new_slopes.push_back(-(slope / (2*slope - 1)));
        }
        
        double new_m0 = -0.5;
        double new_c0 = frequency;

        l = std::lower_bound(zs.begin(), zs.end(), 0.0) - zs.begin();
        double alpha_star = 0.0;

        if (l == 0) {
            alpha_star = (frequency - c0) / (m0 - 0.5);
        } else {
            alpha_star = (frequency + 0.5 * breakpoints[l-1] - cs[l-1]) / (slopes[l-1] - 0.5) + breakpoints[l-1];
        }

        alpha_star = std::max(0.0, alpha_star);
        double new_f0 = this->operator()(alpha_star) - (0.25 * alpha_star * alpha_star + frequency * alpha_star);
        
        PiecewiseQuadraticF result;
        result.f0 = new_f0;
        result.c0 = new_c0;
        result.m0 = new_m0;
        result.breakpoints = new_breakpoints;
        result.slopes = new_slopes;

        return result;
    }
};

void test_piecewisequadraticf() {
    /*
    PiecewiseQuadraticF f1(0.7), f2(0.3), f3(0.9), f4(1.0), f5(0.0);

    PiecewiseQuadraticF f123 = f1 + f2 + f3 + f4 + f5;

    for (double x = -2.0; x <= 2.0; x += 0.2) {
        // std::cout << "x = " << x << " f1(x) + f2(x) + f3(x) = " << f1(x) + f2(x) + f3(x) << std::endl;
        std::cout << "x = " << x << " f123(x) = " << f123(x) << std::endl;
    }

    PiecewiseQuadraticF res = f123.update_representation(0.5);
    std::cout <<" BUILT ASDSADAAAAAAAAAAAAAAAAAa" << std::endl;
    for (double x = -2.0; x <= 2.0; x += 0.2) {
        std::cout << "x = " << x << " res(x) = " << res(x) << std::endl;
    }
    */

    PiecewiseQuadraticF g(0.3);
    auto J1 = g.update_representation(1.0).update_representation(0.5);
    
    // print breakpoints, slopes, and values
    auto J1_cs = J1.get_derivative_intercepts();
    std::cout << J1.f0 << " " << J1.c0 << " " << J1.m0 << std::endl;
    for (size_t i = 0; i < J1.breakpoints.size(); i++) {
        std::cout << "b_" << i << " = " << J1.breakpoints[i] << " c_" << i << " = " << J1_cs[i] << " m_" << i << " = " << J1.slopes[i] << std::endl;
    }

    //auto J2 = J1.update_representation(0.0);
    //auto J2_cs = J2.get_derivative_intercepts();
    //auto J2_vals = J2.evaluate_at_breakpoints();
    //
    //std::cout << J2.f0 << " " << J2.c0 << " " << J2.m0 << std::endl;
    //for (size_t i = 0; i < J2.breakpoints.size(); i++) {
        //////// std::cout << "b_" << i << " = " << J2.breakpoints[i] << " c_" << i << " = " << J2_cs[i] << " m_" << i << " = " << J2.slopes[i] << std::endl;
        //std::cout << "f(b_" << i << ") = " << J2_vals[i] << std::endl;
    //}

    for (double x = -2.0; x <= 2.0; x += 0.2) {
        std::cout << "x = " << x << " J1(x) = " << J1(x) << std::endl;
    }
    
}

#endif
