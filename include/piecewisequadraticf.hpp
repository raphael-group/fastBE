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

    PiecewiseQuadraticF() {}

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

    // runs in O(k) time
    double operator()(double x) const {
        // compute intercepts of the pieces of the derivative, using continuity
        std::vector<double> cs(breakpoints.size());
        cs[0] = c0 + m0 * (breakpoints[0]);
        for (size_t i = 1; i < breakpoints.size(); i++) {
            cs[i] = cs[i - 1] + slopes[i - 1] * (breakpoints[i] - breakpoints[i - 1]);
        }

        // evaluate at all breakpoints
        std::vector<double> values(breakpoints.size() + 1);
        values[0] = f0;
        if (breakpoints[0] >= 0) {
            values[0] += c0 * breakpoints[0] + 0.5 * m0 * breakpoints[0] * breakpoints[0];
        } else {
            values[0] -= c0 * breakpoints[0] + 0.5 * m0 * breakpoints[0] * breakpoints[0];
        }

        for (size_t i = 1; i < breakpoints.size(); i++) {
            values[i] = values[i - 1];
            values[i] += (cs[i - 1] - slopes[i - 1] * breakpoints[i - 1]) * (breakpoints[i] - breakpoints[i - 1]);
            values[i] += 0.5 * slopes[i - 1] * (breakpoints[i] * breakpoints[i] - breakpoints[i - 1] * breakpoints[i - 1]);
        }

        // the above can all be done in O(k) time and precomputed
        // now we find the piece that x is in
        size_t l = std::upper_bound(breakpoints.begin(), breakpoints.end(), x) - breakpoints.begin();
        if (l == 0) {
            return values[0] + (c0 * (x - breakpoints[0]) + 0.5 * m0 * (x * x - breakpoints[0] * breakpoints[0]));
        } else {
            return values[l - 1] + (cs[l - 1] - slopes[l - 1] * breakpoints[l - 1]) * (x - breakpoints[l - 1]) + 0.5 * slopes[l - 1] * (x * x - breakpoints[l - 1] * breakpoints[l - 1]);
        }
    }
};

void test_piecewisequadraticf() {
    PiecewiseQuadraticF f1(0.7), f2(0.3), f3(0.9), f4(1.0), f5(0.0);

    PiecewiseQuadraticF f123 = f1 + f2 + f3 + f4 + f5;

    for (double x = -2.0; x <= 2.0; x += 0.2) {
        // std::cout << "x = " << x << " f1(x) + f2(x) + f3(x) = " << f1(x) + f2(x) + f3(x) << std::endl;
        std::cout << "x = " << x << " f123(x) = " << f123(x) << std::endl;
    }
}

#endif
