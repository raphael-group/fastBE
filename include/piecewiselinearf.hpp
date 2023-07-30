#ifndef PIECEWISELINEARF_HPP
#define PIECEWISELINEARF_HPP

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

class PiecewiseLinearF {
public:
    double intercept;
    std::vector<double> slopes;

    PiecewiseLinearF(const std::vector<double>& slopes, double intercept) 
    : intercept(intercept), slopes(slopes) {
    }

    PiecewiseLinearF(const std::vector<double>& slopes, double intercept, std::size_t capacity) 
      : intercept(intercept), slopes(slopes) {
        this->slopes.reserve(capacity);
    }

    PiecewiseLinearF operator+(const PiecewiseLinearF& other) const {
        std::vector<double> self_slopes = slopes, other_slopes = other.slopes;
        if (self_slopes.size() < other_slopes.size()) {
            self_slopes.resize(other_slopes.size(), self_slopes.back());
        } else {
            other_slopes.resize(self_slopes.size(), other_slopes.back());
        }
        std::transform(self_slopes.begin(), self_slopes.end(), other_slopes.begin(), self_slopes.begin(), std::plus<double>());
        return PiecewiseLinearF(self_slopes, intercept + other.intercept);
    }

    void addInPlace(PiecewiseLinearF&& other) {
        if (other.slopes.size() > slopes.size()) {
            std::swap(slopes, other.slopes);
        }

        for (size_t i = 0; i < other.slopes.size(); ++i) {
            slopes[i] += other.slopes[i];
        }

        for (size_t i = other.slopes.size(); i < slopes.size(); ++i) {
            slopes[i] += other.slopes.back();
        }

        intercept += other.intercept;
    }

    double operator()(double x) const {
        int index = std::min(static_cast<int>(std::floor(x)), static_cast<int>(slopes.size() - 1));
        std::vector<double> intercepts = this->intercepts();
        return slopes[index] * (x - index) + intercepts[index];
    }

    std::vector<double> intercepts() const {
        std::vector<double> intercepts;
        intercepts.resize(slopes.size());
        intercepts[0] = intercept;
        std::partial_sum(slopes.begin(), slopes.end() - 1, intercepts.begin() + 1);
        for (std::size_t i = 1; i < intercepts.size(); ++i) {
            intercepts[i] = intercepts[i] + intercept;
        }
        return intercepts;
    }

    double minimizer() const {
        double minval = intercept, sum = intercept;

        for (std::size_t i = 0; i < slopes.size(); ++i) {
            if (slopes[i] >= 0) break;
            sum += slopes[i];
            minval = std::min(minval, sum);
        }

        return minval;
    }

    void compute_minimizer() {
        int jl = -1;

        for (std::size_t j = 0; j < slopes.size(); j++) {
            if (slopes[j] >= 0) {
                jl = j;
                break;
            }
        }

        if (jl == -1) {
            intercept += slopes[0];
            if (slopes.size() == 1) return;
            slopes.erase(slopes.begin());
            return;
        }

        if (jl == 0) {
            slopes.insert(slopes.begin(), 0.0);
            return;
        }

        intercept += slopes[0];
        slopes.insert(slopes.begin() + jl, 2, 0.0);
        slopes.erase(slopes.begin());
    }
};

PiecewiseLinearF compute_minimizer(const PiecewiseLinearF& f) {
    PiecewiseLinearF g = f;
    g.compute_minimizer();
    return g;
}

#endif
