import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from dataclasses import dataclass

@dataclass
class PiecewiseLinear():
    slopes : np.ndarray
    intercepts : np.ndarray
    intercept : float 

    def __init__(self, slopes: np.ndarray, intercept : float):
        if not isinstance(slopes, np.ndarray):
            slopes = np.array(slopes)

        self.intercept = intercept 
        self.slopes = slopes
        self.intercepts = np.concatenate(([intercept], np.cumsum(self.slopes[:-1]) + intercept))

    def __add__(self, other):
        if len(self.slopes) < len(other.slopes):
            self_slopes = np.concatenate((self.slopes, self.slopes[-1] * np.ones(len(other.slopes) - len(self.slopes))))
            other_slopes = other.slopes
        else:
            self_slopes = self.slopes
            other_slopes = np.concatenate((other.slopes, other.slopes[-1] * np.ones(len(self.slopes) - len(other.slopes))))

        return PiecewiseLinear(self_slopes + other_slopes, self.intercept + other.intercept)

    def __call__(self, x):
        indices = np.minimum(np.floor(x).astype(int), len(self.slopes) - 1)
        return self.slopes[indices] * (x - indices) + self.intercepts[indices]

"""
Given a convex and continuous, piecewise linear function whose breakpoints
are regularly spaced from $1, \ldots, k$, returns the function
$g(\lambda) = \min_{x \in [\lambda - 1, \lambda + 1], x \geq 0} f(x)$
"""
def compute_minimizer(f : PiecewiseLinear) -> PiecewiseLinear:
    if np.all(f.slopes <= 0):
        if len(f.slopes) == 1:
            return PiecewiseLinear(f.slopes, f.intercept + f.slopes[0])

        return PiecewiseLinear(f.slopes[1:], f.intercept + f.slopes[0])

    # jl is the smallest index such that f.slopes[jl] >= 0
    jl = np.argmax(f.slopes >= 0) 
    if jl == 0:
        new_slopes = np.concatenate((np.array([0]), f.slopes))
        return PiecewiseLinear(new_slopes, f.intercept)

    # jh is the largest index such that f.slopes[jh] == 0
    jh = np.argmax(f.slopes == 0)
    if jh == 0: # i.e. no zero slopes
        jh = jl

    new_slopes = np.concatenate((f.slopes[1:jl], np.zeros(jh - jl + 2), f.slopes[jh:]))
    return PiecewiseLinear(new_slopes, f.intercept + f.slopes[0])

"""
Given a convex and continuous, piecewise linear function whose breakpoints
are regularly spaced from $1, \ldots, k$, evaluates the function
$g(\lambda) = \min_{x \in [\lambda - 1, \lambda + 1], x \geq 0} f(x)$
at a given point $\lambda$.
"""
def compute_minimizer_numeric(f, lam, precision=100):
    x = np.linspace(np.maximum(lam - 1, 0), lam + 1, precision)
    return np.min(f(x))

def plot(f : PiecewiseLinear, ax, label):
    x = np.linspace(0, 6, 1000)
    y = f(x)

    sns.lineplot(x=x, y=y, ax=ax, label=label)
    ax.legend(loc='upper right')

def plot_minimizer_numeric(f : PiecewiseLinear, ax, label):
    x = np.linspace(0, 6, 1000)
    y = np.array([compute_minimizer_numeric(f, lam) for lam in x])

    sns.lineplot(x=x, y=y, ax=ax, label=label)
    ax.legend(loc='upper right')

if __name__ == "__main__":
    f1 = PiecewiseLinear([-1.5, -1.0, 0.1, 1, 3], 1)
    f2 = PiecewiseLinear([-3, -2], 1)
    f3 = PiecewiseLinear([3, 4], 1)

    fig, axes = plt.subplots(1, 4, figsize=(15, 7))

    plot(f1, axes[0], "f1")
    plot_minimizer_numeric(f1, axes[0], "g1_approx")
    plot(compute_minimizer(f1), axes[0], "g1_exact")
    print(f1)
    print(compute_minimizer(f1))

    plot(f2, axes[1], "f2")
    plot(compute_minimizer(f2), axes[1], "g2_exact")
    plot_minimizer_numeric(f2, axes[1], "g2_approx")
    print(f2)
    print(compute_minimizer(f2))

    plot(f3, axes[2], "f3")
    plot(compute_minimizer(f3), axes[2], "g3_exact")
    plot_minimizer_numeric(f3, axes[2], "g3_approx")
    print(f3)
    print(compute_minimizer(f3))

    plot(f1 + f2 + f3, axes[3], "f1 + f2 + f3")
    print(f1 + f2 + f3)
    print(compute_minimizer(f1 + f2 + f3))

    plt.show()
