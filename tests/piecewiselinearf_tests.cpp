#include <gtest/gtest.h>
#include <piecewiselinearf.hpp>

// Operator overload for ==
bool operator==(const PiecewiseLinearF& lhs, const PiecewiseLinearF& rhs) {
    // You should define here the condition for PiecewiseLinearF objects to be equal.
    // This is a placeholder code and will not work for your specific case.
    return lhs.slopes == rhs.slopes && lhs.intercept == rhs.intercept;
}

TEST(PiecewiseLinearF, TestComputeMinimizerInPlace) {
    PiecewiseLinearF f1({-1.5, -1.0, 0.1, 1, 3}, 1);
    f1.compute_minimizer();
    ASSERT_EQ(f1, PiecewiseLinearF({-1. ,  0. ,  0. ,  0.1,  1. ,  3.}, -0.5));

    PiecewiseLinearF f2({-3, -2}, 1);
    f2.compute_minimizer();
    ASSERT_EQ(f2, PiecewiseLinearF({-2.}, -2));

    PiecewiseLinearF f3({3, 4}, 1);
    f3.compute_minimizer();
    ASSERT_EQ(f3, PiecewiseLinearF({0, 3., 4}, 1));
}

TEST(PiecewiseLinearF, TestComputeMinimizer) {
    PiecewiseLinearF f1({-1.5, -1.0, 0.1, 1, 3}, 1);
    ASSERT_EQ(compute_minimizer(f1), PiecewiseLinearF({-1. ,  0. ,  0. ,  0.1,  1. ,  3.}, -0.5));

    PiecewiseLinearF f2({-3, -2}, 1);
    ASSERT_EQ(compute_minimizer(f2), PiecewiseLinearF({-2.}, -2));
    
    PiecewiseLinearF f3({3, 4}, 1);
    ASSERT_EQ(compute_minimizer(f3), PiecewiseLinearF({0, 3., 4}, 1));
}

TEST(PiecewiseLinearF, TestAdditionAndComputeMinimizer) {
    PiecewiseLinearF f1({-1.5, -1.0, 0.1, 1, 3}, 1);
    PiecewiseLinearF f2({-3, -2}, 1);
    PiecewiseLinearF f3({3, 4}, 1);

    ASSERT_EQ(f1 + f2 + f3, PiecewiseLinearF({-1.5,  1. ,  2.1,  3. ,  5.}, 3));
    ASSERT_EQ(compute_minimizer(f1 + f2 + f3), PiecewiseLinearF({0. , 0. , 1. , 2.1, 3. , 5.}, 1.5));
}

TEST(PiecewiseLinearF, TestAdditionAndComputeMinimizerInPlace) {
    PiecewiseLinearF f1({-1.5, -1.0, 0.1, 1, 3}, 1);
    PiecewiseLinearF f2({-3, -2}, 1);
    PiecewiseLinearF f3({3, 4}, 1);

    f1.addInPlace(std::move(f2));
    f1.addInPlace(std::move(f3));

    // print out f1 slopes
    for (auto slope : f1.slopes) {
        std::cout << slope << " ";
    }
    std::cout << std::endl;

    ASSERT_EQ(f1, PiecewiseLinearF({-1.5,  1. ,  2.1,  3. ,  5.}, 3));

    f1.compute_minimizer();
    ASSERT_EQ(f1, PiecewiseLinearF({0. , 0. , 1. , 2.1, 3. , 5.}, 1.5));
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
