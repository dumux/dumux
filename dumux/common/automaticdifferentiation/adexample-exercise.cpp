/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.
  License: CC-BY 4.0
*/

#include <cmath>
#include <stdexcept>
#include <iostream>

namespace Ad
{
    class SimpleAd
    {
    public:
        SimpleAd(double val, double deriv)
            : val_(val), deriv_(deriv)
        {}
        double value() const
        {
            return val_;
        }
        double derivative() const
        {
            return deriv_;
        }
        SimpleAd operator+(const SimpleAd& other) const
        {
            return { val_ + other.val_, deriv_ + other.deriv_ };
        }
        SimpleAd operator*(const SimpleAd& other) const
        {
            return { val_ * other.val_, deriv_ * other.val_ + val_ * other.deriv_ };
        }
        SimpleAd operator/(const SimpleAd& other) const
        {
            return { val_ / other.val_, (deriv_ * other.val_ - val_ * other.deriv_)/(other.val_*other.val_) };
        }
    private:
        double val_;
        double deriv_;
    };

    SimpleAd operator+(const SimpleAd& left, const double right)
    {
        return { left.value() + right, left.derivative() };
    }

    SimpleAd operator+(const double left, const SimpleAd& right)
    {
        return right + left;
    }

    SimpleAd operator*(const SimpleAd& left, const double right)
    {
        return { left.value() * right, left.derivative() * right };
    }

    SimpleAd operator*(const double left, const SimpleAd& right)
    {
        return right * left;
    }

    SimpleAd operator/(const SimpleAd& left, const double right)
    {
        return { left.value() / right, left.derivative() / right };
    }

    SimpleAd operator/(const SimpleAd& left, const double right)
    {
        return { left.value() / right, left.derivative() / right };
    }

    SimpleAd operator/(const double left, const SimpleAd& right)
    {
        return right / left;
    }

    SimpleAd sin(const SimpleAd& x)
    {
        return { std::sin(x.value()), std::cos(x.value()) * x.derivative() };
    }

    SimpleAd cos(const SimpleAd& x)
    {
        return { std::cos(x.value()), -std::sin(x.value()) * x.derivative() };
    }
}

template <class F>
double newtonUpdate(F f, double x)
{
    std::cout << "newtonUpdate(), x = " << x << std::endl;
    Ad::SimpleAd adx(x, 1.0);
    auto fx = f(adx);
    return x - fx.value() / fx.derivative();
}

template <class F>
std::pair<double, int> newton(F f, double start)
{
    const double tol = 1e-13;
    const int max_iter = 20;
    int iter = 0;
    double x = start;
    while (fabs(f(x)) > tol && iter < max_iter) {
        x = newtonUpdate(f, x);
        ++iter;
    }
    if (iter == max_iter) {
        throw std::runtime_error("Too many Newton iterations");
    }
    return { x, iter};
}




int main()
{
    // Part 1
    // {
    //     // Define function.
    //     auto f = [](auto x) { return sin(2*x) + x + 1; };

    //     // Solve for a root.
    //     std::cout.precision(17);
    //     double start = 0.0;
    //     auto root = newton(f, start);
    //     std::cout << "Root = " << root.first << "   Iterations = " << root.second << std::endl;

    //     if (root.first == -0.3522884564608728) {
    //         std::cout << "Part 1 success!" << std::endl;
    //     } else {
    //         std::cout << "Part 1 failure!" << std::endl;
    //     }
    // }

    // Uncomment part 2 to do the exercise!

    // Part 2
    {
        // Define function.
        auto f = [](auto x) { return x/(sin(2*x) + 1.5)*cos(x) + 0.5; };

        // Solve for a root.
        std::cout.precision(16);
        double start = 0.0;
        auto root = newton(f, start);
        std::cout << "Root = " << root.first << "   Iterations = " << root.second << std::endl;

        if (root.first == -1.214838293477704) {
            std::cout << "Part 2 success!" << std::endl;
        } else {
            std::cout << "Part 2 failure!" << std::endl;
        }
    }

    // Note that the root found in part 2 is not the one closest to 0!
    // That is a "problem" with Newton's method and has nothing to do with AD though.
}
