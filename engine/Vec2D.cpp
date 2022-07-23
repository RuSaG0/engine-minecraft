//
// Created by Иван Ильин on 10.10.2021.
//

#include <cmath>
#include <cassert>

#include "Vec2D.h"
#include "Consts.h"

Vec2D::Vec2D(const Vec2D &vec) {
    _arr_point[0] = vec.x();
    _arr_point[1] = vec.y();
}

Vec2D::Vec2D(double x, double y) {
    _arr_point[0] = x;
    _arr_point[1] = y;
}

Vec2D::Vec2D(const Vec4D &point4D) {
    _arr_point[0] = point4D.x();
    _arr_point[1] = point4D.y();
}

Vec2D Vec2D::operator-() const {
    return Vec2D(-x(), -y());
}

bool Vec2D::operator==(const Vec2D &vec) const {
    Vec2D diff = *this - vec;
    return diff.sqrAbs() < Consts::EPS;
}

bool Vec2D::operator!=(const Vec2D &vec) const {
    return !(*this == vec);
}

Vec2D Vec2D::operator+(const Vec2D &vec) const {
    return Vec2D(x() + vec.x(), y() + vec.y());
}

Vec2D Vec2D::operator-(const Vec2D &vec) const {
    return Vec2D(x() - vec.x(), y() - vec.y());
}

Vec2D Vec2D::operator*(double number) const {
    return Vec2D(x() * number, y() * number);
}

Vec2D Vec2D::operator/(double number) const {

    if(std::abs(number) > Consts::EPS) {
        return Vec2D(x() / number, y() / number);
    } else {
        throw std::domain_error{"Vec2D::operator/(double number): division by zero"};
    }
}

// Other useful methods
double Vec2D::sqrAbs() const {
    return x() * x() + y() * y();
}

double Vec2D::abs() const {
    return std::sqrt(sqrAbs());
}

Vec2D Vec2D::normalized() const {
    double vecAbs = abs();

    if(vecAbs > Consts::EPS){
        return Vec2D(*this)/abs();
    } else {
        return Vec2D(0);
    }
}

double Vec2D::dot(const Vec2D &vec) const {
    return x() * vec.x() + y() * vec.y();
}

bool Vec2D::isNear(double a, double b) {
    return std::abs(a - b) < Consts::EPS;
}

void Vec2D::test() {
    Vec2D a(1, 2);
    Vec2D b(3, 4);

    // testing copy constructor:
    Vec2D c(a);
    assert(isNear(c.x(), 1) && isNear(c.y(), 2));

    // testing assigment operator (=):
    c = b;
    assert(isNear(c.x(), 3) && isNear(c.y(), 4));

    // testing .x() & .y() methods:
    assert(isNear(a.x(), 1) && isNear(a.y(), 2));
    assert(isNear(b.x(), 3) && isNear(b.y(), 4));

    // testing operator -Vec:
    Vec2D neg = -a;
    assert(isNear(neg.x(), -a.x()) && isNear(neg.y(), -a.y()));

    // testing == & != operators:
    assert(c != a && c == b);

    // testing operators +, -:
    Vec2D summ = a + b;
    assert(isNear(summ.x(), 4) && isNear(summ.y(), 6));
    Vec2D diff = a - b;
    assert(isNear(diff.x(), -2) && isNear(diff.y(), -2));

    // testing scaling operators:
    Vec2D scale1 = a * 2;
    assert(isNear(scale1.x(), 2) && isNear(scale1.y(), 4));
    Vec2D scale2 = a / 2;
    assert(isNear(scale2.x(), 0.5) && isNear(scale2.y(), 1));

    // testing dot product:
    double dot = a.dot(b);
    assert(isNear(dot, 11));

    // testing .abs() & .normalized() methods:
    assert(isNear(b.abs(), 5));
    assert(isNear(b.normalized().abs(), 1));
}