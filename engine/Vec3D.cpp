//
// Created by Иван Ильин on 09.10.2021.
//

#include <cmath>
#include <stdexcept>
#include <cassert>

#include "Vec3D.h"
#include "Consts.h"

Vec3D::Vec3D(const Vec3D &vec) {
    _arr_point[0] = vec.x();
    _arr_point[1] = vec.y();
    _arr_point[2] = vec.z();
}

Vec3D::Vec3D(const Vec4D &point4D) {
    _arr_point[0] = point4D.x();
    _arr_point[1] = point4D.y();
    _arr_point[2] = point4D.z();
}

Vec3D::Vec3D(double x, double y, double z) {
    _arr_point[0] = x;
    _arr_point[1] = y;
    _arr_point[2] = z;
}

Vec3D Vec3D::operator-() const {
    return Vec3D(-x(), -y(), -z());
}

bool Vec3D::operator==(const Vec3D &vec) const {
    return (*this - vec).sqrAbs() < Consts::EPS;
}

bool Vec3D::operator!=(const Vec3D &vec) const {
    return !(*this == vec);
}

// Operations with Vec3D
Vec3D Vec3D::operator+(const Vec3D &vec) const {
    return Vec3D(x() + vec.x(), y() + vec.y(), z() + vec.z());
}

Vec3D Vec3D::operator-(const Vec3D &vec) const {
    return *this + -vec;
}

Vec3D Vec3D::operator*(double number) const {
    return Vec3D(x() * number, y() * number, z() * number);
}

Vec3D Vec3D::operator/(double number) const {
    if (std::abs(number) > Consts::EPS) {
        return *this * (1.0 / number);
    } else {
        throw std::domain_error{"Vec3D::operator/(double number): division by zero"};
    }
}

// Other useful methods
double Vec3D::sqrAbs() const {
    return x() * x() + y() * y() + z() * z();
}

double Vec3D::abs() const {
    return sqrt(sqrAbs());
}

Vec3D Vec3D::normalized() const {
    double vecAbs = sqrAbs();
    if (vecAbs > Consts::EPS) {
        return *this / sqrt(vecAbs);
    } else {
        return Vec3D(1);
    }
}

double Vec3D::dot(const Vec3D &vec) const {
    return vec.x() * x() + vec.y() * y() + vec.z() * z();
}

Vec3D Vec3D::cross(const Vec3D &vec) const {
    return Vec3D(y() * vec.z() - vec.y() * z(),
                 z() * vec.x() - vec.z() * x(),
                 x() * vec.y() - vec.x() * y());
}

Vec4D Vec3D::makePoint4D() const {
    return Vec4D(x(), y(), z(), 1.0);
}

Vec3D Vec3D::Random() {
    return Vec3D((double) rand() / RAND_MAX, (double) rand() / RAND_MAX, (double) rand() / RAND_MAX);
}

bool Vec3D::isNear(double a, double b) {
    return std::abs(a - b) < Consts::EPS;
}

void Vec3D::test() {
    Vec3D a(1, 2, 3);
    Vec3D b(3, 4, 5);

    // testing copy constructor:
    Vec3D c(a);
    assert(isNear(c.x(), 1) && isNear(c.y(), 2) && isNear(c.z(), 3));

    // testing assigment operator (=):
    c = b;
    assert(isNear(c.x(), 3) && isNear(c.y(), 4) && isNear(c.z(), 5));

    // testing .x() & .y() & .z() methods:
    assert(isNear(a.x(), 1) && isNear(a.y(), 2) && isNear(a.z(), 3));
    assert(isNear(b.x(), 3) && isNear(b.y(), 4) && isNear(b.z(), 5));

    // testing operator -Vec:
    Vec3D neg = -a;
    assert(isNear(neg.x(), -a.x()) && isNear(neg.y(), -a.y()) && isNear(neg.z(), -a.z()));

    // testing == & != operators:
    assert(c != a && c == b);

    // testing operators +, -:
    Vec3D summ = a + b;
    assert(isNear(summ.x(), 4) && isNear(summ.y(), 6) && isNear(summ.z(), 8));
    Vec3D diff = a - b;
    assert(isNear(diff.x(), -2) && isNear(diff.y(), -2)  && isNear(diff.z(), -2));

    // testing scaling operators:
    Vec3D scale1 = a * 2;
    assert(isNear(scale1.x(), 2) && isNear(scale1.y(), 4) && isNear(scale1.z(), 6));
    Vec3D scale2 = a / 2;
    assert(isNear(scale2.x(), 0.5) && isNear(scale2.y(), 1) && isNear(scale2.z(), 1.5));

    // testing dot product:
    double dot = a.dot(b);
    assert(isNear(dot, 26));

    // testing cross product:
    Vec3D cross = a.cross(b);
    Vec3D cross_neg = b.cross(a);
    assert(isNear(a.dot(cross), 0) && isNear(b.dot(cross), 0));
    assert(cross == -cross_neg);

    // testing .abs() & .normalized() methods:
    assert(isNear(b.abs(), std::sqrt(50)));
    assert(isNear(b.normalized().abs(), 1));
}
