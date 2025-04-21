#ifndef SCALARFIELD_H
#define SCALARFIELD_H

#include "core.h"
#include <QImage>


class ScalarField2 : public FieldGrid2D
{
public:
    ScalarField2();
    ScalarField2(const Box2& domain, int nx, int ny, double v = 0.0);
    ScalarField2(const Box2& domain, int nx, int ny, const std::vector<double>& v);
    ScalarField2(const Box2& domain, const QImage& image, double a, double b, bool grayscale);
    virtual ~ScalarField2() {};

    double& operator()(int i, int j) { return field[cellId(i, j)]; }
    double& operator()(const Index2& p) { return field[cellId(p.x(), p.y())]; }
    double at(int i, int j) const { return field.at(cellId(i, j)); }
    double at(int i) const { return field.at(i); }
    double at(const Index2& p) const { return field.at(cellId(p.x(), p.y())); }
    double& operator[](int i) { return field[i]; }
    double value(const Vector2& p) const;

    void getRange(double& vmin, double& vmax) const;
    double percentile(double p) const;
    std::vector<double> percentiles(const std::vector<double>& p) const;
    double sum() const;
    double average() const;
    double standardDeviation(double avg) const;

    void fill(double d);
    void smooth(int n);
    void step(double a, double b);
    void threshold(double v);
    void normalize();
    void addGaussian(const Vector2& center, const double& radius, const double& height);

    ScalarField2 setResolution(int nx, int ny) const;

    ScalarField2 gaussianBlur(int r) const;
    ScalarField2 maxFilter(int w) const;
    ScalarField2 minFilter(int w) const;
    ScalarField2 summedAreaTable() const;

    QImage createImage(bool grayscale) const;
    QImage createImage(double a, double b, bool grayscale) const;
    QImage createImage(const ColorPalette& palette) const;
    QImage createImage(double a, double b, const ColorPalette& palette) const;

    ScalarField2& operator+=(const ScalarField2& s);
    ScalarField2& operator+=(const double& d);
    ScalarField2& operator*=(const double& d);
    friend ScalarField2 operator+(const ScalarField2& s1, const ScalarField2& s2);
    friend ScalarField2 operator*(const ScalarField2& s1, const ScalarField2& s2);

    const std::vector<double>& values() const { return field; }
    std::vector<std::pair<double, Index2>> valuesWithIndex() const;

protected:
    std::vector<double> field;

};

#endif // SCALARFIELD_H
