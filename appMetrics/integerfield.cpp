#include "integerfield.h"
#include "scalarfield.h"
#include <algorithm>


IntField2::IntField2() : FieldGrid2D(Box2(), 0, 0)
{
}

IntField2::IntField2(const Box2 &domain, int nx, int ny, int v) : FieldGrid2D(domain, nx, ny)
{
    field = std::vector<int>(nx*ny, v);
}

IntField2::IntField2(const Box2& domain, int nx, int ny, const std::vector<int>& v) : FieldGrid2D(domain, nx, ny)
{
    field = v;
}


/*!
\brief Get the range of the field.
\param a,b Returned minimum and maximum.
*/
void IntField2::getRange(int& a, int& b) const
{
    a = field.at(0);
    b = a;
    for (unsigned int i = 1; i < field.size(); i++) {
        double x = field.at(i);
        if (x < a) {
            a = x;
        }
        else if (x > b) {
            b = x;
        }
    }
}

int IntField2::percentile(double p) const
{
    std::vector<size_t> indices(field.size());
    std::iota(indices.begin(), indices.end(), 0);

    size_t n = std::min(size_t(p * field.size()), field.size()-1);
    std::nth_element(indices.begin(), indices.begin() + n, indices.end(),
        [this](size_t i1, size_t i2) {
            return field[i1] < field[i2];
        });
    return field[indices[n]];
}

/*!
\brief Sets the entire field with a constant value.
\param s Scalar.
*/
void IntField2::fill(int s)
{
    std::fill(field.begin(), field.end(), s);
}

QImage IntField2::createImage(const LookupPalette& palette) const
{
    QImage image(nx, ny, QImage::Format_RGBA8888);
    for (int i = 0; i < image.width(); i++) {
        for (int j = 0; j < image.height(); j++) {
            QColor color = toQColor(palette.getColor(at(i, j)));
            image.setPixel(i, j, color.rgb());
        }
    }
    return image;
}

ScalarField2 IntField2::toScalarField() const
{
    ScalarField2 s(domain, nx, ny);
    for (unsigned int i = 0; i < field.size(); i++) {
        s[i] = field[i];
    }
    return s;
}
