#include "scalarfield.h"
#include <algorithm>

ScalarField2::ScalarField2() : FieldGrid2D(Box2(Vector2(0), Vector2(1)), 0, 0)
{
}

ScalarField2::ScalarField2(const Box2 &domain, int nx, int ny, double v) : FieldGrid2D(domain, nx, ny)
{
    field = std::vector<double>(nx*ny, v);
}

ScalarField2::ScalarField2(const Box2& domain, int nx, int ny, const std::vector<double>& v) : FieldGrid2D(domain, nx, ny)
{
    field = v;
}

ScalarField2::ScalarField2(const Box2& domain, const QImage& image, double a, double b, bool grayscale)
    : FieldGrid2D(domain, image.width(), image.height())
{
    field.resize(nx*ny);

    // Write Heightmap
    for (int i = 0; i < image.width(); i++) {
        for (int j = 0; j < image.height(); j++) {
            double t = 0.0;

            // Grayscale
            if (grayscale) {
                // Grayscale 16 bits
                if (image.format() == QImage::Format_Grayscale16) {
                    QColor thecolor = image.pixelColor(i, j);
                    t = thecolor.blueF();
                }
                // Grayscale 8 bits
                else {
                    QRgb color = image.pixel(i, j);
                    t = double(qGray(color)) / 255.0;
                }
            }
            // Color
            else {
                QRgb color = image.pixel(i, j);
                // Maximum value is 256^3-1
                t = double(qRed(color) << 16 | qGreen(color) << 8 | qBlue(color)) / (16777216.0 - 1.0);
            }

            field[cellId(i, j)] = Math::Lerp(a, b, t);
        }
    }
}

double ScalarField2::value(const Vector2& p) const
{
    double u, v;
    int i = -1, j = -1;
    cellCoords(p, i, j, u, v);

    if (!isValidCell(i, j)) return 0.0;
    int ip1 = std::min(i+1, nx-1);
    int jp1 = std::min(j+1, ny-1);
    return Math::Bilinear(at(i, j), at(ip1, j), at(ip1, jp1), at(i, jp1), u, v);
}

void ScalarField2::getRange(double &vmin, double &vmax) const
{
    vmin = vmax = field.at(0);
    for (unsigned int i = 1; i < field.size(); i++) {
        double x = field.at(i);
        if (x < vmin) vmin = x;
        if (x > vmax) vmax = x;
    }
}

double ScalarField2::percentile(double p) const
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

std::vector<double> ScalarField2::percentiles(const std::vector<double>& percs) const
{
    std::vector<size_t> indices(field.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::vector<double> values;
    for (double p : percs) {
        size_t n = std::min(size_t(p * field.size()), field.size()-1);
        std::nth_element(indices.begin(), indices.begin() + n, indices.end(),
            [this](size_t i1, size_t i2) {
                return field[i1] < field[i2];
            });
        values.push_back(field[indices[n]]);
    }
    return values;
}

double ScalarField2::sum() const
{
    double s = 0.0;
    int size = nx*ny;
    for (int i = 0; i < size; i++) {
        s += field.at(i);
    }
    return s;
}

double ScalarField2::average() const
{
    const int size = nx * ny;
    double s = sum();
    s /= size;
    return s;
}

double ScalarField2::standardDeviation(double m) const
{
    int size = nx * ny;
    double s = 0.0;
    for (int i = 0; i < size; i++) {
        double x = field.at(i) - m;
        s += x * x;
    }
    return sqrt(s / size);
}

void ScalarField2::fill(double d)
{
    std::fill(field.begin(), field.end(), d);
}

void ScalarField2::smooth(int n)
{
    // Smooth the scalar field using a discrete gaussian kernel.
    // The function uses a 3^2 approximation of the Gaussian kernel.
    std::vector<double> smoothed;
    smoothed.resize(nx * ny);
    int k;

    for (int iter = 0; iter < n; iter++) {
        // Smooth center
        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                k = cellId(i, j);
                smoothed[k] = (4.0*at(k) + 2.0*at(k - 1) + 2.0*at(k + 1) + 2.0*at(k - nx) + 2.0*at(k + nx) + at(k - 1 - nx) + at(k + 1 - nx) + at(k - 1 + nx) + at(k + 1 + nx)) / 16.0;
            }
        }

        // Smooth edges
        for (int i = 1; i < nx - 1; i++) {
            k = cellId(i, 0);
            smoothed[k] = (4.0*at(k) + 2.0*at(k - 1) + 2.0*at(k + 1) + 2.0*at(k + nx) + at(k - 1 + nx) + at(k + 1 + nx)) / 12.0;
            k = cellId(i, ny - 1);
            smoothed[k] = (4.0*at(k) + 2.0*at(k - 1) + 2.0*at(k + 1) + 2.0*at(k - nx) + at(k - 1 - nx) + at(k + 1 - nx)) / 12.0;
        }
        for (int j = 1; j < ny - 1; j++) {
            k = cellId(0, j);
            smoothed[k] = (2.0*at(k - nx) + 4.0*at(k) + 2.0*at(k + nx) + at(k + 1 - nx) + 2.0*at(k + 1) + at(k + 1 + nx)) / 12.0;
            k = cellId(nx - 1, j);
            smoothed[k] = (2.0*at(k - nx) + 4.0*at(k) + 2.0*at(k + nx) + at(k - 1 - nx) + 2.0*at(k - 1) + at(k - 1 + nx)) / 12.0;
        }

        // Corners
        k = cellId(0, 0);
        smoothed[k] = (4.0*at(k) + 2.0*at(k + 1) + 2.0*at(k + nx) + 1.0*at(k + nx + 1)) / 9.0;
        k = cellId(nx - 1, 0);
        smoothed[k] = (4.0*at(k) + 2.0*at(k - 1) + 2.0*at(k + nx) + 1.0*at(k + nx - 1)) / 9.0;
        k = cellId(0, ny - 1);
        smoothed[k] = (4.0*at(k) + 2.0*at(k + 1) + 2.0*at(k - nx) + 1.0*at(k - nx + 1)) / 9.0;
        k = cellId(nx - 1, ny - 1);
        smoothed[k] = (4.0*at(k) + 2.0*at(k - 1) + 2.0*at(k - nx) + 1.0*at(k - nx - 1)) / 9.0;
    }
    field = smoothed;
}

/*!
\brief Perform a Gaussian blur.
\param r Radius (in cells) of blur.
*/
ScalarField2 ScalarField2::gaussianBlur(int r) const
{
  // Fixed sigma to half kernel radius
  double sigma = 0.5 * double(r);
  double sigma2 = sigma * sigma;
  double norm = 1.0 / (sqrt(2.0 * Math::Pi) * sigma);

  // Kernel
  std::vector<double> kernel(r + 1);
  for (int i = 0; i <= r; i++)
  {
    kernel[i] = norm * exp((-0.5 * double(i) * double(i)) / sigma2);
  }

  std::vector<double> resx(nx * ny);

  // Filter x direction
  for (int y = 0; y < ny; y++)
  {
    for (int x = 0; x < nx; x++)
    {
      double v = kernel[0] * at(x, y);
      double w = kernel[0];
      for (int k = 1; k <= r; k++)
      {
        if (isValidCell(x + k, y))
        {
          v += kernel[k] * at(x + k, y);
          w += kernel[k];
        }
        if (isValidCell(x - k, y))
        {
          v += kernel[k] * at(x - k, y);
          w += kernel[k];
        }
      }
      resx[cellId(x, y)] = v / w;
    }
  }

  ScalarField2 res(domain, nx, ny);

  // Filter y direction
  for (int x = 0; x < nx; x++)
  {
    for (int y = 0; y < ny; y++)
    {
      double v = kernel[0] * resx[cellId(x, y)];
      double w = kernel[0];
      for (int k = 1; k <= r; k++)
      {
        if (isValidCell(x, y + k))
        {
          v += kernel[k] * resx[cellId(x, y + k)];
          w += kernel[k];
        }
        if (isValidCell(x, y - k))
        {
          v += kernel[k] * resx[cellId(x, y - k)];
          w += kernel[k];
        }
      }
      res(x, y) = v / w;
    }
  }

  return res;
}

ScalarField2 ScalarField2::maxFilter(int w) const
{
  ScalarField2 filtered(getDomain(), nx, ny, 0);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double v = at(i, j);
      for (int di = -w; di <= w; di++)
      {
        for (int dj = -w; dj <= w; dj++)
        {
          if (isValidCell(i + di, j + dj))
          {
            v = std::max(v, at(i + di, j + dj));
          }
        }
      }
      filtered(i, j) = v;
    }
  }

  return filtered;
}


ScalarField2 ScalarField2::minFilter(int w) const
{
  ScalarField2 filtered(getDomain(), nx, ny, 0);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double v = at(i, j);
      for (int di = -w; di <= w; di++)
      {
        for (int dj = -w; dj <= w; dj++)
        {
          if (isValidCell(i + di, j + dj))
          {
            v = std::min(v, at(i + di, j + dj));
          }
        }
      }
      filtered(i, j) = v;
    }
  }

  return filtered;
}


ScalarField2 ScalarField2::summedAreaTable() const
{
    ScalarField2 sat(domain, nx, ny, 0);
    sat[0] = at(0, 0);
    for (int i = 1; i < nx; i++)
      sat[cellId(i, 0)] = at(i, 0) + sat[cellId(i - 1, 0)];
    for (int j = 1; j < ny; j++)
      sat[cellId(0, j)] = at(0, j) + sat[cellId(0, j - 1)];
    for (int i = 1; i < nx; i++)
    {
      for (int j = 1; j < ny; j++)
      {
        sat[cellId(i, j)] = at(i, j)
          + sat[cellId(i - 1, j)]
          + sat[cellId(i, j - 1)]
          - sat[cellId(i - 1, j - 1)];
      }
    }
    return sat;
}


void ScalarField2::step(double a, double b)
{
    for (unsigned int i = 0; i < field.size(); i++) {
        field[i] = Math::LinearStep(field.at(i), a, b);
    }
}

void ScalarField2::threshold(double v)
{
    for (unsigned int i = 0; i < field.size(); i++) {
        field[i] = field.at(i) >= v ? 1 : 0;
    }
}

void ScalarField2::normalize()
{
    double a, b;
    getRange(a, b);
    if (a == b) {
        std::fill(field.begin(), field.end(), 1);
    }
    else {
        for (unsigned int i = 0; i < field.size(); i++) {
            field[i] = (field[i] - a) / (b - a);
        }
    }
}

void ScalarField2::addGaussian(const Vector2& center, const double& radius, const double& height)
{

    Box2 box (center, radius);

    int ia,ib,ja,jb;
    cellCoords(box.getMin(), ia, ja);
    cellCoords(box.getMax(), ib, jb);

    QPoint pa(ia, ja);
    QPoint pb(ib, jb);

    // Rectangle
    QRect area(pa.x(), pa.y(), pb.x() - pa.x() + 1, pb.y() - pa.y() + 1);

    // Limit to domain
    QRect mask(0, 0, nx - 1, ny - 1);
    area = area.intersected(mask);

    // Add to field
    for (int y = area.y(); y <= area.y() + area.height(); y++) {
        for (int x = area.x(); x <= area.x() + area.width(); x++) {
            // Distance between central point and current point
            double u = SquaredNorm(center - domainCoords(x, y));

            if (u < radius * radius) {
                field[cellId(x, y)] += height * Math::CubicSmooth(u, radius * radius);
           }
        }
    }
}

ScalarField2 ScalarField2::setResolution(int x, int y) const
{
    // Sampled scalar field
    ScalarField2 sampled(domain, x, y);

    // Corners
    sampled(0, 0) = at(0, 0);
    sampled(0, y - 1) = at(0, ny - 1);
    sampled(x - 1, 0) = at(nx - 1, 0);
    sampled(x - 1, y - 1) = at(nx - 1, ny - 1);

    // Borders (use linear interpolation)
    for (int i = 1; i < x - 1; i++) {
        double tx = (nx - 1) * (i / double(x - 1));
        int x0 = int(floor(tx));
        int x1 = int(ceil(tx));

        sampled(i, 0) = Math::Lerp(at(x0, 0), at(x1, 0), tx - x0);
        sampled(i, y - 1) = Math::Lerp(at(x0, ny - 1), at(x1, ny - 1), tx - x0);
    }
    for (int j = 1; j < y - 1; j++) {
        double ty = (ny - 1) * (j / double(y - 1));
        int y0 = int(floor(ty));
        int y1 = int(ceil(ty));

        sampled(0, j) = Math::Lerp(at(0, y0), at(0, y1), ty - y0);
        sampled(x - 1, j) = Math::Lerp(at(nx - 1, y0), at(nx - 1, y1), ty - y0);
    }

    // Interior
    for (int i = 1; i < x - 1; i++) {
        for (int j = 1; j < y - 1; j++) {
            sampled(i, j) = value(sampled.domainCoords(i, j));
        }
    }

    return sampled;
}

QImage ScalarField2::createImage(bool grayscale) const
{
    double a, b;
    this->getRange(a, b);
    if (a == b) {
        b = a + 1.0;
    }
    return createImage(a, b, grayscale);
}

QImage ScalarField2::createImage(double a, double b, bool grayscale) const
{
    QImage image(nx, ny, QImage::Format_RGBA8888);
    for (int i = 0; i < image.width(); i++) {
        for (int j = 0; j < image.height(); j++) {
            double x = field.at(cellId(i, j));
            double y = Math::LinearStep(x, a, b);

            QColor color;
            if (grayscale) {
                int c = int(y * 255.0);
                color = QColor(c, c, c);
            }
            else {
                int c = int(y * (256.0 * 256.0 * 256.0 - 1.0));
                int cr = (c >> 16) & 255;
                int cg = (c >> 8) & 255;
                int cb = c & 255;
                color = QColor(cr, cg, cb);
            }
            image.setPixel(i, j, color.rgb());
        }
    }
    return image;
}

QImage ScalarField2::createImage(const ColorPalette& palette) const
{
    double a, b;
    getRange(a, b);
    if (a == b) {
        b = a + 1.0;
    }
    return createImage(a, b, palette);
}

QImage ScalarField2::createImage(double a, double b, const ColorPalette& palette) const
{
    QImage image(nx, ny, QImage::Format_RGBA8888);
    for (int i = 0; i < image.width(); i++) {
        for (int j = 0; j < image.height(); j++) {
            double x = field.at(cellId(i, j));
            double y = Math::LinearStep(x, a, b);
            QColor color = toQColor(palette.getColor(y));
            image.setPixel(i, j, color.rgb());
        }
    }
    return image;
}

ScalarField2& ScalarField2::operator+=(const ScalarField2& s)
{
    for (unsigned int i = 0; i < field.size(); i++) {
        field[i] += s.at(i);
    }
    return *this;
}

ScalarField2 &ScalarField2::operator+=(const double &d)
{
    for (unsigned int i = 0; i < field.size(); i++) {
        field[i] += d;
    }
    return *this;
}

ScalarField2& ScalarField2::operator*=(const double &d)
{
    for (unsigned int i = 0; i < field.size(); i++) {
        field[i] *= d;
    }
    return *this;
}

ScalarField2 operator+(const ScalarField2& s1, const ScalarField2& s2)
{
    ScalarField2 r(s1);
    for (unsigned int i = 0; i < r.field.size(); i++) {
        r.field[i] += s2.at(i);
    }
    return r;
}

ScalarField2 operator*(const ScalarField2& s1, const ScalarField2& s2)
{
    ScalarField2 r(s1);
    for (unsigned int i = 0; i < r.field.size(); i++) {
        r.field[i] *= s2.at(i);
    }
    return r;
}

/*!
\brief Return all the field data points as a list of value-cell.
*/
std::vector<std::pair<double, Index2>> ScalarField2::valuesWithIndex() const
{
    std::vector<std::pair<double, Index2>> pts(nx*ny);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            pts[cellId(i, j)] = std::make_pair(at(i, j), Index2(i, j));
        }
    }
    return pts;
}
