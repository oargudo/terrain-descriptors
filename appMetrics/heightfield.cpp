#include "heightfield.h"

const double HeightField::flat = 1.0e-8;

/*!
\brief Create a flat heightfield.
\param box Rectangle domain of the terrain.
\param nx, ny Samples.
\param v Constant elevation.
*/
HeightField::HeightField(const Box2& box, int nx, int ny, const double& v)
    : ScalarField2(box, nx, ny, v)
{
}

/*!
\brief Create a heightfield.
\param box Rectangle domain of the terrain.
\param nx, ny Samples.
\param v Set of elevation values.
*/
HeightField::HeightField(const Box2& box, int nx, int ny, const std::vector<double>& v)
    : ScalarField2(box, nx, ny, v)
{
}

/*!
\brief Create a heightfield from a scalar field.
This constructor provides implicit conversion.
\param s Scalar field.
*/
HeightField::HeightField(const ScalarField2& s)
    : ScalarField2(s)
{
}

/*!
\brief Create a heightfield from an image.
\param box Rectangle domain of the terrain.
\param image Elevation image.
\param a, b Minimum and maximum elevation range.
\param grayscale Boolean set to false if the image is provided in color.
*/
HeightField::HeightField(const Box2& box, const QImage& image, const double& a, const double& b, bool grayscale)
    : ScalarField2(box, image, a, b, grayscale)
{
}

/*!
\brief Compute the elevation of a given position on the terrain.
\param p Point.
\param triangular Boolean, use triangular interpolation if set to true, bilinear interpolation otherwise.
*/
inline double HeightField::Height(const Vector2& p, bool triangular) const
{
    double u, v;
    int i, j;
    cellCoords(p, i, j, u, v);

    double z = 0.0;
    if (isValidCell(i, j)) {
        if (triangular) {
            if (u > v) {
                z = (1.0 - u) * at(i, j) + (u - v) * at(i + 1, j) + v * at(i + 1, j + 1);
            }
            else {
                z = (1.0 - v) * at(i, j) + u * at(i + 1, j + 1) + (v - u) * at(i, j + 1);
            }
        }
        else {
            z = Math::Bilinear(at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1), u, v);
        }
    }
    return z;
}


/*!
\brief Compute the vertex position on the terrain.
\param p Point.
\param triangular Boolean, use triangular interpolation if set to true, bilinear interpolation otherwise.
*/
Vector3 HeightField::Vertex(const Vector2& p, bool triangular) const
{
    return Vector3(p[0], p[1], Height(p, triangular));
}


/*!
\brief Compute the gradient at a given array vertex.

\param i,j Integer coordinates of the array vertex.
*/
Vector2 HeightField::Gradient(int i, int j) const
{
    Vector2 n;

    // Gradient along x axis
    if (i == 0)
        n[0] = (at(i + 1, j) - at(i, j)) * inverseCellSize[0];
    else if (i == nx - 1)
        n[0] = (at(i, j) - at(i - 1, j)) * inverseCellSize[0];
    else
        n[0] = (at(i + 1, j) - at(i - 1, j)) * 0.5 * inverseCellSize[0];

    // Gradient along y axis
    if (j == 0)
        n[1] = (at(i, j + 1) - at(i, j)) * inverseCellSize[1];
    else if (j == ny - 1)
        n[1] = (at(i, j) - at(i, j - 1)) * inverseCellSize[1];
    else
        n[1] = (at(i, j + 1) - at(i, j - 1)) * 0.5 * inverseCellSize[1];

    return n;
}

Vector2 HeightField::Gradient(const Vector2 &p) const
{
    int i, j;
    double u, v;

    cellCoords(p, i, j, u, v);

    if (!isValidCell(i, j)) return Vector2(0,0);

    Vector2 g00 = Gradient(i, j);
    Vector2 g10, g01, g11;
    if (isValidCell(i+1, j))   g10 = Gradient(i+1, j);
    else                             g10 = g00;
    if (isValidCell(i, j+1))   g01 = Gradient(i, j+1);
    else                             g01 = g00;
    if (isValidCell(i+1, j+1)) g11 = Gradient(i+1, j+1);
    else                             g11 = g00;

    Vector2 g0 = (1-v)*g00 + v*g01;
    Vector2 g1 = (1-v)*g10 + v*g11;
    return (1-u)*g0 + u*g1;
}


/*!
\brief Compute the normal for a given position on the terrain.

Note that this function may be expensive to compute.

\param p Point.
\param triangular Boolean, use triangle normals if set to true, bilinear interpolation of normals at vertices otherwise.
*/
Vector3 HeightField::Normal(const Vector2& p, bool triangular) const
{
    double u, v;
    int i, j;
    cellCoords(p, i, j, u, v);

    // Test position
    if (!isValidCell(i, j))
        return Vector3(0,0,0);

    if (triangular) {
        if (u > v) {
            return Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).Normal();
        }
        else {
            return Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).Normal();
        }
    }
    else {
        return Normalized(Bilinear(Normal(i, j), Normal(i + 1, j), Normal(i + 1, j + 1), Normal(i, j + 1), u, v));
    }
}


/*!
\brief Compute the normal at a given sample.

This function uses the weighted sum (area) of the normals of the
triangles sharing the point on the grid. The returned vector is normalized.

\param i,j Integer coordinates of the sample.
*/
Vector3 HeightField::Normal(int i, int j) const
{
  Vector3 n;
  if (i == 0)
  {
    if (j == 0)
    {
      // Corner: 0/1
      n = Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).AreaNormal() + Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).AreaNormal();
    }
    else if (j == ny - 1)
    {
      // Corner: 5
      n = Triangle(Vertex(i, j), Vertex(i, j - 1), Vertex(i + 1, j)).AreaNormal();
    }
    else
    {
      // Edge: 0/1/5
      n = Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).AreaNormal() + Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).AreaNormal()
        + Triangle(Vertex(i, j), Vertex(i, j - 1), Vertex(i + 1, j)).AreaNormal();
    }
  }
  else if (i == nx - 1)
  {
    if (j == 0)
    {
      // Corner: 2
      n = Triangle(Vertex(i, j), Vertex(i, j + 1), Vertex(i - 1, j)).AreaNormal();

    }
    else if (j == ny - 1)
    {
      // Corner: 3/4
      n = Triangle(Vertex(i, j), Vertex(i - 1, j - 1), Vertex(i, j - 1)).AreaNormal() + Triangle(Vertex(i, j), Vertex(i - 1, j), Vertex(i - 1, j - 1)).AreaNormal();
    }
    else
    {
      // Edge: 2/3/4
      n = Triangle(Vertex(i, j), Vertex(i, j + 1), Vertex(i - 1, j)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j), Vertex(i - 1, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j - 1), Vertex(i, j - 1)).AreaNormal();
    }
  }
  else
  {
    if (j == 0)
    {
      // Edge: 0/1/2
      n = Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i, j + 1), Vertex(i - 1, j)).AreaNormal();
    }
    else if (j == ny - 1)
    {
      // Edge: 3/4/5
      n = Triangle(Vertex(i, j), Vertex(i - 1, j), Vertex(i - 1, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j - 1), Vertex(i, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i, j - 1), Vertex(i + 1, j)).AreaNormal();
    }
    else
    {
      // Face: 0/1/2/3/4/5
      n = Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i, j + 1), Vertex(i - 1, j)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j), Vertex(i - 1, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j - 1), Vertex(i, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i, j - 1), Vertex(i + 1, j)).AreaNormal();
    }
  }
  return Normalized(n);
}
