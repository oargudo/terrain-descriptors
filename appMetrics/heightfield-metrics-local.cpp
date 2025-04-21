#include "heightfield.h"



/*!
\brief Compute the directional slope between two points, i.e. difference in height divided by distance

\param pa,pb start and end cell indices.
*/
inline double HeightField::SlopeDirectional(const Index2& pa, const Index2& pb) const
{
    double h = fabs(at(pa) - at(pb));
    Index2 pab = pb - pa;
    double l = Norm(Vector2(pab.x() * cellSize[0], pab.y() * cellSize[1]));
    return h / l;
}


/*!
\brief Compute the average slope for every position on the terrain.

Simply average the slope in 8 directions.
*/
ScalarField2 HeightField::SlopeAverage() const
{
    ScalarField2 slope(domain, nx, ny);
    double e = Norm(cellSize);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {            
            double s;
            if (i == 0) {
                if (j == 0) {
                    // Corner
                    s = fabs(at(i, j) - at(i + 1, j)) * inverseCellSize[0] + fabs(at(i, j) - at(i + 1, j + 1)) / e + fabs(at(i, j) - at(i, j + 1)) * inverseCellSize[1];
                    s /= 3.0;
                }
                else if (j == ny - 1) {
                    // Corner
                    s = fabs(at(i, j) - at(i + 1, j)) * inverseCellSize[0] + fabs(at(i, j) - at(i + 1, j - 1)) / e + fabs(at(i, j) - at(i, j - 1)) * inverseCellSize[1];
                    s /= 3.0;
                }
                else {
                    // Edge
                    s = fabs(at(i, j) - at(i, j - 1)) * inverseCellSize[1] + fabs(at(i, j) - at(i + 1, j - 1)) / e + fabs(at(i, j) - at(i + 1, j)) * inverseCellSize[0] + fabs(at(i, j) - at(i + 1, j + 1)) / e + fabs(at(i, j) - at(i, j + 1)) * inverseCellSize[1];
                    s /= 5.0;
                }
            }
            else if (i == nx - 1) {
                if (j == 0) {
                    // Corner
                    s = fabs(at(i, j) - at(i - 1, j)) * inverseCellSize[0] + fabs(at(i, j) - at(i - 1, j + 1)) / e + fabs(at(i, j) - at(i, j + 1)) * inverseCellSize[1];
                    s /= 3.0;
                }
                else if (j == ny - 1) {
                    // Corner
                    s = fabs(at(i, j) - at(i - 1, j)) * inverseCellSize[0] + fabs(at(i, j) - at(i - 1, j - 1)) / e + fabs(at(i, j) - at(i, j - 1)) * inverseCellSize[1];
                    s /= 3.0;
                }
                else {
                    // Edge
                    s = fabs(at(i, j) - at(i, j - 1)) * inverseCellSize[1] + fabs(at(i, j) - at(i - 1, j - 1)) / e + fabs(at(i, j) - at(i - 1, j)) * inverseCellSize[0] + fabs(at(i, j) - at(i - 1, j + 1)) / e + fabs(at(i, j) - at(i, j + 1)) * inverseCellSize[1];
                    s /= 5.0;
                }
            }
            else {
                if (j == 0) {
                    // Edge
                    s = fabs(at(i, j) - at(i - 1, j)) * inverseCellSize[0] + fabs(at(i, j) - at(i - 1, j + 1)) / e + fabs(at(i, j) - at(i, j + 1)) * inverseCellSize[1] + fabs(at(i, j) - at(i + 1, j + 1)) / e + fabs(at(i, j) - at(i + 1, j)) * inverseCellSize[0];
                    s /= 5.0;
                }
                else if (j == ny - 1) {
                    // Edge
                    s = fabs(at(i, j) - at(i - 1, j)) * inverseCellSize[0] + fabs(at(i, j) - at(i - 1, j - 1)) / e + fabs(at(i, j) - at(i, j - 1)) * inverseCellSize[1] + fabs(at(i, j) - at(i + 1, j - 1)) / e + fabs(at(i, j) - at(i + 1, j)) * inverseCellSize[0];
                    s /= 5.0;
                }
                else {
                    // Vertex
                    s = fabs(at(i, j) - at(i + 1, j)) * inverseCellSize[0] + fabs(at(i, j) - at(i + 1, j + 1)) / e + fabs(at(i, j) - at(i, j + 1)) * inverseCellSize[1] + fabs(at(i, j) - at(i - 1, j + 1)) / e + fabs(at(i, j) - at(i - 1, j)) * inverseCellSize[0] + fabs(at(i, j) - at(i - 1, j - 1)) / e + fabs(at(i, j) - at(i, j - 1)) * inverseCellSize[1] + fabs(at(i, j) - at(i + 1, j - 1)) / e;
                    s /= 8.0;
                }
            }

            slope(i, j) = s;
        }
    }
    return slope;
}


/*!
\brief Compute the norm of the gradient for every position on the terrain.
*/
ScalarField2 HeightField::GradientNorm() const
{
    // Scalar field of the same size
    ScalarField2 n(domain, nx, ny);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            n(i, j) = Norm(Gradient(i, j));
        }
    }
    return n;
}


/*!
\brief Compute the slope aspect (orientation) for every position on the terrain.

Aspect angle is measured clockwise from North.
*/
ScalarField2 HeightField::Aspect() const
{
    ScalarField2 aspect(domain, nx, ny);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            Vector2 g = Gradient(i, j);
            aspect(i, j) = Math::Pi + std::atan2(g[0], g[1]);
        }
    }
    return aspect;
}


/*!
\brief Compute the laplacian for every position on the terrain.
*/
ScalarField2 HeightField::Laplacian() const
{
    ScalarField2 laplacian(domain, nx, ny);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {

            double lap = 0.0;

            // Divergence along x axis
            if (i == 0) {
                lap += (at(i, j) - 2.0 * at(i + 1, j) + at(i + 2, j)) / (cellSize[0] * cellSize[0]);
            }
            else if (i == nx - 1) {
                lap += (at(i, j) - 2.0 * at(i - 1, j) + at(i - 2, j)) / (cellSize[0] * cellSize[0]);
            }
            else {
                lap += (at(i + 1, j) - 2.0 * at(i, j) + at(i - 1, j)) / (cellSize[0] * cellSize[0]);
            }

            // Divergence along y axis
            if (j == 0) {
                lap += (at(i, j) - 2.0 * at(i, j + 1) + at(i, j + 2)) / (cellSize[1] * cellSize[1]);
            }
            else if (j == ny - 1) {
                lap += (at(i, j) - 2.0 * at(i, j - 1) + at(i, j - 2)) / (cellSize[1] * cellSize[1]);
            }
            else {
                lap += (at(i, j + 1) - 2.0 * at(i, j) + at(i, j - 1)) / (cellSize[1] * cellSize[1]);
            }

            laplacian(i, j) = lap;
        }
    }
    return laplacian;
}


/*!
\brief Compute the fractional laplacian for every position on the terrain.

\param s fractional exponent
\param n half-size (in cells) of integration domain
*/
ScalarField2 HeightField::FractLaplacian(double s, int n) const
{
    ScalarField2 flaplacian(domain, nx, ny);
    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {

            float sum = 0.f;

            // Factorize value at center
            float hxy = 2.f * at(x, y);

            // Area of integrand : not used as we often work on scaled Fractional Laplacian
            // double dxy=cellDiagonal.x*cellDiagonal.y;

            // Integration on half of the domain RxR-* allows to remove the test inside the loop
            for (int i = -n; i <= n; i++) {
                for (int j = -n; j < 0; j++) {
                    double d = pow(double(i*i + j*j), 1 + s);
                    double hp = at(Math::Clamp(x + i, 0, nx-1), std::max(y + j, 0));
                    double hn = at(Math::Clamp(x - i, 0, nx-1), std::min(y - j, ny - 1));
                    sum += (hp - hxy + hn) / d;
                }
            }
            // Complete integration on R- x {0}
            for (int i = -n; i < 0; i++) {
                double hp = at(Math::Clamp(x + i, 0, nx-1), y);
                double hn = at(Math::Clamp(x - i, 0, nx-1), y);
                sum += (hp - hxy + hn) / pow(double(i*i), 1 + s);
            }

            // The constant C is omitted
            flaplacian(x, y) = -sum;
        }
    }
    return flaplacian;
}


/*!
\brief Compute the quadric function that approximates the height field at a given location.

\param cellx,celly cell integer coordinates
\param w window size
\param interpolate if true, returns an interpolating quadric at (cellx, celly)
*/
QuadricSurface HeightField::FitQuadric(int cellx, int celly, int w, bool interpolate) const
{
    double M[6][6] = { { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 },
                       { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 } };
    double r[6] = { 0, 0, 0, 0, 0, 0 };
    int n = w / 2;
    double z0 = at(cellx, celly);

    for (int i = -n; i <= n; i++) {
      for (int j = -n; j <= n; j++) {
        double x = i * cellSize[0];
        double y = j * cellSize[1];
        double z = isValidCell(cellx + i, celly + j) ? at(cellx + i, celly + j) : 0;
        if (interpolate) z -= z0;

        double x2 = x * x;
        double y2 = y * y;

        M[0][0] += x2 * x2;
        M[1][1] += y2 * y2;
        M[0][1] += x2 * y2;
        M[1][0] += x2 * y2;
        M[2][2] += x2 * y2;
        M[3][3] += x2;
        M[4][4] += y2;
        M[5][0] += x2;
        M[5][1] += y2;
        M[0][5] += x2;
        M[0][1] += y2;
        M[5][5] += 1;

        r[0] += z * x2;
        r[1] += z * y2;
        r[2] += z * x * y;
        r[3] += z * x;
        r[4] += z * y;
        r[5] += z;
      }
    }

    double a, b, c, d, e, f;
    if (interpolate) {
      double det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
      a = (M[1][1] * r[0] - M[0][1] * r[1]) / det;
      b = (M[0][0] * r[1] - M[1][0] * r[0]) / det;
      c = r[2] / M[2][2];
      d = r[3] / M[3][3];
      e = r[4] / M[4][4];
      f = 0;
    }
    else {
      Matrix3 mat3(M[0][0], M[0][1], M[0][5], M[1][0], M[1][1], M[1][5], M[5][0], M[5][1], M[5][5]);
      Vector3 sol3 = Inverse(mat3) * Vector3(r[0], r[1], r[5]);
      a = sol3[0];
      b = sol3[1];
      c = r[2] / M[2][2];
      d = r[3] / M[3][3];
      e = r[4] / M[4][4];
      f = sol3[2];
    }

    return QuadricSurface(0, 0, 0, a, b, c, d, e, f);
}


/*!
\brief Compute a curvature measure on the heightfield.

\param curvType the type of curvature
\param w the size of the window. If 0, use the 3x3 fixed-formula curvatures. Otherwise, fits an interpolating quadric on the wxw window.
*/
ScalarField2 HeightField::Curvature(CurvatureType curvType, int w) const
{
    ScalarField2 curv(domain, nx, ny, 0.0);
    const double dx = cellSize[0];
    const double dy = cellSize[1];

    if (w == 0) {
        if (curvType == CurvatureType::MIN || curvType == CurvatureType::MAX) {
            ScalarField2 cGauss = Curvature(CurvatureType::GAUSSIAN);
            ScalarField2 cMean  = Curvature(CurvatureType::MEAN);
            for (int i = 1; i < nx-1; i++) {
                for (int j = 1; j < ny-1; j++) {
                    double k = cGauss.at(i, j);
                    double h = cMean.at(i, j);
                    switch (curvType) {
                    case MIN: curv(i, j) = h - sqrt(h * h - k); break;
                    case MAX: curv(i, j) = h + sqrt(h * h - k); break;
                    default: break;
                    }
                }
            }
        }
        else {
            for (int i = 1; i < nx - 1; i++) {
                for (int j = 1; j < ny - 1; j++) {
                    double z_x = (at(i + 1, j) - at(i - 1, j)) / (2.0 * dx);
                    double z_y = (at(i, j + 1) - at(i, j - 1)) / (2.0 * dy);
                    double z_xx = (at(i + 1, j) + at(i - 1, j) - 2.0 * at(i, j)) / (dx * dx);
                    double z_yy = (at(i, j + 1) + at(i, j - 1) - 2.0 * at(i, j)) / (dy * dy);
                    double z_xy = (at(i + 1, j + 1) + at(i - 1, j - 1) - at(i + 1, j - 1) - at(i - 1, j + 1)) / (4.0 * dx * dy);
                    double p = z_x * z_x + z_y * z_y;
                    double q = p + 1.0;

                    switch (curvType) {
                    case MEAN:      curv(i, j) = ((1 + z_y * z_y) * z_xx - 2 * z_xy * z_x * z_y + (1 + z_x * z_x) * z_yy) / (2.0 * Math::sqrt3_2(q)); break;
                    case GAUSSIAN:  curv(i, j) = (z_xx * z_yy - z_xy * z_xy) / (q * q); break;
                    case PROFILE:   curv(i, j) = p <= 1e-16 ? 0 : (z_xx * z_x * z_x + 2 * z_xy * z_x * z_y + z_yy * z_y * z_y) / (p * Math::sqrt3_2(q)); break;
                    case CONTOUR:   curv(i, j) = p <= 1e-16 ? 0 : (z_xx * z_y * z_y - 2.0 * z_xy * z_x * z_y + z_yy * z_x * z_x) / Math::sqrt3_2(p); break;
                    case TANGENTIAL:curv(i, j) = p <= 1e-16 ? 0 : (z_xx * z_y * z_y - 2.0 * z_xy * z_x * z_y + z_yy * z_x * z_x) / (p * sqrt(q)); break;
                    default:        curv(i, j) = 0; break;
                    }
                }
            }
        }
    }
    else {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                QuadricSurface q = FitQuadric(i, j, w, true);
                double a = q(2, 0);
                double b = q(0, 2);
                double c = q(1, 1);
                double d = q(1, 0);
                double e = q(0, 1);

                switch (curvType) {
                case MIN:       curv(i, j) = a + b - sqrt((a - b) * (a - b) + c * c); break;
                case MAX:       curv(i, j) = a + b + sqrt((a - b) * (a - b) + c * c); break;
                case MEAN:      curv(i, j) = a + b; break;
                case GAUSSIAN:  curv(i, j) = (a + b) * (a + b) - (a - b) * (a - b) + c * c; break;
                case PROFILE:   curv(i, j) = 2.0 * (a * d * d + b * e * e + c * e * d) / ((e * e + d * d) * Math::sqrt3_2(1.0 + d * d + e * e)); break;
                case CONTOUR:   curv(i, j) = 2.0 * (b * d * d + a * e * e - c * d * e) / Math::sqrt3_2(e * e + d * d); break;
                case TANGENTIAL:curv(i, j) = 2.0 * (b * d * d + a * e * e - c * d * e) / (d * d + e * e); break;
                default:          curv(i, j) = 0; break;
                }
            }
        }
    }
    return curv;
}



/*!
\brief Compute the local relief, the difference between maximum and minimum elevation in a neighborhood.

\param w Window size
*/
ScalarField2 HeightField::LocalRelief(int w) const
{
  ScalarField2 lr(domain, nx, ny, 0.0);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double hmin = at(i, j);
      double hmax = hmin;

      IndexArea area = indexAreaFromWindow(i, j, w);
      for (int wi = area.xmin(); wi <= area.xmax(); wi++)
      {
        for (int wj = area.ymin(); wj <= area.ymax(); wj++)
        {
          const double v = at(wi, wj);
          hmin = std::min(hmin, v);
          hmax = std::max(hmax, v);
        }
      }
      lr(i, j) = hmax - hmin;
    }
  }

  return lr;
}


/*!
\brief Compute the local variance

Woodcock and Strahler 1987. The factor of scale in remote sensing.
Dragut et al. 2011. Local variance for multi-scale analysis in geomorphometry.

Local variance is defined as the standard deviation of a window centered around the cell.

Compute the local variance of every cell and take the average over all the domain.

Authors plot the mean local variance with respect to image scale, and show that peaks in this graph
occur at 1/2 to 3/4 the size of objects in the image, so it could be used as a scale indicator.

\param w Window size (authors recommend 3).
*/
ScalarField2 HeightField::LocalVariance(int w) const
{
  ScalarField2 lv(domain, nx, ny, 0.0);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double s = 0;
      double ss = 0;
      int n = 0;

      IndexArea area = indexAreaFromWindow(i, j, w);
      for (int wi = area.xmin(); wi <= area.xmax(); wi++)
      {
        for (int wj = area.ymin(); wj <= area.ymax(); wj++)
        {
          const double v = at(wi, wj);
          s += v;
          ss += v * v;
          n++;
        }
      }

      // clamp to 0 to prevent numerical errors on flat areas that make sqrt negative
      lv(i, j) = sqrt(std::max(0.0, (n * ss - s * s)) / (double(n) * double(n - 1)));
    }
  }

  return lv;
}


/*!
\brief Compute the surface to planar ratio of areas, a measure of terrain ruggedness.

The function measures the area of the 8 triangles formed by a cell and its neighbors, and computes the ratio with the cells planar area.
Larger values indicate steeper terrain, minimum value is 1 for flat regions.

\param w window size. Used to position the 8 neighors at its corners.
*/
ScalarField2 HeightField::AreaRatio(int w) const
{
    ScalarField2 res(domain, nx, ny, 0.0);

    // compute cell areas once
    ScalarField2 cellTrianglesArea(domain, nx, ny, 0);
    for (int i = 0; i < nx - 1; i++)
    {
      for (int j = 0; j < ny - 1; j++)
      {
        // 2x2 vertices
        Vector3 v00 = Vertex(i, j);
        Vector3 v10 = Vertex(i + 1, j);
        Vector3 v11 = Vertex(i + 1, j + 1);
        Vector3 v01 = Vertex(i, j + 1);

        // 3 edges, 2 triangles
        Vector3 e10 = v10 - v00;
        Vector3 e11 = v11 - v00;
        Vector3 e01 = v01 - v00;

        // surface areas
        cellTrianglesArea(i, j) += 0.5 * Norm(cross(e10, e11));
        cellTrianglesArea(i, j) += 0.5 * Norm(cross(e11, e01));
      }
      cellTrianglesArea(i, ny - 1) = cellTrianglesArea(i, ny - 2);
    }
    for (int j = 0; j < ny; j++) {
      cellTrianglesArea(nx - 1, j) = cellTrianglesArea(nx - 2, j);
    }
    double cellPlanarArea = cellSize[0] * cellSize[1];


    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        double areaTris = 0.0;
        double areaPlane = 0.0;

        IndexArea area = indexAreaFromWindow(i, j, w);
        for (int wi = area.xmin(); wi <= area.xmax(); wi++)
        {
          for (int wj = area.ymin(); wj <= area.ymax(); wj++)
          {
            areaTris += cellTrianglesArea(wi, wj);
            areaPlane += cellPlanarArea;
          }
        }

        res(i, j) = areaTris / areaPlane;
      }
    }
    return res;
}


/*!
\brief Compute the TPI, an index that measures how much a cell rises with respect to the mean elevation
of a given radius. It is very similar in spirit to Peakedness, but here we obtain a value in m instead
of the % of cells below it.

This function is widely used (or was) to classify landforms, and it is integrated in many GIS packages.

See Weiss 2001: Topographic position and landforms analysis (poster)
And Wilson and Gallant 2000: Primary Topographic Attributes (chapter 3 in "Terrain Analysis: Principles and Applications")

\param r Radius of the disk (we actually use a squared window).
*/
ScalarField2 HeightField::TopographicPositionIndex(int w) const
{
  ScalarField2 res(domain, nx, ny, 0.0);
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      int count = 0;
      double sum = 0.0;
      IndexArea area = indexAreaFromWindow(i, j, w);
      for (int wi = area.xmin(); wi <= area.xmax(); wi++)
      {
        for (int wj = area.ymin(); wj <= area.ymax(); wj++)
        {
          count++;
          sum += at(wi, wj);
        }
      }
      res(i, j) = at(i, j) - sum / count;
    }
  }
  return res;
}


/*!
\brief Compute the TPI, an index that measures how much a cell rises with respect to the mean elevation
of a given radius. It is very similar in spirit to Peakedness, but here we obtain a value in m instead
of the % of cells below it.

Optimized version for large w using a Summed Area Table.

\param r Radius of the disk (we actually use a squared window).
*/
ScalarField2 HeightField::TopographicPositionIndexSAT(int w) const
{
  int r = w/2;

  // build SAT
  ScalarField2 sat = this->summedAreaTable();

  // compute TPI
  ScalarField2 res(domain, nx, ny, 0.0);
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      int imin = std::max(i - r, 0);
      int jmin = std::max(j - r, 0);
      int imax = std::min(i + r, nx - 1);
      int jmax = std::min(j + r, ny - 1);

      double sum = sat[cellId(imax, jmax)]
        - (imin > 0 ? sat[cellId(imin - 1, jmax)] : 0)
        - (jmin > 0 ? sat[cellId(imax, jmin - 1)] : 0)
        + (imin > 0 && jmin > 0 ? sat[cellId(imin - 1, jmin - 1)] : 0);
      int count  = (imax - imin + 1) * (jmax - jmin + 1);

      res(i, j) = at(i, j) - sum / count;
    }
  }
  return res;
}


/*!
\brief Compute the TRI, an index that measures the ruggedness of a terrain neighborhood by computing the RMS of the
residuals with respect to the center location.

\param r Radius of the disk (we actually use a squared window).
*/
ScalarField2 HeightField::RuggednessIndex(int w) const
{
  ScalarField2 res(domain, nx, ny, 0.0);
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double hp = at(i, j);
      double ss = 0.0;
      IndexArea area = indexAreaFromWindow(i, j, w);
      for (int wi = area.xmin(); wi <= area.xmax(); wi++)
      {
        for (int wj = area.ymin(); wj <= area.ymax(); wj++)
        {
          double dif = at(wi, wj) - hp;
          ss += dif*dif;
        }
      }
      res(i, j) = std::sqrt(ss);
    }
  }
  return res;
}


/*!
\brief Compute the surface roughness.

Lindsay et al 2019. Scale-optimized surface roughness for topographic analysis.

The authors associate surface normals dispersion to roughness, a measure of texture complexity.
In particular, they measure the spherical standard deviation in a window w &times; w,
i.e. the angular spread of the normal directions in the window.
It will be 0 for flat and inclined planes, and increase with surface texture complexity.

They also propose to remove lower frequencies by doing a gaussian low-pass filter on the DEM,
to obtain a signature that depends only on current scale (otherwise, it increases monotonically with w).

\param w Window size.
\param b Spherical deviation. If true, return spherical standard deviation. Otherwise, 1 - |R|/N
*/
ScalarField2 HeightField::SurfaceRoughness(int w, bool dev) const
{
    // Build normal map
    std::vector<Vector3> normals = std::vector<Vector3>(nx * ny);
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        normals[cellId(i, j)] = Normal(i, j);
      }
    }

    // compute roughness
    int r = w / 2;
    ScalarField2 s(domain, nx, ny, 0.0);
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        int imin = std::max(i - r, 0);
        int jmin = std::max(j - r, 0);
        int imax = std::min(i + r, nx - 1);
        int jmax = std::min(j + r, ny - 1);
        int n = (imax - imin + 1) * (jmax - jmin + 1);

        Vector3 vsum(0);
        for (int di = imin; di <= imax; di++) {
          for (int dj = jmin; dj <= jmax; dj++) {
            vsum = vsum + normals[cellId(di, dj)];
          }
        }
        double R = Norm(vsum);

        if (dev)
            s(i, j) = std::sqrt(-2.0 * log(R / n)) * 180.0 / Math::Pi;
        else
            s(i, j) = 1.0 - R / n;
      }
    }

    return s;
}


/*!
\brief Compute the surface roughness.

Lindsay et al 2019. Scale-optimized surface roughness for topographic analysis.

Optimized version for large w using a Summed Area Table.

\param w Window size.
\param b Spherical deviation. If true, return spherical standard deviation. Otherwise, 1 - |R|/N
*/
ScalarField2 HeightField::SurfaceRoughnessSAT(int w, bool dev) const
{
  int r = w / 2;

  // Build normal SAT
  std::vector<Vector3> normalSAT = std::vector<Vector3>(nx * ny);
  normalSAT[0] = Normal(0, 0);
  for (int i = 1; i < nx; i++)
    normalSAT[cellId(i, 0)] = Normal(i, 0) + normalSAT[cellId(i - 1, 0)];
  for (int j = 1; j < ny; j++)
    normalSAT[cellId(0, j)] = Normal(0, j) + normalSAT[cellId(0, j - 1)];
  for (int i = 1; i < nx; i++)
  {
    for (int j = 1; j < ny; j++)
    {
      normalSAT[cellId(i, j)] = Normal(i, j)
        + normalSAT[cellId(i - 1, j)]
        + normalSAT[cellId(i, j - 1)]
        - normalSAT[cellId(i - 1, j - 1)];
    }
  }

  // compute roughness
  ScalarField2 s(domain, nx, ny, 0.0);
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      int imin = std::max(i - r, 0);
      int jmin = std::max(j - r, 0);
      int imax = std::min(i + r, nx - 1);
      int jmax = std::min(j + r, ny - 1);
      int n = (imax - imin + 1) * (jmax - jmin + 1);

      Vector3 nsum = normalSAT[cellId(imax, jmax)]
        - (imin > 0 ? normalSAT[cellId(imin - 1, jmax)] : Vector3(0))
        - (jmin > 0 ? normalSAT[cellId(imax, jmin - 1)] : Vector3(0))
        + (imin > 0 && jmin > 0 ? normalSAT[cellId(imin - 1, jmin - 1)] : Vector3(0));
      double R = Norm(nsum);

      if (dev)
          s(i, j) = std::sqrt(-2.0 * log(R / n)) * 180.0 / Math::Pi;
      else
          s(i, j) = 1.0 - R / n;
    }
  }

  return s;
}


/*!
\brief Compute the hillslope asymmetry map for a given direction

Poulos et al 2012. Hillslope asymmetry maps reveal widespread, multi-scale organization

The authors define hillsope asymmetry as the difference in median slope between slopes of
opposite orientation inside an area.

NOTE: I used average slopes instead of median. First, for a reason of efficiency. Second,
it also makes more sense in my opinion. In any case, results visually very similar.
It is possible to change the behavior in the code, just set MEDIAN_ASYMMETRY = true

\param w Window in which analysis is performed.
\param Direction angle (in degrees, 0 is +X).
\param Tolerance how much we can deviate from this angle.

\author Oscar Argudo
*/
ScalarField2 HeightField::HillslopeAsymmetry(int w, double direction, double tolerance) const
{
  const double MIN_SLOPE = 0.05;
  const bool MEDIAN_ASYMMETRY = false;

  // compute slope and aspect maps
  double dir = Math::DegreeToRadian(direction);
  double tol = Math::DegreeToRadian(tolerance);
  double dirPositive = dir;
  double dirNegative = dir + Math::Pi;
  if (dirNegative > 2*Math::Pi)
  {
    dirNegative -= 2*Math::Pi;
  }

  ScalarField2 slope(domain, nx, ny, 0);
  std::vector<bool> positiveAspect(nx * ny, false);
  std::vector<bool> negativeAspect(nx * ny, false);

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      Vector2 g = Gradient(i, j);
      slope(i, j) = Norm(g);

      // authors skip flat areas (< 5% slope)
      if (slope.at(i, j) > MIN_SLOPE) {
        double aspect = g.Angle() + Math::Pi;

        double dp = aspect - dirPositive;
        if (dp > Math::Pi) dp -= 2 * Math::Pi;
        if (dp < -Math::Pi) dp += 2 * Math::Pi;
        if (std::abs(dp) < tol)
        {
          positiveAspect[slope.cellId(i, j)] = 1;
        }

        double dn = aspect - dirNegative;
        if (dn > Math::Pi) dn -= 2 * Math::Pi;
        if (dn < -Math::Pi) dn += 2 * Math::Pi;
        if (std::abs(dn) < tol)
        {
          negativeAspect[slope.cellId(i, j)] = 1;
        }
      }
    }
  }

  // compute asymmetry
  const int r = std::max(w / 2, 1);

  ScalarField2 ha(domain, nx, ny, 0);

  // the authors always talk about medians...
  if (MEDIAN_ASYMMETRY)
  {
    for (int x = 0; x < nx; x++)
    {
      for (int y = 0; y < ny; y++)
      {
        if (true || slope.at(x, y) > MIN_SLOPE)
        {
          int imin = std::max(0, x - r);
          int jmin = std::max(0, y - r);
          int imax = std::min(x + r, nx - 1);
          int jmax = std::min(y + r, ny - 1);

          std::vector<double> posSlopes;
          posSlopes.reserve(r * r);
          std::vector<double> negSlopes;
          negSlopes.reserve(r * r);

          for (int i = imin; i <= imax; i++)
          {
            for (int j = jmin; j <= jmax; j++)
            {
              if (positiveAspect[slope.cellId(i, j)])
              {
                posSlopes.push_back(slope.at(i, j));
              }
              if (negativeAspect[slope.cellId(i, j)])
              {
                negSlopes.push_back(slope.at(i, j));
              }
            }
          }

          // find medians
          if (posSlopes.size() > 0 && negSlopes.size() > 0)
          {
            std::nth_element(posSlopes.begin(), posSlopes.begin() + posSlopes.size() / 2, posSlopes.end());
            double posMedian = posSlopes[posSlopes.size() / 2];

            std::nth_element(negSlopes.begin(), negSlopes.begin() + negSlopes.size() / 2, negSlopes.end());
            double negMedian = negSlopes[negSlopes.size() / 2];

            // compute logarithmic asymmetry of medians
            // they do it this way to have opposite signs as: log (3 / 1) = - log (1 / 3)
            ha(x, y) = std::log10(posMedian / negMedian);
          }
        }
      }
    }
  }

  // ... although I think means are also valid and more efficient to obtain
  else {

    // compute SAT
    ScalarField2 slopePositiveSAT(domain, nx + 1, ny + 1, 0);
    ScalarField2 slopeNegativeSAT(domain, nx + 1, ny + 1, 0);
    ScalarField2 cntPositiveSAT(domain, nx + 1, ny + 1, 0);
    ScalarField2 cntNegativeSAT(domain, nx + 1, ny + 1, 0);

    for (int i = 1; i <= nx; i++)
    {
      for (int j = 1; j <= ny; j++)
      {
        slopePositiveSAT(i, j) = slope.at(i - 1, j - 1) * positiveAspect[slope.cellId(i - 1, j - 1)]
          + slopePositiveSAT.at(i - 1, j)
          + slopePositiveSAT.at(i, j - 1)
          - slopePositiveSAT.at(i - 1, j - 1);

        slopeNegativeSAT(i, j) = slope.at(i - 1, j - 1) * negativeAspect[slope.cellId(i - 1, j - 1)]
          + slopeNegativeSAT.at(i - 1, j)
          + slopeNegativeSAT.at(i, j - 1)
          - slopeNegativeSAT.at(i - 1, j - 1);

        cntPositiveSAT(i, j) = positiveAspect[slope.cellId(i - 1, j - 1)]
          + cntPositiveSAT.at(i - 1, j)
          + cntPositiveSAT.at(i, j - 1)
          - cntPositiveSAT.at(i - 1, j - 1);

        cntNegativeSAT(i, j) = negativeAspect[slope.cellId(i - 1, j - 1)]
          + cntNegativeSAT.at(i - 1, j)
          + cntNegativeSAT.at(i, j - 1)
          - cntNegativeSAT.at(i - 1, j - 1);
      }
    }

    // compute asymmetry
    for (int x = 0; x < nx; x++)
    {
      for (int y = 0; y < ny; y++)
      {
        if (true || slope.at(x, y) > MIN_SLOPE)
        {
          int imin = std::max(0, x - r);
          int jmin = std::max(0, y - r);
          int imax = std::min(x + r, nx);
          int jmax = std::min(y + r, ny);

          double sumPosSlope = slopePositiveSAT.at(imax, jmax) + slopePositiveSAT.at(imin, jmin)
            - slopePositiveSAT.at(imin, jmax) - slopePositiveSAT.at(imax, jmin);
          double sumNegSlope = slopeNegativeSAT.at(imax, jmax) + slopeNegativeSAT.at(imin, jmin)
            - slopeNegativeSAT.at(imin, jmax) - slopeNegativeSAT.at(imax, jmin);
          int cntPosSlope = cntPositiveSAT.at(imax, jmax) + cntPositiveSAT.at(imin, jmin)
            - cntPositiveSAT.at(imin, jmax) - cntPositiveSAT.at(imax, jmin);
          int cntNegSlope = cntNegativeSAT.at(imax, jmax) + cntNegativeSAT.at(imin, jmin)
            - cntNegativeSAT.at(imin, jmax) - cntNegativeSAT.at(imax, jmin);

          if (cntPosSlope > 0 && cntNegSlope > 0) {
            double avgPosSlope = sumPosSlope / cntPosSlope;
            double avgNegSlope = sumNegSlope / cntNegSlope;
            ha(x, y) = std::log10(avgPosSlope / avgNegSlope);
          }
        }
      }
    }
  }

  return ha;
}



