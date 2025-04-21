#include "heightfield.h"


int DikauWookClass(bool sloping, int signCurvTan, int signCurvPro, int signCurvMin, int signCurvMax)
{
  // sloping area landforms
  if (sloping) {
    return 1 + (signCurvPro + 1) * 3 + signCurvTan + 1;
  }
  // flat area landforms
  else {
    if (signCurvMax < 0) {
      if (signCurvMin < 0) return 10;
    }
    else if (signCurvMax == 0) {
      if (signCurvMin < 0) return 11;
      else if (signCurvMin == 0) return 12;
    }
    else {
      if (signCurvMin < 0) return 13;
      else if (signCurvMin == 0) return 14;
      else return 15;
    }
  }
  return 0;
}


/*!
\brief Classification of terrain features based on the curvatures

Following Dikau and Wood classification matrix based on the tangential, profile, maximum and minimum curvatures.

\param w size of domain w x w in which a quadric is fitted to compute its curvatures.
\param tSlope threshold to differentiate between sloping and flat areas.
\param tCurv absolute value of thre threshold to consider positive, zero or negative curvature.
*/
IntField2 HeightField::LandformsDikauWood(int w, double tSlope, double tCurvs) const
{
    IntField2 landform(getDomain(), nx, ny, 0);
    for (int i = 0; i < nx; i++)
    {
      for (int j = 0; j < ny; j++)
      {
        QuadricSurface q = FitQuadric(i, j, w, true);
        double a = q(2, 0);
        double b = q(0, 2);
        double c = q(1, 1);
        double d = q(1, 0);
        double e = q(0, 1);
        double slope = sqrt((d * d) + (e * e));
        double curvTan = !d && !e ? 0.0 : 2.0 * (b * d * d + a * e * e - c * d * e) / (d * d + e * e);
        double curvPro = !d && !e ? 0.0 : 2.0 * (a * d * d + b * e * e + c * e * d) / ((e * e + d * d) * pow(1.0 + d * d + e * e, 1.5));
        double curvMin = a + b - sqrt((a - b) * (a - b) + c * c);
        double curvMax = a + b + sqrt((a - b) * (a - b) + c * c);
        int signCurvTan = curvTan < -tCurvs ? -1 : (curvTan > tCurvs ? 1 : 0);
        int signCurvPro = curvPro < -tCurvs ? -1 : (curvPro > tCurvs ? 1 : 0);
        int signCurvMin = curvMin < -tCurvs ? -1 : (curvMin > tCurvs ? 1 : 0);
        int signCurvMax = curvMax < -tCurvs ? -1 : (curvMax > tCurvs ? 1 : 0);
        landform(i, j) = DikauWookClass(slope > tSlope, signCurvTan, signCurvPro, signCurvMin, signCurvMax);
      }
    }

    return landform;
}


/*!
\brief Classification of terrain features based on the curvatures on multiple scales

Following Fisher, Wood and Cheng 2004: "Where is Helvellyn? fuzziness of multi-scale landscape morphometry"

\param entropy output scalarfield for the entropy of each classification result.
\param minScale,maxScale minimum and maximum scales to consider.
\param tSlope,tCurv parameters for Dikau-Wood classification on a single scale.
*/
IntField2 HeightField::LandformsFuzzyDW(ScalarField2& entropy, int minScale, int maxScale, double tSlope, double tCurv) const
{
  int NUM_DW_CLASSES = 16;

  int nx = getSizeX();
  int ny = getSizeY();

  // multiscale analysis
  std::vector<std::vector<int> > classes(getNumElements(), std::vector<int>(NUM_DW_CLASSES, 0));
  int numScales = 0;
  for (int w = minScale; w <= maxScale; w += 2) {
    IntField2 landforms = LandformsDikauWood(w, tSlope, tCurv);
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        classes[cellId(i, j)][landforms.at(i, j)]++;
      }
    }
    numScales++;
  }

  // keep the most repeated class
  IntField2 landform(getDomain(), nx, ny, 0);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      int cntMax = 0;
      double e = 0;
      for (int k = 0; k < NUM_DW_CLASSES; k++) {
        int cnt = classes[cellId(i, j)][k] > cntMax;
        if (cnt > cntMax) {
          cntMax = cnt;
          landform(i, j) = k;
        }
        double pk = double(cnt) / double(numScales);
        if (cnt > 0) e += pk * std::log2(pk);
      }
      entropy(i, j) = -e;
    }
  }

  return landform;
}

/*!
\brief Classification of terrain features based on the topographic position index.

See the poster presentation by Weiss 2001: Topographic Position and Landforms Analysis

\param radiusSmall,radiusLarge Radius of the small-scale and large-scale topographic position index.
\param flatSlope threshold below which a landform with no curvatures is considered flat instead of planar slope
*/
IntField2 HeightField::LandformsTPI(double radiusSmall, double radiusLarge, double flatSlope) const
{
  int ws = int(radiusSmall / cellSize[0])*2 + 1;
  int wl = int(radiusLarge / cellSize[0])*2 + 1;

  ScalarField2 tpiSmall = TopographicPositionIndex(ws);
  ScalarField2 tpiLarge = TopographicPositionIndex(wl);
  ScalarField2 slope = Slope();

  double meanS = tpiSmall.average();
  double meanL = tpiLarge.average();
  double stdS = tpiSmall.standardDeviation(meanS);
  double stdL = tpiLarge.standardDeviation(meanL);

  IntField2 c(getDomain(), nx, ny, 0);
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double vs = (tpiSmall.at(i, j) - meanS) / stdS;
      double vl = (tpiLarge.at(i, j) - meanL) / stdL;

      if (vs <= -1 && vl <= -1)         c(i, j) = 1; // V-shape river valleys, deep narrow canyons
      else if (vs <= -1 && abs(vl) < 1) c(i, j) = 2; // lateral midslope incised drainage, local valleys in plains
      else if (vs <= -1 && vl >= 1)     c(i, j) = 3; // upland incised drainage, stream headwaters
      else if (abs(vs) < 1 && vl <= -1) c(i, j) = 4; // U-shape valleys
      else if (abs(vs) < 1 && abs(vl) < 1) {
        if (slope.at(i, j) < flatSlope) c(i, j) = 5; // broad flat areas
        else                            c(i, j) = 6; // broad open slopes
      }
      else if (abs(vs) < 1 && vl >= 1)  c(i, j) = 7; // flat ridge tops, mesa tops
      else if (vs >= 1 && vl <= -1)     c(i, j) = 8; // local ridge/hilltops within broad valleys
      else if (vs >= 1 && abs(vl) < 1)  c(i, j) = 9; // lateral midslope drainage divides, local ridges in plains
      else if (vs >= 1 && vl >= 1)      c(i, j) = 10; // mountain tops, high narrow ridges
    }
  }
  return c;
}

/*!
\brief  Implementation of the peak/valley detection using the black top hat transform.

See Rodriguez et al. 2002, The Black Top Hat function applied to a DEM.

\param w size of the w x w window for the opening and closing operators.
\param tPeak threshold to highlight peaks.
\param tValley threshold to highlight valleys.
*/
ScalarField2 HeightField::BlackWhiteTopHatTransform(int w, double tPeak, double tValley) const
{
  ScalarField2 opening = this->minFilter(w).maxFilter(w);
  ScalarField2 closing = this->maxFilter(w).minFilter(w);

  ScalarField2 wth(getDomain(), nx, ny);
  ScalarField2 bth(getDomain(), nx, ny);
  for (int i = 0; i < nx*ny; i++) {
      wth[i] = at(i) - opening[i];
      bth[i] = closing[i] - at(i);
  }

  ScalarField2 mountValley(getDomain(), nx, ny, 0);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      bool peak = wth.at(i, j) - tPeak > 0;
      bool vall = bth.at(i, j) - tValley > 0;
      if (peak && !vall) mountValley(i, j) = 1;
      if (vall && !peak) mountValley(i, j) = -1;
    }
  }
  return mountValley;
}

/*!
\brief Compute the geomorphon index map.

Geomorphon classify cells into ten categories: flat, peak, ridge, shoulder, hollow, slope, spur, footslope, valley, and pit.

\param maxDist maximum ray traversal distance.
\param flatTangent threshold for considering a ray angle as flat.
*/
IntField2 HeightField::Geomorphons(double maxDist, double flatTangent) const
{
  IntField2 geo(getDomain(), nx, ny, -1);

  // Geomorphon index lookup table
  // FL 0 - flat; PK 1 - peak; RI 2- ridge; SH 3 - shoulder; HL 4 - hollow=convex; SL 5 - slope;
  // SP 6 - spur=concave; FS 7 - footslope; VL 8 - valley; PT 9 - pit.
  static const int index[9][9] =
  {
    { 0, 0, 0, 7, 7, 8, 8, 8, 9},
    { 0, 0, 7, 7, 7, 8, 8, 8,-1},
    { 0, 3, 5, 5, 4, 4, 8,-1,-1},
    { 3, 3, 5, 5, 5, 4,-1,-1,-1},
    { 3, 3, 6, 5, 5,-1,-1,-1,-1},
    { 2, 2, 6, 6,-1,-1,-1,-1,-1},
    { 2, 2, 2,-1,-1,-1,-1,-1,-1},
    { 2, 2,-1,-1,-1,-1,-1,-1,-1},
    { 1,-1,-1,-1,-1,-1,-1,-1,-1},
  };

  const double dt = cellSize[0];

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double zp = at(i, j);

      // Number of pluses and minuses
      int np = 0;
      int nm = 0;
      for (int k = 0; k < 8; k++)
      {
        double t = dt;
        Vector2 dir = Vector2(next[k].x(), next[k].y());
        Vector2 p0 = domainCoords(i, j);
        double slope = -1e12;

        // traverse terrain
        Vector2 p = p0 + t * dir;
        while ((maxDist < 0 || t < maxDist) && isInDomain(p)) {
          double dh = value(p) - zp;
          double s = dh / t;
          slope = std::max(slope, s);
          t += dt;
          p = p0 + t * dir;
        }

        int s = 1 + Math::IntegerSign(slope, flatTangent);

        if (s == 0) { nm++; }
        if (s == 2) { np++; }
      }

      // Geomorphon type
      geo(i, j) = index[nm][np];
    }
  }
  return geo;
}
