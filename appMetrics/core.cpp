#include "core.h"


const Index2 FieldGrid2D::next[8] = { Index2(1, 0), Index2(1, 1), Index2(0, 1), Index2(-1, 1), Index2(-1, 0), Index2(-1, -1), Index2(0, -1), Index2(1, -1) };
const double FieldGrid2D::length[8] = { 1.0, sqrt(2.0), 1.0, sqrt(2.0), 1.0, sqrt(2.0), 1.0, sqrt(2.0) };
const double FieldGrid2D::inverselength[8] = { 1.0, 1.0 / sqrt(2.0), 1.0, 1.0 / sqrt(2.0), 1.0, 1.0 / sqrt(2.0), 1.0, 1.0 / sqrt(2.0) };


Matrix3 Inverse(const Matrix3& A)
{
  double e = A.Determinant();
  Matrix3 inv;
  if (fabs(e) < 1.0e-18)
  {
      for (int i = 0; i < 9; i++) inv.r[i] /= e;
  }
  return inv;
}

bool Box2::Intersect(const Vector2 &s0, const Vector2 &s1, double &tmin, double &tmax)
{
    const double epsilon = 1.0e-5;

    tmin = -1e16;
    tmax = 1e16;

    const Vector2& a = bmin;
    const Vector2& b = bmax;
    Vector2 p = s0;
    Vector2 d = s1 - s0;

    double t;
    // Ox
    if (d[0] < -epsilon) {
        t = (a[0] - p[0]) / d[0];
        if (t < tmin)
            return false;
        if (t <= tmax)
            tmax = t;
        t = (b[0] - p[0]) / d[0];
        if (t >= tmin) {
            if (t > tmax)
                return false;
            tmin = t;
        }
    }
    else if (d[0] > epsilon) {
        t = (b[0] - p[0]) / d[0];
        if (t < tmin)
            return false;
        if (t <= tmax)
            tmax = t;
        t = (a[0] - p[0]) / d[0];
        if (t >= tmin) {
            if (t > tmax)
                return false;
            tmin = t;
        }
    }
    else if (p[0]<a[0] || p[0]>b[0])
        return false;

    // Oy
    if (d[1] < -epsilon) {
        t = (a[1] - p[1]) / d[1];
        if (t < tmin)
            return false;
        if (t <= tmax)
            tmax = t;
        t = (b[1] - p[1]) / d[1];
        if (t >= tmin) {
            if (t > tmax)
                return false;
            tmin = t;
        }
    }
    else if (d[1] > epsilon) {
        t = (b[1] - p[1]) / d[1];
        if (t < tmin)
            return false;
        if (t <= tmax)
            tmax = t;
        t = (a[1] - p[1]) / d[1];
        if (t >= tmin) {
            if (t > tmax)
                return false;
            tmin = t;
        }
    }
    else if (p[1]<a[1] || p[1]>b[1])
        return false;

    return true;
}


bool Box3::Intersect(const Ray& ray, double& tmin, double& tmax) const
{
  const double epsilon = 1.0e-5;

  tmin = -1e16;
  tmax = 1e16;

  Vector3 p = ray.origin();
  Vector3 d = ray.direction();

  double t;
  // Ox
  if (d[0] < -epsilon)
  {
    t = (bmin[0] - p[0]) / d[0];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (bmax[0] - p[0]) / d[0];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (d[0] > epsilon)
  {
    t = (bmax[0] - p[0]) / d[0];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (bmin[0] - p[0]) / d[0];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (p[0]<bmin[0] || p[0]>bmax[0])
    return false;

  // Oy
  if (d[1] < -epsilon)
  {
    t = (bmin[1] - p[1]) / d[1];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (bmax[1] - p[1]) / d[1];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (d[1] > epsilon)
  {
    t = (bmax[1] - p[1]) / d[1];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (bmin[1] - p[1]) / d[1];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (p[1]<bmin[1] || p[1]>bmax[1])
    return false;

  // Oz
  if (d[2] < -epsilon)
  {
    t = (bmin[2] - p[2]) / d[2];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (bmax[2] - p[2]) / d[2];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (d[2] > epsilon)
  {
    t = (bmax[2] - p[2]) / d[2];
    if (t < tmin)
      return false;
    if (t <= tmax)
      tmax = t;
    t = (bmin[2] - p[2]) / d[2];
    if (t >= tmin)
    {
      if (t > tmax)
        return false;
      tmin = t;
    }
  }
  else if (p[2]<bmin[2] || p[2]>bmax[2])
    return false;

  return true;
}



Camera::Camera()
{
    Camera::eye = Vector3(0.0);
    Camera::at = Vector3(0.0, 1.0, 0.0);
    Camera::up = Vector3(0.0, 0.0, 1.0);

    // Near and far planes
    Camera::nearplane = 1.0;
    Camera::farplane = 1000.0;

    // Aperture
    Camera::cah = 0.980;
    Camera::cav = 0.735;
    Camera::fl = 35.0;
}

Camera::Camera(const Vector3& eye, const Vector3& at, const Vector3& up, double near, double far)
{
    Camera::eye = eye;
    Camera::at = at;
    Camera::up = up;

    // Near and far planes
    Camera::nearplane = near;
    Camera::farplane = far;

    // Aperture
    Camera::cah = 0.980;
    Camera::cav = 0.735;
    Camera::fl = 35.0;
}

double Camera::getAngleOfViewH(double, double) const
{
    return 2.0 * atan(cah * 25.4 * 0.5 / fl);
}

double Camera::getAngleOfViewV(double w, double h) const
{
    double avh = getAngleOfViewH(w, h);
    return 2.0 * atan(tan(avh / 2.0) * double(h) / double(w));
}

void Camera::upDownRound(double a)
{
    Vector3 z = at - eye;
    double length = Norm(z);
    z = z/length;
    Vector3 left = Normalized(cross(up, z));

    // Rotate
    z = z * cos(a) + up * sin(a);

    // Update Vector
    up = cross(z, left);
    eye = at - z * length;
}

void Camera::leftRightRound(double a)
{
    Vector3 e = eye - at;
    Vector3 left = cross(up, e);
    e = Vector3(e[0] * cos(a) - e[1] * sin(a), e[0] * sin(a) + e[1] * cos(a), e[2]);
    left = Vector3(left[0] * cos(a) - left[1] * sin(a), left[0] * sin(a) + left[1] * cos(a), 0.0);
    up = Normalized(cross(left, -e));
    eye = at + e;
}

void Camera::backForth(double a, bool moveAt)
{
    Vector3 z = at - eye;
    double length = Norm(z);
    z = z/length;
    eye = eye + a * z;
    if (moveAt) {
        at = at + a * z;
    }
}

void Camera::upDownPlane(double a)
{
    Vector3 z = at - eye;
    double length = Norm(z);
    z = z/length;
    Vector3 left = Normalized(cross(Vector3(0, 0, 1), z));

    eye = eye + a * cross(z, left);
    at = at + a * cross(z, left);
}

void Camera::leftRightPlane(double a)
{
    Vector3 z = at - eye;
    z[2] = 0.0;
    double length = Norm(z);
    z = z/length;
    Vector3 left = Normalized(cross(Vector3(0, 0, 1), z));

    eye = eye + a * left;
    at = at + a * left;
}

Ray Camera::pixelToRay(int px, int py, int w, int h) const
{
    // Get coordinates
    Vector3 view = getViewDir();
    Vector3 horizontal = Normalized(cross(view, up));
    Vector3 vertical = Normalized(cross(horizontal, view));

    double length = 1.0;

    // Convert to radians
    double rad = getAngleOfViewV(w, h);  // fov

    double vLength = tan(rad / 2.0) * length;
    double hLength = vLength * (double(w) / double(h));
    vertical = vertical*vLength;
    horizontal = horizontal*hLength;

    // Translate mouse coordinates so that the origin lies in the center of the view port
    double x = px - w / 2.0;
    double y = h / 2.0 - py;

    // Scale mouse coordinates so that half the view port width and height becomes 1.0
    x /= w / 2.0;
    y /= h / 2.0;

    // Direction is a linear combination to compute intersection of picking ray with view port plane
    return Ray(eye, Normalized(view * length + horizontal * x + vertical * y));
}

Camera Camera::View(const Box3& box)
{
    Vector3 v = 0.5*(box.getMax() - box.getMin());
	v[2] = 0;
	double r = Norm(v);
	v = 2*v;
	v[2] = -r;
    return Camera(box.center() - v, box.center(), Vector3(0.0, 0.0, 1.0), r, 3*r);
}


Vector3 ColorPalette::getColor(double t) const {
    if (colors.size() == 0) return Vector3(1);
    if (colors.size() == 1) return colors[0];

    if (t <= anchors.front()) return colors.front();
    if (t >= anchors.back()) return colors.back();
    for (int i = 0; i < int(colors.size() - 1); i++) {
        if (t < anchors[i+1]) {
            double s = Math::LinearStep(t, anchors[i], anchors[i+1]);
            return (1-s)*colors[i] + s*colors[i+1];
        }
    }

    return Vector3(0);
}

Vector3 LookupPalette::getColor(int i) const {
    if ((i < 0) || (i >= int(c.size()))) return Vector3(0,0,0);
    return c.at(i);
}


/*!
\brief Get the days to J2000.

Only works between 1901 to 2099.
\param h Universal time, in decimal hours.
\param y, m, d Year, month and day.
*/
double Sun::Day(int y, int m, int d, double h)
{
  int luku = -7 * (y + (m + 9) / 12) / 4 + 275 * m / 9 + d;

  luku += (long int)y * 367;
  return (double)luku - 730530.0 + h / 24.0;
}


/*!
\brief Compute the hour-angle.
\param lat, declin Latitude and declination.
*/
double Sun::f0(double lat, double declin)
{
  const double AirRefr = 34.0 / 60.0; // athmospheric refraction degrees //

  double fo, dfo;
  const double SunDia = 0.53;     // Sunradius degrees

  dfo = Math::DegreeToRadian(0.5 * SunDia + AirRefr);
  if (lat < 0.0) dfo = -dfo;	// Southern hemisphere
  fo = tan(declin + dfo) * tan(Math::DegreeToRadian(lat));
  if (fo > 0.99999) fo = 1.0; // to avoid overflow //
  fo = asin(fo) + 0.5*Math::Pi;
  return fo;
}


/*!
Compute the ecliptic longitude of the sun.
*/
double Sun::FNsun(double d, double& L, double& RA, double& delta)
{
  const double rads = Math::toRads;
  const double degs = Math::toDegs;
  double w, M, v, r;
  //   mean longitude of the Sun
  w = 282.9404 + 4.70935E-5 * d;
  M = 356.047 + 0.9856002585 * d;

  // Sun's mean longitude
  L = Math::Angle(Math::DegreeToRadian(w + M));

  //   mean anomaly of the Sun

  double g = Math::Angle(Math::DegreeToRadian(M));

  // eccentricity
  double ecc = 0.016709 - 1.151E-9 * d;

  //   Obliquity of the ecliptic

  double obliq = 23.4393 * rads - 3.563E-7 * rads * d;
  double E = M + degs * ecc * sin(g) * (1.0 + ecc * cos(g));
  E = degs * Math::Angle(E * rads);
  double x = cos(E * rads) - ecc;
  double y = sin(E * rads) * sqrt(1.0 - ecc * ecc);
  r = sqrt(x * x + y * y);
  v = atan2(y, x) * degs;
  // longitude of sun
  double lonsun = v + w;
  lonsun -= 360.0 * (lonsun > 360.0);

  // sun's ecliptic rectangular coordinates
  x = r * cos(lonsun * rads);
  y = r * sin(lonsun * rads);
  double yequat = y * cos(obliq);
  double zequat = y * sin(obliq);
  RA = atan2(yequat, x);
  delta = atan2(zequat, sqrt(x * x + yequat * yequat));
  RA *= degs;

  //   Ecliptic longitude of the Sun

  return Math::Angle(L + 1.915 * rads * sin(g) + 0.02 * rads * sin(2 * g));
};

/*!
\brief Compute azimuth and altitude angles of Sun given a location and datetime.

Azimuth is measured in degrees clockwise from north.
Elevation is measured in degrees up from the horizon.

\param latit,longit Latitude and longitude.
\param year, month, day, hour Explicit.
\param azim, altit Returned azimth and altitude.
*/
void Sun::Angles(const double& latit, const double& longit, int year, int month, int day, const double& hour, double& azim, double& altit)
{
  const double rads = Math::toRads;
  const double degs = Math::toDegs;
  const double tzone = 0.0;
  double UT = hour - tzone;	// universal time
  double jd = Sun::Day(year, month, day, UT);

  double L, RA, delta;

  //   Find the ecliptic longitude of the sun
  FNsun(jd, L, RA, delta);

  //   Obliquity of the ecliptic
  // double obliq = 23.4393 * rads - 3.563E-7 * rads * jd;

  // Sidereal time at Greenwich meridian
  double GMST0 = L * degs / 15.0 + 12.0;	// hours
  double SIDTIME = GMST0 + UT + longit / 15.0;
  // Hour Angle
  double ha = 15.0 * SIDTIME - RA;	// degrees
  ha = Math::Angle(rads * ha);
  double x = cos(ha) * cos(delta);
  double y = sin(ha) * cos(delta);
  double z = sin(delta);
  double rlatit = Math::DegreeToRadian(latit);
  double xhor = x * sin(rlatit) - z * cos(rlatit);
  double yhor = y;
  double zhor = x * cos(rlatit) + z * sin(rlatit);
  // Azimuth angle
  azim = atan2(yhor, xhor) + Math::Pi;
  azim = Math::Angle(azim);

  // Altitude angle
  altit = asin(zhor);
}


std::minstd_rand Random::rng(1337);
std::uniform_real_distribution<double> Random::uniformDist = std::uniform_real_distribution<double>(0.0, 1.0);

