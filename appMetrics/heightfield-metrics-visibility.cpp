#include "heightfield.h"
#include "integerfield.h"


/*!
\brief Compute the Lipschitz constant of an elevation scalar field
*/
double HeightField::LipschitzK() const
{
  double k = 0.0;
  for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
          k = std::max(k, Norm(Gradient(i, j)));
      }
  }
  return k;
}


/*!
\brief Compute the intersection between a ray and the surface of the heightfield.

The algorithm uses a ray marching approach, therefore this function may require
many iterations and the resulting intersection may not be accurate.

\param ray The ray.
\param t Returned distance along the ray.
\param box %Box where intersection will be computed.
\param q Returned intersection point.
\param k Lipschitz constant.
\param length Maximum distance along the ray.
\param epsilon Minimum stepping distance.
*/
bool HeightField::Intersect(const Ray& ray, double& t, Vector3& q, const Box3& box, const double& k, const double& length, const double& epsilon) const
{
    double ta, tb;

    // Check the intersection with the bounding box
    if (!box.Intersect(ray, ta, tb)) return false;
    if (ta<-1.0e12 || tb>+1.0e12)    return false;

    tb = std::min(tb, length);
    t = std::max(ta + epsilon, 0.0);

    // Ray marching
    while (t < tb)
    {
        // Point along the ray
        Vector3 p = ray(t);

        // Heightfield elevation
        double z = value(p);
        double h = p[2] - z;
        if (h < epsilon)
        {
            q = Vector3(p[0], p[1], z);
            return true;
        }
        else
        {
            t += std::max(h / k, epsilon);
        }
    }
    return false;
}


/*!
\brief Viewshed from a given point.

\param viewpoint Observer location coordinates. Note the height is NOT relative to the ground.
\param viewCnt Output parameter, count of viewed locations.
*/
ScalarField2 HeightField::Viewshed(const Vector3& viewpoint, int& viewCnt) const
{
    double tolerance = Norm(cellSize);
    Box3 bbox = getBox();
    double klip = LipschitzK();

    ScalarField2 view(domain, nx, ny, 0);
    viewCnt = 0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            Vector3 q = Vertex(i, j);
            double dist = Norm(q - viewpoint);
            Vector3 dir = (q - viewpoint) / dist;

            Ray ray(viewpoint, dir);
            double t;
            Vector3 pintersect;
            Intersect(ray, t, pintersect, bbox, klip);

            if (t > dist - tolerance) {
                view(i, j) = 1;
                viewCnt++;
            }
        }
    }
    return view;
}


/*!
\brief Total viewshed or visibility index of the terrain.

\param outgoing Set to true/false to compute outgoing or incoming visibility.
\param offset Height of the observer, to be added to the ground height at each viewshed location.
*/
ScalarField2 HeightField::ViewshedTotal(bool outgoing, double offset) const
{
    ScalarField2 view(domain, nx, ny, 0);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            Vector3 vp = Vertex(i, j) + Vector3(0, 0, offset);
            int cnt;
            if (outgoing) {
                Viewshed(vp, cnt);
                view(i, j) = cnt;
            }
            else {
                view += Viewshed(vp, cnt);
            }
        }
    }
    return view;
}


/*!
\brief Random sampled version of Total Viewshed computation.

This function samples numSamples random receivers on the ground for every origin location,
instead of computing the actual viewshed from that location. Resulting visibility indices are
therefore approximations, but there is a reasonable trade-off between results and speed-up.

The result is reported as the percentage of sampled locations.

\param numSamples Number of rays to be sampled per viewshed location.
\param outgoing Set to true/false to compute outgoing or incoming visibility.
\param offset Height of the observer, to be added to the ground height at each viewshed location.
*/
ScalarField2 HeightField::ViewshedTotalSampled(int numSamples, bool outgoing, double offset) const
{
    double tolerance = Norm(cellSize);
    Box3 bbox = getBox();
    double invSamples = 1.0/numSamples;
    double klip = LipschitzK();

    ScalarField2 view(domain, nx, ny, 0);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            Vector3 viewpoint = Vertex(i, j) + Vector3(0, 0, offset);

            for (int k = 0; k < numSamples; k++) {
                int vi = int(Random::Uniform()*nx);
                int vj = int(Random::Uniform()*ny);
                Vector3 q = Vertex(vi, vj);
                double dist = Norm(q - viewpoint);
                Vector3 dir = (q - viewpoint) / dist;

                Ray ray(viewpoint, dir);
                double t;
                Vector3 pintersect;
                Intersect(ray, t, pintersect, bbox, klip);

                if (t > dist - tolerance) {
                    if (outgoing) view(i, j) += invSamples;
                    else view(vi, vj) += invSamples;
                }
            }
        }
    }
    return view;
}



/*!
\brief Terrain openness measure as proposed by [Yokoyama et al. 2002].

Given a terrain position and a ray direction, positive openness is the maximum angle with +Z
that a ray can have such that it is tangent but does not intersect the terrain.
Similarly, negative openness is the maximum angle with -Z.

\param positive True for positive openness, false for negative openness.
\param maxDist Maximum distance for the ray traversal. Negative value for infinite.
\param numDirs Number of sampled directions around the circle.
*/
ScalarField2 HeightField::Openness(bool positive, double maxDist, int numDirs) const
{
  // initialize direction vectors
  std::vector<Vector2> dirs;
  for (int i = 0; i < numDirs; i++)
  {
    double a = 2.0 * Math::Pi * i / double(numDirs);
    dirs.push_back(Vector2(std::cos(a), std::sin(a)));
  }
  const double dt = cellSize[0];

  // for each grid cell
  ScalarField2 openness(domain, nx, ny);
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      // start at p
      Vector2 p0 = domainCoords(i, j);
      double  h0 = at(i, j);

      // result
      double sum = 0;

      // for each direction find openness
      for (int k = 0; k < numDirs; k++)
      {
        double t = dt;
        Vector2 p = p0 + t * dirs[k];
        double tanmax = -1e12;
        double tanmin = 1e12;

        // traverse terrain
        while ((maxDist < 0 || t < maxDist) && domain.isInside(p)) {
          double dh = value(p) - h0;
          double tangent = dh / t;
          tanmax = std::max(tanmax, tangent);
          tanmin = std::min(tanmin, tangent);

          t += dt;
          p = p0 + t * dirs[k];
        }

        if (positive) sum += 0.5 * Math::Pi - std::atan(tanmax);
        else          sum += 0.5 * Math::Pi + std::atan(tanmin);
      }

      // average over all dirs
      openness(i, j) = sum / numDirs;
    }
  }

  return openness;
}



/*!
\brief Compute the direct lighting.

\param u Light direction.
\param cosine Boolean, set to true to use scalar product between normal and light, false to use (1+c)/2.
*/
ScalarField2 HeightField::DirectLight(const Vector3 &u, bool cosine) const
{
    ScalarField2 sf(domain, nx, ny, 1.0);

    if (cosine)
    {
      for (int i = 0; i < nx; i++)
      {
        for (int j = 0; j < ny; j++)
        {
          double light = Normal(i, j) * u;
          light = std::max(0.0, light);
          sf(i, j) = light;
        }
      }
    }
    else
    {
      for (int i = 0; i < nx; i++)
      {
        for (int j = 0; j < ny; j++)
        {
          double light = Normal(i, j) * u;
          light = 0.5 * (1.0 + light);
          sf(i, j) = light;
        }
      }

    }

    return sf;
}


/*!
\brief Compute the self shadowing map.
\param light Light direction.
*/
ScalarField2 HeightField::SelfShadow(const Vector3& light) const
{
  ScalarField2 e(domain, nx, ny, 0);

  Box3 box = getBox();
  Vector3 u = Normalized(light);
  const double epsilon = 0.1;
  const double length = Norm(box.size());

  double k = LipschitzK();
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      Ray ray(Vertex(i, j) + u * epsilon, u);
      Vector3 p;
      double t;
      if (!Intersect(ray, t, p, box, k, length))
      {
        e(i, j) = 1;
      }
    }
  }
  return e;
}


/*!
\brief Compute the average direct lighting from the sun over a year, integrating several instants.

\param latitude Latitude of the terrain, longitude is set to 0.0.
\param daystep step in days used in the year summation loop
\param hourstep step in hours used in the day summation loop
\param dayoffset initial day in [0..365] for year loop
\param houroffset initial hour in [0..24] for the day loop
*/
ScalarField2 HeightField::Sun(const double& latitude, int daystep, int hourstep, int dayoffset, int houroffset) const
{
  ScalarField2 sun(domain, nx, ny);

  const double& longitude = 0.0;

  int n = 0;
  for (int i = dayoffset; i < 365; i += daystep)
  {
    for (int ho = houroffset; ho < 24; ho += hourstep)
    {
      double azim;
      double altit;
      Sun::Angles(latitude, longitude, 2000.0, i / 30.0, double(i % 30), double(ho), azim, altit);
      azim = 0.5 * Math::Pi - azim; // convert to CCW angle starting at 0 in east direction

      const Vector2 east(1.0, 0.0);
      const Vector2 north(0.0, 1.0);
      Vector2 xy = cos(azim) * east + sin(azim) * north;
      Vector3 light = (xy * cos(altit)).toVector3(sin(altit));
      //if (light[2] < 0.0) continue;
      ScalarField2 li = DirectLight(light, true) * SelfShadow(light);
      sun+=li;
      n++;
    }
  }

  sun*=1.0 / double(n);
  return sun;
}


/*!
\brief Compute the accessibility or the sky-view factor.

\param r Radius.
\param n Number of rays.
\param skyView if true, returns the sky-view factor (hemisphere oriented to sky). Otherwise, returns accessibility (hemisphere around terrain normal)
*/
ScalarField2 HeightField::Accessibility(const double& r, int n, bool skyView) const
{
  Random::rng.seed(1337);

  double epsilon = 0.05;
  double lipschitz = LipschitzK();
  Box3 box = getBox();

  ScalarField2 sf(domain, nx, ny, 1.0);

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      Vector3 p = Vertex(i, j) + Vector3(0.0, 0.0, epsilon);
      Vector3 normal = skyView ? Vector3(0, 0, 1) : Normal(i, j);

      int hit = 0;
      for (int k = 0; k < n; k++) {
        Vector3 direction = Random::SampleOnHemisphere(normal);
        Ray ray(p, direction);
        double t;
        Vector3 q;
        if (Intersect(ray, t, q, box, lipschitz, r, epsilon / 2.0)) {
          hit++;
        }
      }
      sf(i, j) = 1.0 - double(hit) / double(n);
    }
  }

  return sf;
}


/*!
\brief Computes the analytical approximation to the sky-view factor.

Uses the formula 1 + cos(atan(slope))/2.0,
see Dozier and Frew 1990: "Rapid calculation of terrain parameters for
radiation modeling from digital elevation data".
*/
ScalarField2 HeightField::SkyViewFactorApproximation() const
{
    ScalarField2 svf = Slope();
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            svf(i, j) = (1.0 + Math::CosAtan(svf.at(i, j))) / 2.0;
        }
    }
    return svf;
}


/*!
\brief Computes the Diurnal Anisotropic Heat Index.

See Böhner and Antonic 2009: "Land-surface parameters specific to topo-climatology".

\param maxAspect orientation (in radian) that has maximum heat index, e.g: 135 degree (SW).
*/
ScalarField2 HeightField::DiurnalAnisotropicHeatIndex(double maxAspect) const
{
    ScalarField2 slope = Slope();
    ScalarField2 aspect = Aspect();
    ScalarField2 dahi(domain, nx, ny, 0);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            dahi(i, j) = cos(maxAspect - aspect.at(i, j)) * atan(atan(slope.at(i, j))); // note: atan(atan()) ?
        }
    }
    return dahi;
}
