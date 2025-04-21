#include "heightfield.h"
#include <queue>


/*!
\brief Compute the peakedness of a terrain, defined as the percentage of cells in a disk that have lower or equal elevation
than the center position.

See Llobera 2001, Building Past Landscape Perception With GIS: Understanding Topographic Prominence.
Although the author calls it <I>prominence</I>, it is not the same definition as used in orometry
for peak prominence.

See Arge et al. 2013: Algorithms for Computing Prominence on Grid Terrains for faster implementations.

\param r Radius of the disk (we actually use a squared window).

\author Oscar Argudo
*/
ScalarField2 HeightField::PeakPercentage(double r) const
{
    ScalarField2 res(domain, nx, ny, 0.0);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double h0 = at(i, j);
            int neighs = 0;
            int lowers = 0;
            IndexArea area = indexAreaFromRadius(i, j, r);
            for (int di = area.xmin(); di <= area.xmax(); di++) {
                for (int dj = area.ymin(); dj <= area.ymax(); dj++) {
                    neighs++;
                    lowers += at(di, dj) < h0 ? 1 : 0;
                }
            }
            res(i, j) = double(lowers) / neighs;
        }
    }
    return res;
}


/*!
\brief Compute the omni-directional relief and steepness.

See Earl and Metzler 2015, Cloud-Capped Towers: Capturing Terrain Characteristics Using Topographic Functionals.

From a theoretical point of view, ORS is the value that results from convergence when we integrate all distances.
ORS monotonically increases as we make this cut-off bigger, since the function we integrate is always positive.
From my experience (and by definition) after only a few km of distance 10-20, the value is already a good approximation
but also, different radius cut-off highlight different interesting properties. If we use a small radius, those points with a very
steep edge nearby or even needle-shaped secondary small peaks will have large values of ORS, more than other bigger but gentler peaks.
When the radius increases, bigger peaks shadow them with larger ORS values.

\param r Radius for cutting the computation at some point.

\author Oscar Argudo
*/
ScalarField2 HeightField::ORS(double r) const
{
  ScalarField2 ors(domain, nx, ny, 0.0);

  double cellW = getCellSize()[0];
  double cellH = getCellSize()[1];
  double cellArea = cellW*cellH;

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      double h0 = at(i, j);
      double value = 0;

      IndexArea area = indexAreaFromRadius(i, j, r);
      for (int di = area.xmin(); di <= area.xmax(); di++)
      {
        for (int dj = area.ymin(); dj <= area.ymax(); dj++)
        {
          double dist = sqrt((di - i) * cellW * (di - i) * cellW +
                             (dj - j) * cellH * (dj - j) * cellH);
          if (dist > r) continue;

          double hij = at(di, dj);

          // points higher than (i,j) do not contribute
          if (hij < h0)
          {
            // slope
            double s = (h0 - hij) / dist;
            // normalized slope, for convergence. Other functions could be used, see original paper
            double atanu = atan(s);
            double f = (4.0 / (Math::Pi * Math::Pi * Math::Pi)) * (2 * s * atanu - log(s * s + 1) - atanu * atanu);
            // integrate. We clamp to 0 to avoid precision errors and negative f
            value += std::max(0.0, f * cellArea);
          }
        }
      }

      ors(i, j) = sqrt(value);
    }
  }
  return ors;
}



/*!
\brief Compute the Jut metric.

See Xu 2022, "Datumless Topography, A Universally Consistent Way to Quantify Relief".

It is defined as the maximum angle-reduced height of p with respect to any other point q.
Jut reaches larger values- more impressiveness- when the terrain rises abruptly in a small area,
such as in the vicinity of cliffs.
Note that this metric is intended to be used taking into account planetary curvature, here
we coded a simplified version assuming the height field domain is a plane.

\param r Radius for cutting the computation at some point.

\author Oscar Argudo
*/
ScalarField2 HeightField::JutPlanar(double r) const
{
    ScalarField2 jut(domain, nx, ny, 0.0);

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            Vector3 p = Vertex(i, j);
            double jmax = 0;

            IndexArea area = indexAreaFromRadius(i, j, r);
            for (int di = area.xmin(); di <= area.xmax(); di++)
            {
                for (int dj = area.ymin(); dj <= area.ymax(); dj++)
                {
                    if (di == i && dj == j) continue;
                    double dist = sqrt((di - i) * cellSize[0] * (di - i) * cellSize[0] +
                                       (dj - j) * cellSize[1] * (dj - j) * cellSize[1]);
                    if (dist > r) continue;

                    Vector3 q = Vertex(di, dj);
                    double z = p[2] - q[2];
                    double d = Norm(p - q);

                    // z' = z * sin(a)
                    // sin(a) = |z|/dist_3d(p, q)
                    double zreduced = z * abs(z) / d;
                    jmax = std::max(jmax, zreduced);
                }
            }

            jut(i, j) = jmax;
        }
    }
    return jut;
}


/*!
\brief Compute the Rut metric.

See Xu 2022, "Datumless Topography, A Universally Consistent Way to Quantify Relief".

It is defined as the maximum angle-reduced height of any point q with respect to p.
Rut measures how sharply or impressively a point dips below its immediate surroundings.
Note that this metric is intended to be used taking into account planetary curvature, here
we coded a simplified version assuming the height field domain is a plane.

\param r Radius for cutting the computation at some point.

\author Oscar Argudo
*/
ScalarField2 HeightField::RutPlanar(double r) const
{
    ScalarField2 rut(domain, nx, ny, 0.0);

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            Vector3 p = Vertex(i, j);
            double rmax = 0;

            IndexArea area = indexAreaFromRadius(i, j, r);
            for (int di = area.xmin(); di <= area.xmax(); di++)
            {
                for (int dj = area.ymin(); dj <= area.ymax(); dj++)
                {
                    if (di == i && dj == j) continue;
                    double dist = sqrt((di - i) * cellSize[0] * (di - i) * cellSize[0] +
                                       (dj - j) * cellSize[1] * (dj - j) * cellSize[1]);
                    if (dist > r) continue;

                    Vector3 q = Vertex(di, dj);
                    double z = q[2] - p[2];
                    double d = Norm(p - q);

                    // z' = z * sin(a)
                    // sin(a) = |z|/dist_3d(p, q)
                    double zreduced = z * abs(z) / d;
                    rmax = std::max(rmax, zreduced);
                }
            }

            rut(i, j) = rmax;
        }
    }
    return rut;
}


double zReduced(double hp, double hq, double planeDist, double planetR, bool converging)
{
    // angle span between the two points on planet surface
    // we assume the distance "on map" is actually the arc length
    double angle = planeDist/planetR;
    Vector2 k_q(std::sin(angle), std::cos(angle));
    Vector2 p = Vector2(0, planetR + hp);
    Vector2 q = (planetR + hq)*k_q;
    Vector2 r_pq = converging ? p - q : q - p;
    double d = Norm(r_pq);
    double z_pq = r_pq[1]; // dot(r_pq, k_p) with kp=(0,0,1);
    double zreduced = z_pq * std::abs(z_pq)/d;  // std::abs(dot(r_pq/d, k_q));
    return zreduced;
}


/*!
\brief Compute the Jut metric.

See Xu 2022, "Datumless Topography, A Universally Consistent Way to Quantify Relief".

It is defined as the maximum angle-reduced height of p with respect to any other point q.
Jut reaches larger values- more impressiveness- when the terrain rises abruptly in a small area,
such as in the vicinity of cliffs.

\param r Radius for cutting the computation at some point,
although planetR limits the domain of jut > 0 for a point.

\param planetR Planet radius, in meters.

\author Oscar Argudo
*/
ScalarField2 HeightField::JutCurved(double r, double planetR) const
{
    ScalarField2 jut(domain, nx, ny, 0.0);

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            double jmax = 0;
            double hp = at(i, j);

            IndexArea area = indexAreaFromRadius(i, j, r);
            for (int di = area.xmin(); di <= area.xmax(); di++)
            {
                for (int dj = area.ymin(); dj <= area.ymax(); dj++)
                {
                    if (di == i && dj == j) continue;

                    double dist = sqrt((di - i) * cellSize[0] * (di - i) * cellSize[0] +
                                       (dj - j) * cellSize[1] * (dj - j) * cellSize[1]);
                    if (dist > r) continue;

                    double zreduced = zReduced(hp, at(di, dj), dist, planetR, true);
                    jmax = std::max(jmax, zreduced);
                }
            }

            jut(i, j) = jmax;
        }
    }
    return jut;
}


/*!
\brief Compute the Rut metric.

See Xu 2022, "Datumless Topography, A Universally Consistent Way to Quantify Relief".

It is defined as the maximum angle-reduced height of any point q with respect to p.
Rut measures how sharply or impressively a point dips below its immediate surroundings.
Note that this metric is intended to be used taking into account planetary curvature, here
we coded a simplified version assuming the height field domain is a plane.

\param r Radius for cutting the computation at some point.

\param planetR Planet radius, in meters.

\author Oscar Argudo
*/
ScalarField2 HeightField::RutCurved(double r, double planetR) const
{
    ScalarField2 rut(domain, nx, ny, 0.0);

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            double rmax = 0;
            double hp = at(i, j);

            IndexArea area = indexAreaFromRadius(i, j, r);
            for (int di = area.xmin(); di <= area.xmax(); di++)
            {
                for (int dj = area.ymin(); dj <= area.ymax(); dj++)
                {
                    if (di == i && dj == j) continue;

                    double dist = sqrt((di - i) * cellSize[0] * (di - i) * cellSize[0] +
                                       (dj - j) * cellSize[1] * (dj - j) * cellSize[1]);
                    if (dist > r) continue;

                    double zreduced = zReduced(hp, at(di, dj), dist, planetR, false);
                    rmax = std::max(rmax, zreduced);
                }
            }

            rut(i, j) = rmax;
        }
    }
    return rut;
}

ScalarField2 HeightField::AngleReducedHeight(double R, int i, int j, bool converging) const
{
    ScalarField2 a(domain, nx, ny, 0);

    double hp = at(i, j);
    for (int di = 0; di < nx; di++) {
        for (int dj = 0; dj < ny; dj++) {
            double dist = sqrt((di - i) * cellSize[0] * (di - i) * cellSize[0] +
                               (dj - j) * cellSize[1] * (dj - j) * cellSize[1]);
            if (dist == 0) continue;
            a(di, dj) = zReduced(hp, at(di, dj), dist, R, converging);
        }
    }

    return a;
}


/*!
\brief Compute the Peakness metric and the peak prototypicality meatric.

See Deng and Wilson 2008, "Multi‐scale and multi‐criteria mapping of mountain peaks as fuzzy entities".

It is defined as the maximum angle-reduced height of any point q with respect to p.
Rut measures how sharply or impressively a point dips below its immediate surroundings.
Note that this metric is intended to be used taking into account planetary curvature, here
we coded a simplified version assuming the height field domain is a plane.

\param prototypicality Output field for prototypicality metric
\param centroidPercentile

\author Oscar Argudo
*/
ScalarField2 HeightField::Peakness(ScalarField2& prototypicality, double centroidPercentile) const
{
    const int NUM_RADII = 3;
    const int RADII[NUM_RADII] = { 5, 20, 50 };
    const double PEAK_THRES[NUM_RADII] = { 20, 80, 80 };
    const double WEIGHT[NUM_RADII] = { 1.0, 1.2, 1.4 };
    const double WEIGHT_SUM = 3.6;
    const int EXTENDED_RADII_FACTOR = 3;
    const int NUM_METRICS = 4;

    ScalarField2 satSlopes = Slope().summedAreaTable();

    prototypicality = ScalarField2(domain, nx, ny, 0.0);
    std::vector<std::vector<int> > peakLocations = std::vector<std::vector<int> >(NUM_RADII);
    std::vector<std::vector<double> > peakScaleProto = std::vector<std::vector<double> >(NUM_RADII);
    std::vector<std::vector<ScalarField2> > scaleMetrics(NUM_RADII, std::vector<ScalarField2>(NUM_METRICS));

    // identify peaks at different scales, compute their properties, and obtain prototypicality
    for (int ri = 0; ri < NUM_RADII; ri++) {

        int r = RADII[ri];
        ScalarField2 radMin = minFilter(r);
        ScalarField2 radMax = maxFilter(r);
        ScalarField2 extMin = minFilter(EXTENDED_RADII_FACTOR * r);
        ScalarField2 extMax = maxFilter(EXTENDED_RADII_FACTOR * r);
        ScalarField2& pRelief = scaleMetrics[ri][0];
        ScalarField2& pSlope = scaleMetrics[ri][1];
        ScalarField2& pRelElev = scaleMetrics[ri][2];
        ScalarField2& pIsolation = scaleMetrics[ri][3];

        // find peaks and compute their four factors
        ScalarField2 isPeak(domain, nx, ny, 0);
        pRelief = ScalarField2(domain, nx, ny, 0);
        pSlope = ScalarField2(domain, nx, ny, 0);
        pRelElev = ScalarField2(domain, nx, ny, 0);
        pIsolation = ScalarField2(domain, nx, ny, 0);

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                int idx = cellId(i, j);
                int imin = std::max(i - r, 0);
                int jmin = std::max(j - r, 0);
                int imax = std::min(i + r, nx - 1);
                int jmax = std::min(j + r, ny - 1);

                // peak if cell has maximum elevation in radius and elevation drop in it is sufficiently large
                isPeak[idx] = radMax.at(idx) <= at(idx) && radMax.at(idx) - radMin.at(idx) > PEAK_THRES[ri] ? 1.0 : 0.0;
                if (isPeak[idx]) {
                    peakLocations[ri].push_back(idx);
                }

                // 1) relief is the difference between minimum and maximum elevation in window
                pRelief[idx] = radMax.at(idx) - radMin.at(idx);

                // 2) slope is the mean slope in window
                pSlope[idx] = satSlopes.at(imin, jmin) + satSlopes.at(imax, jmax)
                  - satSlopes.at(imin, jmax) - satSlopes.at(imax, jmin);
                pSlope[idx] /= (imax - imin) * (jmax - jmin);

                // 3) relative elevation in extended radius window
                pRelElev[idx] = (at(idx) - extMin.at(idx)) / (extMax.at(idx) - extMin.at(idx));
            }
        }

        ScalarField2 satNeighs = isPeak.summedAreaTable();
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                int idx = cellId(i, j);
                int imin = std::max(i - EXTENDED_RADII_FACTOR * r, 0);
                int jmin = std::max(j - EXTENDED_RADII_FACTOR * r, 0);
                int imax = std::min(i + EXTENDED_RADII_FACTOR * r, nx - 1);
                int jmax = std::min(j + EXTENDED_RADII_FACTOR * r, ny - 1);

                // 4) number of peaks in extended radius window, inverted
                pIsolation[idx] = satNeighs.at(imin, jmin) + satNeighs.at(imax, jmax)
                                - satNeighs.at(imin, jmax) - satNeighs.at(imax, jmin);
                pIsolation[idx] = 1.0 / (1.0 + pIsolation[idx]);
            }
        }

        // compute peak prototypicality as the simple average of the four normalized factors
        ScalarField2 proto(domain, nx, ny, 0);
        for (const ScalarField2& s : scaleMetrics[ri]) {
            ScalarField2 normMetric(s);
            normMetric.normalize();
            proto += normMetric;
        }
        proto *= 1.0 / NUM_METRICS;

        // sum peak proto across scales
        proto *= WEIGHT[ri];
        prototypicality += proto;

        // add measure
        for (int ploc : peakLocations[ri]) {
            peakScaleProto[ri].push_back(proto.at(ploc));
        }
    }
    prototypicality *= 1.0 / WEIGHT_SUM;


    // find the topmost peak representatives at each scale
    ScalarField2 peakSimilarity(domain, nx, ny, 0);;

    for (int ri = 0; ri < NUM_RADII; ri++) {

        // get top peaks
        std::priority_queue<std::pair<double, int>> peakProtos;
        for (unsigned int i = 0; i < peakLocations[ri].size(); i++) {
            peakProtos.push(std::make_pair(peakScaleProto[ri][i], peakLocations[ri][i]));
        }

        // compute peak class centroid
        int nCentroid = std::max(1, int(peakLocations[ri].size() * centroidPercentile));
        std::vector<double> peakClassCentre(NUM_METRICS, 0);
        for (int i = 0; i < nCentroid; i++) {
            int loc = peakProtos.top().second;
            peakProtos.pop();
            for (int j = 0; j < NUM_METRICS; j++) {
                peakClassCentre[j] += scaleMetrics[ri][j][loc];
            }
        }
        for (int j = 0; j < NUM_METRICS; j++) {
            peakClassCentre[j] /= nCentroid;
        }

        // variances of each metric
        std::vector<double> metricSum(NUM_METRICS, 0);
        std::vector<double> metricSSum(NUM_METRICS, 0);
        for (int j = 0; j < NUM_METRICS; j++) {
            const ScalarField2& metric = scaleMetrics[ri][j];
            for (int loc : peakLocations[ri]) {
                metricSum[j] += metric.at(loc);
                metricSSum[j] += (metric.at(loc) * metric.at(loc));
            }
        }
        std::vector<double> metricMean(NUM_METRICS, 0);
        std::vector<double> metricVariance(NUM_METRICS, 0);
        for (int j = 0; j < NUM_METRICS; j++) {
            int n = int(peakLocations[ri].size());
            double mean = metricSum[j] / n;
            metricMean[j] = mean;
            metricVariance[j] = (metricSSum[j] / n) - (mean * mean);
        }

        // compute this scale similarity map using diagonal distance to centroid
        ScalarField2 distance(domain, nx, ny, 0);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                double d = 0;
                for (int k = 0; k < NUM_METRICS; k++) {
                    double dif = scaleMetrics[ri][k].at(i, j) - peakClassCentre[k];
                    d += (dif * dif) / metricVariance[k];
                }
                distance(i, j) = std::sqrt(d);
            }
        }

        double minDist, maxDist;
        distance.getRange(minDist, maxDist);
        ScalarField2 simil(domain, nx, ny);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                simil(i, j) = (maxDist - distance.at(i, j)) / (maxDist - minDist);
            }
        }

        // accumulate across scales
        simil *= WEIGHT[ri];
        peakSimilarity += simil;
    }

    peakSimilarity *= 1.0 / WEIGHT_SUM;

    return peakSimilarity;
}
