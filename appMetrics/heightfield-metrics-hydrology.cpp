#include "heightfield.h"


HeightField::Flow HeightField::GetFlow(const Index2& p) const
{
    HeightField::Flow flow;

    double zp = at(p);
    for (int i = 0; i < 8; i++) {
        Index2 q = p + next[i];
        if (!isValidCell(q)) continue;

        double step = at(q) - zp;
        if (step < -HeightField::flat) { // Should be 0.0, but very small negative values can cause crashes
            flow.i[flow.n] = i;
            flow.q[flow.n] = q;
            flow.s[flow.n] = -step * inverselength[i];
            flow.slopesum += flow.s[flow.n];

            // Steepest slope
            if (flow.n == 0) {
                flow.steepest = 0;
            }
            else if (flow.s[flow.n] > flow.s[flow.steepest]) {
                flow.steepest = flow.n;
            }
            flow.n++;
        }
    }

    return flow;
}


IntField2 HeightField::FlowDirectionD8() const
{
    IntField2 dir(getDomain(), nx, ny, -1);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            Flow flow = GetFlow(Index2(i,j));
            if (flow.n > 0) {
                dir(i, j) = flow.i[flow.steepest];
            }
        }
    }
    return dir;
}


ScalarField2 initFlowFieldFromSources(const std::vector<Index2>& sources, const Box2& domain, int nx, int ny)
{
    // no sources, flow from everywhere
    if (sources.size() == 0)
        return ScalarField2(domain, nx, ny, 1);

    // otherwise, mark only source points
    ScalarField2 s(domain, nx, ny, 0);
    for (const Index2& ps : sources) {
        s(ps.x(), ps.y()) = 1;
    }
    return s;
}


ScalarField2 HeightField::StreamAreaD8(const std::vector<Index2>& sources) const
{
    ScalarField2 sa = initFlowFieldFromSources(sources, getDomain(), nx, ny);

    std::vector<std::pair<double, Index2> > pts = valuesWithIndex();
    std::sort(pts.begin(), pts.end());
    for (int i = int(pts.size())-1; i >= 0; i--) {
        Index2 p = pts[i].second;
        Flow flow = GetFlow(p);
        if (flow.n > 0) {
            Index2 q = flow.q[flow.steepest];
            sa(q) += sa.at(p);
        }
    }

    return sa;
}


ScalarField2 HeightField::StreamAreaMFD(const std::vector<Index2>& sources, double gamma) const
{
    ScalarField2 sa = initFlowFieldFromSources(sources, getDomain(), nx, ny);

    std::vector<std::pair<double, Index2> > pts = valuesWithIndex();
    std::sort(pts.begin(), pts.end());
    for (int i = int(pts.size())-1; i >= 0; i--) {
        Index2 p = pts[i].second;
        Flow flow = GetFlow(p);
        if (flow.n > 0) {
            double div = 0;
            for (int j = 0; j < flow.n; j++) {
                div += std::pow(flow.s[j], gamma);
            }
            div = 1.0/div;

            for (int j = 0; j < flow.n; j++) {
                double f = div * std::pow(flow.s[j], gamma);
                Index2 q = flow.q[j];
                sa(q) += f*sa.at(p);
            }
        }
    }

    return sa;
}


ScalarField2 HeightField::StreamAreaDinf(const std::vector<Index2>& sources) const
{
    ScalarField2 sa = initFlowFieldFromSources(sources, getDomain(), nx, ny);

    std::vector<std::pair<double, Index2> > pts = valuesWithIndex();
    std::sort(pts.begin(), pts.end());
    std::vector<Flow> flows(pts.size());
    for (unsigned int i = 0; i < flows.size(); i++) {
        flows[i] = GetFlow(pts[i].second);
    }

    for (int i = int(pts.size())-1; i >= 0; i--) {
        Index2 p = pts[i].second;
        Vector2 grad = Gradient(p.x(), p.y());
        double angle = std::atan2(-grad[1], -grad[0]);
        if (angle < 0) angle += 2*Math::Pi;
        int octan = int(4 * angle / Math::Pi);
        double alpha = 4 * angle / Math::Pi - octan;

        Index2 q1 = p + next[octan];
        Index2 q2 = p + next[(octan + 1) % 8];
        double alpha1 = 1 - alpha;
        double alpha2 = alpha;

        if (isValidCell(q1)) sa(q1) += alpha1 * sa(p);
        if (isValidCell(q2)) sa(q2) += alpha2 * sa(p);
    }

    return sa;
}


ScalarField2 HeightField::StreamAreaRho8(const std::vector<Index2>& sources, int niters) const
{
    ScalarField2 sa(getDomain(), nx, ny, 0);

    std::vector<std::pair<double, Index2> > pts = valuesWithIndex();
    std::sort(pts.begin(), pts.end());
    std::vector<Flow> flows(pts.size());
    for (unsigned int i = 0; i < flows.size(); i++) {
        flows[i] = GetFlow(pts[i].second);
    }

    for (int iter = 0; iter < niters; iter++) {

        ScalarField2 si = initFlowFieldFromSources(sources, getDomain(), nx, ny);

        for (int i = int(pts.size())-1; i >= 0; i--) {
            Index2 p = pts[i].second;
            Flow flow = flows[i];
            if (flow.n > 0) {
                // choose one outflow randomly based on slope
                double srandom = Random::Uniform(0, flow.slopesum);
                double ssum = 0;
                int j = 0;
                while (j < flow.n) {
                    ssum += flow.s[j];
                    if (ssum < srandom) j++;
                    else break;
                }
                // flow there
                Index2 q = flow.q[j];
                si(q) += si.at(p);
            }
        }

        // accumulate all maps
        sa += si;
    }

    return sa;
}


ScalarField2 HeightField::StreamAreaKinematic(const std::vector<Index2>& sources) const
{
    ScalarField2 sa(getDomain(), nx, ny, 0);

    // precompute gradients
    std::vector<Vector2> gradients(getNumElements());
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            gradients[cellId(i, j)] = -Normalized(Gradient(i, j));
        }
    }

    // source points
    std::vector<Index2> sourcePts = sources;
    if (sourcePts.size() == 0) {
        sourcePts.resize(nx*ny);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                sourcePts[cellId(i, j)] = Index2(i,j);
            }
        }
    }

    // compute flow trajectory per source point
    double dt = std::min(cellSize[0], cellSize[1]);
    ScalarField2 si(getDomain(), nx, ny);
    for (const Index2& ps : sourcePts) {

        // clear trajectory buffer
        si.fill(0);

        // init a "rolling ball" at cell
        Index2 cell = ps;
        double hprev = at(cell);
        Vector2 p = domainCoords(cell.x(), cell.y());

        while (isValidCell(cell)) {
            // do not go upstream
            if (hprev < at(cell)) break;
            hprev = at(cell);

            // prevent flow loops
            if (si(cell) > 0) break;
            si(cell) += 1;

            // simplification: use gradient instead of fitting a plane
            Vector2 flowdir = gradients[cellId(cell.x(), cell.y())];

            // flow until we change cell
            Index2 currCell = cell;
            while (currCell == cell) {
                p = p + dt*flowdir;
                cell = cellCoords(p);
            }
        }

        // accumulate
        sa += si;
    }

    return sa;
}


/*!
\brief Compute the wetness index field of the terrain.

Note that the wetness index is defined as ln ( A / s ) where
s is the slope, however the definition does not work well
for small slopes, therefore the implementation defines the wetness
index as ln ( A / ( e + s ) ), where e is 0.001.

\param sa  Stream Area
*/
ScalarField2 HeightField::WetnessIndex(const ScalarField2& sa) const
{
    ScalarField2 wetness(domain, nx, ny, 0);
    ScalarField2 slope  = Slope();
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            wetness(i, j) = std::log(sa.at(i, j) / (1e-3 + slope.at(i,j)));
        }
    }
    return wetness;
}

ScalarField2 HeightField::WetnessIndex() const
{
    return WetnessIndex(StreamAreaMFD());
}


/*!
\brief Compute the stream power field of the terrain.

\param sa  Stream Area
\param m,n Exponent for the stream area and slope, n should be twice as m.
*/
ScalarField2 HeightField::StreamPower(const ScalarField2& sa, double m, double n) const
{
    ScalarField2 power(domain, nx, ny, 0);
    ScalarField2 slope = Slope();
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            power(i, j) = pow(sa.at(i, j), m) * pow(slope.at(i, j), n);
        }
    }
    return power;
}

/*!
\brief Compute the stream power field of the terrain.
This function does not require to pass a Stream Area, as it will compute it
using Multiple Flow Direction (with gamma=1). For other SA methods, or to
reuse an already computed stream area Scalarfield, use the previous function.

\param m,n Exponent for the stream area and slope, n should be twice as m.
*/
ScalarField2 HeightField::StreamPower(double m, double n) const
{
    return StreamPower(StreamAreaMFD(), m, n);
}
