#include "riversnet.h"
#include <queue>


TerrainFlowD8::TerrainFlowD8(const HeightField &h) : hf(h)
{
    flowdir = hf.FlowDirectionD8();
    upArea = hf.StreamAreaD8();
}

TerrainFlowD8::~TerrainFlowD8()
{
}


Index2 TerrainFlowD8::downstream(int i, int j) const
{
    return downstream(Index2(i, j));
}

Index2 TerrainFlowD8::downstream(const Index2& p) const
{
    if (flowdir.at(p) < 0) return p;
    return flowdir.neighborCell(p, flowdir.at(p));
}

std::vector<Index2> TerrainFlowD8::upstream(int i, int j) const
{
    return upstream(Index2(i,j));
}

std::vector<Index2> TerrainFlowD8::upstream(const Index2& p) const
{
    std::vector<Index2> ups;
    for (int i = 0; i < 8; i++) {
        Index2 q = flowdir.neighborCell(p, i);
        if (flowdir.isValidCell(q.x(), q.y()) && flowdir.at(q) == neighborOutdirIntoCell[i]) {
            ups.push_back(q);
        }
    }
    return ups;
}


double TerrainFlowD8::flowConfluenceDistance(const Index2 &p, const Index2 &q) const
{
    Index2 pflow = p;
    Index2 qflow = q;
    double pdist = 0;
    bool qflowEnded = false;
    double hp = hf.at(p);
    double hq = hf.at(q);

    while (pflow != qflow) {

        // P flow above flow from Q, advance path from P
        while (hp >= hq || qflowEnded) {

            // reached a confluence with flow from Q? return traversed distance
            if (pflow == qflow) return pdist;

            // nowhere to go? end of flow from p whether we reached q path or not
            int dirCode = flowdir.at(pflow);
            if (dirCode < 0) return pdist;

            // advance flow from P
            pflow += FieldGrid2D::next[dirCode];
            pdist += FieldGrid2D::length[dirCode];

            // if we leave domain, return
            if (!hf.isValidCell(pflow)) return pdist;
            hp = hf.at(pflow);
        }

        // now advance flow path from Q
        while (!qflowEnded && hq > hp) {
            int dirCode = flowdir.at(qflow);

            // if flow from Q stops before getting lower than P flow, P flow will not reach it
            if (dirCode < 0) qflowEnded = true;
            qflow += FieldGrid2D::next[dirCode];

            // check inside domain
            if (!hf.isValidCell(qflow)) qflowEnded = true;
            else hq = hf.at(qflow);
        }
    }

    return pdist;
}

IntField2 TerrainFlowD8::riverCellsFromChannelThreshold(double cit_s, double cit_t, bool propagate) const
{
    // mark as river cells those that meet the CIT
    IntField2 riverCells(hf.getDomain(), hf.getSizeX(), hf.getSizeY(), 0);
    ScalarField2 s = hf.Slope();
    for (int i = 0; i < riverCells.getSizeX(); i++) {
        for (int j = 0; j < riverCells.getSizeY(); j++) {
            // A * S^s >= t?
            double c = upArea.at(i, j)*std::pow(s(i, j), cit_s);
            if (c >= cit_t) {
                riverCells(i, j) = 1;
            }
        }
    }

    if (propagate) {
        // ensure downstream cells from a river cell are marked as river too
        std::vector<std::pair<double, Index2> > pts = hf.valuesWithIndex();
        std::sort(pts.begin(), pts.end());
        for (int i = int(pts.size())-1; i >= 0; i--) {
            Index2 p = pts[i].second;
            if (riverCells.at(p)) {
                riverCells(downstream(p)) = 1; // downstream returns p if no flow down
            }
        }
    }

    return riverCells;
}

RiverTree *TerrainFlowD8::buildRiverTree(const Index2 &endpoint, const IntField2 &mask) const
{
    RiverTree* river = new RiverTree();

    Index2 p = endpoint;
    bool followUpstream = true;

    while (followUpstream) {

        // add current point to river
        river->cells.push_back(p);

        // get upstream neighbors accounting for mask
        std::vector<Index2> upNodes;
        for (const Index2& q : upstream(p)) {
            if (mask.at(q) > 0) {
                upNodes.push_back(q);
            }
        }

        // source?
        if (upNodes.size() == 0) {
            followUpstream = false;
        }
        // junction?
        else if (upNodes.size() > 1) {
            followUpstream = false;
            // recursively create parent channels
            for (const Index2& q : upNodes) {
                RiverTree* parent = buildRiverTree(q, mask);
                parent->cells.insert(parent->cells.begin(), p);
                river->parents.push_back(parent);
            }
        }
        // just one node upstream, still in river segment
        else {
            p = upNodes[0];
        }
    }

    return river;
}

DrainageBasin TerrainFlowD8::getDrainageBasin(const Index2 &point) const
{
    IntField2 mask(hf.getDomain(), hf.getSizeX(), hf.getSizeY(), 0);

    std::queue<Index2> Q;
    Q.push(point);
    while (!Q.empty()) {
        Index2 p = Q.front();
        Q.pop();
        mask(p) = 1;
        for (const Index2& q : upstream(p)) {
            Q.push(q);
        }
    }

    return DrainageBasin(mask);
}


void RiverTree::computeRiverMetrics(const HeightField &hf, const ScalarField2 &sa)
{
    // recursive call to parents first
    for (RiverTree* r : parents) {
        r->computeRiverMetrics(hf, sa);
    }

    // strahler index
    if (parents.size() > 0) {
        int smin = parents[0]->strahlerIndex;
        int smax = smin;
        for (unsigned int i = 1; i < parents.size(); i++) {
            smin = std::min(smin, parents[i]->strahlerIndex);
            smax = std::max(smax, parents[i]->strahlerIndex);
        }
        strahlerIndex = std::max(smin + 1, smax);
    }

    // profiles of stream area, flown distance, slopes
    nodeStreamArea.resize(cells.size());
    nodeAccumLength.resize(cells.size());
    nodeSlope.resize(cells.size());

    Vector3 p = hf.Vertex(cells[0].x(), cells[0].y());
    nodeStreamArea[0] = sa.at(cells[0]);
    nodeAccumLength[0] = 0;
    nodeSlope[0] = 0;

    if (cells.size() < 2) return;

    for (unsigned int i = 1; i < cells.size(); i++) {
        Vector3 q = hf.Vertex(cells[i].x(), cells[i].y());
        Vector3 vpq = q - p;

        double d = std::sqrt(vpq[0]*vpq[0] + vpq[1]*vpq[1]);
        nodeStreamArea[i] = sa.at(cells[i]);
        nodeAccumLength[i] = nodeAccumLength[i-1] + d;
        nodeSlope[i] = d > 0 ? vpq[2]/d : 0;

        p = q;
    }

    // sinuosity
    double distStraight = Norm(hf.domainCoords(cells.back()) - hf.domainCoords(cells[0]));
    sinuosity = distStraight > 0 ? nodeAccumLength.back()/distStraight : 1;
}

void RiverTree::markRiverCells(IntField2 &mask, int value) const
{
    for (const Index2& p : cells) {
        mask(p.x(), p.y()) = value;
    }
    for (const RiverTree* r : parents) {
        r->markRiverCells(mask);
    }
}

void RiverTree::markRiverLength(ScalarField2 &rivLength) const
{
    for (unsigned int i = 0; i < cells.size(); i++) {
        // at least 0.001 so we have some length reported at the source cells
        rivLength(cells[i].x(), cells[i].y()) = std::max(nodeAccumLength[i], 0.001);
    }
    for (const RiverTree* r : parents) {
        r->markRiverLength(rivLength);
    }
}

void RiverTree::markRiverStreamArea(ScalarField2 &streamArea) const
{
    for (unsigned int i = 0; i < cells.size(); i++) {
        streamArea(cells[i].x(), cells[i].y()) = nodeAccumLength[i];
    }
    for (const RiverTree* r : parents) {
        r->markRiverStreamArea(streamArea);
    }
}


bool DrainageBasin::inPerimeter(int i, int j) const
{
    if (!basinMask.isValidCell(i, j) || !inBasin(i, j)) return false;
    if (!basinMask.isValidCell(i-1, j) || !inBasin(i-1, j)) return true;
    if (!basinMask.isValidCell(i+1, j) || !inBasin(i+1, j)) return true;
    if (!basinMask.isValidCell(i, j-1) || !inBasin(i, j-1)) return true;
    if (!basinMask.isValidCell(i, j+1) || !inBasin(i, j+1)) return true;
    return false;
}

void DrainageBasin::computeBasinMetrics(const HeightField& hf, const ScalarField2& riverLength, const ScalarField2& riverStreamArea)
{
    countInner = 0;
    countPerim = 0;
    countRiver = 0;
    maxStreamLength = 0;
    maxStreamDrainage = 0;

    hf.getRange(hmax, hmin);
    double hsum = 0;

    for (int i = 0; i < basinMask.getSizeX(); i++) {
        for (int j = 0; j < basinMask.getSizeY(); j++) {
            if (inBasin(i, j)) {
                countInner++;
                if (inPerimeter(i, j)) countPerim++;
                if (riverLength.at(i, j) > 0) {
                    countRiver++;
                    maxStreamLength   = std::max(maxStreamLength, riverLength.at(i, j));
                    maxStreamDrainage = std::max(maxStreamDrainage, riverStreamArea.at(i, j));
                }

                double h = hf.at(i, j);
                hsum += h;
                hmin = std::min(h, hmin);
                hmax = std::max(h, hmax);
            }
        }
    }
    havg = hsum/countInner;
}
