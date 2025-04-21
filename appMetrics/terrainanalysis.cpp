#include "terrainanalysis.h"
#include "tree_builder.h"
#include "line_tree.h"
#include "utm_coordinate_system.h"
#include <iostream>
#include <queue>
#include <map>


TerrainAnalysis::TerrainAnalysis(const HeightField &h) : hf(h)
{
}

TerrainAnalysis::~TerrainAnalysis()
{
    if (terrainD8) delete terrainD8;
    if (ridgesPPA) delete ridgesPPA;
    if (divideTree) delete divideTree;
    if (islandTree) delete islandTree;
    if (prunedDivTree) delete prunedDivTree;
    if (prunedIslandTree) delete prunedIslandTree;
    for (RiverTree* r : riverTrees) delete r;
}

TerrainFlowD8 *TerrainAnalysis::getFlowDirs()
{
    if (!terrainD8) {
        terrainD8 = new TerrainFlowD8(hf);
    }
    return terrainD8;
}

std::vector<RiverTree *> TerrainAnalysis::computeRivers(double cit_s, double cit_t)
{
    if (!terrainD8) {
        terrainD8 = new TerrainFlowD8(hf);
    }

    if (!riversComputed || cit_s != riversCIT_s || cit_t != riversCIT_t) {
        for (RiverTree* r : riverTrees) delete r;
        riverTrees.clear();

        riverCells = terrainD8->riverCellsFromChannelThreshold(cit_s, cit_t, true);
        for (int i = 0; i < riverCells.getSizeX(); i++) {
            for (int j = 0; j < riverCells.getSizeY(); j++) {
                if (riverCells.at(i, j) > 0) {
                    if (!terrainD8->hasDownstream(i, j)) {
                        RiverTree* r = terrainD8->buildRiverTree(Index2(i, j), riverCells);
                        r->computeRiverMetrics(hf, terrainD8->streamArea());
                        riverTrees.push_back(r);
                    }
                }
            }
        }

        riversComputed = true;
        riversCIT_s = cit_s;
        riversCIT_t = cit_t;
    }

    return riverTrees;
}

IntField2 TerrainAnalysis::RiversIdMask() const
{
    IntField2 mask(hf.getDomain(), hf.getSizeX(), hf.getSizeY(), 0);
    for (unsigned int i = 0; i < riverTrees.size(); i++) {
        riverTrees[i]->markRiverCells(mask, i+1);
    }
    return mask;
}

ScalarField2 TerrainAnalysis::DistanceToNearestDrainageEuclidean() const
{
    ScalarField2 distRiver(riverCells.getDomain(), riverCells.getSizeX(), riverCells.getSizeY(), 0);
    for (int i = 0; i < distRiver.getSizeX(); i++) {
        for (int j = 0; j < distRiver.getSizeY(); j++) {
            Index2 p(i, j);
            while (riverCells.isValidCell(p) && riverCells.at(p) == 0) {
                if (!terrainD8->hasDownstream(p)) break;
                p = terrainD8->downstream(p);
            }
            double dx = p.x() - i;
            double dy = p.y() - j;
            distRiver(i, j) = std::sqrt(dx * dx + dy * dy);
        }
    }
    return distRiver;
}

ScalarField2 TerrainAnalysis::DistanceToNearestDrainageFlow() const
{
    ScalarField2 distRiver(riverCells.getDomain(), riverCells.getSizeX(), riverCells.getSizeY(), 0);
    for (int i = 0; i < distRiver.getSizeX(); i++) {
        for (int j = 0; j < distRiver.getSizeY(); j++) {
            Index2 p(i, j);
            double dist = 0;
            while (riverCells.isValidCell(p) && riverCells.at(p) == 0) {
                if (!terrainD8->hasDownstream(p)) break;
                Index2 q = terrainD8->downstream(p);
                double dx = q.x() - p.x();
                double dy = q.y() - p.y();
                dist += std::sqrt(dx * dx + dy * dy);
                p = q;
            }
            distRiver(i, j) = dist;
        }
    }
    return distRiver;
}

ScalarField2 TerrainAnalysis::HeightAboveNearestDrainage() const
{
    ScalarField2 height(riverCells.getDomain(), riverCells.getSizeX(), riverCells.getSizeY(), 0);
    for (int i = 0; i < height.getSizeX(); i++) {
        for (int j = 0; j < height.getSizeY(); j++) {
            Index2 p(i, j);
            while (riverCells.isValidCell(p) && riverCells.at(p) == 0) {
                if (!terrainD8->hasDownstream(p)) break;
                p = terrainD8->downstream(p);
            }
            height(i, j) = hf.at(Index2(i,j)) - hf.at(p);
        }
    }
    return height;
}

ScalarField2 TerrainAnalysis::DepthToWater() const
{
    ScalarField2 slope = hf.Slope();

    ScalarField2 dtw(riverCells.getDomain(), riverCells.getSizeX(), riverCells.getSizeY(), 0);
    for (int i = 0; i < dtw.getSizeX(); i++) {
        for (int j = 0; j < dtw.getSizeY(); j++) {
            Index2 p(i, j);
            double d = 0;
            while (riverCells.isValidCell(p) && riverCells.at(p) == 0) {
                if (!terrainD8->hasDownstream(p)) break;
                Index2 q = terrainD8->downstream(p);
                double dx = q.x() - p.x();
                double dy = q.y() - p.y();
                double dist = std::sqrt(dx * dx + dy * dy);
                d += dist * slope.at(q);
                p = q;
            }
            dtw(i, j) = d;
        }
    }
    return dtw;
}


ScalarField2 TerrainAnalysis::MaximumBranchLength()
{
    if (!terrainD8) {
        terrainD8 = new TerrainFlowD8(hf);
    }

    ScalarField2 Bmax(hf.getDomain(), hf.getSizeX(), hf.getSizeY(), 0);
    for (int i = 0; i < Bmax.getSizeX(); i++) {
        for (int j = 0; j < Bmax.getSizeY(); j++) {
            double maxb = 0;
            Index2 p(i, j);
            for (int d = 0; d < 8; d++) {
                Index2 q = Bmax.neighborCell(p, d);
                if (Bmax.isValidCell(q.x(), q.y())) {
                    double b = terrainD8->flowConfluenceDistance(p, q);
                    if (b > maxb) maxb = b;
                }
            }
            Bmax(i, j) = maxb;
        }
    }
    return Bmax;
}




PPA* TerrainAnalysis::computePPA(int w)
{
    if (!ridgesPPA || w != ridgesPPA_w) {
        ridgesPPA_w = w;
        if (ridgesPPA) delete ridgesPPA;
        ridgesPPA = new PPA(hf);
        ridgesPPA->compute(ridgesPPA_w);
    }
    return ridgesPPA;
}


DivideTree* TerrainAnalysis::computeDivideTree(double minProm)
{
    if (!divideTree) {
        Tile* tile = Tile::fromScalarField(hf);
        Box2 domain = hf.getDomain();
        UtmCoordinateSystem coordSystem(1, domain.getMin()[0], domain.getMin()[1],
                domain.getMax()[0], domain.getMax()[1], hf.getCellSize()[0]);
        TreeBuilder* builder = new TreeBuilder(tile, coordSystem);
        divideTree = builder->buildDivideTree();
        delete builder;
        delete tile;

        if (islandTree) delete islandTree;
        islandTree = new IslandTree(*divideTree);
        islandTree->build(true);
        divideTree->deleteRunoffs(); // does not affect divide tree

        updatedRidgeCells = false;
    }

    if (!prunedDivTree || promCutoff != minProm) {
        promCutoff = minProm;
        delete prunedDivTree;
        delete prunedIslandTree;

        prunedDivTree = new DivideTree(*divideTree);
        prunedDivTree->prune(minProm, *islandTree);
        prunedIslandTree = new IslandTree(*prunedDivTree);
        prunedIslandTree->build(true);
        prunedDivTree->deleteRunoffs();

        peaksProminence.clear();
        peaksIsolation.clear();
        ridgeLevelRak2019.clear();
        ridgeLevelScherler2020.clear();
    }

    return prunedDivTree;
}


std::vector<double> TerrainAnalysis::getPeaksProminence()
{
    if (peaksProminence.size() > 0) return peaksProminence;

    const IslandTree& island = *prunedIslandTree;
    peaksProminence = std::vector<double>(island.nodes().size() - 1);
    for (unsigned int i = 0; i < peaksProminence.size(); i++) {
      peaksProminence[i] = island.nodes()[i + 1].prominence;
    }
    return peaksProminence;
}

std::vector<IsolationRecord> TerrainAnalysis::getPeaksIsolation()
{
    if (peaksIsolation.size() > 0) return peaksIsolation;

    // Compute isolation of the peaks in the divide tree
    const DivideTree& dt = *prunedDivTree;
    IsolationFinder isolFinder(hf);
    peaksIsolation = std::vector<IsolationRecord>(dt.peaks().size());

    // Optimization: use a LineTree to bound isolation search
    LineTree lt(dt);
    lt.build();

    for (unsigned int i = 0; i < dt.peaks().size(); i++) {
        const Offsets& peak = dt.peaks()[i].location;
        if (lt.getLineParent(i+1) > 0) {
            const Offsets& line = dt.peaks()[lt.getLineParent(i+1)-1].location;
            int dx = abs(peak.x() - line.x()) + 1;
            int dy = abs(peak.y() - line.y()) + 1;
            // use distance to line parent (higher) as radius for isolation search
            peaksIsolation[i] = isolFinder.findIsolation(Index2(peak.x(), peak.y()), std::max(dx, dy));
        }
        else {
            peaksIsolation[i] = isolFinder.findIsolation(Index2(peak.x(), peak.y()));
        }
    }

    return peaksIsolation;
}

std::vector<std::vector<int> > divtreePeakNeighbors(const DivideTree& dt)
{
    // returns for every peak the list of peak ids directly connected to it
    std::vector<std::vector<int> > peakNeighs(dt.nodes().size() + 1);
    for (unsigned int i = 0; i < dt.nodes().size(); i++) {
        const DivideTree::Node& node = dt.nodes()[i];
        if (node.saddleId != DivideTree::Node::Null && node.saddleId >= 0
            && dt.getSaddle(node.saddleId).type == Saddle::Type::PROM)
        {
            peakNeighs[i].push_back(node.parentId);
            peakNeighs[node.parentId].push_back(i);
        }
    }
    return peakNeighs;
}

std::vector<std::vector<int> > FloydWarshall(const std::vector<std::vector<int> >& neighbors, std::vector<std::vector<int> >& next)
{
    int n = int(neighbors.size());
    const int DIST_INFTY = 10 * n;
    std::vector<std::vector<int> > dist(n, std::vector<int>(n, DIST_INFTY));
    next = std::vector<std::vector<int> >(n, std::vector<int>(n, -1));

    for (int i = 0; i < n; i++) {
        for (int j : neighbors[i]) {
            dist[i][j] = 1;
            next[i][j] = j;
        }
        dist[i][i] = 0;
        next[i][i] = i;
    }

    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (dist[i][j] > dist[i][k] + dist[k][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                    next[i][j] = next[i][k];
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (dist[i][j] >= DIST_INFTY) {
                dist[i][j] = -1;
            }
        }
    }

    return dist;
}


std::vector<int> getPath(const std::vector<std::vector<int> >& next, int i, int j)
{
    std::vector<int> path;

    // no connection
    if (next[i][j] < 0)
        return path;

    // reconstruct path
    int v = i;
    path.push_back(v);
    while (v != j) {
        v = next[v][j];
        path.push_back(v);
    }
    return path;
}


void TerrainAnalysis::computeRidgesLevelRak2019(std::vector<int>& ridgeLevel, int& maxLevel)
{
  if (ridgeLevelRak2019.size() > 0) {
      ridgeLevel = ridgeLevelRak2019;
      maxLevel = maxRidgeLevelRak2019;
      return;
  }

  const DivideTree& dt = *prunedDivTree;

  std::vector<std::vector<int> > peakNeighs = divtreePeakNeighbors(dt);
  std::vector<std::vector<int> > shortestPathNext;
  std::vector<std::vector<int> > shortestPathDist = FloydWarshall(peakNeighs, shortestPathNext);
  int n = int(peakNeighs.size());

  // find graph radius
  int maxDist = 0;
  int radiusI = -1, radiusJ = -1;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (shortestPathDist[i][j] > maxDist) {
        maxDist = shortestPathDist[i][j];
        radiusI = i;
        radiusJ = j;
      }
    }
  }

  // compute levels
  std::vector<int> peakRidgeLevel = std::vector<int>(n, 0);
  int maxPeakRidgeLevel = 0;

  std::queue<std::pair<int, std::vector<int> > > Q;
  Q.push(std::make_pair(1, getPath(shortestPathNext, radiusI, radiusJ)));
  while (!Q.empty()) {

    int level = Q.front().first;
    maxPeakRidgeLevel = std::max(maxPeakRidgeLevel, level);
    std::vector<int> path = Q.front().second;
    Q.pop();

    // set level to all peaks in path
    for (int i : path) {
      peakRidgeLevel[i] = level;
    }
    // recursively find largest paths of branches
    for (int i : path) {
      for (int j : peakNeighs[i]) {
        // if not visited yet and not part of current path
        if (peakRidgeLevel[j] == 0) {

          // find largest sub-branch
          int maxBranch = 0;
          int branchEnd = j;
          for (int k = 0; k < n; k++) {
            // from j but not through the node we came from (i)
            if (shortestPathDist[j][k] > maxBranch && shortestPathNext[j][k] != i) {
              maxBranch = shortestPathDist[j][k];
              branchEnd = k;
            }
          }
          Q.push(std::make_pair(level + 1, getPath(shortestPathNext, j, branchEnd)));
        }
      }
    }
  }

  // assign levels to ridges as maximum level of both peaks
  maxRidgeLevelRak2019 = 0;
  ridgeLevelRak2019.clear();
  ridgeLevelRak2019.reserve(dt.nodes().size());
  for (unsigned int i = 0; i < dt.nodes().size(); i++) {
      const DivideTree::Node& node = dt.nodes()[i];
      if (node.saddleId != DivideTree::Node::Null && node.saddleId >= 0
         && dt.getSaddle(node.saddleId).type == Saddle::Type::PROM)
      {
          int level = std::max(peakRidgeLevel[i], peakRidgeLevel[node.parentId]);
          ridgeLevelRak2019.push_back(level);
          maxRidgeLevelRak2019 = std::max(maxRidgeLevelRak2019, level);
      }
  }

  ridgeLevel = ridgeLevelRak2019;
  maxLevel = maxRidgeLevelRak2019;
}


void TerrainAnalysis::computeRidgesLevelScherler2020(std::vector<int>& ridgeLevel, int& maxRidgeLevel)
{
  if (ridgeLevelScherler2020.size() > 0) {
      ridgeLevel = ridgeLevelScherler2020;
      maxRidgeLevel = maxRidgeLevelScherler2020;
      return;
  }

  const DivideTree& dt = *prunedDivTree;

  // map all ridges to their two peaks and count neighbors
  std::vector<std::map<int, int> > peakNeighbors = std::vector<std::map<int, int> >(dt.peaks().size() + 1);
  for (unsigned int i = 0; i < dt.nodes().size(); i++) {
    const DivideTree::Node& node = dt.nodes()[i];
    if (node.saddleId != DivideTree::Node::Null && node.saddleId >= 0
      && dt.getSaddle(node.saddleId).type == Saddle::Type::PROM)
    {
      peakNeighbors[i][node.parentId] = i;
      peakNeighbors[node.parentId][i] = i;
    }
  }

  std::vector<int> nodesRidgeLevel(dt.nodes().size(), 0);
  maxRidgeLevel = 1;
  std::vector<std::set<int> > incomingPeaks(dt.peaks().size() + 1);
  std::vector<int> outgoingPeak(dt.peaks().size() + 1, -1);

  // init queue with all single-neighbor peaks
  std::queue<std::pair<int, int> > Q;
  for (unsigned int i = 0; i < peakNeighbors.size(); i++) {
    if (peakNeighbors[i].size() == 1) {
      int peak2 = peakNeighbors[i].begin()->first;
      int ridge = peakNeighbors[i].begin()->second;

      Q.push(std::make_pair(i, peak2));
      outgoingPeak[i] = peak2;
      nodesRidgeLevel[ridge] = 1;
    }
  }

  // traversal from the leafs (single-neighbor nodes)
  while (!Q.empty()) {

    std::pair<int, int> p = Q.front();
    Q.pop();

    int peak1 = p.first;
    int peak2 = p.second;

    // add origin to incoming peaks to peak2
    incomingPeaks[peak2].insert(peak1);

    // if visited all but one neighbor in peak2
    if (incomingPeaks[peak2].size() == peakNeighbors[peak2].size() - 1) {

      // compute the out ridge level
      int inLevel = 0;
      for (int p : incomingPeaks[peak2]) {
        int ridge = peakNeighbors[peak2][p];
        inLevel = std::max(inLevel, nodesRidgeLevel[ridge]);
      }
      int outLevel = inLevel + 1;

      // find the missing neighbor
      for (std::map<int, int>::const_iterator it = peakNeighbors[peak2].begin(); it != peakNeighbors[peak2].end(); it++) {
        if (incomingPeaks[peak2].find(it->first) == incomingPeaks[peak2].end()) {
          int peakOut = it->first;
          Q.push(std::make_pair(peak2, peakOut));
          outgoingPeak[peak2] = peakOut;
          nodesRidgeLevel[it->second] = outLevel;
          maxRidgeLevel = std::max(maxRidgeLevel, nodesRidgeLevel[it->second]);
          break;
        }
      }
    }
  }

  // omit null nodes from ridge list
  maxRidgeLevelScherler2020 = 0;
  ridgeLevelScherler2020.clear();
  ridgeLevelScherler2020.reserve(dt.nodes().size());
  for (unsigned int i = 0; i < dt.nodes().size(); i++) {
      const DivideTree::Node& node = dt.nodes()[i];
      if (node.saddleId != DivideTree::Node::Null && node.saddleId >= 0
         && dt.getSaddle(node.saddleId).type == Saddle::Type::PROM)
      {
          int level = nodesRidgeLevel[i];
          ridgeLevelScherler2020.push_back(level);
          maxRidgeLevelScherler2020 = std::max(maxRidgeLevelScherler2020, level);
      }
  }


  ridgeLevel = ridgeLevelScherler2020;
  maxRidgeLevel = maxRidgeLevelScherler2020;
}



std::vector<Index2> Dijkstra(const ScalarField2& costField, const Index2& ini, const Index2& end)
{
    std::vector<int> parent(costField.getNumElements(), -1);
    int idxIni = costField.cellId(ini.x(), ini.y());
    int idxEnd = costField.cellId(end.x(), end.y());

    // dijkstra
    std::priority_queue<std::tuple<double, int, int> > Q;
    Q.push(std::make_tuple(0, idxIni, idxIni));
    while (!Q.empty()) {
        auto q = Q.top();
        Q.pop();
        double cost = -std::get<0>(q);
        int idxCurr = std::get<1>(q);
        int idxPrev = std::get<2>(q);
        Index2 p = costField.idToCell(idxCurr);

        // already seen?
        if (parent[idxCurr] >= 0) continue;
        parent[idxCurr] = idxPrev;

        // done?
        if (idxCurr == idxEnd) break;

        // go to neighbors
        for (int i = 0; i < 8; i++) {
            Index2 pn = costField.neighborCell(p, i);
            if (costField.isValidCell(pn)) {
                int idxNeigh = costField.cellId(pn.x(), pn.y());
                double cneigh = costField.at(idxNeigh);
                Q.push(std::make_tuple(-cost - cneigh, idxNeigh, idxCurr));
            }
        }
    }

    // reconstruct path
    std::vector<int> pathIdx;
    if (parent[idxEnd] >= 0) {
        pathIdx.push_back(idxEnd);
        int idxCurr = idxEnd;
        while (idxCurr != idxIni) {
            idxCurr = parent[idxCurr];
            pathIdx.push_back(idxCurr);
        }
    }
    std::vector<Index2> path(pathIdx.size());
    for (unsigned int i = 0; i < pathIdx.size(); i++) {
        Index2 p = costField.idToCell(pathIdx[i]);
        path[pathIdx.size() - 1 - i] = p;
    }
    return path;
}


IntField2 TerrainAnalysis::computeRidgeCells()
{
    IntField2 mask(hf.getDomain(), hf.getSizeX(), hf.getSizeY(), 0);

    // we will define a cost favoring the divides
    ScalarField2 ridgeCost = MaximumBranchLength();
    double vmin, vmax;
    ridgeCost.getRange(vmin, vmax);
    for (int i = 0; i < ridgeCost.getNumElements(); i++) {
        ridgeCost[i] = 1 + vmax - ridgeCost.at(i);
    }

    // use the divide tree unpruned.
    // it contains all peaks with any prominence, so there will always
    // exist a monotonically increasing path from a saddle to a peak
    for (unsigned int i = 0; i < divideTree->nodes().size(); i++) {
        const DivideTree::Node& node = divideTree->nodes()[i];
        if (node.saddleId != DivideTree::Node::Null && node.saddleId >= 0 &&
            divideTree->getSaddle(node.saddleId).type == Saddle::Type::PROM)
        {
            const Peak& p1 = divideTree->getPeak(i);
            const Peak& p2 = divideTree->getPeak(node.parentId);
            const Saddle& s = divideTree->getSaddle(node.saddleId);

            std::vector<Index2> pathP1 = Dijkstra(ridgeCost,
                                                  Index2(s.location.x(),  s.location.y()),
                                                  Index2(p1.location.x(), p1.location.y()));
            std::vector<Index2> pathP2 = Dijkstra(ridgeCost,
                                                  Index2(s.location.x(),  s.location.y()),
                                                  Index2(p2.location.x(), p2.location.y()));

            for (const Index2& idx : pathP1) mask(idx) = 1;
            for (const Index2& idx : pathP2) mask(idx) = 1;
        }
    }

    return mask;
}



/*!
\brief Compute the divide ridges.
*/
void TerrainAnalysis::getDivideTreeRidges(bool detailedPath, std::vector<std::vector<Index2> >& ridges,
                                          std::vector<std::pair<int,int> >& ridgePeaks, std::vector<int>& ridgeSaddles)
{
    ridges.clear();

    ScalarField2 ridgeCost;
    if (detailedPath) {
        if (!updatedRidgeCells) {
            ridgeCells = computeRidgeCells();
            updatedRidgeCells = true;
        }
        ridgeCost = ScalarField2(ridgeCells.getDomain(), ridgeCells.getSizeX(), ridgeCells.getSizeY());
        for (int i = 0; i < ridgeCost.getNumElements(); i++) {
            ridgeCost[i] = ridgeCells[i] > 0 ? 1 : 1e6;
        }
    }


    for (unsigned int i = 0; i < prunedDivTree->nodes().size(); i++) {
        const DivideTree::Node& node = prunedDivTree->nodes()[i];
        if (node.saddleId != DivideTree::Node::Null && node.saddleId >= 0
           && prunedDivTree->getSaddle(node.saddleId).type == Saddle::Type::PROM)
        {
            const Peak& p1 = prunedDivTree->getPeak(i);
            const Peak& p2 = prunedDivTree->getPeak(node.parentId);
            const Saddle& s = prunedDivTree->getSaddle(node.saddleId);

            ridgePeaks.push_back(std::make_pair(i, node.parentId));
            ridgeSaddles.push_back(node.saddleId);

            if (detailedPath) {
                std::vector<Index2> pathP1 = Dijkstra(ridgeCost,
                                                      Index2(s.location.x(),  s.location.y()),
                                                      Index2(p1.location.x(), p1.location.y()));
                std::vector<Index2> pathP2 = Dijkstra(ridgeCost,
                                                      Index2(s.location.x(),  s.location.y()),
                                                      Index2(p2.location.x(), p2.location.y()));

                std::vector<Index2> ridgeLine;
                ridgeLine.reserve(pathP1.size() + pathP2.size());
                for (unsigned int j = 0; j < pathP1.size(); j++)
                    ridgeLine.push_back(pathP1[pathP1.size() - 1 - j]);
                for (unsigned int j = 0; j < pathP2.size(); j++)
                    ridgeLine.push_back(pathP2[j]);

                ridges.push_back(ridgeLine);

            } else {
                ridges.push_back({ Index2(p1.location.x(), p1.location.y()),
                                   Index2(s.location.x(), s.location.y()),
                                   Index2(p2.location.x(), p2.location.y())
                                 });
            }
        }
    }
}


ScalarField2 TerrainAnalysis::NearestNeighborIndex(double r, double& globalNNI) const
{
    ScalarField2 nni(hf.getDomain(), hf.getSizeX(), hf.getSizeY(), 0);
    if (!prunedDivTree) return nni;

    // peak coords
    int npeaks = int(prunedDivTree->peaks().size());
    std::vector<Vector2> peakCoords(npeaks);
    for (int i = 0; i < npeaks; i++) {
        const Peak& p = prunedDivTree->peaks()[i];
        peakCoords[i] = hf.domainCoords(p.location.x(), p.location.y());
    }

    // for each peak, find its nearest neighbor
    std::vector<double> closestDist(npeaks);
    for (int i = 0; i < npeaks; i++) {
        double d2min = 1e32;
        for (int j = 0; j < npeaks; j++) {
            if (i == j) continue;
            double d2 = SquaredNorm(peakCoords[i] - peakCoords[j]);
            if (d2 < d2min) {
                closestDist[i] = std::sqrt(d2);
                d2min = d2;
            }
        }
    }

    // compute NNI
    double r2 = r*r;
    for (int i = 0; i < hf.getSizeX(); i++) {
        for (int j = 0; j < hf.getSizeY(); j++) {
            Vector2 pij = hf.domainCoords(i, j);
            int n = 0;
            double dsum = 0;
            // get peaks in domain
            for (int k = 0; k < npeaks; k++) {
                double d2 = SquaredNorm(pij - peakCoords[k]);
                if (d2 < r2) {
                    n++;
                    dsum += closestDist[k];
                }
            }
            if (n > 0) {
                //TODO: circle-box overlap
                //QRect overlapBBox = hf.VertexIntegerArea(i, j, r);
                //double areaOverlapBox = overlapBBox.width() * overlapBBox.height() * hf.getCellArea();
                //double approxArea = 0.25 * Math::Pi * areaOverlapBox;
                double minDistAvg = dsum/n;
                double expectDist = 0.5*std::sqrt(Math::Pi*r2/n);
                nni(i, j) = minDistAvg / expectDist;
            }
        }
    }

    // overall NNI
    double dsum = 0;
    for (double d : closestDist) {
        dsum += d;
    }
    globalNNI = (dsum/npeaks)/(0.5*std::sqrt(hf.getNumElements()*hf.getCellArea()/npeaks));

    return nni;
}
