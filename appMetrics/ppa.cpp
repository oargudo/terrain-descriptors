#include "ppa.h"
#include <algorithm>

PPA::PPA(const HeightField& h) : hf(h)
{
  nx = hf.getSizeX();
  ny = hf.getSizeY();
  hf.getRange(hmin, hmax);
}

void PPA::compute(int profileLength)
{
  ridgeCandidate = std::vector<bool>(nx * ny, false);
  isSegment = std::vector<std::vector<bool>>(nx * ny, std::vector<bool>(NUM_DIRS, false));
  segHeight = std::vector<std::vector<double>>(nx * ny, std::vector<double>(NUM_DIRS, 0.0));
  isReliable = std::vector<std::vector<bool>>(nx * ny, std::vector<bool>(NUM_DIRS, false));

  computeCandidatePoints(profileLength);
  connectCandidateSegments();
  checkReliableSegments(false);

  computeRidgesMST();

  const int pruneLeavesIterations = profileLength; // half profile length iterations (recommended by authors)
  const int minLeafBranchSize = profileLength;
  ridgesPrunedMST = RidgesTreeMatrix(ridgesMST);
  for (int i = 0; i < pruneLeavesIterations; i++) {
    //ridgesPrunedMST = pruneRidgeLeaves(ridgesPrunedMST, reliableEdgesMST);
    ridgesPrunedMST = pruneRidgeLeaves(ridgesPrunedMST, EdgeSet()); // do not account for reliables
  }
  ridgesPrunedMST = pruneSmallBranches(ridgesPrunedMST, minLeafBranchSize);
}


std::vector<Index2> PPA::getRidgeCandidates() const {
  std::vector<Index2> pts;
  pts.reserve(nx * ny);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      if (ridgeCandidate[hf.cellId(i, j)]) {
        pts.push_back(Index2(i, j));
      }
    }
  }
  return pts;
}

std::vector<PPA::RidgeSegment> PPA::getSegments() const {
  std::vector<PPA::RidgeSegment> segs;
  segs.reserve(4 * nx * ny);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int d = 0; d < NUM_DIRS; d++) {
        if (isSegment[hf.cellId(i, j)][d]) {
          segs.push_back(PPA::RidgeSegment(
            Index2(i, j),
            Index2(i + PROFILE_DIRS[d][0], j + PROFILE_DIRS[d][1])
          ));
        }
      }
    }
  }
  return segs;
}

std::vector<PPA::RidgeSegment> PPA::getReliableSegments() const {
  std::vector<PPA::RidgeSegment> segs;
  segs.reserve(4 * nx * ny);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int d = 0; d < NUM_DIRS; d++) {
        if (isSegment[hf.cellId(i, j)][d] && isReliable[hf.cellId(i, j)][d]) {
          segs.push_back(PPA::RidgeSegment(
            Index2(i, j),
            Index2(i + PROFILE_DIRS[d][0], j + PROFILE_DIRS[d][1])
          ));
        }
      }
    }
  }
  return segs;
}

std::vector<PPA::RidgeSegment> PPA::getSegmentsInMST() const
{

  std::vector<PPA::RidgeSegment> segs;
  for (int n1 = 0; n1 < nx * ny; n1++)
  {
    for (int n2 : ridgesMST[n1])
    {
      if (n2 < n1) continue; // avoid duplicates

      int i1, i2, j1, j2;
      nodeToIndices(n1, i1, j1);
      nodeToIndices(n2, i2, j2);
      segs.push_back(RidgeSegment(Index2(i1, j1), Index2(i2, j2)));
    }
  }
  return segs;
}

std::vector<PPA::RidgeSegment> PPA::getSegmentsInPrunedMST() const
{

  std::vector<PPA::RidgeSegment> segs;
  for (int n1 = 0; n1 < nx * ny; n1++)
  {
    for (int n2 : ridgesPrunedMST[n1])
    {
      if (n2 < n1) continue; // avoid duplicates

      int i1, i2, j1, j2;
      nodeToIndices(n1, i1, j1);
      nodeToIndices(n2, i2, j2);
      segs.push_back(RidgeSegment(Index2(i1, j1), Index2(i2, j2)));
    }
  }
  return segs;
}

inline void PPA::getDirection(int dir, int& dx, int& dy) const
{
  int i = dir % 8;
  if (i < 0) i += 8;

  if (i < 4)
  {
    dx = PROFILE_DIRS[i][0];
    dy = PROFILE_DIRS[i][1];
  }
  else
  {
    dx = -PROFILE_DIRS[i - 4][0];
    dy = -PROFILE_DIRS[i - 4][1];
  }
}

void PPA::computeCandidatePoints(int profileLength)
{
  int profileHalf = profileLength / 2;

  for (int i = 0; i < hf.getSizeX(); i++)
  {
    for (int j = 0; j < hf.getSizeY(); j++)
    {

      double ctrHeight = hf.at(i, j);

      // check if we have a peak in any of the profile directions
      for (int d = 0; d < NUM_DIRS; d++) {

        // a peak implies finding a lower point along either profile half
        bool lowerLeft = false;
        bool lowerRight = false;
        for (int k = 1; k <= profileHalf; k++) {
          int ileft = i - PROFILE_DIRS[d][0] * k;
          int jleft = j - PROFILE_DIRS[d][1] * k;
          if (hf.isValidCell(ileft, jleft) && hf.at(ileft, jleft) < ctrHeight)
            lowerLeft = true;

          int iright = i + PROFILE_DIRS[d][0] * k;
          int jright = j + PROFILE_DIRS[d][1] * k;
          if (hf.isValidCell(iright, jright) && hf.at(iright, jright) < ctrHeight)
            lowerRight = true;
        }

        // if lower point on both sides, mark as ridge candidate (we can skip other dirs)
        if (lowerLeft && lowerRight) {
          ridgeCandidate[hf.cellId(i, j)] = true;
          break;
        }
      }
    }
  }
}

void PPA::connectCandidateSegments()
{
  for (int i = 0; i < hf.getSizeX(); i++) {
    for (int j = 0; j < hf.getSizeY(); j++) {
      int idx = hf.cellId(i, j);
      // if ridge candidate, look for neighbors
      if (ridgeCandidate[idx]) {
        for (int d = 0; d < NUM_DIRS; d++) {
          int in = i + PROFILE_DIRS[d][0];
          int jn = j + PROFILE_DIRS[d][1];
          if (!hf.isValidCell(in, jn)) continue;

          int idn = hf.cellId(in, jn);
          if (ridgeCandidate[idn]) {
            isSegment[idx][d] = true;
            segHeight[idx][d] = 0.5 * (hf.at(idx) + hf.at(in, jn));

            // if both diagonals in a square are marked, remove lower one
            int idd = hf.cellId(i - 1, j);
            if (d == 3 && i > 0 && isSegment[idd][1]) {
              if (segHeight[idx][3] > segHeight[idd][1])
                isSegment[idd][1] = false;
              else
                isSegment[idx][3] = false;
            }
          }
        }
      }
    }
  }
}

void PPA::checkReliableSegments(bool removeParallels)
{
  for (int i = 1; i < hf.getSizeX() - 1; i++)
  {
    for (int j = 1; j < hf.getSizeY() - 1; j++)
    {
      int idx = hf.cellId(i, j);
      //double h = hf.at(idx);
      for (int d = 0; d < NUM_DIRS; d++)
      {

        if (!isSegment[idx][d]) continue;

        int dx, dy;
        getDirection(d, dx, dy);
        double h = segHeight[idx][d];

        // NS or EW axis
        if (d % 2 == 0) 
        {
          // left-hand parallel segment height
          int dx1, dy1;
          getDirection(d + 2, dx1, dy1);
          int i1 = i + dx1, j1 = j + dy1;
          double h1 = 0.5 * (hf.at(i1, j1) + hf.at(i1 + dx, j1 + dy));

          // right-hand parallel segment height
          int dx2, dy2;
          getDirection(d - 2, dx2, dy2);
          int i2 = i + dx2, j2 = j + dy2;
          double h2 = 0.5 * (hf.at(i2, j2) + hf.at(i2 + dx, j2 + dy));

          // reliable if higher than parallel neighbors
          if (h > h1 && h > h2) {
            isReliable[idx][d] = true;
            if (removeParallels) {
              isSegment[hf.cellId(i1, j1)][d] = false;
              isSegment[hf.cellId(i2, j2)][d] = false;
            }
          }
        }
        // diagonal axis
        else
        {
          // left-hand corner height
          int dx1, dy1;
          getDirection(d + 1, dx1, dy1);
          int i1 = i + dx1, j1 = j + dy1;
          double h1 = hf.at(i1, j1);

          // right-hand corner height
          int dx2, dy2;
          getDirection(d - 1, dx2, dy2);
          int i2 = i + dx2, j2 = j + dy2;
          double h2 = hf.at(i2, j2);

          // reliable if higher than opposite corners at either side
          if (h > h1 && h > h2) {
            isReliable[idx][d] = true;
            if (removeParallels) {
              isSegment[hf.cellId(i1, j1)][d] = false;
              isSegment[hf.cellId(i2, j2)][d] = false;

              // opposite direction segment as well
              // note that [(i,j)][d+4] is encoded in [(i,j)+(dxo,dyo)][d]
              int dxo, dyo;
              getDirection(d + 4, dxo, dyo);
              isSegment[hf.cellId(i1 + dxo, j1 + dyo)][d] = false;
              isSegment[hf.cellId(i2 + dxo, j2 + dyo)][d] = false;
            }
          }
        }
      }
    }
  }
}

std::vector<PPA::RidgeEdge> PPA::getRidgeEdges(bool favourReliables)
{
  // create edges with cost
  std::vector<std::tuple<double, int, int> > edges;

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      int id1 = hf.cellId(i, j);
      for (int d = 0; d < NUM_DIRS; d++) {
        if (isSegment[id1][d]) {
          int id2 = hf.cellId(i + PROFILE_DIRS[d][0], j + PROFILE_DIRS[d][1]);
          double w = segHeight[id1][d];

          // reliable: less cost for this segment
          if (favourReliables && isReliable[id1][d]) {
            w = hmax - w;
          }
          else {
            w = 2 * hmax - w;
          }

          edges.push_back(std::make_tuple(w, id1, id2));
        }
      }
    }
  }

  return edges;
}

// To represent Disjoint Sets 
class DisjointSets
{
public:
  std::vector<int> parent;
  std::vector<int> rank;

  DisjointSets(int n) {
    parent = std::vector<int>(n + 1);
    rank = std::vector<int>(n + 1);
    // Initially, all vertices are in different sets and have rank 0 
    for (int i = 0; i <= n; i++) {
      rank[i] = 0;
      parent[i] = i;
    }
  }

  // Find the parent of a node 'u' and do path compression
  int findParent(int u) {
    // Make the parent of the nodes in the path from u--> parent[u] point to parent[u]
    if (u != parent[u])
      parent[u] = findParent(parent[u]);
    return parent[u];
  }

  // Union by rank 
  void merge(int x, int y) {
    x = findParent(x);
    y = findParent(y);

    // Make tree with smaller height a subtree of the other tree
    if (rank[x] > rank[y])
      parent[y] = x;
    else
      parent[x] = y;

    if (rank[x] == rank[y])
      rank[y]++;
  }
};

void PPA::computeRidgesMST()
{
  // get edges
  std::vector<RidgeEdge> edges = getRidgeEdges(true);

  // MST using Kruskal algorithm
  std::sort(edges.begin(), edges.end());
  DisjointSets ds(nx * ny);
  std::vector<RidgeEdge> mstEdges;

  for (auto e : edges) {

    int p1 = std::get<1>(e);
    int p2 = std::get<2>(e);
    int s1 = ds.findParent(p1);
    int s2 = ds.findParent(p2);

    // Check if the selected edge is creating a cycle or not 
    // A cycle is created if u and v belong to same set
    if (s1 != s2) {
      // if no cycles, then keep this edge
      mstEdges.push_back(e);

      // merge sets
      ds.merge(s1, s2);
    }
  }

  // build sparse matrix
  ridgesMST.clear();
  ridgesMST.resize(nx * ny);
  for (auto e : mstEdges) {
    int n1 = std::get<1>(e);
    int n2 = std::get<2>(e);
    addEdge(ridgesMST, n1, n2);

    double w = std::get<0>(e);
    if (w < hmax) {
      reliableEdgesMST.insert(std::make_pair(n1, n2));
      reliableEdgesMST.insert(std::make_pair(n2, n1));
    }
  }
}

void PPA::addEdge(PPA::RidgesTreeMatrix& matrix, int n1, int n2) {
  matrix[n1].insert(n2);
  matrix[n2].insert(n1);
}

void PPA::removeEdge(PPA::RidgesTreeMatrix& matrix, int n1, int n2) {
  matrix[n1].erase(n2);
  matrix[n2].erase(n1);
}

PPA::RidgesTreeMatrix PPA::pruneRidgeLeaves(const RidgesTreeMatrix& ridges, const EdgeSet& reliable)
{
  // note that we make edits on the copy but always query the original
  // otherwise, unwanted effect happen. For example, imagine we have a
  // tree: A - B - C and want to do do one iteration of leaves removal 
  // We should obtain: B
  // If we modify the original and B is processed after removing A, 
  // it will be considered a leaf and we will obtain an empty tree

  RidgesTreeMatrix pruned(ridges);
  for (int i = 0; i < int(ridges.size()); i++) {
    // is this node a leaf?
    if (ridges[i].size() == 1) {
      int neigh = *ridges[i].begin();
      // unless marked as reliable, remove this segment
      if (reliable.find(std::make_pair(i, neigh)) == reliable.end()) {
        removeEdge(pruned, i, neigh);
      }
    }
  }
  return pruned;
}

PPA::RidgesTreeMatrix PPA::pruneSmallBranches(const RidgesTreeMatrix& ridges, int minBranchLength)
{
  RidgesTreeMatrix pruned(ridges);

  for (int i = 0; i < int(ridges.size()); i++) {
    // is this node a leaf?
    if (ridges[i].size() == 1) {

      // branch until first bifurcation (a node with 3 or more neighbors)
      std::unordered_set<int> branchNodes = { i };
      int curr = *ridges[i].begin(); // we know we have exactly 1 neigh
      while (ridges[curr].size() <= 2 && int(branchNodes.size()) < minBranchLength) {

        branchNodes.insert(curr);

        if (ridges[curr].size() < 2) // e.g. a line graph: A - B - Curr
          break;

        for (int next : ridges[curr]) {
          if (branchNodes.find(next) == branchNodes.end()) {
            curr = next;
          }
        }
      }

      // delete branch if small
      if (int(branchNodes.size()) < minBranchLength) {
        for (int n : branchNodes) {
          for (int neigh : pruned[n]) {
            removeEdge(pruned, n, neigh);
          }
          pruned[n].clear();
        }
      }
    }
  }

  return pruned;
}
