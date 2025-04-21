#ifndef TERRAINANALYSIS_H
#define TERRAINANALYSIS_H

#include "heightfield.h"
#include "riversnet.h"
#include "ppa.h"
#include "divide_tree.h"
#include "island_tree.h"
#include "isolation.h"


class TerrainAnalysis
{
public:
    explicit TerrainAnalysis(const HeightField& h);
    ~TerrainAnalysis();

    TerrainFlowD8 *getFlowDirs();
    std::vector<RiverTree*> computeRivers(double cit_s, double cit_t);
    IntField2 RiversMask() const { return riverCells; }
    IntField2 RiversIdMask() const;
    ScalarField2 DistanceToNearestDrainageEuclidean() const;
    ScalarField2 DistanceToNearestDrainageFlow() const;
    ScalarField2 HeightAboveNearestDrainage() const;
    ScalarField2 DepthToWater() const;

    ScalarField2 MaximumBranchLength();

    PPA* computePPA(int w);

    DivideTree* computeDivideTree(double minProm);
    std::vector<double> getPeaksProminence();
    std::vector<IsolationRecord> getPeaksIsolation();
    void computeRidgesLevelRak2019(std::vector<int>& ridgeLevel, int& maxRidgeLevel);
    void computeRidgesLevelScherler2020(std::vector<int>& ridgeLevel, int& maxRidgeLevel);
    void getDivideTreeRidges(bool detailed, std::vector<std::vector<Index2> >& ridges,
                             std::vector<std::pair<int,int> >& ridgePeaks, std::vector<int>& ridgeSaddles);
    IntField2 computeRidgeCells();

    ScalarField2 NearestNeighborIndex(double rad, double& globalNNI) const;

protected:

    const HeightField& hf;

    TerrainFlowD8* terrainD8 = nullptr;
    IntField2 riverCells;
    std::vector<RiverTree*> riverTrees;
    double riversCIT_s = 0, riversCIT_t = 0;
    bool riversComputed = false;

    PPA* ridgesPPA = nullptr;
    int ridgesPPA_w = 0;

    DivideTree* divideTree = nullptr;
    DivideTree* prunedDivTree = nullptr;
    IslandTree* islandTree = nullptr;
    IslandTree* prunedIslandTree = nullptr;
    int promCutoff = 0;
    IntField2 ridgeCells;
    bool updatedRidgeCells = false;
    std::vector<double> peaksProminence;
    std::vector<IsolationRecord> peaksIsolation;
    std::vector<int> ridgeLevelRak2019, ridgeLevelScherler2020;
    int maxRidgeLevelRak2019 = 0, maxRidgeLevelScherler2020 = 0;
};

#endif // TERRAINANALYSIS_H
