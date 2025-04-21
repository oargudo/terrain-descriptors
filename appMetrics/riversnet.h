#ifndef RIVERSNET_H
#define RIVERSNET_H

#include "integerfield.h"
#include "heightfield.h"
#include <vector>


class RiverTree
{
public:
    RiverTree() {};
    ~RiverTree() {
        for (RiverTree* r : parents) delete r;
    }

    void computeRiverMetrics(const HeightField& hf, const ScalarField2& streamArea);

    void markRiverCells(IntField2& mask, int value = 1) const;
    void markRiverLength(ScalarField2& rivLength) const;
    void markRiverStreamArea(ScalarField2& streamArea) const;

    int getNumParentRivers() const { return int(parents.size()); }
    RiverTree* getParentRiver(int i) const { return parents[i]; }
    const std::vector<RiverTree*>& getParentRivers() const { return parents; }

    int getNumCells() const { return int(cells.size()); }
    const Index2& getCell(int i) const { return cells[i]; }

    const std::vector<Index2>& getStreamCells() const { return cells; }
    const std::vector<double>& getStreamArea() const { return nodeStreamArea; }
    const std::vector<double>& getStreamLength() const { return nodeAccumLength; }
    const std::vector<double>& getStreamSlope() const { return nodeSlope; }
    int getStrahlerIndex()  const { return strahlerIndex; }
    double getSinuosity() const { return sinuosity; }

protected:
    friend class TerrainFlowD8;

    std::vector<RiverTree*> parents;

    std::vector<Index2> cells;
    std::vector<double> nodeStreamArea;
    std::vector<double> nodeAccumLength;
    std::vector<double> nodeSlope;

    int strahlerIndex = 1;
    double sinuosity = 1;
};


class DrainageBasin
{
public:
    DrainageBasin(const IntField2& m) : basinMask(m) {};
    ~DrainageBasin() {};

    bool inBasin(int i, int j) const { return basinMask.at(i, j) > 0; }
    bool inBasin(const Index2& p) const { return inBasin(p.x(), p.y()); }
    bool inPerimeter(int i, int j) const;

    void computeBasinMetrics(const HeightField& hf, const ScalarField2& riverLength, const ScalarField2& riverDrainage);

    double getArea() const { return countInner * basinMask.getCellArea(); };
    double getPerimeter() const { return countPerim * basinMask.getCellSize()[0]; };
    double getCircularity() const { return 4*Math::Pi*getArea() / (getPerimeter() * getPerimeter()); };
    double getCompactness() const { return getPerimeter() / getArea(); };
    double getMeltonRuggednessIndex() const { return (hmax - hmin)/std::sqrt(getArea()); };
    double getHypsommetricIntegral() const { return (havg - hmin)/(hmax - hmin); };
    double getDrainageDensity() const { return countRiver * basinMask.getCellSize()[0] / getArea(); };
    double getReliefRatio() const { return (hmax - hmin)/maxStreamLength; }
    double getFormFactor() const { return maxStreamDrainage / (maxStreamLength * maxStreamLength); };

protected:
    IntField2 basinMask;

    int countInner = 0;
    int countPerim = 0;
    int countRiver = 0;

    double maxStreamLength = 1;
    double maxStreamDrainage = 0;
    double hmax = 0;
    double hmin = 0;
    double havg = 0;
};


class TerrainFlowD8
{
public:
    explicit TerrainFlowD8(const HeightField&);
    ~TerrainFlowD8();

    bool hasDownstream(int i, int j) const { return flowdir.at(i, j) >= 0; }
    bool hasDownstream(const Index2& p) const { return flowdir.at(p) >= 0; }
    Index2 downstream(int i, int j) const;
    Index2 downstream(const Index2& p) const;
    std::vector<Index2> upstream(int i, int j) const;
    std::vector<Index2> upstream(const Index2& p) const;

    double flowConfluenceDistance(const Index2& p, const Index2& q) const;

    ScalarField2 streamArea() const { return upArea; }
    IntField2 riverCellsFromChannelThreshold(double s = 0, double t = 200, bool propagate = true) const;
    RiverTree* buildRiverTree(const Index2& endpoint, const IntField2& mask) const;    
    DrainageBasin getDrainageBasin(const Index2& point) const;

protected:
    const HeightField& hf;
    IntField2 flowdir;
    ScalarField2 upArea;

    static constexpr int neighborOutdirIntoCell[8] = { 4, 5, 6, 7, 0, 1, 2, 3};
};








#endif // RIVERSNET_H
