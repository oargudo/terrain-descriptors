#ifndef _PPA_
#define _PPA_

#include "heightfield.h"
#include <unordered_set>
#include <set>

class PPA {

public:
	PPA(const HeightField& h);
	virtual ~PPA() {}

	void compute(int profileLength = 5);

    typedef std::pair<Index2, Index2> RidgeSegment;
	typedef std::set<std::pair<int, int> > EdgeSet;
	typedef std::vector<std::unordered_set<int> > RidgesTreeMatrix;

    std::vector<Index2>  getRidgeCandidates() const;
	std::vector<RidgeSegment> getSegments() const;
	std::vector<RidgeSegment> getReliableSegments() const;
	std::vector<RidgeSegment> getSegmentsInMST() const;
	std::vector<RidgeSegment> getSegmentsInPrunedMST() const;

	RidgesTreeMatrix getMST() const { return ridgesMST; };
	RidgesTreeMatrix getPrunedMST() const { return ridgesPrunedMST; };

protected:

	static const int NUM_DIRS = 4;
	const int PROFILE_DIRS[NUM_DIRS][2] = {
		{1, 0}, // E
		{1, 1}, // NE
		{0, 1}, // N
		{-1, 1} // NW
	};
	typedef std::tuple<double, int, int> RidgeEdge;

	void getDirection(int dir, int& dx, int& dy) const;
	void nodeToIndices(int nodeId, int& i, int& j) const {
		i = nodeId % nx;
		j = nodeId / nx;
	}

	void computeCandidatePoints(int profileLength = 5);
	void connectCandidateSegments();
	void checkReliableSegments(bool removeParallels);

	std::vector<RidgeEdge> getRidgeEdges(bool favourReliables);
	void computeRidgesMST();

	static void addEdge(RidgesTreeMatrix& tree, int node1, int node2);
	static void removeEdge(RidgesTreeMatrix& tree, int node1, int node2);
	static RidgesTreeMatrix pruneRidgeLeaves(const RidgesTreeMatrix& ridges, const EdgeSet& reliable);
	static RidgesTreeMatrix pruneSmallBranches(const RidgesTreeMatrix& ridges, int minBranchLength);

	
	const HeightField& hf;
	int nx, ny;
	double hmin, hmax;

	std::vector<bool> ridgeCandidate;
	std::vector<std::vector<bool>>   isSegment;
	std::vector<std::vector<double>> segHeight;
	std::vector<std::vector<bool>>   isReliable;

	EdgeSet reliableEdgesMST;
	RidgesTreeMatrix ridgesMST;
	RidgesTreeMatrix ridgesPrunedMST;
};


#endif
