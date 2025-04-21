#ifndef _ISOLATION_FINDER_H_
#define _ISOLATION_FINDER_H_

#include "heightfield.h"
#include <vector>

struct IsolationRecord {
	bool foundHigherGround;
    Index2 closestHigherGround;
	double distance;  // distance to peak in meters

	IsolationRecord()
		: foundHigherGround(false),
		closestHigherGround(0, 0),
		distance(0) {
	}

	IsolationRecord(const IsolationRecord& other)
		: foundHigherGround(other.foundHigherGround),
		closestHigherGround(other.closestHigherGround),
		distance(other.distance) {
	}

	void operator=(const IsolationRecord& other) {
		foundHigherGround = other.foundHigherGround;
		closestHigherGround = other.closestHigherGround;
		distance = other.distance;
	}
};


class IsolationFinder {
private:
	const HeightField& heights; //!< Reference to heightfield
public:
    //! Constructor.
    explicit IsolationFinder(const HeightField& hf) : heights(hf) {}

    IsolationRecord findIsolation(const Index2& peak);
    IsolationRecord findIsolation(const Index2& peak, int w);
};


#endif
