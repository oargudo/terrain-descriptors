#include "isolation.h"

/*!
\brief Find isolation.
\param peak Encoded offsets.
*/ 
IsolationRecord IsolationFinder::findIsolation(const Index2& peak)
{
    int dx = std::max(heights.getSizeX() - peak.x(), peak.x());
    int dy = std::max(heights.getSizeY() - peak.y(), peak.y());
    return findIsolation(peak, std::max(dx, dy));
}

/*!
\brief Find isolation.
\param peak Encoded offsets.
\param w Search window radius.
*/
IsolationRecord IsolationFinder::findIsolation(const Index2& peak, int w)
{
    IsolationRecord record;
    double peakElev = heights.at(peak.x(), peak.y());

    double cx = heights.getCellSize()[0];
    double cy = heights.getCellSize()[1];
    IndexArea area = heights.indexAreaFromHalfWindow(peak, w);

    for (int x = area.xmin(); x <= area.xmax(); x++) {
        for (int y = area.ymin(); y <= area.ymax(); y++) {
            if (heights.at(x, y) > peakElev) {
                double dx = (x - peak.x())*cx;
                double dy = (y - peak.y())*cy;
                double dist = sqrt(double(dx * dx + dy * dy));
                if (!record.foundHigherGround || record.distance > dist) {
                    record.foundHigherGround = true;
                    record.distance = dist;
                    record.closestHigherGround = Index2(x, y);
                }
            }
        }
    }
    return record;
}
