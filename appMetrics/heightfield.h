#ifndef HEIGHTFIELD_H
#define HEIGHTFIELD_H

#include "scalarfield.h"
#include "integerfield.h"

class HeightField : public ScalarField2
{
public:
    HeightField() {}
    HeightField(const ScalarField2&);
    explicit HeightField(const Box2&, const QImage&, const double& = 0.0, const double& = 256.0 * 256.0 - 1.0, bool = true);
    explicit HeightField(const Box2&, int, int, const double& = 0.0);
    explicit HeightField(const Box2&, int, int, const std::vector<double>&);
    ~HeightField() {}

    // Vertices, elevations and bbox
    double  Height(int i, int j) const;
    double  Height(const Vector2& p, bool triangular = true) const;
    Vector3 Vertex(const Vector2& p, bool triangular = true) const;
    Vector3 Vertex(int i, int j) const;
    Box3 getBox() const;    

    // Gradients and normals
    Vector2 Gradient(int i, int j) const;
    Vector2 Gradient(const Vector2& p) const;
    Vector3 Normal(const Vector2&,  bool = true) const;
    Vector3 Normal(int, int) const;
    double  SlopeDirectional(const Index2&, const Index2&) const;

    // Quadric surface approximation
    QuadricSurface FitQuadric(int i, int j, int w, bool interpolate) const;


    /***********************************/
    /*          LOCAL METRICS          */
    /***********************************/

    // Slope
    ScalarField2 GradientNorm() const;
    ScalarField2 SlopeAverage() const;
    ScalarField2 Aspect() const;
    // For code compatibility and readability, define our default Slope method
    ScalarField2 Slope() const { return GradientNorm(); }

    // Laplacian
    ScalarField2 Laplacian() const;
    ScalarField2 FractLaplacian(double, int) const;

    // Curvatures. If w = 0, use approx from 3x3 window. Otherwise, fit quadric
    enum CurvatureType {
        MIN, MAX, MEAN, GAUSSIAN, PROFILE, CONTOUR, TANGENTIAL, N_CURVATURE_TYPES
    };
    ScalarField2 Curvature(CurvatureType type, int w=0) const;

    // Roughness
    ScalarField2 LocalRelief(int w) const;
    ScalarField2 LocalVariance(int w) const;
    ScalarField2 AreaRatio(int w) const;
    ScalarField2 RuggednessIndex(int w) const;
    ScalarField2 TopographicPositionIndex(int w) const;
    ScalarField2 TopographicPositionIndexSAT(int w) const;
    ScalarField2 SurfaceRoughness(int w, bool dev) const;
    ScalarField2 SurfaceRoughnessSAT(int w, bool dev) const;

    // Asymmetry
    ScalarField2 HillslopeAsymmetry(int w, double direction, double tolerance) const;


    /**********************************/
    /*           VISIBILITY           */
    /**********************************/

    // Ray intersections
    double LipschitzK() const;
    bool Intersect(const Ray&, double&, Vector3&, const Box3&, const double&, const double& = 1.0e8, const double& = 1.0e-4) const;

    // Viewsheds
    ScalarField2 Viewshed(const Vector3& viewpoint, int& viewCnt) const;
    ScalarField2 ViewshedTotal(bool outgoing, double offset) const;
    ScalarField2 ViewshedTotalSampled(int numSamples, bool outgoing, double offset) const;

    // Openness
    ScalarField2 Openness(bool positive, double maxDist, int numDirs) const;

    // Lighting and sun
    ScalarField2 DirectLight(const Vector3&, bool = false) const;
    ScalarField2 SelfShadow(const Vector3& light) const;
    ScalarField2 Sun(const double& latitude, int daystep, int hourstep, int dayoffset, int houroffset) const;

    ScalarField2 Accessibility(const double& r, int n, bool skyView = false) const;
    ScalarField2 SkyViewFactor(const double& r, int n) const { return Accessibility(r, n, true); };
    ScalarField2 SkyViewFactorApproximation() const;
    ScalarField2 DiurnalAnisotropicHeatIndex(double maxAspect) const;


    /**********************************/
    /*            LANDFORMS           */
    /**********************************/

    IntField2 LandformsDikauWood(int w, double tSlope, double tCurvs) const;
    IntField2 LandformsFuzzyDW(ScalarField2& entropy, int minScale, int maxScale, double tSlope, double tCurv) const;
    IntField2 LandformsTPI(double radiusSmall, double radiusLarge, double flatSlope = 0.2) const;
    IntField2 Geomorphons(double maxDist, double flatTangent = 0.02) const;
    ScalarField2 BlackWhiteTopHatTransform(int w, double tPeak, double tValley) const;


    /**********************************/
    /*            HYDROLOGY           */
    /**********************************/

    // Modifiers
    int fillDepressions(const double& eps = 0);
    int fillDepressions(const double&, ScalarField2&) const;
    void completeBreach();

    // Flow
    IntField2 FlowDirectionD8() const;

    ScalarField2 StreamAreaD8  (const std::vector<Index2>& sources = std::vector<Index2>()) const;
    ScalarField2 StreamAreaMFD (const std::vector<Index2>& sources = std::vector<Index2>(), double gamma = 1.0) const;
    ScalarField2 StreamAreaDinf(const std::vector<Index2>& sources = std::vector<Index2>()) const;
    ScalarField2 StreamAreaRho8(const std::vector<Index2>& sources = std::vector<Index2>(), int iters = 20) const;
    ScalarField2 StreamAreaKinematic(const std::vector<Index2>& sources = std::vector<Index2>()) const;

    ScalarField2 WetnessIndex() const;
    ScalarField2 WetnessIndex(const ScalarField2& sa) const;
    ScalarField2 StreamPower (double m = 1, double n = 1) const;
    ScalarField2 StreamPower (const ScalarField2& sa, double m = 1, double n = 1) const;


    /**********************************/
    /*            OROMETRY            */
    /**********************************/

    // Orometry
    ScalarField2 PeakPercentage(double r) const;
    ScalarField2 ORS(double r) const;
    ScalarField2 JutPlanar(double r) const;
    ScalarField2 RutPlanar(double r) const;
    ScalarField2 JutCurved(double r, double planetR) const;
    ScalarField2 RutCurved(double r, double planetR) const;
    ScalarField2 AngleReducedHeight(double R, int i, int j, bool converging) const;
    ScalarField2 Peakness(ScalarField2& prototypicality, double centroidPercentile = 0.2) const;


protected:
    static const double flat; //!< Small negative epsilon value used in breaching and flow algorithms.

    struct Flow {
        int    i[8] = { -1 };          //!< Indexes to neighboring cells
        Index2 q[8] = { Index2(0,0) }; //!< Flow points
        double s[8] = { 0.0 };         //!< Slopes
        int    n = 0;                  //!< Num outflow neighbors
        int    steepest = -1;          //!< Steepest outflow (0..n-1)
        double slopesum = 0;           //!< Sum of slopes
    };
    Flow GetFlow(const Index2& p) const;

};

inline Box3 HeightField::getBox() const
{
    double za, zb;
    getRange(za, zb);
    return Box3(getDomain(), za, zb);
}

inline double HeightField::Height(int i, int j) const
{
    return at(i, j);
}

inline Vector3 HeightField::Vertex(int i, int j) const
{
  return Vector3(domain.getMin()[0] + i * cellSize[0], domain.getMin()[1] + j * cellSize[1], at(i, j));
}


#endif // HEIGHTFIELD_H
