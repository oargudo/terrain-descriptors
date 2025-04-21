#ifndef _CORE_H_
#define _CORE_H_

#include <cmath>
#include <vector>
#include <random>
#include <QColor>

namespace Math {
    const double Pi = 3.14159265358979323846;
    const double toDegs = 180.0 / Math::Pi;
    const double toRads = Math::Pi / 180.0;

    inline double sqrt3_2(double a) {
        return std::sqrt(a * a * a);
    }

    inline double CosAtan(double x) {
        return 1.0 / sqrt(1.0 + x * x);
    }

    inline double Mod(double x, double a) {
      return (x >= 0.0) ? fmod(x, a) : a - fmod(-x, a);
    }

    inline double DegreeToRadian(double a) {
        return a * Math::toRads;
    }

    inline double RadianToDegree(double a) {
        return a * Math::toDegs;
    }

    inline double Angle(double a) {
        return Math::Mod(a, 2.0 * Math::Pi);
    }

    inline double LinearStep(double x, double a, double b) {
        if (x < a) return 0;
        if (x > b) return 1;
        return (x - a) / (b - a);
    }

    inline int Clamp(int x, int a, int b) {
        return (x < a ? a : (x > b ? b : x));
    }

    inline double Clamp(double x, double a = 0, double b = 1) {
        return (x < a ? a : (x > b ? b : x));
    }

    template<typename T>
    inline T Lerp(const T& a, const T& b, double t) {
        return a + t * (b - a);
    }

    inline double Bilinear(double a00, double a10, double a11, double a01, double u, double v) {
        return (1 - u)*(1 - v)*a00 + (1 - u)*(v)*a01 + (u)*(1 - v)*a10 + (u)*(v)*a11;
    }

    inline double CubicSmooth(double x, double r) {
        return (1.0 - x / r)*(1.0 - x / r)*(1.0 - x / r);
    }

    inline int Integer(double x) {
        return x > 0.0 ? int(x) : int(x) - 1;
    }

    inline int IntegerSign(double x, double t = 0) {
        if (x < -t) return -1;
        if (x >  t) return 1;
        return 0;
    }

    inline double Mean(const std::vector<double>& v) {
        double s = 0;
        for (double d : v) s += d;
        return s/double(v.size());
    }
}


class Index2 {
public:
    Index2() : i(0), j(0) {};
    Index2(int i, int j) : i(i), j(j) {};
    ~Index2() {};

    int x() const { return i; };
    int y() const { return j; };

    Index2 operator+(const Index2& other) const {
        return Index2(i + other.i, j + other.j);
    }
    Index2 operator-(const Index2& other) const {
        return Index2(i - other.i, j - other.j);
    }
    Index2& operator+=(const Index2& other) {
        i += other.i;
        j += other.j;
        return *this;
    }

    bool operator==(const Index2& other) const {
        return (i == other.i && j == other.j);
    }
    bool operator!=(const Index2& other) const {
        return !(*this == other);
    }
    bool operator<(const Index2& other) const {
        return i < other.i || (i == other.i && j < other.j);
    }

protected:
    int i, j;
};


class IndexArea {
public:
    IndexArea() {
        x0 = x1 = y0 = y1 = 0; };
    IndexArea(int x0, int y0, int x1, int y1) :
        x0(x0), y0(y0), x1(x1), y1(y1) {}
    IndexArea(const Index2& p, int w, int h) :
        x0(p.x()), y0(p.y()), x1(p.x() + w), y1(p.y() + h) {}
    IndexArea(const Index2& p0, const Index2& p1) :
        x0(p0.x()), y0(p0.y()), x1(p1.x()), y1(p1.y()) {}
    ~IndexArea() {}

    int x() const { return x0; }
    int y() const { return y0; }
    //int width() const { return x1 - x0 + 1; }
    //int height() const { return y1 - y0 + 1; }
    int xmin() const { return x0; }
    int xmax() const { return x1; }
    int ymin() const { return y0; }
    int ymax() const { return y1; }

    IndexArea intersected(const IndexArea& other) {
        return IndexArea(std::max(x0, other.x0),
                         std::max(y0, other.y0),
                         std::min(x1, other.x1),
                         std::min(y1, other.y1)
                    );
    }

protected:
    int x0, y0;
    int x1, y1;
};



class Vector3 {
protected:
    double c[3];

public:
    Vector3() {
        c[0] = c[1] = c[2] = 0;
    }
    explicit Vector3(double d) {
        c[0] = c[1] = c[2] = d;
    }
    explicit Vector3(double d0, double d1, double d2) {
        c[0] = d0;
        c[1] = d1;
        c[2] = d2;
    }

    double& operator[] (int i) { return c[i]; };
    double operator[] (int i) const { return c[i]; };
    Vector3 operator-() const { return Vector3(-c[0], -c[1], -c[2]); }

    friend Vector3 operator+(const Vector3& u, const Vector3& v) { return Vector3(u[0]+v[0], u[1]+v[1], u[2]+v[2]); };
    friend Vector3 operator-(const Vector3& u, const Vector3& v) { return Vector3(u[0]-v[0], u[1]-v[1], u[2]-v[2]); };
    friend Vector3 operator*(const Vector3& u, double a) { return Vector3(u[0]*a, u[1]*a, u[2]*a); }
    friend Vector3 operator*(double a, const Vector3& v) { return v * a; }
    friend Vector3 operator/(const Vector3& u, double a) { return Vector3(u[0]/a, u[1]/a, u[2]/a); }
    friend double operator* (const Vector3& u, const Vector3& v)  { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }

    friend double Norm(const Vector3& u) { return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]); }
    friend double SquaredNorm(const Vector3& u) { return u[0]*u[0] + u[1]*u[1] + u[2]*u[2]; }
    friend Vector3 Normalized(const Vector3& u) { return u/Norm(u); }

    friend double dot(const Vector3& u, const Vector3& v) { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }
    friend Vector3 cross(const Vector3& u, const Vector3& v) {
        return Vector3(u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]);
    }

    friend Vector3 Bilinear(const Vector3& a00, const Vector3& a10, const Vector3& a11, const Vector3& a01, const double& u, const double& v) {
        return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
    }
};


class Vector2 {
protected:
    double c[2];

public:
    Vector2() {
        c[0] = c[1] = 0;
    }
    explicit Vector2(double d) {
        c[0] = c[1] = d;
    }
    explicit Vector2(double d0, double d1) {
        c[0] = d0;
        c[1] = d1;
    }
    Vector2(const Vector3& v) {
        c[0] = v[0];
        c[1] = v[1];
    }

    double& operator[] (int i) { return c[i]; };
    double operator[] (int i) const { return c[i]; };

    Vector2 operator- () const { return Vector2(-c[0], -c[1]); };
    friend Vector2 operator+(const Vector2& u, const Vector2& v) { return Vector2(u[0]+v[0], u[1]+v[1]); };
    friend Vector2 operator-(const Vector2& u, const Vector2& v) { return Vector2(u[0]-v[0], u[1]-v[1]); };
    friend Vector2 operator*(const Vector2& u, double a) { return Vector2(u[0]*a, u[1]*a); }
    friend Vector2 operator*(double a, const Vector2& v) { return v * a; }
    friend Vector2 operator/(const Vector2& u, double a) { return Vector2(u[0]/a, u[1]/a); }
    friend double operator* (const Vector2& u, const Vector2& v)  { return u[0]*v[0] + u[1]*v[1]; }

    friend double Norm(const Vector2& u) { return sqrt(u[0]*u[0] + u[1]*u[1]); }
    friend double SquaredNorm(const Vector2& u) { return u[0]*u[0] + u[1]*u[1]; }
    friend Vector2 Normalized(const Vector2& u) { return u/Norm(u); }
    double Angle() const { return std::atan2(c[1], c[0]); }

    Vector3 toVector3(double d) const { return Vector3(c[0], c[1], d); }
};


class Matrix3
{
protected:
  double r[9] = { 1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0 };  //!< The array storing the coefficients of the matrix.

public:
  //! Empty.
  Matrix3() {}
  explicit Matrix3(const double& a00, const double& a01, const double& a02,
                   const double& a10, const double& a11, const double& a12,
                   const double& a20, const double& a21, const double& a22)
  {
    r[0] = a00; r[1] = a01; r[2] = a02;
    r[3] = a10; r[4] = a11; r[5] = a12;
    r[6] = a20; r[7] = a21; r[8] = a22;
  }

  Matrix3 T() const {
    return Matrix3(r[0], r[3], r[6], r[1], r[4], r[7], r[2], r[5], r[8]);
  }

  double Determinant() const {
    return r[0] * r[4] * r[8] + r[1] * r[5] * r[6] + r[2] * r[3] * r[7] - r[2] * r[4] * r[6] - r[1] * r[3] * r[8] - r[0] * r[5] * r[7];
  }

  Matrix3 Adjoint() const {
    return Matrix3(
        r[4] * r[8] - r[7] * r[5], -(r[3] * r[8] - r[6] * r[5]), r[3] * r[7] - r[6] * r[4],
        -(r[1] * r[8] - r[7] * r[2]),  r[0] * r[8] - r[6] * r[2], -(r[0] * r[7] - r[6] * r[1]),
        r[1] * r[5] - r[4] * r[2], -(r[0] * r[5] - r[3] * r[2]), r[0] * r[4] - r[3] * r[1]);
  }

  friend Matrix3 Inverse(const Matrix3&);

  Vector3 operator*(const Vector3& v) const {
    return Vector3(v[0] * r[0] + v[1] * r[3] + v[2] * r[6], v[0] * r[1] + v[1] * r[4] + v[2] * r[7], v[0] * r[2] + v[1] * r[5] + v[2] * r[8]);
  }
};


class Ray
{
protected:
    Vector3 p; // Origin of the ray.
    Vector3 d; // Direction.

public:
    Ray() {}
    explicit Ray(const Vector3& p, const Vector3& d) : p(p), d(d) {}

    Vector3 origin() const { return p; }
    Vector3 direction() const { return d; }

    Vector3 operator()(double t) const { return p + t * d; }

    Ray reflect(const Vector3& p, const Vector3& n) { return Ray(p, n - 2 * n * dot(d, n)); }
};


class Box2 {
protected:
    Vector2 bmin; // min
    Vector2 bmax; // max

public:
    Box2() : bmin(0), bmax(0) {};
    Box2(const double a, const double b) : bmin(-Vector2(a/2, b/2)), bmax(Vector2(a/2,b/2)) {}
    Box2(const Vector2& pmin, const Vector2& pmax) : bmin(pmin), bmax(pmax) {}
    Box2(const Vector2& c, double r) : bmin(c - Vector2(r)), bmax(c + Vector2(r)) {}

    Vector2 getMin() const { return bmin; }
    Vector2 getMax() const { return bmax; }
    Vector2 center() const { return 0.5*(bmin + bmax); }
    double radius() const { return 0.5 * Norm(bmax - bmin); }
    double width() const { return bmax[0] - bmin[0]; }
    double height() const { return bmax[1] - bmin[1]; }

    bool isInside(const Vector2& p) const {
      if ((p[0] < bmin[0]) || (p[0] > bmax[0]) || (p[1] < bmin[1]) || (p[1] > bmax[1]))
        return false;
      return true;
    }
    bool Intersect(const Vector2& s0, const Vector2& s1, double& tmin, double& tmax);
};


class Box3 {
protected:
    Vector3 bmin; // min
    Vector3 bmax; // max

public:
    Box3() : bmin(0), bmax(0) {};
    Box3(const Vector3& pmin, const Vector3& pmax) : bmin(pmin), bmax(pmax) {}
    Box3(const Box2& box2, double zmin, double zmax) {
        bmin = Vector3(box2.getMin()[0], box2.getMin()[1], zmin);
        bmax = Vector3(box2.getMax()[0], box2.getMax()[1], zmax);
    }

    Vector3 getMin() const { return bmin; }
    Vector3 getMax() const { return bmax; }
    Vector3 center() const { return 0.5*(bmin + bmax); }
    Vector3 size() const { return bmax - bmin; }
    double radius() const { return 0.5 * Norm(bmax - bmin); }
    double width() const { return bmax[0] - bmin[0]; }
    double height() const { return bmax[1] - bmin[1]; }
    double depth() const { return bmax[2] - bmin[2]; }

    bool Intersect(const Ray& ray, double& tmin, double& tmax) const;
};


class FieldGrid2D {

public:
    FieldGrid2D() : nx(0), ny(0), domain(Box2()), cellSize(0,0), inverseCellSize(0,0) { }
    FieldGrid2D(const Box2& domain, int nx, int ny) : nx(nx), ny(ny), domain(domain) {
        cellSize = Vector2(domain.width()/(nx-1), domain.height()/(ny-1));
        inverseCellSize = Vector2(1.0/cellSize[0], 1.0/cellSize[1]);
    }
    virtual ~FieldGrid2D() {};

    int getSizeX() const        { return nx; }
    int getSizeY() const        { return ny; }
    int getNumElements() const  { return nx*ny; }
    Vector2 getCellSize() const { return cellSize; }
    double getCellArea() const  { return cellSize[0]*cellSize[1]; }

    constexpr int cellId(int i, int j) const { return i + nx*j; }
    Index2 idToCell(int idx) const { return Index2(idx%nx, idx/nx); }    

    Index2 neighborCell(const Index2& p, int n) const { return p + next[n]; };

    constexpr bool isValidCell(int i, int j) const { return (i >= 0) && (i < nx) && (j >= 0) && (j < ny); }
    bool isValidCell(const Index2& p) const { return isValidCell(p.x(), p.y()); }
    bool isInDomain(const Vector2& p) const {
        return p[0] >= domain.getMin()[0] && p[1] >= domain.getMin()[1] &&
               p[0] <= domain.getMax()[0] && p[1] <= domain.getMax()[1];
    }

    Box2 getDomain() const { return domain; };
    Vector2 domainCoords(int i, int j) const { return domain.getMin() + Vector2(i * cellSize[0], j * cellSize[1]); }
    Vector2 domainCoords(const Index2& idx) const { return domainCoords(idx.x(), idx.y()); }

    void cellCoords(const Vector2& p, int& i, int& j, double& u, double& v) const {
        Vector2 q = p - domain.getMin();
        u = q[0]*inverseCellSize[0];
        v = q[1]*inverseCellSize[1];
        // Integer coordinates
        i = int(u);
        j = int(v);
        // Local coordinates within cell
        u -= i;
        v -= j;
    }
    void cellCoords(const Vector2& p, int& i, int& j) const {
        Vector2 q = p - domain.getMin();
        i = int(q[0]*inverseCellSize[0]);
        j = int(q[1]*inverseCellSize[1]);
    }
    Index2 cellCoords(const Vector2& p) const {
        int i, j;
        cellCoords(p, i, j);
        return Index2(i, j);
    }


    IndexArea indexAreaFromBox(const Box2& box) const
    {
        // Area from indices of box corners
        Index2 pa = cellCoords(box.getMin());
        Index2 pb = cellCoords(box.getMax());
        IndexArea area(pa, pb);
        // Limit to domain
        IndexArea mask(0, 0, nx - 1, ny - 1);
        return area.intersected(mask);
    }

    IndexArea indexAreaFromHalfWindow(const Index2 &ij, int w_2) const
    {
        // Rectangle
        IndexArea area(ij - Index2(w_2,w_2), ij + Index2(w_2, w_2));
        // Limit to domain
        IndexArea mask(0, 0, nx - 1, ny - 1);
        return area.intersected(mask);
    }

    IndexArea indexAreaFromWindow(int i, int j, int w) const
    {
        return indexAreaFromHalfWindow(Index2(i, j), w/2);
    }

    IndexArea indexAreaFromRadius(int i, int j, double r) const
    {
        Vector2 center = domainCoords(i, j);
        Box2 box(center - Vector2(r, r), center + Vector2(r,r));
        return indexAreaFromBox(box);
    }


protected:
    int nx, ny;
    Box2 domain;
    Vector2 cellSize;
    Vector2 inverseCellSize;

public:
    static const Index2 next[8];          //!< Array of points in the 1-ring neighborhood.
    static const double length[8];        //!< Length to the i-th neighbor.
    static const double inverselength[8]; //!< Inverse length.
};


class Triangle
{
protected:
    Vector3 p[3]; //!< Array of vertices.

public:
    Triangle() {}
    explicit Triangle(const Vector3& a, const Vector3& b, const Vector3& c) {
        p[0] = a; p[1] = b; p[2] = c;
    };
    ~Triangle() {}

    Vector3 operator[](int i) const { return p[i]; }

    Vector3 Center() const { return (p[0] + p[1] + p[2]) / 3.0; }
    Vector3 Normal() const { return Normalized(cross(p[1] - p[0], p[2] - p[0])); }
    Vector3 AreaNormal() const { return 0.5 * cross(p[1] - p[0], p[2] - p[0]); }
    double Area() const { return 0.5 * Norm(cross(p[0] - p[1], p[2] - p[0])); }
};


class QuadricSurface {
// z=a22  x<SUP>2</SUP>y<sup>2</sup> + a21 x<sup>2</sup>y +  a21 xy<sup>2</sup> + a20 x<sup>2</sup> + a02 y<sup>2</sup> +  a11 xy +  a10 x + a01 y + a00.
protected:
    double c[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; //!< %Quadric coefficients.
public:
    QuadricSurface() {}
    explicit QuadricSurface(const double& a22, const double& a21, const double& a12, const double& a20, const double& a02,
                            const double& a11, const double& a10, const double& a01, const double& a00)
    {
        c[0] = a22;
        c[1] = a21;
        c[2] = a12;
        c[3] = a20;
        c[4] = a02;
        c[5] = a11;
        c[6] = a10;
        c[7] = a01;
        c[8] = a00;
    }
    double operator()(int i, int j) const { return c[Index(i, j)]; }
protected:
    constexpr int Index(int a, int b) const { return terms[a + 3 * b]; }
    static constexpr int terms[9] = { 8, 6, 3, 7, 5, 2, 4, 1, 0 };
};



class Camera {
protected:
    Vector3 eye;      // Eye
    Vector3 at;       // Look at point
    Vector3 up;       // Up vector
    double cah;       // Camera aperture horizontal
    double cav;       // Camera aperture vertical
    double nearplane; // Near plane
    double farplane;  // Far plane
    double fl;        // Focal length

public:
    Camera();
    Camera(const Vector3& eye, const Vector3& at, const Vector3& up = Vector3(0,0,1), double near = 1.0, double far = 100000.0);

    Vector3 getEye() const { return eye; }
    Vector3 getAt() const { return at; }
    Vector3 getUp() const { return up; }
    Vector3 getViewDir() const { return Normalized(at - eye); }
    double getNearPlane() const { return nearplane; }
    double getFarPlane() const { return farplane; }
    double getAngleOfViewH(double, double) const;
    double getAngleOfViewV(double, double) const;

    void setAt(const Vector3& p) { at = p; up = Vector3(0,0,1); }
    void setEye(const Vector3& p) { eye = p; }
    void setPlanes(double n, double f) { nearplane = n; farplane = f;}

    // Move camera around eye
    void upDownRound(double a);
    void leftRightRound(double a);
    void backForth(double a, bool moveAt = false);

    // Move camera in a plane
    void upDownPlane(double);
    void leftRightPlane(double);

    Ray pixelToRay(int px, int py, int w, int h) const;

    static Camera View(const Box3& box);
};


inline Vector3 fromQColor(const QColor& c) {
    return Vector3(c.red() / 255.0, c.green() / 255.0, c.blue() / 255.0);
}

inline QColor toQColor(const Vector3& c) {
    return QColor(int(255.0 * Math::Clamp(c[0])),
                  int(255.0 * Math::Clamp(c[1])),
                  int(255.0 * Math::Clamp(c[2])),
                  255);
}

class ColorPalette
{
protected:
    std::vector<Vector3> colors;
    std::vector<double> anchors;

public:
    ColorPalette() : colors({Vector3(1)}), anchors({0}) {}
    ColorPalette(const std::vector<Vector3>& c, const std::vector<double>& a) : colors(c), anchors(a) {}

    Vector3 getColor(double u) const;

    static ColorPalette CoolWarm() {
        const Vector3 Cool = Vector3(97, 130, 234) / 255.0;
        const Vector3 White = Vector3(221, 221, 221) / 255.0;
        const Vector3 Warm = Vector3(220, 94, 75) / 255.0;
        static const ColorPalette p = ColorPalette({Cool, White, Warm}, {0, 0.5, 1});
        return p;
    }

    static ColorPalette Reds() {
        static const ColorPalette p = ColorPalette({
            Vector3(255, 245, 240) / 255.0,
            Vector3(254, 224, 210) / 255.0,
            Vector3(252, 187, 161) / 255.0,
            Vector3(252, 146, 114) / 255.0,
            Vector3(251, 106, 74) / 255.0,
            Vector3(239, 59, 44) / 255.0,
            Vector3(203, 24, 29) / 255.0,
            Vector3(153, 0, 13) / 255.0
        }, {
            0.0, 0.1429, 0.2857, 0.4286, 0.5714, 0.7143, 0.8571, 1.0
        });
        return p;
    }

    static ColorPalette Blues() {
        static const ColorPalette p = ColorPalette({
            Vector3(247, 251, 255) / 255.0,
            Vector3(222, 235, 247) / 255.0,
            Vector3(198, 219, 239) / 255.0,
            Vector3(158, 202, 225) / 255.0,
            Vector3(107, 174, 214) / 255.0,
            Vector3(66, 146, 198) / 255.0,
            Vector3(33, 113, 181) / 255.0,
            Vector3(8, 81, 156) / 255.0
        }, {
            0.0, 0.1429, 0.2857, 0.4286, 0.5714, 0.7143, 0.8571, 1.0
        });
        return p;
    }

    static ColorPalette Relief() {
        const std::vector<Vector3> c = { Vector3(0.627, 0.863, 0.411),
                                         Vector3(1.00, 0.90, 0.45),
                                         Vector3(0.659, 0.607, 0.541),
                                         Vector3(0.95, 0.95, 0.95) };
        static const ColorPalette p = ColorPalette(c, {0, 0.375, 0.625, 1.0});
        return p;
    }

};


class LookupPalette
{
protected:
    std::vector<Vector3> c; //!< Set of colors.
public:
    LookupPalette(const std::vector<Vector3>& v) : c(v) {};
    virtual Vector3 getColor(int) const;
};


class Sun
{
public:
    static void Angles(const double&, const double&, int, int, int, const double&, double&, double&);
protected:
    static double Day(int, int, int, double);
    static double f0(double, double);
    static double FNsun(double,double&,double&,double&);
};


class Random {
public:
    static std::minstd_rand rng;
    static std::uniform_real_distribution<double> uniformDist;

    static inline double Uniform(double a = 0, double b = 1) {
        return a + (b-a)*uniformDist(rng);
    }

    static inline Vector3 SampleOnSphere() {
        double z = 2.0 * Uniform() - 1.0;
        double w = sqrt(1.0 - z * z);
        double a = Uniform(0, 2*Math::Pi);
        return Vector3(w*std::cos(a), w*std::sin(a), z);
    }

    static inline Vector3 SampleOnHemisphere(const Vector3& n) {
        Vector3 d = SampleOnSphere();
        if (d*n < 0) {
            d = -d;
        }
        return d;
    }
};



#endif
