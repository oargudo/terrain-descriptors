#ifndef INTEGERFIELD_H
#define INTEGERFIELD_H

#include "core.h"
#include <QImage>

class ScalarField2;

class IntField2 : public FieldGrid2D
{
public:
    IntField2();
    IntField2(const Box2&, int, int, int = 0);
    IntField2(const Box2&, int, int, const std::vector<int>&);
    virtual ~IntField2() {};

    // Access to elements
    int value(int, int) const;
    int at(int) const;
    int at(int, int) const;
    int at(const Index2&) const;
    int& operator()(int, int);
    int& operator()(const Index2&);
    int& operator[](int);

    void getRange(int&, int&) const;
    int percentile(double p) const;

    void fill(int);

    QImage createImage(const LookupPalette&) const;

    ScalarField2 toScalarField() const;

protected:
    std::vector<int> field;
};


/*!
\brief Return the field value at a given array vertex.
\param i,j Integer coordinates of the vertex.
\sa at(int,int)
*/
inline int IntField2::value(int i, int j) const
{
    return field.at(cellId(i, j));
}

/*!
\brief Return the field value at a given array vertex.
\param i,j Integer coordinates of the vertex.
*/
inline int IntField2::at(int i, int j) const
{
  return field.at(cellId(i, j));
}

/*!
\brief Return the field value at a given array vertex.
\param q Point.
*/
inline int IntField2::at(const Index2& q) const
{
  return field.at(cellId(q.x(), q.y()));
}

/*!
\brief Return the field value at a given array vertex.
\param q Point.
*/
inline int& IntField2::operator()(const Index2& q)
{
  return field[cellId(q.x(), q.y())];
}

/*!
\brief Return the field value at a given array vertex.
\param i,j Integer coordinates of the vertex.
*/
inline int& IntField2::operator()(int i, int j)
{
  return field[cellId(i, j)];
}

/*!
\brief Return the data in the field.
\param c Index.
*/
inline int IntField2::at(int c) const
{
  return field.at(c);
}

/*!
\brief Return the data in the field.
\param c Index.
*/
inline int& IntField2::operator[](int c)
{
  return field[c];
}


#endif
