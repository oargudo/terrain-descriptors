#include <vector>
#include <queue>
#include <tuple>
#include <functional>

#include "heightfield.h"

// Neighbour Directions
const int DIRECTIONS = 8;
const int dx[DIRECTIONS] = { -1, -1, 0, 1, 1, 1, 0, -1 };
const int dy[DIRECTIONS] = { 0, -1, -1, -1, 0, 1, 1, 1 };

// State of a vertex during the execution of the Priority-Flood algorithm.
enum State
{
  CLOSED,
  OPEN
};

// Element in the priority queue during the execution of the Priority-Flood algorithm.
struct Element
{
  // Priority
  double elevation;
  // Coordinates
  int x;
  int y;

  Element() : elevation(-1.0), x(-1), y(-1) {}

  Element(int x, int y, double elevation) : elevation(elevation), x(x), y(y) {}

  // Compare elevation, then x, then y
  bool operator>(const Element& rhs) const
  {
    return std::tie(elevation, x, y) > std::tie(rhs.elevation, rhs.x, rhs.y);
  }

  typedef std::queue<Element> Queue;
  typedef std::priority_queue<Element, std::vector<Element>, std::greater<Element> > PriorityQueue;
};

// Init the state of all cells on the edges of the HeightField
void InitBorders(const HeightField* height_field, Element::PriorityQueue& open, std::vector<std::vector<State> >& state)
{
  const int size_x = height_field->getSizeX();
  const int size_y = height_field->getSizeY();

  for (int x = 0; x < size_x; x++)
  {
    open.emplace(x, 0, height_field->at(x, 0));
    state[x][0] = CLOSED;

    open.emplace(x, size_y - 1, height_field->at(x, size_y - 1));
    state[x][size_y - 1] = CLOSED;
  }
  for (int y = 1; y < size_y - 1; y++)
  {
    open.emplace(0, y, height_field->at(0, y));
    state[0][y] = CLOSED;

    open.emplace(size_x - 1, y, height_field->at(size_x - 1, y));
    state[size_x - 1][y] = CLOSED;
  }
}


/*!
\brief Modifies the HeightField to guarantee drainage.

\sa FillDepressions(const double&, ScalarField2&) const
\author Mathieu Gaillard
\param eps The minimal elevation difference between two cells.
\return The number of cells where a significant alteration of the DEM has occurred.
*/
int HeightField::fillDepressions(const double& eps)
{
  // Store the depressions in a separate ScalarField2
  ScalarField2 depression_field(getDomain(), nx, ny, 0.0);
  int nb_false_pit_cell = fillDepressions(eps, depression_field);
  // Add the depressions to this HeightField
  operator+=(depression_field);
  return nb_false_pit_cell;
}


/*!
\brief Modifies the HeightField to guarantee drainage.

Algorithm:
Priority-flood: An optimal depression-filling and watershed-labeling
algorithm for digital elevation models.
Barnes, R., Lehman, C., Mulla, D., 2014.

This version of Priority-Flood starts on the edges of the DEM and then
works its way inwards using a priority queue to determine the lowest cell
which has a path to the edge. The neighbours of this cell are added to the
priority queue if they are higher. If they are lower, then their elevation
is increased by a small amount to ensure that they have a drainage path and
they are added to a "pit" queue which is used to flood pits. Cells which
are higher than a pit being filled are added to the priority queue. In this
way, pits are filled without incurring the expense of the priority queue.

The values to fill the depressions are added to the depression field.

The depression field should have the exact same size as the HeightField.
If the given depression field does not have the exact same size as the height field,
returns the total number of vertices in the heightfield.

To be sure that the height field is never flat, elevations are always increased
by std::nextafter(elevation) + eps. Thus, even if eps == 0, the terrain is never flat.

Warning, because the method is const it is not recommended to use it like that
\code
// Not recommended, although it works.
hf.FillDepressions(0.0, hf);
\endcode

\author Mathieu Gaillard

\param eps The minimal elevation difference between two cells.
\param depression_field A scalar field in which the values to fill depressions are added
\return The number of cells where a significant alteration of the DEM has occurred.
*/
int HeightField::fillDepressions(const double& eps, ScalarField2& depression_field) const
{
  // Heightfield properties
  const int size_x = getSizeX();
  const int size_y = getSizeY();

  // Check ScalarField2 size
  if (depression_field.getSizeX() != size_x || depression_field.getSizeY() != size_y)
  {
    return getNumElements();
  }

  // State of vertices
  std::vector<std::vector<State> > state(size_x, std::vector<State>(size_y, OPEN));

  // A min priority queue of vertices to visit, ordered by elevation
  Element::PriorityQueue open;

  // A queue of vertices to visit in priority, because they have the same elevation as the current element
  Element::Queue pit;

  // True if currently filling a pit
  bool currently_in_pit = false;

  // Elevation of the first vertex of the pit
  double pit_top = 0;

  // The number of cells where a significant alteration of the DEM has occurred
  int nb_false_pit_cell = 0;

  // Init the state of all cells on the edges of the HeightField
  InitBorders(this, open, state);

  // Perform algorithm
  while (!open.empty() || !pit.empty())
  {
    Element curr_vertex;

    // Maintains the total order by ensuring that all cells of equal elevation
    // which may border a depression are treated equally. Otherwise, in priority, 
    // take the next element from the pit queue.
    if ((!open.empty() && !pit.empty() && open.top().elevation == pit.front().elevation) || pit.empty())
    {
      curr_vertex = open.top();
      open.pop();
      currently_in_pit = false;
    }
    else if (!pit.empty())
    {
      curr_vertex = pit.front();
      pit.pop();
      if (currently_in_pit == false)
      {
        currently_in_pit = true;
        pit_top = curr_vertex.elevation;
      }
    }

    // For all neighbors of the current element
    for (int d = 0; d < DIRECTIONS; d++)
    {
      int neighbor_x = curr_vertex.x + dx[d];
      int neighbor_y = curr_vertex.y + dy[d];

      // If the neighbor exists and is not closed
      if (neighbor_x >= 0 && neighbor_x < size_x && neighbor_y >= 0 && neighbor_y < size_y && state[neighbor_x][neighbor_y] == OPEN)
      {
        state[neighbor_x][neighbor_y] = CLOSED;

        // Current elevation plus an epsilon
        double curr_elevation_plus_eps = nextafter(curr_vertex.elevation, std::numeric_limits<double>::infinity()) + eps;

        if (at(neighbor_x, neighbor_y) <= curr_elevation_plus_eps)
        {
          if (currently_in_pit == true && pit_top < at(neighbor_x, neighbor_y))
          {
            // A significant alteration of the DEM has occurred
            // The inside of the pit is now higher than the terrain surrounding it
            nb_false_pit_cell++;
          }

          depression_field(neighbor_x, neighbor_y) += curr_elevation_plus_eps - at(neighbor_x, neighbor_y);
          pit.emplace(neighbor_x, neighbor_y, curr_elevation_plus_eps);
        }
        else
        {
          open.emplace(neighbor_x, neighbor_y, at(neighbor_x, neighbor_y));
        }
      }
    }
  }

  return nb_false_pit_cell;
}
