/*****************************************************************************
 *   Copyright (c) 2024, Hobu, Inc. (info@hobu.co)                           *
 *                                                                           *
 *   All rights reserved.                                                    *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 ****************************************************************************/

#include <stack>

#include "Grid.hpp"

namespace AtlasProcessor
{

using namespace pdal;


void Grid::findShapes(int windowSize)
{
    std::unordered_map<GridIndex, GridCell> cells(m_cells);

    size_t shape = 0;
    while (cells.size())
    {
        auto it = cells.begin();
        GridIndex gi = it->first;
        cells.erase(it);

        Shape s(shape++);
        findShape(s, gi, cells, windowSize);
        m_shapes.push_back(s);
    }

    for (Shape& s : m_shapes)
        for (const GridIndex& gi : s.indices())
            m_cells[gi].setShape(s.id());
}


void Grid::findShape(Shape& s, GridIndex gi,
    std::unordered_map<GridIndex, GridCell>& cells, int windowSize)
{
    std::stack<GridIndex> pending;

    pending.push(gi);

    // This is woefully inefficient.  Not sure I care right now.
    while (pending.size())
    {
        gi = pending.top();
        pending.pop();
        s.push(gi);
        int minx = gi.x() - windowSize;
        int maxx = gi.x() + windowSize;
        int miny = gi.y() - windowSize;
        int maxy = gi.y() + windowSize;

        for (int y = miny; y <= maxy; ++y)
            for (int x = minx; x <= maxx; ++x)
            {
                auto ci = cells.find({x, y});
                if (ci != cells.end())
                {
                    pending.push(ci->first);
                    cells.erase(ci);
                }
            }
    }
}


void Grid::draw(int minx, int miny, int maxx, int maxy, bool shapes)
{
    for (int y = maxy; y >= miny; --y)
    {
        for (int x = minx; x <= maxx; ++x)
        {
            char c;
            if (shapes)
            {
                int s = cell(x, y).shape();
                if (s < 0)
                    c = '.';
                if (s >= 0)
                    c = '0' + (s % 10);

            }
            else
            {
                size_t count = cell(x, y).size();
                if (count == 0)
                    c = '.';
                else if (count > 9)
                    c = '*';
                else
                    c = '0' + count;
            }
            std::cout << c;
        }
        std::cout << " " << y << "\n";
    }
    std::cout << "\n";
}


BOX2D Grid::location(const Shape *s, Point& center, Point& high)
{
    center = {0, 0};
    high = {0, 0};

    BOX2D box;
    size_t count = 0;
    double zmax = 0;
    for (const GridIndex& gi : s->indices())
    {
        GridCell& cell = m_cells[gi];
        for (PointId id : cell.m_points)
        {
            double x = m_view->getFieldAs<double>(Dimension::Id::X, id);
            double y = m_view->getFieldAs<double>(Dimension::Id::Y, id);
            double z = m_view->getFieldAs<double>(Dimension::Id::Z, id);

            box.grow(x, y);
            center.x += x;
            center.y += y;
            center.z += z;
            if (z > zmax)
            {
                high.x = x;
                high.y = y;
            }
            count++;
        }
    }
    center.x /= count;
    center.y /= count;
    center.z /= count;
    return box;
}

/***
//
// GridIter
//

GridIter::GridIter(Grid& grid, pdal::Dimension::Id dim) :
    m_grid(grid), m_pos(0), m_dimOffset(-1)
{
    using namespace pdal::Dimension;

    if (dim != Id::X && dim != Id::Y && dim != Id::Z &&
            dim != Id::Unknown)
//ABELL    throwError("Dimension must be X, Y, Z or Unknown.")
        std::cerr << "Dimension must be X, Y, Z or Unknown.\n";
    m_dimOffset = static_cast<int>(dim) - 1;
}

int32_t GridIter::x() const
{
    return (m_pos % m_grid.ySize()) + m_grid.xOrigin();
}

int32_t GridIter::y() const
{
    return (m_pos / m_grid.ySize()) + m_grid.yOrigin();
}

GridIter& GridIter::operator++()
{
    m_pos++;
    return *this;
}


GridIter GridIter::operator++(int)
{
    GridIter it(*this);
    m_pos++;
    return it;
}

GridIter& GridIter::operator--()
{
    m_pos--;
    return *this;
}

GridIter GridIter::operator--(int)
{
    GridIter it(*this);
    m_pos--;
    return it;
}

GridIter GridIter::operator+(const difference_type& n) const
{
    GridIter it(*this);
    it.m_pos = m_pos + n;
    return it;
}


GridIter GridIter::operator-(const difference_type& n) const
{
    GridIter it(*this);
    it.m_pos = m_pos - n;
    return it;
}

const double& GridIter::operator*() const
{
    static double empty(-9999);

    Eigen::Vector3d *vec = m_grid.getVector(x(), y());
    if (!vec)
        return empty;
    return (*vec)(m_dimOffset);
}

double& GridIter::operator*()
{
    const double& d = const_cast<const GridIter *>(this)->operator *();
    return const_cast<double&>(d);
}

const double *GridIter::operator->() const
{
    static double empty(-9999);

    Eigen::Vector3d *vec = m_grid.getVector(x(), y());
    if (!vec)
        return &empty;
    return &(*vec)(m_dimOffset);
}

double *GridIter::operator->()
{
    const double& d = const_cast<const GridIter *>(this)->operator *();
    double& dd = const_cast<double &>(d);
    return &dd;
}
***/

} // namespace
