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

#pragma once

#include <unordered_map>
#include <stdint.h>

#include "GridIndex.hpp"
#include "Shape.hpp"
#include "Types.hpp"

namespace AtlasProcessor
{

struct GridCell
{
    GridCell() : m_shape(-1)
    {}

    std::vector<pdal::PointId> m_points;
    int m_shape;

    void addPoint(pdal::PointId id)
        { m_points.push_back(id); }
    size_t size() const
        { return m_points.size(); }
    bool empty() const
        { return m_points.size() == 0; }
    void setShape(int shape)
        { m_shape = shape; }
    int shape() const
        { return m_shape; }
};


//ABELL - Should just allocate a raster instead of a map-based structure.
class Grid
{
public:
    Grid(pdal::PointViewPtr v, double len, Point origin) :
        m_view(v), m_len(len), m_origin(origin)
    {
        for (pdal::PointId idx = 0; idx < m_view->size(); ++idx)
            insert(idx);
    }

    void insert(pdal::PointId idx)
    {
        using namespace pdal::Dimension;

        double x = m_view->getFieldAs<double>(Id::X, idx) - m_origin.x;
        double y = m_view->getFieldAs<double>(Id::Y, idx) - m_origin.y;
        int xi = std::floor(x / m_len);
        int yi = std::floor(y / m_len);
        m_cells[GridIndex(xi, yi)].addPoint(idx);
    }

    GridCell& cell(int x, int y)
    {
        auto ci = m_cells.find(GridIndex(x, y));
        if (ci == m_cells.end())
            return m_emptyCell;
        return ci->second;
    }

    void findShapes(int windowSize);
    void findShape(Shape& s, GridIndex gi,
        std::unordered_map<GridIndex, GridCell>& cells, int windowSize);
    void draw(int minx, int miny, int maxx, int maxy, bool drawShapes);
    std::vector<Shape>& shapes()
        { return m_shapes; }
    pdal::BOX2D location(const Shape *, Point& center, Point& high);

private:
    GridCell m_emptyCell;
    pdal::PointViewPtr m_view;
    int m_len;
    Point m_origin;
    std::unordered_map<GridIndex, GridCell> m_cells;
    std::vector<Shape> m_shapes;
};
using GridPtr = std::unique_ptr<Grid>;

/**
class GridIter
{
public:
    GridIter(Grid& grid, pdal::Dimension::Id dim);
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = int32_t;
    using value_type = double;
    using pointer = double*;
    using reference = double&;

private:
    Grid& m_grid;
    uint64_t m_pos;
    int m_dimOffset; // Offset to vector iterator.

    int32_t x() const;
    int32_t y() const;

public:
    GridIter& operator++();
    GridIter operator++(int);
    GridIter operator+(const difference_type& n) const;
    GridIter& operator--();
    GridIter operator--(int);
    GridIter operator-(const difference_type& n) const;
    bool operator==(const GridIter& other)
        { return m_pos == other.m_pos; }
    bool operator!=(const GridIter& other)
        { return !(operator==(other)); }
    const double& operator*() const;
    double& operator*();
    const double *operator->() const;
    double *operator->();
};
**/

} // namespace
