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

#include "Types.hpp"
#include "GridIndex.hpp"

namespace AtlasProcessor
{

struct Shape
{
public:
    Shape(size_t id) : m_id(id),
        m_xmin(std::numeric_limits<int>::max()),
        m_xmax(std::numeric_limits<int>::min()),
        m_ymin(std::numeric_limits<int>::max()),
        m_ymax(std::numeric_limits<int>::min()),
        m_xCenter(0), m_yCenter(0)
    {}

    void push(const GridIndex& gi)
    {
        m_indices.push_back(gi);
        m_xmin = (std::min)(gi.x(), m_xmin);
        m_xmax = (std::max)(gi.x(), m_xmax);
        m_ymin = (std::min)(gi.y(), m_ymin);
        m_ymax = (std::max)(gi.y(), m_ymax);
        m_xCenter += ((gi.x() - m_xCenter) / m_indices.size());
        m_yCenter += ((gi.y() - m_yCenter) / m_indices.size());
    }
    size_t id() const
        { return m_id; }
    const std::vector<GridIndex>& indices() const
        { return m_indices; }
    size_t size() const
        { return m_indices.size(); }
    Coord center() const
        { return {(int)std::round(m_xCenter), (int)std::round(m_yCenter)}; }
    Point exactCenter() const
        { return {m_xCenter, m_yCenter}; }
    int xmin() const
        { return m_xmin; }
    int xmax() const
        { return m_xmax; }
    int ymin() const
        { return m_ymin; }
    int ymax() const
        { return m_ymax; }

private:
    size_t m_id;
    std::vector<GridIndex> m_indices;
    int m_xmin;
    int m_xmax;
    int m_ymin;
    int m_ymax;
    double m_xCenter;
    double m_yCenter;
};
using ShapePair = std::pair<const Shape *, const Shape *>;

} // namespace AtlasProcessor
