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

#include <array>
#include <utility>

#include <pdal/PointView.hpp>

namespace AtlasProcessor
{

const int XCellCount = 200;
const int YCellCount = 100;

struct Point
{
    Point() : x(0), y(0), z(0)
    {}

    Point(double x, double y) : x(x), y(y), z(0)
    {}

    Point(double x, double y, double z) : x(x), y(y), z(z)
    {}

    bool operator==(const Point& p)
        { return p.x == x && p.y == y && p.z == z; }

    bool operator!=(const Point& p)
        { return p.x != x || p.y != y || p.z != z; }

    double x;
    double y;
    double z;
};

inline std::istream& operator>>(std::istream& in, Point& p)
{
    std::string s;
    std::vector<std::string> l;

    in >> s;
    l = pdal::Utils::split(s, ',');
    if (l.size() != 2 && l.size() != 3)
        in.setstate(std::ios_base::failbit);
    else
    {
        for (std::string& t : l)
            pdal::Utils::trim(t);
        p.x = std::stoi(l[0]);
        p.y = std::stoi(l[1]);
        p.z = (l.size() == 3) ? p.z : 0;
    }
    return in;
}

inline std::ostream& operator<<(std::ostream& out, const Point& p)
{
    out << p.x << ", " << p.y << ", " << p.z;
    return out;
}

using Histogram = std::array<pdal::PointViewPtr, 10>;
using Coord = std::pair<int, int>;

} // namespace AtlasProcessor

namespace AP = AtlasProcessor;
