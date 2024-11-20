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

#include "Types.hpp"

namespace AtlasProcessor
{

class Field
{
public:
    Field(size_t width, size_t height) : m_width(width), m_height(height),
        m_x(width * height), m_y(width * height), m_z(width * height),
        m_before(width * height), m_after(width * height), m_valid(width * height),
        m_maxLen2(std::numeric_limits<double>::lowest())
    {}

    const float *xdata() const
    { return m_x.data(); }

    const float *ydata() const
    { return m_y.data(); }

    const float *zdata() const
    { return m_z.data(); }

    const float *bdata() const
    { return m_before.data(); }

    const float *adata() const
    { return m_after.data(); }

    size_t width() const
    { return m_width; }

    size_t height() const
    { return m_height; }

    double maxLen2() const
    { return m_maxLen2; }

    void setOffset(Coord c, Point displacement)
    {
        size_t idx = pos(c);
        if (idx < 0)
            return;

        double len2 = displacement.x * displacement.x + displacement.y * displacement.y;
        m_maxLen2 = std::max(len2, m_maxLen2);
        m_x[idx] = displacement.x;
        m_y[idx] = displacement.y;
        m_z[idx] = displacement.z;
        m_valid[idx] = true;
    }

    void setBeforeCount(Coord c, float count)
    {
        size_t idx = pos(c);
        if (idx < 0)
            return;

        m_before[idx] = count;
    }

    void setAfterCount(Coord c, float count)
    {
        size_t idx = pos(c);
        if (idx < 0)
            return;

        m_after[idx] = count;
    }

    Point offset(Coord c)
    {
        size_t idx = pos(c);
        if (idx < 0)
            return Point();
        return Point(m_x[idx], m_y[idx], m_z[idx]);
    }

    bool valid(Coord c)
    { return pos(c) >= 0; }

    std::pair<bool, Point> initialOffset(Coord c)
    {
        const double Sqrt2Recip = 0.70710678118;
        double x = 0;
        double y = 0;
        double sum = 0;

        auto addit = [&x, &y, &sum, this](Coord c, double rel)
        {
            int idx = pos(c);
            if (idx >= 0 && m_valid[idx])
            {
                x += (m_x[idx] * rel);
                y += (m_y[idx] * rel);
                sum += rel;
            }
        };

        int idx = pos(c);
        if (idx >= 0)
        {
            if (m_valid[idx])
                return { true, { m_x[idx], m_y[idx] } };

            // We weight each full neighbor equally and each corner neighbor
            // by 1/sqrt(2)
            addit(Coord(c.first + 1, c.second), 1);
            addit(Coord(c.first - 1, c.second), 1);
            addit(Coord(c.first, c.second + 1), 1);
            addit(Coord(c.first, c.second - 1), 1);
            addit(Coord(c.first + 1, c.second + 1), Sqrt2Recip);
            addit(Coord(c.first - 1, c.second + 1), Sqrt2Recip);
            addit(Coord(c.first + 1, c.second - 1), Sqrt2Recip);
            addit(Coord(c.first - 1, c.second - 1), Sqrt2Recip);
            if (sum)
            {
                x /= sum;
                y /= sum;
                return { true, { x, y } };
            }
        }
        return { false, { 0, 0 } };
    }

private:
    size_t m_width;
    size_t m_height;

    // Easy output to GDAL
    std::vector<float> m_x;
    std::vector<float> m_y;
    std::vector<float> m_z;
    std::vector<float> m_before;
    std::vector<float> m_after;
    std::vector<bool> m_valid;
    double m_maxLen2;

    int pos(Coord c)
    {
        int x = c.first;
        int y = c.second;

        if (x < 0 || x >= m_width || y < 0 || y >= m_height)
            return -1;
        return (y * m_width) + x;
    }
};

}
