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

#include <fstream>

#include "Types.hpp"

namespace AtlasProcessor
{

class Draw
{
private:
    std::ofstream m_out;
    double m_width;
    double m_height;
    int m_xCells;
    int m_yCells;
    double m_maxLen;
    double m_userWidth;
    double m_userHeight;

public:
    Draw(double width, double height, int xCells, int yCells,
        double maxLen = 7);
    ~Draw();

    void doVector(const Coord& cell, const Point& vec);

private:
    Point center(const Coord& c);
    void label(const Coord& c);
    void dot(const Coord& c);
    void arrow(int len, const Coord& c, const Point& vec);
    void arrowhead(const Coord& c, const Point& vec);
    void arrowhead(const Point& vec, const Point& offset);
};

}
