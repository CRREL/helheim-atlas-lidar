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

#include "Draw.hpp"

namespace AtlasProcessor
{

namespace
{
    int CellSide = 10;
};

Draw::Draw(double width, double height, int xCells, int yCells, double maxLen) :
    m_out("out.svg"),
    m_width(width), m_height(height), m_xCells(xCells), m_yCells(yCells),
    m_maxLen(maxLen),
    m_userWidth(xCells * CellSide), m_userHeight(yCells * CellSide)
{
    m_out << "<svg xmlns=\"http://www.w3.org/2000/svg\"\n"
        "\tversion=\"1.1\"\n"
        "\twidth=\"" << width << "\" height=\"" << m_height << "\" " <<
        "\tviewBox=\"0 0 " << m_userWidth << " " << m_userHeight << "\"> \n";
}

Draw::~Draw()
{
    m_out << "</svg>\n";
}


Point Draw::center(const Coord& c)
{
    Point p;
    p.x = CellSide * (.5 + c.first);

    // Transform y up to y down.
    p.y = CellSide * (.5 + m_yCells - c.second);
    return p;
}


void Draw::doVector(const Coord& cell, const Point& vec)
{
    double mag2 = (vec.x * vec.x) + (vec.y * vec.y);
    double mag = std::sqrt(mag2);

    // Normalize
    Point pHat { vec.x / mag, vec.y / mag };

    int len = std::min<int>(6, (mag / m_maxLen) * 7);
    label(cell);
    switch (len)
    {
    case 0:
        dot(cell);
        break;
    case 1:
        arrowhead(cell, pHat);
        break;
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
        arrow(len - 2, cell, pHat);
        break;
    }
}

void Draw::label(const Coord& c)
{
    if (c.first % 5 || c.second % 5)
        return;

    std::string coord = std::to_string(c.first) + "/" + std::to_string(m_yCells - c.second - 1);
    Point p = center(c);
    m_out << "<g transform=\"translate(" <<
            p.x << ',' << p.y << ")\">\n";
        m_out << "\t<text x=\"0\" y=\"0\">" << coord << "</text>\n";
    m_out << "</g>\n";

}

void Draw::dot(const Coord& c)
{
    Point p = center(c);

    m_out << "<g transform=\"translate(" <<
            p.x << ',' << p.y << ")\">\n";
        m_out << "\t<circle cx=\"0\" cy=\"0\" r=\"1\" /> \n";
    m_out << "</g>\n";
}


void Draw::arrow(int len, const Coord& c, const Point& vec)
{
    Point arrowOffset[] = { {1.5, 0}, {2.25, 0}, {3, 0}, {3.75, 0}, {4.5, 0} };
    Point p = center(c);

    double x1 = -1.5 + (.75 * -len);
    double x2 = 1.5 + (.75 * len);

    m_out << "<g transform=\"translate(" <<
        p.x << ',' << p.y << ")\">\n";
    m_out << "<g transform=\"matrix(" <<
        vec.x << ',' << vec.y << "," <<
        -vec.y << "," << vec.x << ',' << "0,0)\">";
        m_out << "\t<line x1=\"" << x1 << "\" x2=\"" << x2 << "\" "
            "stroke=\"red\" ";
        m_out << "/>\n";
        arrowhead(vec, arrowOffset[len]);
    m_out << "</g>\n";
    m_out << "</g>\n";
}


void Draw::arrowhead(const Coord& c, const Point& vec)
{
    Point p = center(c);

    m_out << "<g transform=\"translate(" <<
        p.x << ',' << p.y << ")\">\n ";
    m_out << "<g transform=\"matrix(" << vec.x << ',' << vec.y <<
        ',' << -vec.y << ',' << vec.x << ",0,0)\">\n ";
        arrowhead(vec, Point{0,0});
    m_out << "</g>\n";
    m_out << "</g>\n";
}


void Draw::arrowhead(const Point& vec, const Point& off)
{
    m_out << "\t<polygon ";
    m_out << "stroke=\"red\" fill=\"red\" ";
    m_out << "points=\"-.6,0 -.8,1.4 .5,0 -.8,-1.4";
    m_out << "\" ";

    m_out << "transform="
        "\"translate(" << off.x << ',' << off.y << ") "
        "\"/>\n";
}

} // namespace AtlasProcessor

