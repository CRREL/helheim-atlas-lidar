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

#include <pdal/PipelineManager.hpp>
#include <pdal/util/ProgramArgs.hpp>

#include "Field.hpp"
#include "Grid.hpp"

namespace AtlasProcessor
{

class Atlas
{
public:
    Atlas();
    void run(const pdal::StringList& s);

private:
    Histogram histogram(pdal::PointViewPtr v, pdal::Dimension::Id dim,
        pdal::PointViewPtr debugView, const std::string& label);
    bool removeFliers(pdal::PointViewPtr& v, const Histogram& hist);
    GridPtr buildGrid(pdal::PointViewPtr v, Point origin);
    void sortShapes(GridPtr& g);
    void dumpShapes(GridPtr& g);
    void dumpSurrounding();
    Coord splitterCoord(const Coord& c) const;
    std::vector<ShapePair> matchShapes(GridPtr& bg, GridPtr& ag);
    Point calculateOffset(GridPtr& bg, GridPtr& ag,
        const std::vector<ShapePair>& shapes);
    void addArgs();
    void load();
    void parse(const pdal::StringList& s);
    bool process(Coord c, Point& offset);
    void processGrid();
    void throwError(const std::string& s);
    void writeSimple();
    void writeTiff(const std::string& filename);
    void writeSvg(const std::string& filename, const pdal::BOX2D& extent);
    void read(const std::string& filename);
    bool shouldDump() const;

    pdal::ProgramArgs m_args;
    std::string m_beforeFilename;
    std::string m_afterFilename;

    // Base shift
    Point m_shift;

    // Splitter Origin
    Point m_origin;

    double m_xShift;
    double m_yShift;
    std::unique_ptr<Grid> m_grid;
    Field m_field;
    Coord m_coord;
    Coord m_dumpij;
    Point m_dumpxy;
    double m_dumpfrac;

    pdal::PipelineManager m_beforeMgr;
    pdal::PipelineManager m_afterMgr;
    const double m_len = 100.0;
};

} // namespace
