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

#include <functional>
#include <cassert>

namespace AtlasProcessor
{

struct GridIndex
{
public:
    GridIndex(int32_t x, int32_t y) : m_key(key(x, y))
    {
        assert(x == this->x());
        assert(y == this->y());
    }

    int32_t x() const
    { return (uint32_t)m_key; }
    int32_t y() const
    { return (uint32_t)(m_key >> 32); }
    uint64_t key() const
    { return m_key; }
    bool operator==(const GridIndex& other) const
    { return m_key == other.m_key; }

private:
    uint64_t key(int32_t x, int32_t y)
    {
        uint32_t ux = (uint32_t)x;
        uint32_t uy = (uint32_t)y;
        return (ux | ((uint64_t)uy << 32));
    }

    uint64_t m_key;
};

} // namespace AtlasProcessor

namespace std
{
    template<>
    struct hash<AtlasProcessor::GridIndex>
    {
        using argument_type = AtlasProcessor::GridIndex;
        using result_type = std::size_t;

        result_type operator()(const AtlasProcessor::GridIndex& val) const
        { return (result_type)val.key(); }
    };
} // namespace std

