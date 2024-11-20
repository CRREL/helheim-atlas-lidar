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
#include "Atlas.hpp"

int main(int argc, const char *argv[])
{
    pdal::StringList slist(argv + 1, argv + argc);

    AtlasProcessor::Atlas atlas;
    atlas.run(slist);
}
