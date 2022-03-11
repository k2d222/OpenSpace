/*****************************************************************************************
 *                                                                                       *
 * OpenSpace                                                                             *
 *                                                                                       *
 * Copyright (c) 2014-2020                                                               *
 *                                                                                       *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this  *
 * software and associated documentation files (the "Software"), to deal in the Software *
 * without restriction, including without limitation the rights to use, copy, modify,    *
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to    *
 * permit persons to whom the Software is furnished to do so, subject to the following   *
 * conditions:                                                                           *
 *                                                                                       *
 * The above copyright notice and this permission notice shall be included in all copies *
 * or substantial portions of the Software.                                              *
 *                                                                                       *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,   *
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A         *
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT    *
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF  *
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE  *
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                         *
 ****************************************************************************************/

#include "fragment.glsl"

in float gs_depthClipSpace;
in vec4 gs_color;
in float gs_component;
in vec2 texCoord;

uniform float opacity;

Fragment getFragment() {
    float radius = length(texCoord);

    if (radius > 1.0)
        discard;

    float maxRadius = 1.0;
    float width =  0.8 / gs_component;
    float minRadius = 1.0 - width;

    float coord = (radius - minRadius) / (maxRadius - minRadius);
    if (coord < 0.0 || coord > 1.0) {
        discard;
    }

    vec4 color = gs_color;
    color.a *= opacity;

    Fragment frag;
    frag.color = color;
    frag.depth = gs_depthClipSpace;
//    frag.gPosition = vec4(-1e32, -1e32, -1e32, 1.0); // From DU billboards shader code
//
    // There is no normal here
    frag.gNormal = vec4(0.0, 0.0, -1.0, 1.0);

    return frag;
}
