/*****************************************************************************************
 *                                                                                       *
 * OpenSpace                                                                             *
 *                                                                                       *
 * Copyright (c) 2014-2022                                                               *
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

#version __CONTEXT__

#include "PowerScaling/powerScalingMath.hglsl"

layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

flat in float vs_component[];
in float vs_depthClipSpace[];
in vec4 vs_color[];

flat out float gs_component;
out float gs_depthClipSpace;
out vec4 gs_positionViewSpace;
out vec4 gs_color;
out vec2 texCoord;

uniform dmat4 modelMatrix;
uniform dmat4 cameraViewProjectionMatrix;
uniform float size;

const vec2 corners[4] = vec2[4](
    vec2(-1.0, -1.0),
    vec2(1.0, -1.0),
    vec2(-1.0, 1.0),
    vec2(1.0, 1.0)
);

void main() {
    vec4 pos = gl_in[0].gl_Position; // model space

    gs_component = vs_component[0];
    gs_depthClipSpace = vs_depthClipSpace[0];
    gs_color = vs_color[0];

    float scale = exp(size) * vs_component[0];

    // Apply projection matrix as last step
    mat4 MVPMatrix = mat4(cameraViewProjectionMatrix * modelMatrix);

    // The added vector needs to be in clip space

    // Lower left
    texCoord = corners[0];
    vec4 offset = vec4(corners[0], 0.0, 0.0);
    gl_Position = z_normalization(MVPMatrix * (pos + scale * offset));
    EmitVertex();

    // Lower right
    texCoord = corners[1];
    offset = vec4(corners[1], 0.0, 0.0);
    gl_Position = z_normalization(MVPMatrix * (pos + scale * offset));
    EmitVertex();

    // Upper left
    texCoord = corners[2];
    offset = vec4(corners[2], 0.0, 0.0);
    gl_Position = z_normalization(MVPMatrix * (pos + scale * offset));
    EmitVertex();

    // Upper right
    texCoord = corners[3];
    offset = vec4(corners[3], 0.0, 0.0);
    gl_Position = z_normalization(MVPMatrix * (pos + scale * offset));
    EmitVertex();

    //gl_PointSize = scale;

    EndPrimitive();
}
