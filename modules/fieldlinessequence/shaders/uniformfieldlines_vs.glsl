/*****************************************************************************************
 *                                                                                       *
 * OpenSpace                                                                             *
 *                                                                                       *
 * Copyright (c) 2014-2021                                                               *
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
layout(location = 0) in vec3 in_position;

out float vs_depth;
out vec4 vs_positionViewSpace;
out vec4 vs_color;


uniform mat4 modelViewTransform;
uniform mat4 projectionTransform;
// General Uniforms that's always needed
uniform vec4      lineColor;

// Uniforms needed for Particle Flow
uniform vec4      flowColor;
uniform float     particleSize;
uniform float     particleSpeed;
uniform float     particleSpacing;
uniform double    time;
uniform bool      usingParticles;



bool isPartOfParticle(double time, int vertexId, float particleSize,
                      float particleSpeed, float particleSpacing)
{
    int modulusResult = int(particleSpeed * time + vertexId) % int(particleSpacing);
    return modulusResult > 0 && modulusResult <= particleSize;
}

void main() {
    bool hasColor = true;

    if (hasColor) {
        bool isParticle = usingParticles && isPartOfParticle(time, gl_VertexID,
                                                                    particleSize,
                                                                    particleSpeed,
                                                                    particleSpacing);
        if (isParticle) {
            vs_color = flowColor;
        }
        else {
            vs_color = lineColor;
        }

    }
    else {
        vs_color = vec4(0);
    }

    vs_positionViewSpace = vec4(modelViewTransform * dvec4(in_position, 1));
    vec4 positionScreenSpace = projectionTransform * vs_positionViewSpace;
    vs_depth = positionScreenSpace.w;
    gl_Position  = positionScreenSpace;

    // Set z to 0 to disable near and far plane, unique handling for perspective in space
    gl_Position.z = 0.f;
}
