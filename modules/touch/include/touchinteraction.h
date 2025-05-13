/*****************************************************************************************
 *                                                                                       *
 * OpenSpace                                                                             *
 *                                                                                       *
 * Copyright (c) 2014-2025                                                               *
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

#ifndef __OPENSPACE_MODULE_TOUCH___TOUCH_INTERACTION___H__
#define __OPENSPACE_MODULE_TOUCH___TOUCH_INTERACTION___H__

#include <modules/touch/include/ringbuffer.h>
#include <openspace/camera/camerapose.h>
#include <openspace/properties/propertyowner.h>
#include <openspace/properties/misc/triggerproperty.h>
#include <openspace/properties/scalar/boolproperty.h>
#include <openspace/properties/scalar/floatproperty.h>
#include <openspace/properties/scalar/doubleproperty.h>
#include <openspace/properties/scalar/intproperty.h>
#include <openspace/properties/vector/ivec2property.h>
#include <openspace/properties/vector/vec4property.h>
#include <openspace/util/touch.h>

// #define TOUCH_DEBUG_PROPERTIES

namespace openspace {

class Camera;
class SceneGraphNode;

class TouchInteraction : public properties::PropertyOwner {
public:
    TouchInteraction();

    void update(const std::vector<TouchInputHolder>& list);
    void reset();

    bool hasNonZeroVelocities() const;

    /**
     * Sets _tap to true, called if tap occured current frame (called from touchmodule).
     */
    void tap();

    void setCamera(Camera* camera);

private:
    /**
     * Function that calculates the new camera state such that it minimizes the L2 error
     * in screenspace between contact points and surface coordinates projected to clip
     * space using LMA.
     */
    bool directControl(const std::vector<TouchInputHolder>& input);
    bool startDirectControl(const std::vector<TouchInputHolder>& input);
    void endDirectControl();

    /**
     * Compute velocity based on double-tap for zooming.
     */
    double computeTapZoomDistance(double zoomGain);

    /**
     * Resets all properties that can be changed in the GUI to default.
     */
    void resetPropertiesToDefault();

    /**
     * Compute the touchInput position on the target node's interaction sphere.
     * Returns a direction from the node's position.
     */
    glm::dvec3 unprojectTouchOnSphere(const TouchInput& input) const;

    /**
     * Compute the unprojected barycenter of a set touchInputs on the target node's
     * interaction sphere.
     */
    glm::dvec3 unprojectTouchesOnSphere(const std::vector<TouchInput>& inputs) const;


    // transforms keep track of the last N transformations to apply velocity when the
    // interaction is ended.
    struct Transforms {
        glm::dquat rotation;
        double scaling;
        double timestamp;
    };
    RingBuffer<Transforms> _lastTransforms;

    std::vector<TouchInput> _startInputs;
    CameraPose _startPose = {};

    Camera* _camera = nullptr;
    const SceneGraphNode* _anchor = nullptr;
    bool _isDirectControlActive = false;

    // Property variables
    properties::BoolProperty _enableOrbit;
    properties::BoolProperty _enablePinchZoom;
    properties::BoolProperty _enableTapZoom;
    properties::BoolProperty _enableRoll;
    properties::BoolProperty  _enablePan;
    properties::BoolProperty _enableDirectManipulation;
    properties::BoolProperty _enableVelocity;

    properties::IntProperty _doubleTapThreshold;
    properties::FloatProperty _rollAngleThreshold;
    properties::FloatProperty _panDistanceThreshold;
    properties::FloatProperty _stationaryThreshold;
    properties::FloatProperty _centroidStillThreshold;
    properties::FloatProperty _directManipulationThreshold;

    properties::FloatProperty _tapZoomFactor;
    properties::FloatProperty _zoomInBoundarySphereMultiplier;
    properties::FloatProperty _zoomOutBoundarySphereMultiplier;
    properties::DoubleProperty _zoomInLimit;
    properties::DoubleProperty _zoomOutLimit;
    properties::IntProperty _velocitySampleCount;
    properties::BoolProperty _useSphericalDisplay;
    properties::TriggerProperty _reset;

#ifdef TOUCH_DEBUG_PROPERTIES
    struct DebugProperties : PropertyOwner {
        DebugProperties();
        properties::StringProperty interactionMode;
        properties::IntProperty nFingers;
        properties::StringProperty interpretedInteraction;
        properties::FloatProperty normalizedCentroidDistance;
        properties::FloatProperty minDiff;
        properties::FloatProperty rollOn;
    } _debugProperties;

    int pinchConsecCt = 0;
    double pinchConsecZoomFactor = 0;
    int stepVelUpdate = 0;
#endif
};

} // openspace namespace

#endif // __OPENSPACE_MODULE_TOUCH___TOUCH_INTERACTION___H__
