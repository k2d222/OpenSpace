/*****************************************************************************************
 *                                                                                       *
 * OpenSpace                                                                             *
 *                                                                                       *
 * Copyright (c) 2014-2024                                                               *
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

#include "openspace/util/touch.h"
#include <glm/common.hpp>
#include <modules/touch/include/touchinteraction.h>

#include <modules/touch/include/directinputsolver.h>
#include <modules/touch/touchmodule.h>
#include <openspace/camera/camera.h>
#include <openspace/engine/globals.h>
#include <openspace/engine/moduleengine.h>
#include <openspace/engine/windowdelegate.h>
#include <openspace/navigation/navigationhandler.h>
#include <openspace/navigation/orbitalnavigator.h>
#include <openspace/query/query.h>
#include <openspace/rendering/renderengine.h>
#include <openspace/scene/scene.h>
#include <openspace/scene/scenegraphnode.h>
#include <openspace/util/keys.h>
#include <openspace/util/time.h>
#include <openspace/util/updatestructures.h>
#include <ghoul/format.h>
#include <ghoul/logging/logmanager.h>
#include <ghoul/misc/invariants.h>
#include <glm/gtx/quaternion.hpp>
#include <cmath>
#include <functional>
#include <numeric>

#ifdef WIN32
#pragma warning (push)
#pragma warning (disable : 4310) // cast truncates constant value
#endif // WIN32

#include <glm/ext.hpp>

#ifdef WIN32
#pragma warning (pop)
#endif // WIN32

namespace {
    constexpr std::string_view _loggerCat = "TouchInteraction";

    constexpr openspace::properties::Property::PropertyInfo UnitTestInfo = {
        "UnitTest",
        "Take a unit test saving the LM data into file",
        "LM - least-squares minimization using Levenberg-Marquardt algorithm."
        "Used to find a new camera state from touch points when doing direct "
        "manipulation.",
        openspace::properties::Property::Visibility::Developer
    };

    constexpr openspace::properties::Property::PropertyInfo DisableZoomInfo = {
        "DisableZoom",
        "Disable zoom navigation",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo DisableRollInfo = {
        "DisableRoll",
        "Disable roll navigation",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo SetDefaultInfo = {
        "SetDefault",
        "Reset all properties to default",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo MaxTapTimeInfo = {
        "MaxTapTime",
        "Max tap delay (in ms) for double tap",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo DecelatesPerSecondInfo = {
        "DeceleratesPerSecond",
        "Number of times velocity is decelerated per second",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo TouchScreenSizeInfo = {
        "TouchScreenSize",
        "Touch Screen size in inches",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo TapZoomFactorInfo = {
        "TapZoomFactor",
        "Scaling distance travelled on tap",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo PinchZoomFactorInfo = {
        "PinchZoomFactor",
        "Scaling distance travelled on pinch",
        "This value is used to reduce the amount of pinching needed. A linear kind of "
        "sensitivity that will alter the pinch-zoom speed.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo RollThresholdInfo = {
        "RollThreshold",
        "Threshold for min angle for roll interpret",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo ZoomSensitivityExpInfo = {
        "ZoomSensitivityExp",
        "Sensitivity of exponential zooming in relation to distance from focus node",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo ZoomSensitivityPropInfo = {
        "ZoomSensitivityProp",
        "Sensitivity of zooming proportional to distance from focus node",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo
        ZoomSensitivityDistanceThresholdInfo =
    {
        "ZoomSensitivityDistanceThreshold",
        "Threshold of distance to target node for whether or not to use exponential "
        "zooming",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo
        ZoomInBoundarySphereMultiplierInfo =
    {
        "ZoomInBoundarySphereMultiplier",
        "Multiplies a node's boundary sphere by this in order to limit zoom in & prevent "
        "surface collision.",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo
        ZoomOutBoundarySphereMultiplierInfo =
    {
        "ZoomOutBoundarySphereMultiplier",
        "Multiplies a node's boundary sphere by this in order to limit zoom out",
        "" // @TODO Missing documentation
    };

    constexpr openspace::properties::Property::PropertyInfo ConstantTimeDecaySecsInfo = {
        "ConstantTimeDecaySecs",
        "Time duration that a pitch/roll/zoom/pan should take to decay to zero (seconds)",
        "",  // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo InputSensitivityInfo = {
        "InputSensitivity",
        "Threshold for interpreting input as still",
        "",  // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo StationaryCentroidInfo = {
        "CentroidStationary",
        "Threshold for stationary centroid",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo PanModeInfo = {
        "PanMode",
        "Allow panning gesture",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo PanDeltaDistanceInfo = {
        "PanDeltaDistance",
        "Delta distance between fingers allowed for interpreting pan interaction",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo FrictionInfo = {
        "Friction",
        "Friction for different interactions (orbit, zoom, roll, pan)",
        "", // @TODO Missing documentation
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo ZoomOutLimitInfo = {
        "ZoomOutLimit",
        "Zoom Out Limit",
        "The maximum distance you are allowed to navigate away from the anchor. "
        "This should always be larger than the zoom in value if you want to be able "
        "to zoom. Defaults to maximum allowed double.",
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo ZoomInLimitInfo = {
        "ZoomInLimit",
        "Zoom In Limit",
        "The minimum distance from the anchor that you are allowed to navigate to. "
        "Its purpose is to limit zooming in on a node. If this value is not set it "
        "defaults to the surface of the current anchor.",
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo
        EnableDirectManipulationInfo =
    {
        "EnableDirectManipulation",
        "Enable direct manipulation",
        "Decides whether the direct manipulation mode should be enabled or not.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo
        DirectManipulationThresholdInfo =
    {
        "DirectManipulationThreshold",
        "Direct manipulation threshold",
        "This threshold affects the distance from the interaction sphere at which the "
        "direct manipulation interaction mode starts being active. The value is given "
        "as a factor times the interaction sphere.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    // Compute coefficient of decay based on current frametime; if frametime has been
    // longer than usual then multiple decay steps may be applied to keep the decay
    // relative to user time
    double computeDecayCoeffFromFrametime(double coeff, int times) {
        if (coeff > 0.00001) {
            return std::pow(coeff, times);
        }
        else {
            return 0.0;
        }
    }
} // namespace

namespace openspace {

TouchInteraction::TouchInteraction()
    : properties::PropertyOwner({ "TouchInteraction", "Touch Interaction" })
    , _unitTest(UnitTestInfo, false)
    , _disableZoom(DisableZoomInfo, false)
    , _disableRoll(DisableRollInfo, false)
    , _reset(SetDefaultInfo)
    , _maxTapTime(MaxTapTimeInfo, 300, 10, 1000)
    , _deceleratesPerSecond(DecelatesPerSecondInfo, 240, 60, 300)
    , _touchScreenSize(TouchScreenSizeInfo, 55.f, 5.5f, 150.f)
    , _tapZoomFactor(TapZoomFactorInfo, 0.2f, 0.f, 0.5f, 0.01f)
    , _pinchZoomFactor(PinchZoomFactorInfo, 0.01f, 0.f, 0.2f)
    , _rollAngleThreshold(RollThresholdInfo, 0.025f, 0.f, 0.05f, 0.001f)
    , _zoomSensitivityExponential(ZoomSensitivityExpInfo, 1.03f, 1.f, 1.1f)
    , _zoomSensitivityProportionalDist(ZoomSensitivityPropInfo, 11.f, 5.f, 50.f)
    , _zoomSensitivityDistanceThreshold(
        ZoomSensitivityDistanceThresholdInfo,
        0.05f,
        0.01f,
        0.25f
    )
    , _zoomInBoundarySphereMultiplier(
        ZoomInBoundarySphereMultiplierInfo,
        1.001f,
        0.01f,
        4e+27f
    )
    , _zoomOutBoundarySphereMultiplier(
        ZoomOutBoundarySphereMultiplierInfo,
        4e+27f,
        1.f,
        4e+27f
    )
    , _zoomInLimit(ZoomInLimitInfo, -1.0, 0.0, 4e+27)
    , _zoomOutLimit(
        ZoomOutLimitInfo,
        4e+27,
        1000.0,
        4e+27
    )
    , _inputStillThreshold(InputSensitivityInfo, 0.0005f, 0.f, 0.001f, 0.0001f)
    // Used to void wrongly interpreted roll interactions
    , _centroidStillThreshold(StationaryCentroidInfo, 0.0018f, 0.f, 0.01f, 0.0001f)
    , _panEnabled(PanModeInfo, false)
    , _interpretPan(PanDeltaDistanceInfo, 0.015f, 0.f, 0.1f)
    , _friction(
        FrictionInfo,
        glm::vec4(0.025f, 0.025f, 0.02f, 0.001f),
        glm::vec4(0.f),
        glm::vec4(0.2f)
    )
    , _constTimeDecay_secs(ConstantTimeDecaySecsInfo, 1.75f, 0.1f, 4.f)
    , _enableDirectManipulation(EnableDirectManipulationInfo, true)
    , _directTouchDistanceThreshold(DirectManipulationThresholdInfo, 5.f, 0.f, 10.f)
    , _pinchInputs({ TouchInput(0, 0, 0.f, 0.f, 0.0), TouchInput(0, 0, 0.f, 0.f, 0.0) })
{
    addProperty(_disableZoom);
    addProperty(_disableRoll);
    addProperty(_unitTest);
    addProperty(_reset);
    addProperty(_maxTapTime);
    addProperty(_deceleratesPerSecond);
    addProperty(_touchScreenSize);
    addProperty(_tapZoomFactor);
    addProperty(_pinchZoomFactor);
    addProperty(_rollAngleThreshold);
    addProperty(_zoomSensitivityExponential);
    addProperty(_zoomSensitivityProportionalDist);
    addProperty(_zoomSensitivityDistanceThreshold);
    addProperty(_zoomInBoundarySphereMultiplier);
    addProperty(_zoomOutBoundarySphereMultiplier);
    addProperty(_zoomInLimit);
    addProperty(_zoomOutLimit);
    addProperty(_constTimeDecay_secs);
    addProperty(_inputStillThreshold);
    addProperty(_centroidStillThreshold);
    addProperty(_panEnabled);
    addProperty(_interpretPan);
    addProperty(_friction);

    addProperty(_enableDirectManipulation);
    addProperty(_directTouchDistanceThreshold);

#ifdef TOUCH_DEBUG_PROPERTIES
    addPropertySubOwner(_debugProperties);
#endif

    _zoomInBoundarySphereMultiplier.setExponent(20.f);
    _zoomOutBoundarySphereMultiplier.setExponent(20.f);
    _zoomInLimit.setExponent(20.f);
    _zoomOutLimit.setExponent(20.f);
    _time = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now().time_since_epoch()
    );

    _reset.onChange([this]() { resetPropertiesToDefault(); });
}

void TouchInteraction::update(const std::vector<TouchInputHolder>& inputs)
{
    size_t numFingers = inputs.size();

#ifdef TOUCH_DEBUG_PROPERTIES
    _debugProperties.nFingers = numFingers;
#endif

    // just ended an interaction by relasing the last finger
    if (numFingers == 0 && _isDirectControlActive) {
        _anchor = nullptr;
        _isDirectControlActive = false;
        endDirectControl();
    }

    // just started an interaction by pressing a finger
    else if(numFingers != 0 && !_isDirectControlActive) {
        _anchor = global::navigationHandler->orbitalNavigator().anchorNode();
        _isDirectControlActive = startDirectControl(inputs);
    }

    // interaction in progress
    else if (_isDirectControlActive) {
        _isDirectControlActive = directControl(inputs);

        if (!_isDirectControlActive) {
            _anchor = nullptr;
            endDirectControl();
        }
    }
}

void TouchInteraction::reset() {
    if (_isDirectControlActive) {
        _anchor = nullptr;
        _isDirectControlActive = false;
        interaction::OrbitalNavigator& orbNav = global::navigationHandler->orbitalNavigator();
        orbNav.resetVelocities();
    }
}

glm::dvec3 TouchInteraction::unprojectTouchOnSphere(const glm::vec2& input) const {
    const glm::dvec3& camPos = _startPose.position;
    const glm::dquat& camQuat = _startPose.rotation;

    // convert touch input in range [0, 1] to NDC in range [-1, 1]
    glm::vec4 inputNdc(input.x * 2.f - 1.f, -(input.y * 2.f - 1.f), -1.f, 1.f);
    glm::mat4 ndcToCam = glm::inverse(_camera->projectionMatrix());
    glm::dvec3 inputWorldSpace = camQuat * glm::dvec3(ndcToCam * inputNdc); // it's not really in world space, the origin is the camera position.
    glm::dvec3 rayDir = glm::normalize(inputWorldSpace);

    // Compute positions on anchor node, by checking if touch input
    // intersect interaction sphere
    double intersectionDist = 0.0;
    const bool intersected = glm::intersectRaySphere(
        glm::dvec3(0.0), // it's more stable to do the computation at camera origin
        rayDir,
        _anchor->worldPosition() - camPos,
        _anchor->interactionSphere() * _anchor->interactionSphere(),
        intersectionDist
    );

    if (intersected) {
        return camPos + rayDir * intersectionDist;
    }
    else {
        return glm::dvec3(0.0);
    }
}

bool TouchInteraction::directControl(const std::vector<TouchInputHolder>& inputs) {
#ifdef TOUCH_DEBUG_PROPERTIES
    LINFO("DirectControl");
#endif

    // check if direct control is possible
    if (inputs.size() == 0 || !_anchor || !_camera) {
        return false;
    }

    const TouchInputHolder& input = inputs[0];

    interaction::OrbitalNavigator& orbNav = global::navigationHandler->orbitalNavigator();

    const glm::vec2 lastInput(input.latestInput().x, input.latestInput().y);
    const glm::dvec3 startWorldSpace = unprojectTouchOnSphere(_firstInput);
    const glm::dvec3 endWorldSpace = unprojectTouchOnSphere(lastInput);
    const glm::dvec3 spherePosition = _anchor->worldPosition();

    // check if raycast failed (outside anchor sphere)
    if (startWorldSpace == glm::dvec3(0.0) || endWorldSpace == glm::dvec3(0.0)) {
        return false;
    }
    
    glm::vec3 startVec = glm::normalize(startWorldSpace - spherePosition);
    glm::vec3 endVec = glm::normalize(endWorldSpace - spherePosition);

    glm::vec3 rotAxis = glm::normalize(glm::cross(startVec, endVec));
    float rotAngle = glm::acos(glm::dot(startVec, endVec));
    glm::dquat rotQuat = glm::angleAxis(-rotAngle, rotAxis); // invert the angle as we want to rotate the camera around the sphere, not the opposite

    // rotAxis is NaN when startVec and endVec are parallel (no rotation needed)
    if (glm::any(glm::isnan(rotAxis))) {
        return true;
    }

    _lastPoses[0] = _lastPoses[1];
    _lastPoses[1].position = rotQuat * (_startPose.position - spherePosition) + spherePosition;
    _lastPoses[1].rotation = rotQuat * _startPose.rotation;
    _camera->setPose(_lastPoses[1]);

    // Mark that a camera interaction happeneddD
    orbNav.updateOnCameraInteraction();

    return true;
}

bool TouchInteraction::startDirectControl(const std::vector<TouchInputHolder>& inputs) {
    if (inputs.size() == 0 || !_camera || !_anchor) {
        return false;
    }

    // TODO: implement Camera::pose() (also used in OrbitalNavigator)
    _startPose = {
        .position = _camera->positionVec3(),
        .rotation = _camera->rotationQuaternion(),
    };

    const TouchInputHolder& input = inputs[0];
    _firstInput = glm::vec2(input.latestInput().x, input.latestInput().y);
    glm::dvec3 proj = unprojectTouchOnSphere(_firstInput);

    // check if touch is inside interaction sphere
    bool isInside = proj != glm::dvec3(0.0);

    if (isInside) {
        interaction::OrbitalNavigator& orbNav = global::navigationHandler->orbitalNavigator();
        orbNav.resetVelocities();
    }

    return isInside;
}

void TouchInteraction::endDirectControl() {
    glm::dquat rotDiff = glm::inverse(_lastPoses[0].rotation) * _lastPoses[1].rotation;
    glm::vec3 axisAngle = glm::axis(rotDiff) * glm::angle(rotDiff);
    glm::vec3 angularVel = glm::inverse(_camera->viewMatrix()) * glm::vec4(axisAngle, 0.f);

    interaction::OrbitalNavigator& orbNav = global::navigationHandler->orbitalNavigator();
    // XXX: idk why y and x axes are swapped here.
    // TODO: compute dt to have correct velocity
    // TODO: ensure _lastPoses were computed
    orbNav.touchStates().setGlobalRotationVelocity(glm::dvec2(angularVel.y, angularVel.x));
}

double TouchInteraction::computeTapZoomDistance(double zoomGain) {
    const SceneGraphNode* anchor =
        global::navigationHandler->orbitalNavigator().anchorNode();

    if (!anchor) {
        return 0.0;
    }

    double dist = glm::distance(_camera->positionVec3(), anchor->worldPosition());
    dist -= anchor->interactionSphere();

    double newVelocity = dist * _tapZoomFactor;
    newVelocity *= std::max(_touchScreenSize.value() * 0.1, 1.0);
    newVelocity *= _zoomSensitivityProportionalDist * zoomGain;

    return newVelocity;
}

// Reset all property values to default
void TouchInteraction::resetPropertiesToDefault() {
    _unitTest.set(false);
    _disableZoom.set(false);
    _disableRoll.set(false);
    _maxTapTime.set(300);
    _deceleratesPerSecond.set(240);
    _touchScreenSize.set(55.f);
    _tapZoomFactor.set(0.2f);
    _pinchZoomFactor.set(0.01f);
    _rollAngleThreshold.set(0.025f);
    _zoomSensitivityExponential.set(1.025f);
    _inputStillThreshold.set(0.0005f);
    _centroidStillThreshold.set(0.0018f);
    _interpretPan.set(0.015f);
    _friction.set(glm::vec4(0.025f, 0.025f, 0.02f, 0.02f));
}

void TouchInteraction::tap() {
    // _tap = true;
    //TODO
}

void TouchInteraction::setCamera(Camera* camera) {
    _camera = camera;
}

void FrameTimeAverage::updateWithNewFrame(double sample) {
    if (sample > 0.0005) {
        _samples[_index++] = sample;
        if (_index >= TotalSamples) {
            _index = 0;
        }
        if (_nSamples < TotalSamples) {
            _nSamples++;
        }
    }
}

double FrameTimeAverage::averageFrameTime() const {
    if (_nSamples == 0) {
        // Just guess at 60fps if no data is available yet
        return 1.0 / 60.0;
    }
    else {
        return std::accumulate(_samples, _samples + _nSamples, 0.0) / (double)(_nSamples);
    }
}

#ifdef TOUCH_DEBUG_PROPERTIES
TouchInteraction::DebugProperties::DebugProperties()
    : properties::PropertyOwner({ "TouchDebugProperties", "Touch Debug Properties"})
    , interactionMode(
        { "interactionMode", "Current interaction mode", "" },
        "Unknown"
    )
    , nFingers(
        {"nFingers", "Number of fingers", ""},
        0, 0, 20
    )
    , interpretedInteraction(
        { "interpretedInteraction", "Interpreted interaction", "" },
        "Unknown"
    )
    , normalizedCentroidDistance(
        { "normalizedCentroidDistance", "Normalized Centroid Distance", "" },
        0.f, 0.f, 0.01f
    )
    , minDiff(
        { "minDiff", "Movement of slowest moving finger", "" },
        0.f, 0.f, 100.f
    )
    , rollOn(
        { "rollOn", "Roll On", "" },
        0.f, 0.f, 100.f
    )
{
    addProperty(interactionMode);
    addProperty(nFingers);
    addProperty(interpretedInteraction);
    addProperty(normalizedCentroidDistance);
    addProperty(minDiff);
    addProperty(rollOn);
}
#endif

} // openspace namespace
