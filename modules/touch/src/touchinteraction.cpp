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

#include <modules/touch/include/touchinteraction.h>

#include <modules/touch/touchmodule.h>
#include <openspace/camera/camera.h>
#include <openspace/camera/camerapose.h>
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
#include <openspace/util/touch.h>
#include <openspace/util/updatestructures.h>
#include <ghoul/format.h>
#include <ghoul/logging/logmanager.h>
#include <ghoul/misc/invariants.h>
#include <glm/geometric.hpp>
#include <glm/common.hpp>
#include <glm/gtx/quaternion.hpp>

#include <cmath>
#include <functional>

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

    constexpr openspace::properties::Property::PropertyInfo EnableOrbitInfo = {
        "EnableOrbit",
        "Enable orbital navigation",
        "Enable orbiting around the anchor by sliding one finger.",
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo EnablePinchZoomInfo = {
        "EnablePinchZoom",
        "Enable pinch zoom navigation",
        "Enable the pinch gesture to zoom in on the anchor.",
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo EnableTapZoomInfo = {
        "EnableTapZoom",
        "Enable double tap zoom navigation",
        "Enable double tapping to zoom in on the anchor.",
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo EnableRollInfo = {
        "EnableRoll",
        "Enable roll navigation",
        "Enable rolling around an axis when twisting two fingers.",
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo EnablePanInfo = {
        "EnablePan",
        "Enable panning navigation",
        "Enable Panning by sliding two fingers.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo
        EnableDirectManipulationInfo =
    {
        "EnableDirectManipulation",
        "Enable direct manipulation",
        "When direct manipulation is enabled, the targed points on the anchor will "
        "'stick' to the fingers. When disabled, finger gestures simply add velocity to "
        "the camera.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo
        EnableVelocityInfo =
    {
        "EnableVelocity",
        "Enable velocity",
        "Compute and apply a velocity to the OrbitalNavigator at the end of the "
        "interaction.",
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo DoubleTapThresholdInfo = {
        "DoubleTapThreshold",
        "Maximum double tap delay (ms)",
        "Maximum delay between taps to interpret a double tap.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo RollAngleThresholdInfo = {
        "RollAngleThreshold",
        "Roll angle threshold",
        "Minimum angle to interpret a gesture as a roll.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo PanDistanceThresholdInfo = {
        "PanDistanceThreshold",
        "Pan distance threshold",
        "Delta distance between fingers allowed for interpreting pan interaction",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo StationaryThresholdInfo = {
        "StationaryThreshold",
        "Stationary threshold",
        "Maximum distance to interpret input as stationary.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo CentroidStillThresholdInfo = {
        "CentroidStillThreshold",
        "Threshold for stationary centroid",
        "Maximum distance to interpret the input centroid as stationary",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo
        DirectManipulationThresholdInfo =
    {
        "DirectManipulationThreshold",
        "Direct manipulation distance threshold",
        "This threshold affects the distance from the interaction sphere at which the "
        "direct manipulation interaction mode starts being active. The value is given "
        "as a factor times the interaction sphere.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo TapZoomFactorInfo = {
        "TapZoomFactor",
        "Tap zoom scaling factor",
        "Factor to scale the distance to the surface on double-tap-zoom.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo
        ZoomInBoundarySphereMultiplierInfo =
    {
        "ZoomInBoundarySphereMultiplier",
        "Zoom in multiplier",
        "Multiplies a node's boundary sphere by this in order to limit zoom in & prevent "
        "surface collision.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo
        ZoomOutBoundarySphereMultiplierInfo =
    {
        "ZoomOutBoundarySphereMultiplier",
        "Zoom out multiplier",
        "Multiplies a node's boundary sphere by this in order to limit zoom out",
    };

    constexpr openspace::properties::Property::PropertyInfo ZoomOutLimitInfo = {
        "ZoomOutLimit",
        "Zoom out limit",
        "The maximum distance you are allowed to navigate away from the anchor. "
        "This should always be larger than the zoom in value if you want to be able "
        "to zoom. Defaults to maximum allowed double.",
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo ZoomInLimitInfo = {
        "ZoomInLimit",
        "Zoom in limit",
        "The minimum distance from the anchor that you are allowed to navigate to. "
        "Its purpose is to limit zooming in on a node. Defaults to the surface of the "
        "current anchor.",
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo
        VelocitySampleCountInfo =
    {
        "VelocitySampleCount",
        "Velocity sample count",
        "Number of touch samples to keep in memory, to compute the average velocity at "
        "the end of interaction. The previous steps will be averaged. It should "
        "correspond to roughly 0.1 seconds worth of sampling, check your touch device "
        "sampling rate.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo
        SphericalDisplayInfo =
    {
        "SphericalDisplay",
        "Spherical display",
        "Enable this for spherical touch displays (e.g. the Pufferfish touch)."
        "The camera will be fixed at the anchor position. Use in conjuntion with the "
        "outside-in cubemap fisheye configuration.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    constexpr openspace::properties::Property::PropertyInfo ResetInfo = {
        "Reset",
        "Reset",
        "Reset all properties to default",
        openspace::properties::Property::Visibility::AdvancedUser
    };


    std::vector<openspace::TouchInput> lastInputs(const std::vector<openspace::TouchInputHolder>& inputs) {
        std::vector<openspace::TouchInput> res;
        res.reserve(inputs.size());

        for (const openspace::TouchInputHolder& input : inputs) {
            res.push_back(input.latestInput());
        }

        return res;
    }
} // namespace

namespace openspace {

TouchInteraction::TouchInteraction()
    : properties::PropertyOwner({ "TouchInteraction", "Touch Interaction" })
    , _lastTransforms(10)
    , _enableOrbit(EnableOrbitInfo, true)
    , _enablePinchZoom(EnablePinchZoomInfo, true)
    , _enableTapZoom(EnableTapZoomInfo, true)
    , _enableRoll(EnableRollInfo, true)
    , _enablePan(EnablePanInfo, false)
    , _enableDirectManipulation(EnableDirectManipulationInfo, true)
    , _enableVelocity(EnableVelocityInfo, true)
    , _doubleTapThreshold(DoubleTapThresholdInfo, 300, 10, 1000)
    , _rollAngleThreshold(RollAngleThresholdInfo, 0.025f, 0.f, 0.05f, 0.001f)
    , _panDistanceThreshold(PanDistanceThresholdInfo, 0.015f, 0.f, 0.1f)
    , _stationaryThreshold(StationaryThresholdInfo, 0.0005f, 0.f, 0.001f, 0.0001f)
    , _centroidStillThreshold(CentroidStillThresholdInfo, 0.0018f, 0.f, 0.01f, 0.0001f)
    , _directManipulationThreshold(DirectManipulationThresholdInfo, 5.f, 0.f, 10.f)
    , _tapZoomFactor(TapZoomFactorInfo, 0.2f, 0.f, 0.5f, 0.01f)
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
    , _zoomInLimit(ZoomInLimitInfo, -1.f, 0.f, 4e+27f)
    , _zoomOutLimit(
        ZoomOutLimitInfo,
        4e+27f,
        1000.f,
        4e+27f
    )
    , _velocitySampleCount(VelocitySampleCountInfo, 10, 2, 100)
    , _useSphericalDisplay(SphericalDisplayInfo, false)
    , _reset(ResetInfo)
{
    addProperty(_enableOrbit);
    addProperty(_enablePinchZoom);
    addProperty(_enableTapZoom);
    addProperty(_enableRoll);
    addProperty(_enablePan);
    addProperty(_enableDirectManipulation);
    addProperty(_enableVelocity);
    addProperty(_doubleTapThreshold);
    addProperty(_rollAngleThreshold);
    addProperty(_panDistanceThreshold);
    addProperty(_stationaryThreshold);
    addProperty(_centroidStillThreshold);
    addProperty(_directManipulationThreshold);
    addProperty(_tapZoomFactor);
    addProperty(_zoomInBoundarySphereMultiplier);
    addProperty(_zoomOutBoundarySphereMultiplier);
    addProperty(_zoomInLimit);
    addProperty(_zoomOutLimit);
    addProperty(_velocitySampleCount);
    addProperty(_useSphericalDisplay);
    addProperty(_reset);

#ifdef TOUCH_DEBUG_PROPERTIES
    addPropertySubOwner(_debugProperties);
#endif

    _zoomInBoundarySphereMultiplier.setExponent(20.f);
    _zoomOutBoundarySphereMultiplier.setExponent(20.f);
    _zoomInLimit.setExponent(20.f);
    _zoomOutLimit.setExponent(20.f);

    _lastTransforms.resize(_velocitySampleCount.value());

    _velocitySampleCount.onChange([this]() {
        _lastTransforms.resize(_velocitySampleCount.value());
    });

    _reset.onChange([this]() {
        resetPropertiesToDefault();
    });
}

void TouchInteraction::update(const std::vector<TouchInputHolder>& inputs)
{
    size_t numFingers = inputs.size();

#ifdef TOUCH_DEBUG_PROPERTIES
    _debugProperties.nFingers = numFingers;
#endif

    if (_useSphericalDisplay) {
        const SceneGraphNode* anchor = global::navigationHandler->orbitalNavigator().anchorNode();
        if (anchor) {
            _camera->setPositionVec3(anchor->worldPosition() + glm::dvec3(0.0, 10.0, 0.0));
        }
    }

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

    // added / removed a finger
    else if (_isDirectControlActive && numFingers != _startInputs.size()) {
        _isDirectControlActive = startDirectControl(inputs);
        if (_isDirectControlActive) {
            _isDirectControlActive = directControl(inputs);

            if (!_isDirectControlActive) {
                _anchor = nullptr;
                endDirectControl();
            }
        }
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

glm::dvec3 TouchInteraction::unprojectTouchOnSphere(const TouchInput& input) const {
    if (_useSphericalDisplay) {
        double yaw = (1.5 - input.x) * glm::two_pi<double>(); // range 0 to 2pi
        double pitch = (input.y - 0.5) * glm::pi<double>(); // range -pi/2 to pi/2
        glm::dvec3 dir = glm::normalize(glm::dvec3(cos(yaw) * cos(pitch), sin(pitch), sin(yaw) * cos(pitch)));
        // NOTE: on spherical displays, that direction is relative to the camera space.
        return dir;
    }

    const glm::dvec3& camPos = _startPose.position;
    const glm::dquat& camRot = _startPose.rotation;

    // convert touch input in range [0, 1] to NDC in range [-1, 1], then unproject
    glm::vec4 inputNdc(input.x * 2.f - 1.f, -(input.y * 2.f - 1.f), -1.f, 1.f);
    glm::mat4 ndcToCam = glm::inverse(_camera->sgctInternal.projectionMatrix());
    // it's not really in world space, the origin is the camera position
    glm::dvec3 inputWorldSpace = camRot * glm::dvec3(ndcToCam * inputNdc);
    glm::dvec3 rayDir = glm::normalize(inputWorldSpace);

    // Compute positions on anchor node, by checking if touch input intersects the
    // interaction sphere
    double intersectionDist = 0.0;
    const bool intersected = glm::intersectRaySphere(
        glm::dvec3(0.0), // it's more stable to do the computation at the camera origin
        rayDir,
        _anchor->worldPosition() - camPos,
        _anchor->interactionSphere() * _anchor->interactionSphere(),
        intersectionDist
    );

    if (intersected) {
        const glm::dvec3 worldPos = camPos + rayDir * intersectionDist;
        const glm::dvec3 worldDir = glm::normalize(worldPos - _anchor->worldPosition());
        return worldDir;
    }
    else {
        return glm::dvec3(0.0);
    }
}

glm::dvec3 TouchInteraction::unprojectTouchesOnSphere(const std::vector<TouchInput>& inputs) const {
    glm::dvec3 acc;

    for (const TouchInput input : inputs) {
        acc += unprojectTouchOnSphere(input);
    }

    return glm::normalize(acc);
}

bool TouchInteraction::directControl(const std::vector<TouchInputHolder>& inputs) {
#ifdef TOUCH_DEBUG_PROPERTIES
    LINFO("DirectControl");
#endif

    // check if direct control is possible
    if (inputs.size() == 0 || !_anchor || !_camera) {
        return false;
    }

    _startPose = _camera->pose();

    // constants
    const std::vector<TouchInput> endInputs = lastInputs(inputs);
    const glm::dvec3 anchorPos = _anchor->worldPosition();
    const glm::dvec3 camAxis = glm::normalize(anchorPos - _startPose.position);

    const glm::dvec3 startP1Dir = unprojectTouchOnSphere(_startInputs.front());
    const glm::dvec3 endP1Dir = unprojectTouchOnSphere(endInputs.front());
    const glm::dvec3 startBarycenterDir = unprojectTouchesOnSphere(_startInputs);
    glm::dvec3 endBarycenterDir = unprojectTouchesOnSphere(endInputs);
    glm::dvec3 startP2Dir = unprojectTouchOnSphere(_startInputs.back());
    glm::dvec3 endP2Dir = unprojectTouchOnSphere(endInputs.back());

    // check if raycast failed (outside anchor sphere)
    if (startP1Dir == glm::dvec3(0.0) || endP1Dir == glm::dvec3(0.0)) {
        return false;
    }

    // transform parameters
    double scaling = 1.0;
    glm::dvec3 orbitAxis = glm::dvec3(1.0, 0.0, 0.0);
    double orbitAngle = 0.0;
    glm::dvec3 rollAxis = glm::dvec3(1.0, 0.0, 0.0);
    double rollAngle = 0.0;

    // compute scaling (first because orbit and roll depend on it)
    if (_enablePinchZoom && inputs.size() >= 2) {
        // this is a simplified scaling. Proper scaling calculation might be impossible to find
        // analytically. I challenge someone to find it, it would be a great contribution.
        double startAngle = glm::distance(startP1Dir, startP2Dir);
        double endAngle = glm::distance(endP1Dir, endP2Dir);
        scaling = endAngle / startAngle;
 
        const glm::dvec3 surface = anchorPos - camAxis * _anchor->interactionSphere(); // we scale the distance from the camera to the surface of the globle by scaling
        _startPose.position = (_startPose.position - surface) / scaling + surface;
        startP2Dir = unprojectTouchOnSphere(_startInputs.back());
        endP2Dir = unprojectTouchOnSphere(endInputs.back());
        endBarycenterDir = unprojectTouchesOnSphere(endInputs);

        // check if raycast failed (outside anchor sphere)
        if (startP2Dir == glm::dvec3(0.0) || endP2Dir == glm::dvec3(0.0)) {
            return false;
        }
    }

    // compute orbit
    if (_enableOrbit) {
        orbitAxis = glm::normalize(glm::cross(startBarycenterDir, endBarycenterDir));
        orbitAngle = glm::acos(glm::dot(startBarycenterDir, endBarycenterDir));

        // failure can happen if startP1Dir and endP1Dir are colinear (no orbit)
        if (glm::any(glm::isnan(orbitAxis))) {
            orbitAxis = glm::dvec3(1.0, 0.0, 0.0);
            orbitAngle = 0.0;
        }
    }

    // compute camera z-axis rotation (roll)
    if (_enableRoll && inputs.size() >= 2) {
        rollAxis = endBarycenterDir;
        // need to project all vectors on the rotation plane to compute the angle from
        // this point of view.
        rollAngle = -glm::orientedAngle(glm::normalize(startP2Dir - startP1Dir), glm::normalize(endP2Dir - endP1Dir), rollAxis);
    }

    glm::dquat rotation =
        glm::angleAxis(-orbitAngle, orbitAxis) *
        glm::angleAxis(rollAngle, rollAxis);

    // compute transforms
    CameraPose endPose;
    endPose.position = rotation * (_startPose.position - anchorPos) + anchorPos;
    endPose.rotation = rotation * _startPose.rotation;

    Transforms transforms;
    transforms.rotation = glm::inverse(_startPose.rotation) * rotation * _startPose.rotation; // transform global to local rotation
    transforms.scaling = scaling - 1.0;
    transforms.timestamp = endInputs.back().timestamp;
    _startInputs = endInputs;

    if (_useSphericalDisplay) {
        // on spherical displays the rotation is applied locally
        endPose.position = anchorPos + glm::dvec3(0.0, 10.0, 0.0);
        endPose.rotation = _startPose.rotation * rotation; 
        transforms.rotation = rotation;
    }

    _camera->setPose(endPose);
    _lastTransforms.push(transforms);

    // Mark that a camera interaction happened
    interaction::OrbitalNavigator& orbNav = global::navigationHandler->orbitalNavigator();
    orbNav.updateOnCameraInteraction();

    return true;
}

bool TouchInteraction::startDirectControl(const std::vector<TouchInputHolder>& inputs) {
    if (inputs.size() == 0 || !_camera || !_anchor) {
        return false;
    }

    _startInputs = lastInputs(inputs);
    glm::dvec3 proj = unprojectTouchesOnSphere(_startInputs);

    // check if touch is inside interaction sphere
    bool isInside = proj != glm::dvec3(0.0);

    if (isInside) {
        interaction::OrbitalNavigator& orbNav = global::navigationHandler->orbitalNavigator();
        orbNav.resetVelocities();
        _lastTransforms.clear();
        _lastTransforms.push(Transforms {
            .rotation = glm::dquat(),
            .scaling = 0.0,
            .timestamp = inputs.back().latestInput().timestamp
        });
    }

    return isInside;
}

void TouchInteraction::endDirectControl() {
    if (!_enableVelocity || _lastTransforms.size() < 2) {
        return;
    }

    double dt = _lastTransforms.back().timestamp - _lastTransforms.front().timestamp;

    if (dt == 0.0) {
        return;
    }

    glm::dquat rotation;
    double scaling = 0.0;

    for (size_t i = 0; i < _lastTransforms.size(); ++i) {
        rotation = _lastTransforms.at(i).rotation * rotation;
        scaling = _lastTransforms.at(i).scaling * scaling;
    }

    if (_useSphericalDisplay) {
        rotation = glm::dquat();
        scaling = 0.0;
        for (size_t i = 0; i < _lastTransforms.size(); ++i) {
            rotation = rotation * _lastTransforms.at(i).rotation;
        }
    }

    glm::dvec3 angularVel = glm::axis(rotation) * glm::angle(rotation) / dt;

    interaction::OrbitalNavigator& orbNav = global::navigationHandler->orbitalNavigator();
    // components x and y are swapped because angularVel refers to rotation around the x
    // and y axes of the camera, while the rotation velocity refers to
    // x=panning and y=tilting. This is mildly confusing.
    // @see OrbitalNavigator::rotateAroundAnchorUp() and translateHorizontally().
    orbNav.touchStates().setGlobalRotationVelocity(glm::dvec2(angularVel.y, angularVel.x));
    orbNav.touchStates().setGlobalRollVelocity(glm::dvec2(angularVel.z, 0.0));
    orbNav.touchStates().setTruckMovementVelocity(glm::dvec2(0.0, scaling));

    if (_useSphericalDisplay) {
        // on spherical displays the rotation is applied locally
        orbNav.touchStates().setGlobalRotationVelocity(glm::dvec2(0.0, 0.0));
        orbNav.touchStates().setGlobalRollVelocity(glm::dvec2(0.0, 0.0));
        orbNav.touchStates().setTruckMovementVelocity(glm::dvec2(0.0, 0.0));
        orbNav.touchStates().setLocalRotationVelocity(glm::dvec2(angularVel.y, angularVel.x));
        orbNav.touchStates().setLocalRollVelocity(glm::dvec2(angularVel.z, 0.0));
    }
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
    // newVelocity *= _zoomSensitivityProportionalDist * zoomGain;

    return newVelocity;
}

// Reset all property values to default
void TouchInteraction::resetPropertiesToDefault() {
    _enableOrbit.set(true);
    _enablePinchZoom.set(true);
    _enableTapZoom.set(true);
    _enableRoll.set(true);
    _enablePan.set(false);
    _enableDirectManipulation.set(true);
    _doubleTapThreshold.set(300);
    _rollAngleThreshold.set(0.025f);
    _panDistanceThreshold.set(0.015f);
    _stationaryThreshold.set(0.0005f);
    _centroidStillThreshold.set(0.0018f);
    _directManipulationThreshold.set(5.f);
    _tapZoomFactor.set(0.2f);
    _zoomInBoundarySphereMultiplier.set(1.001f);
    _zoomOutBoundarySphereMultiplier.set(4e+27f);
    _zoomInLimit.set(-1.f);
    _zoomOutLimit.set(4e+27f);
}

void TouchInteraction::tap() {
    // _tap = true;
    //TODO
}

void TouchInteraction::setCamera(Camera* camera) {
    _camera = camera;

    if (_useSphericalDisplay) {
        const SceneGraphNode* anchor = global::navigationHandler->orbitalNavigator().anchorNode();
        _camera->setPositionVec3(anchor->worldPosition() + glm::dvec3(0.0, 10.0, 0.0));
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
