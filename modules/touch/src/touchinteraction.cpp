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
    // Calculated with two vectors with known diff in length, then
    // projDiffLength/diffLength.
    , _velocity{ glm::dvec2(0.0), 0.0, 0.0, glm::dvec2(0.0) }
    , _sensitivity{ glm::dvec2(0.08, 0.045), 12.0, 2.75, glm::dvec2(0.08, 0.045) }
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

void TouchInteraction::updateStateFromInput(const std::vector<TouchInputHolder>& inputs)
{
    size_t numFingers = inputs.size();

#ifdef TOUCH_DEBUG_PROPERTIES
    _debugProperties.nFingers = numFingers;
#endif

    if (numFingers == 0) {
        // No fingers, no input (note that this function should not even be called then)
        _anchor = nullptr;
        return;
    }

    else if(numFingers == 1 && !_anchor) {
        _anchor = global::navigationHandler->orbitalNavigator().anchorNode();
        _startPose = {
            .position = _camera->positionVec3(),
            .rotation = _camera->rotationQuaternion(),
        };
    }

    directControl(inputs);
    updateNodeSurfacePoints(inputs);
}

glm::dvec3 TouchInteraction::unprojectTouchOnSphere(const TouchInput& input) const {
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

void TouchInteraction::directControl(const std::vector<TouchInputHolder>& inputs) {
#ifdef TOUCH_DEBUG_PROPERTIES
    LINFO("DirectControl");
#endif

    if (inputs.size() != 1 || !_anchor || !_camera) {
        return;
    }

    const TouchInputHolder& input = inputs[0];

    interaction::OrbitalNavigator& orbNav = global::navigationHandler->orbitalNavigator();

    const glm::dvec3 startWorldSpace = unprojectTouchOnSphere(input.firstInput());
    const glm::dvec3 endWorldSpace = unprojectTouchOnSphere(input.latestInput());
    const glm::dvec3 spherePosition = _anchor->worldPosition();

    if (startWorldSpace == glm::dvec3(0.0) || endWorldSpace == glm::dvec3(0.0)) {
        return;
    }
    
    glm::vec3 startVec = glm::normalize(startWorldSpace - spherePosition);
    glm::vec3 endVec = glm::normalize(endWorldSpace - spherePosition);

    glm::vec3 rotAxis = glm::normalize(glm::cross(startVec, endVec));
    float rotAngle = -glm::acos(glm::dot(startVec, endVec));
    glm::dquat rotQuat = glm::angleAxis(rotAngle, rotAxis);

    if (glm::any(glm::isnan(rotAxis))) {
        return;
    }

    _camera->setPose(_startPose);
    _camera->setPositionVec3(rotQuat * (_camera->positionVec3() - spherePosition) + spherePosition);
    _camera->rotate(rotQuat);

    // Mark that a camera interaction happened
    orbNav.updateOnCameraInteraction();

}

void TouchInteraction::updateNodeSurfacePoints(const std::vector<TouchInputHolder>& inputs)
{
    _selectedNodeSurfacePoints.clear();

    const SceneGraphNode* anchor =
        global::navigationHandler->orbitalNavigator().anchorNode();
    SceneGraphNode* node = sceneGraphNode(anchor->identifier());

    // Check if current anchor is valid for direct touch
    TouchModule* module = global::moduleEngine->module<TouchModule>();

    bool isDirectTouchRenderable = node->renderable() &&
        module->isDefaultDirectTouchType(node->renderable()->typeAsString());

    if (!(node->supportsDirectInteraction() || isDirectTouchRenderable)) {
        return;
    }

    glm::dquat camToWorldSpace = _camera->rotationQuaternion();
    glm::dvec3 camPos = _camera->positionVec3();
    std::vector<DirectInputSolver::SelectedBody> surfacePoints;

    for (const TouchInputHolder& inputHolder : inputs) {
        // Normalized -1 to 1 coordinates on screen
        double xCo = 2 * (inputHolder.latestInput().x - 0.5);
        double yCo = -2 * (inputHolder.latestInput().y - 0.5);
        glm::dvec3 cursorInWorldSpace = camToWorldSpace *
            glm::dvec3(glm::inverse(_camera->projectionMatrix()) *
            glm::dvec4(xCo, yCo, -1.0, 1.0));
        glm::dvec3 rayDir = glm::normalize(cursorInWorldSpace);

        size_t id = inputHolder.fingerId();

        // Compute positions on anchor node, by checking if touch input
        // intersect interaction sphere
        double intersectionDist = 0.0;
        const bool intersected = glm::intersectRaySphere(
            camPos,
            rayDir,
            node->worldPosition(),
            node->interactionSphere() * node->interactionSphere(),
            intersectionDist
        );

        if (intersected) {
            glm::dvec3 intersectionPos = camPos + rayDir * intersectionDist;
            glm::dvec3 pointInModelView = glm::inverse(node->worldRotationMatrix()) *
                                            (intersectionPos - node->worldPosition());

            // Note that node is saved as the direct input solver was initially
            // implemented to handle touch contact points on multiple nodes
            surfacePoints.push_back({ id, node, pointInModelView });
        }
    }

    _selectedNodeSurfacePoints = std::move(surfacePoints);
}

double TouchInteraction::computeConstTimeDecayCoefficient(double velocity) {
    constexpr double postDecayVelocityTarget = 1e-6;
    const double stepsToDecay = _constTimeDecay_secs / _frameTimeAvg.averageFrameTime();

    if (stepsToDecay > 0.0 && std::abs(velocity) > postDecayVelocityTarget) {
        return std::pow(postDecayVelocityTarget / std::abs(velocity), 1.0 / stepsToDecay);
    }
    return 1.0;
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

bool TouchInteraction::hasNonZeroVelocities() const {
    glm::dvec2 sum = _velocity.orbit;
    sum += glm::dvec2(_velocity.zoom, 0.0);
    sum += glm::dvec2(_velocity.roll, 0.0);
    sum += _velocity.pan;
    // Epsilon size based on that even if no interaction is happening,
    // there might still be some residual velocity
    return glm::length(sum) > 0.001;
}

// Main update call, calculates the new orientation and position for the camera depending
// on _vel and dt. Called every frame
void TouchInteraction::step(double dt, bool directTouch) {
    return;
    using namespace glm;

    if (!(directTouch || hasNonZeroVelocities())) {
        // No motion => don't update the camera
        return;
    }

    const SceneGraphNode* anchor =
        global::navigationHandler->orbitalNavigator().anchorNode();

    if (anchor && _camera) {
        // Create variables from current state
        dvec3 camPos = _camera->positionVec3();
        const dvec3 centerPos = anchor->worldPosition();

        dvec3 directionToCenter = normalize(centerPos - camPos);
        const dvec3 centerToCamera = camPos - centerPos;
        const dvec3 lookUp = _camera->lookUpVectorWorldSpace();
        const dvec3 camDirection = _camera->viewDirectionWorldSpace();

        // Make a representation of the rotation quaternion with local and global
        // rotations. To avoid problem with lookup in up direction
        const dmat4 lookAtMat = lookAt(
            dvec3(0.0, 0.0, 0.0),
            directionToCenter,
            normalize(camDirection + lookUp)
        );
        dquat globalCamRot = normalize(quat_cast(inverse(lookAtMat)));
        dquat localCamRot = inverse(globalCamRot) * _camera->rotationQuaternion();

        const double interactionSphere = anchor->interactionSphere();

        // Check if camera is within distance for direct manipulation to be applicable
        // if (interactionSphere > 0.0 && _enableDirectManipulation) {
        //     const double distance =
        //         std::max(length(centerToCamera) - interactionSphere, 0.0);
        //     const double maxDistance = interactionSphere * _directTouchDistanceThreshold;
        //     _isWithinDirectTouchDistance = distance <= maxDistance;
        // }
        // else {
        //     _isWithinDirectTouchDistance = false;
        // }

        {
            // Roll
            const dquat camRollRot = angleAxis(_velocity.roll * dt, dvec3(0.0, 0.0, 1.0));
            localCamRot = localCamRot * camRollRot;
        }
        {
            // Panning (local rotation)
            const dvec3 eulerAngles(_velocity.pan.y * dt, _velocity.pan.x * dt, 0.0);
            const dquat rotationDiff = dquat(eulerAngles);
            localCamRot = localCamRot * rotationDiff;
        }
        {
            // Orbit (global rotation)
            const dvec3 eulerAngles(_velocity.orbit.y * dt, _velocity.orbit.x * dt, 0.0);
            const dquat rotationDiffCamSpace = dquat(eulerAngles);

            const dquat rotationDiffWorldSpace = globalCamRot * rotationDiffCamSpace *
                                                 inverse(globalCamRot);
            const dvec3 rotationDiffVec3 = centerToCamera * rotationDiffWorldSpace -
                                           centerToCamera;
            camPos += rotationDiffVec3;

            const dvec3 centerToCam = camPos - centerPos;
            directionToCenter = normalize(-centerToCam);
            const dvec3 lookUpWhenFacingCenter = globalCamRot *
                                           dvec3(_camera->lookUpVectorCameraSpace());

            const dmat4 lookAtMatrix = lookAt(
                dvec3(0.0),
                directionToCenter,
                lookUpWhenFacingCenter
            );
            globalCamRot = normalize(quat_cast(inverse(lookAtMatrix)));
        }
        {
            // Zooming

            // This is a rough estimate of the node surface
            // If nobody has set another zoom in limit, use this as default zoom in bounds
            double zoomInBounds = interactionSphere * _zoomInBoundarySphereMultiplier;
            bool isZoomInLimitSet = (_zoomInLimit.value() >= 0.0);

            if (isZoomInLimitSet && _zoomInLimit.value() < zoomInBounds) {
                // If zoom in limit is less than the estimated node radius we need to
                // make sure we do not get too close to possible height maps
                SurfacePositionHandle posHandle = anchor->calculateSurfacePositionHandle(
                    camPos
                );
                glm::dvec3 centerToActualSurfaceModelSpace =
                    posHandle.centerToReferenceSurface +
                    posHandle.referenceSurfaceOutDirection * posHandle.heightToSurface;
                glm::dvec3 centerToActualSurface = glm::dmat3(anchor->modelTransform()) *
                    centerToActualSurfaceModelSpace;
                const double nodeRadius = length(centerToActualSurface);

                // Because of heightmaps we need to ensure we don't go through the surface
                if (_zoomInLimit.value() < nodeRadius) {
#ifdef TOUCH_DEBUG_PROPERTIES
                    LINFO(std::format(
                        "Zoom In limit should be larger than anchor "
                        "center to surface, setting it to {}", zoomInBounds
                    ));
#endif
                    zoomInBounds = _zoomInLimit.value();
                }
            }

            double zoomOutBounds = std::min(
                interactionSphere *_zoomOutBoundarySphereMultiplier,
                _zoomOutLimit.value()
            );

            // Make sure zoom in limit is not larger than zoom out limit
            if (zoomInBounds > zoomOutBounds) {
               LWARNING(std::format(
                   "Zoom In Limit should be smaller than Zoom Out Limit",
                    zoomOutBounds
               ));
            }
            const double currentPosDistance = length(centerToCamera);

            // Apply the velocity to update camera position
            double zoomVelocity = _velocity.zoom;
            if (!directTouch) {
                const double distanceFromSurface =
                    length(currentPosDistance) - anchor->interactionSphere();
                if (distanceFromSurface > 0.1) {
                    const double ratioOfDistanceToNodeVsSurf =
                        length(currentPosDistance) / distanceFromSurface;
                    if (ratioOfDistanceToNodeVsSurf > _zoomSensitivityDistanceThreshold) {
                        zoomVelocity *= pow(
                            std::abs(distanceFromSurface),
                            static_cast<float>(_zoomSensitivityExponential)
                        );
                    }
                }
                else {
                    zoomVelocity = 1.0;
                }
            }

            const glm::dvec3 zoomDistanceInc = directionToCenter * zoomVelocity * dt;
            const double newPosDistance = length(centerToCamera + zoomDistanceInc);

            // Possible with other navigations performed outside touch interaction
            const bool currentPosViolatingZoomOutLimit =
                (currentPosDistance >= zoomOutBounds);
            const bool willNewPositionViolateZoomOutLimit =
                (newPosDistance >= zoomOutBounds);
            bool willNewPositionViolateZoomInLimit =
                (newPosDistance < zoomInBounds);
            bool willNewPositionViolateDirection =
                (currentPosDistance <= length(zoomDistanceInc));

            if (!willNewPositionViolateZoomInLimit &&
                !willNewPositionViolateDirection &&
                !willNewPositionViolateZoomOutLimit)
            {
                camPos += zoomDistanceInc;
            }
            else if (currentPosViolatingZoomOutLimit) {
#ifdef TOUCH_DEBUG_PROPERTIES
                LINFO(std::format(
                    "You are outside zoom out {} limit, only zoom in allowed",
                    zoomOutBounds
                ));
#endif
                // Only allow zooming in if you are outside the zoom out limit
                if (newPosDistance < currentPosDistance) {
                    camPos += zoomDistanceInc;
                }
            }
            else {
#ifdef TOUCH_DEBUG_PROPERTIES
                LINFO("Zero the zoom velocity close to surface");
#endif
                _velocity.zoom = 0.0;
            }
        }

        decelerate(dt);

        // @TODO (emmbr, 2023-02-08) This is ugly, but for now prevents jittering
        // when zooming in closer than the orbital navigator allows. Long term, we
        // should make the touch interaction tap into the orbitalnavigator and let that
        // do the updating of the camera, instead of handling them separately. Then we
        // would keep them in sync and avoid duplicated camera updating code.
        auto orbitalNavigator = global::navigationHandler->orbitalNavigator();
        camPos = orbitalNavigator.pushToSurfaceOfAnchor(camPos);

        // @TODO (emmbr, 2023-02-08) with the line above, the ZoomInLimit might not be
        // needed anymore. We should make it so that just the limit properties in the
        // OrbitalNavigator is actually needed, and don't have duplicates

        // Update the camera state
        _camera->setPositionVec3(camPos);
        _camera->setRotation(globalCamRot * localCamRot);

        // Mark that a camera interaction happened
        global::navigationHandler->orbitalNavigator().updateOnCameraInteraction();

#ifdef TOUCH_DEBUG_PROPERTIES
        //Show velocity status every N frames
        if (++stepVelUpdate >= 60) {
            stepVelUpdate = 0;
            LINFO(std::format(
                "DistToFocusNode {} stepZoomVelUpdate {}",
                length(centerToCamera), _vel.zoom
            ));
        }
#endif

        // _tap = false;
        // _doubleTap = false;
        // _zoomOutTap = false;
    }
}

// Decelerate velocities, called a set number of times per second to dereference it from
// frame time
// Example:
// Assume: frequency = 0.01, dt = 0.05 (200 fps), _timeSlack = 0.0001
// times = floor((0.05 + 0.0001) / 0.01) = 5
// _timeSlack = 0.0501 % 0.01 = 0.01
void TouchInteraction::decelerate(double dt) {
    _frameTimeAvg.updateWithNewFrame(dt);
    double expectedFrameTime = _frameTimeAvg.averageFrameTime();

    // Number of times velocities should decelerate, depending on chosen frequency and
    // time slack over from last frame
    int times = static_cast<int>((dt + _timeSlack) / expectedFrameTime);
    // Save the new time slack for the next frame
    _timeSlack = fmod((dt + _timeSlack), expectedFrameTime) * expectedFrameTime;

    //Ensure the number of times to apply the decay coefficient is valid
    times = std::min(times, 1);

    _velocity.orbit *= computeDecayCoeffFromFrametime(_constTimeDecayCoeff.orbit, times);
    _velocity.roll  *= computeDecayCoeffFromFrametime(_constTimeDecayCoeff.roll,  times);
    _velocity.pan   *= computeDecayCoeffFromFrametime(_constTimeDecayCoeff.pan,   times);
    _velocity.zoom  *= computeDecayCoeffFromFrametime(_constTimeDecayCoeff.zoom,  times);
}

// Called if all fingers are off the screen
void TouchInteraction::resetAfterInput() {
#ifdef TOUCH_DEBUG_PROPERTIES
    _debugProperties.nFingers = 0;
    _debugProperties.interactionMode = "None";
#endif
    // @TODO (emmbr 2023-02-03) Bring back feature that allows node to spin when
    // the direct manipulaiton finger is let go. Should implement this using the
    // orbitalnavigator's friction values. This also implies passing velocities to
    // the orbitalnavigator, instead of setting the camera directly as is currently
    // done in this class.

    _constTimeDecayCoeff.zoom = computeConstTimeDecayCoefficient(_velocity.zoom);
    _pinchInputs[0].clearInputs();
    _pinchInputs[1].clearInputs();

    _selectedNodeSurfacePoints.clear();
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

void TouchInteraction::resetVelocities() {
    _velocity.orbit = glm::dvec2(0.0);
    _velocity.zoom = 0.0;
    _velocity.roll = 0.0;
    _velocity.pan = glm::dvec2(0.0);
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
