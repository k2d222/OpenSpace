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

#include <modules/base/rendering/screenspaceframebuffer.h>

#include <openspace/documentation/documentation.h>
#include <openspace/engine/globals.h>
#include <openspace/engine/windowdelegate.h>
#include <openspace/rendering/renderengine.h>
#include <ghoul/opengl/framebufferobject.h>
#include <ghoul/opengl/openglstatecache.h>
#include <ghoul/opengl/textureunit.h>

namespace {
    constexpr openspace::properties::Property::PropertyInfo ResolutionInfo = {
        "Resolution",
        "Resolution",
        "This value sets the resolution in pixels of the framebuffer texture."
        "The default resolution is the active window resolution.",
        openspace::properties::Property::Visibility::AdvancedUser
    };

    struct [[codegen::Dictionary(ScreenSpaceFramebuffer)]] Parameters {
        // [[codegen::verbatim(ResolutionInfo.description)]]
        std::optional<glm::ivec2> resolution;
    };
#include "screenspaceframebuffer_codegen.cpp"
} // namespace

namespace openspace {

documentation::Documentation ScreenSpaceFramebuffer::Documentation() {
    return codegen::doc<Parameters>("base_screenspace_framebuffer");
}

ScreenSpaceFramebuffer::ScreenSpaceFramebuffer(const ghoul::Dictionary& dictionary)
    : ScreenSpaceRenderable(dictionary)
    , _resolution(ResolutionInfo, glm::ivec2(0), glm::ivec2(0), glm::ivec2(16384))
{
    documentation::testSpecificationAndThrow(
        Documentation(),
        dictionary,
        "ScreenSpaceFramebuffer"
    );

    int iIdentifier = 0;
    if (_identifier.empty()) {
        iIdentifier = id();

        if (iIdentifier == 0) {
            setIdentifier("ScreenSpaceFramebuffer");
        }
        else {
            setIdentifier("ScreenSpaceFramebuffer" + std::to_string(iIdentifier));
        }
    }

    if (_guiName.empty()) {
        // Adding an extra space to the user-facing name as it looks nicer
        setGuiName("ScreenSpaceFramebuffer " + std::to_string(iIdentifier));
    }

    const Parameters p = codegen::bake<Parameters>(dictionary);

    const glm::uvec2 res = p.resolution.value_or(global::windowDelegate->currentDrawBufferResolution());

    addProperty(_resolution);

    _resolution.onChange([this]() {
        if (_texture) {
            _texture->setDimensions(glm::uvec3(_resolution.value(), 1));
            _texture->uploadTexture();
            _texture->purgeFromRAM();
        }
    });

    _resolution = res;
}

ScreenSpaceFramebuffer::~ScreenSpaceFramebuffer() {}

void ScreenSpaceFramebuffer::initializeGL() {
    ScreenSpaceRenderable::initializeGL();

    createFramebuffer();
}

void ScreenSpaceFramebuffer::deinitializeGL() {
    _framebuffer->activate();
    _framebuffer->detachAll();
    ghoul::opengl::FramebufferObject::deactivate();
    removeAllRenderFunctions();

    ScreenSpaceRenderable::deinitializeGL();
}

void ScreenSpaceFramebuffer::render(const RenderData& renderData) {
    const glm::ivec2& resolution = _resolution.value();

    if (!_renderFunctions.empty()) {
        const GLint defaultFBO = ghoul::opengl::FramebufferObject::getActiveObject();
        _framebuffer->activate();

        std::array<GLint, 4> viewport;
        global::renderEngine->openglStateCache().viewport(viewport.data());

        glViewport(
            static_cast<GLint>(0),
            static_cast<GLint>(0),
            static_cast<GLsizei>(resolution.x),
            static_cast<GLsizei>(resolution.y)
        );

        glClearColor(0.f, 0.f, 0.f, 0.f);
        glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        for (const RenderFunction& renderFunction : _renderFunctions) {
            renderFunction();
        }

        global::renderEngine->openglStateCache().setViewportState(viewport.data());

        glBindFramebuffer(GL_FRAMEBUFFER, defaultFBO);

        _objectSize = resolution;

        const glm::mat4 globalRotation = globalRotationMatrix();
        const glm::mat4 translation = translationMatrix();
        const glm::mat4 localRotation = localRotationMatrix();
        const glm::mat4 scale = glm::scale(
            scaleMatrix(),
            // glm::vec3((1.f / xratio), (1.f / yratio), 1.f)
            glm::vec3(1.f)
        );
        const glm::mat4 modelTransform = globalRotation*translation*localRotation*scale;
        draw(modelTransform, renderData);
    }
}

bool ScreenSpaceFramebuffer::isReady() const {
    return _shader && _texture;
}

void ScreenSpaceFramebuffer::setResolution(glm::uvec2 resolution) {
    _resolution = std::move(resolution);
    if (_texture) {
        _texture->setDimensions(glm::uvec3(_resolution.value(), 1));
        _texture->uploadTexture();
        _texture->purgeFromRAM();
    }
}

void ScreenSpaceFramebuffer::addRenderFunction(RenderFunction renderFunction) {
    _renderFunctions.push_back(std::move(renderFunction));
}

void ScreenSpaceFramebuffer::removeAllRenderFunctions() {
    _renderFunctions.clear();
}

void ScreenSpaceFramebuffer::createFramebuffer() {
    _framebuffer = std::make_unique<ghoul::opengl::FramebufferObject>();
    _framebuffer->activate();
    _texture = std::make_unique<ghoul::opengl::Texture>(
        glm::uvec3(_resolution.value(), 1),
        GL_TEXTURE_2D
    );

    _texture->uploadTexture();
    _texture->setFilter(ghoul::opengl::Texture::FilterMode::Linear);
    _texture->purgeFromRAM();
    _framebuffer->attachTexture(_texture.get(), GL_COLOR_ATTACHMENT0);
    ghoul::opengl::FramebufferObject::deactivate();
}

int ScreenSpaceFramebuffer::id() {
    static int id = 0;
    return id++;
}

void ScreenSpaceFramebuffer::bindTexture() {
    _texture->bind();
}

} //namespace openspace
