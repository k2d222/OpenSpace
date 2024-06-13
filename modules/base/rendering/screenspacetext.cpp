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

#include <modules/base/rendering/screenspacetext.h>

#include <openspace/documentation/documentation.h>
#include <openspace/documentation/verifier.h>
#include <openspace/engine/globals.h>
#include <openspace/rendering/renderengine.h>
#include <openspace/rendering/screenspacerenderable.h>
#include <ghoul/filesystem/filesystem.h>
#include <ghoul/font/font.h>
#include <ghoul/font/fontmanager.h>
#include <ghoul/font/fontrenderer.h>
#include <ghoul/io/texture/texturereader.h>
#include <ghoul/logging/logmanager.h>
#include <ghoul/opengl/openglstatecache.h>
#include <ghoul/opengl/programobject.h>
#include <ghoul/opengl/texture.h>
#include <ghoul/opengl/textureconversion.h>
#include <optional>

namespace {
    constexpr openspace::properties::Property::PropertyInfo ColorInfo = {
        "Color",
        "Color",
        "The text color.",
        openspace::properties::Property::Visibility::NoviceUser
    };

    constexpr openspace::properties::Property::PropertyInfo FontSizeInfo = {
        "FontSize",
        "Font Size",
        "The font size (in points) for the label.",
        openspace::properties::Property::Visibility::User
    };

    constexpr openspace::properties::Property::PropertyInfo TextInfo = {
        "Text",
        "Text",
        "Sets the text to be displayed on this screen space plane.",
        openspace::properties::Property::Visibility::NoviceUser
    };

    // This `ScreenSpaceRenderable` can be used to display text in screen space.
    struct [[codegen::Dictionary(ScreenSpaceText)]] Parameters {
        // [[codegen::verbatim(ColorInfo.description)]]
        std::optional<glm::vec3> color [[codegen::color()]];

        // [[codegen::verbatim(FontSizeInfo.description)]]
        std::optional<float> fontSize;

        // [[codegen::verbatim(TextInfo.description)]]
        std::optional<std::string> text;
    };
#include "screenspacetext_codegen.cpp"
} // namespace

namespace openspace {

documentation::Documentation ScreenSpaceText::Documentation() {
    return codegen::doc<Parameters>("base_screenspace_text");
}

ScreenSpaceText::ScreenSpaceText(const ghoul::Dictionary& dictionary)
    : ScreenSpaceFramebuffer(dictionary)
    , _color(ColorInfo, glm::vec3(1.f), glm::vec3(0.f), glm::vec3(1.f))
    , _fontSize(FontSizeInfo, 50.f, 1.f, 100.f)
    , _text(TextInfo, "")
{
    const Parameters p = codegen::bake<Parameters>(dictionary);

    std::string identifier;
    if (dictionary.hasValue<std::string>(KeyIdentifier)) {
        identifier = dictionary.value<std::string>(KeyIdentifier);
    }
    else {
        identifier = "ScreenSpaceText";
    }
    identifier = makeUniqueIdentifier(identifier);
    setIdentifier(identifier);

    _color = p.color.value_or(_color);
    _color.setViewOption(properties::Property::ViewOptions::Color);
    addProperty(_color);

    _fontSize = p.fontSize.value_or(_fontSize);
    _fontSize.onChange([this]() {
        _font = global::fontManager->font(
            "Mono",
            _fontSize,
            ghoul::fontrendering::FontManager::Outline::Yes,
            ghoul::fontrendering::FontManager::LoadGlyphs::No
        );
    });
    addProperty(_fontSize);

    // @TODO (emmbr, 2021-05-31): Temporarily set as read only, to avoid errors from font
    // rendering/loading
    _fontSize.setReadOnly(true);

    _text = p.text.value_or(_text);
    addProperty(_text);

    _text.onChange(std::bind(&ScreenSpaceText::resizeFramebuffer, this));
}

void ScreenSpaceText::resizeFramebuffer() {
    if (!_font) {
        return;
    }

    const std::string text = _text.value();
    const auto box = _font->boundingBox(text);
    const size_t lines = std::count(text.begin(), text.end(), '\n') + 1;
    const glm::uvec2 res = glm::ceil(box);

    if (res.x > 0 && res.y > 0) {
        _lines = lines;
        setResolution(res);
    }
}

void ScreenSpaceText::renderText() {
    glDepthMask(true);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    float y = float(_lines - 1) * _font->height() + (_font->hasOutline() ? 1.f : 0.f);
    const glm::vec3 transformedPos = glm::vec3(0.f, y, 0.f);

    glm::vec4 textColor = glm::vec4(glm::vec3(_color), 1.f);

    ghoul::fontrendering::FontRenderer::defaultRenderer().render(
        *_font,
        transformedPos,
        _text.value(),
        textColor
    );

    global::renderEngine->openglStateCache().resetBlendState();
    global::renderEngine->openglStateCache().resetDepthState();
}

bool ScreenSpaceText::initializeGL() {
    if (!ScreenSpaceFramebuffer::initializeGL())
        return false;
    

    if (_font == nullptr) {
        _font = global::fontManager->font(
            "Mono",
            _fontSize,
            ghoul::fontrendering::FontManager::Outline::Yes,
            ghoul::fontrendering::FontManager::LoadGlyphs::No
        );
    }

    resizeFramebuffer();

    addRenderFunction(std::bind(&ScreenSpaceText::renderText, this));

    return true;
}

} // namespace openspace
