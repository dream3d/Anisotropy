

#include "AnisotropyGuiPlugin.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AnisotropyGuiPlugin::AnisotropyGuiPlugin()
: AnisotropyPlugin()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AnisotropyGuiPlugin::~AnisotropyGuiPlugin() = default;

#include "Anisotropy/Gui/FilterParameterWidgets/RegisterKnownFilterParameterWidgets.cpp"
