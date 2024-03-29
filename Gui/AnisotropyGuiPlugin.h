#pragma once

#include "Anisotropy/AnisotropyPlugin.h"

class AnisotropyGuiPlugin : public AnisotropyPlugin
{
  Q_OBJECT
  Q_INTERFACES(ISIMPLibPlugin)
  Q_PLUGIN_METADATA(IID "net.bluequartz.dream3d.AnisotropyGuiPlugin")

public:
  AnisotropyGuiPlugin();
   ~AnisotropyGuiPlugin() override;
  
  /**
   * @brief Register all the filters with the FilterWidgetFactory
   */
  void registerFilterWidgets(FilterWidgetManager* fwm) override;
  

public:
  AnisotropyGuiPlugin(const AnisotropyGuiPlugin&) = delete;            // Copy Constructor Not Implemented
  AnisotropyGuiPlugin(AnisotropyGuiPlugin&&) = delete;                 // Move Constructor Not Implemented
  AnisotropyGuiPlugin& operator=(const AnisotropyGuiPlugin&) = delete; // Copy Assignment Not Implemented
  AnisotropyGuiPlugin& operator=(AnisotropyGuiPlugin&&) = delete;      // Move Assignment Not Implemented
};
