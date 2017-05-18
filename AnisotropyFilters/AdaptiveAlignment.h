/* ============================================================================
* Copyright (c) 2016 Czech Academy of Sciences, Institute of Physics,
* Group of Bulk Nanomaterials and Interfaces, http://ams.fzu.cz
*
* Redistribution and use in source and binary forms, with or without modification,
* are permitted provided that the following conditions are met:
*
* Redistributions of source code must retain the above copyright notice, this
* list of conditions and the following disclaimer.
*
* Redistributions in binary form must reproduce the above copyright notice, this
* list of conditions and the following disclaimer in the documentation and/or
* other materials provided with the distribution.
*
* Neither the name of the Czech Academy of Sciences, nor the names of its
* contributors may be used to endorse or promote products derived from this
* software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
* USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* The code contained herein was partially funded by the followig grants:
*    Czech Science Foundation (GA CR), project no. GBP108/12/G043
*    Czech Ministry of Education, Youth and Sports (MSMT), project no. LM2015087
*
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


#ifndef _AdaptiveAlignment_H_
#define _AdaptiveAlignment_H_

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/AbstractFilter.h"
#include "SIMPLib/Common/SIMPLibSetGetMacros.h"

#include "Anisotropy/AnisotropyConstants.h"

/**
* @brief The AdaptiveAlignment class. This class serves as a superclass for other classes
* in the Reconstruction plugin.
*/
class AdaptiveAlignment : public AbstractFilter
{
    Q_OBJECT
  public:
    SIMPL_SHARED_POINTERS(AdaptiveAlignment)
    SIMPL_STATIC_NEW_MACRO(AdaptiveAlignment)
    SIMPL_TYPE_MACRO_SUPER(AdaptiveAlignment, AbstractFilter)

    virtual ~AdaptiveAlignment();

    SIMPL_INSTANCE_STRING_PROPERTY(DataContainerName)

    SIMPL_INSTANCE_STRING_PROPERTY(CellAttributeMatrixName)

    SIMPL_INSTANCE_PROPERTY(bool, WriteAlignmentShifts)
    Q_PROPERTY(bool WriteAlignmentShifts READ getWriteAlignmentShifts WRITE setWriteAlignmentShifts)

    SIMPL_INSTANCE_STRING_PROPERTY(AlignmentShiftFileName)
    Q_PROPERTY(QString AlignmentShiftFileName READ getAlignmentShiftFileName WRITE setAlignmentShiftFileName)

    /////////////// new:

  SIMPL_FILTER_PARAMETER(int, GlobalCorrection)
  Q_PROPERTY(int GlobalCorrection READ getGlobalCorrection WRITE setGlobalCorrection)

    SIMPL_FILTER_PARAMETER(QString, InputPath)
    Q_PROPERTY(QString InputPath READ getInputPath WRITE setInputPath)

    SIMPL_FILTER_PARAMETER(float, ShiftX)
    Q_PROPERTY(float ShiftX READ getShiftX WRITE setShiftX)

    SIMPL_FILTER_PARAMETER(float, ShiftY)
    Q_PROPERTY(float ShiftY READ getShiftY WRITE setShiftY)

    SIMPL_FILTER_PARAMETER(DataArrayPath, ImageDataArrayPath)
    Q_PROPERTY(DataArrayPath ImageDataArrayPath READ getImageDataArrayPath WRITE setImageDataArrayPath)

    SIMPL_FILTER_PARAMETER(QString, NewCellArrayName)
    Q_PROPERTY(QString NewCellArrayName READ getNewCellArrayName WRITE setNewCellArrayName)

    SIMPL_FILTER_PARAMETER(float, MinRadius)
    Q_PROPERTY(float MinRadius READ getMinRadius WRITE setMinRadius)

    SIMPL_FILTER_PARAMETER(float, MaxRadius)
    Q_PROPERTY(float MaxRadius READ getMaxRadius WRITE setMaxRadius)

    SIMPL_FILTER_PARAMETER(int, NumberCircles)
    Q_PROPERTY(int NumberCircles READ getNumberCircles WRITE setNumberCircles)


    /**
    * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
    */
    virtual const QString getCompiledLibraryName();

    /**
  * @brief getBrandingString Returns the branding string for the filter, which is a tag
  * used to denote the filter's association with specific plugins
  * @return Branding string
  */
    virtual const QString getBrandingString();

    /**
  * @brief getFilterVersion Returns a version string for this filter. Default
  * value is an empty string.
  * @return
  */
    virtual const QString getFilterVersion();

    /**
  * @brief newFilterInstance Reimplemented from @see AbstractFilter class
  */
    virtual AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters);

    /**
  * @brief getGroupName Reimplemented from @see AbstractFilter class
  */
    virtual const QString getGroupName();

    /**
  * @brief getSubGroupName Reimplemented from @see AbstractFilter class
  */
    virtual const QString getSubGroupName();

    /**
  * @brief getHumanLabel Reimplemented from @see AbstractFilter class
  */
    virtual const QString getHumanLabel();

    /**
  * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
  */
    virtual void setupFilterParameters();

    /**
  * @brief readFilterParameters Reimplemented from @see AbstractFilter class
  */
    virtual void readFilterParameters(AbstractFilterParametersReader* reader, int index);

    /**
  * @brief execute Reimplemented from @see AbstractFilter class
  */
    virtual void execute();

    /**
  * @brief preflight Reimplemented from @see AbstractFilter class
  */
    virtual void preflight();

  signals:
    /**
  * @brief updateFilterParameters Emitted when the Filter requests all the latest Filter parameters
  * be pushed from a user-facing control (such as a widget)
  * @param filter Filter instance pointer
  */
    void updateFilterParameters(AbstractFilter* filter);

    /**
  * @brief parametersChanged Emitted when any Filter parameter is changed internally
  */
    void parametersChanged();

    /**
  * @brief preflightAboutToExecute Emitted just before calling dataCheck()
  */
    void preflightAboutToExecute();

    /**
  * @brief preflightExecuted Emitted just after calling dataCheck()
  */
    void preflightExecuted();

  protected:
    AdaptiveAlignment();

    /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
    void dataCheck();

    /**
   * @brief Initializes all the private instance variables.
   */
    void initialize();

    /**
  * @brief create_array_for_flattened_image creates a temporary array
  */
    virtual void create_array_for_flattened_image();

    /**
  * @brief delete_array_for_flattened_image deletes a temporary array
  */
    virtual void delete_array_for_flattened_image();

    /**
  * @brief flatten_image for converting an RGB image to grayscale
  */
    virtual void flatten_image();

    /**
  * @brief find_circles for identification of two circles for alignment in SEM images by Hough transform
  */
    virtual bool find_calibrating_circles();

    /**
  * @brief find_rectangles for identification of the green rectangle
  */
    virtual bool find_rectangles();

    /**
  * @brief find_rectangles for identification of the edge of the specimen
  */
    virtual bool find_interface_edges();

    /**
  * @brief estimate_shifts for estimation of the shifts from SEM images
  */
    virtual void estimate_shifts_from_images(std::vector<float>& xneedshifts, std::vector<float>& yneedshifts);

    /**
  * @brief find_shifts Determines the x and y shifts to register a stacked 3D volume
  * @param xshifts Vector of integer shifts in x direction
  * @param yshifts Vector of integer shifts in y direction
  */
    virtual void find_shifts(std::vector<int64_t>& xshifts, std::vector<int64_t>& yshifts, std::vector<float>& xneedshifts, std::vector<float>& yneedshifts);

  private:

    DEFINE_DATAARRAY_VARIABLE(AnisotropyConstants::DefaultPixelType, ImageData)
    DEFINE_DATAARRAY_VARIABLE(AnisotropyConstants::DefaultPixelType, FlatImageData)

    DEFINE_DATAARRAY_VARIABLE(AnisotropyConstants::DefaultPixelType, NewCellArray) // delete this

    std::vector<std::vector<uint64_t>> m_RectangleCorners;
    std::vector<std::vector<float>> m_CalibratingCircles;
    std::vector<uint64_t> m_InterfaceEdges;

    AdaptiveAlignment(const AdaptiveAlignment&); // Copy Constructor Not Implemented
    void operator=(const AdaptiveAlignment&); // Operator '=' Not Implemented
};

#endif /* AdaptiveAlignment_H_ */
