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

#ifndef _SteinerCompact_H_
#define _SteinerCompact_H_

#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

/**
 * @brief The SteinerCompact class. See [Filter documentation](@ref steinercompact) for details.
 */
class SteinerCompact : public AbstractFilter
{
  Q_OBJECT

public:
  SIMPL_SHARED_POINTERS(SteinerCompact)
  SIMPL_STATIC_NEW_MACRO(SteinerCompact)
  SIMPL_TYPE_MACRO_SUPER(SteinerCompact, AbstractFilter)

  virtual ~SteinerCompact();

  SIMPL_INSTANCE_STRING_PROPERTY(DataContainerName)

  SIMPL_INSTANCE_STRING_PROPERTY(CellAttributeMatrixName)

  SIMPL_FILTER_PARAMETER(bool, VtkOutput)
  Q_PROPERTY(bool VtkOutput READ getVtkOutput WRITE setVtkOutput)

  SIMPL_FILTER_PARAMETER(QString, VtkFileName)
  Q_PROPERTY(QString VtkFileName READ getVtkFileName WRITE setVtkFileName)

  SIMPL_FILTER_PARAMETER(bool, TxtOutput)
  Q_PROPERTY(bool TxtOutput READ getTxtOutput WRITE setTxtOutput)

  SIMPL_FILTER_PARAMETER(QString, TxtFileName)
  Q_PROPERTY(QString TxtFileName READ getTxtFileName WRITE setTxtFileName)

  SIMPL_FILTER_PARAMETER(DataArrayPath, FeatureIdsArrayPath)
  Q_PROPERTY(DataArrayPath FeatureIdsArrayPath READ getFeatureIdsArrayPath WRITE setFeatureIdsArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, CellPhasesArrayPath)
  Q_PROPERTY(DataArrayPath CellPhasesArrayPath READ getCellPhasesArrayPath WRITE setCellPhasesArrayPath)

  SIMPL_FILTER_PARAMETER(int, Plane)
  Q_PROPERTY(int Plane READ getPlane WRITE setPlane)

  SIMPL_FILTER_PARAMETER(int, Sites)
  Q_PROPERTY(int Sites READ getSites WRITE setSites)

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
  SteinerCompact();
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  /**
* @brief  rose_of_intersections Counts the numbers of intersections of randomly placed lines at fixed directions with grain boundaries
*/
  virtual void rose_of_intersections(std::vector<std::vector<float>>& ROI);

  /**
* @brief find_intersections Counts the number of intersections of given line with grain boundaries
*/
  virtual void find_intersections(float line[3], int64_t z, int64_t dims[3], float res[3], std::vector<float>& numofintersections);

  /**
* @brief line_rectangle_intersections Finds the intersections of the line (line[0] * x + line[1] * y + line[2] = 0) with the rectangle (x1, y1, x2, y2)
*/
  virtual uint64_t line_rectangle_intersections(float xintersections[2], float yintersections[2], float rectangle[4], float line[3]);

  /**
* @brief find_one_site_vertices Identifies vertices of a site of the Steiner compact
*/
  virtual void find_one_site_vertices(std::vector<float> ROI, int64_t index, float vertices[4], float& length);

  /**
* @brief find_all_vertices Identifies vertices of all sites of the Steiner compact
*/
  virtual void find_all_vertices(std::vector<std::vector<float>>& vertices_x, std::vector<std::vector<float>>& vertices_y, std::vector<std::vector<float>>& radii,
                                 std::vector<std::vector<float>>& ROI);

  /**
* @brief output_vtk Writes the Steiner compact to a vtk file
*/
  virtual void output_vtk(std::vector<std::vector<float>>& vertices_x, std::vector<std::vector<float>>& vertices_y, std::vector<std::vector<float>>& radii, std::vector<std::vector<float>>& ROI);

  /**
* @brief output_txt Writes the Steiner compact to a text file
*/
  virtual void output_txt(std::vector<std::vector<float>>& vertices_x, std::vector<std::vector<float>>& vertices_y, std::vector<std::vector<float>>& ROI);

private:
  DEFINE_DATAARRAY_VARIABLE(int32_t, FeatureIds)
  DEFINE_DATAARRAY_VARIABLE(int32_t, CellPhases)

  SteinerCompact(const SteinerCompact&) = delete; // Copy Constructor Not Implemented
  void operator=(const SteinerCompact&) = delete; // Operator '=' Not Implemented
};

#endif /* _SteinerCompact_H_ */
