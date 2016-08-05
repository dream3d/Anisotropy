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

#include "SteinerCompact.h"

#include <fstream>

#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtCore/QFile>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/SIMPLibVersion.h"
#include "SIMPLib/Math/SIMPLibMath.h"
#include "SIMPLib/Common/ScopedFileMonitor.hpp"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/OutputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"

#include "SIMPLib/Utilities/SIMPLibRandom.h"
#include "SIMPLib/Geometry/ImageGeom.h"

#include "Anisotropy/AnisotropyConstants.h"
#include "Anisotropy/AnisotropyVersion.h"

// Include the MOC generated file for this class
#include "moc_SteinerCompact.cpp"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SteinerCompact::SteinerCompact() :
  AbstractFilter(),
  m_DataContainerName(SIMPL::Defaults::ImageDataContainerName),
  m_VtkOutput(true),
  m_VtkFileName(""),
  m_TxtOutput(false),
  m_TxtFileName(""),
  m_FeatureIdsArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::FeatureIds),
  m_CellPhasesArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::Phases),
  m_Plane(0),
  m_Sites(1),
  m_FeatureIds(NULL),
  m_CellPhases(NULL)
{
  setupFilterParameters();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
SteinerCompact::~SteinerCompact()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SteinerCompact::setupFilterParameters()
{
  FilterParameterVector parameters;
  QStringList linkedProps;
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Section Plane");
    parameter->setPropertyName("Plane");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(SteinerCompact, this, Plane));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(SteinerCompact, this, Plane));
    QVector<QString> choices;
    choices.push_back("XY");
    choices.push_back("XZ");
    choices.push_back("YZ");
    parameter->setChoices(choices);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Number Of Sites");
    parameter->setPropertyName("Sites");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(SteinerCompact, this, Sites));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(SteinerCompact, this, Sites));
    QVector<QString> choices;
    choices.push_back("8");
    choices.push_back("12");
    choices.push_back("16");
    choices.push_back("24");
    choices.push_back("36");
    parameter->setChoices(choices);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);
  }
  //parameters.push_back(SIMPL_NEW_INTEGER_FP("Number Of Sites", Sites, FilterParameter::Parameter, SteinerCompact));
  linkedProps.clear();
  linkedProps << "VtkFileName";
  parameters.push_back(LinkedBooleanFilterParameter::New("Graphical Output As .vtk", "VtkOutput", getVtkOutput(), FilterParameter::Parameter, SIMPL_BIND_SETTER(SteinerCompact, this, VtkOutput), SIMPL_BIND_GETTER(SteinerCompact, this, VtkOutput), linkedProps));
  parameters.push_back(OutputFileFilterParameter::New("Output Vtk File", "VtkFileName", getVtkFileName(), FilterParameter::Parameter, SIMPL_BIND_SETTER(SteinerCompact, this, VtkFileName), SIMPL_BIND_GETTER(SteinerCompact, this, VtkFileName), "*.vtk", "VTK Polydata"));
  linkedProps.clear();
  linkedProps << "TxtFileName";
  parameters.push_back(LinkedBooleanFilterParameter::New("Text Output As .txt", "TxtOutput", getTxtOutput(), FilterParameter::Parameter, SIMPL_BIND_SETTER(SteinerCompact, this, TxtOutput), SIMPL_BIND_GETTER(SteinerCompact, this, TxtOutput), linkedProps));
  parameters.push_back(OutputFileFilterParameter::New("Output Text File", "TxtFileName", getTxtFileName(), FilterParameter::Parameter, SIMPL_BIND_SETTER(SteinerCompact, this, TxtFileName), SIMPL_BIND_GETTER(SteinerCompact, this, TxtFileName), "*.txt", "Text"));
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, SIMPL::AttributeMatrixType::Cell, SIMPL::GeometryType::ImageGeometry);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature Ids", FeatureIdsArrayPath, FilterParameter::RequiredArray, SteinerCompact, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, SIMPL::AttributeMatrixType::Cell, SIMPL::GeometryType::ImageGeometry);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Phases", CellPhasesArrayPath, FilterParameter::RequiredArray, SteinerCompact, req));
  }

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SteinerCompact::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setFeatureIdsArrayPath(reader->readDataArrayPath("FeatureIdsArrayPath", getFeatureIdsArrayPath()));
  setCellPhasesArrayPath(reader->readDataArrayPath("CellPhasesArrayPath", getCellPhasesArrayPath()));
  setPlane(reader->readValue("Plane", getPlane()));
  setSites(reader->readValue("Sites", getSites()));
  setVtkOutput(reader->readValue("VtkOutput", getVtkOutput()));
  setVtkFileName(reader->readString("VtkFileName", getVtkFileName()));
  setTxtOutput(reader->readValue("TxtOutput", getTxtOutput()));
  setTxtFileName(reader->readString("TxtFileName", getTxtFileName()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SteinerCompact::initialize()
{

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SteinerCompact::dataCheck()
{

  QVector<size_t> cDims(1, 1);
  QVector<DataArrayPath> dataArrayPaths;

  if (m_VtkOutput == true)
  {
    if (m_VtkFileName.isEmpty() == true)
    {
      QString ss = QObject::tr("The vtk output file must be set before executing this filter.");
      setErrorCondition(-1);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }

    // Make sure what we are checking is an actual file name and not a directory
    QFileInfo fi(m_VtkFileName);
    if (fi.isDir() == false)
    {
      QDir parentPath = fi.path();
      if (parentPath.exists() == false)
      {
        QString ss = QObject::tr("The directory path for the output file does not exist.");
        notifyWarningMessage(getHumanLabel(), ss, -1);
      }
    }
    else
    {
      QString ss = QObject::tr("The output file path is a path to an existing directory. Please change the path to point to a file");
      setErrorCondition(-1);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }
  }

  if (m_TxtOutput == true)
  {
    if (m_TxtFileName.isEmpty() == true)
    {
      QString ss = QObject::tr("The text output file must be set before executing this filter.");
      setErrorCondition(-1);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }

    // Make sure what we are checking is an actual file name and not a directory
    QFileInfo fi(m_TxtFileName);
    if (fi.isDir() == false)
    {
      QDir parentPath = fi.path();
      if (parentPath.exists() == false)
      {
        QString ss = QObject::tr("The directory path for the output file does not exist.");
        notifyWarningMessage(getHumanLabel(), ss, -1);
      }
    }
    else
    {
      QString ss = QObject::tr("The output file path is a path to an existing directory. Please change the path to point to a file");
      setErrorCondition(-1);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }
  }

  m_FeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter>(this, getFeatureIdsArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if (NULL != m_FeatureIdsPtr.lock().get()) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
  if (getErrorCondition() >= 0) { dataArrayPaths.push_back(getFeatureIdsArrayPath()); }

  m_CellPhasesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter>(this, getCellPhasesArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if (NULL != m_CellPhasesPtr.lock().get()) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  {
    m_CellPhases = m_CellPhasesPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
  if (getErrorCondition() >= 0) { dataArrayPaths.push_back(getCellPhasesArrayPath()); }

  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, dataArrayPaths);

  setErrorCondition(0);

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SteinerCompact::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true); // Set the fact that we are preflighting.
  emit preflightAboutToExecute(); // Emit this signal so that other widgets can do one file update
  emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck(); // Run our DataCheck to make sure everthing is setup correctly
  emit preflightExecuted(); // We are done preflighting this filter
  setInPreflight(false); // Inform the system this filter is NOT in preflight mode anymore.
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// rose of intersections (counts numbers of intersections of randomly placed lines at fixed directions with grain boundaries)
void SteinerCompact::rose_of_intersections(std::vector<std::vector<float>>& ROI)
{
  float line[3];
  float rectangle[4];
  float xintersections[2];
  float yintersections[2];
  float xdif, ydif;
  int64_t progressInt = 0;

  static float rndm;
  static float zero = 0.000001f;

  uint64_t directions = ROI[1].size();

  SIMPL_RANDOMNG_NEW();

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(m_FeatureIdsArrayPath.getDataContainerName());
  size_t udims[3] = { 0, 0, 0 };
  m->getGeometryAs<ImageGeom>()->getDimensions(udims);
  float res[3] = { 0.0f, 0.0f, 0.0f };
  m->getGeometryAs<ImageGeom>()->getResolution(res);

  int64_t dims[3] =
  {
    static_cast<int64_t>(udims[0]),
    static_cast<int64_t>(udims[1]),
    static_cast<int64_t>(udims[2]),
  };

  float maxdim = static_cast<float>(dims[0]) * res[0];
  if (static_cast<float>(dims[1]) * res[1] > maxdim) maxdim = static_cast<float>(dims[1]) * res[1];
  if (static_cast<float>(dims[2]) * res[2] > maxdim) maxdim = static_cast<float>(dims[2]) * res[2];
  float totlength = 1000 * maxdim;

  // rotation of arrays with dimensions and resolution for planes different from xy
  if (m_Plane == 1)	// xz plane
  {
    int rot = dims[1];
    dims[1] = dims[2];
    dims[2] = rot;
    rot = res[1];
    res[1] = res[2];
    res[2] = rot;
  }
  if (m_Plane == 2)	// yz plane
  {
    int rot = dims[0];
    dims[0] = dims[1];
    dims[1] = dims[2];
    dims[2] = rot;
    rot = res[0];
    res[0] = res[1];
    res[1] = res[2];
    res[2] = rot;
  }

  rectangle[0] = 0.0f;
  rectangle[1] = 0.0f;
  rectangle[2] = dims[0] * res[0];
  rectangle[3] = dims[1] * res[1];

  std::vector<float> numofintersections(ROI.size(), 0);
  std::vector<float> meanROI(ROI.size(), 0);

  float length = 0;
  int64_t z = 0;
  float alpha = 0;

  for (int32_t i = 0; i < directions; i++)
  {
    progressInt = (static_cast<float>(i) / static_cast<float>(directions)) * 100.0f;
    QString ss = QObject::tr("Evaluate random intersections || %1% Complete").arg(progressInt);
    notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);

    alpha = (static_cast<float>(i)+0.5f) * SIMPLib::Constants::k_Pi / static_cast<float>(directions);
    //alpha = (static_cast<float>(i)) * SIMPLib::Constants::k_Pi / directions;

    // the line is represented as line[0]*x+line[1]*y+line[2]=0;
    if (fabs(alpha - SIMPLib::Constants::k_Pi / 2) < zero)
    {
      line[0] = 1.0f;
      line[1] = 0.0f;
    }
    else
    {
      line[0] = -tan(alpha);
      line[1] = 1.0f;
    }

    for (int32_t phase = 1; phase < ROI.size(); phase++) ROI[phase][i] = 0;

    length = 0.0f;

    while (length < totlength)
    {
      float random_point_x = static_cast<float>(rg.genrand_res53()) * static_cast<float>(dims[0]) * res[0];
      float random_point_y = static_cast<float>(rg.genrand_res53()) * static_cast<float>(dims[1]) * res[1];
      line[2] = -random_point_x * line[0] - random_point_y * line[1];

      if (line_rectangle_intersections(xintersections, yintersections, rectangle, line) == 2)
      {
        xdif = xintersections[1] - xintersections[0];
        ydif = yintersections[1] - yintersections[0];
        rndm = static_cast<float>(rg.genrand_res53());
        z = static_cast<int64_t>(trunc(rndm * dims[2] * res[2]));
        for (int32_t phase = 1; phase < ROI.size(); phase++)
        {
          find_intersections(line, z, dims, res, numofintersections);
          ROI[phase][i] += numofintersections[phase];
        }
        length += sqrt(xdif * xdif + ydif * ydif);
      }
    }
    for (int32_t phase = 1; phase < ROI.size(); phase++)
    {
      ROI[phase][i] /= length;
      meanROI[phase] += ROI[phase][i] / static_cast<float>(directions);
    }
  }

  for (int32_t phase = 1; phase < ROI.size(); phase++)
  {
    if (fabs(meanROI[phase]) > zero)
    {
      for (int32_t i = 0; i < directions; i++)
      {
        ROI[phase][i] /= meanROI[phase];
      }
    }
  }
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// count the number of intersections of given line with grain boundaries
void SteinerCompact::find_intersections(float line[3], int64_t z, int64_t dims[3], float res[3], std::vector<float>& numofintersections)
{
  std::vector<int64_t> xvoxels;
  std::vector<int64_t> yvoxels;

  std::vector<float> xcoor0;
  std::vector<float> xcoor1;
  std::vector<float> ycoor0;
  std::vector<float> ycoor1;

  size_t numfeatures = m_FeatureIdsPtr.lock()->getNumberOfTuples();
  std::vector<uint64_t> bound;
  float rectangle[4];
  float xintersections[2];
  float yintersections[2];
  float limit1 = 0, limit2 = 0;
  int64_t low = 0, up = 0;
  static float zero = 0.000001f;

  if (fabs(line[1]) > zero)
  {
    for (int64_t x = 0; x < dims[0]; x++)
    {
      rectangle[0] = static_cast<float>(x)* res[0];
      rectangle[2] = static_cast<float>(x + 1) * res[0];
      limit1 = (line[0] * static_cast<float>(x)* res[0] + line[2]) / (-line[1] * res[1]);
      limit2 = (line[0] * static_cast<float>(x + 1)* res[0] + line[2]) / (-line[1] * res[1]);

      low = (int64_t)(trunc(fmin(limit1, limit2)));
      if (low < 0) low = 0;
      up = (int64_t)(ceil(fmax(limit1, limit2)));
      if (up > dims[1] - 1) up = dims[1] - 1;

      if (line[0] <= 0)
        for (int64_t y = low; y <= up; y++)
        {
          rectangle[1] = static_cast<float>(y)* res[1];
          rectangle[3] = static_cast<float>(y + 1) * res[1];
          if (line_rectangle_intersections(xintersections, yintersections, rectangle, line) == 2)
          {
            xvoxels.push_back(x);
            yvoxels.push_back(y);
            xcoor0.push_back(xintersections[0]);
            xcoor1.push_back(xintersections[1]);
            ycoor0.push_back(yintersections[0]);
            ycoor1.push_back(yintersections[1]);
          }
        }
      else if (line[0] > 0)
        for (int64_t y = up; y >= low; y--)
        {
          rectangle[1] = static_cast<float>(y)* res[1];
          rectangle[3] = static_cast<float>(y + 1) * res[1];
          if (line_rectangle_intersections(xintersections, yintersections, rectangle, line) == 2)
          {
            xvoxels.push_back(x);
            yvoxels.push_back(y);
            xcoor0.push_back(xintersections[0]);
            xcoor1.push_back(xintersections[1]);
            ycoor0.push_back(yintersections[0]);
            ycoor1.push_back(yintersections[1]);
          }
        }
    }
  }
  else
  {
    int64_t x = (uint64_t)round(-line[2] / (line[0] * res[0]));
    if (x < 0) x = 0;
    if (x > dims[0] - 1) x = dims[0] - 1;
    rectangle[0] = static_cast<float>(x)* res[0];
    rectangle[2] = static_cast<float>(x + 1) * res[0];
    for (int64_t y = 0; y < dims[1]; y++)
    {
      rectangle[1] = static_cast<float>(y)* res[1];
      rectangle[3] = static_cast<float>(y + 1) * res[1];
      if (line_rectangle_intersections(xintersections, yintersections, rectangle, line) == 2)
      {
        xvoxels.push_back(x);
        yvoxels.push_back(y);
        xcoor0.push_back(xintersections[0]);
        xcoor1.push_back(xintersections[1]);
        ycoor0.push_back(yintersections[0]);
        ycoor1.push_back(yintersections[1]);
      }
    }
  }

  // finally only some intersections are counted to the rose of intersections:
  // 1) intersection between two grains can be counted only once at each direction,
  //    this shall avoid multiple counting of points close to each other at irregular boundaries,
  //    this is utilized by the stack named interfaces
  // 2) grain boundaries are dilated (parameter smooth sets multiple of a voxel for the dilation)
  //    and intersections are identified with the dilated set
  //    this should diminish the influence of voxelation on the shape of Steiner compact
  //    (intersections occur more frequently for specific directions in voxelated data)

  std::vector<uint64_t> interfaces;
  bool penalization;
  float squaredlength;
  uint64_t g0 = 0, g1 = 0, g2 = 0;
  int64_t point0 = 0, point1 = 0, point2 = 0;
  int64_t phase0 = 0, phase1 = 0, phase2 = 0;

  float smooth = 0.25f;

  numofintersections.assign(numofintersections.size(), 0);

  for (int64_t i = 0; i < static_cast<int64_t>(xvoxels.size()) - 1; i++)
  {
    if (m_Plane == 0)
    {
      point0 = z * (dims[0] * dims[1]) + yvoxels[i] * dims[0] + xvoxels[i];
      point1 = z * (dims[0] * dims[1]) + yvoxels[i + 1] * dims[0] + xvoxels[i + 1];
    }
    else if (m_Plane == 1)  // rotation y <-> z
    {
      point0 = yvoxels[i] * (dims[0] * dims[2]) + z * dims[0] + xvoxels[i];
      point1 = yvoxels[i + 1] * (dims[0] * dims[2]) + z * dims[0] + xvoxels[i + 1];
    }
    else if (m_Plane == 2)  // rotation x <-> y and y <-> z
    {
      point0 = yvoxels[i] * (dims[2] * dims[0]) + xvoxels[i] * dims[2] + z;
      point1 = yvoxels[i + 1] * (dims[2] * dims[0]) + xvoxels[i + 1] * dims[2] + z;
    }
    g0 = m_FeatureIds[point0];
    g1 = m_FeatureIds[point1];
    phase0 = m_CellPhases[point0];
    phase1 = m_CellPhases[point1];
    if (g0 != g1 && (interfaces.size() == 0 || (std::find(interfaces.begin(), interfaces.end(), g0*numfeatures + g1) == interfaces.end())))
    {
      penalization = false;
      if (i < static_cast<int64_t>(xvoxels.size()) - 2)
      {
        if (m_Plane == 0)
        {
          point2 = z * (dims[0] * dims[1]) + yvoxels[i + 2] * dims[0] + xvoxels[i + 2];
        }
        else if (m_Plane == 1)  // rotation y <-> z
        {
          point2 = yvoxels[i + 2] * (dims[0] * dims[2]) + z * dims[0] + xvoxels[i + 2];
        }
        else if (m_Plane == 2)  // rotation x <-> y and y <-> z
        {
          point2 = yvoxels[i + 2] * (dims[2] * dims[0]) + xvoxels[i + 2] * dims[2] + z;
        }
        g2 = m_FeatureIds[point2];

        // intersection is not counted if the line just slightly hits a "corner pixel" (corner of L-shaped boundary)
        if (g0 == g2 && llabs(xvoxels[i] - xvoxels[i + 2]) == 1 && llabs(yvoxels[i] - yvoxels[i + 2]) == 1)
        {
          squaredlength = (xcoor0[i + 1] - xcoor1[i + 1]) * (xcoor0[i + 1] - xcoor1[i + 1]) + (ycoor0[i + 1] - ycoor1[i + 1]) * (ycoor0[i + 1] - ycoor1[i + 1]);
          if (squaredlength < (1 - smooth) * (1 - smooth) * (res[0] * res[0] + res[1] * res[1]))
          {
            penalization = true;
          }
        }
      }
      if (penalization == false)
      {
        interfaces.push_back(g0 * numfeatures + g1);
        numofintersections[phase0] += 1.0f;
        numofintersections[phase1] += 1.0f;
      }
    }
    else if (g0 == g1)  // boundary can be identified even here if the line is close enough to the neighbouring point2
    {
      if (fabs(xcoor0[i] - xcoor1[i] - 1) < zero)      // line hits the pixel horizontally
      {
        if (yvoxels[i] > 0 && (ycoor0[i] - yvoxels[i] * res[1] < smooth * res[1] || ycoor1[i] - yvoxels[i] * res[1] < smooth * res[1]))
        {
          if (m_Plane == 0)
          {
            point2 = z * (dims[0] * dims[1]) + (yvoxels[i] - 1) * dims[0] + xvoxels[i + 1];
          }
          else if (m_Plane == 1)  // rotation y <-> z
          {
            point2 = (yvoxels[i] - 1) * (dims[0] * dims[2]) + z * dims[0] + xvoxels[i + 1];
          }
          else if (m_Plane == 2)  // rotation x <-> y and y <-> z
          {
            point2 = (yvoxels[i] - 1) * (dims[2] * dims[0]) + xvoxels[i + 1] * dims[2] + z;

          }
          g2 = m_FeatureIds[point2];
          if (g0 != g2)
          {
            if (interfaces.size() == 0 || (std::find(interfaces.begin(), interfaces.end(), g0*numfeatures + g2) == interfaces.end()))
            {
              phase2 = m_CellPhases[point2];
              interfaces.push_back(g0 * numfeatures + g2);
              numofintersections[phase0] += 1.0f;
              numofintersections[phase2] += 1.0f;
            }
          }
        }
        if (yvoxels[i + 1] < dims[1] - 1 && ((yvoxels[i] + 1) * res[1] - ycoor0[i] < smooth * res[1] || (yvoxels[i] + 1) * res[1] - ycoor1[i] < smooth * res[1]))
        {
          if (m_Plane == 0)
          {
            point2 = z * (dims[0] * dims[1]) + (yvoxels[i] + 1) * dims[0] + xvoxels[i];
          }
          else if (m_Plane == 1)  // rotation y <-> z
          {
            point2 = (yvoxels[i] + 1) * (dims[0] * dims[2]) + z * dims[0] + xvoxels[i];
          }
          else if (m_Plane == 2)  // rotation x <-> y and y <-> z
          {
            point2 = (yvoxels[i] + 1) * (dims[2] * dims[0]) + xvoxels[i] * dims[2] + z;
          }
          g2 = m_FeatureIds[point2];
          if (g0 != g2)
          {
            if (interfaces.size() == 0 || (std::find(interfaces.begin(), interfaces.end(), g0*numfeatures + g2) == interfaces.end()))
            {
              phase2 = m_CellPhases[point2];
              interfaces.push_back(g0 * numfeatures + g2);
              numofintersections[phase0] += 1.0f;
              numofintersections[phase2] += 1.0f;
            }
          }
        }
      }
      else if (fabs(ycoor0[i] - ycoor1[i] - 1) < zero) // line hits the pixel vertically
      {
        if (xvoxels[i] > 0 && (xcoor0[i] - xvoxels[i] * res[0] < smooth * res[0] || xcoor1[i] - xvoxels[i] * res[0] < smooth * res[0]))
        {
          if (m_Plane == 0)
          {
            point2 = z * (dims[0] * dims[1]) + yvoxels[i] * dims[0] + xvoxels[i] - 1;
          }
          else if (m_Plane == 1)  // rotation y <-> z
          {
            point2 = yvoxels[i] * (dims[0] * dims[2]) + z * dims[0] + xvoxels[i] - 1;
          }
          else if (m_Plane == 2)  // rotation x <-> y and y <-> z
          {
            point2 = yvoxels[i] * (dims[2] * dims[0]) + (xvoxels[i] - 1) * dims[2] + z;
          }
          g2 = m_FeatureIds[point2];
          if (g0 != g2)
          {
            if (interfaces.size() == 0 || (std::find(interfaces.begin(), interfaces.end(), g0*numfeatures + g2) == interfaces.end()))
            {
              phase2 = m_CellPhases[point2];
              interfaces.push_back(g0 * numfeatures + g2);
              numofintersections[phase0] += 1.0f;
              numofintersections[phase2] += 1.0f;
            }
          }
        }
        if (xvoxels[i] < dims[0] - 1 && ((xvoxels[i] + 1) * res[0] - xcoor0[i] < smooth * res[0] || (xvoxels[i] + 1) * res[0] - xcoor1[i] < smooth * res[0]))
        {
          if (m_Plane == 0)
          {
            point2 = z * (dims[0] * dims[1]) + yvoxels[i] * dims[0] + xvoxels[i] + 1;
          }
          else if (m_Plane == 1)  // rotation y <-> z
          {
            point2 = yvoxels[i] * (dims[0] * dims[2]) + z * dims[0] + xvoxels[i] + 1;
          }
          else if (m_Plane == 2)  // rotation x <-> y and y <-> z
          {
            point2 = yvoxels[i] * (dims[2] * dims[0]) + (xvoxels[i] + 1) * dims[2] + z;
          }
          g2 = m_FeatureIds[point2];
          if (g0 != g2)
          {
            if (interfaces.size() == 0 || (std::find(interfaces.begin(), interfaces.end(), g0*numfeatures + g2) == interfaces.end()))
            {
              phase2 = m_CellPhases[point2];
              interfaces.push_back(g0 * numfeatures + g2);
              numofintersections[phase0] += 1.0f;
              numofintersections[phase2] += 1.0f;
            }
          }
        }
      }
    }
  }

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// find the intersections of the line (line[0]*x+line[1]*y+line[2]=0) with the rectangle (x1,y1,x2,y2)
uint64_t SteinerCompact::line_rectangle_intersections(float xintersections[2], float yintersections[2], float rectangle[4], float line[3])
{
  uint32_t intersnum = 0;
  static float zero = 0.000001f;

  if (fabs(line[1]) > zero)
  {
    float y = (line[0] * rectangle[0] + line[2]) / (-line[1]);
    if (y > rectangle[1] && y < rectangle[3])
    {
      xintersections[intersnum] = rectangle[0];
      yintersections[intersnum] = y;
      intersnum++;
    }
    y = (line[0] * rectangle[2] + line[2]) / (-line[1]);
    if (y > rectangle[1] && y < rectangle[3])
    {
      xintersections[intersnum] = rectangle[2];
      yintersections[intersnum] = y;
      intersnum++;
    }
  }
  if (fabs(line[0]) > zero)
  {
    float x = (line[1] * rectangle[1] + line[2]) / (-line[0]);
    if (intersnum < 2 && x > rectangle[0] && x < rectangle[2])
    {
      xintersections[intersnum] = x;
      yintersections[intersnum] = rectangle[1];
      intersnum++;
    }
    x = (line[1] * rectangle[3] + line[2]) / (-line[0]);
    if (intersnum < 2 && x > rectangle[0] && x < rectangle[2])
    {
      xintersections[intersnum] = x;
      yintersections[intersnum] = rectangle[3];
      intersnum++;
    }
  }

  return intersnum;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// find coordinates of vertices of given site of the steiner compact
// input: ROI (rose of intersections), index of the site
// output: vertices = (x1, y1, x2, y2) of the site, length of the site
void SteinerCompact::find_one_site_vertices(std::vector<float> ROI, int64_t index, float vertices[4], float& length)
{
  float cmin = std::numeric_limits<float>::max();
  float cmax = -std::numeric_limits<float>::max();
  float ang;
  float radius;
  float coor;
  int64_t size = static_cast<int64_t>(ROI.size());

  std::vector<float> ROI2(size);
  ROI2[0] = ROI[index];

  // ROI is reordered periodically such that index is at zero position
  for (int64_t i = 0; i < size; i++)
  {
    if (index + i < size) ROI2[i] = ROI[index + i];
    else ROI2[i] = ROI[index + i - size];
  }

  // at first, find coordinates for the reordered array which corresponds to a vertical site
  // its x-coordinate is ROI2[0], y-coordinates are found by the method from Benes and Gokhale (2000)
  for (int64_t i = 0; i < size; i++)
  {
    if (i > 0)
    {
      ang = -(static_cast<float>(i)) * SIMPLib::Constants::k_Pi / static_cast<float>(size);
      coor = (ROI2[0] * cos(ang) - ROI2[i]) / sin(ang);
      if (coor < cmin) cmin = coor;
    }
    if (i < size - 1)
    {
      ang = (static_cast<float>(i + 1)) * SIMPLib::Constants::k_Pi / static_cast<float>(size);
      coor = (ROI2[0] * cos(ang) - ROI2[size - 1 - i]) / sin(ang);
      if (coor > cmax) cmax = coor;
    }
  }

  length = cmin - cmax;
  if (length < 0) length = 0;

  // now convert to polar coordinates and rotate the vertices
  radius = sqrt(ROI2[0] * ROI2[0] + cmin * cmin);
  ang = atan2(cmin, ROI2[0]) + 0.5f * SIMPLib::Constants::k_Pi / static_cast<float>(size);
  vertices[0] = radius * cos(ang + SIMPLib::Constants::k_Pi * (static_cast<float>(index) / static_cast<float>(size) + 0.5f));
  vertices[1] = radius * sin(ang + SIMPLib::Constants::k_Pi * (static_cast<float>(index) / static_cast<float>(size) + 0.5f));

  radius = sqrt(ROI2[0] * ROI2[0] + cmax * cmax);
  ang = atan2(cmax, ROI2[0]) + 0.5f * SIMPLib::Constants::k_Pi / static_cast<float>(size);
  vertices[2] = radius * cos(ang + SIMPLib::Constants::k_Pi * (static_cast<float>(index) / static_cast<float>(size) + 0.5f));
  vertices[3] = radius * sin(ang + SIMPLib::Constants::k_Pi * (static_cast<float>(index) / static_cast<float>(size) + 0.5f));
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// find coordinates of all vertices of the Steiner compact
void SteinerCompact::find_all_vertices(std::vector<std::vector<float>>& vertices_x, std::vector<std::vector<float>>& vertices_y, std::vector<std::vector<float>>& radii, std::vector<std::vector<float>>& ROI)
{
  float vertices[4];
  float length = 0.0f;
  int64_t numphases = ROI.size() - 1;
  int64_t numdirections = ROI[1].size();

  vertices_x.resize(ROI.size());
  vertices_y.resize(ROI.size());
  radii.resize(ROI.size());

  // find coordinates of vertices of the Steiner compact from the rose of intersections
  for (int64_t phase = 1; phase <= numphases; phase++)
  {
    for (int64_t direction = 0; direction < numdirections; direction++)
    {
      find_one_site_vertices(ROI[phase], direction, vertices, length);
      if (length > 0)		// use only sites with positive length
      {
        vertices_x[phase].push_back(vertices[0]);
        vertices_y[phase].push_back(vertices[1]);
        radii[phase].push_back(ROI[phase][direction]);
      }
    }
  }
}

void SteinerCompact::output_vtk(std::vector<std::vector<float>>& vertices_x, std::vector<std::vector<float>>& vertices_y, std::vector<std::vector<float>>& radii, std::vector<std::vector<float>>& ROI)
{
  std::ofstream pom("pom.txt");
  pom << m_VtkOutput << std::endl;

  FILE *vtk = NULL;

  // Make sure any directory path is also available as the user may have just typed
  // in a path without actually creating the full path
  QFileInfo fi(m_VtkFileName);
  QString parentPath = fi.path();
  QDir dir;
  if (!dir.mkpath(parentPath))
  {
    QString ss = QObject::tr("Error creating parent path '%1'").arg(parentPath);
    setErrorCondition(-2031000);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    return;
  }

  vtk = fopen(getVtkFileName().toLatin1().data(), "w");
  if (NULL == vtk)
  {
    QString ss = QObject::tr("Error opening output vtk file '%1'\n ").arg(m_VtkFileName);
    setErrorCondition(-2031001);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    return;
  }
  ScopedFileMonitor fMon(vtk);

  size_t numvertices = vertices_x[1].size();
  size_t numphases = ROI.size() - 1;
  size_t numdirections = ROI[1].size();

  // header
  fprintf(vtk, "# vtk DataFile Version 2.0\n");
  fprintf(vtk, "Steiner compact\n");
  fprintf(vtk, "ASCII\n");
  fprintf(vtk, "DATASET POLYDATA\n");

  fprintf(vtk, "\nPOINTS %lld float\n", static_cast<long long int>(2 * numvertices + numphases + numphases * 2 * numdirections) );
  for (size_t phase = 1; phase <= numphases; phase++)
  {
    if (m_Plane == 0)
    {
      for (uint64_t i = 0; i < vertices_x[phase].size(); i++)
      {
        fprintf(vtk, "%f %f %f\n", vertices_x[phase][i], vertices_y[phase][i], static_cast<float>(phase));
      }
      for (uint64_t i = 0; i < vertices_x[phase].size(); i++)
      {
        fprintf(vtk, "%f %f %f\n", -vertices_x[phase][i], -vertices_y[phase][i], static_cast<float>(phase));
      }
      fprintf(vtk, "%f %f %f\n", 0.0f, 0.0f, static_cast<float>(phase));  // origin
    }
    else if (m_Plane == 1)
    {
      for (uint64_t i = 0; i < vertices_x[phase].size(); i++)
      {
        fprintf(vtk, "%f %f %f\n", vertices_x[phase][i], static_cast<float>(phase), vertices_y[phase][i]);
      }
      for (uint64_t i = 0; i < vertices_x[phase].size(); i++)
      {
        fprintf(vtk, "%f %f %f\n", -vertices_x[phase][i], static_cast<float>(phase), -vertices_y[phase][i]);
      }
      fprintf(vtk, "%f %f %f\n", 0.0f, static_cast<float>(phase), 0.0f);  // origin
    }
    else if (m_Plane == 2)
    {
      for (uint64_t i = 0; i < vertices_x[phase].size(); i++)
      {
        fprintf(vtk, "%f %f %f\n", static_cast<float>(phase), vertices_x[phase][i], vertices_y[phase][i]);
      }
      for (uint64_t i = 0; i < vertices_x[phase].size(); i++)
      {
        fprintf(vtk, "%f %f %f\n", static_cast<float>(phase), -vertices_x[phase][i], -vertices_y[phase][i]);
      }
      fprintf(vtk, "%f %f %f\n", static_cast<float>(phase), 0.0f, 0.0f);  // origin
    }
  }

  // points for reference Steiner compact (regular polygon)
  float angle = SIMPLib::Constants::k_Pif / static_cast<float>(numdirections);
  float r = 1.0f / cosf(0.5f * angle);
  float s, c, p;
  for (size_t phase = 1; phase <= numphases; phase++)
  {
    p = static_cast<float>(phase);
    for (uint64_t i = 0; i < 2 * numdirections; i++)
    {
      s = sinf(static_cast<float>(i)* angle);
      c = cosf(static_cast<float>(i)* angle);
      if (m_Plane == 0) fprintf(vtk, "%f %f %f\n", r * c, r * s, p);
      else if (m_Plane == 1) fprintf(vtk, "%f %f %f\n", r * c, p, r * s);
      else if (m_Plane == 2) fprintf(vtk, "%f %f %f\n", p, r * c, r * s);
    }
  }

  // lines of reference Steiner compact(regular polygon)
  int64_t curindex = 0;
  fprintf(vtk, "\nLINES %lld %lld\n", static_cast<long long int>(numphases), static_cast<long long int>(numphases * (2 * numdirections + 2)) );
  for (size_t phase = 1; phase <= numphases; phase++)
  {
    fprintf(vtk, "%lld ", static_cast<long long int>(2 * numdirections + 1));
    for (size_t i = 0; i < 2 * numdirections; i++)
    {
      fprintf(vtk, "%lld ", static_cast<long long int>(2 * numvertices + numphases + curindex + i));
    }
    fprintf(vtk, "%lld\n", static_cast<long long int>(2 * numvertices + numphases + curindex));           // return to 0
    curindex += 2 * numdirections;
  }

  // Steiner compact for each phase
  fprintf(vtk, "\nPOLYGONS %lld %lld\n", static_cast<long long int>(2 * numvertices), static_cast<long long int>(2 * 4 * numvertices));
  curindex = 0;
  size_t origin = 2 * numvertices + numphases - 1;
  for (size_t phase = 1; phase <= numphases; phase++)
  {
    for (size_t i = 0; i < 2 * vertices_x[phase].size() - 1; i++)
    {
      fprintf(vtk, "%d %lld %lld %lld\n", 3, static_cast<long long int>(curindex + i), static_cast<long long int>(curindex + i + 1), static_cast<long long int>(origin));
    }
    fprintf(vtk, "%d %lld %lld %lld\n", 3, static_cast<long long int>(curindex), static_cast<long long int>(curindex + 2 * vertices_x[phase].size() - 1), static_cast<long long int>(origin));  // last triangle
    curindex += 2 * vertices_x[phase].size();
  }

  fprintf(vtk, "\nCELL_DATA %lld\n", static_cast<long long int>(2 * numvertices + numphases));
  fprintf(vtk, "SCALARS data float\n");
  fprintf(vtk, "LOOKUP_TABLE default\n");

  // data for reference Steiner compact (regular polygon)
  for (size_t phase = 1; phase <= numphases; phase++)
  {
    fprintf(vtk, "%f\n", 1.0f);
  }

  // data for Steiner compact for each phase
  for (size_t phase = 1; phase <= numphases; phase++)
  {
    for (int64_t sign = 0; sign <= 1; sign++)
    {
      for (uint64_t i = 0; i < radii[phase].size(); i++)
      {
        fprintf(vtk, "%f ", radii[phase][i]);
      }
    }
    fprintf(vtk, "\n");
  }


}

void SteinerCompact::output_txt(std::vector<std::vector<float>>& vertices_x, std::vector<std::vector<float>>& vertices_y, std::vector<std::vector<float>>& ROI)
{
  FILE *txt = NULL;

  // Make sure any directory path is also available as the user may have just typed
  // in a path without actually creating the full path
  QFileInfo fi(m_TxtFileName);
  QString parentPath = fi.path();
  QDir dir;
  if (!dir.mkpath(parentPath))
  {
    QString ss = QObject::tr("Error creating parent path '%1'").arg(parentPath);
    setErrorCondition(-2031000);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    return;
  }

  txt = fopen(getTxtFileName().toLatin1().data(), "w");
  if (NULL == txt)
  {
    QString ss = QObject::tr("Error opening output txt file '%1'\n ").arg(m_TxtFileName);
    setErrorCondition(-2031001);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    return;
  }
  ScopedFileMonitor fMon(txt);

  size_t numphases = ROI.size() - 1;
  size_t numdirections = ROI[1].size();

  fprintf(txt, "Distances_of_edges_from_origin\n");
  fprintf(txt, "Phase   Angle   Distance\n");
  float stepangle = SIMPLib::Constants::k_Pif / numdirections;
  float angle = 0;
  for (uint64_t phase = 1; phase <= numphases; phase++)
  {
    for (uint64_t symmetry = 0; symmetry <= 1; symmetry++)
    {
      for (uint64_t site = 0; site < numdirections; site++)
      {
        angle = (static_cast<float>(site + symmetry * numdirections) + 0.5f) * stepangle + 0.5f * SIMPLib::Constants::k_Pif;
        angle *= 180.0f / SIMPLib::Constants::k_Pi;
        if (angle >= 360.0f) angle -= 360.0f;
        fprintf(txt, "%lld   %f   %f\n", static_cast<long long int>(phase), angle, ROI[phase][site]);
      }
    }
  }
  fprintf(txt, "\n");

  // Steiner compact for each phase
  fprintf(txt, "Coordinates_of_vertices\n");
  fprintf(txt, "Phase   Vertex   X   Y\n");
  for (size_t phase = 1; phase <= numphases; phase++)
  {
    for (size_t i = 0; i < vertices_x[phase].size(); i++)
    {
      fprintf(txt, "%lld   %lld   %f   %f\n", static_cast<long long int>(phase), static_cast<long long int>(i), vertices_x[phase][i], vertices_y[phase][i]);
    }
    for (uint64_t i = 0; i < vertices_x[phase].size(); i++)
    {
      fprintf(txt, "%lld   %lld   %f   %f\n", static_cast<long long int>(phase), static_cast<long long int>(vertices_x[phase].size() + i), -vertices_x[phase][i], -vertices_y[phase][i]);
    }
  }

}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void SteinerCompact::execute()
{
  setErrorCondition(0);
  dataCheck();
  if(getErrorCondition() < 0) { return; }

  if (getCancel() == true) { return; }

  if (getErrorCondition() < 0)
  {
    QString ss = QObject::tr("Some error message");
    setErrorCondition(-99999999);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    return;
  }

  int32_t maxPhase = 0;
  size_t totalPoints = m_FeatureIdsPtr.lock()->getNumberOfTuples();
  for (int i = 0; i < totalPoints; i++)
  {
    if (m_CellPhases[i] > maxPhase) maxPhase = m_CellPhases[i];
  }

  uint64_t directions;

  switch (m_Sites)
  {
  case 0: directions = 4; break;
  case 1: directions = 6; break;
  case 2: directions = 8; break;
  case 3: directions = 12; break;
  case 4: directions = 18; break;
  }

  std::vector<std::vector<float>> ROI(maxPhase + 1);
  for (int32_t i = 1; i <= maxPhase; i++)
  {
    ROI[i].resize(directions);
  }

  // find number of intersections per unit length of test lines in all directions
  rose_of_intersections(ROI);

  // write the Steiner compact to vtk or txt file
  if (m_VtkOutput == true || m_TxtOutput == true)
  {
    QString ss = QObject::tr("Write Output File(s)");
    notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);
    std::vector<std::vector<float>> draw_vertices_x;
    std::vector<std::vector<float>> draw_vertices_y;
    std::vector<std::vector<float>> radii;
    find_all_vertices(draw_vertices_x, draw_vertices_y, radii, ROI);
    if (m_VtkOutput == true) output_vtk(draw_vertices_x, draw_vertices_y, radii, ROI);
    if (m_TxtOutput == true) output_txt(draw_vertices_x, draw_vertices_y, ROI);
  }

  notifyStatusMessage(getHumanLabel(), "Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer SteinerCompact::newFilterInstance(bool copyFilterParameters)
{
  SteinerCompact::Pointer filter = SteinerCompact::New();
  if(true == copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString SteinerCompact::getCompiledLibraryName()
{
  return AnisotropyConstants::AnisotropyBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString SteinerCompact::getBrandingString()
{
  return "Anisotropy";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString SteinerCompact::getFilterVersion()
{
  QString version;
  QTextStream vStream(&version);
  vStream << SIMPLib::Version::Major() << "." << SIMPLib::Version::Minor() << "." << SIMPLib::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString SteinerCompact::getGroupName()
{
  return SIMPL::FilterGroups::ReconstructionFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString SteinerCompact::getSubGroupName()
{
  return AnisotropyConstants::FilterSubGroups::AnisotropicAlignment;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString SteinerCompact::getHumanLabel()
{
  return "Steiner Compact";
}

