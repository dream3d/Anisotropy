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

#include "AdaptiveAlignment.h"

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/OutputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/SIMPLibVersion.h"
#include "SIMPLib/Utilities/FileSystemPathHelper.h"
#include "SIMPLib/FilterParameters/MultiDataArraySelectionFilterParameter.h"


#include "Anisotropy/AnisotropyConstants.h"
#include "Anisotropy/AnisotropyFilters/ItkBridge.h"

#include "itkHoughTransform2DCirclesImageFilter.h"
#include "itkHoughTransform2DLinesImageFilter.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AdaptiveAlignment::AdaptiveAlignment()
: m_DataContainerName(SIMPL::Defaults::ImageDataContainerName)
, m_CellAttributeMatrixName(SIMPL::Defaults::CellAttributeMatrixName)
, m_WriteAlignmentShifts(false)
, m_AlignmentShiftFileName("")
, m_GlobalCorrection(0)
, m_InputPath("")
, m_ShiftX(0.0f)
, m_ShiftY(0.0f)
, m_ImageDataArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::ImageData)
, m_NewCellArrayName("")
, m_MinRadius(0.0f)
, m_MaxRadius(0.0f)
, m_NumberCircles(0)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AdaptiveAlignment::~AdaptiveAlignment() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignment::setupFilterParameters()
{
  FilterParameterVector parameters;
  QStringList linkedProps("AlignmentShiftFileName");

//  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Write Alignment Shift File", WriteAlignmentShifts, FilterParameter::Parameter, AdaptiveAlignment, linkedProps));
//  parameters.push_back(SIMPL_NEW_OUTPUT_FILE_FP("Alignment File", AlignmentShiftFileName, FilterParameter::Parameter, AdaptiveAlignment, "", "*.txt"));

  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Global Correction");
    parameter->setPropertyName("GlobalCorrection");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(AdaptiveAlignment, this, GlobalCorrection));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(AdaptiveAlignment, this, GlobalCorrection));

    parameter->setDefaultValue(0);

    QVector<QString> choices;
    choices.push_back("None");
    choices.push_back("SEM Images");
    choices.push_back("Own Shifts");

    parameter->setChoices(choices);
    QStringList linkedProps;
    linkedProps << "ImageDataArrayPath"
                << "ShiftX"
                << "ShiftY";
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Parameter);
    parameters.push_back(parameter);

    parameters.push_back(SeparatorFilterParameter::New("Image Data", FilterParameter::RequiredArray));
    {
      DataArraySelectionFilterParameter::RequirementType req =
          DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::UInt8, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
      QVector<QVector<size_t>> cDims;
      cDims.push_back(QVector<size_t>(1, 1));
      cDims.push_back(QVector<size_t>(1, 3));
      cDims.push_back(QVector<size_t>(1, 4));
      req.componentDimensions = cDims;
      parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Image Data", ImageDataArrayPath, FilterParameter::RequiredArray, AdaptiveAlignment, req, 1));
    }

    parameters.push_back(SIMPL_NEW_FLOAT_FP("Total Shift In X-Direction (Microns)", ShiftX, FilterParameter::Parameter, AdaptiveAlignment, 2));
    parameters.push_back(SIMPL_NEW_FLOAT_FP("Total Shift In Y-Direction (Microns)", ShiftY, FilterParameter::Parameter, AdaptiveAlignment, 2));
  }
  {
    MultiDataArraySelectionFilterParameter::RequirementType req;
    parameters.push_back(SIMPL_NEW_MDA_SELECTION_FP("Attribute Arrays to Ignore", IgnoredDataArrayPaths, FilterParameter::Parameter, AdaptiveAlignment, req));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignment::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setAlignmentShiftFileName(reader->readString("AlignmentShiftFileName", getAlignmentShiftFileName()));
  setWriteAlignmentShifts(reader->readValue("WriteAlignmentShifts", getWriteAlignmentShifts()));
  setGlobalCorrection(reader->readValue("UseImages", getGlobalCorrection()));
  setImageDataArrayPath(reader->readDataArrayPath("ImageDataArrayPath", getImageDataArrayPath()));
  setShiftX(reader->readValue("ShiftX", getShiftX()));
  setShiftY(reader->readValue("ShiftY", getShiftY()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignment::initialize()
{
  m_RectangleCorners.clear();
  m_CalibratingCircles.clear();
  m_InterfaceEdges.clear();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignment::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);
  DataArrayPath tempPath;

  ImageGeom::Pointer image = getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom, AbstractFilter>(this, getDataContainerName());
  if(getErrorCondition() < 0)
  {
    return;
  }

  if(image->getXPoints() <= 1 || image->getYPoints() <= 1 || image->getZPoints() <= 1)
  {
    QString ss =
        QObject::tr("The Image Geometry is not 3D and cannot be run through this filter. The dimensions are (%1,%2,%3)").arg(image->getXPoints()).arg(image->getYPoints()).arg(image->getZPoints());
    setErrorCondition(-3010);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  tempPath.update(getDataContainerName(), getCellAttributeMatrixName(), "");
  getDataContainerArray()->getPrereqAttributeMatrixFromPath<AbstractFilter>(this, tempPath, -301);

  if(getWriteAlignmentShifts())
  {
    FileSystemPathHelper::CheckOutputFile(this, "Alignment Shift File Name", getAlignmentShiftFileName(), true);
  }

  if(m_GlobalCorrection == 1)
  {
    int32_t numImageComp = 1;
    QVector<DataArrayPath> imageDataArrayPaths;
    IDataArray::Pointer iDataArray = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getImageDataArrayPath());
    if(getErrorCondition() < 0)
    {
      return;
    }
    if(nullptr != iDataArray.get())
    {
      numImageComp = iDataArray->getNumberOfComponents();
    }
    QVector<size_t> cDims(1, numImageComp);
    m_ImageDataPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<uint8_t>, AbstractFilter>(this, getImageDataArrayPath(),
                                                                                                         cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
    if(nullptr != m_ImageDataPtr.lock())                                                                         /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
    {
      m_ImageData = m_ImageDataPtr.lock()->getPointer(0);
    } /* Now assign the raw pointer to data from the DataArray<T> object */
    if(getErrorCondition() >= 0)
    {
      imageDataArrayPaths.push_back(getImageDataArrayPath());
    }

    getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, imageDataArrayPaths);

    size_t udims1[3] = {0, 0, 0};
    DataContainer::Pointer m1 = getDataContainerArray()->getDataContainer(getDataContainerName());
    std::tie(udims1[0], udims1[1], udims1[2]) = m1->getGeometryAs<ImageGeom>()->getDimensions();

    size_t udims2[3] = {0, 0, 0};
    DataContainer::Pointer m2 = getDataContainerArray()->getDataContainer(getImageDataArrayPath());
    std::tie(udims2[0], udims2[1], udims2[2]) = m2->getGeometryAs<ImageGeom>()->getDimensions();

    if(static_cast<uint64_t>(udims1[2]) != static_cast<uint64_t>(udims2[2]))
    {
      QString ss = QObject::tr("Image Data and Cell Data must have the same Z-dimension.");
      setErrorCondition(-3012);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignment::preflight()
{
  setInPreflight(true);
  emit preflightAboutToExecute();
  emit updateFilterParameters(this);
  dataCheck();
  emit preflightExecuted();
  setInPreflight(false);
}

void AdaptiveAlignment::create_array_for_flattened_image()
{
  QVector<size_t> cDims(1, 1);
  DataArrayPath tempPath;
  tempPath.update(getImageDataArrayPath().getDataContainerName(), getImageDataArrayPath().getAttributeMatrixName(), "tempFlatImageDataName");
  m_FlatImageDataPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<AnisotropyConstants::DefaultPixelType>, AbstractFilter, AnisotropyConstants::DefaultPixelType>(
      this, tempPath, 0, cDims);                 /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_FlatImageDataPtr.lock())       /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_FlatImageData = m_FlatImageDataPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

void AdaptiveAlignment::delete_array_for_flattened_image()
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getImageDataArrayPath());
  AttributeMatrix::Pointer attrMat = m->getAttributeMatrix(getImageDataArrayPath().getAttributeMatrixName());
  attrMat->removeAttributeArray("tempFlatImageDataName");
}

void AdaptiveAlignment::flatten_image()
{
  // flatten image if the number of components is > 1 (RGB or ARGB image)
  int32_t comp = m_ImageDataPtr.lock()->getNumberOfComponents();
  size_t totalPoints = m_ImageDataPtr.lock()->getNumberOfTuples();

  notifyStatusMessage(getHumanLabel(), "Flatten image");
  if(comp > 1)
  {
    float Rfactor = 0.21f;
    float Gfactor = 0.72f;
    float Bfactor = 0.07f;

    for(uint32_t i = 0; i < totalPoints; i++)
    {
      m_FlatImageData[i] = int32_t((m_ImageData[comp * i] * Rfactor) + (m_ImageData[comp * i + 1] * Gfactor) + (m_ImageData[comp * i + 2] * Bfactor));
    }
  }
  else if(comp == 1)
  {
    for(uint32_t i = 0; i < totalPoints; i++)
    {
      m_FlatImageData[i] = m_ImageData[i];
    }
  }
}

bool AdaptiveAlignment::find_calibrating_circles()
{
  QString ss = "";
  bool found = true;

  notifyStatusMessage(getHumanLabel(), "Find circles");
  // DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getImageDataContainerName());
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getImageDataArrayPath());
  QString attrMatName = getImageDataArrayPath().getAttributeMatrixName();

  m_NumberCircles = 2;
  m_MinRadius = 10;
  m_MaxRadius = 60;

  size_t udims[3] = {0, 0, 0};
  std::tie(udims[0], udims[1], udims[2]) = m->getGeometryAs<ImageGeom>()->getDimensions();

  uint64_t dims[3] = {
      static_cast<uint64_t>(udims[0]), static_cast<uint64_t>(udims[1]), static_cast<uint64_t>(udims[2]),
  };

  // wrap raw and processed image data as itk::images
  AnisotropyConstants::DefaultImageType::Pointer inputImage = ITKUtilitiesType::CreateItkWrapperForDataPointer(m, attrMatName, m_FlatImageData);

  // AnisotropyConstants::DefaultImageType::Pointer outputImage = ITKUtilitiesType::CreateItkWrapperForDataPointer(m, attrMatName, m_NewCellArray); // delete this
  // AnisotropyConstants::DefaultSliceType::IndexType localIndex; // delete this
#if defined(ITK_VERSION_MAJOR) && ITK_VERSION_MAJOR == 4
  using HoughTransformFilterType = itk::HoughTransform2DCirclesImageFilter<AnisotropyConstants::DefaultPixelType, AnisotropyConstants::FloatPixelType>;
#else
  using HoughTransformFilterType = itk::HoughTransform2DCirclesImageFilter<AnisotropyConstants::DefaultPixelType, AnisotropyConstants::FloatPixelType, AnisotropyConstants::FloatPixelType>;
#endif
  HoughTransformFilterType::Pointer houghFilter = HoughTransformFilterType::New();
  houghFilter->SetNumberOfCircles(m_NumberCircles);
  houghFilter->SetMinimumRadius(m_MinRadius);
  houghFilter->SetMaximumRadius(m_MaxRadius);
  /*optional parameters, these are the default values
  houghFilter->SetSweepAngle( 0 );
  houghFilter->SetSigmaGradient( 1 );
  houghFilter->SetVariance( 5 );
  houghFilter->SetDiscRadiusRatio( 10 );
  */

  std::vector<uint64_t> icenter(4);
  std::vector<float> fcenter(4);
  std::vector<float> radius(2);
  // loop over slices
  for(int i = 0; i < dims[2]; ++i)
  {

    AnisotropyConstants::DefaultSliceType::Pointer inputSlice = ITKUtilitiesType::ExtractSlice(inputImage, AnisotropyConstants::ZSlice, i);

    // Hough transform is done for the first slice only
    // to roughly identify the area of calibrating circles
    if(i == 0)
    {
      // extract slice and transform
      ss = QObject::tr("Finding Calibrating Circles");
      notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);

      // input slice here
      houghFilter->SetInput(inputSlice);
      houghFilter->Update();
      AnisotropyConstants::FloatSliceType::Pointer localAccumulator = houghFilter->GetOutput();

      // find circles
      HoughTransformFilterType::CirclesListType circles = houghFilter->GetCircles();
      HoughTransformFilterType::CirclesListType::const_iterator itCircles;

      /// circle variables here
      itCircles = circles.begin();
      icenter[0] = static_cast<uint64_t>((*itCircles)->GetObjectToParentTransform()->GetOffset().GetElement(0));
      icenter[1] = static_cast<uint64_t>((*itCircles)->GetObjectToParentTransform()->GetOffset().GetElement(1));
#if defined(ITK_VERSION_MAJOR) && ITK_VERSION_MAJOR == 4
      radius[0] = static_cast<float>((*itCircles)->GetRadius()[0]);
#else
      radius[0] = static_cast<float>((*itCircles)->GetRadiusInObjectSpace()[0]);
#endif
      itCircles++;
      icenter[2] = static_cast<uint64_t>((*itCircles)->GetObjectToParentTransform()->GetOffset().GetElement(0));
      icenter[3] = static_cast<uint64_t>((*itCircles)->GetObjectToParentTransform()->GetOffset().GetElement(1));
#if defined(ITK_VERSION_MAJOR) && ITK_VERSION_MAJOR == 4
      radius[1] = static_cast<float>((*itCircles)->GetRadius()[0]);
#else
      radius[1] = static_cast<float>((*itCircles)->GetRadiusInObjectSpace()[0]);
#endif
      // check that the circles were found correctly:
      // midpoint of the centres is not far from the centre of the image
      if(std::fabs(0.5f * static_cast<float>(icenter[0] + icenter[2]) - 0.5f * static_cast<float>(dims[0])) > 100.0f)
       { found = false; }
      // circles are in lower part of the image
      if(icenter[1] < dims[1] / 2 || icenter[3] < dims[1] / 2)
       { found = false; }
      // y-coordinates of the circles do not differ much
       if((icenter[1] - icenter[3]) > 20)
       {found = false;}
    }

    if(found)
    {
      // exponential smoothing in encapsulating square for finer (float) center
      // this is done for each slice to identify the circles
      ss = QObject::tr("Exponential smoothing");
      notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);

      float weighting_factor = -0.01f;
      float range_factor = 1.5f;
      float weight, totalweight;
      uint64_t xmin, xmax, ymin, ymax;
      for(uint8_t circle = 0; circle <= 1; circle++)
      {
        xmin = icenter[2 * circle] - range_factor * radius[circle];
        xmax = icenter[2 * circle] + range_factor * radius[circle];
        ymin = icenter[2 * circle + 1] - range_factor * radius[circle];
        ymax = icenter[2 * circle + 1] + range_factor * radius[circle];
        if(int64_t(xmin) < 0)
          {xmin = 0;}
        if(xmax > dims[0] - 1)
          {xmax = dims[0] - 1;}
        if(int64_t(ymin) < 0)
          {ymin = 0;}
        if(ymax > dims[1] - 1)
          {ymax = dims[1] - 1;}

        totalweight = 0.0f;
        for(uint64_t x = xmin; x <= xmax; x++)
        {
          for(uint64_t y = ymin; y <= ymax; y++)
          {
            weight = exp(weighting_factor * m_FlatImageData[i * dims[0] * dims[1] + y * dims[0] + x]);
            fcenter[2 * circle] += static_cast<float>(x) * weight;
            fcenter[2 * circle + 1] += static_cast<float>(y) * weight;
            totalweight += weight;
          }
        }
        if(totalweight > 0)
        {
          fcenter[2 * circle] /= totalweight;
          fcenter[2 * circle + 1] /= totalweight;
        }
      }

      // average coordinates of the two circles (midpoint of the segment connecting the two centres)
      m_CalibratingCircles[i][0] = 0.5f * (fcenter[0] + fcenter[2]);
      m_CalibratingCircles[i][1] = 0.5f * (fcenter[1] + fcenter[3]);

      // outFile << m_CalibratingCircles[i][0] << " " << m_CalibratingCircles[i][1] << std::endl;
    }
  }

  return found;
}

bool AdaptiveAlignment::find_rectangles()
{
  bool found = true;

  QString ss = QObject::tr("Finding The Mapped Areas");
  notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getImageDataArrayPath());

  size_t udims[3] = {0, 0, 0};
  std::tie(udims[0], udims[1], udims[2]) = m->getGeometryAs<ImageGeom>()->getDimensions();

  uint64_t dims[3] = {
      static_cast<uint64_t>(udims[0]), static_cast<uint64_t>(udims[1]), static_cast<uint64_t>(udims[2]),
  };

  uint64_t comp = m_ImageDataPtr.lock()->getNumberOfComponents();
  uint64_t index;

  for(uint64_t i = 0; i < dims[2]; i++)
  {
    m_RectangleCorners[i][0] = dims[0] - 1;
    m_RectangleCorners[i][1] = dims[1] - 1;
    m_RectangleCorners[i][2] = 0;
    m_RectangleCorners[i][3] = 0;
    for(uint64_t j = 0; j < dims[1]; j++)
    {
      for(uint64_t k = 0; k < dims[0]; k++)
      {
        index = (i * dims[0] * dims[1]) + (j * dims[0]) + k;
        if(m_ImageData[comp * index] == 0 && m_ImageData[comp * index + 1] == 255 && m_ImageData[comp * index + 2] == 0) // green pixel
        {
          if(k < m_RectangleCorners[i][0])
            {m_RectangleCorners[i][0] = k;}// x-coordinate of upper left corner
          if(j < m_RectangleCorners[i][1])
            {m_RectangleCorners[i][1] = j;} // y-coordinate of upper left corner
          if(k > m_RectangleCorners[i][2])
            {m_RectangleCorners[i][2] = k;} // x-coordinate of bottom right corner
          if(j > m_RectangleCorners[i][3])
            {m_RectangleCorners[i][3] = j;} // y-coordinate of bottom right corner
        }
      }
    }
    if(m_RectangleCorners[i][0] == dims[0] - 1 || m_RectangleCorners[i][1] == dims[1] - 1 || m_RectangleCorners[i][2] == 0 || m_RectangleCorners[i][3] == 0)
    {
      found = false;
    }
  }
  return found;
}

bool AdaptiveAlignment::find_interface_edges()
{
  bool found = true;

  QString ss = QObject::tr("Finding Edge Of The Sample");
  notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getImageDataArrayPath());
  //QString attrMatName = getImageDataArrayPath().getAttributeMatrixName();

  size_t udims[3] = {0, 0, 0};
  std::tie(udims[0], udims[1], udims[2]) = m->getGeometryAs<ImageGeom>()->getDimensions();

  uint64_t dims[3] = {
      static_cast<uint64_t>(udims[0]), static_cast<uint64_t>(udims[1]), static_cast<uint64_t>(udims[2]),
  };

  // loop over slices
  for(int i = 0; i < dims[2]; i++)
  {
    // mean value of each row of pixels between bottom edge of the green rectangle and center of the circle is found
    uint64_t xrange1 = m_RectangleCorners[i][0];
    uint64_t xrange2 = m_RectangleCorners[i][2];
    uint64_t yrange1 = m_RectangleCorners[i][3] + 1;
    uint64_t yrange2 = static_cast<uint64_t>(floor(m_CalibratingCircles[i][1]));

    if(yrange1 >= yrange2)
      {return false;}

    std::vector<float> mean(yrange2 - yrange1 + 1, 0);
    for(uint64_t y = yrange1; y <= yrange2; y++){
      for(uint64_t x = xrange1; x <= xrange2; x++)
      {
        mean[y - yrange1] += m_FlatImageData[i * dims[0] * dims[1] + y * dims[0] + x] / static_cast<float>(mean.size());
      }
    }
    // maximum difference between two rows of vertical distance 2 is found
    float maxdif = 0.0f;
    uint64_t argmaxdif;
    for(uint64_t y = yrange1; y <= yrange2 - 2; y++) {
      if(std::fabs(mean[y - yrange1 + 2] - mean[y - yrange1]) > maxdif)
      {
        maxdif = std::fabs(mean[y - yrange1 + 2] - mean[y - yrange1]);
        argmaxdif = y;
      }
    }
    m_InterfaceEdges[i] = argmaxdif;
  }

  return found;
}

void AdaptiveAlignment::estimate_shifts_from_images(std::vector<float>& xneedshifts, std::vector<float>& yneedshifts)
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getImageDataArrayPath());

  notifyStatusMessage(getHumanLabel(), "Estimate shifts from images");

  size_t udims[3] = {0, 0, 0};
  std::tie(udims[0], udims[1], udims[2]) = m->getGeometryAs<ImageGeom>()->getDimensions();

  uint64_t dims[3] = {
      static_cast<uint64_t>(udims[0]), static_cast<uint64_t>(udims[1]), static_cast<uint64_t>(udims[2]),
  };

  std::vector<float> RectangleTop(2);
  std::vector<float> RectangleBottom(2);
  std::vector<float> Edge(2);

  // linear regression for smoothed position of the top (T) and bottom of green rectangles (B) and edge of the sample (E) according to the image index (X)
  // positions of the interface edge are cleaned from noise by subtracting the y-coordinate of the calibrating circles
  // T = RectangleTop[0] * X + RectangleTop[1]
  // B = RectangleBottom[0] * X + RectangleBottom[1]
  // E = Edge[0] * X + Edge[1]

  float sumX = 0.0f;  // sum(x_i)
  float sumX2 = 0.0f; // sum(x_i^2)

  float sumT = 0.0f;  // sum(y_i)
  float sumXT = 0.0f; // sum(x_i * y_i)

  float sumB = 0.0f;  // sum(y_i)
  float sumXB = 0.0f; // sum(x_i * y_i)

  float sumE = 0.0f;  // sum(z_i)
  float sumXE = 0.0f; // sum(x_i * z_i)

  for(uint64_t i = 0; i < dims[2]; i++)
  {
    sumX += static_cast<float>(i);
    sumX2 += static_cast<float>(i * i);

    sumT += static_cast<float>(m_RectangleCorners[i][1]);
    sumXT += static_cast<float>(i * m_RectangleCorners[i][1]);

    sumB += static_cast<float>(m_RectangleCorners[i][3]);
    sumXB += static_cast<double>(i * m_RectangleCorners[i][3]);

    sumE += static_cast<float>(m_InterfaceEdges[i] - m_CalibratingCircles[i][1]);
    sumXE += static_cast<double>(i * (m_InterfaceEdges[i] - m_CalibratingCircles[i][1]));
  }

  RectangleTop[0] = ((dims[2] * sumXT - sumX * sumT) / (dims[2] * sumX2 - sumX * sumX));
  RectangleTop[1] = ((sumX2 * sumT - sumX * sumXT) / (dims[2] * sumX2 - sumX * sumX));

  RectangleBottom[0] = ((dims[2] * sumXB - sumX * sumB) / (dims[2] * sumX2 - sumX * sumX));
  RectangleBottom[1] = ((sumX2 * sumB - sumX * sumXB) / (dims[2] * sumX2 - sumX * sumX));

  Edge[0] = ((dims[2] * sumXE - sumX * sumE) / (dims[2] * sumX2 - sumX * sumX));
  Edge[1] = ((sumX2 * sumE - sumX * sumXE) / (dims[2] * sumX2 - sumX * sumX));

  // width and height of the rectangle
  // width should be constant, height is obtained as average of the smoothed values
  float rectangle_width = static_cast<float>(m_RectangleCorners[0][2] - m_RectangleCorners[0][0]);
  float rectangle_height = 0;
  for(uint64_t i = 0; i < dims[2]; i++)
  {
    rectangle_height += (i * (RectangleBottom[0] - RectangleTop[0]) + RectangleBottom[1] - RectangleTop[1]) / dims[2];
  }

  float Slope = 0.5f * (RectangleTop[0] + RectangleBottom[0]);

  // xneedshifts: how the calibrating circles moved horizontally
  // yneedshifts: how the vertical distance of the rectangle and the interface edge changed
  // this goes with minus sign because that's how we need to correct the current position of the slices
  for(uint32_t i = 0; i < dims[2] - 1; i++)
  {
    xneedshifts[dims[2] - 2 - i] = -(m_CalibratingCircles[i + 1][0] - m_CalibratingCircles[i][0]);
    yneedshifts[dims[2] - 2 - i] = -(m_CalibratingCircles[i + 1][1] - m_CalibratingCircles[i][1] + Slope - Edge[0]);

    // convert to EBSD-voxel units
    xneedshifts[dims[2] - 2 - i] *= dims[0] / rectangle_width;
    yneedshifts[dims[2] - 2 - i] *= dims[1] / rectangle_height;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignment::find_shifts(std::vector<int64_t>& xshifts, std::vector<int64_t>& yshifts, std::vector<float>& xneedshifts, std::vector<float>& yneedshifts)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void initializeArrayValues(IDataArray::Pointer p, size_t index)
{

  typename DataArray<T>::Pointer ptr = std::dynamic_pointer_cast<DataArray<T>>(p);
  T var = static_cast<T>(0);
  ptr->initializeTuple(index, &var);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignment::execute()
{
  setErrorCondition(0);
  setWarningCondition(0);
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getDataContainerName());

  size_t udims[3] = {0, 0, 0};
  std::tie(udims[0], udims[1], udims[2]) = m->getGeometryAs<ImageGeom>()->getDimensions();

  uint64_t dims[3] = {
      static_cast<uint64_t>(udims[0]), static_cast<uint64_t>(udims[1]), static_cast<uint64_t>(udims[2]),
  };

  if(dims[0] == 0.0f || dims[1] == 0.0f || dims[2] == 0.0f)
  {
    QString ss = QObject::tr("Dimensions were not recognized correctly.");
    setErrorCondition(-99999999);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  float res[3] = {0.0f, 0.0f, 0.0f};
  m->getGeometryAs<ImageGeom>()->getResolution(res);

  if(res[0] == 0.0f || res[1] == 0.0f || res[2] == 0.0f)
  {
    QString ss = QObject::tr("Resolution was not recognized correctly.");
    setErrorCondition(-99999999);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  std::vector<float> xneedshifts;
  std::vector<float> yneedshifts;

  if(m_GlobalCorrection == 2)
  {
    // user-defined shifts related to one pair of consecutive sections and converted to voxels
    float xinitvalue = m_ShiftX / static_cast<float>(dims[2] / res[2]);
    float yinitvalue = m_ShiftY / static_cast<float>(dims[2] / res[2]);
    xneedshifts.resize(1, xinitvalue);
    yneedshifts.resize(1, yinitvalue);
  }

  // estimate shifts between slices from SEM images
  if(m_GlobalCorrection == 1)
  {
    xneedshifts.resize(dims[2] - 1);
    yneedshifts.resize(dims[2] - 1);

    m_CalibratingCircles.resize(dims[2]);
    m_RectangleCorners.resize(dims[2]);
    m_InterfaceEdges.resize(dims[2]);
    for(uint64_t i = 0; i < dims[2]; i++)
    {
      m_CalibratingCircles[i].resize(2);
      m_RectangleCorners[i].resize(4);
    }

    // find corners of green rectangles bounding the EBSD scanned area
    if(!find_rectangles())
    {
      QString ss = QObject::tr("Area of EBSD mapping could not be identified from the input images.");
      setErrorCondition(-1);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
      return;
    }

    // create 3D grayscale image
    create_array_for_flattened_image();
    flatten_image();

    // find calibrating circles in each 2D image
    if(!find_calibrating_circles())
    {
      QString ss = QObject::tr("Calibrating circles could not be identified from the input images.");
      setErrorCondition(-1);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
      return;
    }

    // find interface edge (edge of the specimen between the scanning plane and milling plane)
    if(!find_interface_edges())
    {
      QString ss = QObject::tr("Edge of the sample could not be identified from the input images.");
      setErrorCondition(-1);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
      return;
    }

    // estimate shifts from SEM images
    estimate_shifts_from_images(xneedshifts, yneedshifts);

    // delete greyscale image
    delete_array_for_flattened_image();
  }

  std::vector<int64_t> xshifts(dims[2], 0);
  std::vector<int64_t> yshifts(dims[2], 0);

  // estimate shifts between slices from EBSD images
  find_shifts(xshifts, yshifts, xneedshifts, yneedshifts);

  uint64_t xspot = 0, yspot = 0;
  uint64_t newPosition = 0;
  uint64_t currentPosition = 0;

  QList<QString> voxelArrayNames = m->getAttributeMatrix(getCellAttributeMatrixName())->getAttributeArrayNames();
  for(const auto& dataArrayPath : m_IgnoredDataArrayPaths)
  {
    voxelArrayNames.removeAll(dataArrayPath.getDataArrayName());
  }

  int64_t progIncrement = dims[2] / 100;
  int64_t prog = 1;
  int64_t progressInt = 0;
  uint64_t slice = 0;

  // transfer cell data
  for(uint64_t i = 1; i < dims[2]; i++)
  {
    if(i > prog)
    {

      progressInt = (static_cast<float>(i) / dims[2]) * 100.0f;
      QString ss = QObject::tr("Transferring Cell Data || %1% Complete").arg(progressInt);
      notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);
      prog = prog + progIncrement;
    }
    if(getCancel())
    {
      return;
    }
    slice = (dims[2] - 1) - i;
    for(uint64_t l = 0; l < dims[1]; l++)
    {
      for(uint64_t n = 0; n < dims[0]; n++)
      {
        if(yshifts[i] >= 0)
        {
          yspot = l;
        }
        else if(yshifts[i] < 0)
        {
          yspot = dims[1] - 1 - l;
        }
        if(xshifts[i] >= 0)
        {
          xspot = n;
        }
        else if(xshifts[i] < 0)
        {
          xspot = dims[0] - 1 - n;
        }
        newPosition = (slice * dims[0] * dims[1]) + (yspot * dims[0]) + xspot;
        currentPosition = (slice * dims[0] * dims[1]) + ((yspot + yshifts[i]) * dims[0]) + (xspot + xshifts[i]);
        if(int64_t((yspot + yshifts[i])) >= 0 && (yspot + yshifts[i]) <= dims[1] - 1 && int64_t((xspot + xshifts[i])) >= 0 && (xspot + xshifts[i]) <= dims[0] - 1)
        {
          for(const auto& arrayName : voxelArrayNames)
          {
            IDataArray::Pointer p = m->getAttributeMatrix(getCellAttributeMatrixName())->getAttributeArray(arrayName);
            p->copyTuple(currentPosition, newPosition);
          }
        }
        if(int64_t((yspot + yshifts[i])) < 0 || (yspot + yshifts[i]) > dims[1] - 1 || int64_t((xspot + xshifts[i])) < 0 || (xspot + xshifts[i]) > dims[0] - 1)
        {
          for(const auto& arrayName : voxelArrayNames)
          {
            IDataArray::Pointer p = m->getAttributeMatrix(getCellAttributeMatrixName())->getAttributeArray(arrayName);
            EXECUTE_FUNCTION_TEMPLATE(this, initializeArrayValues, p, p, newPosition)
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer AdaptiveAlignment::newFilterInstance(bool copyFilterParameters) const
{
  AdaptiveAlignment::Pointer filter = AdaptiveAlignment::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignment::getCompiledLibraryName() const
{
  return AnisotropyConstants::AnisotropyBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignment::getBrandingString() const
{
  return "Anisotropy";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignment::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << SIMPLib::Version::Major() << "." << SIMPLib::Version::Minor() << "." << SIMPLib::Version::Patch();
  return version;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignment::getGroupName() const
{
  return SIMPL::FilterGroups::ReconstructionFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid AdaptiveAlignment::getUuid()
{
  return QUuid("{b7394157-f142-59af-902d-7a29ff4f3ccd}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignment::getSubGroupName() const
{
  return AnisotropyConstants::FilterSubGroups::AnisotropicAlignment;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignment::getHumanLabel() const
{
  return "Adaptive Alignment";
}
