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

#include "AdaptiveAlignmentFeature.h"

#include <fstream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/SIMPLibVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AdaptiveAlignmentFeature::AdaptiveAlignmentFeature()
: m_GoodVoxelsArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::Mask)
, m_GoodVoxels(nullptr)
{
  // only setting up the child parameters because the parent constructor has already been called
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AdaptiveAlignmentFeature::~AdaptiveAlignmentFeature() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentFeature::setupFilterParameters()
{
  // getting the current parameters that were set by the parent and adding to it before resetting it
  AdaptiveAlignment::setupFilterParameters();

  FilterParameterVector parameters = getFilterParameters();
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", GoodVoxelsArrayPath, FilterParameter::RequiredArray, AdaptiveAlignmentFeature, req));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentFeature::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  AdaptiveAlignment::readFilterParameters(reader, index);
  reader->openFilterGroup(this, index);
  setGoodVoxelsArrayPath(reader->readDataArrayPath("GoodVoxelsArrayPath", getGoodVoxelsArrayPath()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentFeature::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentFeature::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);

  // Set the DataContainerName and AttributematrixName for the Parent Class (AdaptiveAlignment) to Use.
  // These are checked for validity in the Parent Class dataCheck
  setDataContainerName(m_GoodVoxelsArrayPath.getDataContainerName());
  setCellAttributeMatrixName(m_GoodVoxelsArrayPath.getAttributeMatrixName());

  AdaptiveAlignment::dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  QVector<size_t> cDims(1, 1);
  m_GoodVoxelsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>, AbstractFilter>(this, getGoodVoxelsArrayPath(),
                                                                                                     cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_GoodVoxelsPtr.lock())                                                                      /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_GoodVoxels = m_GoodVoxelsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentFeature::preflight()
{
  setInPreflight(true);
  emit preflightAboutToExecute();
  emit updateFilterParameters(this);
  dataCheck();
  emit preflightExecuted();
  setInPreflight(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentFeature::find_shifts(std::vector<int64_t>& xshifts, std::vector<int64_t>& yshifts, std::vector<float>& xneedshifts, std::vector<float>& yneedshifts)
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getDataContainerName());

  size_t udims[3] = {0, 0, 0};
  std::tie(udims[0], udims[1], udims[2]) = m->getGeometryAs<ImageGeom>()->getDimensions();

  int64_t dims[3] = {
      static_cast<int64_t>(udims[0]), static_cast<int64_t>(udims[1]), static_cast<int64_t>(udims[2]),
  };

  uint64_t maxstoredshifts = 1;
  if(xneedshifts.size() > 0)
    maxstoredshifts = 20;

  float disorientation = 0.0f;

  std::vector<std::vector<int64_t>> newxshift(dims[2]);
  std::vector<std::vector<int64_t>> newyshift(dims[2]);
  std::vector<std::vector<float>> mindisorientation(dims[2]);
  for(size_t a = 1; a < udims[2]; a++)
  {
    newxshift[a].resize(maxstoredshifts, 0);
    newyshift[a].resize(maxstoredshifts, 0);
    mindisorientation[a].resize(maxstoredshifts, std::numeric_limits<float>::max());
  }

  // Allocate a 2D Array which will be reused from slice to slice
  // second dimension is assigned in each cycle separately
  std::vector<std::vector<bool>> misorients(dims[0]);

  const uint64_t halfDim0 = static_cast<uint64_t>(dims[0] * 0.5f);
  const uint64_t halfDim1 = static_cast<uint64_t>(dims[1] * 0.5f);

  int32_t oldxshift = 0;
  int32_t oldyshift = 0;
  float count = 0.0f;
  int64_t slice = 0;
  int64_t refposition = 0;
  int64_t curposition = 0;
  uint64_t progInt = 0;

  for(int64_t iter = 1; iter < dims[2]; iter++)
  {
    progInt = static_cast<uint64_t>(iter * 100 / static_cast<float>(dims[2]));
    QString ss = QObject::tr("Aligning Anisotropic Sections || Determining Shifts || %1% Complete").arg(progInt);

    slice = (dims[2] - 1) - iter;
    oldxshift = -1;
    oldyshift = -1;

    for(size_t i = 0; i < udims[0]; i++)
    {
      misorients[i].assign(dims[1], false);
    }

    while(newxshift[iter][0] != oldxshift || newyshift[iter][0] != oldyshift)
    {
      oldxshift = newxshift[iter][0];
      oldyshift = newyshift[iter][0];

      for(int32_t j = -3; j < 4; j++)
      {
        for(int32_t k = -3; k < 4; k++)
        {
          disorientation = 0.0f;
          count = 0.0f;
          if((llabs(k + oldxshift) < static_cast<long long int>(halfDim0)) && llabs(j + oldyshift) < static_cast<long long int>(halfDim1) &&
             misorients[k + oldxshift + halfDim0][j + oldyshift + halfDim1] == false)
          {
            for(int64_t l = 0; l < dims[1]; l = l + 4)
            {
              for(int64_t n = 0; n < dims[0]; n = n + 4)
              {
                if((l + j + oldyshift) >= 0 && (l + j + oldyshift) < dims[1] && (n + k + oldxshift) >= 0 && (n + k + oldxshift) < dims[0])
                {
                  refposition = ((slice + 1) * dims[0] * dims[1]) + (l * dims[0]) + n;
                  curposition = (slice * dims[0] * dims[1]) + ((l + j + oldyshift) * dims[0]) + (n + k + oldxshift);
                  if(m_GoodVoxels[refposition] != m_GoodVoxels[curposition])
                  {
                    disorientation++;
                  }
                  count++;
                }
                else
                {
                }
              }
            }
            disorientation = disorientation / count;
            misorients[k + oldxshift + halfDim0][j + oldyshift + halfDim1] = true;

            // compare the new shift with currently stored ones
            int64_t s = maxstoredshifts;
            while(s - 1 >= 0 && disorientation < mindisorientation[iter][s - 1])
            {
              s--;
            }

            // new shift is stored with index 's' in the arrays
            if(s < static_cast<int64_t>(maxstoredshifts))
            {
              // lag the shifts already stored
              for(int64_t t = maxstoredshifts - 1; t > s; t--)
              {
                newxshift[iter][t] = newxshift[iter][t - 1];
                newyshift[iter][t] = newyshift[iter][t - 1];
                mindisorientation[iter][t] = mindisorientation[iter][t - 1];
              }
              // store the new shift
              newxshift[iter][s] = k + oldxshift;
              newyshift[iter][s] = j + oldyshift;
              mindisorientation[iter][s] = disorientation;
            }
          }
        }
      }
    }
    xshifts[iter] = xshifts[iter - 1] + newxshift[iter][0];
    yshifts[iter] = yshifts[iter - 1] + newyshift[iter][0];
  }

  std::vector<uint64_t> curindex(dims[2], 0);

  // find corrected shifts
  if(xneedshifts.size() > 0)
  {
    QString ss = QObject::tr("Aligning Anisotropic Sections || Correcting shifts");
    notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);

    std::vector<float> changedisorientation(dims[2], 0);
    std::vector<uint64_t> changeindex(dims[2], 0);
    std::vector<float> changeerror(dims[2], 0);

    std::vector<float> xshiftsest; // cumulative x-shifts estimated from SEM images
    std::vector<float> yshiftsest; // cumulative y-shifts estimated from SEM images

    float curerror = 0.0f;
    float tolerance = 0.0f;

    // evaluate error between current shifts and desired shifts
    if(xneedshifts.size() == 1) // error is computed as misagreement between slopes
    {
      tolerance = 1.0f / static_cast<float>(dims[2] - 1);
      curerror = compute_error1(dims[2], 0, xneedshifts[0], yneedshifts[0], newxshift, newyshift, curindex);
    }
    else if(xneedshifts.size() > 1) // error is computed as misagreement with shifts estimated from SEM images
    {
      tolerance = 1.0f;
      xshiftsest.resize(dims[2], 0);
      yshiftsest.resize(dims[2], 0);
      for(size_t iter = 1; iter < udims[2]; iter++)
      {
        xshiftsest[iter] = xshiftsest[iter - 1] + xneedshifts[iter - 1];
        yshiftsest[iter] = yshiftsest[iter - 1] + yneedshifts[iter - 1];
      }
      curerror = compute_error2(dims[2], 0, xshiftsest, yshiftsest, newxshift, newyshift, curindex);
    }

    // iterative selection of a candidate shift, recomputing of current candidates, evaluation of error
    if(curerror > tolerance)
    {
      float minchangedisorientation = 0;
      float minchangeerror = 0;
      int64_t minchangeindex = 0;
      int64_t minchangeiter = 0;
      float olderror = 0;
      float newerror = 0;
      uint64_t progInt = 0;

      do
      {
        QString ss = QObject::tr("Aligning Anisotropic Sections || Correcting Shifts || Iteration %1").arg(++progInt);
        ;
        notifyStatusMessage(getMessagePrefix(), getHumanLabel(), ss);
        if(getCancel() == true)
        {
          return;
        }

        olderror = curerror;
        for(size_t iter = 1; iter < udims[2]; iter++)
        {
          float newminerror = std::numeric_limits<float>::max();
          float newmindisorientation = std::numeric_limits<float>::max();
          uint64_t newminindex = 0;
          for(uint64_t index = curindex[iter] + 1; index < maxstoredshifts; index++)
          {
            // recompute error for the configuration with this candidate changed
            if(xneedshifts.size() == 1)
            {
              newerror = compute_error1(iter, index, xneedshifts[0], yneedshifts[0], newxshift, newyshift, curindex);
            }
            else if(xneedshifts.size() > 1)
            {
              newerror = compute_error2(iter, index, xshiftsest, yshiftsest, newxshift, newyshift, curindex);
            }

            // compare the new error with the best current error
            if(newerror < curerror && mindisorientation[iter][index] - mindisorientation[iter][0] < newmindisorientation)
            {
              newminerror = newerror;
              newminindex = index;
              newmindisorientation = mindisorientation[iter][index] - mindisorientation[iter][0];
            }
          }
          // assign best error, corresponding index and disorientation value for this slice
          changeerror[iter] = newminerror;
          changeindex[iter] = newminindex;
          changedisorientation[iter] = newmindisorientation;
        }

        // among all slices, find the best candidate (with minimum disorientation change)
        minchangedisorientation = std::numeric_limits<float>::max() - 1;
        minchangeerror = std::numeric_limits<float>::max();
        minchangeindex = 0;
        minchangeiter = 0;
        for(size_t iter = 1; iter < udims[2]; iter++)
        {
          if(changeerror[iter] < curerror && (changedisorientation[iter] < minchangedisorientation || (changedisorientation[iter] == minchangedisorientation &&
                                                                                                       llabs(newxshift[iter][changeindex[iter]]) + llabs(newyshift[iter][changeindex[iter]]) <
                                                                                                           llabs(newxshift[iter][minchangeindex]) + llabs(newyshift[iter][minchangeindex]))))
          {
            minchangeiter = iter;
            minchangeindex = changeindex[iter];
            minchangedisorientation = changedisorientation[iter];
            minchangeerror = changeerror[iter];
          }
        }

        if(minchangeerror < curerror && minchangeerror >= tolerance)
        {
          // assign the best candidate
          changedisorientation[minchangeiter] = minchangedisorientation;
          curindex[minchangeiter] = minchangeindex;
          // reassign current error
          curerror = minchangeerror;
        }

      } while(minchangedisorientation < std::numeric_limits<float>::max() - 1 && curerror < olderror && curerror > tolerance);
    }
  }

  if(getWriteAlignmentShifts() == true)
  {
    std::ofstream outFile;
    outFile.open(getAlignmentShiftFileName().toLatin1().data());
    for(size_t iter = 1; iter < udims[2]; iter++)
    {
      slice = (dims[2] - 1) - iter;
      xshifts[iter] = xshifts[iter - 1] + newxshift[iter][curindex[iter]];
      yshifts[iter] = yshifts[iter - 1] + newyshift[iter][curindex[iter]];
      outFile << slice << "	" << slice + 1 << "	" << newxshift[iter][curindex[iter]] << "	" << newyshift[iter][curindex[iter]] << "	" << xshifts[iter] << "	" << yshifts[iter] << "\n";
    }
    outFile.close();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// find the error if the current shifts are changed by modification provided on slice 'iter'
// where currently selected shift is substituted by the shift stored in position 'index'
float AdaptiveAlignmentFeature::compute_error1(uint64_t iter, uint64_t index, float xneedtrend, float yneedtrend, std::vector<std::vector<int64_t>>& newxshift,
                                               std::vector<std::vector<int64_t>>& newyshift, std::vector<uint64_t>& curindex)
{
  int64_t n = curindex.size() - 1;

  int64_t xshifts = 0;
  int64_t yshifts = 0;

  double sumX = 0.0f;
  double sumX_2 = 0.0f;
  double x_sumY = 0.0f;
  double x_sumXY = 0.0f;
  double y_sumY = 0.0f;
  double y_sumXY = 0.0f;

  for(int64_t i = 1; i <= n; i++)
  {
    if(i != iter) // shifts without modification
    {
      xshifts += newxshift[i][curindex[i]];
      yshifts += newyshift[i][curindex[i]];
    }
    else // modified shift
    {
      xshifts += newxshift[i][index];
      yshifts += newyshift[i][index];
    }
    // all values are divided by n for better numerical stability
    sumX += static_cast<double>(i) / static_cast<double>(n);
    sumX_2 += static_cast<double>(i * i) / static_cast<double>(n * n);
    x_sumY += static_cast<double>(xshifts) / static_cast<double>(n);
    x_sumXY += static_cast<double>(i * xshifts) / static_cast<double>(n * n);
    y_sumY += static_cast<double>(yshifts) / static_cast<double>(n);
    y_sumXY += static_cast<double>(i * yshifts) / static_cast<double>(n * n);
  }

  double xtrend = static_cast<double>((n * x_sumXY - sumX * x_sumY) / (n * sumX_2 - sumX * sumX));
  double ytrend = static_cast<double>((n * y_sumXY - sumX * y_sumY) / (n * sumX_2 - sumX * sumX));

  // error is computed from angular deviation as
  // |xang1 - xang2| + |yang1 - yang2| = |arctg(xtrend) - arctg(xneedtrend)| + |arctg(ytrend) - arctg(yneedtrend)|
  return std::fabs(atan(xtrend) - atan(xneedtrend)) + std::fabs(atan(ytrend) - atan(yneedtrend));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
float AdaptiveAlignmentFeature::compute_error2(uint64_t iter, uint64_t index, std::vector<float>& xshiftsest, std::vector<float>& yshiftsest, std::vector<std::vector<int64_t>>& newxshift,
                                               std::vector<std::vector<int64_t>>& newyshift, std::vector<uint64_t>& curindex)
{
  uint64_t n = curindex.size() - 1;

  int64_t xshifts = 0;
  int64_t yshifts = 0;

  float error = 0;
  float xdif = 0;
  float ydif = 0;
  float divide = static_cast<float>(2 * n);

  for(uint64_t i = 1; i <= n; i++)
  {
    if(i != iter)
    {
      xshifts += newxshift[i][curindex[i]];
      yshifts += newyshift[i][curindex[i]];
    }
    else
    {
      xshifts += newxshift[i][index];
      yshifts += newyshift[i][index];
    }

    xdif = static_cast<float>(xshifts) - xshiftsest[i];
    ydif = static_cast<float>(yshifts) - yshiftsest[i];
    error += (xdif * xdif + ydif * ydif) / divide;
  }
  return error;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void AdaptiveAlignmentFeature::execute()
{
  setErrorCondition(0);
  setWarningCondition(0);
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  AdaptiveAlignment::execute();

  // If there is an error set this to something negative and also set a message
  notifyStatusMessage(getHumanLabel(), "Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer AdaptiveAlignmentFeature::newFilterInstance(bool copyFilterParameters) const
{
  AdaptiveAlignmentFeature::Pointer filter = AdaptiveAlignmentFeature::New();
  if(true == copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentFeature::getCompiledLibraryName() const
{
  return AnisotropyConstants::AnisotropyBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentFeature::getBrandingString() const
{
  return "Anisotropy";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentFeature::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << SIMPLib::Version::Major() << "." << SIMPLib::Version::Minor() << "." << SIMPLib::Version::Patch();
  return version;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentFeature::getGroupName() const
{
  return SIMPL::FilterGroups::ReconstructionFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid AdaptiveAlignmentFeature::getUuid()
{
  return QUuid("{28977b7a-5b88-5145-abcd-d1c933f7d975}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentFeature::getSubGroupName() const
{
  return AnisotropyConstants::FilterSubGroups::AnisotropicAlignment;
  ;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString AdaptiveAlignmentFeature::getHumanLabel() const
{
  return "Adaptive Alignment (Feature)";
}
