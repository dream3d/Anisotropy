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

#include <QtCore/QCoreApplication>
#include <QtCore/QDir>

#include <QtGui/QImage>

#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/FilterParameters/FileListInfoFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/ComparisonInputs.h"
#include "SIMPLib/Filtering/FilterFactory.hpp"
#include "SIMPLib/Filtering/FilterManager.h"
#include "SIMPLib/Filtering/FilterPipeline.h"
#include "SIMPLib/Filtering/QMetaObjectUtilities.h"
#include "SIMPLib/Plugin/ISIMPLibPlugin.h"
#include "SIMPLib/Plugin/SIMPLibPluginLoader.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Utilities/FilePathGenerator.h"
#include "UnitTestSupport.hpp"

#include "AnisotropyTestFileLocations.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void RemoveTestFiles()
{
#if REMOVE_TEST_FILES
  /*
  // remove input files
  QFile::remove(UnitTest::AnisotropyTest::TestInput);
  bool hasMissingFiles = false;
  bool orderAscending = true;
  QVector<QString> fileList = FilePathGenerator::GenerateFileList(UnitTest::AnisotropyTest::TestTifStartIndex,
    UnitTest::AnisotropyTest::TestTifEndIndex, hasMissingFiles, orderAscending,
    UnitTest::TestDataDir, UnitTest::AnisotropyTest::TestTifPrefix,
    UnitTest::AnisotropyTest::TestTifSuffix, UnitTest::AnisotropyTest::TestTifExtension,
    UnitTest::AnisotropyTest::TestTifPaddingDigits);
  if (fileList.size() == 0)
  {
    for (QVector<QString>::iterator filepath = fileList.begin(); filepath != fileList.end(); ++filepath)
    {
      QString imageFName = *filepath;
      QFile::remove(imageFName);
    }
  }
  */

  // remove output files
  QFile::remove(UnitTest::AnisotropyTest::TestOutput1);
  QFile::remove(UnitTest::AnisotropyTest::TestOutput2);
  QFile::remove(UnitTest::AnisotropyTest::TestOutput3);
  QFile::remove(UnitTest::AnisotropyTest::TestOutput4);
  QFile::remove(UnitTest::AnisotropyTest::TestOutput5);
#endif
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int TestFilterAvailability()
{
  {
    // Now instantiate the AdaptiveAlignmentFeature Filter from the FilterManager
    QString filtName = "AdaptiveAlignmentFeature";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);
    if (nullptr == filterFactory.get())
    {
      std::stringstream ss;
      ss << "Unable to initialize the AdaptiveAlignmentFeature while executing the AnisotropyTest.";
      DREAM3D_TEST_THROW_EXCEPTION(ss.str())
    }
  }

  {
    // Now instantiate the AdaptiveAlignmentMisorientation Filter from the FilterManager
    QString filtName = "AdaptiveAlignmentMisorientation";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);
    if (nullptr == filterFactory.get())
    {
      std::stringstream ss;
      ss << "Unable to initialize the AdaptiveAlignmentMisorientation while executing the AnisotropyTest.";
      DREAM3D_TEST_THROW_EXCEPTION(ss.str())
    }
  }

  {
    // Now instantiate the AdaptiveAlignmentMutualInformation Filter from the FilterManager
    QString filtName = "AdaptiveAlignmentMutualInformation";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);
    if (nullptr == filterFactory.get())
    {
      std::stringstream ss;
      ss << "Unable to initialize the AdaptiveAlignmentMutualInformation while executing the AnisotropyTest.";
      DREAM3D_TEST_THROW_EXCEPTION(ss.str())
    }
  }

  {
    // Now instantiate the SteinerCompact Filter from the FilterManager
    QString filtName = "SteinerCompact";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);
    if (nullptr == filterFactory.get())
    {
      std::stringstream ss;
      ss << "Unable to initialize the SteinerCompact while executing the AnisotropyTest.";
      DREAM3D_TEST_THROW_EXCEPTION(ss.str())
    }
  }

  return 0;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void addReadH5EBSDFilter(FilterPipeline::Pointer pipeline)
{
  QString filtName = "ReadH5Ebsd";
  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

  DataContainerArray::Pointer dca = DataContainerArray::New();

  if (nullptr != filterFactory.get())
  {
    AbstractFilter::Pointer filter = filterFactory->create();
    bool propWasSet;
    QVariant var;

    propWasSet = filter->setProperty("InputFile", UnitTest::AnisotropyTest::TestInput);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("ZStartIndex", UnitTest::AnisotropyTest::TestTifStartIndex);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("ZEndIndex", UnitTest::AnisotropyTest::TestTifEndIndex);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    QSet<QString> SelectedArrays = { "Phases", "EulerAngles", "Image Quality" };
    var.setValue(SelectedArrays);
    propWasSet = filter->setProperty("SelectedArrayNames", var);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("UseTransformations", true);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    filter->setDataContainerArray(dca);
    filter->preflight();
    int err = filter->getErrorCondition();
    DREAM3D_REQUIRE_EQUAL(err, 0);

    pipeline->pushBack(filter);
  }
  else
  {
    QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
    DREAM3D_REQUIRE_EQUAL(0, 1)
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void addImportImageStackFilter(FilterPipeline::Pointer pipeline)
{
  QString filtName = "ImportImageStack";
  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

  if (nullptr != filterFactory.get())
  {
    AbstractFilter::Pointer filter = filterFactory->create();
    bool propWasSet;
    QVariant var;

    FileListInfo_t input;
    input.InputPath = UnitTest::TestDataDir  + "/" + UnitTest::AnisotropyTest::TestTifExtension;
    input.StartIndex = UnitTest::AnisotropyTest::TestTifStartIndex;
    input.EndIndex = UnitTest::AnisotropyTest::TestTifEndIndex;
    input.FilePrefix = UnitTest::AnisotropyTest::TestTifPrefix;
    input.FileSuffix = UnitTest::AnisotropyTest::TestTifSuffix;
    input.IncrementIndex = 1;
    input.FileExtension = UnitTest::AnisotropyTest::TestTifExtension;
    input.PaddingDigits = UnitTest::AnisotropyTest::TestTifPaddingDigits;
    input.Ordering = 0;
    var.setValue(input);
    propWasSet = filter->setProperty("InputFileListInfo", var);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("GeometryType", 0);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    FloatVec3_t origin = { 0, 0, 0 };
    var.setValue(origin);
    propWasSet = filter->setProperty("Origin", var);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    FloatVec3_t resolution = { 1, 1, 1 };
    var.setValue(resolution);
    propWasSet = filter->setProperty("Resolution", var);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("DataContainerName", "SEMImageDataContainer");
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)

    pipeline->pushBack(filter);
  }
  else
  {
    QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
    DREAM3D_REQUIRE_EQUAL(0, 1)
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void addMultiThresholdObjectsFilter(FilterPipeline::Pointer pipeline)
{
  QString filtName = "MultiThresholdObjects";
  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

  if (nullptr != filterFactory.get())
  {
    AbstractFilter::Pointer filter = filterFactory->create();
    QVariant var;
    bool propWasSet;

    ComparisonInputs Thresholds;
    Thresholds.addInput("ImageDataContainer", "CellData", "Image Quality", 1, 800.0);

    var.setValue(Thresholds);
    propWasSet = filter->setProperty("SelectedThresholds", var);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("DestinationArrayName", SIMPL::GeneralData::Mask);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)

    pipeline->pushBack(filter);
  }
  else
  {
    QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
    DREAM3D_REQUIRE_EQUAL(0, 1)
  }
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void addConvertOrientationsFilter(FilterPipeline::Pointer pipeline)
{
  QString filtName = "ConvertOrientations";
  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

  if (nullptr != filterFactory.get())
  {
    AbstractFilter::Pointer filter = filterFactory->create();
    bool propWasSet;
    QVariant var;

    propWasSet = filter->setProperty("InputType", 0);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("OutputType", 2);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)

    DataArrayPath path("ImageDataContainer", "CellData", "EulerAngles");
    var.setValue(path);
    propWasSet = filter->setProperty("InputOrientationArrayPath", var);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("OutputOrientationArrayName", "Quats");
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)

    pipeline->pushBack(filter);
  }
  else
  {
    QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
    DREAM3D_REQUIRE_EQUAL(0, 1)
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void addEBSDSegmentFeatures(FilterPipeline::Pointer pipeline)
{
  QString filtName = "EBSDSegmentFeatures";
  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

  if (nullptr != filterFactory.get())
  {
    AbstractFilter::Pointer filter = filterFactory->create();
    bool propWasSet;

    propWasSet = filter->setProperty("MisorientationTolerance", 5.0);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("UseGoodVoxels", false);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)

    pipeline->pushBack(filter);
  }
  else
  {
    QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
    DREAM3D_REQUIRE_EQUAL(0, 1)
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void addAdaptiveAlignmentFeatureFilter(FilterPipeline::Pointer pipeline)
{
  QString filtName = "AdaptiveAlignmentFeature";
  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

  if (nullptr != filterFactory.get())
  {
    AbstractFilter::Pointer filter = filterFactory->create();
    bool propWasSet;

  propWasSet = filter->setProperty("GlobalCorrection", 2);
  DREAM3D_REQUIRE_EQUAL(propWasSet, true)
  propWasSet = filter->setProperty("ShiftX", 0);
  DREAM3D_REQUIRE_EQUAL(propWasSet, true)
  propWasSet = filter->setProperty("ShiftY", 0);
  DREAM3D_REQUIRE_EQUAL(propWasSet, true)
  propWasSet = filter->setProperty("WriteAlignmentShifts", true);
  DREAM3D_REQUIRE_EQUAL(propWasSet, true)
  propWasSet = filter->setProperty("AlignmentShiftFileName", UnitTest::AnisotropyTest::TestOutput1);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)

  pipeline->pushBack(filter);
  }
  else
  {
    QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
    DREAM3D_REQUIRE_EQUAL(0, 1)
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void addAdaptiveAlignmentMisorientationFilter(FilterPipeline::Pointer pipeline)
{
  QString filtName = "AdaptiveAlignmentMisorientation";
  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

  if (nullptr != filterFactory.get())
  {
    AbstractFilter::Pointer filter = filterFactory->create();

    QVariant var;
    bool propWasSet;

    propWasSet = filter->setProperty("GlobalCorrection", 1);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("MisorientationTolerance", 5.0f);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("UseGoodVoxels", false);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)

    DataArrayPath path("SEMImageDataContainer", "CellData", "ImageData");
    var.setValue(path);
    propWasSet = filter->setProperty("ImageDataArrayPath", var);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("WriteAlignmentShifts", true);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("AlignmentShiftFileName", UnitTest::AnisotropyTest::TestOutput2);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)

    pipeline->pushBack(filter);
  }
  else
  {
    QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
    DREAM3D_REQUIRE_EQUAL(0, 1)
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void addAdaptiveAlignmentMutualInformationFilter(FilterPipeline::Pointer pipeline)
{
  QString filtName = "AdaptiveAlignmentMutualInformation";
  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

  if (nullptr != filterFactory.get())
  {
    AbstractFilter::Pointer filter = filterFactory->create();
    bool propWasSet;

    propWasSet = filter->setProperty("GlobalCorrection", 0);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("UseGoodVoxels", false);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("MisorientationTolerance", 5.0f);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("WriteAlignmentShifts", true);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("AlignmentShiftFileName", UnitTest::AnisotropyTest::TestOutput3);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)

    pipeline->pushBack(filter);
  }
  else
  {
    QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
    DREAM3D_REQUIRE_EQUAL(0, 1)
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void addSteinerCompactFilter(FilterPipeline::Pointer pipeline)
{
  QString filtName = "SteinerCompact";
  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer filterFactory = fm->getFactoryFromClassName(filtName);

  if (nullptr != filterFactory.get())
  {
    AbstractFilter::Pointer filter = filterFactory->create();
    bool propWasSet;

    propWasSet = filter->setProperty("Plane", 0);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("Sites", 1);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("VtkOutput", true);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("VtkFileName", UnitTest::AnisotropyTest::TestOutput4);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("TxtOutput", true);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
    propWasSet = filter->setProperty("TxtFileName", UnitTest::AnisotropyTest::TestOutput5);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)

    pipeline->pushBack(filter);
  }
  else
  {
    QString ss = QObject::tr("AnisotropyTest Error creating filter '%1'. Filter was not created/executed. Please notify the developers.").arg(filtName);
    DREAM3D_REQUIRE_EQUAL(0, 1)
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int TestAnisotropy()
{
  FilterPipeline::Pointer pipeline = FilterPipeline::New();

  {
    addReadH5EBSDFilter(pipeline);
    addImportImageStackFilter(pipeline);
    addConvertOrientationsFilter(pipeline);
    addAdaptiveAlignmentMisorientationFilter(pipeline);
    addMultiThresholdObjectsFilter(pipeline);
    addAdaptiveAlignmentFeatureFilter(pipeline);
    addAdaptiveAlignmentMutualInformationFilter(pipeline);
    addEBSDSegmentFeatures(pipeline);
    addSteinerCompactFilter(pipeline);
    pipeline->execute();
    DREAM3D_REQUIRE_EQUAL(pipeline->getErrorCondition(), 0)
  }

  return EXIT_SUCCESS;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void loadFilterPlugins()
{
  // Register all the filters including trying to load those from Plugins
  FilterManager* fm = FilterManager::Instance();
  SIMPLibPluginLoader::LoadPluginFilters(fm);

  // Send progress messages from PipelineBuilder to this object for display
  QMetaObjectUtilities::RegisterMetaTypes();
}


// -----------------------------------------------------------------------------
//  Use test framework
// -----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Instantiate the QCoreApplication that we need to get the current path and load plugins.

  QCoreApplication app(argc, argv);
  QCoreApplication::setOrganizationName("Czech Academy of Sciences");
  QCoreApplication::setOrganizationDomain("http://ams.fzu.cz");
  QCoreApplication::setApplicationName("AnisotropyTest");

  int err = EXIT_SUCCESS;

  DREAM3D_REGISTER_TEST( loadFilterPlugins() );
  DREAM3D_REGISTER_TEST( TestFilterAvailability() );
  DREAM3D_REGISTER_TEST(TestAnisotropy())
  DREAM3D_REGISTER_TEST(RemoveTestFiles())

  PRINT_TEST_SUMMARY();

  return err;
}

