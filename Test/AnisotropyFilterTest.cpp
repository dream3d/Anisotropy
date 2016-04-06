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

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Common/FilterPipeline.h"
#include "SIMPLib/Common/FilterManager.h"
#include "SIMPLib/Common/FilterFactory.hpp"
#include "SIMPLib/Plugin/ISIMPLibPlugin.h"
#include "SIMPLib/Plugin/SIMPLibPluginLoader.h"
#include "SIMPLib/Utilities/UnitTestSupport.hpp"
#include "SIMPLib/Utilities/QMetaObjectUtilities.h"

#include "AnisotropyTestFileLocations.h"


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int TestFilterAvailability()
{
  {
    // Now instantiate the AdaptiveAlignmentFeature Filter from the FilterManager
    QString filtName = "AdaptiveAlignmentFeature";
    FilterManager* fm = FilterManager::Instance();
    IFilterFactory::Pointer filterFactory = fm->getFactoryForFilter(filtName);
    if (NULL == filterFactory.get())
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
	  IFilterFactory::Pointer filterFactory = fm->getFactoryForFilter(filtName);
	  if (NULL == filterFactory.get())
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
	  IFilterFactory::Pointer filterFactory = fm->getFactoryForFilter(filtName);
	  if (NULL == filterFactory.get())
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
	  IFilterFactory::Pointer filterFactory = fm->getFactoryForFilter(filtName);
	  if (NULL == filterFactory.get())
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
void TestAdaptiveAlignmentFeatureFilter()
{
  QString filtName = "AdaptiveAlignmentFeature";
  FilterManager* fm = FilterManager::Instance();
  IFilterFactory::Pointer filterFactory = fm->getFactoryForFilter(filtName);
  
  if (NULL != filterFactory.get())
  {
    // If we get this far, the Factory is good so creating the filter should not fail unless something has
    // horribly gone wrong in which case the system is going to come down quickly after this.
    AbstractFilter::Pointer filter = filterFactory->create();

    QVariant var;
    bool propWasSet;

	DataArrayPath path;
	QVector<size_t> tDims(1, 0);
	QVector<size_t> cDims(1);
	DataContainer::Pointer dc = DataContainer::New("MyDataContainer");
	AttributeMatrix::Pointer am = AttributeMatrix::New(tDims, "MyAttributeMatrix", DREAM3D::AttributeMatrixType::Cell);
	DoubleArrayType::Pointer ar = DoubleArrayType::CreateArray(tDims, cDims, "MyDoubleArray");
	path = DataArrayPath(dc->getName(), am->getName(), ar->getName());

	propWasSet = filter->setProperty("UseImages", true);
	DREAM3D_REQUIRE_EQUAL(propWasSet, true)
		
	propWasSet = filter->setProperty("UserDefinedShifts", false);
	DREAM3D_REQUIRE_EQUAL(propWasSet, true)
	
	propWasSet = filter->setProperty("ShiftX", 0);
	DREAM3D_REQUIRE_EQUAL(propWasSet, true)

	propWasSet = filter->setProperty("ShiftY", 0);
	DREAM3D_REQUIRE_EQUAL(propWasSet, true)
		
	var.setValue(path);
	propWasSet = filter->setProperty("ImageDataArrayPath", var);
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
		
	propWasSet = filter->setProperty("WriteAlignmentShifts", true);
	DREAM3D_REQUIRE_EQUAL(propWasSet, true)

	propWasSet = filter->setProperty("AlignmentShiftFileName", "output");
    DREAM3D_REQUIRE_EQUAL(propWasSet, true)
	
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
void TestAdaptiveAlignmentMisorientationFilter()
{
	QString filtName = "AdaptiveAlignmentMisorientation";
	FilterManager* fm = FilterManager::Instance();
	IFilterFactory::Pointer filterFactory = fm->getFactoryForFilter(filtName);

	if (NULL != filterFactory.get())
	{
		// If we get this far, the Factory is good so creating the filter should not fail unless something has
		// horribly gone wrong in which case the system is going to come down quickly after this.
		AbstractFilter::Pointer filter = filterFactory->create();

		QVariant var;
		bool propWasSet;

		DataArrayPath path;
		QVector<size_t> tDims(1, 0);
		QVector<size_t> cDims(1);
		DataContainer::Pointer dc = DataContainer::New("MyDataContainer");
		AttributeMatrix::Pointer am = AttributeMatrix::New(tDims, "MyAttributeMatrix", DREAM3D::AttributeMatrixType::Cell);
		DoubleArrayType::Pointer ar = DoubleArrayType::CreateArray(tDims, cDims, "MyDoubleArray");
		path = DataArrayPath(dc->getName(), am->getName(), ar->getName());

		propWasSet = filter->setProperty("UseImages", true);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		propWasSet = filter->setProperty("UserDefinedShifts", false);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		propWasSet = filter->setProperty("ShiftX", 0);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		propWasSet = filter->setProperty("ShiftY", 0);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		var.setValue(path);
		propWasSet = filter->setProperty("ImageDataArrayPath", var);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		propWasSet = filter->setProperty("WriteAlignmentShifts", true);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		propWasSet = filter->setProperty("AlignmentShiftFileName", "output");
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

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
void TestAdaptiveAlignmentMutualInformationFilter()
{
	QString filtName = "AdaptiveAlignmentMutualInformation";
	FilterManager* fm = FilterManager::Instance();
	IFilterFactory::Pointer filterFactory = fm->getFactoryForFilter(filtName);

	if (NULL != filterFactory.get())
	{
		// If we get this far, the Factory is good so creating the filter should not fail unless something has
		// horribly gone wrong in which case the system is going to come down quickly after this.
		AbstractFilter::Pointer filter = filterFactory->create();

		QVariant var;
		bool propWasSet;

		DataArrayPath path;
		QVector<size_t> tDims(1, 0);
		QVector<size_t> cDims(1);
		DataContainer::Pointer dc = DataContainer::New("MyDataContainer");
		AttributeMatrix::Pointer am = AttributeMatrix::New(tDims, "MyAttributeMatrix", DREAM3D::AttributeMatrixType::Cell);
		DoubleArrayType::Pointer ar = DoubleArrayType::CreateArray(tDims, cDims, "MyDoubleArray");
		path = DataArrayPath(dc->getName(), am->getName(), ar->getName());

		propWasSet = filter->setProperty("UseImages", true);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		propWasSet = filter->setProperty("UserDefinedShifts", false);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		propWasSet = filter->setProperty("ShiftX", 0);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		propWasSet = filter->setProperty("ShiftY", 0);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		var.setValue(path);
		propWasSet = filter->setProperty("ImageDataArrayPath", var);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		propWasSet = filter->setProperty("WriteAlignmentShifts", true);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		propWasSet = filter->setProperty("AlignmentShiftFileName", "output");
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

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
void TestSteinerCompactFilter()
{
	QString filtName = "SteinerCompact";
	FilterManager* fm = FilterManager::Instance();
	IFilterFactory::Pointer filterFactory = fm->getFactoryForFilter(filtName);

	if (NULL != filterFactory.get())
	{
		// If we get this far, the Factory is good so creating the filter should not fail unless something has
		// horribly gone wrong in which case the system is going to come down quickly after this.
		AbstractFilter::Pointer filter = filterFactory->create();

		QVariant var;
		bool propWasSet;


		propWasSet = filter->setProperty("Plane", 0);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		propWasSet = filter->setProperty("Sites", 12);
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

		propWasSet = filter->setProperty("VtkFileName", "output");
		DREAM3D_REQUIRE_EQUAL(propWasSet, true)

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
    TestAdaptiveAlignmentFeatureFilter();
	TestAdaptiveAlignmentMisorientationFilter();
	TestAdaptiveAlignmentMutualInformationFilter();
	TestSteinerCompactFilter();

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
  QCoreApplication::setOrganizationName("BlueQuartz Software");
  QCoreApplication::setOrganizationDomain("bluequartz.net");
  QCoreApplication::setApplicationName("AnisotropyTest");

  int err = EXIT_SUCCESS;
  
  DREAM3D_REGISTER_TEST( loadFilterPlugins() );
  DREAM3D_REGISTER_TEST( TestFilterAvailability() );

  DREAM3D_REGISTER_TEST( TestAnisotropy())

  PRINT_TEST_SUMMARY();
  
  return err;
}

