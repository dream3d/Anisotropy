#pragma once

#include <QtCore/QString>


#include "itkImage.h"
#include "itkRGBPixel.h"

#include "SIMPLib/DataArrays/DataArray.hpp"

#include "Anisotropy/AnisotropyConfig.h"

namespace AnisotropyConstants
{
  const QString AnisotropyPluginFile("AnisotropyPlugin");
  const QString AnisotropyPluginDisplayName("Anisotropy");
  const QString AnisotropyBaseName("Anisotropy");

  const QString Version("1.0");
  const QString CompatibilityVersion("1.0");

  const QString VendorName("Czech Academy of Sciences, Institute of Physics, Group of Bulk Nanomaterials and Interfaces");
  const QString URL("http://ams.fzu.cz");
  const QString Copyright("(C) 2016 Czech Academy of Sciences, v.v.i.");

  namespace FilterSubGroups
  {
    const QString AnisotropicAlignment("Anisotropic Alignment");
  }

  //define pixels for dream3d variable types
  typedef int8_t Int8PixelType;
  typedef uint8_t UInt8PixelType;
  typedef int16_t Int16PixelType;
  typedef uint16_t UInt16PixelType;
  typedef int32_t Int32PixelType;
  typedef uint32_t UInt32PixelType;
  typedef int64_t Int64PixelType;
  typedef uint64_t UInt64PixelType;
  typedef float FloatPixelType;
  typedef double DoublePixelType;

  //define default pixel type
#if Anisotropy_BitDepth == 8
  typedef UInt8PixelType                     DefaultPixelType;
  typedef DataArray<DefaultPixelType>        DefaultArrayType;
#elif Anisotropy_BitDepth == 16
  typedef UInt16PixelType DefaultPixelType;
  typedef UInt16ArrayType DefaultArrayType;
#elif Anisotropy_BitDepth == 32
  typedef FloatPixelType DefaultPixelType;
  typedef FloatArrayType DefaultArrayType;
#else
  typedef UInt8PixelType DefaultPixelType;
  typedef UInt8ArrayType DefaultArrayType;
#endif

  //multicomponent pixels
  typedef itk::RGBPixel <uint8_t> RGBUInt8PixelType; //ipf color etc
  //typedef itk::RGBAPixel <float> RGBAFloatPixelType; //may be able to handle quats with this?

  //define dimensionality
  const unsigned int SliceDimension = 2;
  const unsigned int ImageDimension = 3;

  //slice directions
  const unsigned int XSlice = 0;
  const unsigned int YSlice = 1;
  const unsigned int ZSlice = 2;

  //define image types
  typedef itk::Image< DefaultPixelType, ImageDimension > DefaultImageType;
  typedef itk::Image< Int8PixelType, ImageDimension > Int8ImageType;
  typedef itk::Image< UInt8PixelType, ImageDimension > UInt8ImageType;
  typedef itk::Image< Int16PixelType, ImageDimension > Int16ImageType;
  typedef itk::Image< UInt16PixelType, ImageDimension > UInt16ImageType;
  typedef itk::Image< Int32PixelType, ImageDimension > Int32ImageType;
  typedef itk::Image< UInt32PixelType, ImageDimension > UInt32ImageType;
  typedef itk::Image< Int64PixelType, ImageDimension > Int64ImageType;
  typedef itk::Image< UInt64PixelType, ImageDimension > UInt64ImageType;
  typedef itk::Image< FloatPixelType, ImageDimension > FloatImageType;
  typedef itk::Image< DoublePixelType, ImageDimension > DoubleImageType;


  typedef itk::Image< DefaultPixelType, SliceDimension > DefaultSliceType;
  typedef itk::Image< Int8PixelType, SliceDimension > Int8SliceType;
  typedef itk::Image< UInt8PixelType, SliceDimension > UInt8SliceType;
  typedef itk::Image< Int16PixelType, SliceDimension > Int16SliceType;
  typedef itk::Image< UInt16PixelType, SliceDimension > UInt16SliceType;
  typedef itk::Image< Int32PixelType, SliceDimension > Int32SliceType;
  typedef itk::Image< UInt32PixelType, SliceDimension > UInt32SliceType;
  typedef itk::Image< Int64PixelType, SliceDimension > Int64SliceType;
  typedef itk::Image< UInt64PixelType, SliceDimension > UInt64SliceType;
  typedef itk::Image< FloatPixelType, SliceDimension > FloatSliceType;
  typedef itk::Image< DoublePixelType, SliceDimension > DoubleSliceType;

}
