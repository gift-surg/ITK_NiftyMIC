/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkOrientedGaussianInterpolateImageFilter_hxx
#define itkOrientedGaussianInterpolateImageFilter_hxx

#include "itkOrientedGaussianInterpolateImageFilter.h"
#include "itkObjectFactory.h"
#include "itkIdentityTransform.h"
#include "itkProgressReporter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageScanlineIterator.h"
#include "itkSpecialCoordinatesImage.h"
#include "itkDefaultConvertPixelTraits.h"

#include "itkAddImageFilter.h"

namespace itk
{
/**
 * Initialize new instance
 */
template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::OrientedGaussianInterpolateImageFilter()
{
  this->m_OutputOrigin.Fill(0.0);
  this->m_OutputSpacing.Fill(1.0);
  this->m_OutputDirection.SetIdentity();

  this->m_UseReferenceImage = false;

  this->m_Size.Fill(0);
  this->m_OutputStartIndex.Fill(0);

  // Pipeline input configuration

  // implicit:
  // #0 "Primary" required

  //  #1 "ReferenceImage" optional
  Self::AddRequiredInputName("ReferenceImage",1);
  Self::RemoveRequiredInputName("ReferenceImage");

  //   "Transform" required ( not numbered )
  Self::AddRequiredInputName("Transform");
  Self::SetTransform(IdentityTransform< TTransformPrecisionType, ImageDimension >::New());

  this->m_Interpolator = dynamic_cast< InterpolatorType * >
    ( LinearInterpolatorType::New().GetPointer() );

  this->m_Extrapolator = ITK_NULLPTR;

  this->m_DefaultPixelValue
    = NumericTraits<PixelType>::ZeroValue( m_DefaultPixelValue );

  this->m_Alpha = 1.0;
  this->m_Sigma.Fill( 1.0 );
  this->m_Covariance.Fill( 0.0 );

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    this->m_Covariance[d*ImageDimension + d] = 1.0;
    }
}

/**
 * Print out a description of self
 *
 * \todo Add details about this class
 */
template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
void
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "DefaultPixelValue: "
     << static_cast< typename NumericTraits< PixelType >::PrintType >
  ( m_DefaultPixelValue )
     << std::endl;
  os << indent << "Size: " << m_Size << std::endl;
  os << indent << "OutputStartIndex: " << m_OutputStartIndex << std::endl;
  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
  os << indent << "OutputDirection: " << m_OutputDirection << std::endl;
  os << indent << "Transform: " << this->GetTransform() << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "Extrapolator: " << m_Extrapolator.GetPointer() << std::endl;
  os << indent << "UseReferenceImage: " << ( m_UseReferenceImage ? "On" : "Off" )
     << std::endl;
}

template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
void
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::ComputeBoundingBox()
{
  if( !this->GetInput() )
    {
    return;
    }

  typename InputImageType::ConstPointer input = this->GetInput();
  typename InputImageType::SpacingType spacing = input->GetSpacing();
  typename InputImageRegionType::SizeType size = input->GetBufferedRegion().GetSize();

  // std::cout << "spacing = " << spacing << std::endl;
  // std::cout << "size = " << size << std::endl;

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    this->m_BoundingBoxStart[d] = -0.5;
    this->m_BoundingBoxEnd[d] = static_cast<RealType>( size[d] ) - 0.5;
    this->m_CutoffDistance[d] = this->m_Sigma[d] * this->m_Alpha / spacing[d];
    }
}

/**
 * Set the output image spacing.
 */
template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
void
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::SetOutputSpacing(const double *spacing)
{
  SpacingType s;
  for(unsigned int i = 0; i < TOutputImage::ImageDimension; ++i)
    {
    s[i] = static_cast< typename SpacingType::ValueType >(spacing[i]);
    }
  this->SetOutputSpacing(s);
}

/**
 * Set the output image origin.
 */
template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
void
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::SetOutputOrigin(const double *origin)
{
  OriginPointType p(origin);

  this->SetOutputOrigin(p);
}

/** Helper method to set the output parameters based on this image */
template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
void
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::SetOutputParametersFromImage(const ImageBaseType *image)
{
  this->SetOutputOrigin ( image->GetOrigin() );
  this->SetOutputSpacing ( image->GetSpacing() );
  this->SetOutputDirection ( image->GetDirection() );
  this->SetOutputStartIndex ( image->GetLargestPossibleRegion().GetIndex() );
  this->SetSize ( image->GetLargestPossibleRegion().GetSize() );
}

/**
 * Cast from interpolotor output to pixel type
 */
template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
typename OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::PixelType
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::CastPixelWithBoundsChecking(const InterpolatorOutputType value,
                              const ComponentType minComponent,
                              const ComponentType maxComponent ) const
{
  const unsigned int nComponents = InterpolatorConvertType::GetNumberOfComponents(value);
  PixelType          outputValue;

  NumericTraits<PixelType>::SetLength( outputValue, nComponents );

  for (unsigned int n=0; n<nComponents; n++)
    {
    ComponentType component = InterpolatorConvertType::GetNthComponent( n, value );

    if ( component < minComponent )
      {
      PixelConvertType::SetNthComponent( n, outputValue, static_cast<PixelComponentType>( minComponent ) );
      }
    else if ( component > maxComponent )
      {
      PixelConvertType::SetNthComponent( n, outputValue, static_cast<PixelComponentType>( maxComponent ) );
      }
    else
      {
      PixelConvertType::SetNthComponent(n, outputValue,
                                        static_cast<PixelComponentType>( component ) );
      }
    }

  return outputValue;
}


/**
 * GenerateData
 */
// template< typename TInputImage,
//           typename TOutputImage,
//           typename TInterpolatorPrecisionType,
//           typename TTransformPrecisionType >
// void
// OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
// ::GenerateData()
// {
//   // Call a method that can be overriden by a subclass to allocate
//   // memory for the filter's outputs
//   this->AllocateOutputs();

//   // Call a method that can be overridden by a subclass to perform
//   // some calculations prior to splitting the main computations into
//   // separate threads
//   this->BeforeThreadedGenerateData();

//   // Compute bounding box for Gaussian exponential
//   this->ComputeBoundingBox();

//   // Set up the multithreaded processing
//   ThreadStruct str;

//   // Get the output pointer
//   typename OutputImageType::Pointer outputPtr = this->GetOutput();

//   // Initialize output
//   outputPtr->FillBuffer( itk::NumericTraits< PixelType >::Zero );

//   // Create additional filter with swaped input and output
//   // Necessary for split input image region in this->ThreaderCallback
//   Pointer foo = OrientedGaussianInterpolateImageFilter<TOutputImage,TInputImage>::New();
//   foo->SetInput(this->GetOutput());
//   foo->SetOutputParametersFromImage(this->GetInput()); //seems not to be necessary
//   foo->GetOutput()->SetRequestedRegion(this->GetInput()->GetLargestPossibleRegion()); //necessary for SplitRequestedRegion

//   str.FilterFoo = foo; // To split input image region in this->ThreaderCallback
//   str.Filter = this;

//   // Get the input pointer and prepare multithreader
//   typename InputImageType::ConstPointer inputPtr = this->GetInput();

//   const ImageRegionSplitterBase * splitter = this->GetImageRegionSplitter();
//   const unsigned int validThreads = splitter->GetNumberOfSplits( inputPtr->GetRequestedRegion(), this->GetNumberOfThreads() );

//   this->GetMultiThreader()->SetNumberOfThreads( validThreads );
//   this->GetMultiThreader()->SetSingleMethod(this->ThreaderCallback, &str);

//   // Create valid threads amount of output images for threads working independently
//   m_OutputPtrThread = std::vector< OutputImagePointer>(validThreads);

//   for (int i = 0; i < validThreads; ++i)
//   {
//       m_OutputPtrThread[i] = OutputImageType::New();
//       m_OutputPtrThread[i]->CopyInformation(outputPtr);
//       m_OutputPtrThread[i]->SetBufferedRegion( outputPtr->GetRequestedRegion() );
//       m_OutputPtrThread[i]->Allocate();
//       m_OutputPtrThread[i]->FillBuffer( itk::NumericTraits< PixelType >::Zero );
//   }

//   // multithread the execution
//   this->GetMultiThreader()->SingleMethodExecute();

//   // Call a method that can be overridden by a subclass to perform
//   // some calculations after all the threads have completed
//   this->AfterThreadedGenerateData();

// }


/**
 * Set up state of filter before multi-threading.
 */
template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
void
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::BeforeThreadedGenerateData()
{
  this->ComputeBoundingBox();

  // std::cout << "m_BoundingBoxStart = " << this->m_BoundingBoxStart << std::endl;
  // std::cout << "m_BoundingBoxEnd = " << this->m_BoundingBoxEnd << std::endl;
  // std::cout << "m_CutoffDistance = " << this->m_CutoffDistance << std::endl;
}


/**
 * ThreadedGenerateData
 */
template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
void
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                       ThreadIdType threadId)
{

  // Get the output pointer
  typename OutputImageType::Pointer outputPtr = this->GetOutput();

  // Get the input pointer
  typename InputImageType::ConstPointer inputPtr = this->GetInput();

  // Get the input transform
  typename TransformType::ConstPointer transformPtr = this->GetTransform();

  // Create an iterator that will walk the input region for this thread.
  typedef ImageRegionIteratorWithIndex< TOutputImage > OutputIterator;
  OutputIterator outIt(outputPtr, outputRegionForThread);

  // Define a few indices that will be used to translate from an input pixel
  // to an output pixel
  PointType outputPoint;         // Coordinates of current output pixel
  PointType inputPoint;          // Coordinates of current input pixel

  ContinuousInputIndexType inputCIndex;
  ContinuousOutputIndexType outputCIndex;

  IndexType inputIndex;
  IndexType outputIndex;

  typedef typename InterpolatorType::OutputType OutputType;
  PixelType pixval;
  OutputType value;

  // Support for progress methods/callbacks
  ProgressReporter progress( this,
                             threadId,
                             outputRegionForThread.GetNumberOfPixels() );

  // Min/max values of the output pixel type AND these values
  // represented as the output type of the interpolator
  const PixelComponentType minValue =  NumericTraits< PixelComponentType >::NonpositiveMin();
  const PixelComponentType maxValue =  NumericTraits< PixelComponentType >::max();
  const ComponentType minOutputValue = static_cast< ComponentType >( minValue );
  const ComponentType maxOutputValue = static_cast< ComponentType >( maxValue );

  // ME: Compute scaled oriented PSF for voxel space
  // Scaling matrix
  const typename InputImageType::SpacingType spacing = inputPtr->GetSpacing();
  itk::Matrix<double,ImageDimension,ImageDimension> S;
  S.Fill(0.0);
  for (int d = 0; d < ImageDimension; ++d)
    {
    S(d,d) = spacing[d];
    }
  // std::cout << "S = \n" << S << std::endl;

  // Scale rotated inverse Gaussian needed for exponential function
  itk::Matrix<double, ImageDimension, ImageDimension> covariance;
  for (int i = 0; i < ImageDimension; i++ )
  {
    for (int j = 0; j < ImageDimension; j++ )
    {
      covariance(i,j) = this->m_Covariance[i*ImageDimension + j];
    }
  }

  itk::Matrix<double, ImageDimension, ImageDimension> CovScaledInv = S * covariance.GetInverse() * S;
  // std::cout << "CovScaledInv = \n" << CovScaledInv << std::endl;
  // std::cout << "m_CutoffDistance = " << m_CutoffDistance << std::endl;
  // std::cout << "m_BoundingBoxStart = " << m_BoundingBoxStart << std::endl;
  // std::cout << "m_BoundingBoxEnd = " << m_BoundingBoxEnd << std::endl;

  // Get region to check whether obtained index lies within
  InputImageRegionType entireInputRegion = inputPtr->GetBufferedRegion();

  // Walk the input region
  outIt.GoToBegin();

  while ( !outIt.IsAtEnd() )
    {
      RealType w = 0.0;
      RealType sum_m = 0.0;
      RealType sum_me = 0.0;

      // Determine the position of the first pixel in the scanline
      outputIndex = outIt.GetIndex();
      outputPtr->TransformIndexToPhysicalPoint(outputIndex, outputPoint);

      // Compute corresponding input pixel position
      inputPoint = transformPtr->TransformPoint(outputPoint);
      inputPtr->TransformPhysicalPointToContinuousIndex(inputPoint, inputCIndex);

      // Check that index is within input region
      if ( entireInputRegion.IsInside(inputCIndex) ){

        // Loop over the voxels in the region identified
        InputImageRegionType inputRegion;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          int boundingBoxSize = static_cast<int>(
            this->m_BoundingBoxEnd[d] - this->m_BoundingBoxStart[d] + 0.5 );  // = size[d]
          int begin = vnl_math_max( 0, static_cast<int>( std::floor( inputCIndex[d] -
            this->m_BoundingBoxStart[d] - this->m_CutoffDistance[d] ) ) );
          int end = vnl_math_min( boundingBoxSize, static_cast<int>( std::ceil(
            inputCIndex[d] - this->m_BoundingBoxStart[d] + this->m_CutoffDistance[d] ) ) );

          inputRegion.SetIndex( d, begin );
          inputRegion.SetSize( d, end - begin );

          // std::cout << "begin = " << begin << ", end = " << end << std::endl;
          // std::cout << "boundingBoxSize = " << boundingBoxSize << std::endl;
          }

        // std::cout << "inputRegion = " << inputRegion << std::endl;

        // ME: Define iterator over chosen region
        ImageRegionConstIteratorWithIndex<InputImageType> inIt( inputPtr, inputRegion );

        // ME: For each voxel of that region do
        for( inIt.GoToBegin(); !inIt.IsAtEnd(); ++inIt )
          {
          typename InputImageType::IndexType index = inIt.GetIndex();
          // std::cout << index;

          w = this->ComputeExponentialFunction(index, inputCIndex, CovScaledInv);
          RealType V = inIt.Get();  // ME: Intensity of current voxel
          sum_me += V * w;        // ME: Add Gaussian weighted intensity
          sum_m += w;             // ME: Record weight sum for subsequent normalization
          }
        RealType rc = sum_me / sum_m;   // ME: Final Gaussian interpolated voxel intensity

        outIt.Set(rc);

      }
      else{
        outIt.Set( this->m_DefaultPixelValue );
      }
      progress.CompletedPixel();
      ++outIt;
    }

  // // DEBUGGING OUTPUT
  // std::cout << "sizeInput = " << inputRegionSize << std::endl;
  // std::cout << "sizeOutput = " << outputRegionSize << std::endl;
  // std::cout << "inputIndex = " << inputIndex << std::endl;
  // std::cout << "inputPoint = " << inputPoint << std::endl;
  // std::cout << "outputPoint = " << outputPoint << std::endl;
  // std::cout << "outputIndex = " << outputIndex << std::endl;
  // std::cout << "outputCIndex = " << outputCIndex << std::endl;

  // ImageAlgorithm::Copy(inputPtr, outputPtr, outputPtr->GetRequestedRegion(),
  //                      outputPtr->GetRequestedRegion() );

}


template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
typename OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::RealType
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::ComputeExponentialFunction(
  IndexType point,
  ContinuousOutputIndexType center,
  itk::Matrix<double, ImageDimension, ImageDimension> CovScaledInv ) const
{
  itk::Vector<double, ImageDimension> tmp;
  itk::Vector<double, ImageDimension> diff;
  RealType result;

  for (int i = 0; i < ImageDimension; ++i)
  {
    diff[i] = point[i] - center[i];
    // printf("diff[%d] = %f\n", i, diff[i]);
  }


  tmp = CovScaledInv*diff;
  result = diff*tmp;

  // std::cout << "Sigma*(point-center) = " << tmp << std::endl;
  // std::cout << "(point-center)'*Sigma*(point-center) = " << result << std::endl;

  result = exp( -0.5*result );
  // std::cout << "exp(.) = " << result << std::endl;

  return result;
}

/**
 * Inform pipeline of necessary input image region
 *
 * Determining the actual input region is non-trivial, especially
 * when we cannot assume anything about the transform being used.
 * So we do the easy thing and request the entire input image.
 */
template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
void
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::GenerateInputRequestedRegion()
{
  // call the superclass's implementation of this method
  Superclass::GenerateInputRequestedRegion();

  if ( !this->GetInput() )
    {
    return;
    }

  // get pointers to the input and output
  InputImagePointer inputPtr  =
    const_cast< TInputImage * >( this->GetInput() );

  // Request the entire input image
  inputPtr->SetRequestedRegionToLargestPossibleRegion();
}

/**
 * Inform pipeline of required output region
 */
template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
void
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::GenerateOutputInformation()
{
  // call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  // get pointers to the input and output
  OutputImageType *outputPtr = this->GetOutput();
  if ( !outputPtr )
    {
    return;
    }

  const ReferenceImageBaseType *referenceImage = this->GetReferenceImage();

  // Set the size of the output region
  if ( m_UseReferenceImage && referenceImage )
    {
    outputPtr->SetLargestPossibleRegion(
      referenceImage->GetLargestPossibleRegion() );
    }
  else
    {
    typename TOutputImage::RegionType outputLargestPossibleRegion;
    outputLargestPossibleRegion.SetSize(m_Size);
    outputLargestPossibleRegion.SetIndex(m_OutputStartIndex);
    outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
    }

  // Set spacing and origin
  if ( m_UseReferenceImage && referenceImage )
    {
    outputPtr->SetSpacing( referenceImage->GetSpacing() );
    outputPtr->SetOrigin( referenceImage->GetOrigin() );
    outputPtr->SetDirection( referenceImage->GetDirection() );
    }
  else
    {
    outputPtr->SetSpacing(m_OutputSpacing);
    outputPtr->SetOrigin(m_OutputOrigin);
    outputPtr->SetDirection(m_OutputDirection);
    }
}

/**
 * Verify if any of the components has been modified.
 */
template< typename TInputImage,
          typename TOutputImage,
          typename TInterpolatorPrecisionType,
          typename TTransformPrecisionType >
ModifiedTimeType
OrientedGaussianInterpolateImageFilter< TInputImage, TOutputImage, TInterpolatorPrecisionType, TTransformPrecisionType >
::GetMTime(void) const
{
  ModifiedTimeType latestTime = Object::GetMTime();

  if ( m_Interpolator )
    {
    if ( latestTime < m_Interpolator->GetMTime() )
      {
      latestTime = m_Interpolator->GetMTime();
      }
    }

  return latestTime;
}

} // end namespace itk

#endif
