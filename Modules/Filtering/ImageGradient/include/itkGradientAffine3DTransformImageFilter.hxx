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
#ifndef itkGradientAffine3DTransformImageFilter_hxx
#define itkGradientAffine3DTransformImageFilter_hxx
#include "itkGradientAffine3DTransformImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkDerivativeOperator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"
#include "itkAffineTransform.h"

namespace itk
{
//
// Constructor
//
template< typename TInputImage, typename TOperatorValueType, typename TOutputValueType , typename TOutputImageType >
GradientAffine3DTransformImageFilter< TInputImage, TOperatorValueType, TOutputValueType, TOutputImageType >
::GradientAffine3DTransformImageFilter()
{
  //   "Transform" required ( not numbered )
  Self::AddRequiredInputName("Transform");
  Self::SetTransform(AffineTransform< TTransformPrecisionType, InputImageDimension >::New());
}

//
// Destructor
//
template< typename TInputImage, typename TOperatorValueType, typename TOutputValueType , typename TOutputImageType >
GradientAffine3DTransformImageFilter< TInputImage, TOperatorValueType, TOutputValueType, TOutputImageType >
::~GradientAffine3DTransformImageFilter()
{}

template< typename TInputImage, typename TOperatorValueType, typename TOutputValueType , typename TOutputImageType >
void
GradientAffine3DTransformImageFilter< TInputImage, TOperatorValueType, TOutputValueType, TOutputImageType >
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

template< typename TInputImage, typename TOperatorValueType, typename TOutputValueType , typename TOutputImageType >
void
GradientAffine3DTransformImageFilter< TInputImage, TOperatorValueType, TOutputValueType, TOutputImageType >
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
  typedef ImageRegionIteratorWithIndex< OutputImageType > OutputIterator;
  OutputIterator outIt(outputPtr, outputRegionForThread);

  IndexType outputIndex;
  PointType outputPoint;         // Coordinates of current output pixel

  CovariantVectorType gradient;
  typename TransformType::JacobianType v;

  // Support for progress methods/callbacks
  ProgressReporter progress( this,
                             threadId,
                             outputRegionForThread.GetNumberOfPixels() );

  // Walk the output region
  outIt.GoToBegin();

  while ( !outIt.IsAtEnd() ) {
    outputIndex = outIt.GetIndex();
    outputPtr->TransformIndexToPhysicalPoint(outputIndex, outputPoint);

    // Compute Jacobian
    transformPtr->ComputeJacobianWithRespectToParameters(outputPoint, v);

    // gradient composed based on row-wise stacking
    for (unsigned int i = 0; i < InputImageDimension; ++i) {
      for (unsigned int j = 0; j < ParametersDimension; ++j) {
        gradient[i*ParametersDimension + j] = v.GetElement(i,j);
      }
    }
    outputPtr->SetPixel(outputIndex, gradient);

    progress.CompletedPixel();
    ++outIt;
  }
}

template< typename TInputImage, typename TOperatorValueType, typename TOutputValueType , typename TOutputImageType >
void
GradientAffine3DTransformImageFilter< TInputImage, TOperatorValueType, TOutputValueType, TOutputImageType >
::GenerateOutputInformation()
{
  // this methods is overloaded so that if the output image is a
  // VectorImage then the correct number of components are set.

  Superclass::GenerateOutputInformation();
  OutputImageType* output = this->GetOutput();

  if ( !output )
    {
    return;
    }
  if ( output->GetNumberOfComponentsPerPixel() != InputImageDimension )
    {
    output->SetNumberOfComponentsPerPixel( InputImageDimension );
    }
}


/**
 * Standard "PrintSelf" method
 */
template< typename TInputImage, typename TOperatorValueType, typename TOutputValueType , typename TOutputImageType >
void
GradientAffine3DTransformImageFilter< TInputImage, TOperatorValueType, TOutputValueType, TOutputImageType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Transform: " << this->GetTransform() << std::endl;

}
} // end namespace itk

#endif
