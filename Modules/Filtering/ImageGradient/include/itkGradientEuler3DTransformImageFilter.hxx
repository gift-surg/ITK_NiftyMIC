/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: GradientEuler3DTransformImageFilter.hxx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) 2017, University College London
  Author:    Michael Ebner, michael.ebner.14@ucl.ac.uk

  Portions of this code are covered under the ITK copyright.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkGradientEuler3DTransformImageFilter_hxx
#define itkGradientEuler3DTransformImageFilter_hxx
#include "itkGradientEuler3DTransformImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkDerivativeOperator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"
#include "itkEuler3DTransform.h"

namespace itk
{
//
// Constructor
//
template< typename TInputImage, typename TOperatorValueType, typename TOutputValueType , typename TOutputImageType >
GradientEuler3DTransformImageFilter< TInputImage, TOperatorValueType, TOutputValueType, TOutputImageType >
::GradientEuler3DTransformImageFilter()
{
  //   "Transform" required ( not numbered )
  Self::AddRequiredInputName("Transform");
  Self::SetTransform(Euler3DTransform< TTransformPrecisionType >::New());
}

//
// Destructor
//
template< typename TInputImage, typename TOperatorValueType, typename TOutputValueType , typename TOutputImageType >
GradientEuler3DTransformImageFilter< TInputImage, TOperatorValueType, TOutputValueType, TOutputImageType >
::~GradientEuler3DTransformImageFilter()
{}

template< typename TInputImage, typename TOperatorValueType, typename TOutputValueType , typename TOutputImageType >
void
GradientEuler3DTransformImageFilter< TInputImage, TOperatorValueType, TOutputValueType, TOutputImageType >
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
GradientEuler3DTransformImageFilter< TInputImage, TOperatorValueType, TOutputValueType, TOutputImageType >
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
GradientEuler3DTransformImageFilter< TInputImage, TOperatorValueType, TOutputValueType, TOutputImageType >
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
GradientEuler3DTransformImageFilter< TInputImage, TOperatorValueType, TOutputValueType, TOutputImageType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Transform: " << this->GetTransform() << std::endl;

}
} // end namespace itk

#endif
