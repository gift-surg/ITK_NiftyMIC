/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkOrientedGaussianInterpolateImageFunction.hxx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkOrientedGaussianInterpolateImageFunction_hxx
#define itkOrientedGaussianInterpolateImageFunction_hxx

#include "itkOrientedGaussianInterpolateImageFunction.h"

#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

/**
 * Constructor
 */
template<typename TImageType, typename TCoordRep>
OrientedGaussianInterpolateImageFunction<TImageType, TCoordRep>
::OrientedGaussianInterpolateImageFunction()
{
  this->m_Alpha = 1.0;
  this->m_Sigma.Fill( 1.0 );
  this->m_Covariance.Fill( 0.0 );

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    this->m_Covariance[d*ImageDimension + d] = 1.0;
    }
}

/**
 * Standard "PrintSelf" method
 */
template<typename TImageType, typename TCoordRep>
void
OrientedGaussianInterpolateImageFunction<TImageType, TCoordRep>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Alpha: " << this->m_Alpha << std::endl;
  os << indent << "Sigma: " << this->m_Sigma << std::endl;
}

template<typename TImageType, typename TCoordRep>
void
OrientedGaussianInterpolateImageFunction<TImageType, TCoordRep>
::ComputeBoundingBox()
{
  if( !this->GetInputImage() )
    {
    return;
    }

  typename InputImageType::ConstPointer input = this->GetInputImage();
  typename InputImageType::SpacingType spacing = input->GetSpacing();
  typename InputImageType::SizeType size = input->GetBufferedRegion().GetSize();

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    this->m_BoundingBoxStart[d] = -0.5;
    this->m_BoundingBoxEnd[d] = static_cast<RealType>( size[d] ) - 0.5;
    this->m_CutoffDistance[d] = this->m_Sigma[d] * this->m_Alpha / spacing[d];
    }
}


// ME: Compute intensity value of one single voxel based on Gaussian interpolation
template<typename TImageType, typename TCoordRep>
typename OrientedGaussianInterpolateImageFunction<TImageType, TCoordRep>
::OutputType
OrientedGaussianInterpolateImageFunction<TImageType, TCoordRep>
::EvaluateAtContinuousIndex( const ContinuousIndexType & cindex, OutputType *grad ) const
{
  // ME: Where is ImageDimension defined?

  RealType sum_me = 0.0;
  RealType sum_m = 0.0;
  ArrayType dsum_me;
  ArrayType dsum_m;
  ArrayType dw;

  dsum_m.Fill( 0.0 );
  dsum_me.Fill( 0.0 );
  dw.Fill( 0.0 );

  // Loop over the voxels in the region identified
  // ME: First, compute region (rectangular bounding box) which is to be considered, i.e. have non-zero Gaussian weights
  ImageRegion<ImageDimension> region;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    int boundingBoxSize = static_cast<int>(
      this->m_BoundingBoxEnd[d] - this->m_BoundingBoxStart[d] + 0.5 );  // = size[d]
    int begin = vnl_math_max( 0, static_cast<int>( std::floor( cindex[d] -
      this->m_BoundingBoxStart[d] - this->m_CutoffDistance[d] ) ) );
    int end = vnl_math_min( boundingBoxSize, static_cast<int>( std::ceil(
      cindex[d] - this->m_BoundingBoxStart[d] + this->m_CutoffDistance[d] ) ) );
    region.SetIndex( d, begin );
    region.SetSize( d, end - begin );
    }

  // ME: Define iterator over chosen region
  ImageRegionConstIteratorWithIndex<InputImageType> It(
    this->GetInputImage(), region );

  // ME: Compute scaled oriented PSF for voxel space
  /* Scaling matrix */
  const typename InputImageType::SpacingType spacing = this->GetInputImage()->GetSpacing();
  itk::Matrix<double,ImageDimension,ImageDimension> S;
  S.Fill(0.0);
  for (int d = 0; d < ImageDimension; ++d)
    {
    S(d,d) = spacing[d];
    }
  // std::cout << "S = \n" << S << std::endl;

  /* Scale rotated inverse Gaussian needed for exponential function */
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

  RealType w = 0.0;

  // ME: For each voxel of that region do
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename InputImageType::IndexType index = It.GetIndex();

    w = this->ComputeExponentialFunction(index, cindex, CovScaledInv);
    RealType V = It.Get();  // ME: Intensity of current voxel
    sum_me += V * w;        // ME: Add Gaussian weighted intensity
    sum_m += w;             // ME: Record weight sum for subsequent normalization
    }

  RealType rc = sum_me / sum_m;   // ME: Final Gaussian interpolated voxel intensity

  return rc;
}


template<typename TImageType, typename TCoordRep>
typename OrientedGaussianInterpolateImageFunction<TImageType, TCoordRep>
::RealType
OrientedGaussianInterpolateImageFunction<TImageType, TCoordRep>
::ComputeExponentialFunction(
  IndexType point,
  ContinuousIndexType center,
  itk::Matrix<double, ImageDimension, ImageDimension> SigmaInverse ) const
{
  itk::Vector<double, ImageDimension> tmp;
  itk::Vector<double, ImageDimension> diff;
  RealType result;

  for (int i = 0; i < ImageDimension; ++i)
  {
    diff[i] = point[i] - center[i];
    // printf("diff[%d] = %f\n", i, diff[i]);
  }


  tmp = SigmaInverse*diff;
  result = diff*tmp;

  // std::cout << "Sigma*(point-center) = " << tmp << std::endl;
  // std::cout << "(point-center)'*Sigma*(point-center) = " << result << std::endl;

  result = exp( -0.5*result );
  // std::cout << "exp(.) = " << result << std::endl;

  return result;
}

} // namespace itk

#endif
