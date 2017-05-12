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
#ifndef itkGradientAffine3DTransformImageFilter_h
#define itkGradientAffine3DTransformImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkCovariantVector.h"
#include "itkImageRegionIterator.h"
#include "itkAffineTransform.h"
#include "itkDataObjectDecorator.h"

namespace itk
{


template <typename TPixelType, unsigned int VImageDimension > class VectorImage;


/** \class GradientAffine3DTransformImageFilter
 * \brief Computes the gradient of an image using directional derivatives.
 *
 * Computes the gradient of an image using directional derivatives.
 * The directional derivative at each pixel location is computed by
 * convolution with a first-order derivative operator.
 *
 * The second template parameter defines the value type used in the
 * derivative operator (defaults to float).  The third template
 * parameter defines the value type used for output image (defaults to
 * float).  The output image is defined as a covariant vector image
 * whose value type is specified as this third template parameter.
 *
 *
 * \sa Image
 * \sa Neighborhood
 * \sa NeighborhoodOperator
 * \sa NeighborhoodIterator
 *
 * \ingroup GradientFilters
 * \ingroup ITKImageGradient
 */
template< typename TInputImage,
          typename TOperatorValueType = float,
          typename TOutputValueType = float,
          typename TOutputImageType = Image< CovariantVector< TOutputValueType,
                                                           TInputImage::ImageDimension*AffineTransform<TOutputValueType, TInputImage::ImageDimension>::ParametersDimension >,
                                          TInputImage::ImageDimension > >
class GradientAffine3DTransformImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImageType >
{
public:
  /** Extract dimension from input image */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(ParametersDimension, unsigned int,
                      AffineTransform<TOutputValueType>::ParametersDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      InputImageDimension * ParametersDimension);

  /** Standard class typedefs. */
  typedef GradientAffine3DTransformImageFilter Self;

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                       InputImageType;
  typedef typename InputImageType::Pointer  InputImagePointer;
  typedef TOutputImageType                  OutputImageType;
  typedef typename OutputImageType::Pointer OutputImagePointer;

  /** Standard class typedefs. */
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GradientAffine3DTransformImageFilter, ImageToImageFilter);

  /** Image typedef support. */
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef TOperatorValueType                  OperatorValueType;
  typedef TOutputValueType                    OutputValueType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  typedef CovariantVector<
    OutputValueType, itkGetStaticConstMacro(OutputImageDimension) >
  CovariantVectorType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  /**
   *  Transform typedef.
   */
  typedef TOutputValueType                                  TTransformPrecisionType;
  typedef Transform< TTransformPrecisionType,
            itkGetStaticConstMacro(InputImageDimension),
            itkGetStaticConstMacro(InputImageDimension) >   TransformType;
  typedef typename TransformType::ConstPointer              TransformPointerType;
  typedef DataObjectDecorator<TransformType>                DecoratedTransformType;
  typedef typename DecoratedTransformType::Pointer          DecoratedTransformPointer;

  /** Image index typedef. */
  typedef typename OutputImageType::IndexType IndexType;

  /** Image point typedef. */
  typedef typename OutputImageType::PointType PointType;


  /** Get/Set the coordinate transformation.
   * Set the coordinate transform to use for resampling.  Note that this must
   * be in physical coordinates and it is the output-to-input transform, NOT
   * the input-to-output transform that you might naively expect.  By default
   * the filter uses an Identity transform. You must provide a different
   * transform here, before attempting to run the filter, if you do not want to
   * use the default Identity transform. */
   itkSetGetDecoratedObjectInputMacro(Transform, TransformType);

  /** GradientAffine3DTransformImageFilter needs a larger input requested region than
   * the output requested region.  As such, GradientAffine3DTransformImageFilter needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   *
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;


#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro( InputConvertibleToOutputCheck,
                   ( Concept::Convertible< InputPixelType, OutputValueType > ) );
  itkConceptMacro( OutputHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< OutputValueType > ) );
  // End concept checking
#endif


protected:
  GradientAffine3DTransformImageFilter();
  virtual ~GradientAffine3DTransformImageFilter();
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  /** GradientAffine3DTransformImageFilter can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData()
   * routine which is called for each processing thread. The output
   * image data is allocated automatically by the superclass prior to
   * calling ThreadedGenerateData().  ThreadedGenerateData can only
   * write to the portion of the output image specified by the
   * parameter "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                            ThreadIdType threadId) ITK_OVERRIDE;

private:
  GradientAffine3DTransformImageFilter(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

  virtual void GenerateOutputInformation() ITK_OVERRIDE;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGradientAffine3DTransformImageFilter.hxx"
#endif

#endif
