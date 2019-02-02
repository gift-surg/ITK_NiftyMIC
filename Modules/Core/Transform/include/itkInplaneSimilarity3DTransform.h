/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: InplaneSimilarity3DTransform.h,v $
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
#ifndef itkInplaneSimilarity3DTransform_h
#define itkInplaneSimilarity3DTransform_h

#include <iostream>
#include "itkVersorRigid3DTransform.h"

namespace itk
{
/** \class InplaneSimilarity3DTransform
 * \brief InplaneSimilarity3DTransform of a vector space (e.g. space coordinates)
 *
 * This transform applies a rotation, translation to the space and isotropic
 * in-plane scaling of x and y voxel spacings. The idea is to use this to
 * extract the right in-plane scaling of the fixed image during registration.
 *
 * The rigid transform parameters and the scaling factor of the in-plane
 * spacing of the fixed image will be estimated so that it matches the moving
 * image.
 *
 * The direction matrix of the fixed image needs to be provided via the
 * SetFixedParameters method, i.e. appended after the center coordinates.
 *
 * The serialization of the optimizable parameters is an array of 7 elements.
 * The first 3 elements are the components of the versor representation
 * of 3D rotation. The next 3 parameters defines the translation in each
 * dimension. The last parameter defines the isotropic in-plane scaling.
 *
 * The serialization of the fixed parameters is an array of 3 elements defining
 * the center of rotation and another 9 elements defining the row-wise
 * representation of the fixed image direction matrix.
 *
 * After registration, update the fixed image spacing with
 * spacing[0] *= scale and spacing[1]*scale. In order to have an aligned fixed
 * image with the moving image, create a VersorRigid3DTransform object and
 * set versors coordinates according to the registration result and set the
 * translation to
 * R*D*Lambda*Dinv*(origin-center) + R*(center-origin) + translation + center
 * with R and center being the rotation matrix and center according to the
 * registration and D and origin the direction matrix and the origin of the
 * fixed image, respectively.
 *
 *
 * \sa VersorRigid3DTransform
 * \ingroup ITKTransform
 */
template<typename TParametersValueType=double>
class InplaneSimilarity3DTransform :
  public VersorRigid3DTransform<TParametersValueType>
{
public:
  /** Standard class typedefs. */
  typedef InplaneSimilarity3DTransform                 Self;
  typedef VersorRigid3DTransform<TParametersValueType> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(InplaneSimilarity3DTransform, VersorRigid3DTransform);

  /** Dimension of parameters. */
  itkStaticConstMacro(SpaceDimension, unsigned int, 3);
  itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(ParametersDimension, unsigned int, 7);

  /** Parameters Type   */
  typedef typename Superclass::ParametersType            ParametersType;
  typedef typename Superclass::FixedParametersType       FixedParametersType;
  typedef typename Superclass::JacobianType              JacobianType;
  typedef typename Superclass::ScalarType                ScalarType;
  typedef typename Superclass::InputPointType            InputPointType;
  typedef typename Superclass::OutputPointType           OutputPointType;
  typedef typename Superclass::InputVectorType           InputVectorType;
  typedef typename Superclass::OutputVectorType          OutputVectorType;
  typedef typename Superclass::InputVnlVectorType        InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType       OutputVnlVectorType;
  typedef typename Superclass::InputCovariantVectorType  InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType OutputCovariantVectorType;
  typedef typename Superclass::MatrixType                MatrixType;
  typedef typename Superclass::InverseMatrixType         InverseMatrixType;
  typedef typename Superclass::CenterType                CenterType;
  typedef typename Superclass::OffsetType                OffsetType;
  typedef typename Superclass::TranslationType           TranslationType;

  /** Versor type. */
  typedef typename Superclass::VersorType VersorType;
  typedef typename Superclass::AxisType   AxisType;
  typedef typename Superclass::AngleType  AngleType;
  typedef          TParametersValueType   ScaleType;

  /** Set the parameters to the IdentityTransform */
  virtual void SetIdentity(void) ITK_OVERRIDE;

  /** Directly set the rotation matrix of the transform.
   *
   * \warning The input matrix must be orthogonal with isotropic scaling
   * to within a specified tolerance, else an exception is thrown.
   *
   * \sa MatrixOffsetTransformBase::SetMatrix() */
  virtual void SetMatrix(const MatrixType & matrix) ITK_OVERRIDE;

  /** Directly set the rotation matrix of the transform.
   *
   * \warning The input matrix must be orthogonal with isotropic scaling
   * to within the specified tolerance, else an exception is thrown.
   *
   * \sa MatrixOffsetTransformBase::SetMatrix() */
  virtual void SetMatrix(const MatrixType & matrix, const TParametersValueType tolerance) ITK_OVERRIDE;

  /** Set the transformation from a container of parameters This is typically
   * used by optimizers.  There are 7 parameters. The first three represent the
   * versor, the next three represent the translation and the last one
   * represents the scaling factor. */
  // ME: Always set SetFixedParameters before SetParameters. Otherwise, no
  // correct update of ComputeMatrixParameters will be done!
  void SetParameters(const ParametersType & parameters) ITK_OVERRIDE;

  virtual const ParametersType & GetParameters(void) const ITK_OVERRIDE;

  /** Set the fixed parameters and update internal transformation. */
  // ME: Always set SetFixedParameters before SetParameters. Otherwise, no
  // correct update of ComputeMatrixParameters will be done!
  virtual void SetFixedParameters(const FixedParametersType &) ITK_OVERRIDE;

  /** Get the Fixed Parameters. */
  virtual const FixedParametersType & GetFixedParameters() const ITK_OVERRIDE;

  /** Set/Get the value of the isotropic scaling factor */
  void SetScale(ScaleType scale);

  itkGetConstReferenceMacro(Scale, ScaleType);

  /** This method computes the Jacobian matrix of the transformation.
   * given point or vector, returning the transformed point or
   * vector. The rank of the Jacobian will also indicate if the
   * transform is invertible at this point. */
  virtual void ComputeJacobianWithRespectToParameters( const InputPointType  & p, JacobianType & jacobian) const ITK_OVERRIDE;

protected:
  InplaneSimilarity3DTransform(const MatrixType & matrix, const OutputVectorType & offset);
  InplaneSimilarity3DTransform(unsigned int paramDim);
  InplaneSimilarity3DTransform();
  ~InplaneSimilarity3DTransform()
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  /** Recomputes the matrix by calling the Superclass::ComputeMatrix() and then
   * applying the scale factor. */
  void ComputeMatrix() ITK_OVERRIDE;

  /** Computes the parameters from an input matrix. */
  void ComputeMatrixParameters() ITK_OVERRIDE;

private:
  InplaneSimilarity3DTransform(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

  MatrixType GetMatrixTransform(const ScaleType &scale) const;

  ScaleType  m_Scale;
  MatrixType m_DirectionMatrix;
  MatrixType m_DirectionMatrixInverse;

}; // class InplaneSimilarity3DTransform
}  // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkInplaneSimilarity3DTransform.hxx"
#endif

#endif /* itkInplaneSimilarity3DTransform_h */
