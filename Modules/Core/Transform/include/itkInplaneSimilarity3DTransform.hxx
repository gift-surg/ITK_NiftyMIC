/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: InplaneSimilarity3DTransform.hxx,v $
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
#ifndef itkInplaneSimilarity3DTransform_hxx
#define itkInplaneSimilarity3DTransform_hxx

#include "itkInplaneSimilarity3DTransform.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_det.h"

namespace itk
{
// Constructor with default arguments
template<typename TParametersValueType>
InplaneSimilarity3DTransform<TParametersValueType>
::InplaneSimilarity3DTransform() :
  Superclass(ParametersDimension),
  m_Scale(1.0)
{
  this->m_DirectionMatrix.SetIdentity();
  this->m_DirectionMatrixInverse.SetIdentity();
}

// Constructor with arguments
template<typename TParametersValueType>
InplaneSimilarity3DTransform<TParametersValueType>
::InplaneSimilarity3DTransform(unsigned int paramDim) :
  Superclass(paramDim),
  m_Scale(1.0)
{
  this->m_DirectionMatrix.SetIdentity();
  this->m_DirectionMatrixInverse.SetIdentity();
}

// Constructor with arguments
template<typename TParametersValueType>
InplaneSimilarity3DTransform<TParametersValueType>
::InplaneSimilarity3DTransform(const MatrixType & matrix, const OutputVectorType & offset) :
  Superclass(matrix, offset),
  m_Scale(1.0)
{
  this->m_DirectionMatrix.SetIdentity();
  this->m_DirectionMatrixInverse.SetIdentity();
}

// / Set the parameters to the IdentityTransform
template<typename TParametersValueType>
void
InplaneSimilarity3DTransform<TParametersValueType>
::SetIdentity(void)
{
  this->Superclass::SetIdentity();
  this->m_Scale = 1.0;
  this->m_DirectionMatrix.SetIdentity();
  this->m_DirectionMatrixInverse.SetIdentity();
}

// Set the scale factor
template<typename TParametersValueType>
void
InplaneSimilarity3DTransform<TParametersValueType>
::SetScale(ScaleType scale)
{
  m_Scale = scale;
  this->ComputeMatrix();
}

// Directly set the matrix
template<typename TParametersValueType>
void
InplaneSimilarity3DTransform<TParametersValueType>
::SetMatrix(const MatrixType & matrix)
{
  const TParametersValueType tolerance = MatrixOrthogonalityTolerance<TParametersValueType>::GetTolerance();
  this->SetMatrix( matrix, tolerance );
}

template<typename TParametersValueType>
void
InplaneSimilarity3DTransform<TParametersValueType>
::SetMatrix(const MatrixType & matrix, const TParametersValueType tolerance)
{
  //
  // Since the matrix should be an orthogonal matrix
  // multiplied by the scale factor, then its determinant
  // must be equal to the cube of the scale factor.
  //
  double det = vnl_det( matrix.GetVnlMatrix() );

  if( det == 0.0 )
    {
    itkExceptionMacro(<< "Attempting to set a matrix with a zero determinant");
    }

  //
  // A negative scale is not acceptable
  // It will imply a reflection of the coordinate system.
  //

  double s = sqrt(det);

  //
  // A negative scale is not acceptable
  // It will imply a reflection of the coordinate system.
  //
  if( s <= 0.0 )
    {
    itkExceptionMacro(<< "Attempting to set a matrix with a negative trace");
    }

  MatrixType testForOrthogonal = matrix;
  const MatrixType matrixTransform = this->GetMatrixTransform(s);

  testForOrthogonal = matrix * matrixTransform.GetInverse();

  if( !this->MatrixIsOrthogonal(testForOrthogonal, tolerance) )
    {
    itkExceptionMacro(<< "Attempting to set a non-orthogonal matrix (after removing scaling)");
    }

  typedef MatrixOffsetTransformBase<TParametersValueType, 3, 3> Baseclass;
  this->Baseclass::SetMatrix(matrix);
}


// Set Parameters
// ME: Always set SetFixedParameters before SetParameters. Otherwise, no
// correct update of ComputeMatrixParameters will be done!
template<typename TParametersValueType>
void
InplaneSimilarity3DTransform<TParametersValueType>
::SetParameters(const ParametersType & parameters)
{
  itkDebugMacro(<< "Setting parameters " << parameters);

  // Save parameters. Needed for proper operation of TransformUpdateParameters.
  if( &parameters != &(this->m_Parameters) )
    {
    this->m_Parameters = parameters;
    }

  // Transfer the versor part

  AxisType axis;

  double norm = parameters[0] * parameters[0];
  axis[0] = parameters[0];
  norm += parameters[1] * parameters[1];
  axis[1] = parameters[1];
  norm += parameters[2] * parameters[2];
  axis[2] = parameters[2];
  if( norm > 0 )
    {
    norm = std::sqrt(norm);
    }

  double epsilon = 1e-10;
  if( norm >= 1.0 - epsilon )
    {
    axis = axis / ( norm + epsilon * norm );
    }
  VersorType newVersor;
  newVersor.Set(axis);
  this->SetVarVersor(newVersor);
  m_Scale = parameters[6]; // must be set before calling ComputeMatrix();
  this->ComputeMatrix();

  itkDebugMacro( << "Versor is now " << this->GetVersor() );

  // Transfer the translation part
  TranslationType newTranslation;
  newTranslation[0] = parameters[3];
  newTranslation[1] = parameters[4];
  newTranslation[2] = parameters[5];
  this->SetVarTranslation(newTranslation);
  this->ComputeOffset();

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();

  // std::cout << "-----------------------------------------------------------" << std::endl;
  // std::cout << "\tParameters = " << this->m_Parameters << std::endl << std::endl;
  // std::cout << "\tDirectionMatrix = " << this->m_DirectionMatrix << std::endl;
  // std::cout << "\tDirectionMatrixInverse = " << this->m_DirectionMatrixInverse << std::endl;
  // std::cout << "\tMatrixTransform(lambda=" << this->m_Scale << ") = " << this->GetMatrixTransform(this->m_Scale) << std::endl;
  // std::cout << "\tMatrix = " << this->GetMatrix() << std::endl;

  itkDebugMacro(<< "After setting parameters ");
}


// ME: Always set SetFixedParameters before SetParameters. Otherwise, no
// correct update of ComputeMatrixParameters will be done!
template<typename TParametersValueType>
void
InplaneSimilarity3DTransform<TParametersValueType>
::SetFixedParameters(const FixedParametersType & fp)
{
  this->m_FixedParameters = fp;
  InputPointType c;
  for( unsigned int i = 0; i < InputSpaceDimension; i++ )
    {
    c[i] = this->m_FixedParameters[i];
    }
  this->SetCenter(c);

  for (unsigned int i = 0; i < OutputSpaceDimension; ++i)
  {
    for (unsigned int j = 0; j < InputSpaceDimension; ++j)
    {
      this->m_DirectionMatrix(i,j) = fp[InputSpaceDimension + i*InputSpaceDimension + j];
    }
  }
  this->m_DirectionMatrixInverse = this->m_DirectionMatrix.GetInverse();
  // std::cout << "DirectionMatrix = " << this->m_DirectionMatrix << std::endl;
}

//
// Get Parameters
//
// Parameters are ordered as:
//
// p[0:2] = right part of the versor (axis times std::sin(t/2))
// p[3:5} = translation components
// p[6:6} = scaling factor (in-plane)
//
template<typename TParametersValueType>
const typename InplaneSimilarity3DTransform<TParametersValueType>::ParametersType
& InplaneSimilarity3DTransform<TParametersValueType>
::GetParameters(void) const
  {
  itkDebugMacro(<< "Getting parameters ");

  this->m_Parameters[0] = this->GetVersor().GetX();
  this->m_Parameters[1] = this->GetVersor().GetY();
  this->m_Parameters[2] = this->GetVersor().GetZ();

  // Transfer the translation
  this->m_Parameters[3] = this->GetTranslation()[0];
  this->m_Parameters[4] = this->GetTranslation()[1];
  this->m_Parameters[5] = this->GetTranslation()[2];

  this->m_Parameters[6] = this->GetScale();

  itkDebugMacro(<< "After getting parameters " << this->m_Parameters);

  return this->m_Parameters;
  }


template<typename TParametersValueType>
const typename InplaneSimilarity3DTransform<TParametersValueType>::FixedParametersType
& InplaneSimilarity3DTransform<TParametersValueType>
::GetFixedParameters(void) const
{
  const InputPointType center = this->GetCenter();

  for( unsigned int i = 0; i < InputSpaceDimension; ++i )
    {
    this->m_FixedParameters[i] = center[i];
    }
  // std::cout<<this->m_FixedParameters << std::endl;
  return this->m_FixedParameters;
}

template<typename TParametersValueType>
void
InplaneSimilarity3DTransform<TParametersValueType>::ComputeJacobianWithRespectToParameters(const InputPointType & p,
                                                                           JacobianType & jacobian) const
{
  typedef typename VersorType::ValueType ValueType;

  // compute derivatives with respect to rotation
  const ValueType vx = this->GetVersor().GetX();
  const ValueType vy = this->GetVersor().GetY();
  const ValueType vz = this->GetVersor().GetZ();
  const ValueType vw = this->GetVersor().GetW();

  jacobian.SetSize( 3, this->GetNumberOfLocalParameters() );
  jacobian.Fill(0.0);

  const MatrixType matrixTransform = this->GetMatrixTransform(m_Scale);

  const InputVectorType pp = matrixTransform * (p - this->GetCenter());

  const double px = pp[0];
  const double py = pp[1];
  const double pz = pp[2];

  const double vxx = vx * vx;
  const double vyy = vy * vy;
  const double vzz = vz * vz;
  const double vww = vw * vw;

  const double vxy = vx * vy;
  const double vxz = vx * vz;
  const double vxw = vx * vw;

  const double vyz = vy * vz;
  const double vyw = vy * vw;
  const double vzw = vz * vw;

  // compute Jacobian with respect to quaternion parameters
  // (Same as superclass, i.e. VersorRigid3DTransform, but with different px, py and pz)
  jacobian[0][0] = 2.0 * ( ( vyw + vxz ) * py + ( vzw - vxy ) * pz )
    / vw;
  jacobian[1][0] = 2.0 * ( ( vyw - vxz ) * px   - 2 * vxw   * py + ( vxx - vww ) * pz )
    / vw;
  jacobian[2][0] = 2.0 * ( ( vzw + vxy ) * px + ( vww - vxx ) * py   - 2 * vxw   * pz )
    / vw;

  jacobian[0][1] = 2.0 * ( -2 * vyw  * px + ( vxw + vyz ) * py + ( vww - vyy ) * pz )
    / vw;
  jacobian[1][1] = 2.0 * ( ( vxw - vyz ) * px                + ( vzw + vxy ) * pz )
    / vw;
  jacobian[2][1] = 2.0 * ( ( vyy - vww ) * px + ( vzw - vxy ) * py   - 2 * vyw   * pz )
    / vw;

  jacobian[0][2] = 2.0 * ( -2 * vzw  * px + ( vzz - vww ) * py + ( vxw - vyz ) * pz )
    / vw;
  jacobian[1][2] = 2.0 * ( ( vww - vzz ) * px   - 2 * vzw   * py + ( vyw + vxz ) * pz )
    / vw;
  jacobian[2][2] = 2.0 * ( ( vxw + vyz ) * px + ( vyw - vxz ) * py )
    / vw;

  // compute Jacobian with respect to the translation parameters
  jacobian[0][3] = 1.0;
  jacobian[1][4] = 1.0;
  jacobian[2][5] = 1.0;

  // compute Jacobian with respect to the scale parameter
  const MatrixType & matrix = this->GetMatrix();
  MatrixType helper;
  for (int i = 0; i < OutputSpaceDimension; ++i)
  {
    for (int j = 0; j < InputSpaceDimension; ++j)
    {
      helper(i,j) = m_DirectionMatrix(i,0)*m_DirectionMatrixInverse(0,j)
                    + m_DirectionMatrix(i,1)*m_DirectionMatrixInverse(1,j);
    }
  }
  const InputVectorType jac_scale = matrix*helper*pp;

  jacobian[0][6] = jac_scale[0];
  jacobian[1][6] = jac_scale[1];
  jacobian[2][6] = jac_scale[2];
}

// Set the scale factor
template<typename TParametersValueType>
void
InplaneSimilarity3DTransform<TParametersValueType>
::ComputeMatrix()
{
  this->Superclass::ComputeMatrix();
  MatrixType newMatrix = this->GetMatrix();
  const MatrixType matrixTransform = this->GetMatrixTransform(m_Scale);
  newMatrix = newMatrix * matrixTransform;
  // std::cout<< "newMatrix = " << newMatrix << std::endl;

  this->SetVarMatrix(newMatrix);
}

/** Compute the matrix */
template<typename TParametersValueType>
void
InplaneSimilarity3DTransform<TParametersValueType>
::ComputeMatrixParameters(void)
{
  MatrixType matrix = this->GetMatrix();
  m_Scale = sqrt( vnl_det( matrix.GetVnlMatrix() ) );
  const MatrixType matrixTransform = this->GetMatrixTransform(m_Scale);

  matrix = matrix * matrixTransform.GetInverse();

  VersorType v;
  v.Set(matrix);
  this->SetVarVersor(v);
}

template<typename TParametersValueType>
typename InplaneSimilarity3DTransform<TParametersValueType>::MatrixType
InplaneSimilarity3DTransform<TParametersValueType>
::GetMatrixTransform(const ScaleType &scale) const
{
  MatrixType Lambda;
  Lambda.SetIdentity();
  Lambda(0,0) = scale;
  Lambda(1,1) = scale;

  // std::cout<< "Lambda = " << Lambda << std::endl;
  // std::cout<< "m_DirectionMatrix = " << m_DirectionMatrix << std::endl;
  // std::cout<< "m_DirectionMatrixInverse = " << m_DirectionMatrixInverse << std::endl;

  return m_DirectionMatrix * Lambda * m_DirectionMatrixInverse;
}


// Print self
template<typename TParametersValueType>
void
InplaneSimilarity3DTransform<TParametersValueType>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Scale (in-plane) = " << m_Scale << std::endl;
}

} // namespace

#endif
