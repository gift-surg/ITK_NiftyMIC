/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkProcrustesRegistrationMetric.txx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) 2000 National Library of Medicine
  All rights reserved.

  See COPYRIGHT.txt for copyright details.

=========================================================================*/
#ifndef _itkProcrustesRegistrationMetric_txx
#define _itkProcrustesRegistrationMetric_txx



#include <itkExceptionObject.h>

namespace itk
{

/**
 * Constructor
 */
template <class TTransform, unsigned int NDimension>
ProcrustesRegistrationMetric<TTransform, NDimension>
::ProcrustesRegistrationMetric()
{
}






/**
 * Compute Performs the evaluation of similarity
 */
template <class TTransform, unsigned int NDimension>
void
ProcrustesRegistrationMetric<TTransform, NDimension>
::Compute( void )
{

  if( m_Reference->Size() != m_Target->Size() )
  {
    ExceptionObject wrongSize;
    wrongSize.SetLocation( "ProcrustesRegistrationMetric Compute" );
    wrongSize.SetDescription( 
                "Reference and Target have different number of points");
    throw wrongSize;
  }
  
  m_MatchMeasure = vnl_vector<double>( m_Reference->Size() * NDimension );

  m_MatchMeasureDerivatives = 
                vnl_matrix<double>( 10, NDimension * m_Reference->Size() );
  
  typedef typename TargetType::Element        PointType;
  
  TargetType::ConstIterator     targetPoint    = m_Target->Begin();  
  ReferenceType::ConstIterator  referencePoint = m_Reference->Begin();

  typedef typename PointType::CoordRepType       CoordinatesType;

  vnl_vector<double>::iterator   similarityMeasure =
                                            m_MatchMeasure.begin();
  
  while(    targetPoint    != m_Target->End()  
         && referencePoint != m_Reference->End() ) 
  {

    PointType transformedPoint = 
                        m_Mapper->Transform( referencePoint.Value() );

    // TODO: This should be converted to using iterators through the
    // coordinates of the point.
    const CoordinatesType * transformed = transformedPoint.GetDataPointer();
    const CoordinatesType * target      = targetPoint.Value().GetDataPointer();

    for(unsigned int i=0; i<NDimension; i++)
    {
      *similarityMeasure = transformed[i] - target[i];
      similarityMeasure++;
    }

    targetPoint++;
    referencePoint++;

  }

}





} // end namespace itk

#endif
