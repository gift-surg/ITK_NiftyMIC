/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkImageIteratorWithIndex.txx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) 2000 National Library of Medicine
  All rights reserved.

  See COPYRIGHT.txt for copyright details.

=========================================================================*/

// #include "itkImageIteratorWithIndex.h"

namespace itk
{



//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
template<class TPixel, unsigned int VImageDimension>
ImageIteratorWithIndex<TPixel, VImageDimension>
::ImageIteratorWithIndex()
{
  m_Position  = 0;
  m_Begin     = 0;
  m_End       = 0;
  m_Remaining = true;
  memset( m_Size, 0, VImageDimension*sizeof(unsigned long) );
}



//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
template<class TPixel, unsigned int VImageDimension>
ImageIteratorWithIndex<TPixel, VImageDimension>
::ImageIteratorWithIndex(const Self& it)
{
  m_Image = it.m_Image;     // copy the smart pointer

  m_PositionIndex     = it.m_PositionIndex;
  m_BeginIndex        = it.m_BeginIndex;
  m_EndIndex          = it.m_EndIndex;
  m_StartBufferIndex  = it.m_StartBufferIndex;

  memcpy(m_Size,        it.m_Size,        VImageDimension*sizeof(unsigned long));
  memcpy(m_OffsetTable, it.m_OffsetTable, (VImageDimension+1)*sizeof(unsigned long));
  
  m_Position    = it.m_Position;
  m_Begin       = it.m_Begin;
  m_End         = it.m_End;
  m_Remaining   = it.m_Remaining;
}



//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
template<class TPixel, unsigned int VImageDimension>
ImageIteratorWithIndex<TPixel, VImageDimension>
::ImageIteratorWithIndex(const SmartPointer<Image> &ptr,
              const Index &start,
              const unsigned long size[VImageDimension])
{
  m_Image = ptr;
  TPixel * m_Buffer   = m_Image->GetBufferPointer();
  m_StartBufferIndex  = m_Image->GetBufferStartIndex();
  m_BeginIndex        = start;
  m_PositionIndex     = m_BeginIndex;
  memcpy(m_Size,        size,             VImageDimension*sizeof(unsigned long));
  memcpy(m_OffsetTable, m_Image->GetOffsetTable(), (VImageDimension+1)*sizeof(unsigned long));

  // Compute the start position
  long offs =  m_Image->ComputeOffset( m_BeginIndex );
  m_Begin = m_Buffer + offs;
  
  // Compute the end offset
  Index pastEnd;
  for (unsigned int i=0; i < VImageDimension; ++i)
    {
    m_EndIndex[i] = m_BeginIndex[i] + m_Size[i];
    pastEnd[i]    = m_BeginIndex[i] + m_Size[i]-1;
    }
  m_End = m_Buffer + m_Image->ComputeOffset( pastEnd );
  m_End++;
  m_Remaining = true;

}
 

//----------------------------------------------------------------------
//    Assignment Operator
//----------------------------------------------------------------------
template<class TPixel, unsigned int VImageDimension>
ImageIteratorWithIndex<TPixel, VImageDimension> &
ImageIteratorWithIndex<TPixel, VImageDimension> 
::operator=(const Self& it)
{
  m_Image = it.m_Image;     // copy the smart pointer

  m_BeginIndex        = it.m_BeginIndex;
  m_EndIndex          = it.m_EndIndex;
  m_StartBufferIndex  = it.m_StartBufferIndex;
  m_PositionIndex     = it.m_PositionIndex;

  memcpy(m_Size,        it.m_Size,        VImageDimension*sizeof(unsigned long));
  memcpy(m_OffsetTable, it.m_OffsetTable, (VImageDimension+1)*sizeof(unsigned long));
  
  m_Position    = it.m_Position;
  m_Begin       = it.m_Begin;
  m_End         = it.m_End;
  m_Remaining   = it.m_Remaining;

  return *this;
} 
  



//----------------------------------------------------------------------------
// Begin() is the first pixel in the region.
//----------------------------------------------------------------------------
template<class TPixel, unsigned int VImageDimension>
void
ImageIteratorWithIndex<TPixel, VImageDimension>
::Begin()
{
  // Set the position at begin

  m_Position       = m_Begin;
  m_PositionIndex  = m_BeginIndex;
  m_Remaining      = true;
}


} // end namespace itk



