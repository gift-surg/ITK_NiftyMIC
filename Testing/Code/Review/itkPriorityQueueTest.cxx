/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkPriorityQueueTest.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <vnl/vnl_random.h>
#include "itkPriorityQueueContainer.h"

using namespace itk;

int itkPriorityQueueTest( int argc, char* argv[] )
{
  (void) argc;
  (void) argv;
  typedef MinPriorityQueueElementWrapper< int, double, int > MinPQElementType;
  typedef MaxPriorityQueueElementWrapper< int, double, int > MaxPQElementType;

  typedef PriorityQueueContainer< MinPQElementType, MinPQElementType, double, int >
    MinPQType;
  MinPQType::Pointer min_priority_queue = MinPQType::New( );

  typedef PriorityQueueContainer< MaxPQElementType, MaxPQElementType, double, int >
      MaxPQType;
  MaxPQType::Pointer max_priority_queue = MaxPQType::New( );

  std::list< double > sequence;
  sequence.push_back( -0.1 );
  sequence.push_back( 0.1 );
  sequence.push_back( 0.4 );
  sequence.push_back( -0.2 );
  sequence.push_back( -0.3 );
  sequence.push_back( 0.3 );
  sequence.push_back( 0.2 );
  sequence.push_back( 0.5 );
  sequence.push_back( -0.6 );
  sequence.push_back( -0.5 );
  sequence.push_back( 0.6 );
  sequence.push_back( 1. );
  sequence.push_back( -1. );

  std::list< double >::const_iterator it = sequence.begin();
  unsigned int i = 0;
  for(; it != sequence.end(); ++it, i++ )
    {
    min_priority_queue->Push( MinPQElementType( i, *it ) );
    max_priority_queue->Push( MaxPQElementType( i, *it ) );
    }

  sequence.sort();
  it = sequence.begin();
  i = sequence.size();

  std::cout <<"Min Priority Queue   ";
  while( !min_priority_queue->Empty() )
    {
    if( min_priority_queue->Peek().m_Priority != *it )
      {
      std::cout <<min_priority_queue->Peek().m_Priority <<" " <<*it <<std::endl;
      return EXIT_FAILURE;
      }
    if( min_priority_queue->Size() != i )
      {
      std::cout <<"Size " <<min_priority_queue->Size() <<" " <<i <<std::endl;
      return EXIT_FAILURE;
      }
    min_priority_queue->Pop();
    it++;
    i--;
    }
  std::cout <<"OK" <<std::endl;

  std::list< double >::const_reverse_iterator r_it = sequence.rbegin();
  i = sequence.size();

  std::cout <<"Max Priority Queue   ";
  while( !max_priority_queue->Empty() )
    {
    if( max_priority_queue->Peek().m_Priority != *r_it )
      {
      std::cout <<max_priority_queue->Peek().m_Priority <<" " <<*r_it <<std::endl;
      return EXIT_FAILURE;
      }
    if( max_priority_queue->Size() != i )
      {
      std::cout <<"Size " <<max_priority_queue->Size() <<" " <<i <<std::endl;
      return EXIT_FAILURE;
      }
    max_priority_queue->Pop();
    r_it++;
    i--;
    }
  std::cout <<"OK" <<std::endl;

  return EXIT_SUCCESS;
}
