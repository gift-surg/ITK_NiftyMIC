/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkVectorFuzzyConnectednessImageFilter.txx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) 2000 National Library of Medicine
  All rights reserved.

  See COPYRIGHT.txt for copyright details.

=========================================================================*/
#ifndef _itkVectorFuzzyConnectednessImageFilter_txx
#define _itkVectorFuzzyConnectednessImageFilter_txx
#include "itkSimpleImageRegionIterator.h"


namespace itk{


/**
 *
 */
template <class TInputImage, class TOutputImage>
VectorFuzzyConnectednessImageFilter<TInputImage,TOutputImage>
::VectorFuzzyConnectednessImageFilter()
{
}

/**
 *
 */
template <class TInputImage, class TOutputImage>
VectorFuzzyConnectednessImageFilter<TInputImage,TOutputImage>
::~VectorFuzzyConnectednessImageFilter()
{
}


/**
 *
 */
template <class TInputImage, class TOutputImage>
void
VectorFuzzyConnectednessImageFilter<TInputImage,TOutputImage>
::SetObjectsMatrix(const MatrixType object_max,const int object_num)
{
	m_ObjectsCovMatrix[object_num] = object_max;
}


/**
 *
 */
template <class TInputImage, class TOutputImage>
void 
VectorFuzzyConnectednessImageFilter<TInputImage,TOutputImage>
::SetObjectsMean(const IntVector mean, const int object_num)
{
	m_ObjectsMean[object_num] = mean;
}

/**
 *
 */
template <class TInputImage, class TOutputImage>
void 
VectorFuzzyConnectednessImageFilter<TInputImage,TOutputImage>
::SetObjectsSeed(const IndexType &seed, const int object_num)
{
	m_ObjectsSeed[object_num].push_front(seed);
}


/**
 *
 */
template <class TInputImage, class TOutputImage>
void 
VectorFuzzyConnectednessImageFilter<TInputImage,TOutputImage>
::Initialization()
{

	m_SpherePointsNum = new int[8+1];
	m_SpherePointsLoc = new OffsetType[8+1];

	m_ObjectsMean = new IntVector[m_Objects];
	m_ObjectsSeed = new ListType[m_Objects];
	m_ObjectsCovMatrix = new MatrixType[m_Objects];
	m_ObjectsMap = new FloatType[m_Objects];
	m_ObjectsMaxDiff = new IntVector[m_Objects];

}


/**
 *
 */
template <class TInputImage, class TOutputImage>
void 
VectorFuzzyConnectednessImageFilter<TInputImage,TOutputImage>
::ScalePrepare()
{

	int i,j,k,l,tti1,tti2;
	double tt1;
	double anisotropy_row,anisotropy_col,anisotropy_slice;
	double const *spaceing;
	int ppptti1[2*(8+5)][2*(8+5)][2*(8+5)];
	IntVector location;


	m_InputImage = this->GetInput();
		
	spaceing = m_InputImage->GetSpacing();
	m_Size = m_InputImage->GetLargestPossibleRegion().GetSize();


	anisotropy_col = spaceing[0];
	anisotropy_row = spaceing[1];
	anisotropy_slice = spaceing[2];

	tt1 = anisotropy_col;
	if(tt1>anisotropy_row)
		tt1 = anisotropy_row;
		
	if(m_Size[2] > 1 && tt1 > anisotropy_slice)
		tt1 = anisotropy_slice;
	
	anisotropy_col = anisotropy_col/tt1;
	anisotropy_row = anisotropy_row/tt1;
	
	if(m_Size[2] > 1)
		anisotropy_slice = anisotropy_slice/tt1;

	for(i=0;i<2*(8+5);i++)
		for(j=0;j<2*(8+5);j++)
			for(l=0;l<2*(8+5);l++)
				ppptti1[i][j][l] = 0;


	for(int i = 0;i<=8;i++)
		m_SpherePointsNum[i] = 0;

  tti1 = 8 + 5;
  if (m_Size[2] > 1)
	{
	  for (k = 0; k <= 8; k++)
		{
		  for (i = -k - 2; i <= k + 2; i++)
			for (j = -k - 2; j <= k + 2; j++)
			  for (l = -k - 2; l <= k + 2; l++)
				if (ppptti1[tti1 + i][tti1 + j][tti1 + l] == 0)
				  {
					tt1 = sqrt(pow(((double) i) * anisotropy_slice, 2.0) +
							pow(((double) j) * anisotropy_row,2.0) + pow(((double) l) * anisotropy_col, 2.0));
					if (tt1 <= ((double) k) + 0.5)
					  {
						m_SpherePointsNum[k] = m_SpherePointsNum[k] + 1;
						ppptti1[tti1 + i][tti1 + j][tti1 + l] = 2;
					  }
				  }

		m_SpherePointsLoc[k].resize(m_SpherePointsNum[k]);

		tti2 = 0;
		for (i = -k - 2; i <= k + 2; i++)
			for (j = -k - 2; j <= k + 2; j++)
			  for (l = -k - 2; l <= k + 2; l++)
				if (ppptti1[(8+5) + i][(8+5) + j][(8+5) + l] == 2)
				  {
					ppptti1[(8+5) + i][(8+5) + j][(8+5) + l] = 1;
					location = l,j,i;
					m_SpherePointsLoc[k][tti2] = location;
					tti2 = tti2+1;
				  }
		}
	}
  else
	{
	  for (k = 0; k <= 8; k++)
		{
		  for (j = -k - 2; j <= k + 2; j++)
			for (l = -k - 2; l <= k + 2; l++)
			  if (ppptti1[tti1][tti1 + j][tti1 + l] == 0)
				{
				  tt1 = sqrt(pow(((double) j) * anisotropy_row, 2.0)
							 + pow(((double) l) * anisotropy_col, 2.0));
				  if (tt1 <= ((double) k) + 0.5)
					{
					  m_SpherePointsNum[k] = m_SpherePointsNum[k] + 1;
					  ppptti1[tti1][tti1 + j][tti1 + l] = 2;
					}
				}

	      m_SpherePointsLoc[k].resize(m_SpherePointsNum[k]);

		  tti2 = 0;
		  for (j = -k - 2; j <= k + 2; j++)
			for (l = -k - 2; l <= k + 2; l++)
			  if (ppptti1[tti1][tti1 + j][tti1 + l] == 2)
				{
				  ppptti1[tti1][tti1 + j][tti1 + l] = 1;
				  location = l,j,0;
				  m_SpherePointsLoc[k][tti2] = location;
				  tti2 = tti2 + 1;
				}
		}
	}

	m_MaskTotal = 0.0;
	for(i = -1;i <= 1;i++)
		for(j = -1; j <= 1;j++)
			for(k = -1;k <= 1;k++)
				m_Mask[i+1][j+1][k+1] = 0.0;
	if(m_Size[2] == 1)
	{
		for(i = 0;i<=0;i++)
			for(j = -1; j <= 1;j++)
				for(k = -1;k <= 1;k++)

				{
					tt1 = pow(anisotropy_col * k, 2.0);
					tt1 = tt1 + pow(anisotropy_row * j, 2.0);
					tt1 = 1 / (1 + tt1);
					m_Mask[i + 1][j + 1][k + 1] = tt1;
					m_MaskTotal = m_MaskTotal + tt1;
				 }
	}
	else
	{
		for(i = -1;i<=1;i++)
			for(j = -1; j <= 1;j++)
				for(k = -1;k <= 1;k++)

				{
					tt1 = pow(anisotropy_col * k, 2.0);
					tt1 = tt1 + pow(anisotropy_row * j, 2.0);
					tt1 = tt1 + pow(anisotropy_row * i, 2.0);
					tt1 = 1 / (1 + tt1);
					m_Mask[i + 1][j + 1][k + 1] = tt1;
					m_MaskTotal = m_MaskTotal + tt1;
				 }
	}

}


template <class TInputImage, class TOutputImage>
void 
VectorFuzzyConnectednessImageFilter<TInputImage,TOutputImage>
::Compute_LookupTable()
{

	const float HistThreshold = 0.90;
	MatrixType      HomogeneityCovarianceMatrix;	
	itk::Vector<double,3>    vectorA,vectorB;

	typedef std::vector<int>        VectorInt;

	VectorInt  *Histogram;
	IntVector            Hist_sum;
	int i,j,k,l,pslices,prow,pcol,tti1;
	IndexType            index1,index2;
	IntVector            value1,value2;
	float     result;   

	Histogram = new VectorInt[3];

	m_InputImage = this->GetInput();
	m_Size = m_InputImage->GetLargestPossibleRegion().GetSize();

	pslices = m_Size[2];
	prow = m_Size[1];
	pcol = m_Size[0];
	
	for(i = 0;i<3;i++)
		m_HomoMaxDiff[i] = 0;

	for(i = 0;i<pslices ;i++)
		for(j = 0;j<prow;j++)
			for(k = 0;k<pcol-1;k++)
			{
				index1[2] = index2[2] = i;
				index1[1] = index2[1] = j;
				index1[0] = k;
				index2[0] = k+1;
				value1 =  m_InputImage->GetPixel(index1);
				value2 =  m_InputImage->GetPixel(index2);

				for(l = 0;l<3;l++)
					if(abs(value1[l] - value2[l]) > m_HomoMaxDiff[l])
						m_HomoMaxDiff[l] = abs(value1[l] - value2[l]);
			}

	for(i = 0;i<pslices ;i++)
		for(j = 0;j<prow-1;j++)
			for(k = 0;k<pcol;k++)
			{
				index1[2] = index2[2] = i;
				index1[1] = j;
				index2[1] = j+1;
				index1[0] = index2[0] = k;

				value1 =  m_InputImage->GetPixel(index1);
				value2 =  m_InputImage->GetPixel(index2);

				for(l = 0;l<3;l++)
					if(abs(value1[l] - value2[l]) > m_HomoMaxDiff[l])
						m_HomoMaxDiff[l] = abs(value1[l] - value2[l]);
			}

	for(i = 0;i<pslices-1;i++)
		for(j = 0;j<prow;j++)
			for(k = 0;k<pcol;k++)
			{
				index1[2] = i;
				index2[2] = i+1;
				index1[1] = index2[1] = j;
				index1[0] = index2[0] = k;

				value1 =  m_InputImage->GetPixel(index1);
				value2 =  m_InputImage->GetPixel(index2);

				for(l = 0;l<3;l++)
					if(abs(value1[l] - value2[l]) > m_HomoMaxDiff[l])
						m_HomoMaxDiff[l] = abs(value1[l] - value2[l]);
			}

	for(i = 0;i<3;i++)
		Histogram[i].resize(m_HomoMaxDiff[i]+1);

	for(i = 0;i<3;i++)
		for(j = 0;j<=m_HomoMaxDiff[i];j++)
			Histogram[i][j] = 0;

	for(i = 0;i<pslices ;i++)
		for(j = 0;j<prow;j++)
			for(k = 0;k<pcol-1;k++)
			{
				index1[2] = index2[2] = i;
				index1[1] = index2[1] = j;
				index1[0] = k;
				index2[0] = k+1;
				value1 =  m_InputImage->GetPixel(index1);
				value2 =  m_InputImage->GetPixel(index2);
									
				for(l = 0;l<3;l++)
				{
					tti1 = abs(value1[l] - value2[l]);
					Histogram[l][tti1] = Histogram[l][tti1] +1;
				}
			}

	for(i = 0;i<pslices ;i++)
		for(j = 0;j<prow-1;j++)
			for(k = 0;k<pcol;k++)
			{
				index1[2] = index2[2] = i;
				index1[1] = j;
				index2[1] = j+1;
				index1[0] = index2[0] = k;

				value1 =  m_InputImage->GetPixel(index1);
				value2 =  m_InputImage->GetPixel(index2);

				for(l = 0;l<3;l++)
				{
					tti1 = abs(value1[l] - value2[l]);
					Histogram[l][tti1] = Histogram[l][tti1] +1;
				}
			}

	for(i = 0;i<pslices-1;i++)
		for(j = 0;j<prow;j++)
			for(k = 0;k<pcol;k++)
			{
				index1[2] = i;
				index2[2] = i+1;
				index1[1] = index2[1] = j;
				index1[0] = index2[0] = k;

				value1 =  m_InputImage->GetPixel(index1);
				value2 =  m_InputImage->GetPixel(index2);

				for(l = 0;l<3;l++)
				{
					tti1 = abs(value1[l] - value2[l]);
					Histogram[l][tti1] = Histogram[l][tti1] +1;
				}
			}

	for( i =0;i<3; i++)
	{
		Hist_sum[i] = 0;
		for(j = 0;j<=m_HomoMaxDiff[i];j++)
			Hist_sum[i] = Hist_sum[i] + Histogram[i][j];
	}


     for(i = 0;i<3;i++)
	 {
	  for(j=0;j<=m_HomoMaxDiff[i];j++)
	    {
	      tti1 = 0;
	      m_FeaturesThreshold[i] = (double)j;
	      for(k=0;k<=j;k++)
			tti1 = tti1+Histogram[i][k];
	      if (((double)tti1 /(double) Hist_sum[i])>=HistThreshold)
			break;
	    }
	  }
    for(i = 0;i<3;i++)
	{
	  tti1 = 1;
	  for(j=0;j<=i-1;j++)
	    tti1 = tti1*m_FeaturesThreshold[j]; 
	  m_PowerValue[i] = tti1;
	}

  /*    To computer the homogeneity covariance matrix     */

  for (int x = 0;x < 3; x++)
    for(int y = x;y < 3; y++)
	{
		HomogeneityCovarianceMatrix[x][y] = 0;
		double tt1 = 0;
		double tt2 = 0;
		int count = 0;
		for(i = 0;i<pslices ;i++)
			for(j = 0;j<prow;j++)
				for(k = 0;k<pcol-1;k++)
				{
					index1[2] = index2[2] = i;
					index1[1] = index2[1] = j;
					index1[0] = k;
					index2[0] = k+1;
					value1 =  m_InputImage->GetPixel(index1);
					value2 =  m_InputImage->GetPixel(index2);
					tt1 = abs(value1[x]-value2[x]) ;
					tt2 = abs(value1[y]-value2[y]) ;
					if((tt1<m_FeaturesThreshold[x]) && (tt2<m_FeaturesThreshold[y]))
					{
						HomogeneityCovarianceMatrix[x][y] += tt1*tt2;	
						count++;
					}			
				}
		for(i = 0;i<pslices ;i++)
			for(j = 0;j<prow-1;j++)
				for(k = 0;k<pcol;k++)
				{
					index1[2] = index2[2] = i;
					index1[1] = j;
					index2[1] = j+1;
					index1[0] = index2[0] = k;

					value1 =  m_InputImage->GetPixel(index1);
					value2 =  m_InputImage->GetPixel(index2);
					tt1 = abs(value1[x]-value2[x]);
					tt2 = abs(value1[y]-value2[y]);

					if(tt1<m_FeaturesThreshold[x] && tt2<m_FeaturesThreshold[y])
					{
						HomogeneityCovarianceMatrix[x][y] += tt1*tt2;	
						count++;
					}			

				}

		for(i = 0;i<pslices-1;i++)
			for(j = 0;j<prow;j++)
				for(k = 0;k<pcol;k++)
				{
					index1[2] = i;
					index2[2] = i+1;
					index1[1] = index2[1] = j;
					index1[0] = index2[0] = k;

					value1 =  m_InputImage->GetPixel(index1);
					value2 =  m_InputImage->GetPixel(index2);
					tt1 = abs(value1[x]-value2[x]);
					tt2 = abs(value1[y]-value2[y]);

					if(tt1<m_FeaturesThreshold[x] && tt2<m_FeaturesThreshold[y])
					{
						HomogeneityCovarianceMatrix[x][y] += tt1*tt2;	
						count++;
					}			

				}
		HomogeneityCovarianceMatrix[x][y] = HomogeneityCovarianceMatrix[x][y] /(double)count;
	}

	for (int x = 0;x < 3; x++)
	  for(int y = 0;y < x; y++)
	    HomogeneityCovarianceMatrix[x][y] = HomogeneityCovarianceMatrix[y][x];
	  
	
	vnl_matrix<double>  HomoInverseMatrix = HomogeneityCovarianceMatrix.GetInverse();	  

	tti1 = 1;
	for(int x=0;x<3;x++)
		tti1 = tti1*m_FeaturesThreshold[x];

	m_HomogeneityMap.resize(tti1);
	m_ScaleMap.resize(tti1);
		
	for(i = 0;i<tti1;i++)
	{
		k = i;
		for(j=0;j<3;j++)
	    {
		  vectorA[j] = (k % m_FeaturesThreshold[j]);
	      k = k/m_FeaturesThreshold[j];
	    }
		
		for(int x = 0;x<3;x++)
		{
			vectorB[x] = 0;
			for(int y = 0;y<3;y++)
				vectorB[x] = vectorB[x] + vectorA[y]*HomoInverseMatrix[x][y];
		}
		result = 0;
		for(int x = 0;x<3;x++)
			result = result + vectorB[x]*vectorA[x];

		m_HomogeneityMap[i] = (float) exp(-0.5*result);
		m_ScaleMap[i] = (float) exp(-0.5*result/9.0);
	}


	/*     compute object-based affinity look-up table  */

	for(l = 0;l<m_Objects;l++)
		for(int x = 0;x<3;x++)
		{
			m_ObjectsMaxDiff[l][x] = 0;
			for(i = 0;i<pslices;i++)
				for(j = 0;j<prow;j++)
					for(k = 0;k<pcol;k++)
					{
						index1[2] = i;
						index1[1] = j;
						index1[0] = k;
						value1 = m_InputImage->GetPixel(index1);
						if(abs(value1[x] - m_ObjectsMean[l][x])> m_ObjectsMaxDiff[l][x])
							m_ObjectsMaxDiff[l][x] = abs(value1[x] - m_ObjectsMean[l][x]);
					}
		}
	
	for(l = 0;l<m_Objects;l++)
		for(int x = 0;x<3;x++)
		 m_ObjectsMaxDiff[l][x] = (int) m_ObjectsMaxDiff[l][x]/3;

	for(l = 0;l<m_Objects;l++)
	{
		tti1 = 1;
		for(int x = 0;x<3;x++)
			tti1 = tti1*(m_ObjectsMaxDiff[l][x]+1);
		m_ObjectsMap[l].resize(tti1);
		 
        vnl_matrix<double> ObjectInverseMatrix = m_ObjectsCovMatrix[l].GetInverse();
		for(i = 0;i<tti1;i++)
	    {
			k = i;
			for(j=0;j<3;j++)
			{
  			     vectorA[j] = (k % (m_ObjectsMaxDiff[l][j]+1));
		         k = k/(m_ObjectsMaxDiff[l][j]+1);
			}

			for(int x = 0;x<3;x++)
			{
				vectorB[x] = 0;
				for(int y = 0;y<3;y++)
					vectorB[x] = vectorB[x] + vectorA[y]*ObjectInverseMatrix[x][y];
			}
			result = 0;
			for(int x = 0;x<3;x++)
				result = result + vectorB[x]*vectorA[x];

	        m_ObjectsMap[l][i] = exp(-0.5*result);
	    }
	}

	for(i = 0;i<3;i++)
		Histogram[i].resize(0);
	
}

template <class TInputImage, class TOutputImage>
void 
VectorFuzzyConnectednessImageFilter<TInputImage,TOutputImage>
::Compute_Scale()
{
	const int Tolerance = 13;
	int i,j,k,x,y,z,pslices,prow,pcol;
	int flag,edge_flag;
	int tti1;
	IndexType  index1,index2;
	IntVector  value1,value2, mean;
	double mask_f[3]; 
	double count_obj,count_nonobj;
	int slice,row,col;
	
	m_InputImage = this->GetInput();
	m_Size = m_InputImage->GetLargestPossibleRegion().GetSize();

	pslices = m_Size[2];
	prow = m_Size[1];
	pcol = m_Size[0];
	
	m_ScaleArray.resize(pslices*prow*pcol);

	for( slice = 0;slice<pslices;slice++)
		for( row = 0;row<prow;row++)
			for( col = 0;col<pcol;col++)
			{

			 index1[2] = slice;
			 index1[1] = row;
			 index1[0] = col;
			 value1 = m_InputImage->GetPixel(index1);

			  flag = 0;
			  edge_flag = 0;
			  
			  for(i = 0;i<3;i++)
				mask_f[i] = 0.0;

			  for (int zz = -1; zz <= 1; zz++)
				for (int yy = -1; yy <= 1; yy++)
				  for (int xx = -1; xx <= 1; xx++)
				  {
					 int  x = xx + col;
					 int  y = yy + row;
					 int  z = zz + slice;

					 index2[2] = z;
					 index2[1] = y;
					 index2[0] = x;	
					 
					 for(i = 0;i<3;i++)
					 {
						if (x >= 0 && y >= 0 && z >= 0 
							  && x < pcol && y < prow && z < pslices)
							  {
								value2 = m_InputImage->GetPixel(index2);
								mask_f[i]  = mask_f[i] + m_Mask[zz + 1][yy + 1][xx + 1]*(double) value2[i];
							  }
						else
							mask_f[i]  = mask_f[i] + m_Mask[zz + 1][yy + 1][xx + 1]*(double) value1[i];
					 }
				  }

			  for(i = 0;i<3;i++)
				mean[i] = (int) (mask_f[i]/m_MaskTotal + 0.5);

			  mean = value1;		
			  for (k = 1; k < 8 && !flag; k++)
			  {
			     count_obj = 0.0;
				 count_nonobj = 0.0;
				 for (i = 0; i < m_SpherePointsNum[k]; i++)
				 {
					x = col + m_SpherePointsLoc[k][i][0];
					y = row + m_SpherePointsLoc[k][i][1];
					z = slice + m_SpherePointsLoc[k][i][2];
					if (x < 0 || x >= pcol)
						x = col;
					if (y < 0 || y >= prow)
						y = row;
				    if (z < 0 || z >= pslices)
						z = slice;
					index1[2] = z;
					index1[1] = y;
					index1[0] = x;
					value1 = m_InputImage->GetPixel(index1);

					tti1 = 0;
					int temp[3];

					for(j=0;j<3;j++)
					{
						temp[j] = abs( value1[j] - mean[j]);
						tti1 = tti1 + temp[j]*m_PowerValue[j];

						if(temp[j]>=m_FeaturesThreshold[j])
							edge_flag = 1;
					}

					if(!edge_flag)
					{
						count_obj = count_obj + m_ScaleMap[tti1];
						count_nonobj = count_nonobj + 1.0 - m_ScaleMap[tti1];
					}
					else
					{
						edge_flag = 0;
						count_obj = count_obj;
						count_nonobj = count_nonobj + 1.0;
					}
				}
				if (100.0 * count_nonobj >= Tolerance * (count_nonobj + count_obj)) 
				{
				  m_ScaleArray[slice*prow*pcol + row * pcol + col] = k;
				  flag = 1;
				}
			  }

			  if (!flag)
				m_ScaleArray[slice*prow*pcol + row * pcol + col] = k;
			}

	m_ScaleMap.resize(0);

}


template <class TInputImage, class TOutputImage>
void 
VectorFuzzyConnectednessImageFilter<TInputImage,TOutputImage>
::Compute_Affinity()
{

	int i,j,k,pslices,prow,pcol;
	IndexType  index1,index2;
	IntVector  value1,value2, mean; 
	int x,y,z,tti1;
	double temp[3],sum_temp[3],sum_total[3];
	double tt1,tt2,inv_k,count;
	double sum_pos,count_pos;
	int col,row,slice,col1,row1,slice1;
	
	m_InputImage = this->GetInput();
	m_Size = m_InputImage->GetLargestPossibleRegion().GetSize();

	pslices = m_Size[2];
	prow = m_Size[1];
	pcol = m_Size[0];
	
	m_Xaffinity.resize(pslices*prow*pcol);
	m_Yaffinity.resize(pslices*prow*pcol);
	m_Zaffinity.resize(pslices*prow*pcol);

	unsigned short *temp2 = new unsigned short[pslices * prow * pcol];
	for( i = 0;i<pslices * prow * pcol;i++)
		temp2[i] = 0;

	/*  compute homogeneity-based affinity  */

     for ( slice = 0; slice < pslices; slice++)
		for ( row = 0; row < prow; row++)
		  for ( col = 0; col < pcol-1; col++)
		    {
			   slice1 = slice;
			   row1 = row;
			   col1 = col + 1;
			  	
			   int scale1 = (int) m_ScaleArray[slice * pcol * prow + row * pcol + col];
			   int scale2 = (int) m_ScaleArray[slice1 * pcol * prow + row1 * pcol + col1];
			   int iscale = (scale1<scale2 ? scale1:scale2);

			  for(i=0;i<3;i++)
			    {
			      sum_temp[i] = 0;
			      sum_total[i] = 0;
			      temp[i] = 0;
			    }
			   sum_pos = 0.0;
			   count_pos = 0.00001;
			 			  			  			  
			  for (k = 0; k < iscale; k++)  {
				  float temp_sum_pos = 0.0;

				  /*   EXPERIMENTS  */
				  tt1 = (double) iscale;
				  tt1 = 0.5 * tt1;
				  double inv_half_scale = -0.5 / pow(tt1, 2.0);
				  inv_k = exp(inv_half_scale * pow((double) k, 2.0));
				 	
				  
				  for(i = 0;i<3;i++)
				    sum_temp[i] = 0;
				  for (i = 0; i < m_SpherePointsNum[k]; i++) {
					  int xx = m_SpherePointsLoc[k][i][0];
					  int yy = m_SpherePointsLoc[k][i][1];
					  int zz = m_SpherePointsLoc[k][i][2];

					  int x1 = col + xx;
					  int x2 = col1 + xx;
					  int y1 = row + yy;
					  int y2 = row1 + yy;
					  int z1 = slice + zz;
					  int z2 = slice1 + zz;

					  if (x1 >= 0 && x2 >= 0 && y1 >= 0 && y2 >= 0
						  && z1 >= 0 && z2 >= 0 && x1 < pcol && x2 < pcol
						  && y1 < prow && y2 < prow && z1 < pslices && z2 < pslices)
					    {
					      tti1 = 0;
					      temp_sum_pos = 0;

	  					  index1[2] = z1;
						  index1[1] = y1;
						  index1[0] = x1;

						  index2[2] = z2;
						  index2[1] = y2;
						  index2[0] = x2;
	
						  value1 = m_InputImage->GetPixel(index1);
						  value2 = m_InputImage->GetPixel(index2);

						  int edge_flag = 0;

					      for(j=0;j<3;j++)
						  {
							temp[j] = (value1[j]-value2[j]);
							if(abs(temp[j]) >= m_FeaturesThreshold[j])
								edge_flag = 1;
							tti1 = tti1 + abs(temp[j])*m_PowerValue[j];
							temp_sum_pos = temp_sum_pos + pow(temp[j],2);
						  }
					      temp_sum_pos = sqrt(temp_sum_pos);

						  if(edge_flag)
							sum_pos = 1/temp_sum_pos;
						  else
						  {
							if(temp_sum_pos != 0)
								sum_pos = (1 - m_HomogeneityMap[tti1])/temp_sum_pos;
							else
								sum_pos = 0;
						  }
					      count_pos = count_pos + inv_k;
					     
					      for(j=0;j<3;j++)
							sum_temp[j] = sum_temp[j] + sum_pos*temp[j];
					    }
				  } 

				  for(j=0;j<3;j++)
				    sum_total[j] = sum_total[j] + inv_k*sum_temp[j];

			  }

			  tt1 = 0.0;
			  for(i = 0; i<3; i++)
			    tt1 = tt1 + pow(sum_total[i],2);
			  tt1  = sqrt(tt1);
			  tt2 =  1 - tt1 / count_pos;
			  
			  temp2[slice * prow * pcol + row * pcol + col] = (unsigned short) (4096 * tt2);
			  	
			  m_Xaffinity[slice * prow * pcol + row * pcol + col]= (unsigned short) (4096 * tt2);
		    }

     for (slice = 0; slice < pslices; slice++)
		for (row = 0; row < prow-1; row++)
		  for (col = 0; col < pcol; col++)
		    {
			  slice1 = slice;
			  row1 = row+1;
			  col1 = col;

			  int scale1 = (int) m_ScaleArray[slice * pcol * prow + row * pcol + col];
			  int scale2 = (int) m_ScaleArray[slice1 * pcol * prow + row1 * pcol + col1];
			  int iscale = (scale1<scale2 ? scale1:scale2);

			  for(i=0;i<3;i++)
			    {
			      sum_temp[i] = 0;
			      sum_total[i] = 0;
			      temp[i] = 0;
			    }
			  sum_pos = 0.0;
			  count_pos = 0.00001;
			  	
			  for (k = 0; k < iscale; k++)  {
				  float temp_sum_pos = 0.0;

				  /*  EXPERIMENTS  */
				  tt1 = (double) iscale;
				  tt1 = 0.5 * tt1;
				  double inv_half_scale = -0.5 / pow(tt1, 2.0);
				  inv_k = exp(inv_half_scale * pow((double) k, 2.0));
				 	
				  for(i = 0;i<3;i++)
				    sum_temp[i] = 0;
				  for (i = 0; i < m_SpherePointsNum[k]; i++) {
					  int xx = m_SpherePointsLoc[k][i][0];
					  int yy = m_SpherePointsLoc[k][i][1];
					  int zz = m_SpherePointsLoc[k][i][2];

					  int x1 = col + xx;
					  int x2 = col1 + xx;
					  int y1 = row + yy;
					  int y2 = row1 + yy;
					  int z1 = slice + zz;
					  int z2 = slice1 + zz;

					  if (x1 >= 0 && x2 >= 0 && y1 >= 0 && y2 >= 0
						  && z1 >= 0 && z2 >= 0 && x1 < pcol && x2 < pcol
						  && y1 < prow && y2 < prow && z1 < pslices && z2 < pslices)
					    {
					      tti1 = 0;
					      temp_sum_pos = 0;

	  					  index1[2] = z1;
						  index1[1] = y1;
						  index1[0] = x1;

						  index2[2] = z2;
						  index2[1] = y2;
						  index2[0] = x2;
	
						  value1 = m_InputImage->GetPixel(index1);
						  value2 = m_InputImage->GetPixel(index2);

						  int edge_flag = 0;

					      for(j=0;j<3;j++)
						  {
							temp[j] = (value1[j]-value2[j]);
							if(abs(temp[j]) >= m_FeaturesThreshold[j])
								edge_flag = 1;
							tti1 = tti1 + abs(temp[j])*m_PowerValue[j];
							temp_sum_pos = temp_sum_pos + pow(temp[j],2);
						  }
					      temp_sum_pos = sqrt(temp_sum_pos);

						  if(edge_flag)
							sum_pos = 1/temp_sum_pos;
						  else
						  {
							if(temp_sum_pos != 0)
								sum_pos = (1 - m_HomogeneityMap[tti1])/temp_sum_pos;
							else
								sum_pos = 0;
						  }
					      count_pos = count_pos + inv_k;
					     
					      for(j=0;j<3;j++)
							sum_temp[j] = sum_temp[j] + sum_pos*temp[j];
					    }
				  }

				  for(j=0;j<3;j++)
				    sum_total[j] = sum_total[j] + inv_k*sum_temp[j];

			  }

			  tt1 = 0.0;
			  for(i = 0; i<3; i++)
			    tt1 = tt1 + pow(sum_total[i],2);
			  tt1  = sqrt(tt1);
			  tt2 =  1 - tt1 / count_pos;
			    
			  m_Yaffinity[slice * prow * pcol + row * pcol + col] = (unsigned short) (4096 * tt2);
		    }

     for (slice = 0; slice < pslices -1; slice++)
		for (row = 0; row < prow; row++)
		  for (col = 0; col < pcol; col++)
		    {
			  slice1 = slice + 1;
			  row1 = row;
			  col1 = col;

			  int scale1 = (int) m_ScaleArray[slice * pcol * prow + row * pcol + col];
			  int scale2 = (int) m_ScaleArray[slice1 * pcol * prow + row1 * pcol + col1];
			  int iscale = (scale1<scale2 ? scale1:scale2);

			  for(i=0;i<3;i++)
			    {
			      sum_temp[i] = 0;
			      sum_total[i] = 0;
			      temp[i] = 0;
			    }
     		  sum_pos = 0.0;
			  count_pos = 0.00001;
			 			  			  			  
			  for (k = 0; k < iscale; k++)  {
				  float temp_sum_pos = 0.0;

				  /* EXPERIMENTS */
				  tt1 = (double) iscale;
				  tt1 = 0.5 * tt1;
				  double inv_half_scale = -0.5 / pow(tt1, 2.0);
				  inv_k = exp(inv_half_scale * pow((double) k, 2.0));
				 	
				  for(i = 0;i<3;i++)
				    sum_temp[i] = 0;
				  for (i = 0; i < m_SpherePointsNum[k]; i++) {
					  int xx = m_SpherePointsLoc[k][i][0];
					  int yy = m_SpherePointsLoc[k][i][1];
					  int zz = m_SpherePointsLoc[k][i][2];

					  int x1 = col + xx;
					  int x2 = col1 + xx;
					  int y1 = row + yy;
					  int y2 = row1 + yy;
					  int z1 = slice + zz;
					  int z2 = slice1 + zz;

					  if (x1 >= 0 && x2 >= 0 && y1 >= 0 && y2 >= 0
						  && z1 >= 0 && z2 >= 0 && x1 < pcol && x2 < pcol
						  && y1 < prow && y2 < prow && z1 < pslices && z2 < pslices)
					    {
					      tti1 = 0;
					      temp_sum_pos = 0;

	  					  index1[2] = z1;
						  index1[1] = y1;
						  index1[0] = x1;

						  index2[2] = z2;
						  index2[1] = y2;
						  index2[0] = x2;
	
						  value1 = m_InputImage->GetPixel(index1);
						  value2 = m_InputImage->GetPixel(index2);

						  int edge_flag = 0;

					      for(j=0;j<3;j++)
						  {
							temp[j] = (value1[j]-value2[j]);
							if(abs(temp[j]) >= m_FeaturesThreshold[j])
								edge_flag = 1;
							tti1 = tti1 + abs(temp[j])*m_PowerValue[j];
							temp_sum_pos = temp_sum_pos + pow(temp[j],2);
						  }
					      temp_sum_pos = sqrt(temp_sum_pos);

						  if(edge_flag)
							sum_pos = 1/temp_sum_pos;
						  else
						  {
							if(temp_sum_pos != 0)
								sum_pos = (1 - m_HomogeneityMap[tti1])/temp_sum_pos;
							else
								sum_pos = 0;
						  }
					      count_pos = count_pos + inv_k;
					     
					      for(j=0;j<3;j++)
							sum_temp[j] = sum_temp[j] + sum_pos*temp[j];
					    }
				  }

				  for(j=0;j<3;j++)
				    sum_total[j] = sum_total[j] + inv_k*sum_temp[j];

			  }

			  tt1 = 0.0;
			  for(i = 0; i<3; i++)
			    tt1 = tt1 + pow(sum_total[i],2);
			  tt1  = sqrt(tt1);
			  tt2 =  1 - tt1 / count_pos;
			    
			  m_Zaffinity[slice * prow * pcol + row * pcol + col]= (unsigned short) (4096 * tt2);
		    }

	  m_HomogeneityMap.resize(0);

	/*   compute object-based affinity */

	  IntVector* ObjectsOffset = new IntVector[m_Objects];

	  for(i = 0;i<m_Objects; i++)
		for(j = 0;j<3; j++)
		{
			ObjectsOffset[i][j] = 1;
			for(k = 0;k<j;k++)
				ObjectsOffset[i][j] = ObjectsOffset[i][j]*(m_ObjectsMaxDiff[i][k]+1);
		}

	  m_Material.resize(pslices*prow*pcol);

	  double *material,largest_material;
	  material = new double[m_Objects];
	    	
	  for (slice = 0; slice < pslices; slice++)
		for (row = 0; row < prow; row++)
		  for (col = 0; col < pcol; col++)
		    {
			  for (int object = 0; object < m_Objects; object++)
				{ 
				  int iscale = m_ScaleArray[slice * prow * pcol + row * pcol + col];

				  for(j=0;j<3;j++)
						temp[j] = 0;
				  count = 0.0;

				  for (k = 0; k < iscale; k++)
					{
					  tt1 = (double) iscale;
					  tt1 = 0.5 * tt1;
					  double inv_half_scale = -0.5 / pow(tt1, 2.0);
					  inv_k = exp(inv_half_scale * pow((double) k, 2.0));
					
					  for (i = 0; i < m_SpherePointsNum[k]; i++)
					    {
						  int xx = m_SpherePointsLoc[k][i][0];
						  int yy = m_SpherePointsLoc[k][i][1];
						  int zz = m_SpherePointsLoc[k][i][2];

					      x = col + xx;
						  y = row + yy;
					      z = slice + zz;

						  index1[2] = z;
						  index1[1] = y;
						  index1[0] = x;



					      if (x >= 0 && y >= 0 && z >= 0 
						  && x < pcol && y < prow && z < pslices)
							{
  							    value1 = m_InputImage->GetPixel(index1);
								for(j=0;j<3;j++)
									temp[j] = temp[j]+inv_k*value1[j];
						  		count = count + inv_k;
							}
					    }

					}

				  int edge_flag = 0;
				  tti1 = 0;
				  for(j=0;j<3;j++)
				    {
				      temp[j] = (int)abs(temp[j]/count - m_ObjectsMean[object][j]);
				      if (temp[j] > m_ObjectsMaxDiff[object][j])
							edge_flag = 1;
				      tti1 = tti1 + temp[j]*ObjectsOffset[object][j];
				    }
				  if(edge_flag ==1)
					material[object] = 0;
				  else
					material[object] = m_ObjectsMap[object][tti1] ;
				  
				}

			  double max_material = material[0];
			  for (int object = 1; object < m_Objects; object++)
				if (max_material < material[object])
				  max_material = material[object];
			  m_Material[slice * prow * pcol + row * pcol + col] = max_material;
			  if(largest_material < max_material)
			    largest_material = max_material;
		    } 

	  delete[] material;

	  for (slice = 0; slice < pslices; slice++)
		for (row = 0; row < prow; row++)
		  for (col = 0; col < pcol; col++)
		    m_Material[slice * prow * pcol + row * pcol + col] = 
		          (double) 4096 * m_Material[slice * prow * pcol + row * pcol + col] / largest_material;

	for(int object = 0;object<m_Objects;object++)
		m_ObjectsMap[object].resize(0);

	m_ScaleArray.resize(0);

	delete[] ObjectsOffset;
	
	/* combine homogeneity_based affinity and object_based affinity */

	for (slice = 0; slice < pslices; slice++)
		for (row = 0; row < prow; row++)
		  for (col = 0; col < pcol - 1; col++)  {
			  slice1 = slice;
			  row1 = row;
			  col1 = col + 1;
			  if(m_Material[slice * prow * pcol + row * pcol + col] < 
			     m_Material[slice1 * prow * pcol + row1 * pcol + col1])
					m_Xaffinity[slice * prow * pcol + row * pcol + col] = 
					 (unsigned short)sqrt(m_Material[slice * prow * pcol + row * pcol + col] 
							* (double)m_Xaffinity[slice*prow * pcol+row*pcol+col]+0.5);
			  else										
					m_Xaffinity[slice * prow * pcol + row * pcol + col] = 
					(unsigned short) sqrt(m_Material[slice1 * prow * pcol + row1 * pcol + col1] 
						   * (double)m_Xaffinity[slice*prow * pcol+row*pcol+col]+0.5);
			}

	  for (slice = 0; slice < pslices; slice++)
		for (row = 0; row < prow-1; row++)
		  for (col = 0; col < pcol; col++)	{
			  slice1 = slice;
			  row1 = row + 1;
			  col1 = col;
			  if(m_Material[slice * prow * pcol + row * pcol + col] < 
			     m_Material[slice1 * prow * pcol + row1 * pcol + col1])
					m_Yaffinity[slice * prow * pcol + row * pcol + col] = 
						(unsigned short)sqrt(m_Material[slice * prow * pcol + row * pcol + col] 
						 * (double)m_Yaffinity[slice*prow * pcol+row*pcol+col]+0.5);
			  else										
					m_Yaffinity[slice * prow * pcol + row * pcol + col] = 
					(unsigned short) sqrt(m_Material[slice1 * prow * pcol + row1 * pcol + col1] 
						 * (double)m_Yaffinity[slice*prow * pcol +row*pcol+col]+0.5);
			}

	  for (slice = 0; slice < pslices-1; slice++)
		for (row = 0; row < prow; row++)
		  for (col = 0; col < pcol; col++)  {
			  slice1 = slice+1;
			  row1 = row;
			  col1 = col;
			  if(m_Material[slice * prow * pcol + row * pcol + col] < 
			     m_Material[slice1 * prow * pcol + row1 * pcol + col1])
				m_Zaffinity[slice * prow * pcol + row * pcol + col] = 
					(unsigned short)sqrt(m_Material[slice * prow * pcol + row * pcol + col] 
								   * (double)m_Zaffinity[slice*prow*pcol + row*pcol+col]+0.5);
			  else										
					m_Zaffinity[slice * prow * pcol + row * pcol + col] = 
						(unsigned short)sqrt(m_Material[slice1 * prow * pcol + row1 * pcol + col1] 
					   * (double)m_Zaffinity[slice*prow * pcol +row*pcol+col]+0.5);
			}

	for (slice = 0; slice < pslices; slice++)
		for (row = 0; row < prow; row++)
			for (col = pcol - 1; col < pcol; col++)
				m_Xaffinity[slice * prow * pcol + row * pcol + col] = 0;

	for (slice = 0; slice < pslices; slice++)
		for (row = prow - 1; row < prow; row++)
			for (col = 0; col < pcol; col++)
				m_Yaffinity[slice * prow * pcol + row * pcol + col] = 0;

	for (slice = pslices - 1; slice < pslices; slice++)
		for (row = 0; row < prow; row++)
			for (col = 0; col < pcol; col++)
				m_Zaffinity[slice * prow * pcol + row * pcol + col] = 0;

	m_Material.resize(0);
}


template <class TInputImage, class TOutputImage>
void VectorFuzzyConnectednessImageFilter<TInputImage,TOutputImage>
::FastTracking(int object_flag)
{

	
	int nbor[6][3] = {
		{ 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 1 },
		{ -1, 0, 0},
		{ 0, -1, 0},
		{ 0, 0, -1}
	};

	ListType           chash[4096+1];
	ListType::iterator  iter;
	int                topIndex;
	int                pslices,prow,pcol;
	IndexType          current,index1,index2;
	unsigned short     affp[6];


	m_InputImage = this->GetInput();
	m_Size = m_InputImage->GetLargestPossibleRegion().GetSize();

	pslices = m_Size[2];
	prow = m_Size[1];
	pcol = m_Size[0];

	for(int i=0 ;i<=4096;i++)
		chash[i].clear();
			
	if(object_flag == 1)
	{
		/* object tracking...*/			
		for(iter = m_ObjectsSeed[m_SelectedObject-1].begin();iter!=m_ObjectsSeed[m_SelectedObject-1].end();iter++)
		{
			current = *iter;
			int x = current[0];
			int y = current[1];
			int z = current[2];
			m_ObjectFuzzyScene->SetPixel(current,4096);
			chash[4096].push_front(current);
		}


		topIndex = 4096;
		while((topIndex>0) && (!chash[topIndex].empty()))
		{
			current = chash[topIndex].back();
			chash[topIndex].pop_back();

			int x = current[0];
			int y = current[1];
			int z = current[2];

			int k = z*prow*pcol + y*pcol + x;

			affp[0] = m_Xaffinity[k];
			affp[1] = m_Yaffinity[k];
			affp[2] = m_Zaffinity[k];
			if(x>0)
				affp[3] = m_Xaffinity[k-1];
			if(y>0)
				affp[4] = m_Yaffinity[k-pcol];
			if(z>0)
				affp[5] = m_Zaffinity[k-prow*pcol];

			while((topIndex > 0) && (chash[topIndex].empty()))
				topIndex--;

			unsigned short pmax,pmin,affn;

			pmax = m_ObjectFuzzyScene->GetPixel(current);

			for (int ei = 0; ei < 6; ei++)
			{
				 int xx = current[0] + nbor[ei][0];
				 int yy = current[1] + nbor[ei][1];
				 int zz = current[2] + nbor[ei][2];

				  if (xx >= 0 && xx < pcol &&
					  yy >= 0 && yy < prow &&
					  zz >= 0 && zz < pslices)
				  {
					affn = affp[ei];
					pmin = (pmax < affn ? pmax: affn);

					index1[0] = xx;
					index1[1] = yy;
					index1[2] = zz;

					unsigned short value = m_ObjectFuzzyScene->GetPixel(index1);

					if (pmin > value)
					{
						if (value == 0)
						{
							chash[pmin].push_front(index1);
							if (pmin>topIndex) topIndex = pmin;
						}
						else
						{

							ListType::iterator iter;

							for( iter = chash[value].begin();iter !=chash[value].end();iter++)
							{
								index2 = *iter;
								if(index2 == index1)
								{
									chash[value].erase(iter);
									break;
								}
							}
							chash[pmin].push_front(index1);
							if (pmin>topIndex) topIndex = pmin;

						}
						m_ObjectFuzzyScene->SetPixel(index1,pmin);
					}
				  }
			}
		}

	}
	else if(object_flag == 0)
	{
		/* background tracking ... */
		for(int i = 0;i<m_Objects;i++)
		{
			if(i != (m_SelectedObject-1))
			{
				for(iter = m_ObjectsSeed[i].begin();iter!=m_ObjectsSeed[i].end();iter++)
				{
					current = *iter;
					int x = current[0];
					int y = current[1];
					int z = current[2];

					m_BackgroundFuzzyScene->SetPixel(current,4096);
					chash[4096].push_front(current);
				}
			}
		}


		topIndex = 4096;
		while((topIndex>0) && (!chash[topIndex].empty()))
		{
			current = chash[topIndex].front();
			chash[topIndex].pop_front();

			while((topIndex > 0) && (chash[topIndex].empty()))		
				topIndex--;

			int x = current[0];
			int y = current[1];
			int z = current[2];

			int k = z*prow*pcol + y*pcol + x;

			affp[0] = m_Xaffinity[k];
			affp[1] = m_Yaffinity[k];
			affp[2] = m_Zaffinity[k];
			if(x>0)
				affp[3] = m_Xaffinity[k-1];
			if(y>0)
				affp[4] = m_Yaffinity[k-pcol];
			if(z>0)
				affp[5] = m_Zaffinity[k-prow*pcol];

			unsigned short pmax,pmin,affn;
			pmax = m_BackgroundFuzzyScene->GetPixel(current);

			for (int ei = 0; ei < 6; ei++)
			{
				 int xx = current[0] + nbor[ei][0];
				 int yy = current[1] + nbor[ei][1];
				 int zz = current[2] + nbor[ei][2];

				  if (xx >= 0 && xx < pcol &&
					  yy >= 0 && yy < prow &&
					  zz >= 0 && zz < pslices)
				  {
					affn = affp[ei];
					pmin = (pmax < affn ? pmax: affn);

					index1[0] = xx;
					index1[1] = yy;
					index1[2] = zz;

					unsigned short value = m_BackgroundFuzzyScene->GetPixel(index1);

					if (pmin > value)
					{
						if (value == 0)
						{
							chash[pmin].push_front(index1);
							if (pmin>topIndex) topIndex = pmin;

						}
						else
						{

							ListType::iterator iter;

							for( iter = chash[value].begin();iter !=chash[value].end();iter++)
							{
								index2 = *iter;
								if(index2 == index1)
								{
									chash[value].erase(iter);
									break;
								}
							}
							chash[pmin].push_front(index1);
							if (pmin>topIndex) topIndex = pmin;

						}
						m_BackgroundFuzzyScene->SetPixel(index1,pmin);

					}
				  }
			}
		}
	}
	
}



template <class TInputImage, class TOutputImage>
void VectorFuzzyConnectednessImageFilter<TInputImage,TOutputImage>
::DoFuzzySegmentation()
{

  int pslices,prow,pcol;

  ScalePrepare();

  Compute_LookupTable();
  Compute_Scale();
  Compute_Affinity();

  m_InputImage = this->GetInput();
  m_SegmentObject = this->GetOutput(); 

  m_Size = m_InputImage->GetLargestPossibleRegion().GetSize();

  pslices = m_Size[2];
  prow = m_Size[1];
  pcol = m_Size[0];

  IndexType index = IndexType::ZeroIndex;

  UShortImage::RegionType region;
  region.SetSize(m_Size);
  region.SetIndex(index);
  m_ObjectFuzzyScene = UShortImage::New();  
  m_ObjectFuzzyScene->SetLargestPossibleRegion( region );
  m_ObjectFuzzyScene->SetBufferedRegion( region );
  m_ObjectFuzzyScene->SetRequestedRegion( region );
  m_ObjectFuzzyScene->Allocate();  

  m_BackgroundFuzzyScene = UShortImage::New();  
  m_BackgroundFuzzyScene->SetLargestPossibleRegion( region );
  m_BackgroundFuzzyScene->SetBufferedRegion( region );
  m_BackgroundFuzzyScene->SetRequestedRegion( region );
  m_BackgroundFuzzyScene->Allocate();

  RegionType region1;
  region1.SetSize(m_Size);
  region1.SetIndex(index);
  m_SegmentObject->SetLargestPossibleRegion( region1 );
  m_SegmentObject->SetBufferedRegion( region1 );
  m_SegmentObject->SetRequestedRegion( region1 );
  m_SegmentObject->Allocate();  

  SimpleImageRegionIterator <UShortImage> it1(this->m_ObjectFuzzyScene,region);
  SimpleImageRegionIterator <UShortImage> it2(this->m_BackgroundFuzzyScene,region);

  SimpleImageRegionIterator <OutputImageType> it3(this->m_SegmentObject,region1);

  it1.Begin();
  it2.Begin();

  while(!it1.IsAtEnd())
  {
	it1.Set(0);
	it2.Set(0);
	++it1;
	++it2;
  }

  int object_flag = 1;
  FastTracking(object_flag);

  int count, old_count;
  	
  if(m_Objects > 1)
  {

	object_flag = 0;
    FastTracking(object_flag);

	count = 0;
	old_count = 0;
	bool flag = 1;
	int iteration = 0; 

	
	while(flag)
	{
		old_count = count;
		count  = 0;
		flag = 0;
		iteration = iteration + 1;

		it1.Begin();
		it2.Begin();
		it3.Begin();

		while(!it1.IsAtEnd())
		{
			if(it1.Get() > it2.Get())
			{
				it3.Set(1);

				count++;

				IndexType current = it1.GetIndex();
				if(current[0]>0)
					m_Xaffinity[current[2]*prow*pcol + current[1]*pcol + current[0]-1] = 0;
				if(current[1]>0)
					m_Yaffinity[current[2]*prow*pcol + (current[1]-1)*pcol + current[0]] = 0;
				if(current[2]>0)
					m_Zaffinity[(current[2]-1)*prow*pcol + current[1]*pcol + current[0]] = 0;
			}
			else
				it3.Set(0);
			++it1;
			++it2;
			++it3;
		}
		if(count>old_count)
		{
			flag = 1;

			it2.Begin();
			while(!it2.IsAtEnd())
			{
				it2.Set(0);
				++it2;
			}
			FastTracking(0);  /* tracking background again */
		}
	}
  }


}

} // end namespace itk

#endif
