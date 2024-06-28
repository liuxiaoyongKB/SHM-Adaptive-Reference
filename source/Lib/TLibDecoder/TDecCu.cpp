/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * Copyright (c) 2010-2016, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     TDecCu.cpp
    \brief    CU decoder class
*/

#include "TDecCu.h"
#include "TLibCommon/TComTU.h"
#include "TLibCommon/TComPrediction.h"
#if SVC_EXTENSION
#include "TDecTop.h"
#endif


//! \ingroup TLibDecoder
//! \{

// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

TDecCu::TDecCu()
{
  m_ppcYuvResi = NULL;
  m_ppcYuvReco = NULL;
  m_ppcCU      = NULL;
}

TDecCu::~TDecCu()
{
}

#if SVC_EXTENSION
Void TDecCu::init(TDecTop** ppcDecTop, TDecEntropy* pcEntropyDecoder, TComTrQuant* pcTrQuant, TComPrediction* pcPrediction, UInt layerId)
{
  m_pcEntropyDecoder  = pcEntropyDecoder;
  m_pcTrQuant         = pcTrQuant;
  m_pcPrediction      = pcPrediction; 
  m_ppcTDecTop = ppcDecTop;
  m_layerId = layerId; 

#if LAYER_CTB
  memcpy(g_auiLayerZscanToRaster[m_layerId], g_auiZscanToRaster, sizeof( g_auiZscanToRaster ) );
  memcpy(g_auiLayerRasterToZscan[m_layerId], g_auiRasterToZscan, sizeof( g_auiRasterToZscan ) );
  memcpy(g_auiLayerRasterToPelX[m_layerId],  g_auiRasterToPelX,  sizeof( g_auiRasterToPelX ) );
  memcpy(g_auiLayerRasterToPelY[m_layerId],  g_auiRasterToPelY,  sizeof( g_auiRasterToPelY ) );
#endif
}
#else
Void TDecCu::init( TDecEntropy* pcEntropyDecoder, TComTrQuant* pcTrQuant, TComPrediction* pcPrediction)
{
  m_pcEntropyDecoder  = pcEntropyDecoder;
  m_pcTrQuant         = pcTrQuant;
  m_pcPrediction      = pcPrediction;
}
#endif

/**
 \param    uiMaxDepth      total number of allowable depth
 \param    uiMaxWidth      largest CU width
 \param    uiMaxHeight     largest CU height
 \param    chromaFormatIDC chroma format
 */
Void TDecCu::create( UInt uiMaxDepth, UInt uiMaxWidth, UInt uiMaxHeight, ChromaFormat chromaFormatIDC )
{
  m_uiMaxDepth = uiMaxDepth+1;  // uiMaxDepth == 4

  m_ppcYuvResi = new TComYuv*[m_uiMaxDepth-1];
  m_ppcYuvReco = new TComYuv*[m_uiMaxDepth-1];
  m_ppcCU      = new TComDataCU*[m_uiMaxDepth-1];

  for ( UInt ui = 0; ui < m_uiMaxDepth-1; ui++ )
  {
    UInt uiNumPartitions = 1<<( ( m_uiMaxDepth - ui - 1 )<<1 );
    UInt uiWidth  = uiMaxWidth  >> ui;
    UInt uiHeight = uiMaxHeight >> ui;

    // The following arrays (m_ppcYuvResi, m_ppcYuvReco and m_ppcCU) are only required for CU depths
    // although data is allocated for all possible depths of the CU/TU tree except the last.
    // Since the TU tree will always include at least one additional depth greater than the CU tree,
    // there will be enough entries for these arrays.
    // (Section 7.4.3.2: "The CVS shall not contain data that result in (Log2MinTrafoSize) MinTbLog2SizeY
    //                    greater than or equal to MinCbLog2SizeY")
    // TODO: tidy the array allocation given the above comment.

    m_ppcYuvResi[ui] = new TComYuv;    m_ppcYuvResi[ui]->create( uiWidth, uiHeight, chromaFormatIDC );
    m_ppcYuvReco[ui] = new TComYuv;    m_ppcYuvReco[ui]->create( uiWidth, uiHeight, chromaFormatIDC );
    m_ppcCU     [ui] = new TComDataCU; m_ppcCU     [ui]->create( chromaFormatIDC, uiNumPartitions, uiWidth, uiHeight, true, uiMaxWidth >> (m_uiMaxDepth - 1) );
  }

  m_bDecodeDQP = false;
  m_IsChromaQpAdjCoded = false;

  // initialize partition order.
  UInt* piTmp = &g_auiZscanToRaster[0];
  initZscanToRaster(m_uiMaxDepth, 1, 0, piTmp);
  initRasterToZscan( uiMaxWidth, uiMaxHeight, m_uiMaxDepth );

  // initialize conversion matrix from partition index to pel
  initRasterToPelXY( uiMaxWidth, uiMaxHeight, m_uiMaxDepth );
}

Void TDecCu::destroy()
{
  for ( UInt ui = 0; ui < m_uiMaxDepth-1; ui++ )
  {
    m_ppcYuvResi[ui]->destroy(); delete m_ppcYuvResi[ui]; m_ppcYuvResi[ui] = NULL;
    m_ppcYuvReco[ui]->destroy(); delete m_ppcYuvReco[ui]; m_ppcYuvReco[ui] = NULL;
    m_ppcCU     [ui]->destroy(); delete m_ppcCU     [ui]; m_ppcCU     [ui] = NULL;
  }

  delete [] m_ppcYuvResi; m_ppcYuvResi = NULL;
  delete [] m_ppcYuvReco; m_ppcYuvReco = NULL;
  delete [] m_ppcCU     ; m_ppcCU      = NULL;
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

/** 
 Parse a CTU.  解析CTU的信息
 \param    pCtu                      [in/out] pointer to CTU data structure
 \param    isLastCtuOfSliceSegment   [out]    true, if last CTU of the slice segment
 */
Void TDecCu::decodeCtu( TComDataCU* pCtu, Bool& isLastCtuOfSliceSegment )
{
  if ( pCtu->getSlice()->getPPS()->getUseDQP() )
  {
    setdQPFlag(true);
  }

  if ( pCtu->getSlice()->getUseChromaQpAdj() )
  {
    setIsChromaQpAdjCoded(true);
  }
  
  // start from the top level CU
  xDecodeCU( pCtu, 0, 0, isLastCtuOfSliceSegment);
}

/** 
 Decoding process for a CTU.
 \param    pCtu                      [in/out] pointer to CTU data structure
 */
Void TDecCu::decompressCtu( TComDataCU* pCtu )
{
  xDecompressCU( pCtu, 0,  0 );
}

// ====================================================================================================================
// Protected member functions
// ====================================================================================================================

//! decode end-of-slice flag
Bool TDecCu::xDecodeSliceEnd( TComDataCU* pcCU, UInt uiAbsPartIdx )
{
  UInt uiIsLastCtuOfSliceSegment;

  if (pcCU->isLastSubCUOfCtu(uiAbsPartIdx))
  {
    m_pcEntropyDecoder->decodeTerminatingBit( uiIsLastCtuOfSliceSegment );
  }
  else
  {
    uiIsLastCtuOfSliceSegment=0;
  }

  return uiIsLastCtuOfSliceSegment>0;
}

//! decode CU block recursively
Void TDecCu::xDecodeCU( TComDataCU*const pcCU, const UInt uiAbsPartIdx, const UInt uiDepth, Bool &isLastCtuOfSliceSegment)
{
  TComPic* pcPic        = pcCU->getPic();
  const TComSPS &sps    = pcPic->getPicSym()->getSPS();
  const TComPPS &pps    = pcPic->getPicSym()->getPPS();
  const UInt maxCuWidth = sps.getMaxCUWidth();
  const UInt maxCuHeight= sps.getMaxCUHeight();
  UInt uiCurNumParts    = pcPic->getNumPartitionsInCtu() >> (uiDepth<<1);  // 当前深度CU包含的4x4_CU（基本编码单元）的数量
  UInt uiQNumParts      = uiCurNumParts>>2;  // 4个子节点：下一级深度CU包含的4x4_CU（基本编码单元）的数量


  Bool bBoundary = false;
  UInt uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
  UInt uiRPelX   = uiLPelX + (maxCuWidth>>uiDepth)  - 1;
  UInt uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
  UInt uiBPelY   = uiTPelY + (maxCuHeight>>uiDepth) - 1;

  if( ( uiRPelX < sps.getPicWidthInLumaSamples() ) && ( uiBPelY < sps.getPicHeightInLumaSamples() ) )
  {
    m_pcEntropyDecoder->decodeSplitFlag( pcCU, uiAbsPartIdx, uiDepth );  // 解码划分标志（设置当前CU的深度信息uiDepth）
  }
  else
  {
    bBoundary = true;
  }

  /*
  // 调试
  fstream outDepth("CTU_Depth.txt", ios::out | ios::app);
  UChar RasterDepthCTU[256] = { 0 };
  for ( Int i = 0; i < pcPic->getNumPartitionsInCtu(); i++)
  {
	  RasterDepthCTU[i] = pcCU->getDepth(g_auiRasterToZscan[i]);
	  outDepth << (Int)(RasterDepthCTU[i]) << ' ';
	  if ((i + 1) % 16 == 0)
	  {
		  outDepth << '\n';
	  }
  }
  outDepth.close();
  */

  if( ( ( uiDepth < pcCU->getDepth( uiAbsPartIdx ) ) && ( uiDepth < sps.getLog2DiffMaxMinCodingBlockSize() ) ) || bBoundary )
  {
    UInt uiIdx = uiAbsPartIdx;
    if( uiDepth == pps.getMaxCuDQPDepth() && pps.getUseDQP())
    {
      setdQPFlag(true);
      pcCU->setQPSubParts( pcCU->getRefQP(uiAbsPartIdx), uiAbsPartIdx, uiDepth ); // set QP to default QP
    }

    if( uiDepth == pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth() && pcCU->getSlice()->getUseChromaQpAdj() )
    {
      setIsChromaQpAdjCoded(true);
    }

    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++ )
    {
      uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiIdx] ];
      uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiIdx] ];

      if ( !isLastCtuOfSliceSegment && ( uiLPelX < sps.getPicWidthInLumaSamples() ) && ( uiTPelY < sps.getPicHeightInLumaSamples() ) )
      {
        xDecodeCU( pcCU, uiIdx, uiDepth+1, isLastCtuOfSliceSegment );
      }
      else
      {
        pcCU->setOutsideCUPart( uiIdx, uiDepth+1 );
      }

      uiIdx += uiQNumParts;  // 同一划分深度中四叉树的下一个节点的Z扫描索引
    }
    if( uiDepth == pps.getMaxCuDQPDepth() && pps.getUseDQP())
    {
      if ( getdQPFlag() )
      {
        UInt uiQPSrcPartIdx = uiAbsPartIdx;
        pcCU->setQPSubParts( pcCU->getRefQP( uiQPSrcPartIdx ), uiAbsPartIdx, uiDepth ); // set QP to default QP
      }
    }
    return;
  }

  if( uiDepth <= pps.getMaxCuDQPDepth() && pps.getUseDQP())
  {
    setdQPFlag(true);
    pcCU->setQPSubParts( pcCU->getRefQP(uiAbsPartIdx), uiAbsPartIdx, uiDepth ); // set QP to default QP
  }

  if( uiDepth <= pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth() && pcCU->getSlice()->getUseChromaQpAdj() )
  {
    setIsChromaQpAdjCoded(true);
  }

  if (pps.getTransquantBypassEnabledFlag())
  {
    m_pcEntropyDecoder->decodeCUTransquantBypassFlag( pcCU, uiAbsPartIdx, uiDepth );
  }

  // decode CU mode and the partition size
  if( !pcCU->getSlice()->isIntra())
  {
    m_pcEntropyDecoder->decodeSkipFlag( pcCU, uiAbsPartIdx, uiDepth );  // 解码skip标志
  }
 
#if SVC_EXTENSION
  // Check CU skip for higher layer IRAP skip flag
  if( pcCU->getSlice()->getVPS()->getHigherLayerIrapSkipFlag() && pcCU->getSlice()->getVPS()->getSingleLayerForNonIrapFlag() && pcCU->getPic()->getLayerId() > 0 )
  {
    Bool lowerLayerExist = false;
    for( Int i = 0; i < pcCU->getPic()->getLayerId(); i++ )
    {
      if(pcCU->getSlice()->getBaseColPic(pcCU->getSlice()->getInterLayerPredLayerIdc(i)))
      {
        lowerLayerExist = true;
      }
    }

    if( lowerLayerExist && !pcCU->isSkipped(uiAbsPartIdx) )
    {
      printf( "Warning: CU is not skipped with enabled higher layer IRAP skip flag\n" );
    }
  }
#endif
 
  if( pcCU->isSkipped(uiAbsPartIdx) )  // Skip模式
  {
    m_ppcCU[uiDepth]->copyInterPredInfoFrom( pcCU, uiAbsPartIdx, REF_PIC_LIST_0 );
    m_ppcCU[uiDepth]->copyInterPredInfoFrom( pcCU, uiAbsPartIdx, REF_PIC_LIST_1 );
    TComMvField cMvFieldNeighbours[MRG_MAX_NUM_CANDS << 1]; // double length for mv of both lists
    UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
    Int numValidMergeCand = 0;
    for( UInt ui = 0; ui < m_ppcCU[uiDepth]->getSlice()->getMaxNumMergeCand(); ++ui )
    {
      uhInterDirNeighbours[ui] = 0;
    }
    m_pcEntropyDecoder->decodeMergeIndex( pcCU, 0, uiAbsPartIdx, uiDepth );  // 解码Merge索引
    UInt uiMergeIndex = pcCU->getMergeIndex(uiAbsPartIdx);
    m_ppcCU[uiDepth]->getInterMergeCandidates( 0, 0, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand, uiMergeIndex );
    pcCU->setInterDirSubParts( uhInterDirNeighbours[uiMergeIndex], uiAbsPartIdx, 0, uiDepth );

    TComMv cTmpMv( 0, 0 );
    for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
    {
      if ( pcCU->getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
      {
        pcCU->setMVPIdxSubParts( 0, RefPicList( uiRefListIdx ), uiAbsPartIdx, 0, uiDepth);
        pcCU->setMVPNumSubParts( 0, RefPicList( uiRefListIdx ), uiAbsPartIdx, 0, uiDepth);
        pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->setAllMvd( cTmpMv, SIZE_2Nx2N, uiAbsPartIdx, uiDepth );
        pcCU->getCUMvField( RefPicList( uiRefListIdx ) )->setAllMvField( cMvFieldNeighbours[ 2*uiMergeIndex + uiRefListIdx ], SIZE_2Nx2N, uiAbsPartIdx, uiDepth );
      }
    }
    xFinishDecodeCU( pcCU, uiAbsPartIdx, uiDepth, isLastCtuOfSliceSegment );
    return;
  }

  m_pcEntropyDecoder->decodePredMode( pcCU, uiAbsPartIdx, uiDepth );
  m_pcEntropyDecoder->decodePartSize( pcCU, uiAbsPartIdx, uiDepth );

  if (pcCU->isIntra( uiAbsPartIdx ) && pcCU->getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N )
  {
    m_pcEntropyDecoder->decodeIPCMInfo( pcCU, uiAbsPartIdx, uiDepth );

    if(pcCU->getIPCMFlag(uiAbsPartIdx))
    {
      xFinishDecodeCU( pcCU, uiAbsPartIdx, uiDepth, isLastCtuOfSliceSegment );
      return;
    }
  }

  // prediction mode ( Intra : direction mode, Inter : Mv, reference idx )  通过熵解码获取预测信息
  m_pcEntropyDecoder->decodePredInfo( pcCU, uiAbsPartIdx, uiDepth, m_ppcCU[uiDepth]);

  // Coefficient decoding  变换系数解码
  Bool bCodeDQP = getdQPFlag();
  Bool isChromaQpAdjCoded = getIsChromaQpAdjCoded();
  m_pcEntropyDecoder->decodeCoeff( pcCU, uiAbsPartIdx, uiDepth, bCodeDQP, isChromaQpAdjCoded );
  setIsChromaQpAdjCoded( isChromaQpAdjCoded );
  setdQPFlag( bCodeDQP );
  xFinishDecodeCU( pcCU, uiAbsPartIdx, uiDepth, isLastCtuOfSliceSegment );
}

Void TDecCu::xFinishDecodeCU( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth, Bool &isLastCtuOfSliceSegment)
{
  if(  pcCU->getSlice()->getPPS()->getUseDQP())
  {
    pcCU->setQPSubParts( getdQPFlag()?pcCU->getRefQP(uiAbsPartIdx):pcCU->getCodedQP(), uiAbsPartIdx, uiDepth ); // set QP
  }

  if (pcCU->getSlice()->getUseChromaQpAdj() && !getIsChromaQpAdjCoded())
  {
    pcCU->setChromaQpAdjSubParts( pcCU->getCodedChromaQpAdj(), uiAbsPartIdx, uiDepth ); // set QP
  }

  isLastCtuOfSliceSegment = xDecodeSliceEnd( pcCU, uiAbsPartIdx );
}

Void TDecCu::xDecompressCU( TComDataCU* pCtu, UInt uiAbsPartIdx,  UInt uiDepth )
{
  TComPic* pcPic = pCtu->getPic();
  TComSlice * pcSlice = pCtu->getSlice();
  const TComSPS &sps=*(pcSlice->getSPS());

  Bool bBoundary = false;
  UInt uiLPelX   = pCtu->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
  UInt uiRPelX   = uiLPelX + (sps.getMaxCUWidth()>>uiDepth)  - 1;
  UInt uiTPelY   = pCtu->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
  UInt uiBPelY   = uiTPelY + (sps.getMaxCUHeight()>>uiDepth) - 1;

  if( ( uiRPelX >= sps.getPicWidthInLumaSamples() ) || ( uiBPelY >= sps.getPicHeightInLumaSamples() ) )
  {
    bBoundary = true;
  }

  if( ( ( uiDepth < pCtu->getDepth( uiAbsPartIdx ) ) && ( uiDepth < sps.getLog2DiffMaxMinCodingBlockSize() ) ) || bBoundary )
  {
    UInt uiNextDepth = uiDepth + 1;  // 下一级深度
    UInt uiQNumParts = pCtu->getTotalNumPart() >> (uiNextDepth<<1);  // 下一级深度CU包含的4x4_CU（基本编码单元）的数量
    UInt uiIdx = uiAbsPartIdx;
    for ( UInt uiPartIdx = 0; uiPartIdx < 4; uiPartIdx++ )
    {
      uiLPelX = pCtu->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiIdx] ];
      uiTPelY = pCtu->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiIdx] ];

      if( ( uiLPelX < sps.getPicWidthInLumaSamples() ) && ( uiTPelY < sps.getPicHeightInLumaSamples() ) )
      {
        xDecompressCU(pCtu, uiIdx, uiNextDepth );
      }

      uiIdx += uiQNumParts;  // 同一划分深度中四叉树的下一个节点的Z扫描索引
    }
    return;
  }

  // Residual reconstruction
  m_ppcYuvResi[uiDepth]->clear();  // YUV的残差

  m_ppcCU[uiDepth]->copySubCU( pCtu, uiAbsPartIdx );  // 把pCtu的信息复制给四叉树划分得到的子CU，m_ppcCU[uiDepth]是当前深度的子CU

  switch( m_ppcCU[uiDepth]->getPredictionMode(0) )
  {
    case MODE_INTER:
      xReconInter( m_ppcCU[uiDepth], uiDepth );  // 获得使用帧间预测模式的CU的重建YUV
      break;
    case MODE_INTRA:
      xReconIntraQT( m_ppcCU[uiDepth], uiDepth );  // 获得使用帧内预测模式的CU的重建YUV
      break;
    default:
      assert(0);
      break;
  }

#if DEBUG_STRING
  const PredMode predMode=m_ppcCU[uiDepth]->getPredictionMode(0);
  if (DebugOptionList::DebugString_Structure.getInt()&DebugStringGetPredModeMask(predMode))
  {
    PartSize eSize=m_ppcCU[uiDepth]->getPartitionSize(0);
    std::ostream &ss(std::cout);

    ss <<"###: " << (predMode==MODE_INTRA?"Intra   ":"Inter   ") << partSizeToString[eSize] << " CU at " << m_ppcCU[uiDepth]->getCUPelX() << ", " << m_ppcCU[uiDepth]->getCUPelY() << " width=" << UInt(m_ppcCU[uiDepth]->getWidth(0)) << std::endl;
  }
#endif

  if ( m_ppcCU[uiDepth]->isLosslessCoded(0) && (m_ppcCU[uiDepth]->getIPCMFlag(0) == false))
  {
    xFillPCMBuffer(m_ppcCU[uiDepth], uiDepth);
  }

  xCopyToPic( m_ppcCU[uiDepth], pcPic, uiAbsPartIdx, uiDepth );  // 把CU的重建YUV复制到图像
}

Void TDecCu::xReconInter( TComDataCU* pcCU, UInt uiDepth )
{

  // 运动补偿获得YUV的预测值，m_ppcYuvReco[uiDepth]存储当前CU的YUV的预测值
  // inter prediction
  m_pcPrediction->motionCompensation( pcCU, m_ppcYuvReco[uiDepth] );

#if DEBUG_STRING
  const Int debugPredModeMask=DebugStringGetPredModeMask(MODE_INTER);
  if (DebugOptionList::DebugString_Pred.getInt()&debugPredModeMask)
  {
    printBlockToStream(std::cout, "###inter-pred: ", *(m_ppcYuvReco[uiDepth]));
  }
#endif

  // inter recon
  xDecodeInterTexture( pcCU, uiDepth );  // 解码YUV的残差

#if DEBUG_STRING
  if (DebugOptionList::DebugString_Resi.getInt()&debugPredModeMask)
  {
    printBlockToStream(std::cout, "###inter-resi: ", *(m_ppcYuvResi[uiDepth]));
  }
#endif

  // clip for only non-zero cbp case
  if  ( pcCU->getQtRootCbf( 0) )  // 获得非全零CU的重建值
  {
    m_ppcYuvReco[uiDepth]->addClip( m_ppcYuvReco[uiDepth], m_ppcYuvResi[uiDepth], 0, pcCU->getWidth( 0 ), pcCU->getSlice()->getSPS()->getBitDepths() );
  }
  else  // CBF == 0，当前CU的变换参数都为0，残差也为0，YUV预测值就是最终的YUV重建值
  {
    m_ppcYuvReco[uiDepth]->copyPartToPartYuv( m_ppcYuvReco[uiDepth],0, pcCU->getWidth( 0 ),pcCU->getHeight( 0 ));
  }
#if DEBUG_STRING
  if (DebugOptionList::DebugString_Reco.getInt()&debugPredModeMask)
  {
    printBlockToStream(std::cout, "###inter-reco: ", *(m_ppcYuvReco[uiDepth]));
  }
#endif

}

#if DecodedELImage_Restoration
#if ImageRestorationRefinement
// 重建丢失EL图像中的帧内PU的YUV
Void TDecCu::xReconIntraQTLostPic(TComDataCU* pcCU, UInt uiZorderIdx, UInt uiDepth)
{
	// Residual reconstruction
	m_ppcYuvResi[uiDepth]->clear();  // YUV的残差

	if (pcCU->getIPCMFlag(0))
	{
		xReconPCM(pcCU, uiDepth);
		return;
	}
	const UInt numChType = pcCU->getPic()->getChromaFormat() != CHROMA_400 ? 2 : 1;
	for (UInt chType = CHANNEL_TYPE_LUMA; chType < numChType; chType++)
	{
		const ChannelType chanType = ChannelType(chType);
		const Bool NxNPUHas4Parts = ::isChroma(chanType) ? enable4ChromaPUsInIntraNxNCU(pcCU->getPic()->getChromaFormat()) : true;
		const UInt uiInitTrDepth = (pcCU->getPartitionSize(0) != SIZE_2Nx2N && NxNPUHas4Parts ? 1 : 0);

		TComTURecurse tuRecurseCU(pcCU, 0);
		TComTURecurse tuRecurseWithPU(tuRecurseCU, false, (uiInitTrDepth == 0) ? TComTU::DONT_SPLIT : TComTU::QUAD_SPLIT);

		do
		{
			xIntraRecQT(m_ppcYuvReco[uiDepth], m_ppcYuvReco[uiDepth], m_ppcYuvResi[uiDepth], chanType, tuRecurseWithPU);
		} while (tuRecurseWithPU.nextSection(tuRecurseCU));
	}

	xCopyToPic(pcCU, pcCU->getPic(), uiZorderIdx, uiDepth);  // 把CU的重建YUV和残差YUV复制到丢失图像
}

// 使用BL图像和同位参考图像的重建YUV的差重建丢失EL图像中的PU的YUV
Void TDecCu::ReconLostPicFromBLPic(TComDataCU* pcCU, TComDataCU* pcBLCU, UInt uiZorderIdx, UInt uiDepth)
{
	// Residual reconstruction
	m_ppcYuvResi[uiDepth]->clear();

	// 运动补偿获得YUV的预测值，m_ppcYuvReco[uiDepth]存储当前CU的YUV的预测值
	// inter prediction
	//m_pcPrediction->motionCompensation(pcCU, m_ppcYuvReco[uiDepth]);

	// clip for only non-zero cbf case
	if (pcCU->getQtRootCbf(0))  // 获得非全零CU的重建值
	{
		Int         iWidth;
		Int         iHeight;
		UInt        uiPartAddr;
		// 给当前CU中的每个PU复制同位参考图像的重建YUV的差
		for (Int iPartIdx = 0; iPartIdx < pcCU->getNumPartitions(); iPartIdx++)
		{
			pcCU->getPartIndexAndSize(iPartIdx, uiPartAddr, iWidth, iHeight);
			Int iRefIdx = pcCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx(uiPartAddr);    assert(iRefIdx >= 0);

			TComMv     cMv = pcCU->getCUMvField(REF_PIC_LIST_0)->getMv(uiPartAddr);  // 使用BL图像中PU的MV
			pcCU->clipMv(cMv);

			TComPicYuv *refPicResi = pcCU->getSlice()->getRefPic(REF_PIC_LIST_0, iRefIdx)->getPicYuvResi();  // 同位参考图像的重建YUV的差
			TComPicYuv *PicRecoBL = pcBLCU->getPic()->getPicYuvRec();  // 同位BL图像的重建YUV

			for (UInt comp = COMPONENT_Y; comp < m_ppcYuvReco[uiDepth]->getNumberValidComponents(); comp++)  // 分别复制3个分量的残差
			{
				const ComponentID compID = ComponentID(comp);
				UInt cxWidth = iWidth >> refPicResi->getComponentScaleX(compID);
				UInt cxHeight = iHeight >> refPicResi->getComponentScaleY(compID);
				Int     refStride = refPicResi->getStride(compID);
				Int     dstStride = m_ppcYuvResi[uiDepth]->getStride(compID);
				Int shiftHor = (2 + refPicResi->getComponentScaleX(compID));
				Int shiftVer = (2 + refPicResi->getComponentScaleY(compID));

				Int     refOffset = (cMv.getHor() >> shiftHor) + (cMv.getVer() >> shiftVer) * refStride;

				Pel*    refYuvResi = refPicResi->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr) + refOffset;  // 同位参考图像中PU的重建YUV的差的原点像素
				Pel*    PicYuvRecoBL = PicRecoBL->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr);  // BL图像中同位PU的重建YUV的原点像素

				Pel*    dst = m_ppcYuvResi[uiDepth]->getAddr(compID, uiPartAddr);  // 当前PU残差值的原点像素
				Pel*    dstReco = m_ppcYuvReco[uiDepth]->getAddr(compID, uiPartAddr);  // 当前PU重建值的原点像素

				for (Int y = cxHeight; y != 0; y--)
				{
					::memcpy(dst, refYuvResi, sizeof(Pel)*cxWidth);  // 复制参考图像中的重建YUV的差
					//::memset(dst, 0, sizeof(Pel)*cxWidth);
					dst += dstStride;
					refYuvResi += refStride;

					::memcpy(dstReco, PicYuvRecoBL, sizeof(Pel)*cxWidth);  // 复制BL图像中的重建YUV
					dstReco += dstStride;
					PicYuvRecoBL += refStride;
				}

			}

			/*for (UInt comp = COMPONENT_Y; comp < m_ppcYuvReco[uiDepth]->getNumberValidComponents(); comp++)  // 分别复制3个分量的残差
			{
				const ComponentID compID = ComponentID(comp);
				UInt cxWidth = iWidth >> refPicResi->getComponentScaleX(compID);
				UInt cxHeight = iHeight >> refPicResi->getComponentScaleY(compID);
				Int     refStride = refPicResi->getStride(compID);
				Int     dstStride = m_ppcYuvResi[uiDepth]->getStride(compID);
				Int shiftHor = (2 + refPicResi->getComponentScaleX(compID));
				Int shiftVer = (2 + refPicResi->getComponentScaleY(compID));

				Int     refOffset = (cMv.getHor() >> shiftHor) + (cMv.getVer() >> shiftVer) * refStride;

				Pel*    refYuvResi = refPicResi->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr) + refOffset;  // 同位参考图像中PU的重建YUV的差的原点像素
				Pel*    PicYuvRecoBL = PicRecoBL->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr);  // BL图像中同位PU的重建YUV的原点像素

				Pel*    dst = m_ppcYuvResi[uiDepth]->getAddr(compID, uiPartAddr);  // 当前PU残差值的原点像素
				Pel*    dstReco = m_ppcYuvReco[uiDepth]->getAddr(compID, uiPartAddr);  // 当前PU重建值的原点像素

				// 预测值的Cb Cr首地址
				Pel*    refYuvPredCb = m_ppcYuvReco[uiDepth]->getAddr(COMPONENT_Cb, uiPartAddr);  // 同位参考图像中PU的重建YUV的差的原点像素
				Pel*    refYuvPredCr = m_ppcYuvReco[uiDepth]->getAddr(COMPONENT_Cr, uiPartAddr);  // BL图像中同位PU的重建YUV的原点像素
				Int     refStrideCbCr = m_ppcYuvReco[uiDepth]->getStride(COMPONENT_Cb);

				// BL图像的Cb Cr首地址
				Pel*    PicYuvRecoBLCb = PicRecoBL->getAddr(COMPONENT_Cb, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr);  // BL图像中同位PU的重建YUV的原点像素
				Pel*    PicYuvRecoBLCr = PicRecoBL->getAddr(COMPONENT_Cr, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr);  // BL图像中同位PU的重建YUV的原点像素
				Int     RecoBLStrideCbCr = PicRecoBL->getStride(COMPONENT_Cb);

				if (compID == COMPONENT_Y)
				{
					for (Int y = 0; y < cxHeight; y += 2)
					{
						for (Int x = 0; x < cxWidth; x += 2)
						{
							// 亮度分量的YUV残差
							dst[x] = ((PicYuvRecoBLCb[x >> 1] - refYuvPredCb[x >> 1]) + (PicYuvRecoBLCr[x >> 1] - refYuvPredCr[x >> 1])) >> 1;
							dst[x + 1] = ((PicYuvRecoBLCb[x >> 1] - refYuvPredCb[x >> 1]) + (PicYuvRecoBLCr[x >> 1] - refYuvPredCr[x >> 1])) >> 1;
							dst[x + dstStride] = ((PicYuvRecoBLCb[x >> 1] - refYuvPredCb[x >> 1]) + (PicYuvRecoBLCr[x >> 1] - refYuvPredCr[x >> 1])) >> 1;
							dst[x + dstStride + 1] = ((PicYuvRecoBLCb[x >> 1] - refYuvPredCb[x >> 1]) + (PicYuvRecoBLCr[x >> 1] - refYuvPredCr[x >> 1])) >> 1;
						}
						dst += dstStride;
						dst += dstStride;
						PicYuvRecoBLCb += RecoBLStrideCbCr;
						PicYuvRecoBLCr += RecoBLStrideCbCr;
						refYuvPredCb += refStrideCbCr;
						refYuvPredCr += refStrideCbCr;
					}
				}
				else
				{
					for (Int y = cxHeight; y != 0; y--)
					{
						// 色度分量YUV残差为0
						::memset(dst, 0, sizeof(Pel)*cxWidth);
						dst += dstStride;
						//refYuvResi += refStride;

						::memcpy(dstReco, PicYuvRecoBL, sizeof(Pel)*cxWidth);  // 复制BL图像中的重建YUV
						dstReco += dstStride;
						PicYuvRecoBL += refStride;

					}
				}
			}*/

		}

		m_ppcYuvReco[uiDepth]->addClip(m_ppcYuvReco[uiDepth], m_ppcYuvResi[uiDepth], 0, pcCU->getWidth(0), pcCU->getSlice()->getSPS()->getBitDepths());
	}
	else  // CBF == 0，当前CU的变换参数都为0，残差也为0，YUV预测值就是最终的YUV重建值（skip模式）
	{
		// 运动补偿获得YUV的预测值，m_ppcYuvReco[uiDepth]存储当前CU的YUV的预测值
        // inter prediction
		m_pcPrediction->motionCompensation(pcCU, m_ppcYuvReco[uiDepth]);
	}
	
	xCopyToPic(pcCU, pcCU->getPic(), uiZorderIdx, uiDepth);  // 把CU的重建YUV和残差YUV复制到丢失图像
}



// 使用前一帧EL图像重建丢失EL图像中的PU的YUV
Void TDecCu::ReconLostPicFromLastELPic(TComDataCU* pcCU, TComDataCU* pcLastELCU, UInt uiZorderIdx, UInt uiDepth)
{
	// Residual reconstruction
	m_ppcYuvResi[uiDepth]->clear();

	// 运动补偿获得YUV的预测值，m_ppcYuvReco[uiDepth]存储当前CU的YUV的预测值
	// inter prediction
	m_pcPrediction->motionCompensation(pcCU, m_ppcYuvReco[uiDepth]);

	// inter recon
	//xDecodeInterTexture(pcCU, uiDepth);  // 解码YUV的残差

	// clip for only non-zero cbf case
	if (pcCU->getQtRootCbf(0))  // 获得非全零CU的重建值
	{
		
		Int         iWidth;
		Int         iHeight;
		UInt        uiPartAddr;
		// 给当前CU中的每个PU复制前一帧EL图像的残差
		for (Int iPartIdx = 0; iPartIdx < pcCU->getNumPartitions(); iPartIdx++)
		{
			pcCU->getPartIndexAndSize(iPartIdx, uiPartAddr, iWidth, iHeight);
			// Int iRefIdx = pcCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx(uiPartAddr);    assert(iRefIdx >= 0);

			// TComMv     cMv = pcCU->getCUMvField(REF_PIC_LIST_0)->getMv(uiPartAddr);
			// pcCU->clipMv(cMv);

			TComPicYuv *refPicYuvResi = pcLastELCU->getPic()->getPicYuvResi();  // 前一帧EL图像的残差值
			// TComPicYuv *refPicYuvReco = pcLastELCU->getPic()->getPicYuvRec();  // 前一帧EL图像的重建值

			for (UInt comp = COMPONENT_Y; comp < m_ppcYuvReco[uiDepth]->getNumberValidComponents(); comp++)  // 分别复制3个分量的残差
			{
				const ComponentID compID = ComponentID(comp);
				UInt cxWidth = iWidth >> refPicYuvResi->getComponentScaleX(compID);
				UInt cxHeight = iHeight >> refPicYuvResi->getComponentScaleY(compID);
				Int     refStride = refPicYuvResi->getStride(compID);
				Int     dstStride = m_ppcYuvResi[uiDepth]->getStride(compID);
				Int shiftHor = (2 + refPicYuvResi->getComponentScaleX(compID));
				Int shiftVer = (2 + refPicYuvResi->getComponentScaleY(compID));

				// Int     refOffset = (cMv.getHor() >> shiftHor) + (cMv.getVer() >> shiftVer) * refStride;

				// Pel*    refYuvResi = refPicYuvResi->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr) + refOffset;  // 参考PU的残差值的原点像素
				// Pel*    refYuvReco = refPicYuvReco->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr) + refOffset;  // 参考PU的重建值的原点像素
				Pel*    refYuvResi = refPicYuvResi->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr);  // 同位PU的残差值的原点像素

				Pel*    dst = m_ppcYuvResi[uiDepth]->getAddr(compID, uiPartAddr);  // 当前PU残差值的原点像素
				//Pel*    dstReco = m_ppcYuvReco[uiDepth]->getAddr(compID, uiPartAddr);  // 当前PU重建值的原点像素

				for (Int y = cxHeight; y != 0; y--)
				{
					::memcpy(dst, refYuvResi, sizeof(Pel)*cxWidth);  // 复制残差值
					dst += dstStride;
					refYuvResi += refStride;

					//::memcpy(dstReco, refYuvReco, sizeof(Pel)*cxWidth);  // 复制重建值
					//dstReco += dstStride;
					//refYuvReco += refStride;
				}

			}

		}

		m_ppcYuvReco[uiDepth]->addClip(m_ppcYuvReco[uiDepth], m_ppcYuvResi[uiDepth], 0, pcCU->getWidth(0), pcCU->getSlice()->getSPS()->getBitDepths());
	}
	else  // CBF == 0，当前CU的变换参数都为0，残差也为0，YUV预测值就是最终的YUV重建值
	{
		m_ppcYuvReco[uiDepth]->copyPartToPartYuv(m_ppcYuvReco[uiDepth], 0, pcCU->getWidth(0), pcCU->getHeight(0));
	}

	xCopyToPic(pcCU, pcCU->getPic(), uiZorderIdx, uiDepth);  // 把CU的重建YUV和残差YUV复制到丢失图像
}

// 重建丢失EL图像中的帧间PU的YUV
Void TDecCu::ReconInterLostPic(TComDataCU* pcCU, UInt uiZorderIdx, UInt uiDepth)
{
	// Residual reconstruction
	m_ppcYuvResi[uiDepth]->clear();

	// 运动补偿获得YUV的预测值，m_ppcYuvReco[uiDepth]存储当前CU的YUV的预测值
	// inter prediction
	m_pcPrediction->motionCompensation(pcCU, m_ppcYuvReco[uiDepth]);

	// inter recon
	// xDecodeInterTexture(pcCU, uiDepth);  // 解码YUV的残差

    // clip for only non-zero cbf case
	if (pcCU->getQtRootCbf(0))  // 获得非全零CU的重建值
	{

		Int         iWidth;
		Int         iHeight;
		UInt        uiPartAddr;
		for (Int iPartIdx = 0; iPartIdx < pcCU->getNumPartitions(); iPartIdx++)  // 给当前CU中的每个PU复制残差
		{
			pcCU->getPartIndexAndSize(iPartIdx, uiPartAddr, iWidth, iHeight);
			Int iRefIdx = pcCU->getCUMvField(REF_PIC_LIST_0)->getRefIdx(uiPartAddr);    assert(iRefIdx >= 0);

			// 在解码端，如果当前EL图像丢失且前一帧EL图像没有丢失，因为BL图像的参考列表改变了，所以需要把BL图像中CU的参考索引转换为EL图像中CU的参考索引
			if (pcCU->getDecodedFlag() && pcCU->getPic()->getLayerId() > 0 && pcCU->getPic()->getSlice(0)->getAddParam()->getPicLossFlag(1)[pcCU->getPic()->getPOC()] == 0)
			{
				if (pcCU->getPic()->getPOC() > 0 && pcCU->getPic()->getSlice(0)->getAddParam()->getPicLossFlag(1)[pcCU->getPic()->getPOC() - 1] != 0)
				{
					Int iBLRefIdxToELRefIdx[4] = { 0, 0, 1, 2 };
					iRefIdx = iBLRefIdxToELRefIdx[iRefIdx];
				}
			}

			TComMv     cMv = pcCU->getCUMvField(REF_PIC_LIST_0)->getMv(uiPartAddr);
			pcCU->clipMv(cMv);

			TComPicYuv *refPicYuvResi = pcCU->getSlice()->getRefPic(REF_PIC_LIST_0, iRefIdx)->getPicYuvResi();  // 参考图像的残差值
			//TComPicYuv *refPicYuvReco = pcCU->getSlice()->getRefPic(REF_PIC_LIST_0, iRefIdx)->getPicYuvRec();  // 参考图像的重建值

			for (UInt comp = COMPONENT_Y; comp < m_ppcYuvReco[uiDepth]->getNumberValidComponents(); comp++)  // 分别复制3个分量的残差
			{
				const ComponentID compID = ComponentID(comp);
				UInt cxWidth = iWidth >> refPicYuvResi->getComponentScaleX(compID);
				UInt cxHeight = iHeight >> refPicYuvResi->getComponentScaleY(compID);
				Int     refStride = refPicYuvResi->getStride(compID);
				Int     dstStride = m_ppcYuvResi[uiDepth]->getStride(compID);
				Int shiftHor = (2 + refPicYuvResi->getComponentScaleX(compID));
				Int shiftVer = (2 + refPicYuvResi->getComponentScaleY(compID));

				Int     refOffset = (cMv.getHor() >> shiftHor) + (cMv.getVer() >> shiftVer) * refStride;

				Pel*    refYuvResi = refPicYuvResi->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr) + refOffset;  // 参考PU的残差值的原点像素
				//Pel*    refYuvReco = refPicYuvReco->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr) + refOffset;  // 参考PU的重建值的原点像素

		        /*// 同位BL图像中PU的残差值的原点像素
				const Pel*    ReplacedYuvResi = pcCU->getPic()->getPicYuvResi()->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr);

				for (Int y = cxHeight - 1; y >= 0; y--)
				{
					for (Int x = cxWidth - 1; x >= 0; x--)
					{
						refYuvResi[x] = refYuvResi[x] + ReplacedYuvResi[x];
					}
					refYuvResi += refStride;
					ReplacedYuvResi += refStride;
				}*/

				Pel*    dst = m_ppcYuvResi[uiDepth]->getAddr(compID, uiPartAddr);  // 当前PU残差值的原点像素
				//Pel*    dstReco = m_ppcYuvReco[uiDepth]->getAddr(compID, uiPartAddr);  // 当前PU重建值的原点像素

				for (Int y = cxHeight; y != 0; y--)
				{
					::memcpy(dst, refYuvResi, sizeof(Pel)*cxWidth);  // 复制残差值
					dst += dstStride;
					refYuvResi += refStride;
					/*
					::memcpy(dstReco, refYuvReco, sizeof(Pel)*cxWidth);  // 复制重建值
					dstReco += dstStride;
					refYuvReco += refStride;
					*/
				}

			}

		}

		m_ppcYuvReco[uiDepth]->addClip(m_ppcYuvReco[uiDepth], m_ppcYuvResi[uiDepth], 0, pcCU->getWidth(0), pcCU->getSlice()->getSPS()->getBitDepths());
	}
	else  // CBF == 0，当前CU的变换参数都为0，残差也为0，YUV预测值就是最终的YUV重建值
	{
		m_ppcYuvReco[uiDepth]->copyPartToPartYuv(m_ppcYuvReco[uiDepth], 0, pcCU->getWidth(0), pcCU->getHeight(0));
	}
	
	xCopyToPic(pcCU, pcCU->getPic(), uiZorderIdx, uiDepth);  // 把CU的重建YUV和残差YUV复制到丢失图像
}
#endif
#endif

Void
TDecCu::xIntraRecBlk(       TComYuv*    pcRecoYuv,
                            TComYuv*    pcPredYuv,
                            TComYuv*    pcResiYuv,
                      const ComponentID compID,
                            TComTU     &rTu)
{
  if (!rTu.ProcessComponentSection(compID))
  {
    return;
  }
  const Bool       bIsLuma = isLuma(compID);


  TComDataCU *pcCU = rTu.getCU();
  const TComSPS &sps=*(pcCU->getSlice()->getSPS());
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU();

  const TComRectangle &tuRect  =rTu.getRect(compID);
  const UInt uiWidth           = tuRect.width;
  const UInt uiHeight          = tuRect.height;
  const UInt uiStride          = pcRecoYuv->getStride (compID);
        Pel* piPred            = pcPredYuv->getAddr( compID, uiAbsPartIdx );
  const ChromaFormat chFmt     = rTu.GetChromaFormat();

  if (uiWidth != uiHeight)
  {
    //------------------------------------------------

    //split at current level if dividing into square sub-TUs

    TComTURecurse subTURecurse(rTu, false, TComTU::VERTICAL_SPLIT, true, compID);

    //recurse further
    do
    {
      xIntraRecBlk(pcRecoYuv, pcPredYuv, pcResiYuv, compID, subTURecurse);
    } while (subTURecurse.nextSection(rTu));

    //------------------------------------------------

    return;
  }

  const UInt uiChPredMode  = pcCU->getIntraDir( toChannelType(compID), uiAbsPartIdx );
  const UInt partsPerMinCU = 1<<(2*(sps.getMaxTotalCUDepth() - sps.getLog2DiffMaxMinCodingBlockSize()));
  const UInt uiChCodedMode = (uiChPredMode==DM_CHROMA_IDX && !bIsLuma) ? pcCU->getIntraDir(CHANNEL_TYPE_LUMA, getChromasCorrespondingPULumaIdx(uiAbsPartIdx, chFmt, partsPerMinCU)) : uiChPredMode;
  const UInt uiChFinalMode = ((chFmt == CHROMA_422)       && !bIsLuma) ? g_chroma422IntraAngleMappingTable[uiChCodedMode] : uiChCodedMode;

  //===== init availability pattern =====
  const Bool bUseFilteredPredictions=TComPrediction::filteringIntraReferenceSamples(compID, uiChFinalMode, uiWidth, uiHeight, chFmt, pcCU->getSlice()->getSPS()->getSpsRangeExtension().getIntraSmoothingDisabledFlag());

#if DEBUG_STRING
  std::ostream &ss(std::cout);
#endif

  DEBUG_STRING_NEW(sTemp)
  m_pcPrediction->initIntraPatternChType( rTu, compID, bUseFilteredPredictions  DEBUG_STRING_PASS_INTO(sTemp) );


  //===== get prediction signal =====

  m_pcPrediction->predIntraAng( compID,   uiChFinalMode, 0 /* Decoder does not have an original image */, 0, piPred, uiStride, rTu, bUseFilteredPredictions );

#if DEBUG_STRING
  ss << sTemp;
#endif

  //===== inverse transform =====
  Pel*      piResi            = pcResiYuv->getAddr( compID, uiAbsPartIdx );
  TCoeff*   pcCoeff           = pcCU->getCoeff(compID) + rTu.getCoefficientOffset(compID);//( uiNumCoeffInc * uiAbsPartIdx );

  const QpParam cQP(*pcCU, compID);


  DEBUG_STRING_NEW(sDebug);
#if DEBUG_STRING
  const Int debugPredModeMask=DebugStringGetPredModeMask(MODE_INTRA);
  std::string *psDebug=(DebugOptionList::DebugString_InvTran.getInt()&debugPredModeMask) ? &sDebug : 0;
#endif

  if (pcCU->getCbf(uiAbsPartIdx, compID, rTu.GetTransformDepthRel()) != 0)
  {
    m_pcTrQuant->invTransformNxN( rTu, compID, piResi, uiStride, pcCoeff, cQP DEBUG_STRING_PASS_INTO(psDebug) );
  }
  else
  {
    for (UInt y = 0; y < uiHeight; y++)
    {
      for (UInt x = 0; x < uiWidth; x++)
      {
        piResi[(y * uiStride) + x] = 0;
      }
    }
  }

#if DEBUG_STRING
  if (psDebug)
  {
    ss << (*psDebug);
  }
#endif

  //===== reconstruction =====
  const UInt uiRecIPredStride  = pcCU->getPic()->getPicYuvRec()->getStride(compID);

  const Bool useCrossComponentPrediction = isChroma(compID) && (pcCU->getCrossComponentPredictionAlpha(uiAbsPartIdx, compID) != 0);
  const Pel* pResiLuma  = pcResiYuv->getAddr( COMPONENT_Y, uiAbsPartIdx );
  const Int  strideLuma = pcResiYuv->getStride( COMPONENT_Y );

        Pel* pPred      = piPred;
        Pel* pResi      = piResi;
        Pel* pReco      = pcRecoYuv->getAddr( compID, uiAbsPartIdx );
        Pel* pRecIPred  = pcCU->getPic()->getPicYuvRec()->getAddr( compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiAbsPartIdx );


#if DEBUG_STRING
  const Bool bDebugPred=((DebugOptionList::DebugString_Pred.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
  const Bool bDebugResi=((DebugOptionList::DebugString_Resi.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
  const Bool bDebugReco=((DebugOptionList::DebugString_Reco.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
  if (bDebugPred || bDebugResi || bDebugReco)
  {
    ss << "###: " << "CompID: " << compID << " pred mode (ch/fin): " << uiChPredMode << "/" << uiChFinalMode << " absPartIdx: " << rTu.GetAbsPartIdxTU() << std::endl;
  }
#endif

  const Int clipbd = sps.getBitDepth(toChannelType(compID));
#if O0043_BEST_EFFORT_DECODING
  const Int bitDepthDelta = sps.getStreamBitDepth(toChannelType(compID)) - clipbd;
#endif

  if( useCrossComponentPrediction )
  {
    TComTrQuant::crossComponentPrediction( rTu, compID, pResiLuma, piResi, piResi, uiWidth, uiHeight, strideLuma, uiStride, uiStride, true );
  }

  for( UInt uiY = 0; uiY < uiHeight; uiY++ )
  {
#if DEBUG_STRING
    if (bDebugPred || bDebugResi || bDebugReco)
    {
      ss << "###: ";
    }

    if (bDebugPred)
    {
      ss << " - pred: ";
      for( UInt uiX = 0; uiX < uiWidth; uiX++ )
      {
        ss << pPred[ uiX ] << ", ";
      }
    }
    if (bDebugResi)
    {
      ss << " - resi: ";
    }
#endif

    for( UInt uiX = 0; uiX < uiWidth; uiX++ )
    {
#if DEBUG_STRING
      if (bDebugResi)
      {
        ss << pResi[ uiX ] << ", ";
      }
#endif
#if O0043_BEST_EFFORT_DECODING
      pReco    [ uiX ] = ClipBD( rightShiftEvenRounding<Pel>(pPred[ uiX ] + pResi[ uiX ], bitDepthDelta), clipbd );
#else
      pReco    [ uiX ] = ClipBD( pPred[ uiX ] + pResi[ uiX ], clipbd );
#endif
      pRecIPred[ uiX ] = pReco[ uiX ];
    }
#if DEBUG_STRING
    if (bDebugReco)
    {
      ss << " - reco: ";
      for( UInt uiX = 0; uiX < uiWidth; uiX++ )
      {
        ss << pReco[ uiX ] << ", ";
      }
    }

    if (bDebugPred || bDebugResi || bDebugReco)
    {
      ss << "\n";
    }
#endif
    pPred     += uiStride;
    pResi     += uiStride;
    pReco     += uiStride;
    pRecIPred += uiRecIPredStride;
  }
}


Void
TDecCu::xReconIntraQT( TComDataCU* pcCU, UInt uiDepth )
{
  if (pcCU->getIPCMFlag(0))
  {
    xReconPCM( pcCU, uiDepth );
    return;
  }
  const UInt numChType = pcCU->getPic()->getChromaFormat()!=CHROMA_400 ? 2 : 1;
  for (UInt chType=CHANNEL_TYPE_LUMA; chType<numChType; chType++)
  {
    const ChannelType chanType=ChannelType(chType);
    const Bool NxNPUHas4Parts = ::isChroma(chanType) ? enable4ChromaPUsInIntraNxNCU(pcCU->getPic()->getChromaFormat()) : true;
    const UInt uiInitTrDepth = ( pcCU->getPartitionSize(0) != SIZE_2Nx2N && NxNPUHas4Parts ? 1 : 0 );

    TComTURecurse tuRecurseCU(pcCU, 0);
    TComTURecurse tuRecurseWithPU(tuRecurseCU, false, (uiInitTrDepth==0)?TComTU::DONT_SPLIT : TComTU::QUAD_SPLIT);

    do
    {
      xIntraRecQT( m_ppcYuvReco[uiDepth], m_ppcYuvReco[uiDepth], m_ppcYuvResi[uiDepth], chanType, tuRecurseWithPU );
    } while (tuRecurseWithPU.nextSection(tuRecurseCU));
  }
}



/** Function for deriving reconstructed PU/CU chroma samples with QTree structure
 * \param pcRecoYuv pointer to reconstructed sample arrays
 * \param pcPredYuv pointer to prediction sample arrays
 * \param pcResiYuv pointer to residue sample arrays
 * \param chType    texture channel type (luma/chroma)
 * \param rTu       reference to transform data
 *
 \ This function derives reconstructed PU/CU chroma samples with QTree recursive structure
 */

Void
TDecCu::xIntraRecQT(TComYuv*    pcRecoYuv,
                    TComYuv*    pcPredYuv,
                    TComYuv*    pcResiYuv,
                    const ChannelType chType,
                    TComTU     &rTu)
{
  UInt uiTrDepth    = rTu.GetTransformDepthRel();
  TComDataCU *pcCU  = rTu.getCU();
  UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
  UInt uiTrMode     = pcCU->getTransformIdx( uiAbsPartIdx );
  if( uiTrMode == uiTrDepth )
  {
    if (isLuma(chType))
    {
      xIntraRecBlk( pcRecoYuv, pcPredYuv, pcResiYuv, COMPONENT_Y,  rTu );
    }
    else
    {
      const UInt numValidComp=getNumberValidComponents(rTu.GetChromaFormat());
      for(UInt compID=COMPONENT_Cb; compID<numValidComp; compID++)
      {
        xIntraRecBlk( pcRecoYuv, pcPredYuv, pcResiYuv, ComponentID(compID), rTu );
      }
    }
  }
  else
  {
    TComTURecurse tuRecurseChild(rTu, false);
    do
    {
      xIntraRecQT( pcRecoYuv, pcPredYuv, pcResiYuv, chType, tuRecurseChild );
    } while (tuRecurseChild.nextSection(rTu));
  }
}

Void TDecCu::xCopyToPic( TComDataCU* pcCU, TComPic* pcPic, UInt uiZorderIdx, UInt uiDepth )
{
  UInt uiCtuRsAddr = pcCU->getCtuRsAddr();

  m_ppcYuvReco[uiDepth]->copyToPicYuv( pcPic->getPicYuvRec (), uiCtuRsAddr, uiZorderIdx );  // 把CU的重建YUV复制到图像

#if DecodedELImage_Restoration
#if ImageRestorationRefinement

  // 把CU的重建YUV存储到图像的YUV残差变量中
  // m_ppcYuvReco[uiDepth]->copyToPicYuv(pcPic->getPicYuvResi(), uiCtuRsAddr, uiZorderIdx);

  
  // 存储图像的残差值
  if (pcCU->getQtRootCbf(0))
  {
	  m_ppcYuvResi[uiDepth]->copyToPicYuv(pcPic->getPicYuvResi(), uiCtuRsAddr, uiZorderIdx);  // 把CU的YUV残差复制到图像
  }
  else  // 没有残差时（Skip模式），把残差值设为0
  {
	  TComPicYuv* pcPicYuvDst = pcPic->getPicYuvResi();
	  for (Int comp = 0; comp < m_ppcYuvReco[uiDepth]->getNumberValidComponents(); comp++)
	  {
		  const ComponentID compID = ComponentID(comp);
		  const Int iWidth = m_ppcYuvReco[uiDepth]->getWidth(compID);
		  const Int iHeight = m_ppcYuvReco[uiDepth]->getHeight(compID);

		  const Pel* pSrc = m_ppcYuvResi[uiDepth]->getAddr(compID, 0, iWidth);  // 当前CU的原点像素
		  Pel* pDst = pcPicYuvDst->getAddr(compID, uiCtuRsAddr, uiZorderIdx);

		  const UInt  iSrcStride = m_ppcYuvReco[uiDepth]->getStride(compID);
		  const UInt  iDstStride = pcPicYuvDst->getStride(compID);

		  for (Int y = iHeight; y != 0; y--)
		  {   
			  ::memset(pDst, 0, sizeof(Pel)*iWidth);  // 把残差值设为0
			  //::memcpy(pDst, 0, sizeof(Pel)*iWidth);  // 把残差值设为0
			  pDst += iDstStride;
		  }
	  }
  }
  
#endif
#endif

  return;
}

Void TDecCu::xDecodeInterTexture ( TComDataCU* pcCU, UInt uiDepth )
{

  TComTURecurse tuRecur(pcCU, 0, uiDepth);

  for(UInt ch=0; ch<pcCU->getPic()->getNumberValidComponents(); ch++)
  {
    const ComponentID compID=ComponentID(ch);
    DEBUG_STRING_OUTPUT(std::cout, debug_reorder_data_inter_token[compID])

    m_pcTrQuant->invRecurTransformNxN ( compID, m_ppcYuvResi[uiDepth], tuRecur );  // TU进行递归反变换获得YUV的残差
  }

  DEBUG_STRING_OUTPUT(std::cout, debug_reorder_data_inter_token[MAX_NUM_COMPONENT])
}

/** Function for deriving reconstructed luma/chroma samples of a PCM mode CU.
 * \param pcCU pointer to current CU
 * \param uiPartIdx part index
 * \param piPCM pointer to PCM code arrays
 * \param piReco pointer to reconstructed sample arrays
 * \param uiStride stride of reconstructed sample arrays
 * \param uiWidth CU width
 * \param uiHeight CU height
 * \param compID colour component ID
 * \returns Void
 */
Void TDecCu::xDecodePCMTexture( TComDataCU* pcCU, const UInt uiPartIdx, const Pel *piPCM, Pel* piReco, const UInt uiStride, const UInt uiWidth, const UInt uiHeight, const ComponentID compID)
{
        Pel* piPicReco         = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu()+uiPartIdx);
  const UInt uiPicStride       = pcCU->getPic()->getPicYuvRec()->getStride(compID);
  const TComSPS &sps           = *(pcCU->getSlice()->getSPS());
  const UInt uiPcmLeftShiftBit = sps.getBitDepth(toChannelType(compID)) - sps.getPCMBitDepth(toChannelType(compID));

  for(UInt uiY = 0; uiY < uiHeight; uiY++ )
  {
    for(UInt uiX = 0; uiX < uiWidth; uiX++ )
    {
      piReco[uiX] = (piPCM[uiX] << uiPcmLeftShiftBit);
      piPicReco[uiX] = piReco[uiX];
    }
    piPCM += uiWidth;
    piReco += uiStride;
    piPicReco += uiPicStride;
  }
}

/** Function for reconstructing a PCM mode CU.
 * \param pcCU pointer to current CU
 * \param uiDepth CU Depth
 * \returns Void
 */
Void TDecCu::xReconPCM( TComDataCU* pcCU, UInt uiDepth )
{
  const UInt maxCuWidth     = pcCU->getSlice()->getSPS()->getMaxCUWidth();
  const UInt maxCuHeight    = pcCU->getSlice()->getSPS()->getMaxCUHeight();
  for (UInt ch=0; ch < pcCU->getPic()->getNumberValidComponents(); ch++)
  {
    const ComponentID compID = ComponentID(ch);
    const UInt width  = (maxCuWidth >>(uiDepth+m_ppcYuvResi[uiDepth]->getComponentScaleX(compID)));
    const UInt height = (maxCuHeight>>(uiDepth+m_ppcYuvResi[uiDepth]->getComponentScaleY(compID)));
    const UInt stride = m_ppcYuvResi[uiDepth]->getStride(compID);
    Pel * pPCMChannel = pcCU->getPCMSample(compID);
    Pel * pRecChannel = m_ppcYuvReco[uiDepth]->getAddr(compID);
    xDecodePCMTexture( pcCU, 0, pPCMChannel, pRecChannel, stride, width, height, compID );
  }
}

/** Function for filling the PCM buffer of a CU using its reconstructed sample array
 * \param pCU   pointer to current CU
 * \param depth CU Depth
 */
Void TDecCu::xFillPCMBuffer(TComDataCU* pCU, UInt depth)
{
  const ChromaFormat format = pCU->getPic()->getChromaFormat();
  const UInt numValidComp   = getNumberValidComponents(format);
  const UInt maxCuWidth     = pCU->getSlice()->getSPS()->getMaxCUWidth();
  const UInt maxCuHeight    = pCU->getSlice()->getSPS()->getMaxCUHeight();

  for (UInt componentIndex = 0; componentIndex < numValidComp; componentIndex++)
  {
    const ComponentID component = ComponentID(componentIndex);

    const UInt width  = maxCuWidth  >> (depth + getComponentScaleX(component, format));
    const UInt height = maxCuHeight >> (depth + getComponentScaleY(component, format));

    Pel *source      = m_ppcYuvReco[depth]->getAddr(component, 0, width);
    Pel *destination = pCU->getPCMSample(component);

    const UInt sourceStride = m_ppcYuvReco[depth]->getStride(component);

    for (Int line = 0; line < height; line++)
    {
      for (Int column = 0; column < width; column++)
      {
        destination[column] = source[column];
      }

      source      += sourceStride;
      destination += width;
    }
  }
}

//! \}
