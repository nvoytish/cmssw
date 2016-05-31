/** 
 * \file CSCSegAlgRU.cc 
 * 
 *  \authors V.Palichik & N.Voytishin   
 *  \some functions and structure taken from SK algo by M.Sani 
 */ 
	  
#include "CSCSegAlgoRU.h" 
#include "Geometry/CSCGeometry/interface/CSCLayer.h" 
// For clhep Matrix::solve 
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h" 
#include "DataFormats/GeometryVector/interface/GlobalPoint.h" 
 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 
 
#include <algorithm> 
#include <cmath> 
#include <iostream> 
#include <string> 
 
CSCSegAlgoRU::CSCSegAlgoRU(const edm::ParameterSet& ps) : CSCSegmentAlgorithm(ps), 
							  myName("CSCSegAlgoRU") { 
		 
  doCollisions = ps.getParameter<bool>("doCollisions"); 
  chi2_str_   = ps.getParameter<double>("chi2_str_"); 
  chi2Norm_2D_   = ps.getParameter<double>("chi2Norm_2D_");    
  dRMax       = ps.getParameter<double>("dRMax"); 
  dPhiMax        = ps.getParameter<double>("dPhiMax"); 
  dRIntMax   = ps.getParameter<double>("dRIntMax"); 
  dPhiIntMax    = ps.getParameter<double>("dPhiIntMax"); 
  chi2Max        = ps.getParameter<double>("chi2Max"); 
  wideSeg        = ps.getParameter<double>("wideSeg"); 
  minLayersApart = ps.getParameter<int>("minLayersApart"); 
		 
  LogDebug("CSC") << myName << " has algorithm cuts set to: \n" 
		  << "--------------------------------------------------------------------\n" 
		  << "dRMax     = " << dRMax << '\n' 
		  << "dPhiMax      = " << dPhiMax << '\n' 
		  << "dRIntMax = " << dRIntMax << '\n' 
		  << "dPhiIntMax  = " << dPhiIntMax << '\n' 
		  << "chi2Max      = " << chi2Max << '\n' 
		  << "wideSeg      = " << wideSeg << '\n' 
		  << "minLayersApart = " << minLayersApart << std::endl; 
  //reset the thresholds for non-collision data 
  if(!doCollisions){ 
    dRMax = 2.0; 
    dPhiMax = 2*dPhiMax; 
    dRIntMax = 2*dRIntMax; 
    dPhiIntMax = 2*dPhiIntMax; 
    chi2Norm_2D_ = 5*chi2Norm_2D_; 
    chi2_str_ = 100; 
    chi2Max = 2*chi2Max; 
  } 
} 
 
std::vector<CSCSegment> CSCSegAlgoRU::run(const CSCChamber* aChamber, const ChamberHitContainer& rechits){ 
  theChamber = aChamber;  
  return buildSegments(rechits);  
} 
 
std::vector<CSCSegment> CSCSegAlgoRU::buildSegments(const ChamberHitContainer& urechits) { 
	  
  ChamberHitContainer rechits = urechits; 
  LayerIndex layerIndex(rechits.size()); 
  int recHits_per_layer[6] = {0,0,0,0,0,0}; 

  //skip events with high multiplicity of hits 
  if (rechits.size()>150){ 
    return std::vector<CSCSegment>(); 
  } 
 
  int iadd = 0; 
 
  for(unsigned int i = 0; i < rechits.size(); i++) { 
    recHits_per_layer[rechits[i]->cscDetId().layer()-1]++;//count rh per chamber 
    layerIndex[i] = rechits[i]->cscDetId().layer(); 
  } 
 
  double z1 = theChamber->layer(1)->position().z(); 
  double z6 = theChamber->layer(6)->position().z(); 
	   
  if ( z1 > 0. ) { 
    if ( z1 > z6 ) {  
      reverse(layerIndex.begin(), layerIndex.end()); 
      reverse(rechits.begin(), rechits.end()); 
    }     
  } 
  else if ( z1 < 0. ) { 
    if ( z1 < z6 ) { 
      reverse(layerIndex.begin(), layerIndex.end()); 
      reverse(rechits.begin(), rechits.end()); 
    }     
  } 
	   
  if (rechits.size() < 2) { 
    return std::vector<CSCSegment>();  
  } 
		 
  //-----coments to be reviewed and updated!!! 
 
  // We have at least 2 hits. We intend to try all possible pairs of hits to start  
  // segment building. 'All possible' means each hit lies on different layers in the chamber. 
  // after all same size segs are build we get rid of the overcrossed segments using the chi2 criteria 
  // the hits from the segs that are left are marked as used and are not added to segs in future iterations 
 
  // the hits from 3p segs are marked as used separately in order to try to assamble them in longer segments  
  // in case there is a second pass    
	   
  // Choose first hit (as close to IP as possible) h1 and second hit 
  // (as far from IP as possible) h2 To do this we iterate over hits 
  // in the chamber by layer - pick two layers.  Then we 
  // iterate over hits within each of these layers and pick h1 and h2 
  // these.  If they are 'close enough' we build an empty 
  // segment.  Then try adding hits to this segment. 
	   
  // Initialize flags that a given hit has been allocated to a segment 
  BoolContainer used(rechits.size(), false); 
  BoolContainer used3p(rechits.size(), false); 
 	   
  // Define buffer for segments we build  
  std::vector<CSCSegment> segments; 
 
	   
  ChamberHitContainerCIt ib = rechits.begin(); 
  ChamberHitContainerCIt ie = rechits.end(); 
	   
  // Possibly allow 2 passes, second widening scale factor for cuts 
  windowScale = 1.; // scale factor for cuts 
	   
  int npass = (wideSeg > 1.)? 2 : 1; 
	   
  for (int ipass = 0; ipass < npass; ++ipass) { 
 
    if(windowScale >1.){ 
      iadd = 1; 
      strip_iadd = 2; 
      chi2D_iadd = 2; 
    } 
    for(int n_seg_min = 6; n_seg_min > 2 + iadd; --n_seg_min){ 
 
      //define buffer for overcrossed seg 
 
      BoolContainer common_used(rechits.size(),false); 
      std::array<BoolContainer, 120> common_used_it = {}; 
      for (unsigned int i = 0; i < common_used_it.size(); i++) { 
	common_used_it[i] = common_used; 
      } 
 
      ChamberHitContainer best_proto_segment[120]; 
      float min_chi[120] = {9999}; 
      int common_it = 0; 
      bool first_proto_segment = true; 
      for (ChamberHitContainerCIt i1 = ib; i1 != ie; ++i1) { 
	bool segok = false; 
	//skip if rh is used and the layer that has big rh multiplicity(>25RHs) 
	if(used[i1-ib] || recHits_per_layer[int(layerIndex[i1-ib])-1]>25) continue; 
		 
	int layer1 = layerIndex[i1-ib];  
	const CSCRecHit2D* h1 = *i1; 
 
	for (ChamberHitContainerCIt i2 = ie-1; i2 != i1; --i2) { 
 
	  if(used[i2-ib] || recHits_per_layer[int(layerIndex[i2-ib])-1]>25) continue; 
			 
	  int layer2 = layerIndex[i2-ib];  
 
	  if((abs(layer2 - layer1) + 1) < n_seg_min) break;//decrease n_seg_min 
		  
	  const CSCRecHit2D* h2 = *i2; 
		  
	  if (areHitsCloseInR(h1, h2) && areHitsCloseInGlobalPhi(h1, h2)) { 
			   
	    proto_segment.clear(); 
			  
	    if (!addHit(h1, layer1))continue; 
	    if (!addHit(h2, layer2))continue; 
			   
	    tryAddingHitsToSegment(rechits, used, layerIndex, i1, i2);  
	    // calculate error matrix... 
	    AlgebraicSymMatrix errors = calculateError(); 
 
	    // but reorder components to match what's required by TrackingRecHit interface 
	    // i.e. slopes first, then positions 
 
	    flipErrors( errors ); 
	    fitSlopes(); 
	    fillChiSquared(); 
	    segok = isSegmentGood(rechits); 
	    if (segok) { 
	      //run optimal point rejection method if the segment is long enough and has a big chi2   
	      if(abs(layer2-layer1) + 1 > n_seg_min) baseline(n_seg_min); 
 
	      // Copy the proto_segment and its properties inside a CSCSegment. 
	      // Then fill the segment vector.. 
	      // calculate error matrix... 
	      errors = calculateError(); 
 
	      // but reorder components to match what's required by TrackingRecHit interface 
	      // i.e. slopes first, then positions 
 
	      flipErrors( errors ); 
	      fitSlopes(); 
			  
	      fillChiSquared(); 
 
	      //calculate chi2_str 
 
 
	      if((theChi2/(2*proto_segment.size() - 4)) > chi2Norm_2D_*chi2D_iadd){ 
		proto_segment.clear(); 
	      } 
	      if (!proto_segment.empty()) { 
		// calculate error matrix... 
		AlgebraicSymMatrix errors = calculateError(); 
 
		// but reorder components to match what's required by TrackingRecHit interface 
		// i.e. slopes first, then positions 
 
		flipErrors( errors ); 
		fitSlopes(); 
		fillChiSquared(); 
		updateParameters(); 
					 
		//add same-size overcrossed protosegments to the collection 
		if(first_proto_segment){ 
		  flagHitsAsUsed(rechits, common_used_it[0]); 
		  min_chi[0] = (theChi2/(2*proto_segment.size() - 4)); 
		  best_proto_segment[0] = proto_segment; 
		  first_proto_segment = false; 
		}else{  //for the rest of found proto_segments  
		  common_it++; 
		  flagHitsAsUsed(rechits, common_used_it[common_it]); 
		  min_chi[common_it] = (theChi2/(2*proto_segment.size() - 4)); 
		  best_proto_segment[common_it] = proto_segment; 
		   
		  ChamberHitContainerCIt hi, iu, ik; 
		  int iter = common_it; 
		  for(iu = ib; iu != ie; ++iu) { 
		    for(hi = proto_segment.begin(); hi != proto_segment.end(); ++hi) { 
		     
		      if(*hi == *iu) { 
		 
			int merge_nr = -1; 
			for(int k = 0; k < iter+1; k++){ 
			  if(common_used_it[k][iu-ib] == true){ 
			    if(merge_nr != -1){ 
			      //merge the k and merge_nr collections of flaged hits into the merge_nr collection and unmark the k collection hits                                                    
			      for(ik = ib; ik != ie; ++ik) { 
				if(common_used_it[k][ik-ib] == true){ 
				  common_used_it[merge_nr][ik-ib] = true; 
				  common_used_it[k][ik-ib] = false; 
				} 
			      } 
			      //change best_protoseg if min_chi_2 is smaller                                                                                                                       
			      if(min_chi[k] < min_chi[merge_nr]){ 
				min_chi[merge_nr] = min_chi[k]; 
				best_proto_segment[merge_nr] = best_proto_segment[k]; 
				best_proto_segment[k].clear(); 
				min_chi[k] = 9999; 
			      } 
			      common_it--; 
			    } 
			    else{ 
			      merge_nr = k; 
			    } 
 
			  }//if(common_used[k][iu-ib] == true)                                                                                                                      
			}//for k                                                                                                                                                                  
		      }//if 
		    }//for proto_seg                                                                                                                                                                 
		  }//for rec_hits   
		} 
	      } 
	    } 
	  }  //   h1 & h2 close 
			 
	  if (segok)  
	    break; 
	}  //  i2 
      }  //  i1 
	    
      //add the reconstructed segments 
 
 
      for(int j = 0;j < common_it+1; j++){ 
	proto_segment = best_proto_segment[j]; 
	best_proto_segment[j].clear(); 
	// calculate error matrix...                                                                                                                                                  
	AlgebraicSymMatrix errors = calculateError();                                                                                                                                 
	// but reorder components to match what's required by TrackingRecHit interface                                                                                                
	// i.e. slopes first, then positions                                                                                                                                          
	flipErrors( errors );                                                                                                                                                         
	fitSlopes();  
	fillChiSquared(); 
	updateParameters();                                                                                                                                                           
	//SKIP empty proto-segments                                                                                                                                                                                                                                                                                                                                    
	if(proto_segment.size() == 0) continue;                                                                                                                                                                                                                                                                                                                          
	CSCSegment temp(proto_segment, theOrigin, theDirection, errors, theChi2);                                                                                                      
 
	segments.push_back(temp);                         
 
	//if the segment has 3 hits flag them as used in a particular way 
 
	if(proto_segment.size() == 3){ 
	  flagHitsAsUsed(rechits, used3p); 
	} 
	else{ 
	  flagHitsAsUsed(rechits, used);   
 
	} 
	proto_segment.clear(); 
      }  
 
    }//for n_seg_min 
 
    std::vector<CSCSegment>::iterator it =segments.begin(); 
    bool good_segs = false; 
    while(it != segments.end()) { 
      if ((*it).nRecHits() > 3){ 
	good_segs = true; 
	break; 
      } 
      ++it;	 
    }     
 
	    
    if (good_segs) break;  // only change window if not enough good segments were found (bool can be changed to int if a >0 number of good segs is required) 
		 
    // Increase cut windows by factor of wideSeg only for collisions 
    if(!doCollisions) break;     
    windowScale = wideSeg; 
 
  }  //  ipass 
  int used_rh = 0; 
  for (ChamberHitContainerCIt i1 = ib; i1 != ie; ++i1) { 
    if(used[i1-ib])used_rh++; 
  }  
 
 
 
  //cycle for displaced mu 
  windowScale = 1.; // scale factor for cuts 
  if(doCollisions && int(rechits.size()-used_rh)>2){//check if there are enough recHits left to build a segment 
    doCollisions = false; 
 
    dRMax = 2.0; 
    dPhiMax = 2*dPhiMax; 
    dRIntMax = 2*dRIntMax; 
    dPhiIntMax = 2*dPhiIntMax; 
    chi2Norm_2D_ = 5*chi2Norm_2D_; 
    chi2_str_ = 100; 
    chi2Max = 2*chi2Max; 
 
		 
    for(int n_seg_min = 6; n_seg_min > 2 + iadd; --n_seg_min){ 
      BoolContainer common_used(rechits.size(),false); 
      std::array<BoolContainer, 120> common_used_it = {}; 
      for (unsigned int i = 0; i < common_used_it.size(); i++) { 
	common_used_it[i] = common_used; 
      } 
      ChamberHitContainer best_proto_segment[120]; 
      float min_chi[120] = {9999}; 
      int common_it = 0; 
      bool first_proto_segment = true; 
      for (ChamberHitContainerCIt i1 = ib; i1 != ie; ++i1) { 
	bool segok = false; 
	//skip if rh is used and the layer tat has big rh multiplicity(>20RHs) 
	if(used[i1-ib] || recHits_per_layer[int(layerIndex[i1-ib])-1]>25) continue; 
 
	int layer1 = layerIndex[i1-ib];  
	const CSCRecHit2D* h1 = *i1; 
 
	for (ChamberHitContainerCIt i2 = ie-1; i2 != i1; --i2) { 
 
	  if(used[i2-ib] || recHits_per_layer[int(layerIndex[i2-ib])-1]>25) continue; 
	   
	  int layer2 = layerIndex[i2-ib]; 
 
	  if((abs(layer2 - layer1) + 1) < n_seg_min) break;//decrease n_seg_min 
		  
	  const CSCRecHit2D* h2 = *i2; 
	  if (areHitsCloseInR(h1, h2) && areHitsCloseInGlobalPhi(h1, h2)) { 
	    proto_segment.clear(); 
	    if (!addHit(h1, layer1))continue; 
	    if (!addHit(h2, layer2))continue; 
			   
	    tryAddingHitsToSegment(rechits, used, layerIndex, i1, i2);  
	    // calculate error matrix... 
	    AlgebraicSymMatrix errors = calculateError(); 
 
	    // but reorder components to match what's required by TrackingRecHit interface 
	    // i.e. slopes first, then positions 
 
	    flipErrors( errors ); 
	    fitSlopes(); 
	    fillChiSquared(); 
	    //	  float iadd_chi2 = 1; 
	    // Check no. of hits on segment, and if enough flag them as used 
	    // and store the segment 
	    segok = isSegmentGood(rechits); 
	    if (segok) { 
	      if(abs(layer2-layer1) + 1 > n_seg_min) baseline(n_seg_min); 
	      // Copy the proto_segment and its properties inside a CSCSegment. 
	      // Then fill the segment vector.. 
	      // calculate error matrix... 
	      errors = calculateError(); 
 
	      // but reorder components to match what's required by TrackingRecHit interface 
	      // i.e. slopes first, then positions 
 
	      flipErrors( errors ); 
	      fitSlopes(); 
	      fillChiSquared(); 
 
 
	      if((theChi2/(2*proto_segment.size() - 4)) > chi2Norm_2D_*chi2D_iadd){ 
		proto_segment.clear(); 
	      } 
	      if (!proto_segment.empty()) { 
		// calculate error matrix... 
		AlgebraicSymMatrix errors = calculateError(); 
 
		// but reorder components to match what's required by TrackingRecHit interface 
		// i.e. slopes first, then positions 
 
		flipErrors( errors ); 
		fitSlopes(); 
		fillChiSquared(); 
		updateParameters(); 
				  
		//add same-size overcrossed protosegments to the collection 
		if(first_proto_segment){ 
		  flagHitsAsUsed(rechits, common_used_it[0]); 
		  min_chi[0] = (theChi2/(2*proto_segment.size() - 4)); 
		  best_proto_segment[0] = proto_segment; 
		  first_proto_segment = false; 
		}else{  //for the rest of found proto_segments  
		  common_it++; 
		  flagHitsAsUsed(rechits, common_used_it[common_it]); 
		   
		  min_chi[common_it] = (theChi2/(2*proto_segment.size() - 4)); 
		  best_proto_segment[common_it] = proto_segment; 
		  ChamberHitContainerCIt hi, iu, ik; 
		  int iter = common_it; 
		  for(iu = ib; iu != ie; ++iu) { 
		    for(hi = proto_segment.begin(); hi != proto_segment.end(); ++hi) { 
		      if(*hi == *iu) { 
			int merge_nr = -1; 
			for(int k = 0; k < iter+1; k++){ 
			  if(common_used_it[k][iu-ib] == true){ 
			    if(merge_nr != -1){ 
			      //merge the k and merge_nr collections of flaged hits into the merge_nr collection and unmark the k collection hits                                   
			      for(ik = ib; ik != ie; ++ik) { 
				if(common_used_it[k][ik-ib] == true){ 
				  common_used_it[merge_nr][ik-ib] = true; 
				  common_used_it[k][ik-ib] = false; 
				} 
			      } 
			      //change best_protoseg if min_chi_2 is smaller                                                                                                                       
			      if(min_chi[k] < min_chi[merge_nr]){ 
				min_chi[merge_nr] = min_chi[k]; 
				best_proto_segment[merge_nr] = best_proto_segment[k]; 
				best_proto_segment[k].clear(); 
				min_chi[k] = 9999; 
			      } 
			      common_it--; 
			    } 
			    else{ 
			      merge_nr = k; 
			    } 
			  }//if(common_used[k][iu-ib] == true)                                                                                                                                     
			}//for k                                                                                                                                                                    
		      }//if                                                                                                                                                                        
		    }//for proto_seg                                                                                                                                                                 
		  }//for rec_hits   
		} 
	      } 
	    } 
	  }  //   h1 & h2 close 
			 
	  if (segok)  
	    break; 
	}  //  i2 
      }  //  i1 
	    
      //add the reconstructed segments 
 
      for(int j = 0;j < common_it+1; j++){ 
	proto_segment = best_proto_segment[j]; 
	best_proto_segment[j].clear(); 
	// calculate error matrix...                                                                                                                                                  
	AlgebraicSymMatrix errors = calculateError();                                                                                                                                 
	// but reorder components to match what's required by TrackingRecHit interface                                                                                                
	// i.e. slopes first, then positions                                                                                                                                          
	flipErrors( errors );                                                                                                                                                         
	fitSlopes();  
	fillChiSquared(); 
	updateParameters();                                                                                                                                                           
	//SKIP empty proto-segments                                                                                                                                                                                                                                                                                                                                    
	if(proto_segment.size() == 0) continue;                                                                                                                                                                                                                                                                                                                        
	CSCSegment temp(proto_segment, theOrigin, theDirection, errors, theChi2); 
	segments.push_back(temp);                         
	//if the segment has 3 hits flag them as used in a particular way 
	if(proto_segment.size() == 3){ 
	  flagHitsAsUsed(rechits, used3p); 
	} 
	else{ 
	  flagHitsAsUsed(rechits, used);   
	} 
	proto_segment.clear(); 
      }  
    }//for n_seg_min 
		 
 
    //reset the docollisions flag for the next chamber 
    doCollisions = true; 
    dRMax = 2.0; 
    dPhiMax = dPhiMax/2; 
    dRIntMax = dRIntMax/2; 
    dPhiIntMax = dPhiIntMax/2; 
    chi2Norm_2D_ = chi2Norm_2D_/5; 
    chi2_str_ = 100; 
    chi2Max = chi2Max/2; 
 
 
  }//if(!docollisions - end cycle for displaced mu 
 
  //get rid of enchansed 3p segments 
  std::vector<CSCSegment>::iterator it =segments.begin(); 
  while(it != segments.end()) { 
    if((*it).nRecHits() == 3){ 
      bool found_common = false; 
      std::vector<CSCRecHit2D> theseRH = (*it).specificRecHits(); 
      for (ChamberHitContainerCIt i1 = ib; i1 != ie; ++i1) { 
	if(used[i1-ib] && used3p[i1-ib]){ 
	 
	  const CSCRecHit2D* sh1 = *i1; 
		   
	  CSCDetId id = (CSCDetId)(*sh1).cscDetId(); 
	  int sh1layer = id.layer(); 
	  int RH_centerid     =  sh1->nStrips()/2; 
	  int RH_centerStrip =  sh1->channels(RH_centerid); 
	  int RH_wg = sh1->hitWire(); 
	  std::vector<CSCRecHit2D>::const_iterator sh; 
	  for(sh = theseRH.begin(); sh != theseRH.end(); ++sh){ 
	    //find segment hit coord 
	    CSCDetId idRH = (CSCDetId)(*sh).cscDetId(); 
	    int shlayer = idRH.layer(); 
	    int SegRH_centerid     =  sh->nStrips()/2; 
	    int SegRH_centerStrip =  sh->channels(SegRH_centerid); 
	    int SegRH_wg = sh->hitWire(); 
	    if(sh1layer == shlayer && SegRH_centerStrip == RH_centerStrip && SegRH_wg == RH_wg){ 
	      //remove the enchansed 3p segment 
	      segments.erase(it,(it+1)); 
	      found_common = true; 
	      break; 
	    } 
	  }//theserh 
	} 
	if(found_common) break;//current seg has already been erased 
      }//camber hits 
      if(!found_common)++it;  
    }//its a 3p seg 
    else{ 
      ++it;//go to the next seg 
    }  
  }//while 
  	   
  // Give the segments to the CSCChamber 
  return segments; 
}//build segments 
 
 
void CSCSegAlgoRU::tryAddingHitsToSegment(const ChamberHitContainer& rechits,  
					  const BoolContainer& used, const LayerIndex& layerIndex, 
					  const ChamberHitContainerCIt i1, const ChamberHitContainerCIt i2) { 
 
  // Iterate over the layers with hits in the chamber 
  // Skip the layers containing the segment endpoints 
  // Test each hit on the other layers to see if it is near the segment 
  // If it is, see whether there is already a hit on the segment from the same layer 
  //    - if so, and there are more than 2 hits on the segment, copy the segment, 
  //      replace the old hit with the new hit. If the new segment chi2 is better 
  //      then replace the original segment with the new one (by swap) 
  //    - if not, copy the segment, add the hit. If the new chi2/dof is still satisfactory 
  //      then replace the original segment with the new one (by swap) 
	   
  ChamberHitContainerCIt ib = rechits.begin(); 
  ChamberHitContainerCIt ie = rechits.end(); 
	   
  for (ChamberHitContainerCIt i = ib; i != ie; ++i) { 
 
    //continue; 
    if(layerIndex[i1-ib]<layerIndex[i2-ib]){ 
      if (layerIndex[i-ib] <= layerIndex[i1-ib] || layerIndex[i-ib] >= layerIndex[i2-ib] || i  == i1 || i == i2 || used[i-ib]){  
	if ( i  == i1 || i == i2 || used[i-ib]) 
	  continue;  
      } 
    } 
    else{ 
      if (layerIndex[i-ib] >= layerIndex[i1-ib] || layerIndex[i-ib] <= layerIndex[i2-ib] || i  == i1 || i == i2 || used[i-ib]){ 
	if ( i  == i1 || i == i2 || used[i-ib])                                                                                                              
	  continue; 
      }  
 
    } 
    int layer = layerIndex[i-ib];  
    const CSCRecHit2D* h = *i; 
    if (isHitNearSegment(h)) { 
      // Don't consider alternate hits on layers holding the two starting points 
      if (hasHitOnLayer(layer)) { 	   
	if (proto_segment.size() <= 2)continue; 
	compareProtoSegment(h, layer); 
      }  
      else{ 
	increaseProtoSegment(h, layer, chi2D_iadd);  
      } 
    }   // h & seg close 
  }   // i 
} 
 
bool CSCSegAlgoRU::areHitsCloseInR(const CSCRecHit2D* h1, const CSCRecHit2D* h2) const { 
  float maxWG_width[10] = {0, 0, 4.1, 5.69, 2.49, 5.06, 2.49, 5.06, 1.87, 5.06}; 
  int iStn = stationNumber(h1); 
		 
  //find maxWG_width for ME11 (tilt = 29deg) 
  int wg_num = h2->hitWire(); 
  if(iStn == 0 || iStn == 1){ 
    if (wg_num == 1){  
      maxWG_width[0] = 9.25; 
      maxWG_width[1] = 9.25;	     
    } 
    if (wg_num > 1 && wg_num < 48){  
      maxWG_width[0] = 3.14; 
      maxWG_width[1] = 3.14; 
    } 
    if (wg_num == 48){  
      maxWG_width[0] = 10.75; 
      maxWG_width[1] = 10.75; 
    } 
  } 
  const CSCLayer* l1 = theChamber->layer(h1->cscDetId().layer()); 
  GlobalPoint gp1 = l1->toGlobal(h1->localPosition());	 
  const CSCLayer* l2 = theChamber->layer(h2->cscDetId().layer()); 
  GlobalPoint gp2 = l2->toGlobal(h2->localPosition());	 
  //find z to understand the direction 
  float h1z = gp1.z(); 
  float h2z = gp2.z(); 
	   
  //switch off the IP check for non collisions case 
  if (!doCollisions){ 
    h1z = 1; 
    h2z = 1; 
  } 
  if (gp2.perp() > ((gp1.perp() - dRMax*maxWG_width[iStn])*h2z)/h1z && gp2.perp() < ((gp1.perp() + dRMax*maxWG_width[iStn])*h2z)/h1z){ 
    return true; 
  } 
  else{ 
    return false;    
  } 
} 
 
bool CSCSegAlgoRU::areHitsCloseInGlobalPhi(const CSCRecHit2D* h1, const CSCRecHit2D* h2) const { 
	   
  float strip_width[10] = {0.003878509, 0.002958185, 0.002327105, 0.00152552, 0.00465421, 0.002327105, 0.00465421, 0.002327105, 0.00465421, 0.002327105};//in rad 
 
  const CSCLayer* l1 = theChamber->layer(h1->cscDetId().layer()); 
  GlobalPoint gp1 = l1->toGlobal(h1->localPosition());	 
  const CSCLayer* l2 = theChamber->layer(h2->cscDetId().layer()); 
  GlobalPoint gp2 = l2->toGlobal(h2->localPosition());	 
  float err_stpos_h1 = h1->errorWithinStrip(); 
  float err_stpos_h2 = h2->errorWithinStrip(); 
	  
  int iStn = stationNumber(h1); 
	   
  float dphi_incr = 0; 
  if(err_stpos_h1>0.25*strip_width[iStn] || err_stpos_h2>0.25*strip_width[iStn])dphi_incr = 0.5*strip_width[iStn]; 
  float h1p = gp1.phi(); 
  float h2p = gp2.phi(); 
  float dphi12 = h1p - h2p; 
	   
  // Into range [-pi, pi) (phi() returns values in this range) 
  if (dphi12 < -M_PI)  
    dphi12 += 2.*M_PI;   
  if (dphi12 >  M_PI)  
    dphi12 -= 2.*M_PI; 
 
  return (fabs(dphi12) < (dPhiMax*strip_iadd+dphi_incr))? true:false;  // +v 
} 
 
bool CSCSegAlgoRU::isHitNearSegment(const CSCRecHit2D* h) const { 
	   
  // Is hit near segment?  
  // Requires deltaphi and deltaR within ranges specified in parameter set. 
  // Note that to make intuitive cuts on delta(phi) one must work in 
  // phi range (-pi, +pi] not [0, 2pi) 
 
  float strip_width[10] = {0.003878509, 0.002958185, 0.002327105, 0.00152552, 0.00465421, 0.002327105, 0.00465421, 0.002327105, 0.00465421, 0.002327105};//in rad 
 
  const CSCLayer* l1 = theChamber->layer((*(proto_segment.begin()))->cscDetId().layer()); 
  GlobalPoint gp1 = l1->toGlobal((*(proto_segment.begin()))->localPosition()); 
  const CSCLayer* l2 = theChamber->layer((*(proto_segment.begin()+1))->cscDetId().layer()); 
  GlobalPoint gp2 = l2->toGlobal((*(proto_segment.begin()+1))->localPosition()); 
  float err_stpos_h1 = (*(proto_segment.begin()))->errorWithinStrip(); 
  float err_stpos_h2 = (*(proto_segment.begin()+1))->errorWithinStrip(); 
 
  const CSCLayer* l = theChamber->layer(h->cscDetId().layer()); 
  GlobalPoint hp = l->toGlobal(h->localPosition()); 
  float err_stpos_h = h->errorWithinStrip(); 
 
  float hphi = hp.phi();          // in (-pi, +pi] 
  if (hphi < 0.) 
    hphi += 2.*M_PI;            // into range [0, 2pi) 
  float sphi = phiAtZ(hp.z());    // in [0, 2*pi) 
  float phidif = sphi-hphi; 
  if (phidif < 0.) 
    phidif += 2.*M_PI;          // into range [0, 2pi) 
  if (phidif > M_PI) 
    phidif -= 2.*M_PI;          // into range (-pi, pi] 
  
  
  CLHEP::HepMatrix  r_glob(6,1); 
     r_glob(1,1) = -99.99; 
     r_glob(2,1) = -99.99; 
     r_glob(3,1) = -99.99; 
     r_glob(4,1) = -99.99; 
     r_glob(5,1) = -99.99; 
     r_glob(6,1) = -99.99; 

       int iStn = stationNumber(h);
       float dphi_incr = 0; 
	 float pos_str = 1; 

     //increase dPhi cut if the hit is on the edge of the strip                                                                                                                                 
     float stpos = (*h).positionWithinStrip(); 
     bool centr_str = false; 
     if(iStn != 0 && iStn != 1){ 
       if (stpos > -0.25 && stpos < 0.25) centr_str = true; 
     } 
     if(err_stpos_h1<0.25*strip_width[iStn] || err_stpos_h2<0.25*strip_width[iStn] || err_stpos_h < 0.25*strip_width[iStn]){ 
      dphi_incr = 0.5*strip_width[iStn]; 
       }else{if(centr_str) pos_str = 1.3;} 
     
  r_glob((*(proto_segment.begin()))->cscDetId().layer(),1) = gp1.perp(); 
  r_glob((*(proto_segment.begin()+1))->cscDetId().layer(),1) = gp2.perp(); 
 
   
  float R =  hp.perp(); 
  int layer = h->cscDetId().layer(); 
  float r_interpolated = fit_r_phi(r_glob,layer); 
  float dr = fabs(r_interpolated - R); 
  
 
  float maxWG_width[10] = {0, 0, 4.1, 5.69, 2.49, 5.06, 2.49, 5.06, 1.87, 5.06}; 
 
   //find maxWG_width for ME11 (tilt = 29deg) 
    int wg_num = h->hitWire(); 
    if(iStn == 0 || iStn == 1){ 
      if (wg_num == 1){  
	maxWG_width[0] = 9.25; 
	maxWG_width[1] = 9.25;	     
      } 
      if (wg_num > 1 && wg_num < 48){  
	maxWG_width[0] = 3.14; 
	maxWG_width[1] = 3.14; 
      } 
      if (wg_num == 48){  
	maxWG_width[0] = 10.75; 
	maxWG_width[1] = 10.75; 
      } 
    } 
   
  return (fabs(phidif) <  dPhiIntMax*strip_iadd*pos_str+dphi_incr && fabs(dr) <  dRIntMax*maxWG_width[iStn])? true:false; 
  
} 
 
float CSCSegAlgoRU::phiAtZ(float z) const { 
   
  // Returns a phi in [ 0, 2*pi ) 
  const CSCLayer* l1 = theChamber->layer((*(proto_segment.begin()))->cscDetId().layer()); 
  GlobalPoint gp = l1->toGlobal(theOrigin);	 
  GlobalVector gv = l1->toGlobal(theDirection);	 
  float x = gp.x() + (gv.x()/gv.z())*(z - gp.z()); 
  float y = gp.y() + (gv.y()/gv.z())*(z - gp.z()); 
  float phi = atan2(y, x); 
  if (phi < 0.f )  
    phi += 2. * M_PI; 
   
  return phi ; 
} 
 
 
bool CSCSegAlgoRU::isSegmentGood(const ChamberHitContainer& rechitsInChamber) const { 
	   
  // If the chamber has 20 hits or fewer, require at least 3 hits on segment 
  // If the chamber has >20 hits require at least 4 hits 
  //@@ THESE VALUES SHOULD BECOME PARAMETERS? 
  bool ok = false; 
	   
  unsigned int iadd = ( rechitsInChamber.size() > 20)?  1 : 0;   
	   
  if (windowScale > 1.) 
    iadd = 1; 
	   
  if (proto_segment.size() >= 3+iadd) 
    ok = true; 
  return ok; 
} 
 
float CSCSegAlgoRU::hitPosInStrips(const CSCRecHit2D* h) const { 
  CSCDetId idRH = (CSCDetId)(*h).cscDetId(); 
  int kRing    = idRH.ring(); 
  int kStation = idRH.station(); 
  int kLayer = idRH.layer(); 
  // Find the strip containing this hit 
  int centerid     =  h->nStrips()/2; 
  int centerStrip =  h->channels(centerid); 
  float stpos = (*h).positionWithinStrip(); 
  float sp = -1;  
  // Take into account half-strip staggering of layers (ME1/1 has no staggering) 
  if (kStation == 1 && (kRing == 1 || kRing == 4)) sp = stpos + centerStrip; 
  else{ 
    if (kLayer == 1 || kLayer == 3 || kLayer == 5) sp = stpos + centerStrip; 
    if (kLayer == 2 || kLayer == 4 || kLayer == 6) sp = stpos - 0.5 + centerStrip; 
  } 
  return sp;   
} 
 
void CSCSegAlgoRU::flagHitsAsUsed(const ChamberHitContainer& rechitsInChamber,  
				  BoolContainer& used ) const { 
	   
  // Flag hits on segment as used 
  ChamberHitContainerCIt ib = rechitsInChamber.begin(); 
  ChamberHitContainerCIt hi, iu; 
	   
  for(hi = proto_segment.begin(); hi != proto_segment.end(); ++hi) { 
    for(iu = ib; iu != rechitsInChamber.end(); ++iu) { 
      if(*hi == *iu) 
	used[iu-ib] = true; 
    } 
  } 
} 
 
bool CSCSegAlgoRU::addHit(const CSCRecHit2D* aHit, int layer) { 
	   
  // Return true if hit was added successfully  
  // (and then parameters are updated). 
  // Return false if there is already a hit on the same layer, or insert failed. 
	   
  ChamberHitContainer::const_iterator it; 
	   
  for(it = proto_segment.begin(); it != proto_segment.end(); it++) 
    if (((*it)->cscDetId().layer() == layer) && (aHit != (*it))) 
      return false;  
	   
  proto_segment.push_back(aHit); 
  updateParameters(); 
	   
  return true; 
} 
 
void CSCSegAlgoRU::updateParameters() { 
	   
  // Note that we need local position of a RecHit w.r.t. the CHAMBER 
  // and the RecHit itself only knows its local position w.r.t. 
  // the LAYER, so need to explicitly transform to global position. 
	   
  //  no. of hits in the RecHitsOnSegment 
  //  By construction this is the no. of layers with hits 
  //  since we allow just one hit per layer in a segment. 
	   
  int nh = proto_segment.size(); 
	   
  // First hit added to a segment must always fail here 
  if (nh < 2)  
    return; 
	   
  if (nh == 2) { 
		 
    // Once we have two hits we can calculate a straight line  
    // (or rather, a straight line for each projection in xz and yz.) 
    ChamberHitContainer::const_iterator ih = proto_segment.begin(); 
    int il1 = (*ih)->cscDetId().layer(); 
    const CSCRecHit2D& h1 = (**ih); 
    ++ih;     
    int il2 = (*ih)->cscDetId().layer(); 
    const CSCRecHit2D& h2 = (**ih); 
		 
    //@@ Skip if on same layer, but should be impossible 
    if (il1 == il2)  
      return; 
		 
    const CSCLayer* layer1 = theChamber->layer(il1); 
    const CSCLayer* layer2 = theChamber->layer(il2); 
		 
    GlobalPoint h1glopos = layer1->toGlobal(h1.localPosition()); 
    GlobalPoint h2glopos = layer2->toGlobal(h2.localPosition()); 
		 
    // localPosition is position of hit wrt layer (so local z = 0) 
    theOrigin = h1.localPosition(); 
		 
    // We want hit wrt chamber (and local z will be != 0) 
    LocalPoint h1pos = theChamber->toLocal(h1glopos);   
    LocalPoint h2pos = theChamber->toLocal(h2glopos);   
		 
    float dz = h2pos.z()-h1pos.z(); 
    uz = (h2pos.x()-h1pos.x())/dz ; 
    vz = (h2pos.y()-h1pos.y())/dz ; 
		 
    theChi2 = 0.; 
  } 
  else if (nh > 2) { 
		 
    // When we have more than two hits then we can fit projections to straight lines 
    fitSlopes();   
    fillChiSquared(); 
  } // end of 'if' testing no. of hits 
	   
  fillLocalDirection();  
} 
 
//function used for wire group and strip width selection  
int CSCSegAlgoRU::stationNumber(const CSCRecHit2D* currHit) const{ 
  CSCDetId id = (CSCDetId)(*currHit).cscDetId(); 
  int thering = id.ring(); 
  int thestation = id.station(); 
  int iStn = -1; 
  if( thestation == 1 && thering == 4) iStn = 0; 
  if( thestation == 1 && thering == 1) iStn = 1; 
  if( thestation == 1 && thering == 2) iStn = 2; 
  if( thestation == 1 && thering == 3) iStn = 3; 
  if( thestation == 2 && thering == 1) iStn = 4; 
  if( thestation == 2 && thering == 2) iStn = 5; 
  if( thestation == 3 && thering == 1) iStn = 6; 
  if( thestation == 3 && thering == 2) iStn = 7; 
  if( thestation == 4 && thering == 1) iStn = 8; 
  if( thestation == 4 && thering == 2) iStn = 9; 
  return iStn; 
} 
 
 
float CSCSegAlgoRU::fit_r_phi(CLHEP::HepMatrix points, int layer) const{ 
  //find R or Phi on the given layer using the given points for the interpolation 
  float Sx  = 0; 
  float Sy  = 0; 
  float Sxx = 0; 
  float Sxy = 0; 
  for (int i=1;i<7;i++){ 
    if (points(i,1) < -10) continue; 
    Sy = Sy + (points(i,1)); 
    Sx = Sx + i; 
    Sxx = Sxx + (i*i); 
    Sxy = Sxy + ((i)*points(i,1)); 
  } 
  float delta = 2*Sxx - Sx*Sx; 
  float intercept = (Sxx*Sy - Sx*Sxy)/delta; 
  float slope = (2*Sxy - Sx*Sy)/delta;                                                                                                                        
  return (intercept + slope*layer); 
 
} 
 
void CSCSegAlgoRU::baseline(int n_seg_min){ 
	 
  int nhits      = proto_segment.size(); 
  ChamberHitContainer::const_iterator iRH_worst; 
		   
  CLHEP::HepMatrix sp(6,1); 
		 
  sp(1,1) = -999; 
  sp(2,1) = -999; 
  sp(3,1) = -999; 
  sp(4,1) = -999; 
  sp(5,1) = -999; 
  sp(6,1) = -999; 
 
  CLHEP::HepMatrix se(6,1);  
  ChamberHitContainer buffer1; 
  buffer1.clear(); 
 
  unsigned int init_size = proto_segment.size(); 
  while (buffer1.size()< init_size){ 
    ChamberHitContainer::iterator min; 
    int min_layer = 10; 
    for(ChamberHitContainer::iterator k = proto_segment.begin(); k != proto_segment.end(); k++){       
      const CSCRecHit2D* iRHk = *k;  
      CSCDetId idRHk = (CSCDetId)(*iRHk).cscDetId(); 
      int kLayer   = idRHk.layer(); 
      if(kLayer < min_layer){ 
	min_layer = kLayer; 
	min = k; 
      } 
    } 
    buffer1.push_back(*min); 
    proto_segment.erase(min); 
  }//while 
  proto_segment.clear();  
		 
  for (ChamberHitContainer::const_iterator cand = buffer1.begin(); cand != buffer1.end(); cand++) { 
    proto_segment.push_back(*cand); 
  } 
	 
  for(ChamberHitContainer::const_iterator iRH = proto_segment.begin(); iRH != proto_segment.end(); iRH++){       
    hitPosInStrips(*iRH);		 
  } 
	 
 
  float chi2_str; 
  fitX(sp, se, chi2_str); 
		  
  //----------------------------------------------------- 
  // Optimal point rejection method 
  //----------------------------------------------------- 
 
 
  float minSum = 1000; 
  int i1b = 0; 
  int i2b = 0;  
  int iworst = -1;  
  int bad_layer = -1; 
  ChamberHitContainer::iterator rh_to_be_deleted_1; 
  ChamberHitContainer::iterator rh_to_be_deleted_2; 
  if ( nhits > n_seg_min && (chi2_str/(nhits-2)) > chi2_str_*chi2D_iadd){ 
    for (ChamberHitContainer::iterator i1 = proto_segment.begin(); i1 != proto_segment.end();++i1) { 
      ++i1b; 
      const CSCRecHit2D* i1_1 = *i1;  
      CSCDetId idRH1 = (CSCDetId)(*i1_1).cscDetId(); 
      int z1 = idRH1.layer(); 
      i2b = i1b; 
      for (ChamberHitContainer::iterator i2 = i1+1; i2 != proto_segment.end(); ++i2) { 
	++i2b;  
	const CSCRecHit2D* i2_1 = *i2; 
	CSCDetId idRH2 = (CSCDetId)(*i2_1).cscDetId(); 
	int z2 = idRH2.layer(); 
	int irej = 0; 
			 
	for ( ChamberHitContainer::iterator ir = proto_segment.begin(); ir != proto_segment.end(); ++ir) { 
	  ++irej;  
			 
	  if (ir == i1 || ir == i2) continue;  
	  float dsum = 0; 
	  int hit_nr = 0; 
	  const CSCRecHit2D* ir_1 = *ir; 
	  CSCDetId idRH = (CSCDetId)(*ir_1).cscDetId(); 
	  int worst_layer = idRH.layer(); 
	  for (ChamberHitContainer::iterator i = proto_segment.begin(); i != proto_segment.end(); ++i) {  
	    ++hit_nr;   
	    const CSCRecHit2D* i_1 = *i; 
	    if (i == i1 || i == i2 || i == ir) continue;  
	    float slope = (sp(z2,1)-sp(z1,1))/(z2-z1); 
	    float intersept = sp(z1,1) - slope*z1; 
	    CSCDetId idRH = (CSCDetId)(*i_1).cscDetId(); 
	    int z = idRH.layer(); 
	    float di = fabs(sp(z,1) - intersept - slope*z); 
	    dsum = dsum + di; 
					  
	  }//i 
	  if (dsum < minSum){ 
	    minSum = dsum;  
	    bad_layer = worst_layer;		   
	    iworst = irej;  
	    rh_to_be_deleted_1 = ir; 
	  } 
	  
			 
	}//ir  
		    
			 
      }//i2 
    }//i1 
    fitX_ir(sp, se, bad_layer, chi2_str); 
		 
  }//if chi2prob<1.0e-4 
 
					 
  //find worst from n-1 hits 
 
  int iworst2 = -1; 
  int bad_layer2 = -1; 
  if (iworst > -1 && (nhits-1) > n_seg_min && (chi2_str/(nhits-3)) > chi2_str_*chi2D_iadd){ 
    iworst = -1; 
    float minSum = 1000; 
    int i1b = 0; 
    int i2b = 0;  
    for (ChamberHitContainer::iterator i1 = proto_segment.begin(); i1 != proto_segment.end();++i1) { 
      ++i1b;  
      const CSCRecHit2D* i1_1 = *i1; 
      CSCDetId idRH1 = (CSCDetId)(*i1_1).cscDetId(); 
      int z1 = idRH1.layer(); 
      i2b = i1b; 
      for ( ChamberHitContainer::iterator i2 = i1+1; i2 != proto_segment.end(); ++i2) { 
	++i2b;  
	const CSCRecHit2D* i2_1 = *i2; 

	CSCDetId idRH2 = (CSCDetId)(*i2_1).cscDetId(); 
	int z2 = idRH2.layer(); 
	int irej = 0; 
			 
	for ( ChamberHitContainer::iterator ir = proto_segment.begin(); ir != proto_segment.end(); ++ir) { 
			  
	  ++irej;   
	  int irej2 = 0; 
	  if (ir == i1 || ir == i2 ) continue;  
	  const CSCRecHit2D* ir_1 = *ir; 
	  CSCDetId idRH = (CSCDetId)(*ir_1).cscDetId(); 
	  int worst_layer = idRH.layer(); 
	  for (  ChamberHitContainer::iterator ir2 = proto_segment.begin(); ir2 != proto_segment.end(); ++ir2) { 
			  
	    ++irej2;   
	    if (ir2 == i1 || ir2 == i2 || ir2 ==ir ) continue;  
	    float dsum = 0; 
	    int hit_nr = 0; 
	    const CSCRecHit2D* ir2_1 = *ir2; 
	    CSCDetId idRH = (CSCDetId)(*ir2_1).cscDetId(); 
	    int worst_layer2 = idRH.layer(); 
	    for (  ChamberHitContainer::iterator i = proto_segment.begin(); i != proto_segment.end(); ++i) {  
	      ++hit_nr;  
	      const CSCRecHit2D* i_1 = *i; 
	      if (i == i1 || i == i2 || i == ir|| i == ir2 ) continue;  
	      float slope = (sp(z2,1)-sp(z1,1))/(z2-z1); 
	      float intersept = sp(z1,1) - slope*z1; 
	      CSCDetId idRH = (CSCDetId)(*i_1).cscDetId(); 
	      int z = idRH.layer(); 
	      float di = fabs(sp(z,1) - intersept - slope*z); 
	      dsum = dsum + di; 
					  
	    }//i 
	    if (dsum < minSum){ 
	      minSum = dsum;  
	      iworst2 = irej2;  
	      iworst = irej; 
	      bad_layer = worst_layer; 
	      bad_layer2 = worst_layer2;  
	      rh_to_be_deleted_1 = ir; 
	      rh_to_be_deleted_2 = ir2; 
	    } 
	  }//ir2  
	}//ir  
      }//i2 
    }//i1	 
 
    fitX_ir2(sp, se, bad_layer ,bad_layer2, chi2_str); 
  }//if prob(n-1)<e-4 
			 
  //---------------------------------- 
  //erase bad_hits 
  //---------------------------------- 
		 
 
  if( iworst2-1 >= 0 && iworst2 <= int(proto_segment.size())  ) { 
    proto_segment.erase( rh_to_be_deleted_2); 
  } 
		   
  if( iworst-1 >= 0 && iworst <= int(proto_segment.size())  ){ 
    proto_segment.erase(rh_to_be_deleted_1); 
  } 
} 
 
float CSCSegAlgoRU::fit_sp(CLHEP::HepMatrix points, CLHEP::HepMatrix errors, int layer){ 
 
  float S   = 0; 
  float Sx  = 0; 
  float Sy  = 0; 
  float Sxx = 0; 
  float Sxy = 0; 
  float sigma2 = 0; 
  for (int i=1;i<7;i++){ 
    if (points(i,1) < 0) continue; 
    sigma2 = errors(i,1)*errors(i,1); 
    S = S + (1/sigma2); 
    Sy = Sy + (points(i,1)/sigma2); 
    Sx = Sx + ((i)/sigma2); 
    Sxx = Sxx + (i*i)/sigma2; 
    Sxy = Sxy + (((i)*points(i,1))/sigma2); 
  } 
  float delta = S*Sxx - Sx*Sx; 
  float intercept = (Sxx*Sy - Sx*Sxy)/delta; 
  float slope = (S*Sxy - Sx*Sy)/delta; 
  return (intercept + slope*layer); 
 
} 
 
 
float  CSCSegAlgoRU::fitX(CLHEP::HepMatrix points, CLHEP::HepMatrix errors, float &chi2_str){ 
 
  float S   = 0; 
  float Sx  = 0; 
  float Sy  = 0; 
  float Sxx = 0; 
  float Sxy = 0; 
  float sigma2 = 0; 
 
  for (int i=1;i<7;i++){ 
    if (points(i,1) < 0) continue; 
    sigma2 = errors(i,1)*errors(i,1); 
    float i1 = i - 3.5; 
    S = S + (1/sigma2); 
    Sy = Sy + (points(i,1)/sigma2); 
    Sx = Sx + ((i1)/sigma2); 
    Sxx = Sxx + (i1*i1)/sigma2; 
    Sxy = Sxy + (((i1)*points(i,1))/sigma2); 
  } 
	   
  float delta = S*Sxx - Sx*Sx; 
  float intercept = (Sxx*Sy - Sx*Sxy)/delta; 
  float slope = (S*Sxy - Sx*Sy)/delta; 
 
  float chi_str = 0; 
  chi2_str = 0; 
 
  // calculate chi2_str 
  for (int i=1;i<7;i++){ 
    if (points(i,1) < 0) continue; 
    chi_str = (points(i,1) - intercept - slope*(i-3.5))/(errors(i,1)); 
    chi2_str = chi2_str + chi_str*chi_str;  
  } 
  return (intercept + slope*0); 
 
} 
 
float CSCSegAlgoRU::fitX_ir(CLHEP::HepMatrix points, CLHEP::HepMatrix errors, int ir, float &chi2_str){ 
 
  float S   = 0; 
  float Sx  = 0; 
  float Sy  = 0; 
  float Sxx = 0; 
  float Sxy = 0; 
  float sigma2 = 0; 
 
  for (int i=1;i<7;i++){ 
    if (i == ir || points(i,1) < 0) continue; 
    sigma2 = errors(i,1)*errors(i,1); 
    float i1 = i - 3.5; 
    S = S + (1/sigma2); 
    Sy = Sy + (points(i,1)/sigma2); 
    Sx = Sx + ((i1)/sigma2); 
    Sxx = Sxx + (i1*i1)/sigma2; 
    Sxy = Sxy + (((i1)*points(i,1))/sigma2); 
  } 
	   
  float delta = S*Sxx - Sx*Sx; 
  float intercept = (Sxx*Sy - Sx*Sxy)/delta; 
  float slope = (S*Sxy - Sx*Sy)/delta; 
 
  float chi_str = 0; 
  chi2_str = 0; 
 
  // calculate chi2_str 
  for (int i=1;i<7;i++){ 
    if (i == ir || points(i,1) < 0) continue; 
    chi_str = (points(i,1) - intercept - slope*(i-3.5))/(errors(i,1)); 
    chi2_str = chi2_str + chi_str*chi_str;  
  } 
 
  return (intercept + slope*0); 
 
} 
float CSCSegAlgoRU::fitX_ir2(CLHEP::HepMatrix points, CLHEP::HepMatrix errors, int ir, int ir2, float &chi2_str){ 
 
 
  float S   = 0; 
  float Sx  = 0; 
  float Sy  = 0; 
  float Sxx = 0; 
  float Sxy = 0; 
  float sigma2 = 0; 
 
  for (int i=1;i<7;i++){ 
    if (i == ir || i == ir2 || points(i,1) < 0) continue; 
    sigma2 = errors(i,1)*errors(i,1); 
    float i1 = i - 3.5; 
    S = S + (1/sigma2); 
    Sy = Sy + (points(i,1)/sigma2); 
    Sx = Sx + ((i1)/sigma2); 
    Sxx = Sxx + (i1*i1)/sigma2; 
    Sxy = Sxy + (((i1)*points(i,1))/sigma2); 
  } 
	   
  float delta = S*Sxx - Sx*Sx; 
  float intercept = (Sxx*Sy - Sx*Sxy)/delta; 
  float slope = (S*Sxy - Sx*Sy)/delta; 
 
  float chi_str = 0; 
  chi2_str = 0; 
 
 
  // calculate chi2_str 
  for (int i=1;i<7;i++){ 
    if (i == ir || i == ir2 || points(i,1) < 0) continue; 
    chi_str = (points(i,1) - intercept - slope*(i-3.5))/(errors(i,1)); 
    chi2_str = chi2_str + chi_str*chi_str;  
  } 
 
  return (intercept + slope*0); 
 
} 
 
void CSCSegAlgoRU::fitSlopes() { 
	   
  // Update parameters of fit 
  // ptc 13-Aug-02: This does a linear least-squares fit 
  // to the hits associated with the segment, in the z projection. 
	  
  // In principle perhaps one could fit the strip and wire 
  // measurements (u, v respectively), to 
  // u = u0 + uz * z 
  // v = v0 + vz * z 
  // where u0, uz, v0, vz are the parameters resulting from the fit. 
  // But what is actually done is fit to the local x, y coordinates  
  // of the RecHits. However the strip measurement controls the precision 
  // of x, and the wire measurement controls that of y. 
  // Precision in local coordinate: 
  //       u (strip, sigma~200um), v (wire, sigma~1cm) 
	   
  // I have verified that this code agrees with the formulation given 
  // on p246-247 of 'Data analysis techniques for high-energy physics 
  // experiments' by Bock, Grote, Notz & Regler, and that on p111-113 
  // of 'Statistics' by Barlow. 
	   
  // Formulate the matrix equation representing the least-squares fit 
  // We have a vector of measurements m, which is a 2n x 1 dim matrix 
  // The transpose mT is (u1, v1, u2, v2, ..., un, vn) 
  // where ui is the strip-associated measurement and vi is the 
  // wire-associated measurement for a given RecHit i. 
  // The fit is to 
  // u = u0 + uz * z 
  // v = v0 + vz * z 
  // where u0, uz, v0, vz are the parameters to be obtained from the fit. 
  // These are contained in a vector p which is a 4x1 dim matrix, and 
  // its transpose pT is (u0, v0, uz, vz). Note the ordering! 
  // The covariance matrix for each pair of measurements is 2 x 2 and 
  // the inverse of this is the error matrix E. 
  // The error matrix for the whole set of n measurements is a diagonal 
  // matrix with diagonal elements the individual 2 x 2 error matrices 
  // (because the inverse of a diagonal matrix is a diagonal matrix 
  // with each element the inverse of the original.) 
	   
  // The matrix 'matrix' in method 'CSCSegment::weightMatrix()' is this  
  // block-diagonal overall covariance matrix. It is inverted to the  
  // block-diagonal error matrix right before it is returned. 
	   
  // Use the matrix A defined by 
  //    1   0   z1  0 
  //    0   1   0   z1 
  //    1   0   z2  0 
  //    0   1   0   z2 
  //    ..  ..  ..  .. 
  //    1   0   zn  0 
  //    0   1   0   zn 
	   
  // The matrix A is returned by 'CSCSegment::derivativeMatrix()'. 
	   
  // Then the normal equations are encapsulated in the matrix equation 
  // 
  //    (AT E A)p = (AT E)m 
  // 
  // where AT is the transpose of A. 
  // We'll call the combined matrix on the LHS, M, and that on the RHS, B: 
  //     M p = B 
	   
  // We solve this for the parameter vector, p. 
  // The elements of M and B then involve sums over the hits 
	   
  // The error matrix of the parameters is obtained by  
  // (AT E A)^-1 calculated in 'calculateError()'. 
	   
  // The 4 values in p can be accessed from 'CSCSegment::parameters()' 
  // in the order uz, vz, u0, v0. 
	   
  // NOTE 1 
  // Do the #hits = 2 case separately. 
  // (I hope they're not on the same layer! They should not be, by construction.) 
	   
  // NOTE 2 
  // We need local position of a RecHit w.r.t. the CHAMBER 
  // and the RecHit itself only knows its local position w.r.t. 
  // the LAYER, so we must explicitly transform global position. 
	   
  CLHEP::HepMatrix M(4,4,0); 
  CLHEP::HepVector B(4,0); 
	   
  ChamberHitContainer::const_iterator ih = proto_segment.begin(); 
	   
  for (ih = proto_segment.begin(); ih != proto_segment.end(); ++ih) { 
		 
    const CSCRecHit2D& hit = (**ih); 
    const CSCLayer* layer = theChamber->layer(hit.cscDetId().layer()); 
    GlobalPoint gp = layer->toGlobal(hit.localPosition()); 
    LocalPoint  lp  = theChamber->toLocal(gp);  
		 
    // ptc: Local position of hit w.r.t. chamber 
    double u = lp.x(); 
    double v = lp.y(); 
    double z = lp.z(); 
		 
    // ptc: Covariance matrix of local errors  
    CLHEP::HepMatrix IC(2,2); 
    IC(1,1) = hit.localPositionError().xx(); 
    IC(1,2) = hit.localPositionError().xy(); 
    IC(2,1) = IC(1,2); // since Cov is symmetric 
    IC(2,2) = hit.localPositionError().yy(); 
		 
    // ptc: Invert covariance matrix (and trap if it fails!) 

    int ierr = 0; 
    IC.invert(ierr); // inverts in place 
   		 
    // ptc: Note that IC is no longer Cov but Cov^-1 
    M(1,1) += IC(1,1); 
    M(1,2) += IC(1,2); 
    M(1,3) += IC(1,1) * z; 
    M(1,4) += IC(1,2) * z; 
    B(1) += u * IC(1,1) + v * IC(1,2); 
		 
    M(2,1) += IC(2,1); 
    M(2,2) += IC(2,2); 
    M(2,3) += IC(2,1) * z; 
    M(2,4) += IC(2,2) * z; 
    B(2) += u * IC(2,1) + v * IC(2,2); 
		 
    M(3,1) += IC(1,1) * z; 
    M(3,2) += IC(1,2) * z; 
    M(3,3) += IC(1,1) * z * z; 
    M(3,4) += IC(1,2) * z * z; 
    B(3) += ( u * IC(1,1) + v * IC(1,2) ) * z; 
		 
    M(4,1) += IC(2,1) * z; 
    M(4,2) += IC(2,2) * z; 
    M(4,3) += IC(2,1) * z * z; 
    M(4,4) += IC(2,2) * z * z; 
    B(4) += ( u * IC(2,1) + v * IC(2,2) ) * z; 
  } 
	   
  // Solve the matrix equation using CLHEP's 'solve' 
  //@@ ptc: CAN solve FAIL?? UNCLEAR FROM (LACK OF) CLHEP DOC 
  CLHEP::HepVector p = solve(M, B); 
	   
  // Update member variables uz, vz, theOrigin 
  // Note it has local z = 0 
  theOrigin = LocalPoint(p(1), p(2), 0.); 
  uz = p(3); 
  vz = p(4); 
} 
 
void CSCSegAlgoRU::fillChiSquared() { 
	   
  // The chi-squared is (m-Ap)T E (m-Ap) 
  // where T denotes transpose. 
  // This collapses to a simple sum over contributions from each 
  // pair of measurements. 
  float u0 = theOrigin.x(); 
  float v0 = theOrigin.y(); 
  double chsq = 0.; 
	   
  ChamberHitContainer::const_iterator ih; 
  for (ih = proto_segment.begin(); ih != proto_segment.end(); ++ih) { 
		 
    const CSCRecHit2D& hit = (**ih); 
    const CSCLayer* layer = theChamber->layer(hit.cscDetId().layer()); 
    GlobalPoint gp = layer->toGlobal(hit.localPosition()); 
    LocalPoint lp = theChamber->toLocal(gp);  // FIX !! 
		 
    double hu = lp.x(); 
    double hv = lp.y(); 
    double hz = lp.z(); 
		 
    double du = u0 + uz * hz - hu; 
    double dv = v0 + vz * hz - hv; 
		 
    CLHEP::HepMatrix IC(2,2); 
    IC(1,1) = hit.localPositionError().xx(); 
    IC(1,2) = hit.localPositionError().xy(); 
    IC(2,1) = IC(1,2); 
    IC(2,2) = hit.localPositionError().yy(); 
		 
    // Invert covariance matrix 
    int ierr = 0; 
    IC.invert(ierr); 
   		 
    chsq += du*du*IC(1,1) + 2.*du*dv*IC(1,2) + dv*dv*IC(2,2); 
  } 
  theChi2 = chsq; 
} 
 
void CSCSegAlgoRU::fillLocalDirection() { 
  // Always enforce direction of segment to point from IP outwards 
  // (Incorrect for particles not coming from IP, of course.) 
	   
  double dxdz = uz; 
  double dydz = vz; 
  double dz = 1./sqrt(1. + dxdz*dxdz + dydz*dydz); 
  double dx = dz*dxdz; 
  double dy = dz*dydz; 
  LocalVector localDir(dx,dy,dz); 
	   
  // localDir may need sign flip to ensure it points outward from IP 
  // ptc: Examine its direction and origin in global z: to point outward 
  // the localDir should always have same sign as global z... 
	   
  double globalZpos = ( theChamber->toGlobal( theOrigin ) ).z(); 
  double globalZdir = ( theChamber->toGlobal( localDir ) ).z(); 
  double directionSign = globalZpos * globalZdir; 
	   
  theDirection = (directionSign * localDir).unit(); 
	   
} 
 
bool CSCSegAlgoRU::hasHitOnLayer(int layer) const { 
	   
  // Is there is already a hit on this layer? 
  ChamberHitContainerCIt it; 
	   
  for(it = proto_segment.begin(); it != proto_segment.end(); it++) 
    if ((*it)->cscDetId().layer() == layer) 
      return true;  
	   
  return false; 
} 
 
bool CSCSegAlgoRU::replaceHit(const CSCRecHit2D* h, int layer) { 
	   
  // replace a hit from a layer  
  ChamberHitContainer::iterator it; 
  for (it = proto_segment.begin(); it != proto_segment.end();) { 
    if ((*it)->cscDetId().layer() == layer) 
      it = proto_segment.erase(it); 
    else 
      ++it;    
  } 
	   
  return addHit(h, layer);				     
} 
 
void CSCSegAlgoRU::compareProtoSegment(const CSCRecHit2D* h, int layer) { 
	   
  // compare the chi2 of two segments 
  double oldChi2 = theChi2; 
  LocalPoint oldOrigin = theOrigin; 
  LocalVector oldDirection = theDirection; 
  ChamberHitContainer oldSegment = proto_segment; 
	   
  bool ok = replaceHit(h, layer); 
	   
  if ((theChi2 > oldChi2) || (!ok)) { 
    proto_segment = oldSegment; 
    theChi2 = oldChi2; 
    theOrigin = oldOrigin; 
    theDirection = oldDirection; 
  } 
} 
 
void CSCSegAlgoRU::increaseProtoSegment(const CSCRecHit2D* h, int layer, int chi2_factor) { 
	   
  double oldChi2 = theChi2; 
  LocalPoint oldOrigin = theOrigin; 
  LocalVector oldDirection = theDirection; 
  ChamberHitContainer oldSegment = proto_segment; 
	   
  bool ok = addHit(h, layer); 
	   
  int ndf = 2*proto_segment.size() - 4; 
  
  if ( !ok || ((ndf > 0) && (theChi2/ndf > chi2Max*chi2_factor)) ) { 
    proto_segment = oldSegment; 
    theChi2 = oldChi2; 
    theOrigin = oldOrigin; 
    theDirection = oldDirection; 
  }			 
}		 
 
CLHEP::HepMatrix CSCSegAlgoRU::derivativeMatrix() const { 
	   
  ChamberHitContainer::const_iterator it; 
  int nhits = proto_segment.size(); 
  CLHEP::HepMatrix matrix(2*nhits, 4); 
  int row = 0; 
	   
  for(it = proto_segment.begin(); it != proto_segment.end(); ++it) { 
		 
    const CSCRecHit2D& hit = (**it); 
    const CSCLayer* layer = theChamber->layer(hit.cscDetId().layer()); 
    GlobalPoint gp = layer->toGlobal(hit.localPosition());    	 
    LocalPoint lp = theChamber->toLocal(gp);  
    float z = lp.z(); 
    ++row; 
    matrix(row, 1) = 1.; 
    matrix(row, 3) = z; 
    ++row; 
    matrix(row, 2) = 1.; 
    matrix(row, 4) = z; 
  } 
  return matrix; 
} 
 
 
AlgebraicSymMatrix CSCSegAlgoRU::weightMatrix() const { 
	   
  std::vector<const CSCRecHit2D*>::const_iterator it; 
  int nhits = proto_segment.size(); 
  AlgebraicSymMatrix matrix(2*nhits, 0); 
  int row = 0; 
	   
  for (it = proto_segment.begin(); it != proto_segment.end(); ++it) { 
		 
    const CSCRecHit2D& hit = (**it); 
    ++row; 
    matrix(row, row)   = hit.localPositionError().xx(); 
    matrix(row, row+1) = hit.localPositionError().xy(); 
    ++row; 
    matrix(row, row-1) = hit.localPositionError().xy(); 
    matrix(row, row)   = hit.localPositionError().yy(); 
  } 
  int ierr; 
  matrix.invert(ierr); 
  return matrix; 
} 
 
AlgebraicSymMatrix CSCSegAlgoRU::calculateError() const { 
	   
  AlgebraicSymMatrix weights = weightMatrix(); 
  AlgebraicMatrix A = derivativeMatrix(); 
	   
  // (AT W A)^-1 
  // from http://www.phys.ufl.edu/~avery/fitting.html, part I 
  int ierr; 
  AlgebraicSymMatrix result = weights.similarityT(A); 
  result.invert(ierr); 
	   
  // blithely assuming the inverting never fails... 
  return result; 
} 
 
void CSCSegAlgoRU::flipErrors( AlgebraicSymMatrix& a ) const { 
 
  // The CSCSegment needs the error matrix re-arranged 
 
  AlgebraicSymMatrix hold( a ); 
 
  // errors on slopes into upper left 
  a(1,1) = hold(3,3); 
  a(1,2) = hold(3,4); 
  a(2,1) = hold(4,3); 
  a(2,2) = hold(4,4); 
	   
  // errors on positions into lower right 
  a(3,3) = hold(1,1); 
  a(3,4) = hold(1,2); 
  a(4,3) = hold(2,1); 
  a(4,4) = hold(2,2); 
 
  // off-diagonal elements remain unchanged 
 
} 
