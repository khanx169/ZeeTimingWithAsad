void func1( int sc1, int sc2, int eta1_min, int phi1_min, int eta1_max, int phi1_max, int eta2_min, int phi2_min, int eta2_max, int phi2_max )
{
	int eventCounter = 0;
	for (int entry = 0 ; (entry < nEntries && eventCounter < numEvents_); ++entry)
	{
		chain->GetEntry (entry) ;
		// Keep the event?
		bool keepEvent = includeEvent(treeVars_.l1ActiveTriggers,
			treeVars_.l1NActiveTriggers,trigIncludeVector,trigExcludeVector)
				&& includeEvent(treeVars_.l1ActiveTechTriggers,
					treeVars_.l1NActiveTechTriggers,ttrigIncludeVector,ttrigExcludeVector);
		if(!keepEvent)
			continue;
		 // do analysis if the run is in the desired range
		if( treeVars_.runId<minRun_ || maxRun_<treeVars_.runId) continue;
		// do analysis if the LS is in the desired range
		if( treeVars_.lumiSection<minLS_ || maxLS_<treeVars_.lumiSection) continue;
		bool verticesAreOnlyNextToNominalIP;
		int count=0;
		for(int v=0; v<treeVars_.nVertices; v++ )
		{ if (fabs(treeVars_.vtxZ[0])<15) count++; }
		if ( treeVars_.nVertices >0 && count==treeVars_.nVertices ) verticesAreOnlyNextToNominalIP = true;
		else verticesAreOnlyNextToNominalIP = false;
		// --vertex: require vertex@IP (1), veto it (2) or either (0, or unset)
		if (flagOneVertex_ ==1 && (!verticesAreOnlyNextToNominalIP) ) continue;
		if (flagOneVertex_ ==2 && (verticesAreOnlyNextToNominalIP) ) continue;
		
		// if evet being actually processed, increment counter of analyzed events
		eventCounter++;
		speak_=false;
		if (entry<10 || entry%10000==0) speak_=true;
		if (speak_) std::cout << "\n\n------> reading entry " << entry << "\tLS: " << treeVars_.lumiSection << " <------\n" ;
		if (speak_) std::cout << " found " << treeVars_.nSuperClusters << " superclusters" << std::endl ;
		if (speak_) std::cout << " found " << treeVars_.nClusters << " basic clusters" << std::endl ;

		
		//New Code 
		
		float et1 = treeVars_.superClusterRawEnergy[sc1]/cosh( treeVars_.superClusterEta[sc1] );
		if (et1 < eta1_min || eta1 > eta1_max) continue;  //.................... CHECK IF ETA1 IS WITHIN THE SPECIFIED RANGE ...........................\\

		float phi1 = treeVars_.superClusterPhi[sc1];
		if (phi1 < phi1_min || phi1 > phi1_max) continue  //.................... CHECK IF PHI1 IS WITHIN THE SPECIFIED RANGE ...........................\\

		math::PtEtaPhiELorentzVectorD el1(et1 ,
						treeVars_.superClusterEta[sc1],
						treeVars_.superClusterPhi[sc1],
						treeVars_.superClusterRawEnergy[sc1] ); 

		
		float et2 = treeVars_.superClusterRawEnergy[sc2]/cosh( treeVars_.superClusterEta[sc2] );
		if (et2 < eta2_min || eta2 > eta2_max) continue;  //.................... CHECK IF ETA2 IS WITHIN THE SPECIFIED RANGE ...........................\\

		float phi2 = treeVars_.superClusterPhi[sc1];
		if (phi2 < phi2_min || phi2 > phi2_max) continue  //.................... CHECK IF PHI2 IS WITHIN THE SPECIFIED RANGE ...........................\\

		math::PtEtaPhiELorentzVectorD el2(et2 ,
						treeVars_.superClusterEta[sc2],
						treeVars_.superClusterPhi[sc2],
						treeVars_.superClusterRawEnergy[sc2] );

		
		float dvertex = pow(treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterVertexZ[sc2],2);
		dvertex = sqrt(dvertex);

		math::PtEtaPhiELorentzVectorD diEle = el1;
		diEle += el2;

		// ////////////////////////
		mass_ ->Fill(diEle.M());
		dZvertices_->Fill(dvertex);
		Zvertices_->Fill( (treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterVertexZ[sc2])/2 );
		nVertices_->Fill(treeVars_.nVertices);
		
		// require invariant mass
		if( fabs( diEle.M() - 91 ) > 40 ) continue;
		// require two electrons from the same vertex
		if ( dvertex > 0.1 ) continue;


		


		// at this stage I have a suitable di-electron system for time studies
		float tmpEne=-9999;
		// loop on BC and match to sc1 ===============
		int bc1=-1;
		for (int bc=0; bc<treeVars_.nClusters; bc++){
		if ( (pow(treeVars_.superClusterEta[sc1]-treeVars_.clusterEta[bc],2)+ pow(treeVars_.superClusterPhi[sc1]-treeVars_.clusterPhi[bc],2) ) < 0.02
		&& treeVars_.clusterEnergy[bc]>tmpEne) {
		tmpEne=treeVars_.clusterEnergy[bc];
		bc1=bc;
		}// end - if good bc candidate
		}// end - loop over BC
		tmpEne=-9999;
		// loop on BC and match to sc2 ==============
		int bc2=-1;
		for (int bc=0; bc<treeVars_.nClusters; bc++){
		if ( pow(treeVars_.superClusterEta[sc2]-treeVars_.clusterEta[bc],2)+ pow(treeVars_.superClusterPhi[sc2]-treeVars_.clusterPhi[bc],2) < 0.02
		&& treeVars_.clusterEnergy[bc]>tmpEne) {
		tmpEne=treeVars_.clusterEnergy[bc];
		bc2=bc;
		}// end - if good bc candidate
		}// end - loop over BC
		// protect in case of no matching
		if(bc1==-1 || bc2==-1) continue;
		if(0) {
		std::cout << "\n\nsc1 : " << treeVars_.superClusterEta[sc1] << " " << treeVars_.superClusterPhi[sc1] << " " << treeVars_.superClusterRawEnergy[sc1] << std::endl;
		std::cout << "bc1 : " << treeVars_.clusterEta[bc1] << " " << treeVars_.clusterPhi[bc1] << " " << treeVars_.clusterEnergy[bc1] << "\n"<< std::endl;
		std::cout << "sc2 : " << treeVars_.superClusterEta[sc2] << " " << treeVars_.superClusterPhi[sc2] << " " << treeVars_.superClusterRawEnergy[sc2] << std::endl;
		std::cout << "bc2 : " << treeVars_.clusterEta[bc2] << " " << treeVars_.clusterPhi[bc2] << " " << treeVars_.clusterEnergy[bc2] << std::endl;
		}
		ClusterTime bcTime1 = timeAndUncertSingleCluster(bc1,treeVars_);
		ClusterTime bcTime2 = timeAndUncertSingleCluster(bc2,treeVars_);
		if(! (bcTime1.isvalid && bcTime2.isvalid) ) continue;
		// fill the structures which hold all the plots
		plotsECALECAL.fill(sc1,sc2, bc1,bc2);
		if ( fabs(treeVars_.clusterEta[bc1])<1.4 && fabs(treeVars_.clusterEta[bc2])<1.4 ){
		plotsEBEB.fill(sc1,sc2, bc1,bc2);
		float energyRatio1 = treeVars_.xtalInBCEnergy[bc1][bcTime1.seed];
		if(bcTime1.second>-1) {energyRatio1 /= treeVars_.xtalInBCEnergy[bc1][bcTime1.second]; }
		else { energyRatio1 /= 99999; }
		float energyRatio2 = treeVars_.xtalInBCEnergy[bc2][bcTime2.seed];
		if(bcTime2.second>-1) {energyRatio2 /= treeVars_.xtalInBCEnergy[bc2][bcTime2.second]; }
		else { energyRatio2 /= 99999; }
		float minRatio = 0.7; float maxRatio = 1.3;
		if(minRatio<energyRatio1 && minRatio<energyRatio2 && energyRatio1<maxRatio && energyRatio2<maxRatio) plotsEBEBequalShare.fill(sc1,sc2, bc1,bc2);
		minRatio = 2; maxRatio = 10;
		if(minRatio<energyRatio1 && minRatio<energyRatio2 && energyRatio1<maxRatio && energyRatio2<maxRatio) plotsEBEBunevenShare.fill(sc1,sc2, bc1,bc2);
		}// if EBEB, and subcases
		else if ( fabs(treeVars_.clusterEta[bc1])>1.5 && fabs(treeVars_.clusterEta[bc2])>1.5 ) plotsEEEE.fill(sc1,sc2, bc1,bc2);
		else if ( (fabs(treeVars_.clusterEta[bc1])<1.4 && fabs(treeVars_.clusterEta[bc2])>1.5) ||
		(fabs(treeVars_.clusterEta[bc1])>1.5 && fabs(treeVars_.clusterEta[bc2])<1.4) ) plotsEBEE.fill(sc1,sc2, bc1,bc2);
		// if I've found a pair of supercluster, bail out of loop to repeat using twice the same supercluster
		break;


		
				
	}
}
