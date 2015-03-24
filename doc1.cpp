//! main program
int main (int argc, char** argv)
{
// First parse arguments
parseArguments(argc, argv);
if (listOfFiles_.size()==0){
std::cout << "\tno input file found" << std::endl;
return(1);
}
else{
std::cout << "\tfound " << listOfFiles_.size() << " input files: " << std::endl;
for(std::vector<std::string>::const_iterator file_itr=listOfFiles_.begin(); file_itr!=listOfFiles_.end(); file_itr++){
std::cout << "\t" << (*file_itr) << std::endl;
}
}



// Tree construction
// FIX should turn this string into a configurable
TChain * chain = new TChain ("EcalTimeAnalysis") ; // ntuple producer in CMSSW CVS
//TChain * chain = new TChain ("EcalTimePi0Analysis") ; // ntuple producer in UserCode/UMN space
std::vector<std::string>::const_iterator file_itr;
for(file_itr=listOfFiles_.begin(); file_itr!=listOfFiles_.end(); file_itr++){
chain->Add( (*file_itr).c_str() );
}
int nEntries = chain->GetEntries () ;
if (numEvents_==-1) numEvents_ = nEntries;
std::cout << "\n\tFOUND " << listOfFiles_.size() << " input files" << std::endl ;
std::cout << "\n\tFOUND " << nEntries << " events" << std::endl ;
std::cout << "\tWILL run on: " << numEvents_ << " events" << std::endl;
std::cout << "\tOutput file: " << outputRootName_ << std::endl;
std::cout << "\tminAOverSigma: " << minAmpliOverSigma_ << std::endl;
std::cout << "\teTGammaMinEB: " << eTGammaMinEB_ << std::endl;
std::cout << "\ts4s9GammaMinEB: " << s4s9GammaMinEB_ << std::endl;
std::cout << "\teTPi0MinEB: " << eTPi0MinEB_ << std::endl;
std::cout << "\teTGammaMinEE: " << eTGammaMinEE_ << std::endl;
std::cout << "\ts4s9GammaMinEE: " << s4s9GammaMinEE_ << std::endl;
std::cout << "\teTPi0MinEE: " << eTPi0MinEE_ << std::endl;
std::cout << "\tminRun: " << minRun_ << std::endl;
std::cout << "\tmaxRun: " << maxRun_ << std::endl;
std::cout << "\tminLS: " << minLS_ << std::endl;
std::cout << "\tmaxLS: " << maxLS_ << std::endl;
setBranchAddresses (chain, treeVars_);
// setting up the TFileService in the ServiceRegistry;
edmplugin::PluginManager::Config config;
edmplugin::PluginManager::configure(edmplugin::standard::config());
std::vector<edm::ParameterSet> psets;
edm::ParameterSet pSet;
pSet.addParameter("@service_type",std::string("TFileService"));
pSet.addParameter("fileName",std::string("TimePerf-plots.root")); // this is the file TFileService will write into
psets.push_back(pSet);
static edm::ServiceToken services(edm::ServiceRegistry::createSet(psets));
static edm::ServiceRegistry::Operate operate(services);
edm::Service<TFileService> fs;
TFileDirectory subDirECALECAL=fs->mkdir("ECALECAL");
HistSet plotsECALECAL;
plotsECALECAL.book(subDirECALECAL,std::string("ECALECAL"));
TFileDirectory subDirEBEB=fs->mkdir("EBEB");
HistSet plotsEBEB;
plotsEBEB.book(subDirEBEB,std::string("EBEB"));
TFileDirectory subDirEEEE=fs->mkdir("EEEE");
HistSet plotsEEEE;
plotsEEEE.book(subDirEEEE,std::string("EEEE"));
TFileDirectory subDirEBEE=fs->mkdir("EBEE");
HistSet plotsEBEE;
plotsEBEE.book(subDirEBEE,std::string("EBEE"));
TFileDirectory subDirEBEBequalShare=fs->mkdir("EBEBequalShare");
HistSet plotsEBEBequalShare;
plotsEBEBequalShare.book(subDirEBEBequalShare,std::string("EBEBequalShare"));
TFileDirectory subDirEBEBunevenShare=fs->mkdir("EBEBunevenShare");
HistSet plotsEBEBunevenShare;
plotsEBEBunevenShare.book(subDirEBEBunevenShare,std::string("EBEBunevenShare"));
timeCorrector theCorr;
std::cout << "\ncreated object theCorr to be used for timeVsAmpliCorrections" << std::endl;
std::cout << "\ninitializing theCorr" << std::endl;
//theCorr.initEB( std::string("EBmod4") );
//theCorr.initEE( std::string("EElow") );
theCorr.initEB( "EB" );
theCorr.initEE( "EE" );
//Initialize output root file
//saving_ = new TFile(outputRootName_.c_str (),"recreate");
// Initialize the histograms
TFileDirectory subDirGeneral=fs->mkdir("General");
initializeHists(subDirGeneral);
int eventCounter = 0;




/////////////////////////////////////////////////////
// Main loop over entries
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
///////////////////////////////////////////////////////////////////////
// outer loop on supercluster
for (int sc1=0; sc1<treeVars_.nSuperClusters; sc1++){
float et1 = treeVars_.superClusterRawEnergy[sc1]/cosh( treeVars_.superClusterEta[sc1] );
if (et1<20) continue;
math::PtEtaPhiELorentzVectorD el1(et1 ,
treeVars_.superClusterEta[sc1],
treeVars_.superClusterPhi[sc1],
treeVars_.superClusterRawEnergy[sc1] );
///////////////////////////////////////////////////////////////////////
// inner loop on supercluster
for (int sc2=(sc1+1); sc2<treeVars_.nSuperClusters; sc2++){
float et2 = treeVars_.superClusterRawEnergy[sc2]/cosh( treeVars_.superClusterEta[sc2] );
if (et2<20) continue;
math::PtEtaPhiELorentzVectorD el2(et2 ,
treeVars_.superClusterEta[sc2],
treeVars_.superClusterPhi[sc2],
treeVars_.superClusterRawEnergy[sc2] );
// there seems to be a problem with vertexing - since nearly none of the electrons have the same vertex... CHECK!
float dvertex = pow(treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterVertexZ[sc2],2);
//dvertex += pow(treeVars_.superClusterVertexY[sc1]-treeVars_.superClusterVertexY[sc2],2);
//dvertex += pow(treeVars_.superClusterVertexX[sc1]-treeVars_.superClusterVertexX[sc2],2);
dvertex = sqrt(dvertex);
math::PtEtaPhiELorentzVectorD diEle = el1;
diEle += el2;
// ////////////////////////
mass_ ->Fill(diEle.M());
dZvertices_->Fill(dvertex);
Zvertices_->Fill( (treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterVertexZ[sc2])/2 );
nVertices_->Fill(treeVars_.nVertices);
// require invariant mass
//if( fabs( diEle.M() - 91 ) > 20 ) continue;
if( fabs( diEle.M() - 91 ) > 40 ) continue;
// require two electrons from the same vertex
//if ( dvertex > 0.01 ) continue;
if ( dvertex > 0.1 ) continue;
if(0) std::cout << "di-electron system mass: " << diEle.M() << " vertex distance: " << dvertex << std::endl;
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
}// end loop sc2
}// end loop sc1
} // end of loop over entries
delete chain ;
return 0 ;
}
