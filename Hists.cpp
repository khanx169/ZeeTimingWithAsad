struct HistSet{
//book histogram set w/ common suffix inside the provided TFileDirectory
//void book(edm::Service<TFileService>& td,const std::string&);
void book(TFileDirectory subDir,const std::string&);
// fill all histos of the set with the two electron candidates
void fill(int sc1, int sc2, int cl1, int cl2);


TH1F* eta1_;
TH1F* phi1_;
TH1F* eta2_;
TH1F* phi2_;
TH1F* e1pt_;
TH1F* e2pt_;
TH1F* e1eta_;
TH1F* e2eta_;
TH1F* e1phi_;
TH1F* e2phi_;
TH1F* e1E_;
TH1F* e2E_;

} theHists;




void HistSet::book(TFileDirectory subDir, const std::string& post) {

	eta1_ =(TH1F*) subDir.make<TH1F>("eta","eta",100,eta_max,eta_min);
	phi1_ =(TH1F*) subDir.make<TH1F>("phi","phi",100,phi_max,phi_min);
	.
	.
	.
}


void HistSet::fill(int sc1, int sc2, int bc1, int bc2 ){

	float et1 = treeVars_.superClusterRawEnergy[sc1]/cosh( treeVars_.superClusterEta[sc1] );
	float phi1 = treeVars_.superClusterPhi[sc1];

	math::PtEtaPhiELorentzVectorD el1(et1 ,
						treeVars_.superClusterEta[sc1],
						treeVars_.superClusterPhi[sc1],
						treeVars_.superClusterRawEnergy[sc1] );


	float et2 = treeVars_.superClusterRawEnergy[sc2]/cosh( treeVars_.superClusterEta[sc2] );
	float phi2 = treeVars_.superClusterPhi[sc1];

	math::PtEtaPhiELorentzVectorD el2(et2 ,
						treeVars_.superClusterEta[sc2],
						treeVars_.superClusterPhi[sc2],
						treeVars_.superClusterRawEnergy[sc2] );

								
	math::PtEtaPhiELorentzVectorD diEle = el1;
	diEle += el2;

	

	eta1_		->Fill(et1);
	phi1_		->Fill(phi1);

	e1E_		->Fill(treeVars_.superClusterRawEnergy[sc1])
	.
	.
	.

}s
	

	








	
