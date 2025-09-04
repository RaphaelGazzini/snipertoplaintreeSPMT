#include "SNiPERToPlainTreeSPMT.h"
#include "TOF.h"

#include "EvtNavigator/NavBuffer.h"
#include "EvtNavigator/EvtNavHelper.h"
#include "SniperKernel/AlgFactory.h"

#include "Event/CdVertexRecHeader.h"
#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdSpmtCalibHeader.h"
#include "Event/CdLpmtElecTruthHeader.h"
#include "Event/CdSpmtElecTruthHeader.h"
#include "Event/CdWaveformHeader.h"
#include "Event/WpCalibHeader.h"
#include "Event/WpCalibEvt.h"

#include "Event/SimHeader.h"

#include "Event/CdLpmtElecEvt.h"
#include "Event/CdSpmtElecEvt.h"
#include "Event/CdLpmtElecHeader.h"
#include "Event/CdSpmtElecHeader.h"

#include "Event/CdTriggerHeader.h"
#include "Event/CdTriggerEvt.h"
#include "Event/WpTriggerHeader.h"
#include "Event/WpTriggerEvt.h"
#include "Event/TtTriggerHeader.h"
#include "Event/MMTriggerHeader.h"

#include "Event/CdVertexRecHeader.h"

#include "Identifier/Identifier.h"
#include "Identifier/IDService.h"
#include "Identifier/CdID.h"

#include "RootWriter/RootWriter.h"
#include "BufferMemMgr/IDataMemMgr.h"

#include "SpmtElecConfigSvc/SpmtElecConfigSvc.h"

#include "OECTagSvc/OECTagSvc.h"
#include "OECTagID/OECTagID.h"

#include "Geometry/IPMTParamSvc.h"

#include "Event/OecHeader.h"
#include "Event/OecEvt.h"

#include "TTree.h"
//#include "Geometry/IRecGeomSvc.hh"
#include <fstream>

DECLARE_ALGORITHM(SNiPERToPlainTreeSPMT);

SNiPERToPlainTreeSPMT::SNiPERToPlainTreeSPMT(const std::string& name)
: AlgBase(name),
  m_iEvt(0),
  m_buf(0),
  m_spmtSvc(0),
  m_tagsvc(0)
{
  declProp("TxtFileName", m_txtFileName);
  declProp("enableElec", saveElec=false);
  declProp("interface", interface = -40000);
  declProp("enableSim", saveSim = false);
  declProp("savePMTInfo", savePMT = false);
}

bool SNiPERToPlainTreeSPMT::initialize()
{
	//----------------------------------------------------------------------------
	const std::string compilation_date = __DATE__;
	const std::string compilation_time = __TIME__;
	std::cout <<"##################################################################"<<std::endl
		<<"The source file was compiled on " << compilation_date<< " at " << compilation_time <<std::endl
		<<"##################################################################"<<std::endl;
	//----------------------------------------------------------------------------

	// =======================================================================
	// Loading PMT positions
	// =======================================================================

	if(savePMT){
		PMTFile = new TFile("PmtInfo.root", "RECREATE");
		PMTTree = new TTree("PMTInfo", "PMT Tree");
		PMTTree->Branch("CopyNo", &m_CopyNo, "CopyNo/I");
		PMTTree->Branch("PMTX", &m_PMTX, "PMTX/D");
		PMTTree->Branch("PMTY", &m_PMTY, "PMTY/D");
		PMTTree->Branch("PMTZ", &m_PMTZ, "PMTZ/D");
		PMTTree->Branch("CircleNo", &m_Circle, "CircleNo/I");
		PMTTree->Branch("Type", &m_PMTType);
		PMTTree->Branch("Pole", &m_Pole);
	}

	IDService* idServ = IDService::getIdServ();
	idServ->init();

	SniperPtr<IPMTParamSvc> pmtsvc(getParent(), "PMTParamSvc");

	if (ALL_LPMT_pos.size()==0 && pmtsvc.valid()) {
        TotalLPMT = pmtsvc->get_NTotal_CD_LPMT();

        std::cout << " PMT Information " << std::endl;

		std::cout << "LPMT" <<std::endl;
        for (unsigned int ith = 0; ith < TotalLPMT; ith++)
        {
            TVector3 all_pmtCenter(pmtsvc->getPMTX(ith), pmtsvc->getPMTY(ith), pmtsvc->getPMTZ(ith));
            ALL_LPMT_pos.push_back(all_pmtCenter);

			if(savePMT){
				m_CopyNo = ith;
				m_PMTX = pmtsvc->getPMTX(ith);
				m_PMTY = pmtsvc->getPMTY(ith);
				m_PMTZ = pmtsvc->getPMTZ(ith);

				Identifier id = idServ->copyNo2Id(ith);
				int type = CdID::pmtType(id);
				m_Circle = CdID::circleNo(id);
				if(type == 1) m_PMTType = "Hamamatsu";
				else if(type == 2) m_PMTType = "NNVT";
				if(CdID::northOrSouth(id) == 0) m_Pole = "North";
				else m_Pole = "South";

				PMTTree->Fill();
			}
        }
	}

	if (ALL_SPMT_pos.size()==0 && pmtsvc.valid()) {
		TotalSPMT = pmtsvc->get_NTotal_CD_SPMT();

		std::cout << "SPMT" <<std::endl;
        for (unsigned int ith = 0; ith < TotalSPMT; ith++)
        {
            TVector3 all_pmtCenter(pmtsvc->getPMTX(ith+20000), pmtsvc->getPMTY(ith+20000), pmtsvc->getPMTZ(ith+20000));
            ALL_SPMT_pos.push_back(all_pmtCenter);
			if(savePMT){
				m_CopyNo = ith+20000;
				m_PMTX = pmtsvc->getPMTX(ith+20000);
				m_PMTY = pmtsvc->getPMTY(ith+20000);
				m_PMTZ = pmtsvc->getPMTZ(ith+20000);
				m_PMTType = "SPMT";

				Identifier id = idServ->copyNo2Id(ith+20000);
				int type = CdID::pmtType(id);
				m_Circle = CdID::circleNo(id);
				if(CdID::northOrSouth(id) == 0) m_Pole = "North";
				else m_Pole = "South";

				PMTTree->Fill();
			}
        }
	}

	if(savePMT){
		PMTFile->cd();
		PMTTree->Write();
		PMTFile->Close();
	}


	// =======================================================================
    // GET EVENT
    // =======================================================================


	LogDebug << "initializing" << std::endl;
	std::cout<<"36"<<std::endl;
	gDirectory->pwd();

	SniperDataPtr<JM::NavBuffer>  navBuf(getParent()/*getRoot()*/,"/Event");
	std::cout<<"Parent="<<getParent()<<"Root="<<getRoot()<<std::endl;//" Parent Name"<<getParentName()<<std::endl;
	if ( navBuf.invalid() ) {
	    LogError << "cannot get the NavBuffer @ /Event" << std::endl;
	    return false;
	}
	m_buf = navBuf.data();
	gDirectory->pwd();
	SniperPtr<SpmtElecConfigSvc> svc(*getRoot(), "SpmtElecConfigSvc");
	if (svc.invalid()) {
	    LogError << "can't find service SpmtElecConfigSvc" << std::endl;
	    return false;
	}
	m_spmtSvc = svc.data();

	//get OECConfig service
	SniperPtr<OECTagSvc> tagsvc(getParent(),"OECTagSvc");
	if( tagsvc.invalid()) {
	    LogError << "Unable to locate tagsvc" << std::endl;
	    return false;
	}
	m_tagsvc = tagsvc.data();
	//get predefined tags (as example BiPo212)
	i_pBiPo212pair= m_tagsvc->getpTag("BiPo212Pair");
	i_dBiPo212pair= m_tagsvc->getdTag("BiPo212Pair");

	i_pBiPo214pair= m_tagsvc->getpTag("BiPo214Pair");
	i_dBiPo214pair= m_tagsvc->getdTag("BiPo214Pair");

	i_pIBD = m_tagsvc->getpTag("InverseBetaDecay");
	i_dIBD = m_tagsvc->getdTag("InverseBetaDecay");

	i_SingleSE = m_tagsvc->getpTag("Single_SE");
	i_SingleME = m_tagsvc->getpTag("Single_ME");
	i_SingleLE = m_tagsvc->getpTag("Single_LE");

	i_pSpallNeutron = m_tagsvc->getpTag("SpallationNeutron_WF");
	i_dSpallNeutron = m_tagsvc->getdTag("SpallationNeutron_WF");

	i_pB12 = m_tagsvc->getpTag("B12");
	i_dB12 = m_tagsvc->getpTag("B12");
		strict_dIBD = 0x01000000 | 0x04000000 | 0x00000200 | 0x00000010;

	Book_tree();
	/*
	std::cout<<"file "<<m_txtFileName<<std::endl;
	std::string txtFile = "fileVtx_RUN.5162.JUNODAQ.LSFilling.ds-2.global_trigger.20250415222246.137.txt"; // ou récupéré dynamiquement

	if (!LoadCalibValues(m_txtFileName)) {
			return false; // ou voir comment gérer l'erreur
		}
	*/
  	return true;
}

bool SNiPERToPlainTreeSPMT::execute()
{
	std::cout << "executing: " << m_iEvt++
		<< std::endl;

	gDirectory->pwd();

//      JM::EvtNavigator* navig = 0;
	JM::SimEvt* simevent = 0;
	JM::CdVertexRecEvt* recevent = 0;
	JM::CdLpmtCalibEvt* calibeventLPMT = 0;
	JM::CdSpmtCalibEvt* calibeventSPMT = 0;
	JM::WpCalibEvt* calibeventWP = 0;
	JM::CdWaveformEvt* elecevent = 0;
	JM::CdLpmtElecTruthEvt *trutheventLPMT = 0;
	JM::CdSpmtElecTruthEvt *trutheventSPMT = 0;
	JM::CdLpmtElecEvt *eventLPMT = 0;
	JM::CdSpmtElecEvt *eventSPMT = 0;
	JM::CdTriggerEvt *triggerevent =0;
	JM::WpTriggerEvt *Wptriggerevent =0;
	JM::OecEvt *oecevt = 0;
	// Get the events of different stages
	// calculate block charge and time

	auto nav = m_buf->curEvt();
	gDirectory->pwd();
	const auto& paths = nav->getPath();
	const auto& refs = nav->getRef();

	LogInfo << "Detector type is  " <<nav->getDetectorType()<<std::endl;
	LogInfo << "Start to Explore SmartRef: " << std::endl;
	LogInfo << "Size of paths: " << paths.size() << std::endl;
	LogInfo << "Size of refs: " << refs.size() << std::endl;

	for (size_t i = 0; i < paths.size(); ++i) {
        LogInfo << refs[i]<<" -> ref: " << std::endl;
        const std::string& path = paths[i];
        JM::SmartRef* ref = refs[i];
        JM::EventObject* evtobj = ref->GetObject();

        LogInfo << " path: " << path
            << " ref->entry(): " << ref->entry()
            << " evtobj: " << evtobj;

        if (path=="/Event/Sim") {
            auto hdr = dynamic_cast<JM::SimHeader*>(evtobj);
            LogInfo <<i<<" SimHeader: " << hdr;
        }
        LogInfo << std::endl;
	}

	auto simheader = JM::getHeaderObject<JM::SimHeader>(nav);
	if(simheader){
		simevent = (JM::SimEvt*)simheader->event();
		LogInfo << "SimEvent Read in: " << simevent << std::endl;
		LogInfo << "SimEvent Track: " << simevent->getTracksVec().size() << std::endl;
		LogInfo << "SimEvent Hits: " << simevent->getCDHitsVec().size() << std::endl;
	}

	auto recheader = JM::getHeaderObject<JM::CdVertexRecHeader>(nav);
	if (recheader) {
	  recevent = recheader->event();
	  LogInfo << "RecEvent Read in: " << recevent << std::endl;
	}

	auto calibheaderLPMT = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
	if (calibheaderLPMT) {
	  calibeventLPMT = calibheaderLPMT->event();
	  LogInfo << "CalibEventLPMT Read in: " << calibeventLPMT << std::endl;
	}

	auto calibheaderSPMT = JM::getHeaderObject<JM::CdSpmtCalibHeader>(nav);
	if (calibheaderSPMT) {
	  calibeventSPMT = calibheaderSPMT->event();
	  LogInfo << "CalibEventSPMT Read in: " << calibeventSPMT << std::endl;
	}
	auto calibheaderWP = JM::getHeaderObject<JM::WpCalibHeader>(nav);
	if (calibheaderWP) {
	  calibeventWP = calibheaderWP->event();
	  LogInfo << "CalibEventWP Read in: " << calibeventWP << std::endl;
	}

	auto elecheaderLPMT = JM::getHeaderObject<JM::CdLpmtElecHeader>(nav);
	if(elecheaderLPMT){eventLPMT = elecheaderLPMT->event();
	  LogInfo << "ElecEventLPMT Read in: " << eventLPMT << std::endl;
	}
	auto elecheaderSPMT = JM::getHeaderObject<JM::CdSpmtElecHeader>(nav);
	if(elecheaderSPMT){eventSPMT = elecheaderSPMT->event();
	  LogInfo << "ElecEventSPMT Read in: " << eventSPMT << std::endl;
	}

	auto triggerheader = JM::getHeaderObject<JM::CdTriggerHeader>(nav);
	if(triggerheader){triggerevent = triggerheader->event();
	  std::cout<<" CD TriggerEvent Read in: " << triggerevent <<std::endl;
	}

	auto Wptriggerheader = JM::getHeaderObject<JM::WpTriggerHeader>(nav);
	if(Wptriggerheader){Wptriggerevent = Wptriggerheader->event();
	  std::cout<<" WP TriggerEvent Read in: " << Wptriggerevent <<std::endl;
	}

        auto wptriggerhdr = JM::getHeaderObject<JM::WpTriggerHeader>(nav);
        if (wptriggerhdr) {
	  std::cout<<" WP TriggerEvent"<<std::endl;
        }

        auto mmtriggerhdr = JM::getHeaderObject<JM::MMTriggerHeader>(nav);
        if (mmtriggerhdr) {
	  std::cout<<" MM TriggerEvent" <<std::endl;
        }

	auto OEChdr = JM::getHeaderObject<JM::OecHeader>(nav);
	if(OEChdr){
	  std::cout<<" OEC Header"<<std::endl;
	  oecevt = dynamic_cast<JM::OecEvt*>(OEChdr->event("JM::OecEvt"));
	}

	clearAllTrees(); //Clearing all trees variables

	if(oecevt)
	{
        m_oecevt_E = oecevt->getEnergy();
        m_oecevt_totCharge = oecevt->getTotalCharge();
        m_oecevt_X = oecevt->getVertexX();
        m_oecevt_Y = oecevt->getVertexY();
        m_oecevt_Z = oecevt->getVertexZ();
        m_mytag = oecevt->getTag();
        m_EventTag.clear();
        if((m_mytag & i_pBiPo212pair)== i_pBiPo212pair) m_EventTag.push_back(1);
        if((m_mytag & i_dBiPo212pair)== i_dBiPo212pair) m_EventTag.push_back(2);
        if((m_mytag & i_pBiPo214pair)== i_pBiPo214pair) m_EventTag.push_back(3);
        if((m_mytag & i_dBiPo214pair)== i_dBiPo214pair) m_EventTag.push_back(4);
        if((m_mytag & i_pIBD)== i_pIBD) m_EventTag.push_back(5);
        if((m_mytag & i_dIBD)== i_dIBD) m_EventTag.push_back(6);
        if((m_mytag & i_SingleSE) == i_SingleSE) m_EventTag.push_back(7);
        if((m_mytag & i_SingleME) == i_SingleME) m_EventTag.push_back(8);
        if((m_mytag & i_SingleLE) == i_SingleLE) m_EventTag.push_back(9);
        if((m_mytag & i_pSpallNeutron) == i_pSpallNeutron) m_EventTag.push_back(10);
        if((m_mytag & i_dSpallNeutron) == i_dSpallNeutron) m_EventTag.push_back(11);
        if((m_mytag & i_pB12) == i_pB12)  m_EventTag.push_back(12);
        if((m_mytag & i_dB12) == i_dB12) m_EventTag.push_back(13);
        if(m_tagsvc->isMuon(oecevt)) m_EventTag.push_back(14);
        if(m_tagsvc->isWPMuon(oecevt)) m_EventTag.push_back(15);
	}

	PMT_R = 35.4; //m
	LS_R = 17.7; //m
	RfrIndxLS = 1.5;
	RfrIndxWR = 1.355;

	c = 299792458.0; //m/s

	IDService* idservice = IDService::getIdServ(); //I think this is needed to retrieve the gcu number
	idservice->init();

	// char TriggerName[100];
	if(triggerevent)
	{
	    const auto& type = triggerevent->triggerType();
	    const auto& pmtFired = triggerevent->nHitMultiplicity();
	    //const auto& volID = triggerevent->volumeId();
	    const auto& trigTime = triggerevent->triggerTime();

	    //std::cout<< "Trigger type size " <<type.size()<<std::endl;
	    //std::cout<< "Triggered PMT size " <<pmtFired<<std::endl;
	    //std::cout<<"Trigger time "<<trigTime.GetNanoSec()<<std::endl;
	    m_Trigger = trigTime.GetNanoSec();
	    //std::cout<< "volume size " <<volID.size()<<std::endl;
	    //std::cout<< "Trigger Time size " <<trigTime.size()<<std::endl;
	    m_TriggerSize = type.size();
	    // if (!type.empty() && !type[0].empty()) {
	    //   strncpy(m_TriggerName, type[0].c_str(), sizeof(m_TriggerName) - 1);
	    //   m_TriggerName[sizeof(m_TriggerName) - 1] = '\0'; // Sécurité pour éviter les dépassements
	    //   //m_TriggerName[0] = type[0][0]; // Premier caractère de la première string
	    //   // m_TriggerName[1] = '\0';       // Terminateur nul
	    // } else {
	    //   m_TriggerName[0] = '\0'; // Valeur par défaut si le vecteur est vide
	    // }

	    std::cout<<"type size "<<type.size()<<std::endl;
		if(type.size() != 0){
			for(auto it = 0; it<type.size(); it++){
				std::cout<<"Trigger type = "<<type[it]<<std::endl;
				m_TriggerName.push_back(type[it]);
			}
		}
	}

	const auto& timestamp = nav->TimeStamp();
	int RunNumber = nav->RunID();
	uint32_t EventNumber = nav->EventID();
	//std::cout<<"RunNumber "<<RunNumber<<" EventNumber "<<EventNumber<<std::endl;
	m_iRun = RunNumber;
	m_EvtID = EventNumber;
	m_AssembleID = nav->AssembleID();
	unsigned long long tempTS =  timestamp.GetSec()*1e9 + timestamp.GetNanoSec() ;
	uint64_t TimeRef_2 =  timestamp.GetSec()*1000000000ULL + timestamp.GetNanoSec() ;
	unsigned long long TimeStampLPMT = (tempTS&0xFFFFFFFF00000000) ;

	m_TimeStamp = TimeRef_2;


	if(simevent && saveSim){
		const auto& tracks = simevent->getTracksVec();
		const auto& pmthit = simevent->getCDHitsVec();
		for (auto it = tracks.begin(); it != tracks.end(); it++){
			m_PDG.push_back((*it)->getPDGID());
			m_Edep.push_back((*it)->getEdep());
			m_Qedep.push_back((*it)->getQEdep());
			if(it == tracks.begin()){
				m_Vtx = (*it)->getInitX();
				m_Vty = (*it)->getInitY();
				m_Vtz = (*it)->getInitZ();
			}
		}
		for (auto it = pmthit.begin(); it != pmthit.end(); it++){
			m_TotSimPE += (*it)->getNPE();
		}
		
		m_ntuple->Fill();
	}



	m_WPnHit=0;
	m_WPTrigTime=0;
	if(Wptriggerevent)
	{
	    const auto& WPnHit = Wptriggerevent->nHitMultiplicity();
	    const auto& WPTrigTime = Wptriggerevent->triggerTime();
	    m_WPnHit = WPnHit;
	    m_WPTrigTime = WPTrigTime.GetNanoSec();
	}

	if(eventLPMT && saveElec)
	{
	    //std::cout<<"eventLPMT"<<std::endl;
	    LogInfo <<"Processing JM:CdLpmtElecEvt "<<std::endl;
	    //m_DateTimeLPMTElec = trigTime64bits;
	    //		cout<<"m_DateTimeLPMTElec="<<m_DateTimeLPMTElec<<endl;
	    //m_TriggerTimeLPMTElec = trigTime48bits;
	    m_TimeStampLPMTElec = timestamp.GetSec()*1e9 + timestamp.GetNanoSec();
	    const auto& elecLPMT = eventLPMT->channelData();
	    //		cout<<"Number of fired LPMT channel in the event "<<elecLPMT.size()<<endl;
	    LogInfo <<"Number of fired LPMT channel in the event "<<elecLPMT.size()<<std::endl;
	    m_NbFiredChannelLPMTElec = elecLPMT.size();
	    for(auto it = elecLPMT.begin();it!=elecLPMT.end();it++)
	    {
			m_GCUNumberLPMTElec = idservice->id2GCU(Identifier(it->first));
			m_ChargeLPMTElec.reserve( ((it->second)->charge()).size() );
			m_HitTimeLPMTElec.reserve( ((it->second)->time()).size() );

			m_CopyNoLPMTElec = idservice->id2CopyNo(Identifier(it->first));
			m_PmtIDLPMTElec = (it->first);
			m_HitTimeLPMTElec.assign( ((it->second)->time()).begin(), ((it->second)->time()).end() );
			m_ChargeLPMTElec.assign( ((it->second)->charge()).begin(), ((it->second)->charge()).end()  );
			m_NbHitLPMTElec+=m_HitTimeLPMTElec.size();

			m_ntuple4->Fill();
	    }
	}

	if(eventSPMT && saveElec)
	{
	    //std::cout<<"eventSPMT"<<std::endl;
	    const auto& elecSPMT = eventSPMT->channelData();
//std::cout<<"SPMT channel Data Map size:"<<elecSPMT.size()<<std::endl;
		m_NbHitSPMTElec=elecSPMT.size();//=m_HitTimeElec.size();
	    for(auto it = elecSPMT.begin();it!=elecSPMT.end();it++)
	      //for(int spmtIdIdx=0;spmtIdIdx<326000;spmtIdIdx++)
	      {
		m_GCUNumber = idservice->id2GCU(Identifier(it->first));
		int channelID = idservice->id2Channel(Identifier(it->first));
		Identifier id = static_cast<Identifier>(it->first);
		int pmtid = idservice->id2CopyNo(id);
		//std::cout<<"gcuid = "<<" channelID="<<channelID<<std::endl;
		//std::cout<<"it->first "<<int(it->first)<<" "<<pmtid<<std::endl;
		m_charge.reserve((it->second)->charge().size()); // Reserving a vector large enough avoid reallocation when push_back hits the last pre-allocated memory case
		m_blockCounter.reserve((it->second)->blockCounter().size());
		m_channelCounter.reserve((it->second)->channelCounter().size());
		m_BECTime.reserve((it->second)->becTime().size());
		m_coarseTime.reserve((it->second)->coarseTime().size());
		m_CTOverflow.reserve((it->second)->ctOverflow().size());
		m_fineTime.reserve((it->second)->fineTime().size());
		m_channelFlag.reserve((it->second)->channelFlag().size());

		m_PmtIdElec = pmtid;
		m_blockID = ((it->second)->blockID());
		m_channelNumber = ((it->second)->channelNumber());
		m_channelID = channelID;
		m_charge.clear();
		m_blockCounter.clear();
		m_channelCounter.clear();
		m_BECTime.clear();
		m_coarseTime.clear();
		m_CTOverflow.clear();
		m_fineTime.clear();
		m_channelFlag.clear();
		for(int l=0;l<((it->second)->charge()).size();l++)
		{
		    if((it->second)->isHighGain(l))
		    	m_gain.push_back(1);
		    else
		    	m_gain.push_back(0);

		    m_charge.push_back( ((it->second)->getCharge(l)));
		    m_blockCounter.push_back( ((it->second)->getBlockCounter(l)));
		    m_channelCounter.push_back( ((it->second)->getchannelCounter(l)));
		    m_BECTime.push_back( ((it->second)->getBECTime(l)));
		    m_coarseTime.push_back( ((it->second)->getCoarseTime(l)));
		    m_CTOverflow.push_back( ((it->second)->getCTOverflow(l)));
		    m_fineTime.push_back( ((it->second)->getFineTime(l)));
		    m_channelFlag.push_back( ((it->second)->isHighGain(l)));
		    long long int TimeSPMT = ((it->second)->getBECTime(l))*pow(2,24)*8 + 25*(((it->second)->getCoarseTime(l)) + ( ((it->second)->getCTOverflow(l))*pow(2,26) ) - ((it->second)->getFineTime(l))/1024 )  ;
		}

		m_ntuple3->Fill();
		//			m_NbHitSPMTElec+=m_HitTimeElec.size();
	      }
	}

	double ChargeTot=0;
	double ChargeTotUp=0;
	std::vector<int> tempPmtIds;
	std::vector<double> tempHitTimes;
	std::vector<double> tempCharges;
	std::vector<int> tempCircleNo;
	std::vector<std::string> tempType;

	double VtxReco=m_oecevt_X;
	double VtyReco=m_oecevt_Y;
	double VtzReco=m_oecevt_Z;
	if(calibeventLPMT){

		std::cout<<"INIT calibevent"<<std::endl;
		const auto& chhlistLPMT = calibheaderLPMT->event()->calibPMTCol();
		// std::cout<<"NLPMT Hitted "<<chhlistLPMT.size()<<std::endl;
		for (auto chit = chhlistLPMT.begin(); chit!=chhlistLPMT.end(); ++chit){
			auto calib = *chit;
			unsigned int pmtId = calib->pmtId();
			Identifier id = Identifier(pmtId);
			//PmtGeom *pG = m_cdGeom->getPmt(id);
			int TruePM=CdID::module(id);
			int CircleNo = CdID::circleNo(id);
			int PmtType = CdID::pmtType(id);

			SniperPtr<IPMTParamSvc> pmtsvc(getParent(), "PMTParamSvc");
			double ZZ = pmtsvc->getPMTX(TruePM);
			//Identifier id = static_cast<Identifier>(it->first);
			//int pmtid = idservice->id2CopyNo(id);
			//std::cout<<"TruePM "<<TruePM<<std::endl;
			m_NbHitLPMTCalib+=calib->size();

			// std::cout<<"DANS CALIB "<<m_iEvt<<" "<<VtxReco<<" "<<VtyReco<<" "<<VtzReco<<std::endl;
			for(unsigned int j=0;j<calib->size();j++){
				tempPmtIds.push_back(TruePM);
				tempCharges.push_back(calib->charge(j));
				tempCircleNo.push_back(CircleNo);
				if(PmtType == 1) tempType.push_back("Hamamatsu");
				else if(PmtType == 2) tempType.push_back("NNVT");

				double timeTOF2 = 0.0;
				//timeTOF2 = calib->time(j) - ComputeLTOF(TruePM, VtxReco, VtyReco, VtzReco);
				timeTOF2 = calib->time(j);
				tempHitTimes.push_back(timeTOF2);
				//tempHitTimes.push_back(calib->time(j));
				//m_PmtIdCalib.push_back(TruePM);
				//m_HitTimeCalib.push_back(calib->time(j));
				//m_ChargeCalib.push_back(calib->charge(j));
				ChargeTot+=calib->charge(j);
				if(ZZ>14000) ChargeTotUp+=calib->charge(j);
			}
		}
		m_ChargeTotLPMT = ChargeTot;
	}
	 //std::cout<<"ChargeTot "<<ChargeTot<<std::endl;
//	 if (ChargeTot > 1000 && ChargeTot < 60000)
	{
		m_PmtIdCalib.insert(m_PmtIdCalib.end(), tempPmtIds.begin(), tempPmtIds.end());
		m_HitTimeCalib.insert(m_HitTimeCalib.end(), tempHitTimes.begin(), tempHitTimes.end());
		m_ChargeCalib.insert(m_ChargeCalib.end(), tempCharges.begin(), tempCharges.end());
		m_CircleNo.insert(m_CircleNo.end(), tempCircleNo.begin(), tempCircleNo.end());
		m_PmtType.insert(m_PmtType.end(), tempType.begin(), tempType.end());
	}

	 ChargeTot=0;
	if(calibeventSPMT)
	{
		const auto& chhlistSPMT = calibheaderSPMT->event()->calibPMTCol();
		for (auto cchit = chhlistSPMT.begin(); cchit!=chhlistSPMT.end(); ++cchit) {
			auto calibSPMT = *cchit;

			IDService* idServ = IDService::getIdServ();
			idServ->init();

			unsigned int pmtId = calibSPMT->pmtId();
			Identifier id = Identifier(pmtId);
			int TruePM=CdID::module(id);
			int CircleNo = CdID::circleNo(id);
			TruePM+=20000-17612;

			m_NbHitSPMTCalib+=calibSPMT->size();
			for(unsigned int j=0;j<calibSPMT->size();j++)
			{
				m_PmtIdCalib.push_back(TruePM);
				double timeTOF2 = 0.0;
				//timeTOF2 = calibSPMT->time(j) - ComputeSTOF(TruePM, VtxReco, VtyReco, VtzReco);
				timeTOF2 = calibSPMT->time(j);
				m_HitTimeCalib.push_back(timeTOF2);
				//m_HitTimeCalib.push_back(calibSPMT->time(j));
				double CalibSPMTBackADC = 124.*0.48*calibSPMT->charge(j)+76.;
				m_ChargeCalib.push_back( CalibSPMTBackADC);
				ChargeTot+=CalibSPMTBackADC;

				m_CircleNo.push_back(CircleNo);
				m_PmtType.push_back("SPMT");
				//m_ChargeCalib.push_back(calibSPMT->charge(j));
				//ChargeTot+=calibSPMT->charge(j);
			}
		}
	}
	 m_ChargeTotSPMT = ChargeTot;
	 //	std::cout<<"CHARGE SPMT "<<m_ChargeTotSPMT<<std::endl;

	ChargeTot=0;
	if(calibeventWP)
	{
	    const auto& chhlistWP = calibheaderWP->event()->calibPMTCol();
	    for (auto wphit = chhlistWP.begin(); wphit!=chhlistWP.end(); ++wphit) {
			auto calibWP = *wphit;
			m_NbHitWPCalib+=calibWP->size();
			for(unsigned int j=0;j<calibWP->size();j++){
				ChargeTot+=calibWP->charge(j);
		 	}
	    }
	}
	 m_ChargeTotWP = ChargeTot;
	 //std::cout<<"CHARGE WP "<<m_ChargeTotWP<<std::endl;

	if (recevent)
	{
		// std::cout << "ADDING REC INFO TO PLAIN" << std::endl;
		const auto& recvertices = recevent->vertices();
		LogInfo << " CdVertexRecEvt: " << std::endl;
		LogInfo << " - number of vertices: " << recevent->nVertices() << std::endl;

		for(auto vertex: recvertices)
		{
			m_TotalPERec = vertex->peSum();
			m_NFiredPMT = vertex->nfiredpmts();
			m_RecE = vertex->energy();
			m_RecX = vertex->x();
			m_RecY = vertex->y();
			m_RecZ = vertex->z();
			m_T0 = vertex->t0();
		}

		if (calibeventLPMT || calibeventSPMT){

			TVector3 vertex(m_RecX, m_RecY, m_RecZ);
			for (size_t i = 0; i < m_PmtIdCalib.size(); i++)
			{
				int pmtid = m_PmtIdCalib.at(i);
				if(pmtid < 17612){
					TOFCalculator TOF(vertex, ALL_LPMT_pos.at(pmtid), interface); //interface is set to -40 m by default
					double timeTOF = m_HitTimeCalib.at(i) - TOF.CalTOF();
					m_HitTimeCalibTOF.push_back(timeTOF);
				}
				else if (pmtid >= 20000){
					TOFCalculator TOF(vertex, ALL_SPMT_pos.at(pmtid - 20000), interface);
					double timeTOF = m_HitTimeCalib.at(i) - TOF.CalTOF();
					m_HitTimeCalibTOF.push_back(timeTOF);
				}
			}

		}
	}
	if(calibeventLPMT || calibeventSPMT || calibeventWP || oecevt) //If there is a trigger fill calib tree
	{
		std::cout << "Fill Calib " << std::endl;
		m_ntuple1->Fill();
	}
	if(oecevt)
	{
		std::cout<<"Fill OEC "<<std::endl;
		m_ntuple5->Fill();
	}
	if(recevent)
	{
		std::cout<<"Fill Rec evt "<<std::endl;
		m_ntuple6->Fill();
	}

	 return true;

}

bool SNiPERToPlainTreeSPMT::Book_tree()
{
	SniperPtr<RootWriter> svc(*getRoot(),"RootWriter");

	if(saveSim){
		m_ntuple = svc->bookTree(*m_par, "Data/Sim", "SimTree");
		m_ntuple->Branch("EntryNb", &m_iEvt, "EntryNb/I");
		m_ntuple->Branch("TotalPE", &m_TotSimPE, "TotalPE/I");
		m_ntuple->Branch("Vtx", &m_Vtx, "Vtx/D");
		m_ntuple->Branch("Vty", &m_Vty, "Vty/D");
		m_ntuple->Branch("Vtz", &m_Vtz, "Vtz/D");
		m_ntuple->Branch("PDGID", &m_PDG);
		m_ntuple->Branch("Edep", &m_Edep);
		m_ntuple->Branch("Qedep", &m_Qedep);
	}
	m_ntuple1 = svc->bookTree(*m_par,"Data/Calib", "CalibTree");
	m_ntuple1->Branch("EntryNb", &m_iEvt, "EntryNb/I");
	m_ntuple1->Branch("RunNb", &m_iRun,"RunNb/I");
	m_ntuple1->Branch("EvtID",&m_EvtID,"EvtID/I");
	m_ntuple1->Branch("AssID",&m_AssembleID,"AssID/I");
	m_ntuple1->Branch("TriggerTime",&m_Trigger,"TriggerTime/l");
	m_ntuple1->Branch("TriggerName", &m_TriggerName);
	m_ntuple1->Branch("TimeStamp",&m_TimeStamp,"TimeStamp/l");
	m_ntuple1->Branch("ChargeTotLPMT",&m_ChargeTotLPMT,"ChargeTotLPMT/D");
	m_ntuple1->Branch("ChargeTotSPMT",&m_ChargeTotSPMT,"ChargeTotSPMT/D");
	m_ntuple1->Branch("ChargeTotWP",&m_ChargeTotWP,"ChargeTotWP/D");
	m_ntuple1->Branch("NbHitLPMTCalib", &m_NbHitLPMTCalib, "NbHitLPMTCalib/I");
	m_ntuple1->Branch("NbHitSPMTCalib", &m_NbHitSPMTCalib, "NbHitSPMTCalib/I");
	m_ntuple1->Branch("NbHitWPCalib", &m_NbHitWPCalib, "NbHitWPCalib/I");
	m_ntuple1->Branch("PmtIDCalib", &m_PmtIdCalib);
	m_ntuple1->Branch("PmtCircle", &m_CircleNo);
	m_ntuple1->Branch("PmtType", &m_PmtType);
	m_ntuple1->Branch("HitTimeCalib", &m_HitTimeCalib);
	m_ntuple1->Branch("HitTimeCalibTOF", &m_HitTimeCalibTOF);
	m_ntuple1->Branch("ChargeCalib", &m_ChargeCalib);

	if(saveElec){
		m_ntuple3 = svc->bookTree(*m_par,"Data/SPMTElec", "ElecTreeSPMT");
		m_ntuple3->Branch("EntryNb", &m_iEvt, "EntryNb/I");
		m_ntuple3->Branch("NbHiSPMTElec", &m_NbHitSPMTElec, "NbHitSPMTElec/I");
		m_ntuple3->Branch("GCUID", &m_GCUNumber,"GCUID/I");
		m_ntuple3->Branch("PmtIDElec", &m_PmtIdElec,"PmtIDElec/I");
		m_ntuple3->Branch("blockID", &m_blockID,"blockID/I");
		m_ntuple3->Branch("channelNumber", &m_channelNumber,"channelNumber/I");
		m_ntuple3->Branch("channelID", &m_channelID,"channelID/I");
		m_ntuple3->Branch("blockCounter", &m_blockCounter);
		m_ntuple3->Branch("gain",&m_gain);
		m_ntuple3->Branch("channelCounter", &m_channelCounter);
		m_ntuple3->Branch("BECTime", &m_BECTime);
		m_ntuple3->Branch("coarseTime", &m_coarseTime);
		m_ntuple3->Branch("CTOverflow", &m_CTOverflow);
		m_ntuple3->Branch("fineTime", &m_fineTime);
		m_ntuple3->Branch("chargeSPMTElec", &m_charge);
		m_ntuple3->Branch("channelFlag", &m_channelFlag);
		m_ntuple3->Branch("NbHitSPMTElec", &m_NbHitSPMTElec, "NbHitSPMTElec/I");

		// m_ntuple4 = svc->bookTree(*m_par,"Data/LPMTElec", "ElecTreeLPMT");
		// m_ntuple4->Branch("EntryNb", &m_iEvt, "EntryNb/I");
		// m_ntuple4->Branch("TriggerTime",&m_Trigger,"TriggerTime/l");
		// m_ntuple4->Branch("TimeStamp",&m_TimeStamp,"TimeStamp/l");
		// m_ntuple4->Branch("NbFiredChannelLPMTElec",&m_NbFiredChannelLPMTElec,"NbFiredChannelLPMTElec/I");
		// m_ntuple4->Branch("GCUNumberLPMTElec", &m_GCUNumberLPMTElec,"GCUNumberLPMTElec/I");
		// m_ntuple4->Branch("PmtIDLPMTElec", &m_PmtIDLPMTElec,"PmtIDLPMTElec/I");
		// m_ntuple4->Branch("CopyNoLPMTElec",&m_CopyNoLPMTElec,"CopyNoLPMTElec/I");
		// m_ntuple4->Branch("HitTimeLPMTElec", &m_HitTimeLPMTElec);
		// m_ntuple4->Branch("ChargeLPMTElec", &m_ChargeLPMTElec);
		// m_ntuple4->Branch("NbHitLPMTElec", &m_NbHitLPMTElec, "NbHitLPMTElec/I");
	}

	m_ntuple5 = svc->bookTree(*m_par,"Data/EventInfo", "EventInfo");
	m_ntuple5->Branch("EntryNb", &m_iEvt, "EntryNb/I");
	m_ntuple5->Branch("RunNb", &m_iRun,"Run/I");
	m_ntuple5->Branch("EvtID",&m_EvtID,"EvtID/I");
	m_ntuple5->Branch("AssID",&m_AssembleID,"AssID/I");
	m_ntuple5->Branch("TriggerTime",&m_TriggerTime,"TriggerTime/l");
	m_ntuple5->Branch("TriggerSize",&m_TriggerSize,"TriggerSize/I");
	m_ntuple5->Branch("TriggerName",&m_TriggerName);
	m_ntuple5->Branch("TimeStamp",&m_TimeStamp,"TimeStamp/l");
	m_ntuple5->Branch("OECEnergy",&m_oecevt_E,"OECEnergy/D");
	m_ntuple5->Branch("OECTotCharge",&m_oecevt_totCharge,"OECTotCharge/D");
	m_ntuple5->Branch("OECX",&m_oecevt_X,"OECX/D");
	m_ntuple5->Branch("OECY",&m_oecevt_Y,"OECY/D");
	m_ntuple5->Branch("OECZ",&m_oecevt_Z,"OECZ/D");
    m_ntuple5->Branch("OECTag",&m_mytag,"OECTag/I");
	m_ntuple5->Branch("TagEvent",&m_EventTag);
	m_ntuple5->Branch("WPnHit",&m_WPnHit,"WPnHit/I");
	m_ntuple5->Branch("WPTrigTime",&m_WPTrigTime,"WPTrigTime/l");

	m_ntuple6 = svc->bookTree(*m_par, "Data/Reco", "Reconstruction Tree");
	m_ntuple6->Branch("EntryNb", &m_iEvt, "EntryNb/I");
	m_ntuple6->Branch("NFiredPMT", &m_NFiredPMT, "NFiredPMT/I");
	m_ntuple6->Branch("TotalPERec", &m_TotalPERec);
	m_ntuple6->Branch("RecEnergy", &m_RecE, "RecEnergy/D");
	m_ntuple6->Branch("Recx", &m_RecX, "Recx/D");
	m_ntuple6->Branch("Recy", &m_RecY, "Recy/D");
	m_ntuple6->Branch("Recz", &m_RecZ, "Recz/D");
	m_ntuple6->Branch("RecT0", &m_T0, "RecT0/D");



	return true;
}
void SNiPERToPlainTreeSPMT::clearAllTrees()
{
	m_TotSimPE = 0;
	m_Vtx = 0.0;
	m_Vty = 0.0;
	m_Vtz = 0.0;
	m_PDG.clear();
	m_Edep.clear();
	m_Qedep.clear();

	m_ChargeTotLPMT=0;
	m_ChargeTotUpLPMT=0;
	m_ChargeTotSPMT=0;
	m_ChargeTotWP=0;
	m_TimeStampInNanoSec=0;
	m_NbHitLPMTCalib=0;
	m_NbHitSPMTCalib=0;
	m_NbHitWPCalib=0;
	m_PmtIdCalib.clear();
	m_HitTimeCalib.clear();
	m_HitTimeCalibTOF.clear();
	m_ChargeCalib.clear();
	m_CircleNo.clear();
	m_PmtType.clear();
	m_NbHitLPMTElec=0;
	m_NbHitSPMTElec=0;
	m_GCUNumber=0;
	m_PmtIdElec=0;
	m_blockID=0;
	m_channelNumber=0;
	m_channelID=0;
	m_gain.clear();
	m_blockCounter.clear();
	m_channelCounter.clear();
	m_BECTime.clear();
	m_coarseTime.clear();
	m_CTOverflow.clear();
	m_fineTime.clear();
	m_charge.clear();
	m_channelFlag.clear();
	m_TimeStampLPMTElec=0;
	m_NbFiredChannelLPMTElec=0; //Number of different channels fired in the event
	m_GCUNumberLPMTElec=0; //GCU Number for LPMT
	m_CopyNoLPMTElec=0;
	m_PmtIDLPMTElec=0;    // PmtID
	m_HitTimeLPMTElec.clear(); // vector of Hit Time of LPMT
	m_ChargeLPMTElec.clear(); // vector of charge of LPMT
	m_TriggerTime=0;
	m_TriggerSize=0;
	//   m_TriggerName[100];
	m_TriggerName.clear();
	m_oecevt_E=0;
	m_oecevt_totCharge=0;
	m_oecevt_X=0;
	m_oecevt_Y=0;
	m_oecevt_Z=0;
	m_mytag=0;
	m_EventTag.clear();
	m_WPnHit=0;
	m_WPTrigTime=0;
	m_TotalPERec=0;
	m_NFiredPMT=0;
	m_RecE=0;
	m_RecX=0;
	m_RecY=0;
	m_RecZ=0;
	m_T0=0;

	PMT_R=0;
	LS_R=0;

	RfrIndxLS=0;
	RfrIndxWR=0;
	c=0;
}
bool SNiPERToPlainTreeSPMT::LoadCalibValues(const std::string& txtFileName) {
    std::ifstream infile(txtFileName);
    if (!infile.is_open()) {
        std::cerr << "Erreur : impossible d’ouvrir le fichier " << txtFileName << std::endl;
        return false;
    }

    int evtID;
    double E1, E2,E3;
    while (infile >> evtID >> E1 >> E2 >> E3) {
      fCalibValues[evtID] = CalibValues{E1,E2,E3};
      std::cout<<"lecture file "<<evtID<<" E1: "<<E1<<" "<<E2<<" "<<E3<<std::endl;
    }

    infile.close();
    return true;
}


bool SNiPERToPlainTreeSPMT::finalize()
{
	LogDebug << "finalizing" << std::endl;
	return true;
}


// Can be removed if TOFCalculator class is kept
double SNiPERToPlainTreeSPMT::ComputeLTOF(double pmtid, double evtx, double evty, double evtz){
  double pmt_pos_x = ALL_LPMT_pos.at(pmtid).X();
  double pmt_pos_y = ALL_LPMT_pos.at(pmtid).Y();
  double pmt_pos_z = ALL_LPMT_pos.at(pmtid).Z();

  double dx = (pmt_pos_x - evtx);
  double dy = (pmt_pos_y - evty);
  double dz = (pmt_pos_z - evtz);

  double Evt = sqrt(evtx*evtx + evty*evty + evtz*evtz);
  double Dist = sqrt(dx*dx + dy*dy + dz*dz);
  double costheta = (Dist*Dist + PMT_R*PMT_R*1e6 - Evt*Evt)/(2.*Dist*PMT_R*1e3); //Al Kashi
  double LengthWater = 1e3*PMT_R*costheta - 1e3*sqrt(PMT_R*costheta*PMT_R*costheta - PMT_R*PMT_R + LS_R*LS_R);

  return RfrIndxLS*(Dist-LengthWater)*1e6/c + RfrIndxWR*LengthWater*1e6/c;
}

double SNiPERToPlainTreeSPMT::ComputeSTOF(double pmtid, double evtx, double evty, double evtz){
  pmtid = pmtid - 20000;
  double pmt_pos_x = ALL_SPMT_pos.at(pmtid).X();
  double pmt_pos_y = ALL_SPMT_pos.at(pmtid).Y();
  double pmt_pos_z = ALL_SPMT_pos.at(pmtid).Z();

  double dx = (pmt_pos_x - evtx);
  double dy = (pmt_pos_y - evty);
  double dz = (pmt_pos_z - evtz);

  double Evt = sqrt(evtx*evtx + evty*evty + evtz*evtz);
  double Dist = sqrt(dx*dx + dy*dy + dz*dz);
  double costheta = (Dist*Dist + PMT_R*PMT_R*1e6 - Evt*Evt)/(2.*Dist*PMT_R*1e3); //Al Kashi
  double LengthWater = 1e3*PMT_R*costheta - 1e3*sqrt(PMT_R*costheta*PMT_R*costheta - PMT_R*PMT_R + LS_R*LS_R);

  return RfrIndxLS*(Dist-LengthWater)*1e6/c + RfrIndxWR*LengthWater*1e6/c;
}
