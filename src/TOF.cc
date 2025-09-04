#include "TOF.h"

TOFCalculator::TOFCalculator(TVector3 Vertex, TVector3 PMT, double interface)
{
    params.vertex = Vertex;
    params.pmt = PMT;
    params.z_interface = interface;
    // std::cout << "LS/Water interface: " << params.z_interface << std::endl;
    params.n_LS = 1.49;
    params.n_Water = 1.355;

    c = 299792458.0; //m/s
    PMT_R = 35.4; //m
	LS_R = 17.7; //m
}

TOFCalculator::~TOFCalculator()
{
}


double TOFCalculator::TOFFunction(const double* xy){
    TVector3 r(xy[0], xy[1], params.z_interface);

    double dLS = (r - params.vertex).Mag();
    double dWater = (params.pmt - r ).Mag();

    double t1 = dLS/(c*1e-6 / params.n_LS); //ns
    double t2 = dWater/(c*1e-6 / params.n_Water);

    return t1 + t2;
}

double TOFCalculator::TOFMinimizer(){

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    if(!min){
        std::cerr << "Error: cannot create minimizer Minuit2" << std::endl;
        return 0;
    }

    ROOT::Math::Functor fcn([this](const double* x) {
    return this->TOFFunction(x);
    }, 2);

    min->SetMaxFunctionCalls(100000);
    min->SetMaxIterations(10000);
    min->SetTolerance(0.001);
    min->SetPrintLevel(1);

    min->SetFunction(fcn);

    min->SetVariable(0, "x", 0.5 * (params.vertex.X() + params.pmt.X()), 0.01); // first guess in the middle
    min->SetVariable(1, "y", 0.5 * (params.vertex.Y() + params.pmt.Y()), 0.01);

    min->Minimize();

    const double* xs = min->X();
    TVector3 interface(xs[0], xs[1], params.z_interface);

    double tof = TOFFunction(xs);

    std::cout << "Refraction point: " << interface.X() << ", " << interface.Y() << ", " << interface.Z() << std::endl;
    std::cout << "Time of flight: " << tof  << "ns" << std::endl;

    delete min;
    return tof;

}

double TOFCalculator::CalTOF(){

    double dx = (params.pmt.X() - params.vertex.X());
    double dy = (params.pmt.Y() - params.vertex.Y());
    double dz = (params.pmt.Z() - params.vertex.Z());

    double tof = 0.0;

    if(params.pmt.Z() < params.z_interface){ // find position on the interface that gives the fastest traveling path
        tof = TOFMinimizer();
    }
    else{ // https://arxiv.org/pdf/1803.09394
        double Evt = sqrt(params.vertex.X()*params.vertex.X() + params.vertex.Y()*params.vertex.Y() + params.vertex.Z()*params.vertex.Z());
	    double Dist = sqrt(dx*dx + dy*dy + dz*dz);
	    double costheta = (Dist*Dist + PMT_R*PMT_R*1e6 - Evt*Evt)/(2.*Dist*PMT_R*1e3); //Al Kashi
	    double LengthWater = 1e3*PMT_R*costheta - 1e3*sqrt(PMT_R*costheta*PMT_R*costheta - PMT_R*PMT_R + LS_R*LS_R);

        tof = params.n_LS*(Dist-LengthWater)*1e6/c + params.n_Water*LengthWater*1e6/c;
    }

    return tof;
}


