#include <PeleLMeX.H>
#include <iostream>
#include <fstream>
using namespace amrex;

void
PeleLM::Regler(double& T_center, double& T_coflow, double& V_mean, double& V_coflow, double& phi){
//    V_mean = V_mean * 1.1 
}
void
PeleLM::Feedback()
{
    double T_center       = PeleLM::prob_parm->T_center;
    double T_coflow       = PeleLM::prob_parm->T_coflow;
    double V_mean         = PeleLM::prob_parm->V_mean;
    double V_coflow       = PeleLM::prob_parm->V_coflow;
    double phi            = PeleLM::prob_parm->phi;

    Regler(T_center,T_coflow,V_mean,V_coflow,phi);

    std::ofstream controlFile("Control.inp");
    if (!controlFile.is_open()) {
        std::cerr << "Error: Control.inp can't be opened." << std::endl;
    }

    controlFile << "prob.T_center ="       << T_center       << "\n";
    controlFile << "prob.T_coflow ="       << T_coflow       << "\n";
    controlFile << "prob.V_mean ="         << V_mean         << "\n";
    controlFile << "prob.V_coflow ="       << V_coflow       << "\n";
    controlFile << "prob.phi ="            << phi            << "\n";

    controlFile.close();

}