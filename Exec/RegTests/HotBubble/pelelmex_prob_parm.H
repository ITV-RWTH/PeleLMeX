#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real T_mean = 300.0_rt;
  amrex::Real T_bubble = 600.0_rt;
  amrex::Real bubble_rad = 0.005_rt;
  amrex::Real bubble_y0 = 0.01;
  int bubble_is_mix = 0;
  int is_sym = 0;
  amrex::Real T_center = 298.0_rt;
  amrex::Real T_coflow = 298.0_rt;
  amrex::Real T_wall = 298.0_rt;
  amrex::Real V_mean = 1.0_rt;
  amrex::Real V_coflow = 1.0_rt;
  amrex::Real phi = 1.0_rt;
  amrex::Real slot_width = 0.01_rt;
  amrex::Real wall_thickness =  slot_width / 10.0_rt;
  amrex::Real wall_height = slot_width;
};
#endif
