#include "fretdynamics.h"
#include "ticktock.h"
#include "tools.h"


FRETDynamics::FRETDynamics(FluorophorManager* fluorophors, std::string object_name):
    FluorophorDynamics(fluorophors, object_name)
{
    parent="";
}

FRETDynamics::~FRETDynamics()
{
}

void FRETDynamics::read_config_internal(jkINIParser2& parser) {
    FluorophorDynamics::read_config_internal(parser);

    parent=tolower(strstrip(parser.getAsString("parent", parent)));
    n_fluorophores=1;
}

void FRETDynamics::init() {
    FluorophorDynamics* p=get_parent();
    if (p) {
        p->ensure_dynamics_is_hooked(this);
        handle_parent_walker_count_changed(p->get_walker_count(), n_fluorophores);
        sim_time=p->get_sim_time();
        endoftrajectory=p->end_of_trajectory_reached();
        set_sim_timestep(p->get_sim_timestep());
    } else {
        change_walker_count(0,n_fluorophores);
        sim_time=0;
        endoftrajectory=false;
        set_sim_timestep(sim_timestep);
    }

    init_additional_walkers();

    store_step_protocol();

}

void FRETDynamics::handle_parent_walker_count_changed(unsigned long N_walker, unsigned long N_fluorophores) {
    FluorophorDynamics* p=get_parent();
    if (p) {
        change_walker_count(N_walker, n_fluorophores);
    } else {
        change_walker_count(0,n_fluorophores);
    }

    init_additional_walkers();

    propagate(false);
}

void FRETDynamics::propagate(bool boundary_check) {
    FluorophorDynamics* p=get_parent();
    if (p) {
        sim_time=p->get_sim_time();
        endoftrajectory=p->end_of_trajectory_reached();
        FluorophorDynamics::walkerState* ws=p->get_walker_state();
        for (unsigned long i=0; i<walker_count; i++) {
            walker_state[i].x=ws[i].x;
            walker_state[i].y=ws[i].y;
            walker_state[i].z=ws[i].z;
            walker_state[i].x0=ws[i].x0;
            walker_state[i].y0=ws[i].y0;
            walker_state[i].z0=ws[i].z0;
            walker_state[i].ix=ws[i].ix;
            walker_state[i].iy=ws[i].iy;
            walker_state[i].iz=ws[i].iz;
            walker_state[i].ix0=ws[i].ix0;
            walker_state[i].iy0=ws[i].iy0;
            walker_state[i].iz0=ws[i].iz0;

            propagate_photophysics(i);
        }
        store_step_protocol();
    }
}

std::string FRETDynamics::report() {
    std::string s=FluorophorDynamics::report();
    s+="parent = "+parent+"\n";
    return s;
}


FluorophorDynamics* FRETDynamics::get_parent() const {
    if (dynmap.count(parent)>0) return dynmap[parent];
    return NULL;
}

bool FRETDynamics::depends_on(const FluorophorDynamics* other) const {
    return other==get_parent();
}
