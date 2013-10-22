#include "childdynamics.h"
#include "ticktock.h"
#include "tools.h"


ChildDynamics::ChildDynamics(FluorophorManager* fluorophors, std::string object_name):
    FluorophorDynamics(fluorophors, object_name)
{
    parent="";
    initial_walker_visible=true;
    dont_copy_photophysics=false;
}

ChildDynamics::~ChildDynamics()
{
}

void ChildDynamics::read_config_internal(jkINIParser2& parser) {
    FluorophorDynamics::read_config_internal(parser);

    parent=tolower(strstrip(parser.getAsString("parent", parent)));
    initial_walker_visible=parser.getAsBool("initial_walker_visible", initial_walker_visible);
    dont_copy_photophysics=parser.getAsBool("dont_copy_photophysics", dont_copy_photophysics);
}

void ChildDynamics::init() {
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

void ChildDynamics::handle_parent_walker_count_changed(unsigned long N_walker, unsigned long N_fluorophores) {
    FluorophorDynamics* p=get_parent();
    if (p) {
        change_walker_count(N_walker, n_fluorophores);
    } else {
        change_walker_count(0,n_fluorophores);
    }

    init_additional_walkers();

    propagate(false);
}

void ChildDynamics::propagate(bool boundary_check) {
    FluorophorDynamics* p=get_parent();
    if (p) {
        sim_time=p->get_sim_time();
        endoftrajectory=p->end_of_trajectory_reached();
        FluorophorDynamics::walkerState* ws=p->get_walker_state();
        if (dont_copy_photophysics) {
            FluorophorDynamics* p=get_parent();
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
        } else {
            for (unsigned long i=0; i<walker_count; i++) {
                walker_state[i]=ws[i];
                propagate_photophysics(i);
            }
        }
    }
}

std::string ChildDynamics::report() {
    std::string s=FluorophorDynamics::report();
    s+="parent = "+parent+"\n";
    s+="initial_walker_visible = "+booltostr(initial_walker_visible)+"\n";
    s+="dont_copy_photophysics = "+booltostr(dont_copy_photophysics)+"\n";
    return s;
}


FluorophorDynamics* ChildDynamics::get_parent() const {
    if (dynmap.count(parent)>0) return dynmap[parent];
    return NULL;
}

bool ChildDynamics::depends_on(const FluorophorDynamics* other) const {
    return other==get_parent();
}

unsigned long ChildDynamics::get_visible_walker_count()  {
    if (!initial_walker_visible) return (walker_count)*(n_fluorophores-1);
    return walker_count*n_fluorophores;
}

double ChildDynamics::get_visible_walker_sigma_times_qfl(unsigned long i) {
    if (!initial_walker_visible) return FluorophorDynamics::get_walker_sigma_times_qfl(i+walker_count);
    return FluorophorDynamics::get_walker_sigma_times_qfl(i);

}

FluorophorDynamics::walkerState* ChildDynamics::get_visible_walker_state() {
    if (!initial_walker_visible) {
        return &(walker_state[walker_count]);
    }
    return walker_state;
}
