/*
    Copyright (c) 2008-2015 Jan W. Krieger (<jan@jkrieger.de>), German Cancer Research Center + IWR, University Heidelberg

    This software is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


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
