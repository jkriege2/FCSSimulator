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

#include "nulldynamics.h"
#include "ticktock.h"
#include "tools.h"


NullDynamics::NullDynamics(FluorophorManager* fluorophors, std::string object_name):
    FluorophorDynamics(fluorophors, object_name)
{
    has_walker=false;
    x=y=z=0;
}

NullDynamics::~NullDynamics()
{
}

void NullDynamics::read_config_internal(jkINIParser2& parser) {
    FluorophorDynamics::read_config_internal(parser);

    has_walker=parser.getAsBool("has_walker", has_walker);
    x=parser.getAsDouble("walker_x", x);
    y=parser.getAsDouble("walker_y", y);
    z=parser.getAsDouble("walker_z", z);
}

void NullDynamics::init() {
    FluorophorDynamics::init();
    if (has_walker) change_walker_count(1,n_fluorophores);
    else change_walker_count(0,n_fluorophores);
    sim_time=0;
    endoftrajectory=false;
    set_sim_timestep(sim_timestep);

    for (unsigned long i=0; i<walker_count; i++) {
        init_walker(i, x, y, z);
    }
}



void NullDynamics::propagate(bool boundary_check) {
    FluorophorDynamics::propagate(boundary_check);
    for (unsigned long i=0; i<walker_count; i++) {
        walker_state[i].x=x;
        walker_state[i].y=y;
        walker_state[i].z=z;
    }
}


