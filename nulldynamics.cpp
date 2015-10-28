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
}

NullDynamics::~NullDynamics()
{
}

void NullDynamics::read_config_internal(jkINIParser2& parser) {
    FluorophorDynamics::read_config_internal(parser);

}

void NullDynamics::init() {
    FluorophorDynamics::init();
    change_walker_count(0,n_fluorophores);
    sim_time=0;
    endoftrajectory=false;
    set_sim_timestep(sim_timestep);

}



void NullDynamics::propagate(bool boundary_check) {
    FluorophorDynamics::propagate(boundary_check);

}


