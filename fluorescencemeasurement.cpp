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


#include "fluorescencemeasurement.h"

FluorescenceMeasurement::FluorescenceMeasurement(FluorophorManager* fluorophors, std::string objectname):
    TickTock()
{
    this->fluorophors=fluorophors;
    this->object_name=objectname;
    sim_time=0;
    sim_timestep=1e-6;
    description=objectname;
    object_number=0;
    group="";
    supergroup="";
     // init GSL random number generator
    gsl_rng_env_setup();
    rng_type = gsl_rng_taus;
    rng = gsl_rng_alloc (rng_type);
    gsl_rng_set(rng, gsl_rng_get(global_rng));
}

FluorescenceMeasurement::~FluorescenceMeasurement()
{
    gsl_rng_free(rng);
    clearDynamics();
}

void FluorescenceMeasurement::read_config_internal(jkINIParser2& parser) {
}

void FluorescenceMeasurement::read_config(jkINIParser2& parser, std::string group, std::string supergroup) {
    this->group=group;
    this->supergroup=supergroup;
    basename=parser.getSetAsString("simulation.basename", "");
    std::string rng=tolower(parser.getSetAsString("simulation.rng", "taus"));
    duration=parser.getAsDouble("simulation.duration", 1.0);
    sim_timestep=parser.getAsDouble("simulation.timestep", 1e-6);
    if (rng=="mt19937") {
        rng_type=gsl_rng_mt19937;
    } else if (rng=="ranlxs0") {
        rng_type=gsl_rng_ranlxs0;
    } else if (rng=="ranlxs1") {
        rng_type=gsl_rng_ranlxs1;
    } else if (rng=="ranlxs2") {
        rng_type=gsl_rng_ranlxs2;
    } else if (rng=="ranlxd1") {
        rng_type=gsl_rng_ranlxd1;
    } else if (rng=="ranlxd2") {
        rng_type=gsl_rng_ranlxd2;
    } else if (rng=="ranlux") {
        rng_type=gsl_rng_ranlux;
    } else if (rng=="ranlux389") {
        rng_type=gsl_rng_ranlux389;
    } else if (rng=="cmrg") {
        rng_type=gsl_rng_cmrg;
    } else if (rng=="mrg") {
        rng_type=gsl_rng_mrg;
    } else if (rng=="taus") {
        rng_type=gsl_rng_taus;
    } else if (rng=="taus2") {
        rng_type=gsl_rng_taus2;
    } else if (rng=="gfsr4") {
        rng_type=gsl_rng_gfsr4;
    } else if (rng=="rand") {
        rng_type=gsl_rng_rand;
    } else if (rng=="bsd") {
        rng_type=gsl_rng_random_bsd;
    } else if (rng=="libc5") {
        rng_type=gsl_rng_random_libc5;
    } else if (rng=="glibc2") {
        rng_type=gsl_rng_random_glibc2;
    } else if (rng=="rand48") {
        rng_type=gsl_rng_rand48;
    } else if (rng=="ranf") {
        rng_type=gsl_rng_ranf;
    } else if (rng=="ranmar") {
        rng_type=gsl_rng_ranmar;
    } else if (rng=="r250") {
        rng_type=gsl_rng_r250;
    } else if (rng=="tt800") {
        rng_type=gsl_rng_tt800;
    } else if (rng=="minstd") {
        rng_type=gsl_rng_minstd;
    } else if (rng=="knuthran2") {
        rng_type=gsl_rng_knuthran2;
    } else if (rng=="knuthran2002") {
        rng_type=gsl_rng_knuthran2002;
    } else if (rng=="knuthran") {
        rng_type=gsl_rng_knuthran;
    }
    if (this->rng!=NULL) gsl_rng_free(this->rng);
    this->rng = gsl_rng_alloc(rng_type);
    gsl_rng_set(this->rng, gsl_rng_get(global_rng));


    for (int i=0; i<2; i++) {
        if (i==0) {
            if (supergroup.size()<=0) break;
            parser.enterGroup(supergroup);
        } else if (i==1) {
            parser.enterGroup(group);
        }

        if (parser.exists("rng_seed")) {
            gsl_rng_set(this->rng, parser.getAsInt("rng_seed", gsl_rng_get(global_rng)));
        }

        read_config_internal(parser);
        description=parser.getAsString("description", description);

        parser.leaveGroup();
    }
    std::cout<<object_name<<" RNG-TEST: "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<"\n";

    //init();
}


std::string FluorescenceMeasurement::report() {
    std::string s="";
    s+="object_name = "+object_name+"\n";
    s+="description = "+description+"\n";
    s+="runtime = "+floattostr(runtime)+" secs\n";
    s+="rng = "+std::string(gsl_rng_name(rng))+"\n";
    s+="data from: ";
    for (size_t i=0; i<dyn.size(); i++) {
        if (i>0) s+=", ";
        s+=dyn[i]->get_object_name();
    }
    s+="\n";
    return s;
}

void FluorescenceMeasurement::init() {
    sim_time=0;
    std::cout<<object_name<<" RNG-TEST: "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<"\n";
}

void FluorescenceMeasurement::propagate() {
    sim_time=sim_time+sim_timestep;
    //std::cout<<">>> meas   sim_time = "<<sim_time<<"\n";

}


void FluorescenceMeasurement::finalize_sim() {

}


void FluorescenceMeasurement::save_results() {
    save();
}


bool FluorescenceMeasurement::depends_on(const FluorescenceMeasurement* other) const {
    return false;
}

std::string FluorescenceMeasurement::dot_get_properties()  {
    std::string s="";
    s+="description = "+description+"<BR/>";
    s+="rng = "+std::string(gsl_rng_name(rng))+"<BR/>";
    return s;
}

std::string FluorescenceMeasurement::dot_get_node(bool simple)  {
    std::string s;
    s+="     "+get_group()+" [ label =<<FONT POINT-SIZE=\"10\"><B>"+get_group()+": "+object_name+"</B></FONT><BR/><BR/><FONT POINT-SIZE=\"7\">";
    if (!simple) s+=dot_get_properties();
    s+="</FONT><BR/><BR/>> color=red4 fontcolor=red4 fillcolor=gray75 ];\n";
    return s;
}
std::string FluorescenceMeasurement::dot_get_links()  {
    std::string s;
    for (size_t i=0; i<dyn.size(); i++) {
        s+="     "+dyn[i]->get_group()+" -> "+get_group()+"  [ color=red4 fontcolor = \"red\" label = \"trajectories\" ];\n";
    }
    s+="\n";
    return s;
}
