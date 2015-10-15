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


#include "fluorophordynamics.h"

//#include <boost/thread/mutex.hpp>

#include "datatable.h"
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include "fluorescencemeasurement.h"

double RelativeAbsorbanceReader::get_relative_absorbance_for(FluorophorDynamics* dyn, int i) {
    const FluorophorDynamics::walkerState* dynn=dyn->get_visible_walker_state();
    if (i<0 || i>=dyn->get_visible_walker_count()) {
        return 0;
    }
    const FluorophorDynamics::walkerState* walker=&(dynn[i]);
    return get_relative_absorbance_for(dyn, i, walker->x, walker->y, walker->z);
}

FluorophorDynamics::walkerState::walkerState(FluorophorDynamics::walkerState& other) {
    operator=(other);
}

FluorophorDynamics::walkerState& FluorophorDynamics::walkerState::operator=(FluorophorDynamics::walkerState& other) {
    exists=other.exists;
    time=other.time;
    x=other.x;
    y=other.y;
    z=other.z;
    x0=other.x0;
    y0=other.y0;
    z0=other.z0;

    ix=other.ix;
    iy=other.iy;
    iz=other.iz;
    ix0=other.ix0;
    iy0=other.iy0;
    iz0=other.iz0;
    was_just_reset=other.was_just_reset;

    qm_state=other.qm_state;
    memcpy(sigma_abs, other.sigma_abs, N_FLUORESCENT_STATES*sizeof(double));
    memcpy(q_fluor, other.q_fluor, N_FLUORESCENT_STATES*sizeof(double));
    memcpy(photophysics_transition, other.photophysics_transition, N_FLUORESCENT_STATES*N_FLUORESCENT_STATES*sizeof(double));
    used_qm_states=other.used_qm_states;
    p_x=other.p_x;
    p_y=other.p_y;
    p_z=other.p_z;
    type=other.type;
    spectrum=other.spectrum;

    return *this;
}

FluorophorDynamics::FluorophorDynamics(FluorophorManager* fluorophors, std::string object_name)
{
    this->fluorophors=fluorophors;
    this->object_name=object_name;
    protocol_trajectories=0;
    sim_time=0;
    sim_timestep=1e-6;
    init_fluorophor="atto488";
    group="";
    supergroup="";

    use_two_walkerstates=false;
    depletion_propability=0;

    endoftrajectory=false;

    absorbance_reader="";
    photophysics_absorbance_dependent=false;
    photophysics_absorbance_factor=1;

     // init GSL random number generator
    gsl_rng_env_setup();
    rng_type = gsl_rng_taus;
    rng = gsl_rng_alloc (rng_type);
    gsl_rng_set(rng, gsl_rng_get(global_rng));

    heatup_steps=0;

    init_p_x=1;
    init_p_y=0;
    init_p_z=0;
    init_qm_state=0;
    init_type=0;
    init_spectrum=-1;
    init_used_qm_states=1;
    reset_qmstate_at_simboxborder=false;

    for (int j=0; j<N_FLUORESCENT_STATES; j++) init_sigma_abs[j]=2.2e-20;
    for (int j=0; j<N_FLUORESCENT_STATES; j++) init_q_fluor[j]=0.1;
    for (int j=0; j<N_FLUORESCENT_STATES*N_FLUORESCENT_STATES; j++) init_photophysics_transition[j]=0;

    init_spectrum=-1;
    use_photophysics=true;
    protocol_timestep_count=-1;

    additional_walker_photophysics=true;
    additional_walker_position_mode=SamePosition;
    additional_walker_sphere_radius=0.1;
    additional_walker_off_if_main_off=false;
    n_fluorophores=1;

    store_walker_statistics=false;
    walker_statistics_averageduration=100e-3;
    walker_statistics.clear();
    walker_statistics.push_back(FluorophorDynamics::walker_statistics_entry());

    walker_state=NULL;
    walker_state_other=NULL;
    walker_dx=NULL;
    walker_dy=NULL;
    walker_dz=NULL;
    this->sim_x=14;
    this->sim_y=14;
    this->sim_z=4;
    this->sim_radius=10;
    this->volume_shape=Box;
    set_c_fluor(0.1);
    set_sim_timestep(1e-6);
    init();
}




FluorophorDynamics::FluorophorDynamics(FluorophorManager* fluorophors, double sim_x, double sim_y, double sim_z, double c_fluor, std::string object_name)
{
    this->fluorophors=fluorophors;
    this->object_name=object_name;
    protocol_trajectories=0;
    endoftrajectory=false;
    sim_time=0;
    sim_timestep=1e-6;
    use_two_walkerstates=false;
    init_fluorophor="atto488";
    group="";
    supergroup="";

     // init GSL random number generator
    gsl_rng_env_setup();
    rng_type = gsl_rng_taus;
    rng = gsl_rng_alloc (rng_type);
    gsl_rng_set(rng, gsl_rng_get(global_rng));
    depletion_propability=0;
    reset_qmstate_at_simboxborder=false;

    absorbance_reader="";
    photophysics_absorbance_dependent=false;
    photophysics_absorbance_factor=1;


    init_p_x=1;
    init_p_y=0;
    init_p_z=0;
    init_qm_state=0;
    init_type=0;
    init_spectrum=-1;
    init_used_qm_states=1;
    heatup_steps=0;

    for (int j=0; j<N_FLUORESCENT_STATES; j++) init_sigma_abs[j]=2.2e-20;
    for (int j=0; j<N_FLUORESCENT_STATES; j++) init_q_fluor[j]=0.1;
    for (int j=0; j<N_FLUORESCENT_STATES*N_FLUORESCENT_STATES; j++) init_photophysics_transition[j]=0;
    use_photophysics=true;
    protocol_timestep_count=-1;
    additional_walker_photophysics=true;
    additional_walker_position_mode=SamePosition;
    additional_walker_sphere_radius=0.1;
    additional_walker_off_if_main_off=false;
    n_fluorophores=1;

    store_walker_statistics=false;
    walker_statistics_averageduration=100e-3;
    walker_statistics.clear();
    walker_statistics.push_back(FluorophorDynamics::walker_statistics_entry());

    walker_state=NULL;
    walker_state_other=NULL;
    walker_dx=NULL;
    walker_dy=NULL;
    walker_dz=NULL;
    this->sim_x=sim_x;
    this->sim_y=sim_y;
    this->sim_z=sim_z;
    this->sim_radius=sim_z;
    this->volume_shape=Box;
    set_c_fluor(c_fluor);
    set_sim_timestep(1e-6);
    trajectories.clear();
    init();
}

FluorophorDynamics::FluorophorDynamics(FluorophorManager* fluorophors, double sim_radius, double c_fluor, std::string object_name)
{
    this->fluorophors=fluorophors;
    this->object_name=object_name;
    protocol_trajectories=0;
    endoftrajectory=false;
    sim_time=0;
    sim_timestep=1e-6;
    use_two_walkerstates=false;
    depletion_propability=0;
    init_fluorophor="atto488";
    group="";
    supergroup="";
    reset_qmstate_at_simboxborder=false;

    absorbance_reader="";
    photophysics_absorbance_dependent=false;
    photophysics_absorbance_factor=1;


     // init GSL random number generator
    gsl_rng_env_setup();
    rng_type = gsl_rng_taus;
    rng = gsl_rng_alloc (rng_type);
    gsl_rng_set(rng, gsl_rng_get(global_rng));

    store_walker_statistics=false;
    walker_statistics_averageduration=100e-3;
    walker_statistics.clear();
    walker_statistics.push_back(FluorophorDynamics::walker_statistics_entry());

    init_p_x=1;
    init_p_y=0;
    init_p_z=0;
    init_qm_state=0;
    init_type=0;
    init_spectrum=-1;
    init_used_qm_states=1;
    heatup_steps=0;

    for (int j=0; j<N_FLUORESCENT_STATES; j++) init_sigma_abs[j]=2.2e-20;
    for (int j=0; j<N_FLUORESCENT_STATES; j++) init_q_fluor[j]=0.1;
    for (int j=0; j<N_FLUORESCENT_STATES*N_FLUORESCENT_STATES; j++) init_photophysics_transition[j]=0;
    use_photophysics=true;
    protocol_timestep_count=-1;
    additional_walker_photophysics=true;
    additional_walker_position_mode=SamePosition;
    additional_walker_sphere_radius=0.1;
    additional_walker_off_if_main_off=false;
    n_fluorophores=1;

    walker_state=NULL;
    walker_state_other=NULL;
    walker_dx=NULL;
    walker_dy=NULL;
    walker_dz=NULL;
    this->sim_x=5;
    this->sim_y=5;
    this->sim_z=5;
    this->sim_radius=sim_radius;
    this->volume_shape=Ball;
    set_c_fluor(c_fluor);
    set_sim_timestep(1e-6);
    trajectories.clear();
    init();
}


FluorophorDynamics::~FluorophorDynamics()
{
    if (walker_state!=NULL) free(walker_state);
    if (use_two_walkerstates && walker_state_other!=NULL) free(walker_state_other);
    gsl_rng_free(rng);
}


unsigned long FluorophorDynamics::calc_walker_count() {
    //std::cout<<object_name<<": calc_walker_count() ...  c_fluor="<<c_fluor<<"  sim_radius="<<sim_radius<<"  sim_x="<<sim_x<<"  sim_y="<<sim_y<<"  sim_z="<<sim_z;
    unsigned long r=0;
    if (volume_shape==Box) {
        r= (unsigned long)round(c_fluor*1e-9*6.022e23*sim_x*1e-5*sim_y*1e-5*sim_z*1e-5);
    } else if (volume_shape==Ball) {
        r= (unsigned long)round(c_fluor*1e-9*6.022e23*4.0*M_PI/3.0*gsl_pow_3(sim_radius*1e-5));
    }

    //std::cout<<"  => result="<<r<<std::endl;
    return r;
}

void FluorophorDynamics::read_config_internal(jkINIParser2& parser) {
    std::string ivshape="sphere";
    if (volume_shape==Box) ivshape="box";
    std::string vshape=tolower(parser.getSetAsString("volume_shape", ivshape));
    n_fluorophores=parser.getSetAsInt("n_fluorophores", n_fluorophores);
    additional_walker_photophysics=parser.getSetAsBool("additional_walker_photophysics", additional_walker_photophysics);
    additional_walker_off_if_main_off=parser.getSetAsBool("additional_walker_off_if_main_off", additional_walker_off_if_main_off);
    additional_walker_sphere_radius=parser.getSetAsDouble("additional_walker_sphere_radius", additional_walker_sphere_radius);

    std::string awpm="same";
    if (additional_walker_position_mode==InSphere) awpm="in_sphere";
    awpm=tolower(parser.getSetAsString("additional_walker_position_mode", awpm));
    if (awpm=="same") additional_walker_position_mode=SamePosition;
    if (awpm=="in_sphere") additional_walker_position_mode=InSphere;

    if (vshape=="box") {
        set_sim_box(parser.getSetAsDouble("sim_x", sim_x), parser.getSetAsDouble("sim_y", sim_y), parser.getSetAsDouble("sim_z", sim_z));
    } else if (vshape=="sphere" || vshape=="ball") {
        set_sim_sphere(parser.getSetAsDouble("sim_radius", sim_radius));
    }
    set_c_fluor(parser.getSetAsDouble("c_fluor", c_fluor));


    init_p_x=parser.getSetAsDouble("init_p_x", init_p_x);
    init_p_y=parser.getSetAsDouble("init_p_y", init_p_y);
    init_p_z=parser.getSetAsDouble("init_p_z", init_p_z);
    //std::string init_fluorophor="atto488";
    //std::cout<<"fluorophor: "<<parser.getAsString("init_fluorophor")<<std::endl;
    if (parser.exists("init_fluorophor")) {
        init_fluorophor=tolower(parser.getSetAsString("init_fluorophor", init_fluorophor));
        std::cout<<"WARNING: found init_fluorophor='"<<init_fluorophor<<"' in group '"<<parser.getGroupName()<<"'!\n           THIS WILL RESET ALL init_q_fluor_X and init_sigma_abs_X DEFINED IN A PARENT CLASS TO ITS DEFAULTS FROM THE FLUOROPHORE DATASET!\n\n";
        if (!fluorophors->fluorophorExists(init_fluorophor)) {
            throw FluorophorException(format("didn't find fluorophor %s in database", init_fluorophor.c_str()));
        }
        for (int j=0; j<N_FLUORESCENT_STATES; j++) {
            init_sigma_abs[j]=fluorophors->getFluorophorData(init_fluorophor).sigma_abs;
            //std::cout<<"fluorophor: "<<init_fluorophor<<"   init_sigma_abs["<<j<<"]="<<init_sigma_abs[j]<<std::endl;
            init_q_fluor[j]=fluorophors->getFluorophorData(init_fluorophor).fluorescence_efficiency;
        }
        /*
        init_tau_fl=fluorophors->getFluorophorData(init_fluorophor).fluorescence_lifetime;
        init_bleaching_propability=fluorophors->getFluorophorData(init_fluorophor).bleaching_propability;
        init_triplet_lifetime=fluorophors->getFluorophorData(init_fluorophor).triplet_lifetime;
        init_triplet_propability=fluorophors->getFluorophorData(init_fluorophor).triplet_propability;
        std::cout<<init_fluorophor<<": "<<init_q_fluor<<", "<<init_tau_fl<<", "<<init_sigma_abs<<std::endl;*/
    }
    //std::cout<<"fluorophor: "<<init_fluorophor<<std::endl;
    init_qm_state=parser.getSetAsInt("init_qm_state", init_qm_state);
    init_type=parser.getSetAsInt("init_type", init_type);

    //init_used_qm_states=1;

    for (int j=0; j<N_FLUORESCENT_STATES; j++) {
        std::string prop1="init_sigma_abs_"+inttostr(j);
        std::string prop2="init_q_fluor_"+inttostr(j);
        //std::cout<<parser.getGroupName()<<prop1<<"("<<parser.exists(prop1)<<") -> init_sigma_abs["<<j<<"]="<<init_sigma_abs[j]<<"     "<<parser.getGroupName()<<prop2<<"("<<parser.exists(prop2)<<") -> init_q_fluor["<<j<<"]="<<init_q_fluor[j]<<std::endl;
        if ((parser.exists(prop1) || parser.exists(prop2)) && (j+1>init_used_qm_states)) init_used_qm_states=j+1;
        init_sigma_abs[j]=parser.getSetAsDouble(prop1, parser.getAsDouble("init_sigma_abs", init_sigma_abs[j]));
        init_q_fluor[j]=parser.getSetAsDouble(prop2, parser.getAsDouble("init_q_fluor", init_q_fluor[j]));
        //std::cout<<"init_sigma_abs["<<j<<"]="<<init_sigma_abs[j]<<std::endl;
        //std::cout<<"init_used_qm_states="<<init_used_qm_states<<std::endl;
        //std::cout<<parser.getGroupName()<<prop1<<"("<<parser.exists(prop1)<<") -> init_sigma_abs["<<j<<"]="<<init_sigma_abs[j]<<"     "<<parser.getGroupName()<<prop2<<"("<<parser.exists(prop2)<<") -> init_q_fluor["<<j<<"]="<<init_q_fluor[j]<<std::endl;
    }

    for (int i=0; i<N_FLUORESCENT_STATES; i++) {
        bool has_ii=false;
        for (int f=0; f<N_FLUORESCENT_STATES; f++) {
            std::string prop="init_photophysics_transition_"+inttostr(i)+"_"+inttostr(f);
            if ((i==f) && parser.exists(prop)) has_ii=true;
            if (parser.exists(prop) && (i+1>init_used_qm_states)) init_used_qm_states=i+1;
            if (parser.exists(prop) && (f+1>init_used_qm_states)) init_used_qm_states=f+1;
            init_photophysics_transition[i*N_FLUORESCENT_STATES+f]=parser.getSetAsDouble(prop, init_photophysics_transition[i*N_FLUORESCENT_STATES+f]);
            //std::cout<<"init_used_qm_states="<<init_used_qm_states<<std::endl;
        }
        if (!has_ii) {
            std::string prop="init_photophysics_transition_"+inttostr(i)+"_"+inttostr(i);
            std::cout<<"calculating probability "<<prop<<": ";
            double sum=0;
            for (int f=0; f<N_FLUORESCENT_STATES; f++) {
                if (i!=f) sum=sum+init_photophysics_transition[i*N_FLUORESCENT_STATES+f];
            }
            if (sum<=1.0) {
                init_photophysics_transition[i*N_FLUORESCENT_STATES+i]=1.0-sum;
                std::cout<<1.0-sum<<std::endl;
            } else {
                std::cout<<1.0-sum<<std::endl;
                throw FluorophorException(format("probability to leave state %d larger than one!!!", i));
            }
        }
    }
    init_used_qm_states=parser.getSetAsInt("init_used_qm_states", init_used_qm_states);

    use_photophysics=parser.getSetAsBool("use_photophysics", use_photophysics);
    reset_qmstate_at_simboxborder=parser.getSetAsBool("reset_qmstate_at_simboxborder", reset_qmstate_at_simboxborder);
    depletion_propability=parser.getSetAsDouble("depletion_propability", depletion_propability);


    photophysics_absorbance_dependent=parser.getSetAsBool("photophysics_absorbance_dependent", photophysics_absorbance_dependent);
    photophysics_absorbance_factor=parser.getSetAsDouble("photophysics_absorbance_factor", photophysics_absorbance_factor);
    absorbance_reader=parser.getSetAsString("absorbance_reader", absorbance_reader);




    std::string spec=init_fluorophor;
    if (parser.exists("init_spectrum")) {
        spec=tolower(parser.getSetAsString("init_spectrum", spec));
        //init_spectrum=-1;
    }
    std::cout<<"spec="<<spec<<"   init_spectrum="<<init_spectrum<<std::endl;
    if (fluorophors->getFluorophorData(spec).spectrum!=-1) {
        init_spectrum=fluorophors->getFluorophorData(spec).spectrum;
    } else {
        std::cout<<std::endl<<std::endl<<"didn't find spectrum for "<<spec<<std::endl;
    }

    //std::cout<<spec<<" ("<<init_spectrum<<"): "<<init_q_fluor[0]<<", "<<init_sigma_abs[0]<<<<init_q_fluor[1]<<", "<<init_sigma_abs[1]<<std::endl;
    protocol_trajectories=parser.getSetAsInt("protocol_trajectories", protocol_trajectories);
    protocol_timestep_count=parser.getSetAsInt("protocol_timestep_count", protocol_timestep_count);;

    heatup_steps=parser.getSetAsInt("heatup_steps", heatup_steps);

    store_walker_statistics=parser.getSetAsBool("store_walker_statistics", store_walker_statistics);;
    walker_statistics_averageduration=parser.getSetAsDouble("walker_statistics_averageduration", walker_statistics_averageduration);
    walker_statistics.clear();
    walker_statistics.push_back(walker_statistics_entry());

}

void FluorophorDynamics::read_config(jkINIParser2& parser, std::string group, std::string supergroup) {
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

        parser.leaveGroup();
    }

    std::cout<<object_name<<" RNG-TEST: "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<"\n";

    //init();
}

void FluorophorDynamics::change_walker_count(unsigned long N_walker, unsigned long N_fluorophores) {
    std::cout<<object_name<<": change_walker_count ... N_walker="<<N_walker<<"   N_fluorophores="<<N_fluorophores<<"\n";
    if (walker_state!=NULL) {
        free(walker_state);
    }
    if (use_two_walkerstates && walker_state_other!=NULL) {
        free(walker_state_other);
    }
    if (walker_dx) free(walker_dx);
    if (walker_dy) free(walker_dy);
    if (walker_dz) free(walker_dz);
    walker_dx=NULL;
    walker_dy=NULL;
    walker_dz=NULL;
    walker_state=NULL;
    walker_state_other=NULL;
    if (N_walker>0) {
        walker_state=(walkerState*)calloc(N_walker*N_fluorophores, sizeof(walkerState));
        walker_dx=(double*)calloc(N_walker*N_fluorophores, sizeof(double));
        walker_dy=(double*)calloc(N_walker*N_fluorophores, sizeof(double));
        walker_dz=(double*)calloc(N_walker*N_fluorophores, sizeof(double));
        if (use_two_walkerstates) {
            walker_state_other=(walkerState*)calloc(N_walker*N_fluorophores, sizeof(walkerState));
        } else {
            walker_state_other=walker_state;
        }
        if (walker_state==NULL || walker_dx==NULL ||  walker_dy==NULL ||  walker_dz==NULL ) {
            throw FluorophorException(format("could not reserve enough memory walker_state=0x%X N_walker=%lu N_fluorophores=%lu sizeof(walkerState)=%lu", walker_state, N_walker, N_fluorophores, sizeof(walkerState)));
        }

    }
    walker_count=N_walker;
    n_fluorophores=N_fluorophores;
    init_walkers();

    if (notify_when_walkercount_changes.size()>0) {
        for (size_t i=0; i<notify_when_walkercount_changes.size(); i++) {
                notify_when_walkercount_changes[i]->handle_parent_walker_count_changed(walker_count, n_fluorophores);
        }
    }
}

void FluorophorDynamics::init_walkers() {
    std::cout<<object_name<<": initing all walkers ... walker_count="<<walker_count<<"   n_fluorophores="<<n_fluorophores<<"\n";
    for (unsigned long i=0; i<walker_count*n_fluorophores; i++) {
        init_walker(i);
    }
    if (n_fluorophores>1) {
        if (additional_walker_position_mode==SamePosition) {
            std::cout<<object_name<<": initializing additional walkser in same position\n";
            for (unsigned long i=0; i<walker_count*n_fluorophores; i++) {
                walker_dx[i]=0;
                walker_dy[i]=0;
                walker_dz[i]=0;
            }
        } else if (additional_walker_position_mode==InSphere) {
            std::cout<<object_name<<": initializing additional walkser in sphere: walker_count="<<walker_count<<" n_fluorophores="<<n_fluorophores<<" additional_walker_sphere_radius="<<additional_walker_sphere_radius<<"\n";
            for (unsigned long i=0; i<walker_count*n_fluorophores; i++) {
                double r2=0;
                do {
                    walker_dx[i]=gsl_ran_flat(rng, -additional_walker_sphere_radius, additional_walker_sphere_radius);
                    walker_dy[i]=gsl_ran_flat(rng, -additional_walker_sphere_radius, additional_walker_sphere_radius);
                    walker_dz[i]=gsl_ran_flat(rng, -additional_walker_sphere_radius, additional_walker_sphere_radius);
                    r2 = walker_dx[i]*walker_dx[i] + walker_dy[i]*walker_dy[i] + walker_dz[i]*walker_dz[i];
                } while (r2>additional_walker_sphere_radius*additional_walker_sphere_radius);
            }
        }
    }
    init_additional_walkers();
}

void FluorophorDynamics::init_additional_walkers() {
    std::cout<<object_name<<": init_additional_walkers()   n_fluorophores="<<n_fluorophores<<"  walker_count="<<walker_count<<"\n";
    if (n_fluorophores>1 && walker_count>0) {
        for (unsigned long i=0; i<walker_count; i++) {
            for (unsigned long j=1; j<n_fluorophores; j++) {
                init_walker(i+j*walker_count, walker_state[i].x, walker_state[i].y, walker_state[i].z);
                //std::cout<<object_name<<": w["<<i+j*walker_count<<"]: qeff="<<walker_state[i+j*walker_count].q_fluor[walker_state[i+j*walker_count].qm_state]<<"  state="<<walker_state[i+j*walker_count].qm_state<<"\n";
            }
        }
        propagate_additional_walkers();
    }
}

void FluorophorDynamics::propagate_additional_walkers() {
    if (n_fluorophores>1 && walker_count>0) {
        if (additional_walker_position_mode==SamePosition) {
            for (long i=0; i<walker_count; i++) {
                for (long j=1; j<n_fluorophores; j++) {
                    walker_state[i+j*walker_count].time=walker_state[i].time;
                    walker_state[i+j*walker_count].p_z=walker_state[i].p_z;
                    walker_state[i+j*walker_count].p_x=walker_state[i].p_x;
                    walker_state[i+j*walker_count].p_y=walker_state[i].p_y;
                    walker_state[i+j*walker_count].ix=walker_state[i].ix;
                    walker_state[i+j*walker_count].iy=walker_state[i].iy;
                    walker_state[i+j*walker_count].iz=walker_state[i].iz;
                    walker_state[i+j*walker_count].x=walker_state[i].x;
                    walker_state[i+j*walker_count].y=walker_state[i].y;
                    walker_state[i+j*walker_count].z=walker_state[i].z;
                    walker_state[i+j*walker_count].was_just_reset=walker_state[i].was_just_reset;
                    if (additional_walker_off_if_main_off && !walker_state[i].exists) walker_state[i+j*walker_count].exists=false;
                }
            }
        } else if (additional_walker_position_mode==InSphere) {
            for (long i=0; i<walker_count; i++) {
                for (long j=1; j<n_fluorophores; j++) {
                    walker_state[i+j*walker_count].time=walker_state[i].time;
                    walker_state[i+j*walker_count].p_z=walker_state[i].p_z;
                    walker_state[i+j*walker_count].p_x=walker_state[i].p_x;
                    walker_state[i+j*walker_count].p_y=walker_state[i].p_y;
                    walker_state[i+j*walker_count].ix=walker_state[i].ix;
                    walker_state[i+j*walker_count].iy=walker_state[i].iy;
                    walker_state[i+j*walker_count].iz=walker_state[i].iz;
                    walker_state[i+j*walker_count].x=walker_state[i].x+walker_dx[i+j*walker_count];
                    walker_state[i+j*walker_count].y=walker_state[i].y+walker_dy[i+j*walker_count];
                    walker_state[i+j*walker_count].z=walker_state[i].z+walker_dz[i+j*walker_count];
                    walker_state[i+j*walker_count].was_just_reset=walker_state[i].was_just_reset;
                    if (additional_walker_off_if_main_off && !walker_state[i].exists) walker_state[i+j*walker_count].exists=false;
                }
            }
        }

        if (additional_walker_photophysics && use_photophysics) {
            for (long i=0; i<walker_count; i++) {
                for (long j=1; j<n_fluorophores; j++) {
                    propagate_photophysics(i+j*walker_count);
                }
            }
        }
        /*for (long i=0; i<3; i++) {
            std::cout<<i<<": x="<<walker_state[i].x<<"   "<<object_name;
            for (long j=1; j<n_fluorophores; j++) {
                std::cout<<"\n   dx"<<j<<"="<<walker_state[i+j*walker_count].x-walker_state[i].x<<" q"<<j<<"="<<walker_state[i+j*walker_count].q_fluor[walker_state[i+j*walker_count].qm_state]<<" s"<<j<<"="<<walker_state[i+j*walker_count].qm_state;
            }
            std::cout<<"\n";
        }*/
    }
}

void FluorophorDynamics::init_walker(unsigned long i, double x, double y, double z) {
    walker_state[i].time=0;
    walker_state[i].exists=true;
    walker_state[i].was_just_reset=false;
    walker_state[i].x=x;
    walker_state[i].y=y;
    walker_state[i].z=z;
    walker_state[i].ix=x;
    walker_state[i].iy=y;
    walker_state[i].iz=z;
    walker_state[i].x0=x;
    walker_state[i].y0=y;
    walker_state[i].z0=z;
    walker_state[i].ix0=x;
    walker_state[i].iy0=y;
    walker_state[i].iz0=z;
    walker_state[i].p_x=init_p_x;
    walker_state[i].p_y=init_p_y;
    walker_state[i].p_z=init_p_z;
    for (int j=0; j<N_FLUORESCENT_STATES; j++) walker_state[i].sigma_abs[j]=init_sigma_abs[j];
    for (int j=0; j<N_FLUORESCENT_STATES; j++) walker_state[i].q_fluor[j]=init_q_fluor[j];
    for (int j=0; j<N_FLUORESCENT_STATES*N_FLUORESCENT_STATES; j++) walker_state[i].photophysics_transition[j]=init_photophysics_transition[j];
    walker_state[i].qm_state=init_qm_state;
    walker_state[i].used_qm_states=init_used_qm_states;
    walker_state[i].type=init_type;
    walker_state[i].spectrum=init_spectrum;
    //walker_state[i].user_data=NULL;
    fluorophors->load_spectrum(init_spectrum);

    //std::cout<<"init_walker: init_spectrum="<<init_spectrum<<" init_fluorophor="<<init_fluorophor<<std::endl;
}

void FluorophorDynamics::init(){

    std::cout<<object_name<<" RNG-TEST: "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<", "<<gsl_rng_get(this->rng)<<"\n";

    sim_time=0;
    endoftrajectory=false;

    // ensure that all calculated variables are set properly.
    // set_sim_timestep also allocates memory for the walker states
    set_c_fluor(c_fluor);
    set_sim_timestep(sim_timestep);

    // now we initialize the walkers at random possitions inside the simulation volume
    if (volume_shape==Box) {
        for (unsigned long i=0; i<walker_count; i++) {
            double x=sim_x*gsl_rng_uniform(rng);
            double y=sim_y*gsl_rng_uniform(rng);
            double z=sim_z*gsl_rng_uniform(rng);
            init_walker(i, x, y, z);
        }
    } else if (volume_shape==Ball) {
        for (unsigned long i=0; i<walker_count; i++) {
            register double x,y,z;
            do {
                x=sim_x*gsl_rng_uniform(rng);
                y=sim_y*gsl_rng_uniform(rng);
                z=sim_z*gsl_rng_uniform(rng);
            } while (gsl_pow_2(x)+gsl_pow_2(y)+gsl_pow_2(z)>gsl_pow_2(sim_radius));
            init_walker(i,x,y,z);
        }
    }
    init_additional_walkers();

    //sim_time=0;
    char fn[1024];
    std::cout<<"protocol_trajectories="<<protocol_trajectories<<std::endl;
    if (protocol_trajectories>0) {
        trajectoryFile=(FILE**)malloc(sizeof(FILE*)*protocol_trajectories);
        for (unsigned int i=0; i<protocol_trajectories; i++) {
            sprintf(fn, "%s%straj%03d.plt", basename.c_str(), object_name.c_str(), i);
            std::cout<<"writing '"<<fn<<"' ... ";
            FILE* f=fopen(fn, "w");
            sprintf(fn, "%s%straj%03d.dat", basename.c_str(), object_name.c_str(), i);
            fprintf(f, "unset multiplot\n");
            fprintf(f, "reset\n");
            fprintf(f, "set title \"3D Translation Trajectory\"\n");
            fprintf(f, "set xlabel \"x\"\n");
            fprintf(f, "set ylabel \"y\"\n");
            fprintf(f, "set zlabel \"z\"\n");
            fprintf(f, "splot \"%s\" using 2:3:4 notitle with lines\n", extract_file_name(fn).c_str());
            fprintf(f, "pause -1\n");
            fprintf(f, "set title \"Rotation Trajectory on r=1 sphere\"\n");
            fprintf(f, "set dummy u,v\n");
            fprintf(f, "set xlabel \"x\"\n");
            fprintf(f, "set ylabel \"y\"\n");
            fprintf(f, "set zlabel \"z\"\n");
            fprintf(f, "set angles degrees\n");
            fprintf(f, "set parametric\n");
            fprintf(f, "set view 60, 136, 1.22, 1.26\n");
            fprintf(f, "set urange [ -90.0000 : 90.0000 ] noreverse nowriteback\n");
            fprintf(f, "set vrange [ 0.00000 : 360.000 ] noreverse nowriteback\n");
            fprintf(f, "splot cos(u)*cos(v),cos(u)*sin(v),sin(u) notitle with lines lt 13,\\\n");
            fprintf(f, "      \"%s\" using 11:12:13 notitle with lines\n", extract_file_name(fn).c_str());
            fprintf(f, "pause -1\n");
            fprintf(f, "unset parametric\n");
            fclose(f);

            sprintf(fn, "%s%straj%03d.dat", basename.c_str(), object_name.c_str(), i);
            std::cout<<"done!\nopening '"<<fn<<"' ... ";
            trajectoryFile[i]=fopen(fn, "w");
            fprintf(trajectoryFile[i], "# t [seconds], x [microns], y [microns], z [microns], qm_state, type, spectrum, sigma_abs [meters^2], q_fluor [0..1], tau_fl [secs], p_x, p_y, p_z\n");
            std::cout<<"done!\n";
        }
    }

    walker_statistics_nextsavetime=sim_time+walker_statistics_averageduration;
    // create new step protocol
    walker_statistics.push_back(walker_statistics_entry());
    walker_statistics[walker_statistics.size()-1].time=sim_time+walker_statistics_averageduration/2.0;

    store_step_protocol();

    if (heatup_steps>0) {
        std::cout<<"heating up "<<object_name<<" ("<<heatup_steps<<" steps)... ";
        for (int i=0; i<heatup_steps; i++) {
            propagate();
            propagate_additional_walkers();
        }
        sim_time=0;
        std::cout<<" DONE!\n";
    }

}

void FluorophorDynamics::propagate(bool boundary_check){


    sim_time=sim_time+sim_timestep;

    //std::cout<<">>> dyn    sim_time = "<<sim_time<<"\n";
}

RelativeAbsorbanceReader* FluorophorDynamics::get_abs_reader() const {
    if (measmap.count(absorbance_reader)>0) {
        FluorescenceMeasurement* fm=measmap[absorbance_reader];
        return dynamic_cast<RelativeAbsorbanceReader*>(fm);
    } else {
        return NULL;
    }
}
void FluorophorDynamics::propagate_photophysics(int i) {
    //    abs_reader=NULL;
    //photophysics_absorbance_dependent=false;
    //photophysics_absorbance_factor=1;
    register double alpha=1.0;
    if (photophysics_absorbance_dependent) {
        RelativeAbsorbanceReader* absorbance_reader=get_abs_reader();
        if(absorbance_reader) {
            alpha=absorbance_reader->get_relative_absorbance_for(this, i)*photophysics_absorbance_factor;
        }
    }
    propagate_photophysics_scaled(i, alpha);
}

void FluorophorDynamics::propagate_photophysics_scaled(int i, double scale) {
    if (!use_photophysics) return;
    if (walker_state[i].used_qm_states<=1) return;
    const int& state=walker_state[i].qm_state;
    register double p=0;
    const double r=gsl_rng_uniform(rng);
    register int maxs=mmin(walker_state[i].used_qm_states,N_FLUORESCENT_STATES);
    // omit the diagonal-element
    register int f=0;
    while (f<maxs) {
        if (f!=state) {
            p=p+walker_state[i].photophysics_transition[state*N_FLUORESCENT_STATES+f]*scale;
            if (r <= p) {
                walker_state[i].qm_state=f;
                return;
            }
        }
        f++;
    }

    walker_state[i].qm_state=state;
}

void FluorophorDynamics::test_photophysics(unsigned int steps, unsigned int walkers) {
    std::cout<<"init photophysics-test ... \n";
    unsigned int nw=walkers;
    unsigned int plot_walkers=5;
    init();
    change_walker_count(nw, 1);
    std::cout<<"setup photophysics-test ... \n";
    int used_states=0;
    for (unsigned int i=0; i<nw; i++) {
        walker_state[i].x=0;
        walker_state[i].y=0;
        walker_state[i].z=0;
        walker_state[i].p_x=0;
        walker_state[i].p_y=0;
        walker_state[i].p_z=1;
        walker_state[i].qm_state=0;
        used_states=mmax(used_states, walker_state[i].used_qm_states);
    }

    //sim_time=0;
    std::cout<<"opening output file ... basename="<<basename<<"  object_name="<<object_name<<"\n";
    FILE* f1=fopen((basename+object_name+"photophysicstest_traj.dat").c_str(), "w");
    int nhist=100;
    double state_count[N_FLUORESCENT_STATES];
    for (int i=0; i<N_FLUORESCENT_STATES; i++) {
        state_count[i]=0;
    }

    std::cout<<"starting photophysics-test propagation ";
    for (unsigned int i=0; i<steps; i++) {
        for (unsigned int w=0; w<nw; w++) {
            propagate_photophysics_scaled(w, 1);
        }
        //propagate(true);
        for (unsigned int w=0; w<nw; w++) {
            state_count[walker_state[w].qm_state]++;
        }
        fprintf(f1, "%20.10lf", (double)(i+1)*sim_timestep*1e6);
        for (unsigned int w=0; w<mmin(plot_walkers,nw); w++) {
            fprintf(f1, ", %d", walker_state[w].qm_state);
        }
        fprintf(f1, "\n");
        if (i%100==0) std::cout<<".";
    }
    fclose(f1);
    std::cout<<" ready!"<<std::endl;
    std::cout<<"writing output-files ..."<<std::endl;

    used_states=0;
    for (int i=0; i<N_FLUORESCENT_STATES; i++) {
        state_count[i]=state_count[i]/double(steps*nw);
        if (state_count[i]>0.0) used_states=i+1;
    }
    std::cout<<"   - histogram"<<std::endl;
    FILE* f2=fopen((basename+object_name+"photophysicstest_hist.dat").c_str(), "w");
    for (int i=0; i<N_FLUORESCENT_STATES; i++) {
        fprintf(f2, "%d, %20.10lf\n", i, state_count[i]);
    }
    fclose(f2);
    std::cout<<"   - plot"<<std::endl;
    FILE* f=fopen((basename+object_name+"photophysicstest_plot.plt").c_str(), "w");
    fprintf(f, "unset multiplot\n");
    fprintf(f, "reset\n");
    for (int plt=0; plt<2; plt++) {
        std::cout<<"     * plot plt="<<plt<<std::endl;
        if (plt==0) {
            fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
            fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"photophysicstest_plot.pdf").c_str());
        } else if (plt==1) {
            fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
            fprintf(f, "set output\n");
        }
        fprintf(f, "set multiplot layout 3,1\n");
        fprintf(f, "set title \"State Trajectory\"\n");
        fprintf(f, "set xlabel \"simulation time [microseconds]\"\n");
        fprintf(f, "set ylabel \"QM state\"\n");
        fprintf(f, "plot ");
        unsigned int pw=mmin(plot_walkers,nw);
        std::cout<<"     * plot-walkers="<<pw<<std::endl;
        for (unsigned int i=0; i<pw; i++) {
            if (i>0) fprintf(f, ", ");
            fprintf(f, "\"%s\" using 1:(%d+($%d)) title \"walker %d\" with steps", extract_file_name(basename+object_name+"photophysicstest_traj.dat").c_str(), i*(used_states+1), 2+i, i+1);
        }
        fprintf(f, "\n");
        fprintf(f, "plot [0:1000]");
        for (unsigned int i=0; i<pw; i++) {
            if (i>0) fprintf(f, ", ");
            fprintf(f, "\"%s\" using 1:(%d+($%d)) title \"walker %d\" with steps", extract_file_name(basename+object_name+"photophysicstest_traj.dat").c_str(), (used_states+1)*i, 2+i, i+1);
        }
        fprintf(f, "\n");
        fprintf(f, "plot [0:100]");
        for (unsigned int i=0; i<pw; i++) {
            if (i>0) fprintf(f, ", ");
            fprintf(f, "\"%s\" using 1:(%d+($%d)) title \"walker %d\" with steps", extract_file_name(basename+object_name+"photophysicstest_traj.dat").c_str(), (used_states+1)*i, 2+i, i+1);
        }
        fprintf(f, "\n");
        fprintf(f, "unset multiplot\n");
        if (plt==1) fprintf(f, "pause -1\n");
        fprintf(f, "set title \"state histogram\"\n");
        fprintf(f, "set xlabel \"QM state\"\n");
        fprintf(f, "set ylabel \"frequency\"\n");
        fprintf(f, "set boxwidth 0.9 relative\n");
        fprintf(f, "set style fill solid 1.0\n");

        fprintf(f, "plot [0:%d] \"%s\" using 1:2 notitle with boxes\n", used_states, extract_file_name(basename+object_name+"photophysicstest_hist.dat").c_str());
        if (plt==1) fprintf(f, "pause -1\n");
    }
    fclose(f);

    f=fopen((basename+object_name+"photophysicstest_config.txt").c_str(), "w");
    fprintf(f, "%s\n", report().c_str());
    fclose(f);
    std::cout<<"test done!\n";
}

std::string FluorophorDynamics::report() {
    std::string s="";
    s+="object_name = "+object_name+"\n";
    s+="rng = "+std::string(gsl_rng_name(rng))+"\n";
    s+="sim_timestep = "+floattostr(sim_timestep/1e-6)+" microsec = "+floattostr(sim_timestep/1e-9)+" nanosec\n";
    if (volume_shape==0) {
        s+="sim_volume = "+floattostr(sim_x*sim_y*sim_z)+" femto litres\n";
        s+="           = box: x=0.."+floattostr(sim_x)+" microns,     y=0.."+floattostr(sim_y)+" microns,     z=0.."+floattostr(sim_z)+" microns\n";
    } else if (volume_shape==1) {
        s+="sim_volume = "+floattostr(4.0*M_PI/3.0*gsl_pow_3(sim_radius))+" femto litres\n";
        s+="           = sphere, r="+floattostr(sim_radius)+" microns\n";
    }
    s+="walker_count = "+inttostr(walker_count)+"\n";
    s+="visible_walker_count = "+inttostr(get_visible_walker_count())+"\n";
    s+="c_fluor = "+floattostr(c_fluor)+" nMolar\n";
    s+="n_fluorophores = "+inttostr(n_fluorophores)+"\n";
    s+="additional_walker_photophysics = "+booltostr(additional_walker_photophysics)+"\n";
    s+="additional_walker_off_if_main_off = "+booltostr(additional_walker_off_if_main_off)+"\n";
    if (additional_walker_position_mode==InSphere) {
        s+="additional_walker_mode = in_sphere\n";
        s+="additional_walker_sphere_radius = "+floattostr(additional_walker_sphere_radius)+" micron\n";
    } else if (additional_walker_position_mode==SamePosition) {
        s+="additional_walker_mode = same_position\n";
    }

    s+="init_dipole_vec = ("+floattostr(init_p_x)+", "+floattostr(init_p_y)+", "+floattostr(init_p_z)+")\n";
    s+="init_spectrum = "+inttostr(init_spectrum);
    if (init_spectrum==-1) s+=" [none]\n";
    else s+=" ["+extract_file_name(fluorophors->getSpectrumFilename(init_spectrum))+"]\n";
    s+="init_q_fluor = "+doublevectortostr(init_q_fluor, N_FLUORESCENT_STATES)+" \n";
    s+="init_qm_state = "+inttostr(init_qm_state)+"\n";
    s+="init_type = "+inttostr(init_type)+"\n";
    double init_sigma_absnm2[N_FLUORESCENT_STATES];
    for (int i=0; i<N_FLUORESCENT_STATES; i++) init_sigma_absnm2[i]=init_sigma_abs[i]*1e18;
    s+="init_sigma_abs = "+doublevectortostr(init_sigma_abs, N_FLUORESCENT_STATES)+" meters^2\n";
    s+="               = "+doublevectortostr(init_sigma_absnm2, N_FLUORESCENT_STATES)+" nanometers^2\n";
    if (!use_photophysics) {
        s+="no photophysics simulation!\n";
    } else {
        s+="simulation using photophysics with these constants:\n";
        s+="init_used_qm_states = "+inttostr(init_used_qm_states)+"\n";
        for (int i=0; i<N_FLUORESCENT_STATES; i++) {
            s+="P["+inttostr(i)+"->f] = "+doublevectortostr(&init_photophysics_transition[i*N_FLUORESCENT_STATES], N_FLUORESCENT_STATES)+"\n";
        }
        s+="\n";
        //s+=doublearraytostr(init_photophysics_transition, N_FLUORESCENT_STATES,N_FLUORESCENT_STATES,false)+"\n";
        s+="reset_qmstate_at_simboxborder = "+booltostr(reset_qmstate_at_simboxborder)+"\n";
        if (photophysics_absorbance_dependent) {
            s+="photophysics_absorbance_factor = "+floattostr(photophysics_absorbance_factor)+"\n\n";
        }
    }
    s+="depletion_propability = "+floattostr(depletion_propability)+"\n\n";
    s+="heatup_steps = "+inttostr(heatup_steps)+"\n\n";

    s+="store_walker_statistics = "+booltostr(store_walker_statistics)+"\n";
    if (store_walker_statistics) {
        s+="walker_statistics_averageduration = "+floattostr(walker_statistics_averageduration)+" seconds\n";
        s+="walker_statistics_saved_steps = "+inttostr(walker_statistics.size())+"\n";
    }
    s+="\n";

    return s;
}

void FluorophorDynamics::save_trajectories() {
    //FILE* f;
    //char fn[1024];
    if (protocol_trajectories>0) {
        for (unsigned int i=0; i<protocol_trajectories; i++) {
            fclose(trajectoryFile[i]);
        }
        free(trajectoryFile);
    }
};

void FluorophorDynamics::store_step_protocol() {
    //std::cout<<"store step protocol\n";
    for (unsigned int i=0; i<protocol_trajectories; i++) {
        walkerState ws=walker_state[i];
        if (ws.time<protocol_timestep_count || protocol_timestep_count<0)
            fprintf(trajectoryFile[i], "%lg, %lg, %lg, %lg, %d, %d, %d, %lg, %lg, %lg\n", ws.time*sim_timestep, ws.x, ws.y, ws.z, ws.qm_state, ws.type, ws.spectrum, ws.p_x, ws.p_y, ws.p_z);
    }

    if (store_walker_statistics) {
        if (sim_time>=walker_statistics_nextsavetime) {

            // normalize current step
            uint64_t cur=walker_statistics.size()-1;
            walker_statistics[cur].average_brightness=walker_statistics[cur].average_brightness/double(walker_statistics[cur].average_steps);
            walker_statistics[cur].count_all=walker_statistics[cur].count_all/double(walker_statistics[cur].average_steps);
            walker_statistics[cur].count_existing=walker_statistics[cur].count_existing/double(walker_statistics[cur].average_steps);
            walker_statistics[cur].posx=walker_statistics[cur].posx/double(walker_statistics[cur].average_steps);
            walker_statistics[cur].posy=walker_statistics[cur].posy/double(walker_statistics[cur].average_steps);
            walker_statistics[cur].posz=walker_statistics[cur].posz/double(walker_statistics[cur].average_steps);
            for (int i=0; i<N_FLUORESCENT_STATES; i++) walker_statistics[cur].state_distribution[i]=walker_statistics[cur].state_distribution[i]/double(walker_statistics[cur].average_steps);


            // create new step protocol
            walker_statistics.push_back(walker_statistics_entry());
            walker_statistics[walker_statistics.size()-1].time=sim_time+walker_statistics_averageduration/2.0;
            walker_statistics_nextsavetime=walker_statistics_nextsavetime+walker_statistics_averageduration;
        } else {
            uint64_t cur=walker_statistics.size()-1;

            unsigned int wc=get_visible_walker_count();
            walkerState* ws=get_visible_walker_state();
            uint64_t viswalk=0;
            double average_brightness=0;
            double px=0, py=0, pz=0;
            double state_distribution[N_FLUORESCENT_STATES];
            for (int i=0; i<N_FLUORESCENT_STATES; i++) state_distribution[i]=0;
            for (unsigned int i=0; i<wc; i++) {
                if (ws[i].exists) {
                    viswalk++;
                    average_brightness=average_brightness+ws[i].sigma_abs[ws[i].qm_state]*ws[i].q_fluor[ws[i].qm_state];
                    px=px+ws[i].x;
                    py=py+ws[i].y;
                    pz=pz+ws[i].z;
                    state_distribution[ws[i].qm_state]++;
                }
            }
            walker_statistics[cur].average_steps++;
            walker_statistics[cur].average_brightness=walker_statistics[cur].average_brightness+average_brightness/double(viswalk);
            walker_statistics[cur].posx=walker_statistics[cur].posx+px/double(viswalk);
            walker_statistics[cur].posy=walker_statistics[cur].posy+py/double(viswalk);
            walker_statistics[cur].posz=walker_statistics[cur].posz+pz/double(viswalk);
            walker_statistics[cur].count_all=walker_statistics[cur].count_all+wc;
            walker_statistics[cur].count_existing=walker_statistics[cur].count_existing+viswalk;
            for (int i=0; i<N_FLUORESCENT_STATES; i++) walker_statistics[cur].state_distribution[i]=walker_statistics[cur].state_distribution[i]+state_distribution[i]/double(viswalk);
        }
    }
    //std::cout<<"store step protocol done\n";
}

/*void FluorophorDynamics::store_step_protocol(int w_start, int w_end) {
    for (unsigned int i=0; i<protocol_trajectories; i++) {
        if (i>=w_start && i<=w_end) {
            walkerState ws=walker_state[i];
            if (ws.time<protocol_timestep_count || protocol_timestep_count<0)
                fprintf(trajectoryFile[i], "%lg, %lg, %lg, %lg, %d, %d, %d, %lg, %lg, %lg, %lg, %lg, %lg\n", ws.time*sim_timestep, ws.x, ws.y, ws.z, ws.qm_state, ws.type, ws.spectrum, ws.sigma_abs, ws.q_fluor, ws.tau_fl, ws.p_x, ws.p_y, ws.p_z);
        }
    }
}*/




void FluorophorDynamics::load_all_used_spectra() {
    if (walker_count>0) for (unsigned long i=0; i<walker_count; i++) {
        fluorophors->load_spectrum(walker_state[i].spectrum);
    }
}





void FluorophorDynamics::perform_boundary_check(unsigned long i) {
    register double nx=walker_state[i].x;
    register double ny=walker_state[i].y;
    register double nz=walker_state[i].z;
    walker_state[i].was_just_reset=false;
    if (walker_state[i].exists) {
        if (volume_shape==0) {
            if (   (nx<0) || (nx>sim_x)
                || (ny<0) || (ny>sim_y)
                || (nz<0) || (nz>sim_z) ) {

                if (depletion_propability<=0 || gsl_ran_flat(rng,0,1)>depletion_propability) {
                    //std::cout<<"initializing new walker ... depletion_propability="<<depletion_propability<<"\n";


                    // first choose one face of the simulation volume and then set the walker
                    // to any position on the face ... also shift a bit inwards
                    char face=gsl_rng_uniform_int(rng, 6)+1;
                    double x=0,y=0,z=0;
                    switch(face) {
                        case 1:
                            //x-y-plane at z=0
                            x=gsl_ran_flat(rng, 0, sim_x);
                            y=gsl_ran_flat(rng, 0, sim_y);
                            z=0;
                            break;
                        case 2:
                            //x-y-plane at z=sim_z
                            x=gsl_ran_flat(rng, 0, sim_x);
                            y=gsl_ran_flat(rng, 0, sim_y);
                            z=sim_z;
                            break;
                        case 3:
                            //x-z-plane at y=0
                            x=gsl_ran_flat(rng, 0, sim_x);
                            y=0;
                            z=gsl_ran_flat(rng, 0, sim_z);
                            break;
                        case 4:
                            //x-z-plane at y=sim_y
                            x=gsl_ran_flat(rng, 0, sim_x);
                            y=sim_y;
                            z=gsl_ran_flat(rng, 0, sim_z);
                            break;
                        case 5:
                            //z-y-plane at x=0
                            x=0;
                            y=gsl_ran_flat(rng, 0, sim_y);
                            z=gsl_ran_flat(rng, 0, sim_z);
                            break;
                        case 6:
                            //z-y-plane at x=sim_x
                            x=sim_x;
                            y=gsl_ran_flat(rng, 0, sim_y);
                            z=gsl_ran_flat(rng, 0, sim_z);
                            break;
                    }
                    if (reset_qmstate_at_simboxborder)  {
                        init_walker(i, x, y, z);
                    }   else {
                        walker_state[i].x=x;
                        walker_state[i].y=y;
                        walker_state[i].z=z;
                    }
                    walker_state[i].was_just_reset=true;
                    walker_state[i].time=0;
                    walker_state[i].x0=walker_state[i].x;
                    walker_state[i].y0=walker_state[i].y;
                    walker_state[i].z0=walker_state[i].z;
                } else {
                    //std::cout<<"killing walker ... depletion_propability="<<depletion_propability<<"\n";
                    walker_state[i].exists=false;
                }
            }
        } else if (volume_shape==1) {
            if (gsl_pow_2(nx)+gsl_pow_2(ny)+gsl_pow_2(nz)>gsl_pow_2(sim_radius)) {
                if (depletion_propability<=0 || gsl_ran_flat(rng,0,1)>depletion_propability) {
                    //std::cout<<"initializing new walker ... depletion_propability="<<depletion_propability<<"\n";
                    gsl_ran_dir_3d(rng, &nx, &ny, &nz);
                    double x=sim_radius*nx;
                    double y=sim_radius*ny;
                    double z=sim_radius*nz;
                    if (reset_qmstate_at_simboxborder)  {
                        init_walker(i, x, y, z);
                    }   else {
                        walker_state[i].x=x;
                        walker_state[i].y=y;
                        walker_state[i].z=z;

                    }
                    walker_state[i].time=0;
                    walker_state[i].was_just_reset=true;
                    walker_state[i].x0=walker_state[i].x;
                    walker_state[i].y0=walker_state[i].y;
                    walker_state[i].z0=walker_state[i].z;
                } else {
                    //std::cout<<"killing walker ... depletion_propability="<<depletion_propability<<"\n";
                    walker_state[i].exists=false;
                }
            }
        }
    }
}

void FluorophorDynamics::set_sim_timestep(double value) {
    sim_timestep=value;
};

void FluorophorDynamics::set_c_fluor(double value) {
    c_fluor=value;
    change_walker_count(calc_walker_count(), n_fluorophores);

};

void FluorophorDynamics::set_sim_box(double vx, double vy, double vz) {
    volume_shape=Box;
    sim_x=vx;
    sim_y=vy;
    sim_z=vz;
    change_walker_count(calc_walker_count(), n_fluorophores);
};

void FluorophorDynamics::set_sim_sphere(double rad) {
    volume_shape=Ball;
    sim_radius=rad;
    unsigned long wc=calc_walker_count();
    //std::cout<<object_name<<": set_sim_sphere(rad="<<rad<<"):  walker_count="<<wc<<"  n_fluorophores="<<n_fluorophores<<std::endl;
    change_walker_count(wc, n_fluorophores);
};

void FluorophorDynamics::finalize_sim() {
};

double FluorophorDynamics::estimate_runtime() {
    return 0;
}

void FluorophorDynamics::save() {
}

bool FluorophorDynamics::depends_on(const FluorophorDynamics* other) const {
    return false;
}

void FluorophorDynamics::save_results() {
    save();

    FILE* f;
    char fn[1024];

    if (additional_walker_position_mode==InSphere) {
        sprintf(fn, "%s%sinspherefluorophores.dat", basename.c_str(), object_name.c_str());
        std::cout<<"writing '"<<fn<<"' ...";
        f=fopen(fn, "w");

        int density_width=100;
        uint64_t* countxy=(uint64_t*)calloc(density_width*density_width,sizeof(uint64_t));
        uint64_t* countxz=(uint64_t*)calloc(density_width*density_width,sizeof(uint64_t));
        uint64_t* countyz=(uint64_t*)calloc(density_width*density_width,sizeof(uint64_t));
        uint64_t countN=0;
        uint64_t countNout=0;

        for (long i=0; i<walker_count; i++) {
            for (long j=1; j<n_fluorophores; j++) {
                int x=(walker_dx[i+j*walker_count]+additional_walker_sphere_radius)/(2.0*additional_walker_sphere_radius)*density_width;
                int y=(walker_dy[i+j*walker_count]+additional_walker_sphere_radius)/(2.0*additional_walker_sphere_radius)*density_width;
                int z=(walker_dz[i+j*walker_count]+additional_walker_sphere_radius)/(2.0*additional_walker_sphere_radius)*density_width;
                if (x>=0 && y>=0 && z>=0 && x<density_width && y<density_width && z<density_width) {
                    countxy[y*density_width+x]++;
                    countxz[z*density_width+x]++;
                    countyz[z*density_width+y]++;
                    countN++;
                } else {
                    countNout++;
                }
            }
        }
        for (int y=0; y<density_width; y++) {
            for (int x=0; x<density_width; x++) {
                fprintf(f, "%s ", inttostr(countxy[y*density_width+x]).c_str());
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n\n");
        for (int y=0; y<density_width; y++) {
            for (int x=0; x<density_width; x++) {
                fprintf(f, "%s ", inttostr(countxz[y*density_width+x]).c_str());
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n\n");
        for (int y=0; y<density_width; y++) {
            for (int x=0; x<density_width; x++) {
                fprintf(f, "%s ", inttostr(countyz[y*density_width+x]).c_str());
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n\n");
        int example_walker=mmin(walker_count, 5);
        for (long i=0; i<mmin(walker_count, example_walker); i++) {
            fprintf(f, "0 0 0\n");
            for (long j=1; j<n_fluorophores; j++) {
                fprintf(f, "%15.10lf %15.10lf %15.10lf\n", walker_dx[i+j*walker_count], walker_dy[i+j*walker_count], walker_dz[i+j*walker_count]);
            }
            fprintf(f, "\n\n");
        }
        fclose(f);
        std::cout<<" done!\n";
        std::string fnsphere=fn;

        free(countxy);
        free(countyz);
        free(countxz);


        sprintf(fn, "%s%sinspherefluorophores.plt", basename.c_str(), object_name.c_str());
        std::cout<<"writing '"<<fn<<"' ...";
        f=fopen(fn, "w");
        fprintf(f, "sphere_radius=%lf\n", additional_walker_sphere_radius);
        fprintf(f, "density_width=%s\n", inttostr(density_width).c_str());
        fprintf(f, "countN=%s\n", inttostr(countN).c_str());
        fprintf(f, "countNout=%s\n", inttostr(countNout).c_str());
        fprintf(f, "deltax=%lf\n", 2.0*additional_walker_sphere_radius/double(density_width));
        fprintf(f, "deltay=%lf\n", 2.0*additional_walker_sphere_radius/double(density_width));
        fprintf(f, "deltaz=%lf\n", 2.0*additional_walker_sphere_radius/double(density_width));
        for (int plt=0; plt<2; plt++) {
            if (plt==0) {
                fprintf(f, "set terminal pdfcairo color solid font \"%s, 5\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"inspherefluorophores.pdf").c_str());
            } else if (plt==1) {
                fprintf(f, "set terminal wxt font \"%s, 5\"\n", GNUPLOT_FONT);
                fprintf(f, "set output\n");
            }

            fprintf(f, "set size noratio\n");
            fprintf(f, "set multiplot layout 1,3\n");
            fprintf(f, "set xlabel \"position x [micron]\"\n");
            fprintf(f, "set ylabel \"position y [micron]\"\n");
            fprintf(f, "set size ratio deltay/deltax\n");
            fprintf(f, "set title sprintf(\"fluorophore density in sphere: xy-projection: %s  (outside: %%d)\", countNout)\n", object_name.c_str());
            fprintf(f, "plot [0:density_width*deltax] [0:density_width*deltay] \"%s\" using (($1)*deltax):(($2)*deltay):(($3)/countN) matrix index 0 notitle with image"
            "\n", extract_file_name(fnsphere).c_str());
            fprintf(f, "set title sprintf(\"fluorophore density in sphere: xz-projection: %s  (outside: %%d)\", countNout)\n", object_name.c_str());
            fprintf(f, "set xlabel \"position x [micron]\"\n");
            fprintf(f, "set ylabel \"position z [micron]\"\n");
            fprintf(f, "set size ratio deltaz/deltax\n");
            fprintf(f, "plot [0:density_width*deltax] [0:density_width*deltaz] \"%s\" using (($1)*deltax):(($2)*deltaz):(($3)/countN) matrix index 1 notitle with image"
            "\n", extract_file_name(fnsphere).c_str());
            fprintf(f, "set title sprintf(\"fluorophore density in sphere: yz-projection: %s  (outside: %%d)\", countNout)\n", object_name.c_str());
            fprintf(f, "set xlabel \"position y [micron]\"\n");
            fprintf(f, "set ylabel \"position z [micron]\"\n");
            fprintf(f, "set size ratio deltaz/deltay\n");
            fprintf(f, "plot [0:density_width*deltay] [0:density_width*deltaz] \"%s\" using (($1)*deltay):(($2)*deltaz):(($3)/countN) matrix index 2 notitle with image"
            "\n", extract_file_name(fnsphere).c_str());
            fprintf(f, "unset multiplot\n");
            if (plt==1) fprintf(f, "pause -1\n");

            for (int ex=0; ex<example_walker; ex++) {
                fprintf(f, "set size noratio\n");
                fprintf(f, "set multiplot layout 2,2\n");
                fprintf(f, "set xlabel \"position x [micron]\"\n");
                fprintf(f, "set ylabel \"position y [micron]\"\n");
                fprintf(f, "set size ratio 1\n");
                fprintf(f, "set title \"fluorophore in walker %d: xy-projection: %s\"\n", ex, object_name.c_str());
                fprintf(f, "plot [-1.0*sphere_radius:sphere_radius] [-1.0*sphere_radius:sphere_radius] \"%s\" using 1:2 index %d notitle with points"
                "\n", extract_file_name(fnsphere).c_str(), ex+3);
                fprintf(f, "set xlabel \"position x [micron]\"\n");
                fprintf(f, "set ylabel \"position z [micron]\"\n");
                fprintf(f, "set title \"fluorophore in walker %d: xz-projection: %s\"\n", ex, object_name.c_str());
                fprintf(f, "set size ratio 1\n");
                fprintf(f, "plot [-1.0*sphere_radius:sphere_radius] [-1.0*sphere_radius:sphere_radius] \"%s\" using 1:3 index %d notitle with points"
                "\n", extract_file_name(fnsphere).c_str(), ex+3);
                fprintf(f, "set xlabel \"position y [micron]\"\n");
                fprintf(f, "set ylabel \"position z [micron]\"\n");
                fprintf(f, "set title \"fluorophore in walker %d: yz-projection: %s\"\n", ex, object_name.c_str());
                fprintf(f, "set size ratio 1\n");
                fprintf(f, "plot [-1.0*sphere_radius:sphere_radius] [-1.0*sphere_radius:sphere_radius] \"%s\" using 2:3 index %d notitle with points"
                "\n", extract_file_name(fnsphere).c_str(), ex+3);
                fprintf(f, "set xlabel \"position x [micron]\"\n");
                fprintf(f, "set ylabel \"position y [micron]\"\n");
                fprintf(f, "set zlabel \"position z [micron]\"\n");
                fprintf(f, "set title \"fluorophore in walker %d: 3D-projection: %s\"\n", ex, object_name.c_str());
                fprintf(f, "set size ratio 1\n");
                fprintf(f, "splot [-1.0*sphere_radius:sphere_radius] [-1.0*sphere_radius:sphere_radius] [-1.0*sphere_radius:sphere_radius] \"%s\" using 1:2:3 index %d notitle with points"
                "\n", extract_file_name(fnsphere).c_str(), ex+3);
                fprintf(f, "unset multiplot\n");
                if (plt==1) fprintf(f, "pause -1\n");
            }
        }


        fclose(f);
        std::cout<<" done!\n";
    }



    if (store_walker_statistics) {
        sprintf(fn, "%s%swalkerstatistics.dat", basename.c_str(), object_name.c_str());
        std::cout<<"writing '"<<fn<<"' ...";
        f=fopen(fn, "w");

        if (walker_statistics.size()>2) {
            for (uint64_t i=2; i<walker_statistics.size()-2; i++) {
                fprintf(f, "%15.10lf %lld %lld %15.10lf %lld %15.10lf %15.10lf %15.10lf", walker_statistics[i].time,
                                                            (int64_t)walker_statistics[i].count_all,
                                                            (int64_t)walker_statistics[i].count_existing,
                                                            walker_statistics[i].average_brightness,
                                                            (int64_t)walker_statistics[i].average_steps,
                                                            walker_statistics[i].posx,
                                                            walker_statistics[i].posy,
                                                            walker_statistics[i].posz
                                                        );
                for (long j=0; j<N_FLUORESCENT_STATES; j++) {
                    fprintf(f, " %15.10lf", walker_statistics[i].state_distribution[j]);
                }
                fprintf(f, "\n");
            }
        }
        fprintf(f, "\n");
        fclose(f);
        std::cout<<" done!\n";
        std::string fnwalkerstat=fn;



        sprintf(fn, "%s%swalkerstatistics.plt", basename.c_str(), object_name.c_str());
        std::cout<<"writing '"<<fn<<"' ...";
        f=fopen(fn, "w");
        for (int plt=0; plt<2; plt++) {
            if (plt==0) {
                fprintf(f, "set terminal pdfcairo color solid font \"%s, 5\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"walkerstatistics.pdf").c_str());
            } else if (plt==1) {
                fprintf(f, "set terminal wxt font \"%s, 5\"\n", GNUPLOT_FONT);
                fprintf(f, "set output\n");
            }

            fprintf(f, "set size noratio\n");
            fprintf(f, "set multiplot layout 2,1\n");
            fprintf(f, "set xlabel \"time t [seconds]\"\n");
            fprintf(f, "set ylabel \"average walker count\"\n");
            fprintf(f, "set title \"walker count: %s, average over %lfms\"\n", object_name.c_str(), walker_statistics_averageduration*1000.0);
            fprintf(f, "plot \"%s\" using 1:2 title \"all walkers\" with linespoints,  \"%s\" using 1:3 title \"visible walkers\" with linespoints\n"
            "\n", extract_file_name(fnwalkerstat).c_str(), extract_file_name(fnwalkerstat).c_str());

            fprintf(f, "set xlabel \"time t [seconds]\"\n");
            fprintf(f, "set ylabel \"average brightness (q*sigma)\"\n");
            fprintf(f, "set title \"walker count: %s, average over %lfms\"\n", object_name.c_str(), walker_statistics_averageduration*1000.0);
            fprintf(f, "plot \"%s\" using 1:4 notitle  with linespoints\n"
            "\n", extract_file_name(fnwalkerstat).c_str());

            fprintf(f, "unset multiplot\n");
            if (plt==1) fprintf(f, "pause -1\n\n");

            fprintf(f, "set xlabel \"time t [seconds]\"\n");
            fprintf(f, "set ylabel \"position [micron]\"\n");
            fprintf(f, "set title \"center of mass: %s, average over %lfms\"\n", object_name.c_str(), walker_statistics_averageduration*1000.0);
            fprintf(f, "plot \"%s\" using 1:6 title \"x-position\"  with linespoints, \"%s\" using 1:7 title \"y-position\"  with linespoints, \"%s\" using 1:8 title \"z-position\"  with linespoints\n", extract_file_name(fnwalkerstat).c_str(), extract_file_name(fnwalkerstat).c_str(), extract_file_name(fnwalkerstat).c_str());
            if (plt==1) fprintf(f, "pause -1\n\n");

            fprintf(f, "set xlabel \"time t [seconds]\"\n");
            fprintf(f, "set ylabel \"fraction of particles in QM state\"\n");
            fprintf(f, "set title \"QM state occupation: %s, average over %lfms\"\n", object_name.c_str(), walker_statistics_averageduration*1000.0);
            fprintf(f, "plot [][-0.1:1.1]");
            for (long j=0; j<N_FLUORESCENT_STATES; j++) {
                if (j>0) fprintf(f, ",  \\\n");
                                fprintf(f, "   \"%s\" using 1:%ld title \"qm_state %ld\"  with lines", extract_file_name(fnwalkerstat).c_str(), j+9, j+1);
            }
            fprintf(f, "\n\n");
            if (plt==1) fprintf(f, "pause -1\n\n");

        }


        fclose(f);
        std::cout<<" done!\n";
    }
}

void FluorophorDynamics::handle_parent_walker_count_changed(unsigned long N_walker, unsigned long N_fluorophores) {
}

void FluorophorDynamics::ensure_dynamics_is_hooked(FluorophorDynamics* other) {
    if (count(notify_when_walkercount_changes.begin(), notify_when_walkercount_changes.end(), other)<=0) {
        notify_when_walkercount_changes.push_back(other);
    }
}

unsigned long FluorophorDynamics::get_walker_count() {
     return walker_count;
}

unsigned long FluorophorDynamics::get_visible_walker_count() {
    return walker_count*n_fluorophores;
}


double FluorophorDynamics::get_walker_sigma_times_qfl(unsigned long i) {
    register double s=0;
    register int state=walker_state[i].qm_state;
    if ((state>=0)&&(state<N_FLUORESCENT_STATES)) {
        s=walker_state[i].sigma_abs[state]*walker_state[i].q_fluor[state];
    }
    return s;
}

double FluorophorDynamics::get_visible_walker_sigma_times_qfl(unsigned long i) {
    return get_walker_sigma_times_qfl(i);
}

FluorophorDynamics::walkerState* FluorophorDynamics::get_visible_walker_state() {
    return get_walker_state();
}
