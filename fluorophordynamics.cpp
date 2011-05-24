#include "fluorophordynamics.h"

//#include <boost/thread/mutex.hpp>

#include "../../../LIB/trunk/datatable.h"
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>

FluorophorDynamics::FluorophorDynamics(FluorophorManager* fluorophors, std::string object_name)
{
    this->fluorophors=fluorophors;
    this->object_name=object_name;
    protocol_trajectories=0;
    sim_time=0;
    sim_timestep=1e-6;

    use_two_walkerstates=false;

    endoftrajectory=false;

     // init GSL random number generator
    gsl_rng_env_setup();
    rng_type = gsl_rng_taus2;
    rng = gsl_rng_alloc (rng_type);
    gsl_rng_set(rng, time(0));

    init_p_x=1;
    init_p_y=0;
    init_p_z=0;
    init_qm_state=0;
    init_type=0;
    init_spectrum=-1;
    init_used_qm_states=1;

    for (int j=0; j<N_FLUORESCENT_STATES; j++) init_sigma_abs[j]=2.2e-20;
    for (int j=0; j<N_FLUORESCENT_STATES; j++) init_q_fluor[j]=0.1;
    for (int j=0; j<N_FLUORESCENT_STATES*N_FLUORESCENT_STATES; j++) init_photophysics_transition[j]=0;

    init_spectrum=-1;
    use_photophysics=true;
    protocol_timestep_count=-1;

    walker_state=NULL;
    walker_state_other=NULL;
    this->sim_x=14;
    this->sim_y=14;
    this->sim_z=4;
    this->sim_radius=10;
    this->volume_shape=Box;
    set_c_fluor(50);
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

     // init GSL random number generator
    gsl_rng_env_setup();
    rng_type = gsl_rng_taus2;
    rng = gsl_rng_alloc (rng_type);
    gsl_rng_set(rng, time(0));

    init_p_x=1;
    init_p_y=0;
    init_p_z=0;
    init_qm_state=0;
    init_type=0;
    init_spectrum=-1;
    init_used_qm_states=1;

    for (int j=0; j<N_FLUORESCENT_STATES; j++) init_sigma_abs[j]=2.2e-20;
    for (int j=0; j<N_FLUORESCENT_STATES; j++) init_q_fluor[j]=0.1;
    for (int j=0; j<N_FLUORESCENT_STATES*N_FLUORESCENT_STATES; j++) init_photophysics_transition[j]=0;
    use_photophysics=true;
    protocol_timestep_count=-1;

    walker_state=NULL;
    walker_state_other=NULL;
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

     // init GSL random number generator
    gsl_rng_env_setup();
    rng_type = gsl_rng_taus2;
    rng = gsl_rng_alloc (rng_type);
    gsl_rng_set(rng, time(0));

    init_p_x=1;
    init_p_y=0;
    init_p_z=0;
    init_qm_state=0;
    init_type=0;
    init_spectrum=-1;
    init_used_qm_states=1;

    for (int j=0; j<N_FLUORESCENT_STATES; j++) init_sigma_abs[j]=2.2e-20;
    for (int j=0; j<N_FLUORESCENT_STATES; j++) init_q_fluor[j]=0.1;
    for (int j=0; j<N_FLUORESCENT_STATES*N_FLUORESCENT_STATES; j++) init_photophysics_transition[j]=0;
    use_photophysics=true;
    protocol_timestep_count=-1;

    walker_state=NULL;
    walker_state_other=NULL;
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

void FluorophorDynamics::read_config_internal(jkINIParser2& parser) {
    std::string ivshape="sphere";
    if (volume_shape==Box) ivshape="box";
    std::string vshape=tolower(parser.getSetAsString("volume_shape", ivshape));
    if (vshape=="box") {
        set_sim_box(parser.getSetAsDouble("sim_x", sim_x), parser.getSetAsDouble("sim_y", sim_y), parser.getSetAsDouble("sim_z", sim_z));
    } else if (vshape=="sphere" || vshape=="ball") {
        set_sim_sphere(parser.getSetAsDouble("sim_radius", sim_radius));
    }
    set_c_fluor(parser.getSetAsDouble("c_fluor", c_fluor));


    init_p_x=parser.getSetAsDouble("init_p_x", init_p_x);
    init_p_y=parser.getSetAsDouble("init_p_y", init_p_y);
    init_p_z=parser.getSetAsDouble("init_p_z", init_p_z);
    std::string spec="rho6g";
    if (parser.exists("init_fluorophor")) {
        spec=tolower(parser.getSetAsString("init_fluorophor", spec));
        if (!fluorophors->fluorophorExists(spec)) {
            throw FluorophorException(format("didn't find fluorophor %s in database", spec.c_str()));
        }
        for (int j=0; j<N_FLUORESCENT_STATES; j++) {
            init_sigma_abs[j]=fluorophors->getFluorophorData(spec).sigma_abs;
            init_q_fluor[j]=fluorophors->getFluorophorData(spec).fluorescence_efficiency;
        }
        /*
        init_tau_fl=fluorophors->getFluorophorData(spec).fluorescence_lifetime;
        init_bleaching_propability=fluorophors->getFluorophorData(spec).bleaching_propability;
        init_triplet_lifetime=fluorophors->getFluorophorData(spec).triplet_lifetime;
        init_triplet_propability=fluorophors->getFluorophorData(spec).triplet_propability;
        std::cout<<spec<<": "<<init_q_fluor<<", "<<init_tau_fl<<", "<<init_sigma_abs<<std::endl;*/
    }
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

    if (parser.exists("init_spectrum")) {
        spec=tolower(parser.getSetAsString("init_spectrum", spec));
        //init_spectrum=-1;
        if (fluorophors->getFluorophorData(spec).spectrum!=-1) {
            init_spectrum=fluorophors->getFluorophorData(spec).spectrum;
        } else {
            std::cout<<std::endl<<std::endl<<"didn't find spectrum for "<<spec<<std::endl;
        }
    }
    //std::cout<<spec<<" ("<<init_spectrum<<"): "<<init_q_fluor<<", "<<init_tau_fl<<", "<<init_sigma_abs<<std::endl;
    protocol_trajectories=parser.getSetAsInt("protocol_trajectories", protocol_trajectories);
    protocol_timestep_count=parser.getSetAsInt("protocol_timestep_count", protocol_timestep_count);;
}

void FluorophorDynamics::read_config(jkINIParser2& parser, std::string group, std::string supergroup) {
    basename=parser.getSetAsString("simulation.basename", "");
    std::string rng=tolower(parser.getSetAsString("simulation.rng", "taus2"));
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
    gsl_rng_set(this->rng, time(0));


    for (int i=0; i<2; i++) {
        if (i==0) {
            if (supergroup.size()<=0) break;
            parser.enterGroup(supergroup);
        } else if (i==1) {
            parser.enterGroup(group);
        }

        read_config_internal(parser);

        parser.leaveGroup();
    }

    init();
}

void FluorophorDynamics::change_walker_count(unsigned long N_walker){
    if (walker_state!=NULL) {
        free(walker_state);
    }
    if (use_two_walkerstates && walker_state_other!=NULL) {
        free(walker_state_other);
    }
    walker_state=NULL;
    walker_state_other=NULL;
    if (N_walker>0) {
        walker_state=(walkerState*)calloc(N_walker, sizeof(walkerState));
        if (use_two_walkerstates) {
            walker_state_other=(walkerState*)calloc(N_walker, sizeof(walkerState));
        } else {
            walker_state_other=walker_state;
        }
    }
    walker_count=N_walker;
}

void FluorophorDynamics::init_walker(unsigned long i, double x, double y, double z) {
    walker_state[i].time=0;
    walker_state[i].exists=true;
    walker_state[i].x=x;
    walker_state[i].y=y;
    walker_state[i].z=z;
    walker_state[i].x0=x;
    walker_state[i].y0=y;
    walker_state[i].z0=z;
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
    walker_state[i].user_data=NULL;
    fluorophors->load_spectrum(init_spectrum);
}

void FluorophorDynamics::init(){
    sim_time=0;
    endoftrajectory=false;

    // ensure that all calculated variables are set properly.
    // set_sim_timestep also allocates memory for the walker states
    set_c_fluor(c_fluor);
    set_sim_timestep(sim_timestep);

    // now we initialize the walkers at random possitions inside the simulation volume
    if (volume_shape==Box) {
        for (unsigned long i=0; i<walker_count; i++) {
            init_walker(i, sim_x*gsl_rng_uniform(rng), sim_y*gsl_rng_uniform(rng), sim_z*gsl_rng_uniform(rng));
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

    store_step_protocol();

}

void FluorophorDynamics::propagate(bool boundary_check){
    sim_time=sim_time+sim_timestep;
    //std::cout<<">>> dyn    sim_time = "<<sim_time<<"\n";
}

void FluorophorDynamics::propagate_photophysics(int i) {
    if (!use_photophysics) return;
    if (walker_state[i].used_qm_states<=1) return;
    register int state=walker_state[i].qm_state;
    register double p=0;
    register double r=gsl_rng_uniform(rng);
    for (register int f=0; f<mmin(walker_state[i].used_qm_states,N_FLUORESCENT_STATES); f++) {
        p=p+walker_state[i].photophysics_transition[state*N_FLUORESCENT_STATES+f];
        if (p >= r) {
            walker_state[i].qm_state=f;
            break;
        }
    }

    /*switch(walker_state[i].qm_state) {
        case -2: return; break; // may not leave bleached state
        case -1: { // in triplet state
                   if (gsl_rng_uniform(rng)<sim_timestep/walker_state[i].triplet_lifetime) {
                       walker_state[i].qm_state=0; // return to ground state
                   }
                   return;
                 } break;
        case 0: { // in ground state
                   if (gsl_rng_uniform(rng)<walker_state[i].bleaching_propability) {
                       walker_state[i].qm_state=-2; // return to ground state
                   } else if (gsl_rng_uniform(rng)<walker_state[i].triplet_propability) {
                       walker_state[i].qm_state=-1; // return to ground state
                   }
                   return;
                 } break;
    }*/
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

    s+="c_fluor = "+floattostr(c_fluor)+" nMolar\n";
    s+="init_dipole_vec = ("+floattostr(init_p_x)+", "+floattostr(init_p_y)+", "+floattostr(init_p_z)+")\n";
    s+="init_spectrum = "+inttostr(init_spectrum);
    if (init_spectrum==-1) s+=" [none]\n";
    else s+=" ["+extract_file_name(fluorophors->getSpectrumFilename(init_spectrum))+"]\n";
    s+="init_q_fluor = "+doublevectortostr(init_q_fluor, N_FLUORESCENT_STATES)+" \n";
    s+="init_qm_state = "+inttostr(init_qm_state)+"\n";
    s+="init_type = "+inttostr(init_type)+"\n";
    s+="init_sigma_abs = "+doublevectortostr(init_sigma_abs, N_FLUORESCENT_STATES)+" meters^2\n";
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
    }


    return s;
}

void FluorophorDynamics::save_trajectories() {
    FILE* f;
    char fn[1024];
    if (protocol_trajectories>0) {
        for (unsigned int i=0; i<protocol_trajectories; i++) {
            fclose(trajectoryFile[i]);
        }
        free(trajectoryFile);
    }
};

void FluorophorDynamics::store_step_protocol() {
    for (unsigned int i=0; i<protocol_trajectories; i++) {
        walkerState ws=walker_state[i];
        if (ws.time<protocol_timestep_count || protocol_timestep_count<0)
            fprintf(trajectoryFile[i], "%lg, %lg, %lg, %lg, %d, %d, %d, %lg, %lg, %lg\n", ws.time*sim_timestep, ws.x, ws.y, ws.z, ws.qm_state, ws.type, ws.spectrum, ws.p_x, ws.p_y, ws.p_z);
    }
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


FluorophorDynamics::walkerState* FluorophorDynamics::copy_walker_state(FluorophorDynamics::walkerState* start) {
    FluorophorDynamics::walkerState* next=start+walker_count;
    memcpy(start, walker_state, walker_count*sizeof(FluorophorDynamics::walkerState));
    return next;
}


void FluorophorDynamics::load_all_used_spectra() {
    if (walker_count>0) for (unsigned long i=0; i<walker_count; i++) {
        fluorophors->load_spectrum(walker_state[i].spectrum);
    }
}



