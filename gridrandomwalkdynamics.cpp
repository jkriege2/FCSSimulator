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


#include "gridrandomwalkdynamics.h"

GridRandomWalkDynamics::GridRandomWalkDynamics(FluorophorManager* fluorophors, std::string object_name):
    FluorophorDynamics(fluorophors, object_name)
{
    save_msd_every_n_timesteps=-1;
    msd=NULL;
    msd2=NULL;
    msdg=NULL;
    msd2g=NULL;
    msd_count=NULL;
    obstacles=NULL;
    obstacle_size=0;
    msd_size=10000;
    gridsize_X=0;
    gridsize_Y=0;
    gridsize_Z=0;
    obstacle_fraction=0;
    //set_grid_constant(0.01);
    set_diff_coeff(1000);
    init();
}

GridRandomWalkDynamics::~GridRandomWalkDynamics() {
    if (msd) {
        free(msd);
        free(msd2);
        free(msdg);
        free(msd2g);
        free(msd_count);
    }
    if (obstacles) free(obstacles);
}

void GridRandomWalkDynamics::recalcJumpProp() {
    double tau=grid_constant*grid_constant/4.0/diff_coeff;
    jump_propability=sim_timestep/tau;
    //jump_propability=sqrt(6.0*diff_coeff*sim_timestep)
}

void GridRandomWalkDynamics::recalcGridSize() {
    grid_constant=sqrt(6.0*diff_coeff*sim_timestep);

    gridsize_X=ceil(sim_x/grid_constant);
    gridsize_Y=ceil(sim_y/grid_constant);
    gridsize_Z=ceil(sim_z/grid_constant);

    recalcJumpProp();
}

void GridRandomWalkDynamics::set_grid_constant(double gconst) {
    grid_constant=gconst;

    volume_shape= FluorophorDynamics::Box;
    gridsize_X=ceil(sim_x/gconst);
    gridsize_Y=ceil(sim_y/gconst);
    gridsize_Z=ceil(sim_z/gconst);

    recalcJumpProp();
}

void GridRandomWalkDynamics::set_diff_coeff(double value) {
    diff_coeff=value;
    //sigma_jump[i]=sqrt(2.0*diff_coeff[i]*sim_timestep);
    recalcJumpProp();
    recalcGridSize();
};


void GridRandomWalkDynamics::set_sim_timestep(double value) {
    FluorophorDynamics::set_sim_timestep(value);
    recalcJumpProp();
    /*for (int i=0; i<DCOUNT; i++) {
        sigma_jump[i]=sqrt(2.0*diff_coeff[i]*sim_timestep);
    }*/
    recalcGridSize();

};


void GridRandomWalkDynamics::init_walker(unsigned long i, double x, double y, double z) {
    walker_state[i].time=0;
    walker_state[i].exists=true;
    walker_state[i].x=x;
    walker_state[i].y=y;
    walker_state[i].z=z;
    walker_state[i].ix=floor(x/grid_constant);
    walker_state[i].iy=floor(y/grid_constant);
    walker_state[i].iz=floor(z/grid_constant);
    walker_state[i].x0=x;
    walker_state[i].y0=y;
    walker_state[i].z0=z;
    walker_state[i].ix0=walker_state[i].ix;
    walker_state[i].iy0=walker_state[i].iy;
    walker_state[i].iz0=walker_state[i].iz;
    walker_state[i].p_x=init_p_x;
    walker_state[i].p_y=init_p_y;
    walker_state[i].p_z=init_p_z;
    //for (int j=0; j<N_FLUORESCENT_STATES; j++) walker_state[i].sigma_abs[j]=init_sigma_abs[j];
    //for (int j=0; j<N_FLUORESCENT_STATES; j++) walker_state[i].q_fluor[j]=init_q_fluor[j];
    //for (int j=0; j<N_FLUORESCENT_STATES*N_FLUORESCENT_STATES; j++) walker_state[i].photophysics_transition[j]=init_photophysics_transition[j];
    walker_state[i].qm_state=init_qm_state;
    walker_state[i].used_qm_states=init_used_qm_states;
    walker_state[i].type=init_type;
    walker_state[i].spectrum=init_spectrum;
    //walker_state[i].user_data=NULL;
    fluorophors->load_spectrum(init_spectrum);
}

void GridRandomWalkDynamics::read_config_internal(jkINIParser2& parser) {
    FluorophorDynamics::read_config_internal(parser);
    set_diff_coeff(parser.getSetAsDouble("diff_coeff", diff_coeff));
    save_msd_every_n_timesteps=parser.getSetAsInt("save_msd_every_n_timesteps", save_msd_every_n_timesteps);
    msd_size=parser.getSetAsInt("msd_size", msd_size);
    obstacle_fraction=parser.getSetAsDouble("obstacle_fraction", obstacle_fraction);
    //set_grid_constant(parser.getSetAsDouble("grid_constant", grid_constant));
    //init();
}


void GridRandomWalkDynamics::init(){
    FluorophorDynamics::init();
    if (msd!=NULL) {
        free(msd);
        free(msd2);
        free(msdg);
        free(msd2g);
        free(msd_count);
    }
    if (save_msd_every_n_timesteps>0) {
        msd=(double*)calloc(msd_size, sizeof(double));
        msd2=(double*)calloc(msd_size, sizeof(double));
        msdg=(double*)calloc(msd_size, sizeof(double));
        msd2g=(double*)calloc(msd_size, sizeof(double));
        msd_count=(uint64_t*)calloc(msd_size, sizeof(uint64_t));
        for (int i=0; i<msd_size; i++) {
            msd[i]=0;
            msd2[i]=0;
            msdg[i]=0;
            msd2g[i]=0;
            msd_count[i]=0;
        }
    }
    init_obstacles();
}

void GridRandomWalkDynamics::init_obstacles() {
    if (obstacles) free(obstacles);
    int64_t cells=(gridsize_X+3)*(gridsize_Y+3)*(gridsize_Z+3);
    std::cout<<"allocating "<<bytestostr(cells)<<" for obstacle grid ("<<gridsize_X+3<<" * "<<gridsize_Y+3<<" * "<<gridsize_Z+3<<") ... ";
    obstacles=(uint8_t*)calloc(cells, sizeof(uint8_t));
    obstacle_size=cells;
    if (obstacles) std::cout<<"OK!\n";
    else std::cout<<"ERROR!\n";
    if (!obstacles) {
        throw FluorophorException("could not allocate "+bytestostr(cells)+" of memory for obstacle grid!");
    }
    if (obstacle_fraction>0 && obstacle_fraction<1.0) {
        numobstacles=0;
        std::cout<<"placing obstacles on a  "<<gridsize_X+3<<" * "<<gridsize_Y+3<<" * "<<gridsize_Z+3<<"  grid (fraction = "+floattostr(obstacle_fraction)+") ... ";
        for (int64_t i=0; i<obstacle_size; i++) {
            obstacles[i]=0;
            if (gsl_rng_uniform(rng)<=obstacle_fraction) {
                obstacles[i]=1;
                numobstacles++;
            }
        }
        std::cout<<"DONE!\n";
        std::cout<<"  placed "<<numobstacles<<" obstacles (real fraction: "<<double(numobstacles)/double(cells)<<"\n";
    }
}


void GridRandomWalkDynamics::set_sim_box(double vx, double vy, double vz) {
    FluorophorDynamics::set_sim_box(vx, vy, vz);

    gridsize_X=ceil(sim_x/grid_constant);
    gridsize_Y=ceil(sim_y/grid_constant);
    gridsize_Z=ceil(sim_z/grid_constant);

    recalcJumpProp();

    recalcGridSize();
}

void GridRandomWalkDynamics::set_sim_sphere(double rad) {
    set_sim_box(rad*2.0, rad*2.0, rad*2.0);
}


void GridRandomWalkDynamics::propagate(bool boundary_check){
    FluorophorDynamics::propagate(boundary_check);

    if (volume_shape != FluorophorDynamics::Box)
        throw FluorophorException("random walks on grids need a box-shaped simulation volume (other options are not allowd)");

    //static walkerState oldstate;
    // now wepropagate every walker
    for (register unsigned long i=0; i<walker_count; i++) {
        if (walker_state[i].exists) {
            walker_state[i].time++;

            /*if (jump_propability>=1 || gsl_rng_uniform(rng)>jump_propability) */
            {
                // translational diffusion
                register int32_t x0=walker_state[i].ix;
                register int32_t y0=walker_state[i].iy;
                register int32_t z0=walker_state[i].iz;

                register int32_t nx,ny,nz;

                const int dx[6]={-1, 1, 0, 0, 0, 0};
                const int dy[6]={0, 0, -1, 1, 0, 0};
                const int dz[6]={0, 0, 0, 0, -1, 1};

                int direction=gsl_rng_uniform_int (rng, 6);
                nx=x0+dx[direction];
                ny=y0+dy[direction];
                nz=z0+dz[direction];


                // update walker position
                int64_t idx=nz*(gridsize_X+3)*(gridsize_Y+3)+ny*(gridsize_X+3)+nx;
                if (idx<obstacle_size && idx>0 && obstacles[idx]==0) {
                    walker_state[i].ix=nx;
                    walker_state[i].iy=ny;
                    walker_state[i].iz=nz;
                    walker_state[i].x=nx*grid_constant;
                    walker_state[i].y=ny*grid_constant;
                    walker_state[i].z=nz*grid_constant;
                }
            }
            // check whether walker left the volume and if so: delete the walker and add
            // anew one at some position on the border of the simulation volume
            if (boundary_check) {
                //std::cout<<"boundary ";
                perform_boundary_check(i);
            }

            register double nx=walker_state[i].x;
            register double ny=walker_state[i].y;
            register double nz=walker_state[i].z;

            //std::cout<<"msd ";

            if ((save_msd_every_n_timesteps>0)&&((int64_t)(walker_state[i].time%save_msd_every_n_timesteps)==(save_msd_every_n_timesteps-1))) {
                double m=gsl_pow_2(nx-walker_state[i].x0)+gsl_pow_2(ny-walker_state[i].y0)+gsl_pow_2(nz-walker_state[i].z0);
                double mg=gsl_pow_2(nx-walker_state[i].x0)/gsl_pow_2(grid_constant)+gsl_pow_2(ny-walker_state[i].y0)/gsl_pow_2(grid_constant)+gsl_pow_2(nz-walker_state[i].z0)/gsl_pow_2(grid_constant);
                int idx=walker_state[i].time/save_msd_every_n_timesteps;
                if ((idx>=0)&&(idx<msd_size)) {
                    msd[idx] += m;
                    msd2[idx] += m*m;
                    msdg[idx] += mg;
                    msd2g[idx] += mg*mg;
                    msd_count[idx]++;
                }
            }



            //std::cout<<"photophys ";

            // photophysics
            propagate_photophysics(i);
            //std::cout<<"prop "<<i<<std::endl;


        }
    }
    store_step_protocol();

    if ((save_msd_every_n_timesteps>0) && (sim_time>duration-sim_timestep)) {
        // store MSD curve
        char fn[1024];
        sprintf(fn, "%s%smsd.dat", basename.c_str(), object_name.c_str());
        std::cout<<"writing MSD file '"<<fn<<"' ... ";
        FILE* f=fopen(fn, "w");
        for (int i=0; i<msd_size; i++) {
            double tau=(double)(i)*(double)sim_timestep*(double)save_msd_every_n_timesteps;
            double taug=(double)(i)*(double)save_msd_every_n_timesteps;
            double cnt=(double)msd_count[i];
            double mean=msd[i]/cnt;
            double stddev=sqrt(msd2[i]/cnt-mean*mean);
            double meang=msdg[i]/cnt;
            double stddevg=sqrt(msd2g[i]/cnt-meang*meang);
            if (msd_count[i]>0) fprintf(f, "%lg, %lg, %lg, %lg, %lg, %lg, %lg\n", tau, mean, stddev, cnt, taug, meang, stddevg);
        }
        fclose(f);
        std::cout<<"DONE!\n";

        sprintf(fn, "%s%smsd.plt", basename.c_str(), object_name.c_str());
        std::cout<<"writing MSD plot file '"<<fn<<"' ... ";
        f=fopen(fn, "w");
        fprintf(f, "unset multiplot\n");
        fprintf(f, "reset\n");
        fprintf(f, "msd_fit(tau) = 6 * Diff * tau\n");
        fprintf(f, "Diff=%lf;\n", diff_coeff);
        fprintf(f, "fit msd_fit(x) \"%s\" using 1:2 via Diff\n", extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "msda_fit(tau) = 6 * Gamma * (tau**alpha)\n");
        fprintf(f, "Gamma=%lf;\n", diff_coeff);
        fprintf(f, "alpha=1;\n");
        fprintf(f, "fit msd_fit(x) \"%s\" using 1:2 via Diff\n", extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "fit msda_fit(x) \"%s\" using 1:2 via Gamma, alpha\n", extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "set title \"sigma_translation^2 against time\"\n");
        fprintf(f, "set xlabel \"time t [s]\"\n");
        fprintf(f, "set ylabel \"sigma^2 [microns^2]\"\n");
        fprintf(f, "msd(tau)=6 * %lf * tau\n", diff_coeff);
        fprintf(f, "set style fill transparent solid 0.5 noborder\n");
        fprintf(f, "plot \"%s\" using 1:(($2)-($3)):(($2)+($3)) notitle with filledcurves, \"%s\" using 1:2 title \"simulation result\" with points, "
                   "msd(x) title \"theory\" with lines, "
                   "msd_fit(x) title sprintf(\"fit D=%%f\",Diff) with lines, "
                   "msda_fit(x) title sprintf(\"fit Gamma=%%f alpha=%%f\",Gamma,alpha) with lines\n",
                   extract_file_name(basename+object_name+"msd.dat").c_str(), extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "\n");
        fprintf(f, "pause -1\n");
        fprintf(f, "set logscale xy\n");
        fprintf(f, "set title \"sigma_translation^2 against time\"\n");
        fprintf(f, "set xlabel \"time t [s]\"\n");
        fprintf(f, "set ylabel \"sigma^2 [microns^2]\"\n");
        fprintf(f, "msd(tau)=6 * %lf * tau\n", diff_coeff);
        fprintf(f, "set style fill transparent solid 0.5 noborder\n");
        fprintf(f, "plot \"%s\" using 1:(($2)-($3)):(($2)+($3)) notitle with filledcurves, "
                   "\"%s\" using 1:2 title \"simulation result\" with points, "
                   "msd(x) title \"theory\" with lines, "
                   "msd_fit(x) title sprintf(\"fit D=%%f\",Diff) with lines, "
                   "msda_fit(x) title sprintf(\"fit Gamma=%%f alpha=%%f\",Gamma,alpha) with lines\n", extract_file_name(basename+object_name+"msd.dat").c_str(), extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "\n");
        fprintf(f, "pause -1\n");
        fprintf(f, "unset logscale xy\n");
        fprintf(f, "set title \"particle averaged over against time\"\n");
        fprintf(f, "set xlabel \"time t [s]\"\n");
        fprintf(f, "set ylabel \"particles\"\n");
        fprintf(f, "plot \"%s\" using 1:4 title \"simulation result\" with linespoints\n", extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "\n");
        fprintf(f, "pause -1\n");
        fclose(f);
        std::cout<<"DONE!\n";

        sprintf(fn, "%s%smsdsimple.plt", basename.c_str(), object_name.c_str());
        std::cout<<"writing simpleMSD plot file '"<<fn<<"' ... ";
        f=fopen(fn, "w");
        fprintf(f, "unset multiplot\n");
        fprintf(f, "reset\n");
        fprintf(f, "set title \"sigma_translation^2 against time\"\n");
        fprintf(f, "set xlabel \"time t [s]\"\n");
        fprintf(f, "set ylabel \"sigma^2 [microns^2]\"\n");
        fprintf(f, "msd(tau)=6 * %lf * tau\n", diff_coeff);
        fprintf(f, "set style fill transparent solid 0.5 noborder\n");
        fprintf(f, "plot \"%s\" using 1:(($2)-($3)):(($2)+($3)) notitle with filledcurves, \"%s\" using 1:2 title \"simulation result\" with points, "
                   "msd(x) title \"theory\" with lines\n",
                   extract_file_name(basename+object_name+"msd.dat").c_str(), extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "\n");
        fprintf(f, "pause -1\n");
        fprintf(f, "set logscale xy\n");
        fprintf(f, "set title \"sigma_translation^2 against time\"\n");
        fprintf(f, "set xlabel \"time t [s]\"\n");
        fprintf(f, "set ylabel \"sigma^2 [microns^2]\"\n");
        fprintf(f, "msd(tau)=6 * %lf * tau\n", diff_coeff);
        fprintf(f, "set style fill transparent solid 0.5 noborder\n");
        fprintf(f, "plot \"%s\" using 1:(($2)-($3)):(($2)+($3)) notitle with filledcurves, "
                   "\"%s\" using 1:2 title \"simulation result\" with points, "
                   "msd(x) title \"theory\" with lines\n", extract_file_name(basename+object_name+"msd.dat").c_str(), extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "\n");
        fprintf(f, "pause -1\n");
        fprintf(f, "unset logscale xy\n");
        fprintf(f, "set title \"particle averaged over against time\"\n");
        fprintf(f, "set xlabel \"time t [s]\"\n");
        fprintf(f, "set ylabel \"particles\"\n");
        fprintf(f, "plot \"%s\" using 1:4 title \"simulation result\" with linespoints\n", extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "\n");
        fprintf(f, "pause -1\n");
        fclose(f);
        std::cout<<"DONE!\n";

        sprintf(fn, "%s%smsdgrid.plt", basename.c_str(), object_name.c_str());
        std::cout<<"writing grid MSD plot file '"<<fn<<"' ... ";
        f=fopen(fn, "w");
        fprintf(f, "unset multiplot\n");
        fprintf(f, "reset\n");
        fprintf(f, "msd_fit(tau) = 6 * Diff * tau\n");
        fprintf(f, "Diff=%lf;\n", diff_coeff/grid_constant/grid_constant);
        fprintf(f, "fit msd_fit(x) \"%s\" using 5:6 via Diff\n", extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "msda_fit(tau) = 6 * Gamma * (tau**alpha)\n");
        fprintf(f, "Gamma=%lf;\n", diff_coeff/grid_constant/grid_constant*sim_timestep);
        fprintf(f, "alpha=1;\n");
        fprintf(f, "fit msd_fit(x) \"%s\" using 5:6 via Diff\n", extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "fit msda_fit(x) \"%s\" using 5:6 via Gamma, alpha\n", extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "set title \"sigma_translation^2 against time\"\n");
        fprintf(f, "set xlabel \"time t [timesteps]\"\n");
        fprintf(f, "set ylabel \"sigma^2 [gridconst^2]\"\n");
        fprintf(f, "msd(tau)=6 * %lf * tau\n", diff_coeff/grid_constant/grid_constant*sim_timestep);
        fprintf(f, "set style fill transparent solid 0.5 noborder\n");
        fprintf(f, "plot \"%s\" using 5:(($6)-($7)):(($6)+($7)) notitle with filledcurves, \"%s\" using 5:6 title \"simulation result\" with points, "
                   "msd(x) title \"theory\" with lines, "
                   "msd_fit(x) title sprintf(\"fit D=%%f\",Diff) with lines, "
                   "msda_fit(x) title sprintf(\"fit Gamma=%%f alpha=%%f\",Gamma,alpha) with lines\n",
                   extract_file_name(basename+object_name+"msd.dat").c_str(), extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "\n");
        fprintf(f, "pause -1\n");
        fprintf(f, "set logscale xy\n");
        fprintf(f, "set title \"sigma_translation^2 against time\"\n");
        fprintf(f, "set xlabel \"time t [timesteps]\"\n");
        fprintf(f, "set ylabel \"sigma^2 [gridconst^2]\"\n");
        fprintf(f, "msd(tau)=6 * %lf * tau\n", diff_coeff/grid_constant/grid_constant*sim_timestep);
        fprintf(f, "set style fill transparent solid 0.5 noborder\n");
        fprintf(f, "plot \"%s\" using 5:(($6)-($7)):(($6)+($7)) notitle with filledcurves, "
                   "\"%s\" using 5:6 title \"simulation result\" with points, "
                   "msd(x) title \"theory\" with lines, "
                   "msd_fit(x) title sprintf(\"fit D=%%f\",Diff) with lines, "
                   "msda_fit(x) title sprintf(\"fit Gamma=%%f alpha=%%f\",Gamma,alpha) with lines\n", extract_file_name(basename+object_name+"msd.dat").c_str(), extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "\n");
        fprintf(f, "pause -1\n");
        fprintf(f, "unset logscale xy\n");
        fprintf(f, "set title \"particle averaged over against time\"\n");
        fprintf(f, "set xlabel \"time t [timesteps]\"\n");
        fprintf(f, "set ylabel \"particles\"\n");
        fprintf(f, "plot \"%s\" using 5:4 title \"simulation result\" with linespoints\n", extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "\n");
        fprintf(f, "pause -1\n");
        fclose(f);
        std::cout<<"DONE!\n";


        sprintf(fn, "%s%smsdgridsimple.plt", basename.c_str(), object_name.c_str());
        std::cout<<"writing grid MSD plot file '"<<fn<<"' ... ";
        f=fopen(fn, "w");
        fprintf(f, "unset multiplot\n");
        fprintf(f, "reset\n");
        fprintf(f, "set title \"sigma_translation^2 against time\"\n");
        fprintf(f, "set xlabel \"time t [timesteps]\"\n");
        fprintf(f, "set ylabel \"sigma^2 [gridconst^2]\"\n");
        fprintf(f, "msd(tau)=6 * %lf * tau\n", diff_coeff/grid_constant/grid_constant*sim_timestep);
        fprintf(f, "set style fill transparent solid 0.5 noborder\n");
        fprintf(f, "plot \"%s\" using 5:(($6)-($7)):(($6)+($7)) notitle with filledcurves, \"%s\" using 5:6 title \"simulation result\" with points, "
                   "msd(x) title \"theory\" with lines, "
                   "msd_fit(x) title sprintf(\"fit D=%%f\",Diff) with lines, "
                   "msda_fit(x) title sprintf(\"fit Gamma=%%f alpha=%%f\",Gamma,alpha) with lines\n",
                   extract_file_name(basename+object_name+"msd.dat").c_str(), extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "\n");
        fprintf(f, "pause -1\n");
        fprintf(f, "set logscale xy\n");
        fprintf(f, "set title \"sigma_translation^2 against time\"\n");
        fprintf(f, "set xlabel \"time t [timesteps]\"\n");
        fprintf(f, "set ylabel \"sigma^2 [gridconst^2]\"\n");
        fprintf(f, "msd(tau)=6 * %lf * tau\n", diff_coeff/grid_constant/grid_constant*sim_timestep);
        fprintf(f, "set style fill transparent solid 0.5 noborder\n");
        fprintf(f, "plot \"%s\" using 5:(($6)-($7)):(($6)+($7)) notitle with filledcurves, "
                   "\"%s\" using 5:6 title \"simulation result\" with points, "
                   "msd(x) title \"theory\" with lines, "
                   "msd_fit(x) title sprintf(\"fit D=%%f\",Diff) with lines, "
                   "msda_fit(x) title sprintf(\"fit Gamma=%%f alpha=%%f\",Gamma,alpha) with lines\n", extract_file_name(basename+object_name+"msd.dat").c_str(), extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "\n");
        fprintf(f, "pause -1\n");
        fprintf(f, "unset logscale xy\n");
        fprintf(f, "set title \"particle averaged over against time\"\n");
        fprintf(f, "set xlabel \"time t [timesteps]\"\n");
        fprintf(f, "set ylabel \"particles\"\n");
        fprintf(f, "plot \"%s\" using 5:4 title \"simulation result\" with linespoints\n", extract_file_name(basename+object_name+"msd.dat").c_str());
        fprintf(f, "\n");
        fprintf(f, "pause -1\n");
        fclose(f);
        std::cout<<"DONE!\n";
    }
    //std::cout<<"propagate done\n";
}

std::string GridRandomWalkDynamics::report() {
    std::string s=FluorophorDynamics::report();
    s+="gridsize_X = "+inttostr(gridsize_X)+"\n";
    s+="gridsize_Y = "+inttostr(gridsize_Y)+"\n";
    s+="gridsize_Z = "+inttostr(gridsize_Z)+"\n";    s+="diff_coeff = "+floattostr(diff_coeff)+" micron^2/s\n";
    //s+="jump_propability = "+floattostr(jump_propability)+"\n";
    s+="grid_constant = "+floattostr(grid_constant)+" micron\n";
    s+="  => diff_coeff_from_grid = "+floattostr(grid_constant*grid_constant/6.0/sim_timestep)+" micron^2/s\n";
    s+="set_obstacle_fraction = "+floattostr(obstacle_fraction)+"\n";
    s+="numobstacles = "+inttostr(numobstacles)+"\n";
    s+="real_obstacle_fraction = "+floattostr(double(numobstacles)/double(obstacle_size))+"\n";

    if (save_msd_every_n_timesteps>0) s+="saving MSD every "+inttostr(save_msd_every_n_timesteps)+" timesteps with "+inttostr(msd_size)+" items\n";

    return s;
}

std::string GridRandomWalkDynamics::dot_get_properties() {
    std::string s=FluorophorDynamics::dot_get_properties();
    s+="gridsize_X = "+inttostr(gridsize_X)+"<BR/>";
    s+="gridsize_Y = "+inttostr(gridsize_Y)+"<BR/>";
    s+="gridsize_Z = "+inttostr(gridsize_Z)+"<BR/>";
    s+="diff_coeff = "+floattostr(diff_coeff)+" &mu;m&sup2;/s<BR/>";
    //s+="jump_propability = "+floattostr(jump_propability)+"\n";
    s+="grid_constant = "+floattostr(grid_constant)+" &mu;m<BR/>";
    s+="  => diff_coeff_from_grid = "+floattostr(grid_constant*grid_constant/6.0/sim_timestep)+" &mu;m&sup2;/s<BR/>";
    s+="set_obstacle_fraction = "+floattostr(obstacle_fraction)+"<BR/>";
    s+="numobstacles = "+inttostr(numobstacles)+"<BR/>";
    s+="real_obstacle_fraction = "+floattostr(double(numobstacles)/double(obstacle_size))+"<BR/>";

    return s;
}

void GridRandomWalkDynamics::perform_boundary_check(unsigned long i) {
    register int32_t nx=walker_state[i].ix;
    register int32_t ny=walker_state[i].iy;
    register int32_t nz=walker_state[i].iz;

    if (walker_state[i].exists) {
        if (   (nx<0) || (nx>=gridsize_X)
            || (ny<0) || (ny>=gridsize_Y)
            || (nz<0) || (nz>=gridsize_Z) ) {

            if (depletion_propability<=0 || gsl_ran_flat(rng,0,1)>depletion_propability) {
                //std::cout<<"initializing new walker ... depletion_propability="<<depletion_propability<<"\n";


                // first choose one face of the simulation volume and then set the walker
                // to any position on the face ... also shift a bit inwards
                char face=gsl_rng_uniform_int(rng, 6)+1;
                switch(face) {
                    case 1:
                        //x-y-plane at z=0
                        walker_state[i].ix=gsl_rng_uniform_int(rng, gridsize_X);
                        walker_state[i].iy=gsl_rng_uniform_int(rng, gridsize_Y);
                        walker_state[i].iz=0;
                        break;
                    case 2:
                        //x-y-plane at z=gridsize_Z
                        walker_state[i].ix=gsl_rng_uniform_int(rng, gridsize_X);
                        walker_state[i].iy=gsl_rng_uniform_int(rng, gridsize_Y);
                        walker_state[i].iz=gridsize_Z;
                        break;
                    case 3:
                        //x-z-plane at y=0
                        walker_state[i].ix=gsl_rng_uniform_int(rng, gridsize_X);
                        walker_state[i].iy=0;
                        walker_state[i].iz=gsl_rng_uniform_int(rng, gridsize_Z);
                        break;
                    case 4:
                        //x-z-plane at y=gridsize_Y
                        walker_state[i].ix=gsl_rng_uniform_int(rng, gridsize_X);
                        walker_state[i].iy=gridsize_Y;
                        walker_state[i].iz=gsl_rng_uniform_int(rng, gridsize_Z);
                        break;
                    case 5:
                        //z-y-plane at x=0
                        walker_state[i].ix=0;
                        walker_state[i].iy=gsl_rng_uniform_int(rng, gridsize_Y);
                        walker_state[i].iz=gsl_rng_uniform_int(rng, gridsize_Z);
                        break;
                    case 6:
                        //z-y-plane at x=gridsize_X
                        walker_state[i].ix=gridsize_X;
                        walker_state[i].iy=gsl_rng_uniform_int(rng, gridsize_Y);
                        walker_state[i].iz=gsl_rng_uniform_int(rng, gridsize_Z);
                        break;
                }
                walker_state[i].time=0;
                walker_state[i].ix0=walker_state[i].ix;
                walker_state[i].iy0=walker_state[i].iy;
                walker_state[i].iz0=walker_state[i].iz;
                walker_state[i].x=walker_state[i].ix*grid_constant;
                walker_state[i].y=walker_state[i].iy*grid_constant;
                walker_state[i].z=walker_state[i].iz*grid_constant;
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
