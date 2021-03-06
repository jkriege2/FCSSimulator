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


#include "browniandynamics.h"

#include <gsl/gsl_histogram.h>

BrownianDynamics::BrownianDynamics(FluorophorManager* fluorophors, std::string object_name):
    FluorophorDynamics(fluorophors, object_name)
{
    save_msd_every_n_timesteps=-1;
    msd=NULL;
    velocity[0]=velocity[1]=velocity[2]=0;
    msd2=NULL;
    msd_count=NULL;
    msd_size=10000;
    set_diff_coeff(10);
    diff_rot=100;
    use_rotational_diffusion=false;
    n_fluorophores_local=1;
    init();
}

BrownianDynamics::BrownianDynamics(FluorophorManager* fluorophors, double sim_x, double sim_y, double sim_z, double c_fluor, std::string object_name):
    FluorophorDynamics(fluorophors, sim_x, sim_y, sim_z, c_fluor, object_name)
{
    save_msd_every_n_timesteps=-1;
    msd=NULL;
    velocity[0]=velocity[1]=velocity[2]=0;
    msd2=NULL;
    msd_size=10000;
    msd_count=NULL;
    set_diff_coeff(10);
    diff_rot=100;
    use_rotational_diffusion=false;
    n_fluorophores_local=1;
    init();
}

BrownianDynamics::BrownianDynamics(FluorophorManager* fluorophors, double sim_radius, double c_fluor, std::string object_name):
    FluorophorDynamics(fluorophors, sim_radius, c_fluor, object_name)
{
    save_msd_every_n_timesteps=-1;
    msd=NULL;
    velocity[0]=velocity[1]=velocity[2]=0;
    msd2=NULL;
    msd_size=10000;
    msd_count=NULL;
    set_diff_coeff(20);
    n_fluorophores_local=1;
    diff_rot=100;
    use_rotational_diffusion=false;
    init();
}

BrownianDynamics::~BrownianDynamics()
{
    if (msd) {
        free(msd);
        free(msd2);
        free(msd_count);
    }
}

unsigned long BrownianDynamics::calc_walker_count() {
    unsigned long c=FluorophorDynamics::calc_walker_count();
    return c*n_fluorophores_local;
}

void BrownianDynamics::read_config_internal(jkINIParser2& parser) {
    n_fluorophores_local=parser.getSetAsInt("n_fluorophores_local", n_fluorophores_local);
    FluorophorDynamics::read_config_internal(parser);
    //diffarea_x0=parser.getAsDouble("diffarea_x0", this->sim_x/2.0);
    //set_diff_coeff(parser.getAsDouble("diff_coeff", diff_coeff[0]), parser.getAsDouble("diff_coeff1", diff_coeff[1]));
    set_diff_coeff(parser.getSetAsDouble("diff_coeff", parser.getSetAsDouble("diff_coeff0", diff_coeff)));

    set_diffrot_coeff(parser.getSetAsDouble("diff_rot", diff_rot/M_PI/M_PI*180.0*180.0)/180.0/180.0*M_PI*M_PI);
    use_rotational_diffusion=parser.getSetAsBool("use_rotational_diffusion", use_rotational_diffusion);
    save_msd_every_n_timesteps=parser.getSetAsInt("save_msd_every_n_timesteps", save_msd_every_n_timesteps);
    msd_size=parser.getSetAsInt("msd_size", msd_size);
    velocity[0]=parser.getSetAsDouble("vx", parser.getSetAsDouble("v0", velocity[0]));
    velocity[1]=parser.getSetAsDouble("vy", parser.getSetAsDouble("v1", velocity[1]));
    velocity[2]=parser.getSetAsDouble("vz", parser.getSetAsDouble("v2", velocity[2]));
    init();
}


void BrownianDynamics::init(){
    FluorophorDynamics::init();
    if (msd!=NULL) {
        free(msd);
        free(msd2);
        free(msd_count);
    }
    if (save_msd_every_n_timesteps>0) {
        msd=(double*)calloc(msd_size, sizeof(double));
        msd2=(double*)calloc(msd_size, sizeof(double));
        msd_count=(uint64_t*)calloc(msd_size, sizeof(uint64_t));
        for (int i=0; i<msd_size; i++) {
            msd[i]=0;
            msd2[i]=0;
            msd_count[i]=0;
        }
    }
}

void BrownianDynamics::propagate(bool boundary_check){
    FluorophorDynamics::propagate(boundary_check);
    sigma_jump=sqrt(2.0*diff_coeff*sim_timestep);
    //static walkerState oldstate;
    // now wepropagate every walker
    for (register unsigned long i=0; i<walker_count; i++) {
        if (walker_state[i].exists) {
            if (i%n_fluorophores_local==0) {
                walker_state[i].time++;

                // translational diffusion
                register double x0=walker_state[i].x;
                register double y0=walker_state[i].y;
                register double z0=walker_state[i].z;

                register double nx,ny,nz;

                //double s=(x0>=diffarea_x0)?sigma_jump1:sigma_jump;
                register double s=sigma_jump;

                //if (x0>=diffarea_x0) s=sigma_jump1;
                nx=x0+gsl_ran_gaussian_ziggurat(rng, s)+velocity[0]*sim_timestep;
                ny=y0+gsl_ran_gaussian_ziggurat(rng, s)+velocity[1]*sim_timestep;
                nz=z0+gsl_ran_gaussian_ziggurat(rng, s)+velocity[2]*sim_timestep;

                // update walker position
                walker_state[i].x=nx;
                walker_state[i].y=ny;
                walker_state[i].z=nz;

                if ((save_msd_every_n_timesteps>0)&&msd&&msd2&&msd_count&&((int64_t)(walker_state[i].time%save_msd_every_n_timesteps)==(save_msd_every_n_timesteps-1))) {
                    double m=gsl_pow_2(nx-walker_state[i].x0)+gsl_pow_2(ny-walker_state[i].y0)+gsl_pow_2(nz-walker_state[i].z0);
                    int idx=walker_state[i].time/save_msd_every_n_timesteps;
                    //std::cout<<save_msd_every_n_timesteps<<", "<<msd<<", "<<msd_size<<", "<<idx<<"\n";
                    if ((idx>=0)&&(idx<msd_size)) {
                        msd[idx] += m;
                        msd2[idx] += m*m;
                        msd_count[idx]++;
                    }
                }

                // check whether walker left the volume and if so: delete the walker and add
                // anew one at some position on the border of the simulation volume
                if (boundary_check) {
                    perform_boundary_check(i);
                }



                // rotational diffusion
                if (use_rotational_diffusion) {
                    register double mu1=walker_state[i].p_x;
                    register double mu2=walker_state[i].p_y;
                    register double mu3=walker_state[i].p_z;
                    //double mu;

                    double ns1, ns2, ns3, nss1, nss2, nss3, ns, nss;
                    if (mu3!=0) {
                        ns1=1;
                        ns2=2;
                        ns3=-(mu1+2.0*mu2)/mu3;
                        ns=sqrt(1+4+ns3*ns3);
                    } else if (mu2!=0) {
                        ns1=1;
                        ns2=-(mu1+2.0*mu3)/mu2;
                        ns3=2;
                        ns=sqrt(1+4+ns2*ns2);
                    } else {
                        ns1=-(mu2+2.0*mu3)/mu1;
                        ns2=1;
                        ns3=1;
                        ns=sqrt(1+4+ns1*ns1);
                    }
                    ns1=ns1/ns; ns2=ns2/ns; ns3=ns3/ns;

                    nss1=mu2*ns3-mu3*ns2;
                    nss2=mu3*ns1-mu1*ns3;
                    nss3=mu1*ns2-mu2*ns1;
                    nss=sqrt(nss1*nss1+nss2*nss2+nss3*nss3);
                    nss1=nss1/nss; nss2=nss2/nss; nss3=nss3/nss;

                    double phi=gsl_ran_flat(rng, -M_PI/2.0, M_PI/2.0);
                    double dalpha=gsl_ran_gaussian_ziggurat(rng, sigma_rotjump);
                    /*while (fabs(dalpha)>=M_PI/2.0) {
                        dalpha=gsl_ran_gaussian_ziggurat(rng, sigma_rotjump);
                    }*/
                    double dmu1, dmu2, dmu3;//, dmu;
                    dmu1=cos(phi)*ns1+sin(phi)*nss1;
                    dmu2=cos(phi)*ns2+sin(phi)*nss2;
                    dmu3=cos(phi)*ns3+sin(phi)*nss3;
                    //dmu=sqrt(dmu1*dmu1+dmu2*dmu2+dmu3*dmu3);
                    /*double x=sqrt(1.0/gsl_pow_2(cos(dalpha))-1.0)/dmu;

                    if (dalpha<0) {
                        x=-1.0*x;
                    }
                    mu1=mu1+x*dmu1;
                    mu2=mu2+x*dmu2;
                    mu3=mu3+x*dmu3;
                    mu=sqrt(mu1*mu1+mu2*mu2+mu3*mu3);
                    mu1=mu1/mu; mu2=mu2/mu; mu3=mu3/mu;*/

                    // normalized vector (mu x dmu), i.e. rotation axis
                    nss1=mu2*dmu3-mu3*dmu2;
                    nss2=mu3*dmu1-mu1*dmu3;
                    nss3=mu1*dmu2-mu2*dmu1;
                    nss=sqrt(nss1*nss1+nss2*nss2+nss3*nss3);
                    nss1=nss1/nss; nss2=nss2/nss; nss3=nss3/nss;

                    double ca=cos(dalpha);
                    double sa=sin(dalpha);
                    double R[3][3]= {{ca+nss1*nss1*(1.0-ca),      nss1*nss2*(1.0-ca)-nss3*sa, nss1*nss3*(1.0-ca)+nss2*sa},
                                     {nss2*nss1*(1.0-ca)+nss3*sa, ca+nss2*nss2*(1.0-ca),      nss2*nss3*(1.0-ca)-nss1*sa},
                                     {nss3*nss1*(1.0-ca)-nss2*sa, nss3*nss2*(1.0-ca)+nss1*sa, ca+nss3*nss3*(1.0-ca)}};


                    walker_state[i].p_x=R[0][0]*mu1+R[0][1]*mu2+R[0][2]*mu3;
                    walker_state[i].p_y=R[1][0]*mu1+R[1][1]*mu2+R[1][2]*mu3;
                    walker_state[i].p_z=R[2][0]*mu1+R[2][1]*mu2+R[2][2]*mu3;
                    double l=sqrt(walker_state[i].p_x*walker_state[i].p_x+walker_state[i].p_y*walker_state[i].p_y+walker_state[i].p_z*walker_state[i].p_z);
                    walker_state[i].p_x=walker_state[i].p_x/l;
                    walker_state[i].p_y=walker_state[i].p_y/l;
                    walker_state[i].p_z=walker_state[i].p_z/l;

                }

                // photophysics
                propagate_photophysics(i);
                //std::cout<<"prop "<<i<<std::endl;

            } else {
                walker_state[i].x=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].x;
                walker_state[i].y=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].y;
                walker_state[i].z=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].z;
                walker_state[i].x0=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].x0;
                walker_state[i].y0=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].y0;
                walker_state[i].z0=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].z0;
                walker_state[i].ix=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].ix;
                walker_state[i].iy=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].iy;
                walker_state[i].iz=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].iz;
                walker_state[i].ix0=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].ix0;
                walker_state[i].iy0=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].iy0;
                walker_state[i].iz0=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].iz0;
                walker_state[i].time=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].time;
                walker_state[i].p_x=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].p_x;
                walker_state[i].p_y=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].p_y;
                walker_state[i].p_z=walker_state[n_fluorophores_local*(i/n_fluorophores_local)].p_z;

                propagate_photophysics(i);
                //std::cout<<"copy "<<i<<" from "<<n_fluorophores_local*(i/n_fluorophores_local)<<std::endl;
            }
        }
    }
    store_step_protocol();

    if (!heating_up && (save_msd_every_n_timesteps>0) && (sim_time>duration-sim_timestep) && msd_count && msd && msd2 && msd_size>0) {
        // store MSD curve
        char fn[1024];
        sprintf(fn, "%s%smsd.dat", basename.c_str(), object_name.c_str());
        std::cout<<"writing MSD file '"<<fn<<"' ... msd_size="<<msd_size<<" ... ";
        FILE* f=fopen(fn, "w");
        for (int i=0; i<msd_size; i++) {
            double tau=(double)(i)*(double)sim_timestep*(double)save_msd_every_n_timesteps;
            double cnt=(double)msd_count[i];
            double mean=msd[i]/cnt;
            double stddev=sqrt(msd2[i]/cnt-mean*mean);
            if (msd_count[i]>0) fprintf(f, "%lg, %lg, %lg, %lg\n", tau, mean, stddev, cnt);
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
        for (int plt=0; plt<2; plt++) {
            if (plt==0) {
                fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"msd.pdf").c_str());
            } else if (plt==1) {
                fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
                fprintf(f, "set output\n");
            }
            fprintf(f, "set title \"sigma_translation^2 against time\"\n");
            fprintf(f, "set xlabel \"time t [s]\"\n");
            fprintf(f, "set ylabel \"sigma^2 [microns^2]\"\n");
            fprintf(f, "msd(tau)=6 * %lf * tau\n", diff_coeff);
            fprintf(f, "set style fill transparent solid 0.5 noborder\n");
            fprintf(f, "plot \"%s\" using 1:(($2)-($3)):(($2)+($3)) notitle with filledcurves, \"%s\" using 1:2 title \"simulation result\" with points, msd(x) title \"theory\" with lines, msd_fit(x) title sprintf(\"fit D=%%.2f\",Diff) with lines\n", extract_file_name(basename+object_name+"msd.dat").c_str(), extract_file_name(basename+object_name+"msd.dat").c_str());
            fprintf(f, "\n");
            if (plt==1) fprintf(f, "pause -1\n");
            fprintf(f, "set logscale xy\n");
            fprintf(f, "set title \"sigma_translation^2 against time\"\n");
            fprintf(f, "set xlabel \"time t [s]\"\n");
            fprintf(f, "set ylabel \"sigma^2 [microns^2]\"\n");
            fprintf(f, "msd(tau)=6 * %lf * tau\n", diff_coeff);
            fprintf(f, "set style fill transparent solid 0.5 noborder\n");
            fprintf(f, "plot \"%s\" using 1:(($2)-($3)):(($2)+($3)) notitle with filledcurves, \"%s\" using 1:2 title \"simulation result\" with points, msd(x) title \"theory\" with lines, msd_fit(x) title sprintf(\"fit D=%%.2f\",Diff) with lines\n", extract_file_name(basename+object_name+"msd.dat").c_str(), extract_file_name(basename+object_name+"msd.dat").c_str());
            fprintf(f, "\n");
            if (plt==1) fprintf(f, "pause -1\n");
            fprintf(f, "unset logscale xy\n");
            fprintf(f, "set title \"particle averaged over against time\"\n");
            fprintf(f, "set xlabel \"time t [s]\"\n");
            fprintf(f, "set ylabel \"particles\"\n");
            fprintf(f, "plot \"%s\" using 1:4 title \"simulation result\" with linespoints\n", extract_file_name(basename+object_name+"msd.dat").c_str());
            fprintf(f, "\n");
            if (plt==1) fprintf(f, "pause -1\n");
        }
        fclose(f);
        std::cout<<"DONE!\n";
    }
}

std::string BrownianDynamics::report() {
    std::string s=FluorophorDynamics::report();
    s+="n_fluorophores_local = "+inttostr(n_fluorophores_local)+" fluorophore(s) per particle/walker\n";
    s+="diff_coeff = "+floattostr(diff_coeff)+" micron^2/s\n";
    s+="sigma_jump = MSD(sim_timestep) = "+floattostr(sigma_jump)+" microns\n";
    s+="MSD(0.1 msec) = "+floattostr(sqrt(6.0*diff_coeff*0.1e-3))+" microns\n";
    s+="MSD(1 msec)   = "+floattostr(sqrt(6.0*diff_coeff*1e-3))+" microns\n";
    s+="MSD(5 msec)   = "+floattostr(sqrt(6.0*diff_coeff*5e-3))+" microns\n";
    s+="MSD(10 msec)  = "+floattostr(sqrt(6.0*diff_coeff*10e-3))+" microns\n";
    s+="velocity = ( "+floattostr(velocity[0])+", "+floattostr(velocity[1])+", "+floattostr(velocity[2])+" ) micron/s\n";
    if (use_rotational_diffusion) {
        s+="diff_rot = "+floattostr(diff_rot/M_PI*180.0/M_PI*180.0)+" deg^2/s\n";
        s+="sigma_rotjump = "+floattostr(sigma_rotjump/M_PI*180.0)+" deg\n";
    } else {
        s+="no rotational diffusion\n";
    }
    if (save_msd_every_n_timesteps>0) s+="saving MSD every "+inttostr(save_msd_every_n_timesteps)+" timesteps with "+inttostr(msd_size)+" items\n";

    return s;
}

std::string BrownianDynamics::dot_get_properties() {
    std::string s=FluorophorDynamics::dot_get_properties();
    s+="n_fluorophores_local = "+inttostr(n_fluorophores_local)+" fl./walker<BR/>";
    s+="diff_coeff = "+floattostr(diff_coeff)+" &mu;m&sup2;/s<BR/>";
    s+="velocity = ( "+floattostr(velocity[0])+", "+floattostr(velocity[1])+", "+floattostr(velocity[2])+" ) &mu;m/s<BR/>";
    if (use_rotational_diffusion) {
        s+="diff_rot = "+floattostr(diff_rot/M_PI*180.0/M_PI*180.0)+" deg^2/s<BR/>";
    } else {
        s+="no rotation<BR/>";
    }
    s+="<BR/>";
    return s;
}

void BrownianDynamics::test(unsigned int steps, unsigned int walkers) {
    FluorophorDynamics::test(steps, walkers);
    std::cout<<"TEST: init test ... ";
    unsigned int nw=walkers;
    use_rotational_diffusion=true;
    init();
    change_walker_count(nw, 1);
    std::cout<<"TEST: heating up ... \n";
    run_heatup();
    std::cout<<"TEST: setup test ... \n";
    for (unsigned int i=0; i<nw; i++) {
        walker_state[i].x=0;
        walker_state[i].y=0;
        walker_state[i].z=0;
        walker_state[i].p_x=0;
        walker_state[i].p_y=0;
        walker_state[i].p_z=1;
    }
    //sim_time=0;
    std::cout<<"TEST: opening output file ... \n";
    FILE* f=fopen((basename+object_name+"test_x2.dat").c_str(), "w");
    FILE* f1=fopen((basename+object_name+"test_traj.dat").c_str(), "w");
    int nhist=100;
    gsl_histogram * h_theta = gsl_histogram_alloc (nhist+1);
    gsl_histogram_set_ranges_uniform (h_theta, 0, 180);
    gsl_histogram * h_phi = gsl_histogram_alloc(nhist+1);
    gsl_histogram_set_ranges_uniform(h_phi, 0, 360);
    std::cout<<"TEST: starting test propagation ";
    for (unsigned int i=0; i<steps; i++) {
        propagate(false);
        double xs=0;
        double dalpha=0;
        for (unsigned int w=0; w<nw; w++) {
            double theta=(M_PI/2.0-atan(walker_state[w].p_z/sqrt(gsl_pow_2(walker_state[w].p_x)+gsl_pow_2(walker_state[w].p_y))));
            double phi=acos(walker_state[w].p_x/sqrt(gsl_pow_2(walker_state[w].p_x)+gsl_pow_2(walker_state[w].p_y)));
            if (walker_state[w].p_y<0) {
                phi=2*M_PI-phi;
            }
            phi=phi/M_PI*180.0;
            theta=theta/M_PI*180.0;
            gsl_histogram_increment(h_theta, theta);
            gsl_histogram_increment(h_phi, phi);
            xs=xs+gsl_pow_2(walker_state[w].x)+gsl_pow_2(walker_state[w].y)+gsl_pow_2(walker_state[w].z);
            dalpha=dalpha+gsl_pow_2(acos(walker_state[w].p_z));
        }
        fprintf(f, "%20.10lf,      %20.10lf,  %20.10lf,      %20.10lf,  %20.10lf\n", (double)(i+1)*sim_timestep*1e6, xs/(double)nw, 6.0*diff_coeff*(double)((double)(i+1)*sim_timestep), dalpha/(double)nw/M_PI/M_PI*180.0*180.0, 2.0*diff_rot*(double)((double)(i+1)*sim_timestep)/M_PI/M_PI*180.0*180.0);
        fprintf(f1, "%20.10lf,   %20.10lf,  %20.10lf,     %20.10lf,   %20.10lf,  %20.10lf\n", walker_state[0].x, walker_state[0].y, walker_state[0].z, walker_state[0].p_x, walker_state[0].p_y, walker_state[0].p_z);
        if (i%100==0) std::cout<<".";
    }
    std::cout<<" ready!"<<std::endl;
    std::cout<<"TEST: writing test results ... "<<std::endl;
    fclose(f);
    fclose(f1);
    FILE* f2=fopen((basename+object_name+"test_hist_rot.dat").c_str(), "w");
    for (int i=0; i<nhist; i++) {
        double mi, ma, phi, theta;
        gsl_histogram_get_range(h_phi, i, &mi, &ma);
        phi=(ma+mi)/2.0;
        gsl_histogram_get_range(h_theta, i, &mi, &ma);
        theta=(ma+mi)/2.0;
        fprintf(f2, "%20.10lf, %20.10lf,    %20.10lf, %20.10lf\n", theta, gsl_histogram_get(h_theta, i), phi, gsl_histogram_get(h_phi, i));
    }
    fclose(f2);
    f=fopen((basename+object_name+"test_plot.plt").c_str(), "w");
    fprintf(f, "unset multiplot\n");
    fprintf(f, "reset\n");
    for (int plt=0; plt<2; plt++) {
        if (plt==0) {
            fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
            fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"test_plot.pdf").c_str());
        } else if (plt==1) {
            fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
            fprintf(f, "set output\n");
        }
        fprintf(f, "set title \"3D Translation Trajectory\"\n");
        fprintf(f, "set xlabel \"x\"\n");
        fprintf(f, "set ylabel \"y\"\n");
        fprintf(f, "set zlabel \"z\"\n");
        fprintf(f, "splot \"%s\" using 1:2:3 notitle with lines\n", extract_file_name(basename+object_name+"test_traj.dat").c_str());
        if (plt==1) fprintf(f, "pause -1\n");
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
        fprintf(f, "      \"%s\" using 4:5:6 notitle with lines\n", extract_file_name(basename+object_name+"test_traj.dat").c_str());
        if (plt==1) fprintf(f, "pause -1\n");
        fprintf(f, "unset parametric\n");
        fprintf(f, "set title \"sigma_translation^2 against time\"\n");
        fprintf(f, "set xlabel \"time t\"\n");
        fprintf(f, "set ylabel \"sigma^2\"\n");
        fprintf(f, "plot \"%s\" using 1:2 title \"simulation result\" with points, \"%s\" using 1:3 title \"theory\" with lines\n", extract_file_name(basename+object_name+"test_x2.dat").c_str(), extract_file_name(basename+object_name+"test_x2.dat").c_str());
        fprintf(f, "\n");
        if (plt==1) fprintf(f, "pause -1\n");
        fprintf(f, "set title \"sigma_rotation^2 against time\"\n");
        fprintf(f, "set xlabel \"time t\"\n");
        fprintf(f, "set ylabel \"sigma_alpha^2\"\n");
        fprintf(f, "plot \"%s\" using 1:4 title \"simulation result\" with points, \"%s\" using 1:5 title \"theory\" with lines, 90*90 notitle\n", extract_file_name(basename+object_name+"test_x2.dat").c_str(), extract_file_name(basename+object_name+"test_x2.dat").c_str());
        fprintf(f, "\n");
        if (plt==1) fprintf(f, "pause -1\n");
        fprintf(f, "set multiplot layout 2,1\n");
        fprintf(f, "set title \"theta histogram for rotational diffusion\"\n");
        fprintf(f, "set xlabel \"theta [deg]\"\n");
        fprintf(f, "set ylabel \"count\"\n");
        fprintf(f, "plot \"%s\" using 1:2 title \"simulation result\" with steps\n", extract_file_name(basename+object_name+"test_hist_rot.dat").c_str());
        fprintf(f, "set title \"phi histogram for rotational diffusion\"\n");
        fprintf(f, "set xlabel \"phi [deg]\"\n");
        fprintf(f, "set ylabel \"count\"\n");
        fprintf(f, "plot \"%s\" using 3:4 title \"simulation result\" with steps\n", extract_file_name(basename+object_name+"test_hist_rot.dat").c_str());
        fprintf(f, "unset multiplot\n");
        if (plt==1) fprintf(f, "pause -1\n");
    }
    fclose(f);
    gsl_histogram_free (h_theta);
    gsl_histogram_free (h_phi);

    f=fopen((basename+object_name+"test_config.txt").c_str(), "w");
    fprintf(f, "%s\n", report().c_str());
    fclose(f);
    std::cout<<"TEST: test DONE!"<<std::endl;
}

