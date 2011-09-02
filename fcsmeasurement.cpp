#include "fcsmeasurement.h"

#include <unistd.h>

FCSMeasurement::FCSMeasurement(FluorophorManager* fluorophors, std::string objectname):
    FluorescenceMeasurement(fluorophors, objectname)
{
    save_binning=false;
    online_correlation=false;
    save_timeseries=false;

    psf_region_factor=3;

    expsf_r0=0.5;
    expsf_z0=0.5*6;
    detpsf_r0=0.53;
    detpsf_z0=0.53*6;
    I0=200/0.25e-12;

    // Rho 6G
    lambda_ex=527;

    q_det=0.1;

    img_z0=2;
    img_x0=2;
    img_y0=2;

    e_x=1;
    e_y=0;
    e_z=0;
    e_pol_fraction=1.0;

    d_x=1;
    d_y=0;
    d_z=0;

    polarised_detection=false;
    polarised_excitation=false;


    corr_taumin=sim_timestep;
    S=8;
    m=2;
    P=16;
    correlator_type=0;
    duration=1.0;
    timesteps=(unsigned long long)ceil(duration/corr_taumin);

    save_binning=false;
    save_binning_time=corr_taumin*10.0;

    correlator=NULL;
    corrjanb=NULL;
    corr=NULL;
    corr_tau=NULL;
    timeseries=NULL;
    binned_timeseries=NULL;
    slots=0;

    correlation_runtime=0;

}

FCSMeasurement::~FCSMeasurement()
{
    if (correlator!=NULL) delete correlator; correlator=NULL;
    if (corrjanb!=NULL) delete corrjanb; corrjanb=NULL;
    if (timeseries!=NULL) free(timeseries); timeseries=NULL;
    if (binned_timeseries!=NULL) free(binned_timeseries); binned_timeseries=NULL;
}

void FCSMeasurement::clear() {
    if (correlator!=NULL) delete correlator; correlator=NULL;
    if (corrjanb!=NULL) delete corrjanb; corrjanb=NULL;
    if (timeseries!=NULL) free(timeseries); timeseries=NULL;
    if (binned_timeseries!=NULL) free(binned_timeseries); binned_timeseries=NULL;
    slots=0;
    corr=NULL;
    corr_tau=NULL;
}

void FCSMeasurement::read_config_internal(jkINIParser2& parser) {
    FluorescenceMeasurement::read_config_internal(parser);
    online_correlation=parser.getSetAsBool("online_correlation", online_correlation);

    psf_region_factor=parser.getAsDouble("psf_region_factor", psf_region_factor);

    expsf_r0=parser.getSetAsDouble("expsf_r0", expsf_r0);
    expsf_z0=parser.getSetAsDouble("expsf_z0", expsf_z0);

    detpsf_r0=parser.getSetAsDouble("detpsf_r0", detpsf_r0);
    detpsf_z0=parser.getSetAsDouble("detpsf_z0", detpsf_z0);

    corr_taumin=parser.getSetAsDouble("corr_taumin", corr_taumin);
    //std::cout<<"  >>>> read corr_taumin = "<<corr_taumin<<std::endl;
    S=parser.getSetAsInt("corr_S", S);
    m=parser.getSetAsInt("corr_m", m);
    P=parser.getSetAsInt("corr_P", P);
    correlator_type=parser.getSetAsInt("corr_type", correlator_type);

    if (parser.exists("I0")) {
        I0=parser.getSetAsDouble("I0", I0);
    } else {
        I0=parser.getSetAsDouble("P0", I0*(M_PI*gsl_pow_2(2.0*1e-6*expsf_r0)))/(M_PI*gsl_pow_2(2.0*1e-6*expsf_r0));
    }
    lambda_ex=parser.getSetAsDouble("lambda_ex", lambda_ex);
    q_det=parser.getSetAsDouble("q_det", q_det);

    img_z0=parser.getSetAsDouble("img_z0", img_z0);
    img_x0=parser.getSetAsDouble("img_x0", img_x0);
    img_y0=parser.getSetAsDouble("img_y0", img_y0);

    ex_z0=parser.getSetAsDouble("ex_z0", img_z0);
    ex_x0=parser.getSetAsDouble("ex_x0", img_x0);
    ex_y0=parser.getSetAsDouble("ex_y0", img_y0);

    e_x=parser.getSetAsDouble("e_x", e_x);
    e_y=parser.getSetAsDouble("e_y", e_y);
    e_z=parser.getSetAsDouble("e_z", e_z);
    e_pol_fraction=parser.getSetAsDouble("e_pol_fraction", e_pol_fraction);

    double le=sqrt(e_x*e_x+e_y*e_y+e_z*e_z);
    e_x/=le; e_y/=le; e_z/=le;

    d_x=parser.getSetAsDouble("d_x", d_x);
    d_y=parser.getSetAsDouble("d_y", d_y);
    d_z=parser.getSetAsDouble("d_z", d_z);

    le=sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
    d_x/=le; d_y/=le; d_z/=le;

    polarised_detection=parser.getSetAsBool("polarised_detection", polarised_detection);
    polarised_excitation=parser.getSetAsBool("polarised_excitation", polarised_excitation);

    save_binning=parser.getSetAsBool("save_binning", save_binning);
    save_binning_time=parser.getSetAsDouble("save_binning_time", save_binning_time);
    save_timeseries=parser.getSetAsBool("save_timeseries", save_timeseries);


    init();
}


void FCSMeasurement::init(){
    FluorescenceMeasurement::init();
    correlation_runtime=0;

    clear();
    timesteps=(unsigned long long)ceil(duration/sim_timestep);

    // when in online_correlation mode we do NOT build up a timeseries.
    // if the user selects to get a binned_timeseries we create a biined_timeseries
    // array. In non-online_correlation mode the binned series may be calculated from
    // the raw timeseries, so there is no need for a binned_timeseries array!
    if (!online_correlation)  timeseries=(uint16_t*)calloc(timesteps, sizeof(uint16_t));
    if (save_binning && online_correlation) {
        int b=round(save_binning_time/corr_taumin);
        binned_timeseries=(uint32_t*)calloc((timesteps)/b+10, sizeof(uint32_t));
    }
    slots=0;
    corr=NULL;
    corr_tau=NULL;

    //unsigned int correlation_slots=P+(S-1)*(P-P/m);
    //corr=(double*)calloc(correlation_slots, sizeof(double));
    //corr_tau=(double*)calloc(correlation_slots, sizeof(double));
    correlator=new MultiTauCorrelator<double, double>(S, m, P, corr_taumin);
    corrjanb=new correlatorjb<double, double>(S, P);
    Ephoton=6.626e-34*2.99e8/(lambda_ex*1e-9);

    if (corr_taumin<0.99*sim_timestep) {
        throw MeasurementException("corr_taumin may not be smaller than simulation timestep\ncorr_taumin="+floattostr(corr_taumin)+"   sim_timestep="+floattostr(sim_timestep));
    }

    double r=corr_taumin/sim_timestep;
    if (fabs(round(r)-r)>0.0000001) {
        throw MeasurementException("corr_taumin has to be an integer multiple of sim_timestep!!!\ncorr_taumin="+floattostr(corr_taumin)+"   sim_timestep="+floattostr(sim_timestep)+"   r="+floattostr(r));
    } else corr_taumin=round(r)*sim_timestep;

    for (size_t i=0; i<dyn.size(); i++) {
        dyn[i]->load_all_used_spectra();
    }
    endCurrentStep=sim_time+corr_taumin;
    nphot_sum=0;
    current_timestep=0;
    display_temp=0;
    bin_sum=0;
    bin_counter=0;
    bin_i=0;
}

void FCSMeasurement::propagate(){
    tick();
    FluorescenceMeasurement::propagate();
    run_fcs_simulation();
    if (sim_time>=duration) {
        // calculate correlation function
        printf("correlating ...\n");
        PublicTickTock tim;
        tim.tick();
        long* taus=(long*)malloc(S*P*sizeof(long));
        if (!online_correlation) {
            if (correlator_type==0) {
                for (int i=0; i<timesteps; i++) {
                    correlator->correlate_step(timeseries[i]);
                }
            } else if (correlator_type==1) {
                for (int i=0; i<timesteps; i++) {
                    corrjanb->run(timeseries[i], timeseries[i]);
                }
            } else if (correlator_type==2) {
                statisticsAutocorrelateCreateMultiTau(taus, S, m, P);
                corr=(double*)malloc(S*P*sizeof(double));
                statisticsAutocorrelateMultiTauSymmetric(corr, timeseries, timesteps, taus, S*P);
            }
        }
        if (correlator_type==0) {
            correlator->normalize();
            corr=correlator->getCor();
            corr_tau=correlator->getCorTau();
            slots=correlator->getSlots();
        } else if (correlator_type==1) {
            slots=S*P;
            double** corr1=corrjanb->get_array_G();
            corr=(double*)malloc(slots*sizeof(double));
            corr_tau=(double*)malloc(slots*sizeof(double));
            for (int i=0; i<slots; i++) {
                corr_tau[i]=corr1[0][i]*corr_taumin;
                corr[i]=corr1[1][i];
            }
        } else if (correlator_type==2) {
            slots=S*P;
            corr_tau=(double*)malloc(slots*sizeof(double));
            for (int i=0; i<slots; i++) {
                corr_tau[i]=(double)taus[i]*corr_taumin;
            }
        }

        free(taus);

        tim.tock();
        correlation_runtime+=tim.get_duration();
        printf(" done!\n");
    }
    tock();
    runtime=runtime+get_duration();
}



void FCSMeasurement::run_fcs_simulation(){
    // number of fluorescence photons per molecule and sim_timestep, but you have to multiply by sigma_abs and q_fluor (walker dependent!)
    double n00=I0*1e-6/Ephoton*sim_timestep;
    int bin_r=round(save_binning_time/corr_taumin); // binning ration

    // first we go through alle the fluorophors and integrate their contribution
    for (size_t d=0; d<dyn.size(); d++) { // go through all dynamics objects that provide data for this measurement object
        FluorophorDynamics::walkerState* dynn=dyn[d]->get_walker_state();
        unsigned long wc=dyn[d]->get_walker_count();
        if (!dyn[d]->end_of_trajectory_reached()) for (unsigned long i=0; i<wc; i++) { // iterate through all walkers in the d-th dynamics object
            if (dynn[i].exists) {
                double x0,y0,z0;
                //dynn->get_walker_position(i, &x0, &y0, &z0);
                x0=dynn[i].x;
                y0=dynn[i].y;
                z0=dynn[i].z;
                double dx=x0-img_x0;
                double dz=z0-img_z0;
                double dy=y0-img_y0;
                double edx=x0-ex_x0;
                double edz=z0-ex_z0;
                double edy=y0-ex_y0;
                double edxs=edx*edx;
                double edys=edy*edy;
                if ((fabs(edz)<(psf_region_factor*detpsf_z0)) && ((edxs+edys)<gsl_pow_2(psf_region_factor*detpsf_r0))) {
                    double dxs=dx*dx;
                    double dys=dy*dy;
                    double n0=n00*dyn[d]->get_walker_sigma_times_qfl(i);
                    n0=n0*fluorophors->get_spectral_absorbance(dynn[i].spectrum, lambda_ex);
                    //dynn->get_walker_spectral_absorbance(i, lambda_ex);
                    double dpx, dpy, dpz;
                    dpx=dynn[i].p_x;
                    dpy=dynn[i].p_y;
                    dpz=dynn[i].p_z;
                    if (polarised_excitation) {
                        n0=n0*(1.0-e_pol_fraction+e_pol_fraction*gsl_pow_2(dpx*e_x+dpy*e_y+dpz*e_z));
                    }
                    n0=n0*exp(-2.0*(edxs+edys)/gsl_pow_2(expsf_r0)-2.0*gsl_pow_2(edz)/gsl_pow_2(expsf_z0));
                    if (polarised_detection) {
                        n0=n0*gsl_pow_2(dpx*d_x+dpy*d_y+dpz*d_z);
                    }
                    nphot_sum=nphot_sum+n0*q_det*exp(-2.0*(dxs+dys)/gsl_pow_2(detpsf_r0)-2.0*gsl_pow_2(dz)/gsl_pow_2(detpsf_z0));
                    //std::cout<<"nphot_sum="<<nphot_sum<<"\n";
                }
            }
        }
    }

    // if the current integration step ended, we may add a new value to the correlator
    //std::cout<<"sim_time="<<sim_time<<"   endCurrentStep="<<endCurrentStep<<"\n";
    if (sim_time>=endCurrentStep) {
        endCurrentStep=sim_time+corr_taumin;
        register uint16_t N=gsl_ran_poisson(rng, nphot_sum);
        //std::cout<<"N="<<N<<std::endl;
        if (online_correlation) {
            if (correlator_type==0) {
                correlator->correlate_step(N);
            } else if (correlator_type==1) {
                corrjanb->run(N,N);
            }
        } else {
            timeseries[current_timestep]=N;
        }
        if (save_binning && online_correlation) {
            bin_sum=bin_sum+N;
            bin_counter++;
            if (bin_counter==bin_r) {
                bin_counter=0;
                binned_timeseries[bin_i]=bin_sum;
                bin_sum=0;
                bin_i++;
            }
        }
        if (current_timestep%(timesteps/1000)==0) {
            std::cout<<format("%4.1lf", (double)current_timestep/(timesteps)*100.0)<<"%:   "<<display_temp<<std::endl;
            display_temp=0;
        } else {
            display_temp+=N;
        }
        current_timestep++;
        nphot_sum=0;
    }
}


void FCSMeasurement::save() {
    FILE* f;
    char fn[255];


    sprintf(fn, "%s%scorr.dat", basename.c_str(), object_name.c_str());
    printf("writing '%s' ...", fn);
    f=fopen(fn, "w");
    unsigned long long istart=1;
    if (correlator_type==2) istart=0;

    for (unsigned long long i=istart; i<slots; i++) {
        if (correlator_type==1) fprintf(f, "%15.10lf, %15.10lf\n", corr_tau[i], corr[i]);
        else fprintf(f, "%15.10lf, %15.10lf\n", corr_tau[i], corr[i]);
    }
    fclose(f);
    printf(" done!\n");
    std::string corrfn=fn;

    sprintf(fn, "%s%scorrplot.plt", basename.c_str(), object_name.c_str());
    printf("writing '%s' ...", fn);
    f=fopen(fn, "w");
    double psf_r0=1.0/sqrt(1.0/detpsf_r0/detpsf_r0+1.0/expsf_r0/expsf_r0);
    fprintf(f, "g(t)=1.0+1.0/N/(1.0+t/tauD)/sqrt(1.0+t/gamma/gamma/tauD)\n");
    fprintf(f, "N=1\n");
    fprintf(f, "tauD=100e-6\n");
    fprintf(f, "gamma=1\n");
    fprintf(f, "wxy=%lf\n", psf_r0);
    fprintf(f, "fit g(x) \"%s\" via N, tauD,gamma\n", extract_file_name(corrfn).c_str());

    fprintf(f, "set logscale x\n");
    fprintf(f, "plot \"%s\" title \"simulation data\" with points, g(x) title sprintf(\"fit N=%%f, tauD=%%f microSecs, gamma=%%f, D=%%f micron^2/s\",N, tauD*1e6, gamma, wxy*wxy/4.0/tauD)\n", extract_file_name(corrfn).c_str());
    fprintf(f, "pause -1\n");
    fclose(f);
    printf(" done!\n");

    if (!online_correlation) {
        if (save_timeseries) {
            sprintf(fn, "%s%sts.dat", basename.c_str(), object_name.c_str());
            printf("writing '%s' ...", fn);
            f=fopen(fn, "w");
            double t=0;
            for (unsigned long long i=0; i<timesteps; i++) {
                fprintf(f, "%15.10lf, %d\n", t, timeseries[i]);
                t=t+corr_taumin;
            }
            fclose(f);
            std::string tsfn=fn;
            printf(" done!\n");

            sprintf(fn, "%s%stsplot.plt", basename.c_str(), object_name.c_str());
            printf("writing '%s' ...", fn);
            f=fopen(fn, "w");
            fprintf(f, "set xlabel \"time [seconds]\"\n");
            fprintf(f, "set ylabel \"photon count [photons/%lfsec]\"\n", corr_taumin);
            fprintf(f, "plot \"%s\" with steps\n", extract_file_name(tsfn).c_str());
            fprintf(f, "pause -1\n");
            fprintf(f, "set xlabel \"time [seconds]\"\n");
            fprintf(f, "set ylabel \"photon count [Hz]\"\n");
            fprintf(f, "plot \"%s\" using 1:(($2)/%lf) with steps\n", extract_file_name(tsfn).c_str(), corr_taumin);
            fprintf(f, "pause -1\n");
            fclose(f);
            printf(" done!\n");
        }
    }

    if (save_binning) {
        if (online_correlation) {
            // in online-correlation-mode no complete timeseries is built up, only a binned one, so we save that
            sprintf(fn, "%s%sbts.dat", basename.c_str(), object_name.c_str());
            printf("writing '%s' ...", fn);
            f=fopen(fn, "w");
            double t=0;
            int b=round(save_binning_time/corr_taumin);
            for (unsigned long long i=0; i<timesteps/b-1; i++) {
                fprintf(f, "%15.10lf, %u\n", t, binned_timeseries[i]);
                t=t+corr_taumin*(double)b;
            }
            fclose(f);
            std::string tsfn=fn;
            printf(" done!\n");

            sprintf(fn, "%s%sbtsplot.plt", basename.c_str(), object_name.c_str());
            printf("writing '%s' ...", fn);
            f=fopen(fn, "w");
            fprintf(f, "set xlabel \"time [seconds]\"\n");
            fprintf(f, "set ylabel \"photon count [photons/%lfsec]\"\n", corr_taumin*b);
            fprintf(f, "plot \"%s\" with steps\n", extract_file_name(tsfn).c_str());
            fprintf(f, "pause -1\n");
            fprintf(f, "set xlabel \"time [seconds]\"\n");
            fprintf(f, "set ylabel \"photon count [Hz]\"\n");
            fprintf(f, "plot \"%s\" using 1:(($2)/%lf) with steps\n", extract_file_name(tsfn).c_str(), corr_taumin*b);
            fprintf(f, "pause -1\n");
            fclose(f);
            printf(" done!\n");
        } else {
            // if we are not in online-correlation mode we may simply calculate the binned timeseries from the
            // complete timeseries
            sprintf(fn, "%s%sbts.dat", basename.c_str(), object_name.c_str());
            printf("writing '%s' ...", fn);
            f=fopen(fn, "w");
            double t=0;
            int b=round(save_binning_time/corr_taumin);
            for (unsigned long long i=0; i<timesteps-b; i=i+b) {
                register uint32_t ts=0;
                for (int j=0; j<b; j++) {
                    ts=ts+timeseries[i+j];
                }
                fprintf(f, "%15.10lf, %lu\n", t, ts);
                t=t+corr_taumin;
            }
            fclose(f);
            std::string tsfn=fn;
            printf(" done!\n");

            sprintf(fn, "%s%sbtsplot.plt", basename.c_str(), object_name.c_str());
            printf("writing '%s' ...", fn);
            f=fopen(fn, "w");
            fprintf(f, "set xlabel \"time [seconds]\"\n");
            fprintf(f, "set ylabel \"photon count [photons/%lfsec]\"\n", corr_taumin*b);
            fprintf(f, "plot \"%s\" with steps\n", extract_file_name(tsfn).c_str());
            fprintf(f, "pause -1\n");
            fprintf(f, "set xlabel \"time [seconds]\"\n");
            fprintf(f, "set ylabel \"photon count [Hz]\"\n");
            fprintf(f, "plot \"%s\" using 1:(($2)/%lf) with steps\n", extract_file_name(tsfn).c_str(), corr_taumin*b);
            fprintf(f, "pause -1\n");
            fclose(f);
            printf(" done!\n");
        }
    }



}

std::string FCSMeasurement::report(){
    std::string s=FluorescenceMeasurement::report();
    s+="pos_laser     = ["+floattostr(ex_x0)+", "+floattostr(ex_y0)+", "+floattostr(ex_z0)+"] um\n";
    s+="pos_detection = ["+floattostr(img_x0)+", "+floattostr(img_y0)+", "+floattostr(img_z0)+"] um\n";
    s+="distance_laser_detector = ["+floattostr(img_x0-ex_x0)+", "+floattostr(img_y0-ex_y0)+", "+floattostr(img_z0-ex_z0)+"] um\n";
    s+="                        = "+floattostr(sqrt(gsl_pow_2(img_x0-ex_x0)+gsl_pow_2(img_y0-ex_y0)+gsl_pow_2(img_z0-ex_z0)))+" um\n";
    s+="expsf_r0 = "+floattostr(expsf_r0)+" microns\n";
    s+="expsf_z0 = "+floattostr(expsf_z0)+" microns\n";
    s+="detpsf_r0 = "+floattostr(detpsf_r0)+" microns\n";
    s+="detpsf_z0 = "+floattostr(detpsf_z0)+" microns\n";
    double psf_r0=1.0/sqrt(1.0/detpsf_r0/detpsf_r0+1.0/expsf_r0/expsf_r0);
    s+="psf_system = sqrt(1/expsf_r0^2 + 1/detpsf_r0^2) = "+floattostr(psf_r0)+" microns\n";

    double Veff=pow(M_PI, 1.5)*(psf_r0*psf_r0*psf_r0*detpsf_z0/detpsf_r0);
    s+="focus_volume (Veff) = pi^(3/2) * psf_system^3 * detpsf_z0/detpsf_r0 = "+floattostr(Veff)+" femto litre\n";
    double sum=0;
    if (dyn.size()>0) {
        for (size_t i=0; i<dyn.size(); i++) {
            s+="<"+dyn[i]->get_object_name()+" particles in Veff> = "+floattostr(Veff*dyn[i]->get_c_fluor()*6.022e23*1e-9/1e15)+"\n";
            sum=sum+Veff*dyn[i]->get_c_fluor()*6.022e23*1e-9/1e15;
        }
        s+="  SUM <particles in Veff> = "+floattostr(sum)+"\n";
    }
    s+="psf_region_factor = "+floattostr(psf_region_factor)+"\n";
    if (polarised_excitation) {
        s+="polarised excitation: pe = ["+floattostr(e_x)+", "+floattostr(e_y)+", "+floattostr(e_z)+"]\n";
        s+="polarised excitation: e_pol_fraction ="+floattostr(e_pol_fraction*100.0)+"% \n";
    } else {
        s+="non-polarised excitation\n";
    }
    if (polarised_detection) {
        s+="polarised detection: pe = ["+floattostr(d_x)+", "+floattostr(d_y)+", "+floattostr(d_z)+"]\n";
    } else {
        s+="non-polarised detection\n";
    }
    s+="lambda_ex = "+floattostr(lambda_ex)+" nm\n";
    s+="EPhoton_ex = "+floattostr(Ephoton/1.602176487e-19/1e-3)+" meV\n";
    s+="I0 = "+floattostr(I0)+" uW/m^2  =  "+floattostr(I0/1e12)+" uW/micron^2  =  "+floattostr(I0/1e6)+" uW/mm^2\n";
    s+="P0 [on focus, i.e. on A=pi*(2*expsf_r0)^2] = "+floattostr(I0*(M_PI*gsl_pow_2(2.0*expsf_r0*1e-6)))+" uW\n";
    for (size_t i=0; i<dyn.size(); i++) {
        s+="  using "+dyn[i]->get_object_name()+" data:\n";
        s+="    init_spectrum = "+floattostr(dyn[i]->get_init_spectrum())+"\n";
        s+="    sigma_abs = "+format("%g", dyn[i]->get_init_sigma_abs(0),-1, false, 1e-30)+" m^2\n";
        s+="    q_fluor = "+floattostr(dyn[i]->get_init_q_fluor(0)*100)+" %\n";
        s+="    abs_photons/(molecule*s)  = "+floattostr(I0*1e-6/Ephoton*dyn[i]->get_init_sigma_abs(0))+"\n";
        s+="    fluor_photons/(molecule*s) = "+floattostr(dyn[i]->get_init_q_fluor(0)*I0*1e-6/Ephoton*dyn[i]->get_init_sigma_abs(0))+"\n";
        s+="    det_photons/(molecule*s) = "+floattostr(q_det*dyn[i]->get_init_q_fluor(0)*I0*1e-6/Ephoton*dyn[i]->get_init_sigma_abs(0))+"\n";
        s+="    absorbtion @ lambda_ex = "+floattostr(fluorophors->get_spectral_absorbance(dyn[i]->get_init_spectrum(), lambda_ex)*100.0)+" %\n";
    }
    s+="duration = "+floattostr(duration*1e3)+" msecs\n";
    s+="=> timesteps = "+inttostr(timesteps)+"     à  timestep-duration = "+floattostr(sim_timestep*1e9)+" nsecs\n";
    s+="correlator:  S(#corr)="+inttostr(S)+"   P(#chan/dec)="+inttostr(P)+"    m(binRatio)="+inttostr(m)+"\n";
    s+="             corr_tau_min= "+floattostr_fmt(corr_taumin*1e9, "%lg")+" ns = "+floattostr(corr_taumin/sim_timestep)+" simulation timesteps\n";
    s+="             correlator_type= "+inttostr(correlator_type);
    if (correlator_type==0) s+=" (corr_jk)";
    if (correlator_type==1) s+=" (corr_jb)";
    if (correlator_type==2) s+=" (corr_direct)";
    s+="\n";

    double mem=0;
    if (!online_correlation)  mem=mem+timesteps*sizeof(uint16_t);
    if (save_binning && online_correlation) {
        int b=round(save_binning_time/corr_taumin);
        mem=mem+(timesteps+1)/b*sizeof(uint32_t);
    }

    mem=mem+correlator->calc_mem_consumption();
    s+="estimated_memory_consumtion_detection = "+floattostr(mem/1024.0/1024.0)+" MBytes\n";
    sum=0;
    for (size_t i=0; i<dyn.size(); i++) {
        s+="estimated_memory_consumtion_dynamics("+dyn[i]->get_object_name()+") = "+floattostr(dyn[i]->calc_mem_consumption()/1024.0/1024.0)+" MBytes\n";
        sum=sum+dyn[i]->calc_mem_consumption();
    }
    s+="estimated_memory_consumtion_all = "+floattostr((mem+sum)/1024.0/1024.0)+" MBytes\n";
    s+="runtime correlation = "+floattostr(correlation_runtime)+" secs\n";

    return s;
}
