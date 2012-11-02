#include "fcsmeasurement.h"

#include <unistd.h>

FCSMeasurement::FCSMeasurement(FluorophorManager* fluorophors, std::string objectname):
    FluorescenceMeasurement(fluorophors, objectname)
{
    save_binning=false;
    online_correlation=false;
    save_timeseries=false;

    psf_region_factor=3;
    det_wavelength_max=-1;
    det_wavelength_min=-1;
    background_rate=0;
    offset_rate=0;
    offset_std=0;
    offset_correction=0;

    ill_distribution=0;
    det_distribution=0;
    pixel_size=0.4;

    detector_type=0;
    lindet_bits=14;
    lindet_gain=10;
    lindet_var_factor=10;

    expsf_r0=0.5;
    expsf_z0=0.5*6;
    detpsf_r0=0.53;
    detpsf_z0=0.53*6;
    I0=200/0.25e-12;

    max_photons=0xFFFFFF;

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
    timeseries_size=0;
    binned_timeseries=NULL;
    binned_timeseries_size=0;
    slots=0;

    correlation_runtime=0;

}

FCSMeasurement::~FCSMeasurement()
{
    if (correlator!=NULL) delete correlator; correlator=NULL;
    if (corrjanb!=NULL) delete corrjanb; corrjanb=NULL;
    if (timeseries!=NULL) free(timeseries); timeseries=NULL;; timeseries_size=0;
    if (binned_timeseries!=NULL) free(binned_timeseries); binned_timeseries=NULL; binned_timeseries_size=0;
}

void FCSMeasurement::clear() {
    if (correlator!=NULL) delete correlator; correlator=NULL;
    if (corrjanb!=NULL) delete corrjanb; corrjanb=NULL;
    if (timeseries!=NULL) free(timeseries); timeseries=NULL; timeseries_size=0;
    if (binned_timeseries!=NULL) free(binned_timeseries); binned_timeseries=NULL; binned_timeseries_size=0;
    slots=0;
    corr=NULL;
    corr_tau=NULL;
}

void FCSMeasurement::read_config_internal(jkINIParser2& parser) {
    FluorescenceMeasurement::read_config_internal(parser);
    online_correlation=parser.getSetAsBool("online_correlation", online_correlation);

    max_photons=parser.getAsInt("max_photons", max_photons);
    det_wavelength_min=parser.getAsDouble("det_wavelength_min", det_wavelength_min);
    det_wavelength_max=parser.getAsDouble("det_wavelength_max", det_wavelength_max);

    background_rate=parser.getAsDouble("background_rate", background_rate);
    offset_rate=parser.getAsDouble("offset_rate", offset_rate);
    offset_std=parser.getAsDouble("offset_std", offset_std);
    offset_correction=parser.getAsDouble("offset_correction", offset_correction);
    psf_region_factor=parser.getAsDouble("psf_region_factor", psf_region_factor);

    ill_distribution=str_to_ill_distribution(parser.getSetAsString("ill_distribution", ill_distribution_to_str(ill_distribution)));
    det_distribution=str_to_det_distribution(parser.getSetAsString("det_distribution", ill_distribution_to_str(det_distribution)));

    pixel_size=parser.getSetAsDouble("pixel_size", pixel_size);
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

    std::string dt=tolower(parser.getAsString("detector_type", "photon_counting"));
    if (dt=="0" || dt=="photon" || dt=="photon_counting") detector_type=0;
    if (dt=="1" || dt=="linear" || dt=="emccd" || dt=="ccd" || dt=="camera") detector_type=1;
    lindet_bits=parser.getSetAsInt("lindet_bits", lindet_bits);
    lindet_gain=parser.getSetAsDouble("lindet_gain", lindet_gain);
    lindet_var_factor=parser.getSetAsDouble("lindet_var_factor", lindet_var_factor);
    init();
}


void FCSMeasurement::init(){
    FluorescenceMeasurement::init();
    correlation_runtime=0;

    clear();
    timesteps=(unsigned long long)ceil(duration/corr_taumin);

    // when in online_correlation mode we do NOT build up a timeseries.
    // if the user selects to get a binned_timeseries we create a biined_timeseries
    // array. In non-online_correlation mode the binned series may be calculated from
    // the raw timeseries, so there is no need for a binned_timeseries array!
    if (!online_correlation)  {
      timeseries=(uint32_t*)calloc((unsigned long long)ceil(duration/corr_taumin)+5, sizeof(uint32_t));
      timeseries_size=(unsigned long long)ceil(duration/corr_taumin)+5;
    }
    if (save_binning && online_correlation) {
        int b=round(save_binning_time/corr_taumin);
        binned_timeseries=(uint32_t*)calloc((timesteps)/b+10, sizeof(uint32_t));
	binned_timeseries_size=(timesteps)/b+10;
    }
    slots=0;
    corr=NULL;
    corr_tau=NULL;

    //unsigned int correlation_slots=P+(S-1)*(P-P/m);
    //corr=(double*)calloc(correlation_slots, sizeof(double));
    //corr_tau=(double*)calloc(correlation_slots, sizeof(double));
    correlator=new MultiTauCorrelator<double, double>(S, m, P, corr_taumin);
    corrjanb=new correlatorjb<double, double>(S, P, double(0.0));
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
                statisticsAutocorrelateMultiTauSymmetric<uint32_t,uint64_t>(corr, timeseries, timesteps, taus, S*P);
            } else if (correlator_type==3) {
                statisticsAutocorrelateCreateMultiTau(taus, S, m, P);
                corr=(double*)malloc(S*P*sizeof(double));
                statisticsAutocorrelateMultiTauAvgSymmetric<uint32_t,uint64_t,uint64_t>(corr, timeseries, timesteps, S, m, P, 1);
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
        } else if (correlator_type==2 || correlator_type==3) {
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
                    switch (ill_distribution) {
                        case 0:
                            n0=n0*exp(-2.0*(edxs+edys)/gsl_pow_2(expsf_r0)-2.0*gsl_pow_2(edz)/gsl_pow_2(expsf_z0));
                            break;
                        case 1:
                            n0=n0*exp(-2.0*gsl_pow_2(edz)/gsl_pow_2(expsf_z0));
                            break;
                    }
                    if (polarised_detection) {
                        n0=n0*gsl_pow_2(dpx*d_x+dpy*d_y+dpz*d_z);
                    }
                    if (det_wavelength_min>0 && det_wavelength_max>0) {
                        n0=n0*fluorophors->get_spectral_fluorescence(dynn[i].spectrum, det_wavelength_min, det_wavelength_max);
                    }
                    double sq2=sqrt(2.0);
                    switch (det_distribution) {
                        case 0:
                            nphot_sum=nphot_sum+n0*q_det*exp(-2.0*(dxs+dys)/gsl_pow_2(detpsf_r0)-2.0*gsl_pow_2(dz)/gsl_pow_2(detpsf_z0));
                            break;
                        case 1:
                            nphot_sum=nphot_sum+n0*q_det*exp(gsl_pow_2(dz)/gsl_pow_2(detpsf_z0))*(erf((pixel_size-2.0*dx)/(sq2*gsl_pow_2(detpsf_r0)))+erf((pixel_size+2.0*dy)/(gsl_pow_2(detpsf_r0)*sq2)))/(2.0*erf(pixel_size/(sq2*gsl_pow_2(detpsf_r0))));
                            break;
                    }

                    //std::cout<<"nphot_sum="<<nphot_sum<<"\n";
                }
            }
        }
    }

    // if the current integration step ended, we may add a new value to the correlator
    //std::cout<<"sim_time="<<sim_time<<"   endCurrentStep="<<endCurrentStep<<"\n";
    if (sim_time>=endCurrentStep) {
        endCurrentStep=sim_time+corr_taumin;
        register int32_t N=0;
        if (detector_type==0) N=gsl_ran_poisson(rng, nphot_sum);
        else {
            double d=gsl_ran_gaussian_ziggurat(rng, sqrt(nphot_sum*lindet_gain*lindet_var_factor))+nphot_sum*lindet_gain;
            const uint32_t Nmax=gsl_pow_int(2,lindet_bits)-1;
            if (d>Nmax) d=Nmax;
            if (d<0) d=0;
            N=floor(d);
        }
        if (background_rate>0) {
            N=N+gsl_ran_poisson(rng, background_rate*corr_taumin);
        }
        if (offset_rate>0 && offset_std>0) {
            N=N+gsl_ran_gaussian_ziggurat(rng, offset_std)+offset_rate;
        }
        N=N-offset_correction;
        if (N<0) N=0;
        N=mmin(N, max_photons);
        //std::cout<<"nphot_sum="<<nphot_sum<<"    N="<<N<<std::endl;
        if (online_correlation) {
            if (correlator_type==0) {
                correlator->correlate_step(N);
            } else if (correlator_type==1) {
                corrjanb->run(N,N);
            }
        } else {
            if (current_timestep<timeseries_size) timeseries[current_timestep]=N;
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
        if (current_timestep%mmax(1,(timesteps/1000))==0) {
            display_temp+=N;
            std::cout<<format("%4.1lf", ((double)current_timestep*corr_taumin/corr_taumin)/(timesteps)*100.0)<<"%:   "<<display_temp<<std::endl;
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
    if (correlator_type==1) istart=0;
    if (correlator_type==2) istart=0;
    if (correlator_type==3) istart=0;

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
    double psf_z0=1.0/sqrt(1.0/detpsf_z0/detpsf_z0+1.0/expsf_z0/expsf_z0);
    fprintf(f, "g(t,N,tauD,gamma,alpha)=1.0+1.0/N/(1.0+((t/tauD)**alpha))/sqrt(1.0+((t/tauD)**alpha)/gamma/gamma)\n");
    fprintf(f, "N=1\n");
    fprintf(f, "tauD=100e-6\n");
    fprintf(f, "gamma=%lf\n", psf_z0/psf_r0);
    fprintf(f, "Nf=1\n");
    fprintf(f, "tauDf=100e-6\n");
    fprintf(f, "gammaf=%lf\n", psf_z0/psf_r0);
    fprintf(f, "Na=1\n");
    fprintf(f, "tauDa=100e-6\n");
    fprintf(f, "gammaa=%lf\n", psf_z0/psf_r0);
    fprintf(f, "alphaa=1\n");
    fprintf(f, "wxy=%lf\n", psf_r0);
    fprintf(f, "fit g(x, N, tauD, gamma, 1) \"%s\" via N, tauD,gamma\n", extract_file_name(corrfn).c_str());
    fprintf(f, "fit g(x, Nf, tauDf, gammaf, 1) \"%s\" via Nf, tauDf\n", extract_file_name(corrfn).c_str());
    fprintf(f, "fit g(x, Na, tauDa, gammaa, alphaa) \"%s\" via Na, tauDa, alphaa\n", extract_file_name(corrfn).c_str());

    fprintf(f, "Veffa=pi**(3/2)*wxy*wxy*wxy*gammaa\n");
    fprintf(f, "Vefff=pi**(3/2)*wxy*wxy*wxy*gammaf\n");
    fprintf(f, "Veff=pi**(3/2)*wxy*wxy*wxy*gamma\n");
    for (int plt=0; plt<2; plt++) {
        if (plt==0) {
            fprintf(f, "set terminal pdfcairo color solid font \"Arial, 7\" linewidth 2 size 20cm,15cm\n");
            fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"corrplot.pdf").c_str());
        } else if (plt==1) {
            fprintf(f, "set terminal wxt font \"Arial, 8\"\n");
            fprintf(f, "set output\n");
        }
        fprintf(f, "set logscale x\n");
        fprintf(f, "set title \"object description: %s\"\n", description.c_str());
        fprintf(f, "plot \"%s\" title \"simulation data\" with points, "
                   "g(x,N,tauD,gamma,1) title sprintf(\"fit N=%%.3f, tauD=%%.3f uS, gamma=%%.3f, D=%%.3f um^2/s, c=%%.3f nM\",N, tauD*1e6, gamma, wxy*wxy/4.0/tauD, N/6.022e14/Veff), "
                   "g(x,Nf,tauDf,gammaf,1) title sprintf(\"fit N=%%.3f, tauD=%%.3f uS, gamma=%%.3f, D=%%.3f um^2/s, c=%%.3f nM\",Nf, tauDf*1e6, gammaf, wxy*wxy/4.0/tauDf, Nf/6.022e14/Vefff), "
                   "g(x,Na,tauDa,gammaa,alphaa) title sprintf(\"fit N=%%.3f, tauD=%%.3f uS, gamma=%%.3f, alpha=%%.3f, D=%%.3f um^2/s, c=%%.3f nM\",Na, tauDa*1e6, gammaa,alphaa,  wxy*wxy/4.0/tauDa, Na/6.022e14/Veffa)"
                   "\n", extract_file_name(corrfn).c_str());
    }
    fprintf(f, "pause -1\n");
    fclose(f);
    printf(" done!\n");

    sprintf(fn, "%s%scorrplot_simple.plt", basename.c_str(), object_name.c_str());
    printf("writing '%s' ...", fn);
    f=fopen(fn, "w");
    for (int plt=0; plt<2; plt++) {
        if (plt==0) {
            fprintf(f, "set terminal pdfcairo color solid font \"Arial, 7\" linewidth 2 size 20cm,15cm\n");
            fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"corrplot_simple.pdf").c_str());
        } else if (plt==1) {
            fprintf(f, "set terminal wxt font \"Arial, 8\"\n");
            fprintf(f, "set output\n");
        }

        fprintf(f, "set logscale x\n");
        fprintf(f, "set title \"object description: %s\"\n", description.c_str());
        fprintf(f, "plot \"%s\" title \"simulation data\" with points\n", extract_file_name(corrfn).c_str());
    }
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
            for (int plt=0; plt<2; plt++) {
                if (plt==0) {
                    fprintf(f, "set terminal pdfcairo color solid font \"Arial, 7\" linewidth 2 size 20cm,15cm\n");
                    fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"tsplot.pdf").c_str());
                } else if (plt==1) {
                    fprintf(f, "set terminal wxt font \"Arial, 8\"\n");
                    fprintf(f, "set output\n");
                }
                fprintf(f, "set xlabel \"time [seconds]\"\n");
                fprintf(f, "set ylabel \"photon count [photons/%lfsec]\"\n", corr_taumin);
                fprintf(f, "set title \"object description: %s\"\n", description.c_str());
                fprintf(f, "plot \"%s\" with steps\n", extract_file_name(tsfn).c_str());
                if (plt==1) fprintf(f, "pause -1\n");
                fprintf(f, "set xlabel \"time [seconds]\"\n");
                fprintf(f, "set ylabel \"photon count [Hz]\"\n");
                fprintf(f, "set title \"object description: %s\"\n", description.c_str());
                fprintf(f, "plot \"%s\" using 1:(($2)/%lf) with steps\n", extract_file_name(tsfn).c_str(), corr_taumin);
                if (plt==1) fprintf(f, "pause -1\n");
            }
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
                t=t+save_binning_time;
            }
            fclose(f);
            std::string tsfn=fn;
            printf(" done!\n");

            sprintf(fn, "%s%sbtsplot.plt", basename.c_str(), object_name.c_str());
            printf("writing '%s' ...", fn);
            f=fopen(fn, "w");
            for (int plt=0; plt<2; plt++) {
                if (plt==0) {
                    fprintf(f, "set terminal pdfcairo color solid font \"Arial, 7\" linewidth 2 size 20cm,15cm\n");
                    fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"btsplot.pdf").c_str());
                } else if (plt==1) {
                    fprintf(f, "set terminal wxt font \"Arial, 8\"\n");
                    fprintf(f, "set output\n");
                }
                fprintf(f, "set xlabel \"time [seconds]\"\n");
                fprintf(f, "set ylabel \"photon count [photons/%lfsec]\"\n", corr_taumin*b);
                fprintf(f, "set title \"object description: %s\"\n", description.c_str());
                fprintf(f, "plot \"%s\" with steps\n", extract_file_name(tsfn).c_str());
                if (plt==1) fprintf(f, "pause -1\n");
                fprintf(f, "set xlabel \"time [seconds]\"\n");
                fprintf(f, "set ylabel \"photon count [Hz]\"\n");
                fprintf(f, "set title \"object description: %s\"\n", description.c_str());
                fprintf(f, "plot \"%s\" using 1:(($2)/%lf) with steps\n", extract_file_name(tsfn).c_str(), corr_taumin*b);
                if (plt==1) fprintf(f, "pause -1\n");
            }
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
                register long unsigned int ts=0;
                for (int j=0; j<b; j++) {
                    ts=ts+timeseries[i+j];
                }
                fprintf(f, "%15.10lf, %lu\n", t, ts);
                //fprintf(stdout, "%15.10lf, %lu\n", t, ts);
                t=t+corr_taumin;
            }
            fclose(f);
            std::string tsfn=fn;
            printf(" done!\n");

            sprintf(fn, "%s%sbtsplot.plt", basename.c_str(), object_name.c_str());
            printf("writing '%s' ...", fn);
            f=fopen(fn, "w");
            for (int plt=0; plt<2; plt++) {
                if (plt==0) {
                    fprintf(f, "set terminal pdfcairo color solid font \"Arial, 7\" linewidth 2 size 20cm,15cm\n");
                    fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"btsplot.pdf").c_str());
                } else if (plt==1) {
                    fprintf(f, "set terminal wxt font \"Arial, 8\"\n");
                    fprintf(f, "set output\n");
                }
                fprintf(f, "set xlabel \"time [seconds]\"\n");
                fprintf(f, "set ylabel \"photon count [photons/%lfsec]\"\n", corr_taumin*b);
                fprintf(f, "set title \"object description: %s\"\n", description.c_str());
                fprintf(f, "plot \"%s\" with steps\n", extract_file_name(tsfn).c_str());
                if (plt==1) fprintf(f, "pause -1\n");
                fprintf(f, "set xlabel \"time [seconds]\"\n");
                fprintf(f, "set ylabel \"photon count [Hz]\"\n");
                fprintf(f, "set title \"object description: %s\"\n", description.c_str());
                fprintf(f, "plot \"%s\" using 1:(($2)/%lf) with steps\n", extract_file_name(tsfn).c_str(), corr_taumin*b);
                if (plt==1) fprintf(f, "pause -1\n");
            }
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
    s+="pixel_size = "+floattostr(pixel_size)+" microns\n";
    double psf_r0=1.0/sqrt(1.0/detpsf_r0/detpsf_r0+1.0/expsf_r0/expsf_r0);
    s+="confocal_psf_system = sqrt(1/expsf_r0^2 + 1/detpsf_r0^2) = "+floattostr(psf_r0)+" microns\n";

    double Veff=pow(M_PI, 1.5)*(psf_r0*psf_r0*psf_r0*detpsf_z0/detpsf_r0);
    s+="confocal_focus_volume (Veff) = pi^(3/2) * psf_system^3 * detpsf_z0/detpsf_r0 = "+floattostr(Veff)+" femto litre\n";

    s+="illumination_distribution = "+ill_distribution_to_str(ill_distribution)+"\n";
    s+="detection_distribution = "+det_distribution_to_str(det_distribution)+"\n";
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
    s+="max_photons = "+inttostr(max_photons)+"\n";

    if (detector_type==1) {
        s+="detector_type = linear\n";
        s+="detector_resolution = "+inttostr(lindet_bits)+"bits  =>  range = 0..."+inttostr(pow(2,lindet_bits))+"\n";
        s+="detector_gain = "+floattostr(lindet_gain)+"\n";
        s+="detector_intensity_vs_variance = "+floattostr(lindet_var_factor)+"\n";
    } else {
        s+="detector_type = photon_counting\n";
    }

    s+="lambda_ex = "+floattostr(lambda_ex)+" nm\n";
    if (det_wavelength_min>0 && det_wavelength_max>0) {
        s+="det_wavelength_min = "+floattostr(det_wavelength_min)+" nm\n";
        s+="det_wavelength_max = "+floattostr(det_wavelength_max)+" nm\n";
    } else {
        s+="detecting all photons (no \"fluorescence filter\")\n";
    }
    if (background_rate>0) s+="background_rate = "+floattostr(background_rate)+" Hz (poissonian distribution)\n";
    else s+="no background photons\n";
    if (offset_rate>0 && offset_std>0) {
        s+="offset = "+floattostr(offset_rate/sim_timestep)+"Hz = "+floattostr(offset_rate)+"photons/detectionstep\n";
        s+="offset_std = "+floattostr(offset_std)+" photons   (gaussian distribution)\n";
    } else s+="no offset photons\n";
    if (offset_correction!=0) s+="offset_correction = "+floattostr(offset_correction)+" photons/detectionstep\n";
    else s+="no offset correction\n";
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
    if (correlator_type==3) s+=" (corr_directavg)";
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

 std::string FCSMeasurement::ill_distribution_to_str(int i) const {
    if (i==0) return "gaussian";
    if (i==1) return "gaussian_spim";
    return "unknown";
 }

 std::string FCSMeasurement::det_distribution_to_str(int i) const {
    if (i==0) return "gaussian";
    if (i==1) return "square_pixel";
    return "unknown";
 }


int FCSMeasurement::str_to_ill_distribution(std::string i) const {
    std::string s=tolower(i);
    if (s=="gaussian") return 0;
    if (s=="gaussian_spim") return 1;
    return strtol(s.c_str(), NULL, 10);
 }

int FCSMeasurement::str_to_det_distribution(std::string i) const {
    std::string s=tolower(i);
    if (s=="gaussian") return 0;
    if (s=="square_pixel") return 1;
    return strtol(s.c_str(), NULL, 10);
 }
