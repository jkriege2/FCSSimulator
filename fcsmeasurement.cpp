#include "fcsmeasurement.h"

#include <unistd.h>

FCSMeasurement::FCSMeasurement(FluorophorManager* fluorophors, std::string objectname):
    FluorescenceMeasurement(fluorophors, objectname)
{
    save_binning=false;
    online_correlation=false;
    save_timeseries=false;
    timeseries_ended=false;
    lastNPhotons=0;
    fccs_partner="";

    gaussbeam_pixel_normalization=1;

    psf_rz_image_rwidth=500;
    psf_rz_image_zwidth=1000;

    psf_rz_image_rresolution=0.01;
    psf_rz_image_zresolution=0.02;

    psf_rz_image_detection.resize(psf_rz_image_rwidth, psf_rz_image_zwidth);
    psf_rz_image_detection.setAll(0);
    psf_rz_image_illumination.resize(psf_rz_image_rwidth, psf_rz_image_zwidth);
    psf_rz_image_illumination.setAll(0);

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
    pixel_size_integrationdelta=0.04;

    detector_type=0;
    lindet_bits=14;
    lindet_gain=10;
    lindet_var_factor=10;
    lindet_readnoise=0;

    expsf_r0=0.5;
    expsf_z0=0.5*6;
    expsf_r02=0.6;
    expsf_z02=0.6*6;
    detpsf_r0=0.53;
    detpsf_z0=0.53*6;
    I0=200/0.25e-12;
    I02=200/0.25e-12;

    ex_x0=0;
    ex_y0=0;
    ex_z0=0;
    ex_x02=0;
    ex_y02=0;
    ex_z02=0;

    psfplot_xmax=2;
    psfplot_ymax=2;
    psfplot_zmax=10;

    max_photons=0xFFFFFFF;
    min_photons=0;

    // Rho 6G
    lambda_ex=488;
    lambda_ex2=0;

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

    save_arrivaltimes=false;
    arrivaltimes_onlyonce=true;

    plot_with.clear();

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

    correlator_fccs=NULL;
    corrjanb_fccs=NULL;
    corr_fccs=NULL;

    corr_tau=NULL;
    timeseries=NULL;
    timeseries_size=0;
    binned_timeseries=NULL;
    binned_timeseries_size=0;
    slots=0;

    correlation_runtime=0;

    ndettest_max=300;
    ndettest_step=1;

}

FCSMeasurement::~FCSMeasurement()
{
    if (correlator!=NULL) delete correlator; correlator=NULL;
    if (corrjanb!=NULL) delete corrjanb; corrjanb=NULL;
    if (correlator_fccs!=NULL) delete correlator_fccs; correlator_fccs=NULL;
    if (corrjanb_fccs!=NULL) delete corrjanb_fccs; corrjanb_fccs=NULL;
    if (timeseries!=NULL) free(timeseries); timeseries=NULL;; timeseries_size=0;
    if (binned_timeseries!=NULL) free(binned_timeseries); binned_timeseries=NULL; binned_timeseries_size=0;
}

void FCSMeasurement::clear() {
    if (correlator!=NULL) delete correlator; correlator=NULL;
    if (corrjanb!=NULL) delete corrjanb; corrjanb=NULL;
    if (correlator_fccs!=NULL) delete correlator_fccs; correlator_fccs=NULL;
    if (corrjanb_fccs!=NULL) delete corrjanb_fccs; corrjanb_fccs=NULL;
    if (timeseries!=NULL) free(timeseries); timeseries=NULL; timeseries_size=0;
    if (binned_timeseries!=NULL) free(binned_timeseries); binned_timeseries=NULL; binned_timeseries_size=0;
    slots=0;
    corr=NULL;
    corr_fccs=NULL;
    corr_tau=NULL;
    arrival_times.clear();
}

void FCSMeasurement::read_config_internal(jkINIParser2& parser) {
    FluorescenceMeasurement::read_config_internal(parser);
    online_correlation=parser.getSetAsBool("online_correlation", online_correlation);

    plot_with=tokenize_string(parser.getSetAsString("plot_with", fromStringVector(plot_with, ",")), ",");

    //plot_with_more.clear();
    for (int i=2; i<=9; i++) {
        if (parser.exists("plot_with"+inttostr(i))) {
            plot_with_more[i]=tokenize_string(parser.getSetAsString("plot_with"+inttostr(i), ""), ",");
        }
    }


    max_photons=parser.getAsInt("max_photons", max_photons);
    min_photons=parser.getAsInt("min_photons", min_photons);
    det_wavelength_min=parser.getAsDouble("det_wavelength_min", det_wavelength_min);
    det_wavelength_max=parser.getAsDouble("det_wavelength_max", det_wavelength_max);

    background_rate=parser.getAsDouble("background_rate", background_rate);
    offset_rate=parser.getAsDouble("offset_rate", offset_rate);
    offset_std=parser.getAsDouble("offset_std", offset_std);
    offset_correction=parser.getAsDouble("offset_correction", offset_correction);
    psf_region_factor=parser.getAsDouble("psf_region_factor", psf_region_factor);

    ill_distribution=str_to_ill_distribution(parser.getSetAsString("ill_distribution", ill_distribution_to_str(ill_distribution)));
    det_distribution=str_to_det_distribution(parser.getSetAsString("det_distribution", det_distribution_to_str(det_distribution)));

    pixel_size=parser.getSetAsDouble("pixel_size", pixel_size);
    pixel_size_integrationdelta=parser.getSetAsDouble("pixel_size_integrationdelta", pixel_size_integrationdelta);
    expsf_r0=parser.getSetAsDouble("expsf_r0", expsf_r0);
    expsf_z0=parser.getSetAsDouble("expsf_z0", expsf_z0);
    expsf_r02=parser.getSetAsDouble("expsf_r02", expsf_r02);
    expsf_z02=parser.getSetAsDouble("expsf_z02", expsf_z02);

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
    if (parser.exists("I02")) {
        I02=parser.getSetAsDouble("I02", I0);
    } else {
        I02=parser.getSetAsDouble("P02", I02*(M_PI*gsl_pow_2(2.0*1e-6*expsf_r02)))/(M_PI*gsl_pow_2(2.0*1e-6*expsf_r02));
    }
    lambda_ex=parser.getSetAsDouble("lambda_ex", lambda_ex);
    lambda_ex2=parser.getSetAsDouble("lambda_ex2", lambda_ex2);
    q_det=parser.getSetAsDouble("q_det", q_det);

    img_z0=parser.getSetAsDouble("img_z0", img_z0);
    img_x0=parser.getSetAsDouble("img_x0", img_x0);
    img_y0=parser.getSetAsDouble("img_y0", img_y0);

    ex_z0=parser.getSetAsDouble("ex_z0", img_z0);
    ex_x0=parser.getSetAsDouble("ex_x0", img_x0);
    ex_y0=parser.getSetAsDouble("ex_y0", img_y0);

    ex_z02=parser.getSetAsDouble("ex_z02", ex_z0);
    ex_x02=parser.getSetAsDouble("ex_x02", ex_x0);
    ex_y02=parser.getSetAsDouble("ex_y02", ex_y0);

    e_x=parser.getSetAsDouble("e_x", e_x);
    e_y=parser.getSetAsDouble("e_y", e_y);
    e_z=parser.getSetAsDouble("e_z", e_z);
    e_pol_fraction=parser.getSetAsDouble("e_pol_fraction", e_pol_fraction);

    double le=sqrt(e_x*e_x+e_y*e_y+e_z*e_z);
    e_x/=le; e_y/=le; e_z/=le;

    d_x=parser.getSetAsDouble("d_x", d_x);
    d_y=parser.getSetAsDouble("d_y", d_y);
    d_z=parser.getSetAsDouble("d_z", d_z);

    psfplot_xmax=parser.getSetAsDouble("psfplot_xmax", psfplot_xmax);
    psfplot_ymax=parser.getSetAsDouble("psfplot_ymax", psfplot_ymax);
    psfplot_zmax=parser.getSetAsDouble("psfplot_zmax", psfplot_zmax);

    le=sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
    d_x/=le; d_y/=le; d_z/=le;

    polarised_detection=parser.getSetAsBool("polarised_detection", polarised_detection);
    polarised_excitation=parser.getSetAsBool("polarised_excitation", polarised_excitation);
    save_arrivaltimes=parser.getSetAsBool("save_arrivaltimes", save_arrivaltimes);
    arrivaltimes_onlyonce=parser.getSetAsBool("arrivaltimes_onlyonce", arrivaltimes_onlyonce);

    save_binning=parser.getSetAsBool("save_binning", save_binning);
    save_binning_time=parser.getSetAsDouble("save_binning_time", save_binning_time);
    save_timeseries=parser.getSetAsBool("save_timeseries", save_timeseries);

    std::string dt=tolower(parser.getAsString("detector_type", "photon_counting"));
    if (dt=="0" || dt=="photon" || dt=="photon_counting") detector_type=0;
    if (dt=="1" || dt=="linear" || dt=="emccd" || dt=="ccd" || dt=="camera") detector_type=1;
    lindet_bits=parser.getSetAsInt("lindet_bits", lindet_bits);
    lindet_gain=parser.getSetAsDouble("lindet_gain", lindet_gain);
    lindet_var_factor=parser.getSetAsDouble("lindet_var_factor", lindet_var_factor);
    lindet_readnoise=parser.getSetAsDouble("lindet_readnoise", lindet_readnoise);
    fccs_partner=parser.getSetAsString("fccs_partner", fccs_partner);
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
      timeseries=(int32_t*)calloc((unsigned long long)ceil(duration/corr_taumin)+5, sizeof(int32_t));
      timeseries_size=(unsigned long long)ceil(duration/corr_taumin)+5;
    }
    if (save_binning && online_correlation) {
        int b=round(save_binning_time/corr_taumin);
        binned_timeseries=(int32_t*)calloc((timesteps)/b+10, sizeof(int32_t));
        binned_timeseries_size=(timesteps)/b+10;
    }
    slots=0;
    corr=NULL;
    corr_fccs=NULL;
    corr_tau=NULL;
    photoncounter=0;

    //unsigned int correlation_slots=P+(S-1)*(P-P/m);
    //corr=(double*)calloc(correlation_slots, sizeof(double));
    //corr_tau=(double*)calloc(correlation_slots, sizeof(double));
    correlator=new MultiTauCorrelator<double, double>(S, m, P, corr_taumin);
    corrjanb=new correlatorjb<double, double>(S, P, double(0.0));
    correlator_fccs=new MultiTauCorrelator<double, double>(S, m, P, corr_taumin);
    corrjanb_fccs=new correlatorjb<double, double>(S, P, double(0.0));
    Ephoton=6.626e-34*2.99e8/(lambda_ex*1e-9);
    if (lambda_ex2>0) Ephoton2=6.626e-34*2.99e8/(lambda_ex2*1e-9); else Ephoton2=0;

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

    psf_rz_image_detection.resize(psf_rz_image_rwidth, psf_rz_image_zwidth);
    psf_rz_image_detection.setAll(0);
    psf_rz_image_illumination.resize(psf_rz_image_rwidth, psf_rz_image_zwidth);
    psf_rz_image_illumination.setAll(0);
    gaussbeam_pixel_normalization=1;

    if (det_distribution<3) {
        for (uint32_t zz=0; zz<psf_rz_image_zwidth; zz++) {
            const double z=(double(zz)-double(psf_rz_image_zwidth/2))*psf_rz_image_zresolution;
            for (uint32_t rr=0; rr<psf_rz_image_rwidth; rr++) {
                const double r=double(rr)*psf_rz_image_rresolution;
                psf_rz_image_detection.setAt(rr,zz,detectionEfficiency(r,0,z));
                //if (zz==psf_rz_image_zwidth/2) std::cout<<"("<<zz<<", "<<rr<<") -> ("<<z<<", "<<r<<") = "<<detectionEfficiency(r,0,z)<<", "<<illuminationEfficiency(r,0,z)<<"\n";
            }
        }
    } else if (det_distribution==3) {
        for (uint32_t zz=0; zz<psf_rz_image_zwidth; zz++) {
            const double z=(double(zz)-double(psf_rz_image_zwidth/2))*psf_rz_image_zresolution;
            for (uint32_t rr=0; rr<psf_rz_image_rwidth; rr++) {
                const double r=double(rr)*psf_rz_image_rresolution;
                const double eff=gsl_pow_2(detpsf_r0/gaussbeam_w(z, detpsf_z0, detpsf_r0))*exp(-2.0*gsl_pow_2(r)/gsl_pow_2(gaussbeam_w(z, detpsf_z0, detpsf_r0)));
                psf_rz_image_detection.setAt(rr,zz,eff);
            }
        }

        gaussbeam_pixel_normalization=square_integrate(0, 0, 0, pixel_size, pixel_size_integrationdelta, psf_rz_image_detection, psf_rz_image_rresolution, psf_rz_image_zresolution);
    }

    for (uint32_t zz=0; zz<psf_rz_image_zwidth; zz++) {
        const double z=(double(zz)-double(psf_rz_image_zwidth)/2.0)*psf_rz_image_zresolution;
        for (uint32_t rr=0; rr<psf_rz_image_rwidth; rr++) {
            const double r=double(rr)*psf_rz_image_rresolution;
            psf_rz_image_illumination.setAt(rr, zz, illuminationEfficiency(r,0,z,expsf_r0,expsf_z0));
        }
    }


    estimate_psf_integrals();
    time(&start_time);
}

void FCSMeasurement::propagate(){
    tick();
    FluorescenceMeasurement::propagate();
    run_fcs_simulation();

    tock();
    runtime=runtime+get_duration();
}

void FCSMeasurement::finalize_sim(){
    tick();
    // calculate correlation function
    std::cout<<"correlating ... dur="<<duration<<"\n";
    PublicTickTock tim;
    tim.tick();
    long* taus=(long*)malloc(S*P*sizeof(long));
    FCSMeasurement* partner=get_fccs_partner_object();
    if (!online_correlation) {
        if (correlator_type==0) {
            for (unsigned long long i=0; i<timesteps; i++) {
                correlator->correlate_step(timeseries[i]);
            }
        } else if (correlator_type==1) {
            for (unsigned long long i=0; i<timesteps; i++) {
                corrjanb->run(timeseries[i], timeseries[i]);
            }
        } else if (correlator_type==2) {
            statisticsAutocorrelateCreateMultiTau(taus, S, m, P);
            corr=(double*)malloc(S*P*sizeof(double));
            statisticsAutocorrelateMultiTauSymmetric<int32_t,int64_t>(corr, timeseries, timesteps, taus, S*P);
        } else if (correlator_type==3) {
            statisticsAutocorrelateCreateMultiTau(taus, S, m, P);
            corr=(double*)malloc(S*P*sizeof(double));
            statisticsAutocorrelateMultiTauAvgSymmetric<int32_t,int64_t,int64_t>(corr, timeseries, timesteps, S, m, P, 1);
        }

        if (partner) {
            if (timesteps==partner->timesteps && corr_taumin==partner->corr_taumin) {
                if (correlator_type==0) {
                    for (unsigned long long i=0; i<timesteps; i++) {
                        correlator_fccs->crosscorrelate_step(timeseries[i], partner->timeseries[i]);
                    }
                } else if (correlator_type==1) {
                    for (unsigned long long i=0; i<timesteps; i++) {
                        corrjanb_fccs->run(timeseries[i], partner->timeseries[i]);
                    }
                } else if (correlator_type==2) {
                    corr_fccs=(double*)malloc(S*P*sizeof(double));
                    statisticsCrosscorrelateMultiTauSymmetric<int32_t,int64_t>(corr_fccs, timeseries, partner->timeseries, timesteps, taus, S*P);
                } else if (correlator_type==3) {
                    corr_fccs=(double*)malloc(S*P*sizeof(double));
                    statisticsCrosscorrelateMultiTauAvgSymmetric<int32_t,int64_t,int64_t>(corr_fccs, timeseries, partner->timeseries, timesteps, S, m, P, 1);
                }
            } else {
                throw MeasurementException("The FCCS partner object recorded with different settings! Crosscorreation not possible!");
            }
        }
    }
    if (correlator_type==0) {
        correlator->normalize();
        corr=correlator->getCor();
        corr_tau=correlator->getCorTau();
        if (partner) {
            correlator_fccs->crossnormalize();
            corr_fccs=correlator_fccs->getCor();
        }
        slots=correlator->getSlots();
    } else if (correlator_type==1) {
        slots=S*P;
        double** corr1=corrjanb->get_array_G();
        corr=(double*)malloc(slots*sizeof(double));
        corr_tau=(double*)malloc(slots*sizeof(double));
        for (unsigned int i=0; i<slots; i++) {
            corr_tau[i]=corr1[0][i]*corr_taumin;
            corr[i]=corr1[1][i];
        }
        if (partner) {
            corr1=corrjanb_fccs->get_array_G();
            corr_fccs=(double*)malloc(slots*sizeof(double));
            for (unsigned int i=0; i<slots; i++) {
                corr_fccs[i]=corr1[1][i];
            }
        }
    } else if (correlator_type==2 || correlator_type==3) {
        slots=S*P;
        corr_tau=(double*)malloc(slots*sizeof(double));
        for (unsigned int i=0; i<slots; i++) {
            corr_tau[i]=(double)taus[i]*corr_taumin;
        }
    }

    free(taus);

    tim.tock();
    correlation_runtime+=tim.get_duration();
    std::cout<<" done!\n";
    tock();
    runtime=runtime+get_duration();
}

void FCSMeasurement::estimate_psf_integrals() {
    double maxs=0;
    psf_rz_image_detection_integral_max=1;
    for (uint32_t z=0; z<psf_rz_image_detection.height(); z++) {
        double sum=0;
        double zz=(double(z)-double(psf_rz_image_detection.height()/2))*psf_rz_image_zresolution;
        for (uint32_t r=0; r<psf_rz_image_detection.width(); r++) {
            double rr=M_PI*gsl_pow_2(double(r+1)*psf_rz_image_rresolution)-M_PI*gsl_pow_2(double(r)*psf_rz_image_rresolution);
            sum=sum+rr*psf_rz_image_detection.get(r, z);
            //std::cout<<r<<", "<<psf_rz_image_detection(r, z)<<"|  ";
        }
        //sum=square_integrate(0,0,zz,2*psf_rz_image_detection.width()*psf_rz_image_rresolution, psf_rz_image_rresolution, psf_rz_image_detection, psf_rz_image_rresolution, psf_rz_image_zresolution);
        if (sum>maxs) maxs=sum;
        //std::cout<<zz<<"("<<z<<"): "<<sum<<"\n";
    }
    psf_rz_image_detection_integral_max=maxs;
}


double FCSMeasurement::square_integrate(double dx, double dy, double dz, double a, double da, const JKImage<double>& image, double image_dr, double image_dz) const {
    int apoints=ceil(a/da);
    if (apoints<10) apoints=10;
    double daa=a/double(apoints);
    double sum=0;
    double zz=double(image.height()/2)+dz/image_dz;
    long long zzi=round(zz);
    if (zzi<0) zzi=0;
    if (zz>=0 && zz<=image.height()) {
        for (int i=0; i<apoints*apoints; i++) {
            const double x=double(i%apoints)*daa+daa/2.0-a/2.0;
            const double y=double(i/apoints)*daa+daa/2.0-a/2.0;
            const double r=sqrt(gsl_pow_2(x-dx)+gsl_pow_2(y-dy))/image_dr;
            const double res=image.bilinear_interpolate<double>(r, zz);
            sum=sum+res;
            //if (i==5) std::cout<<r<<", "<<zz<<":  "<<res<<"\n";
        }
        sum=sum*daa*daa/psf_rz_image_detection_integral_max;
    } else {
        sum=0;
    }
    return sum;
}

double FCSMeasurement::detectionEfficiency(double dx, double dy, double dz) const {
    double eff=0;
    const double sq2=sqrt(2.0);
    const double r=sqrt(dx*dx+dy*dy);
    switch (det_distribution) {
        case 0:
            eff=exp(-2.0*(gsl_pow_2(dx)+gsl_pow_2(dy))/gsl_pow_2(detpsf_r0)-2.0*gsl_pow_2(dz)/gsl_pow_2(detpsf_z0));
            break;
        case 1:
            eff=exp(-2.0*gsl_pow_2(dz)/gsl_pow_2(detpsf_z0))
                *(erf((pixel_size-2.0*dx)/(sq2*detpsf_r0))+erf((pixel_size+2.0*dx)/(detpsf_r0*sq2)))/(2.0*erf(pixel_size/(sq2*detpsf_r0)))
                *(erf((pixel_size-2.0*dy)/(sq2*detpsf_r0))+erf((pixel_size+2.0*dy)/(detpsf_r0*sq2)))/(2.0*erf(pixel_size/(sq2*detpsf_r0)));
            break;
        case 2:
            eff=gsl_pow_2(detpsf_r0/gaussbeam_w(dz, detpsf_z0, detpsf_r0))*exp(-2.0*gsl_pow_2(r)/gsl_pow_2(gaussbeam_w(dz, detpsf_z0, detpsf_r0)));
            break;
        case 3:
            eff=square_integrate(dx, dy, dz, pixel_size, pixel_size_integrationdelta, psf_rz_image_detection, psf_rz_image_rresolution, psf_rz_image_zresolution)/gaussbeam_pixel_normalization;
            break;
    }
    return eff;
}

double FCSMeasurement::illuminationEfficiency(double dx, double dy, double dz, double expsf_r0, double expsf_z0) const {
    double eff=0;
    double zz=0;
    switch (ill_distribution) {
        case 0:
            eff=exp(-2.0*(gsl_pow_2(dx)+gsl_pow_2(dy))/gsl_pow_2(expsf_r0)-2.0*gsl_pow_2(dz)/gsl_pow_2(expsf_z0));
            break;
        case 1:
            eff=exp(-2.0*gsl_pow_2(dz)/gsl_pow_2(expsf_z0));
            break;
        case 2:
            zz=dz/expsf_z0;
            eff=gsl_pow_2(gsl_sf_sinc(zz));
            break;
    }
    return eff;
}

void FCSMeasurement::run_fcs_simulation(){
    // number of fluorescence photons per molecule and sim_timestep, but you have to multiply by sigma_abs and q_fluor (walker dependent!)
    double n00=I0*1e-6/Ephoton*sim_timestep;
    double n02=I02*1e-6/Ephoton2*sim_timestep;
    int bin_r=round(save_binning_time/corr_taumin); // binning ration
    std::string walker_cnt="";
    // first we go through alle the fluorophors and integrate their contribution
    for (size_t d=0; d<dyn.size(); d++) { // go through all dynamics objects that provide data for this measurement object
        FluorophorDynamics::walkerState* dynn=dyn[d]->get_visible_walker_state();
        unsigned long wc=dyn[d]->get_visible_walker_count();
        if (walker_cnt.size()>0) walker_cnt=walker_cnt+", ";
        walker_cnt=walker_cnt+dyn[d]->get_object_name()+":"+inttostr(wc);
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
                double n0sum=0;
                double dxs=dx*dx;
                double dys=dy*dy;



                // first excitation focus
                double edx=x0-ex_x0;
                double edz=z0-ex_z0;
                double edy=y0-ex_y0;
                if ((fabs(dz)<(psf_region_factor*detpsf_z0)) && ((dxs+dys)<gsl_pow_2(psf_region_factor*detpsf_r0))) {
                    double n0=dyn[d]->get_visible_walker_sigma_times_qfl(i);
                    n0=n0*n00*fluorophors->get_spectral_absorbance(dynn[i].spectrum, lambda_ex);
                    n0=n0*illuminationEfficiency(edx, edy, edz, expsf_r0, expsf_z0);

                    n0sum=n0sum+n0;
                    //std::cout<<"nphot_sum="<<nphot_sum<<"\n";


                    // second excitation focus
                    if (lambda_ex2>0) {
                        edx=x0-ex_x02;
                        edz=z0-ex_z02;
                        edy=y0-ex_y02;
                        double n0=dyn[d]->get_visible_walker_sigma_times_qfl(i);
                        n0=n0*n02*fluorophors->get_spectral_absorbance(dynn[i].spectrum, lambda_ex2);
                        n0=n0*illuminationEfficiency(edx, edy, edz, expsf_r02, expsf_z02);

                        n0sum=n0sum+n0;
                        //std::cout<<"nphot_sum="<<nphot_sum<<"\n";
                    }


                }

                // excitation polarisation
                double dpx, dpy, dpz;
                dpx=dynn[i].p_x;
                dpy=dynn[i].p_y;
                dpz=dynn[i].p_z;
                if (polarised_excitation) {
                    n0sum=n0sum*(1.0-e_pol_fraction+e_pol_fraction*gsl_pow_2(dpx*e_x+dpy*e_y+dpz*e_z));
                }

                // and the detection!
                if (polarised_detection) {
                    n0sum=n0sum*gsl_pow_2(dpx*d_x+dpy*d_y+dpz*d_z);
                }
                if (det_wavelength_min>0 && det_wavelength_max>0) {
                    n0sum=n0sum*fluorophors->get_spectral_fluorescence(dynn[i].spectrum, det_wavelength_min, det_wavelength_max);
                }
                nphot_sum=nphot_sum+n0sum*q_det*detectionEfficiency(dx, dy, dz);

            }
        }
    }



    // if the current integration step ended, we may add a new value to the correlator
    //std::cout<<"sim_time="<<sim_time<<"   endCurrentStep="<<endCurrentStep<<"\n";
    if (sim_time>=endCurrentStep) {
        endCurrentStep=sim_time+corr_taumin;
        register int32_t N=getDetectedPhotons(nphot_sum);

        photoncounter=photoncounter+N;

        if (save_arrivaltimes && N>0) {
            if (arrivaltimes_onlyonce) arrival_times.push_back(sim_time);
            else {
                for (int i=0; i<N; i++) arrival_times.push_back(sim_time);
            }
        }
        //std::cout<<"nphot_sum="<<nphot_sum<<"    N="<<N<<std::endl;
        if (online_correlation) {
            if (correlator_type==0) {
                correlator->correlate_step(N);
            } else if (correlator_type==1) {
                corrjanb->run(N,N);
            }
            lastNPhotons=N;
            timeseries_ended=false;

            FCSMeasurement* partner=get_fccs_partner_object();
            if (partner) {
                if (timesteps==partner->timesteps && corr_taumin==partner->corr_taumin) {
                    if (correlator_type==0) {
                        correlator_fccs->crosscorrelate_step(N, partner->lastNPhotons);
                    } else if (correlator_type==1) {
                        corrjanb_fccs->run(N,partner->lastNPhotons);
                    }
                } else {
                    throw MeasurementException("The FCCS partner object recorded with different settings! Crosscorreation not possible!");
                }
            }
        } else {
            if (current_timestep<timeseries_size) {
                timeseries[current_timestep]=N;
                timeseries_ended=false;
            } else {
                timeseries_ended=true;
            }
        }
        if (save_binning && online_correlation) {
            bin_sum=bin_sum+N;
            bin_counter++;
            if (int64_t(bin_counter)==bin_r) {
                bin_counter=0;
                binned_timeseries[bin_i]=bin_sum;
                bin_sum=0;
                bin_i++;
            }
        }
        if (current_timestep%mmax(1,(timesteps/1000))==0) {
            display_temp+=N;
            time_t now;
            time(&now);
            double seconds=difftime(now, start_time);
            double percent=((double)current_timestep)/(timesteps)*100.0;
            double eta=seconds/percent*100.0-seconds;
            int hrs=floor(eta/3600.0);
            int mins=floor((eta-double(hrs*3600))/60.0);
            double secs=floor(eta-double(hrs*3600)-double(mins*60));
            //std::cout<<format("%4.1lf", percent)<<"% (ETA: "+format("%8.2lf", eta/60.0)+"min):   "<<repeated_string("     ", object_number-1)<<display_temp<<std::endl;
            std::cout<<format("%4.1lf", percent)<<"% (ETA: "+inttostr(hrs)+":"+inttostr(mins)+":"+inttostr(secs)+"):   "<<repeated_string("     ", object_number-1)<<display_temp<<"  wc=["<<walker_cnt<<"]"<<std::endl;
            display_temp=0;
        } else {
            display_temp+=N;
        }
        current_timestep++;
        nphot_sum=0;
    }
}

int32_t FCSMeasurement::getDetectedPhotons(double nphot_sum) const {
    register int32_t N=0;
    if (detector_type==0) {
        N=gsl_ran_poisson(rng, nphot_sum);
    } else {
        double d=gsl_ran_gaussian(rng, sqrt(nphot_sum*lindet_gain*lindet_gain*lindet_var_factor+lindet_readnoise*lindet_readnoise))+nphot_sum*lindet_gain;
        N=round(d);
        if (d>max_photons) N=max_photons;
    }
    //if (detector_type==1) std::cout<<"#      nphot_sum="<<nphot_sum<<"  lindet_bits="<<lindet_bits"\n";
    //if (detector_type==1) std::cout<<"*      "<<N<<"\n";
    if (background_rate>0) {
        N=N+gsl_ran_poisson(rng, background_rate*corr_taumin);
    }
    //if (detector_type==1) std::cout<<"**     "<<N<<"\n";
    if (offset_rate>0 && offset_std>0) {
        double o=gsl_ran_gaussian(rng, offset_std)+offset_rate;
        N=N+round(o);
    } else if (offset_rate>0 && offset_std<=0) {
        N=N+round(offset_rate);
    }
    //if (detector_type==1) std::cout<<"***    "<<N<<"\n";
    N=N-offset_correction;
    //if (detector_type==1) std::cout<<"****   "<<N<<"\n";

    if (detector_type==1) {
        const int32_t Nmax=gsl_pow_int(2,lindet_bits)-1;
        if (N>Nmax) N=Nmax;
    }
    //if (detector_type==1) std::cout<<"*****  "<<N<<"\n";
    if (N<min_photons) N=min_photons;
    //if (detector_type==1) std::cout<<"****** "<<N<<"\n";
    if (N>max_photons) N=max_photons;
    //if (detector_type==1) std::cout<<"*******"<<N<<"\n";
    return N;
}

void FCSMeasurement::save() {
    FILE* f;
    char fn[255], fno[255];


    sprintf(fn, "%s%scorr.dat", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    f=fopen(fn, "w");
    unsigned long long istart=1;
    if (correlator_type==1) istart=0;
    if (correlator_type==2) istart=0;
    if (correlator_type==3) istart=0;

    for (unsigned long long i=istart; i<slots; i++) {
        if (correlator_type==1) fprintf(f, "%15.10lf, %15.10lf\n", corr_tau[i], corr[i]);
        else                    fprintf(f, "%15.10lf, %15.10lf\n", corr_tau[i], corr[i]);
    }
    fclose(f);
    std::cout<<" done!\n";
    std::string corrfn=fn;

    sprintf(fn, "%s%scorrplot.plt", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
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

    fprintf(f, "Veffa=sqrt(pi*pi*pi)*wxy*wxy*wxy*1e-15*gammaa\n");
    fprintf(f, "Vefff=sqrt(pi*pi*pi)*wxy*wxy*wxy*1e-15*gammaf\n");
    fprintf(f, "Veff=sqrt(pi*pi*pi)*wxy*wxy*wxy*1e-15*gamma\n");
    for (int plt=0; plt<2; plt++) {
        if (plt==0) {
            fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
            fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"corrplot.pdf").c_str());
        } else if (plt==1) {
            fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
            fprintf(f, "set output\n");
        }
        fprintf(f, "set logscale x\n");
        fprintf(f, "set title \"object description: %s\"\n", description.c_str());
        fprintf(f, "plot \"%s\" title \"simulation data\" with points, "
                   "g(x,N,tauD,gamma,1) title sprintf(\"fit N=%%.3f, tauD=%%.3f uS, gamma=%%.3f, D=%%.3f um^2/s, c=%%.3f nM\",N, tauD*1e6, gamma, wxy*wxy/4.0/tauD, N/6.022e14/Veff), "
                   "g(x,Nf,tauDf,gammaf,1) title sprintf(\"fit N=%%.3f, tauD=%%.3f uS, gamma=%%.3f, D=%%.3f um^2/s, c=%%.3f nM\",Nf, tauDf*1e6, gammaf, wxy*wxy/4.0/tauDf, Nf/6.022e14/Vefff)"
                   "\n", extract_file_name(corrfn).c_str());
        if (plt==1) fprintf(f, "pause -1\n");
        fprintf(f, "Na=Nf\n");
        fprintf(f, "tauDa=tauDf\n");
        fprintf(f, "fit g(x, Na, tauDa, gammaa, alphaa) \"%s\" via Na, tauDa, alphaa\n", extract_file_name(corrfn).c_str());
        fprintf(f, "set logscale x\n");
        fprintf(f, "set title \"object description: %s\"\n", description.c_str());
        fprintf(f, "plot \"%s\" title \"simulation data\" with points, "
                   "g(x,N,tauD,gamma,1) title sprintf(\"fit N=%%.3f, tauD=%%.3f uS, gamma=%%.3f, D=%%.3f um^2/s, c=%%.3f nM\",N, tauD*1e6, gamma, wxy*wxy/4.0/tauD, N/6.022e14/Veff), "
                   "g(x,Nf,tauDf,gammaf,1) title sprintf(\"fit N=%%.3f, tauD=%%.3f uS, gamma=%%.3f, D=%%.3f um^2/s, c=%%.3f nM\",Nf, tauDf*1e6, gammaf, wxy*wxy/4.0/tauDf, Nf/6.022e14/Vefff), "
                   "g(x,Na,tauDa,gammaa,alphaa) title sprintf(\"fit N=%%.3f, tauD=%%.3f uS, gamma=%%.3f, alpha=%%.3f, D=%%.3f um^2/s, c=%%.3f nM\",Na, tauDa*1e6, gammaa,alphaa,  wxy*wxy/4.0/tauDa, Na/6.022e14/Veffa)"
                   "\n", extract_file_name(corrfn).c_str());
        if (plt==1) fprintf(f, "pause -1\n");
    }

    fclose(f);
    std::cout<<" done!\n";

    sprintf(fn, "%s%scorrplot_simple.plt", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    f=fopen(fn, "w");
    for (int plt=0; plt<2; plt++) {
        if (plt==0) {
            fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
            fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"corrplot_simple.pdf").c_str());
        } else if (plt==1) {
            fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
            fprintf(f, "set output\n");
        }

        fprintf(f, "set logscale x\n");
        fprintf(f, "set title \"object description: %s\"\n", description.c_str());
        fprintf(f, "plot \"%s\" title \"simulation data\" with points\n", extract_file_name(corrfn).c_str());
    }
    fprintf(f, "pause -1\n");
    fclose(f);
    std::cout<<" done!\n";



    FCSMeasurement* partner=get_fccs_partner_object();
    if (partner) {
        if (timesteps==partner->timesteps && corr_taumin==partner->corr_taumin) {
            char fn2[255], fn1[255];
            sprintf(fn, "%s%scrosscorr.dat", basename.c_str(), object_name.c_str());
            sprintf(fn1, "%s%scorr.dat", basename.c_str(), object_name.c_str());
            sprintf(fn2, "%s%scorr.dat", basename.c_str(), partner->get_object_name().c_str());
            std::cout<<"writing '"<<fn<<"' ...";
            f=fopen(fn, "w");
            unsigned long long istart=1;
            if (correlator_type==1) istart=0;
            if (correlator_type==2) istart=0;
            if (correlator_type==3) istart=0;

            for (unsigned long long i=istart; i<slots; i++) {
                if (correlator_type==1) fprintf(f, "%15.10lf, %15.10lf\n", corr_tau[i], corr_fccs[i]);
                else fprintf(f, "%15.10lf, %15.10lf\n", corr_tau[i], corr_fccs[i]);
            }
            fclose(f);
            std::cout<<" done!\n";
            std::string ccffn=fn;
            std::string acf1fn=fn1;
            std::string acf2fn=fn2;

            sprintf(fn, "%s%scrosscorrplot_simple.plt", basename.c_str(), object_name.c_str());
            std::cout<<"writing '"<<fn<<"' ...";
            f=fopen(fn, "w");
            for (int plt=0; plt<2; plt++) {
                if (plt==0) {
                    fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                    fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"crosscorrplot_simple.pdf").c_str());
                } else if (plt==1) {
                    fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
                    fprintf(f, "set output\n");
                }

                fprintf(f, "set logscale x\n");
                fprintf(f, "set title \"object description: %s\"\n", description.c_str());
                fprintf(f, "plot \"%s\" title \"ACF1\" with lines lc rgb \"dark-green\",  \"%s\" title \"ACF2\" with lines lc rgb \"dark-red\",  \"%s\" title \"CCF\" with lines lc rgb \"dark-blue\"\n", extract_file_name(acf1fn).c_str(), extract_file_name(acf2fn).c_str(), extract_file_name(ccffn).c_str());
            }
            fprintf(f, "pause -1\n");
            fclose(f);
            std::cout<<" done!\n";


        } else {
            throw MeasurementException("The FCCS partner object recorded with different settings! Crosscorreation not possible!");
        }
    }



    sprintf(fn, "%s%sdettest.dat", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    f=fopen(fn, "w");
    long ndettests=1000;
    for (  double n_phot=0; n_phot<=ndettest_max; n_phot=n_phot+ndettest_step) {
        int32_t* Nm=(int32_t*)calloc(ndettests, sizeof(int32_t));
        for (int i=0; i<ndettests; i++) {
            Nm[i]=getDetectedPhotons(n_phot);
        }

        fprintf(f, "%15.10lf, %15.10lf, %15.10lf\n", n_phot, statisticsAverage(Nm, ndettests), statisticsVariance(Nm, ndettests));
        free(Nm);
    }
    fclose(f);
    std::cout<<" done!\n";
    std::string dettestfn=fn;



    sprintf(fn, "%s%sdettest.plt", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    f=fopen(fn, "w");
    fprintf(f, "reset\n");
    fprintf(f, "f(x,offset,a)=offset+a*x\n");
    fprintf(f, "offset=0\n");
    fprintf(f, "a=1\n");
    fprintf(f, "stats \"%s\" using 2:3 name 'DETTEST'\n", extract_file_name(dettestfn).c_str());
    fprintf(f, "offset=DETTEST_intercept\n");
    fprintf(f, "a=DETTEST_slope\n");
    fprintf(f, "fit f(x,offset,a) \"%s\" using 2:3 via offset, a\n", extract_file_name(dettestfn).c_str());
    for (int plt=0; plt<2; plt++) {
        if (plt==0) {
            fprintf(f, "set terminal pdfcairo color solid font \"%s, 8\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
            fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"dettest.pdf").c_str());
        } else if (plt==1) {
            fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
            fprintf(f, "set output\n");
        }

        fprintf(f, "set size noratio\n");

        fprintf(f, "set xlabel \"average signal I\"\n");
        fprintf(f, "set ylabel \"signal variance sigma_I^2\"\n");
        fprintf(f, "set title \"detector signal vs. variance: %s\"\n", description.c_str());
        fprintf(f, "plot \"%s\" using 2:3 title 'data' with points"
        ", f(x, offset, a) title sprintf('fit %%f + %%f * I', offset, a) with lines"
        "\n", extract_file_name(dettestfn).c_str());
        if (plt==1) fprintf(f, "pause -1\n");

    }

    fclose(f);
    std::cout<<" done!\n";



    sprintf(fn, "%s%sdetbackgroundtest.dat", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    f=fopen(fn, "w");
    int32_t* Nm=(int32_t*)calloc(ndettests, sizeof(int32_t));
    int32_t* Nm05=(int32_t*)calloc(ndettests, sizeof(int32_t));
    int32_t* Nm1=(int32_t*)calloc(ndettests, sizeof(int32_t));
    int32_t* Nm2=(int32_t*)calloc(ndettests, sizeof(int32_t));
    int32_t* Nm3=(int32_t*)calloc(ndettests, sizeof(int32_t));
    double IperParticle=q_det*dyn[0]->get_init_q_fluor(0)*I0*1e-6/Ephoton*dyn[0]->get_init_sigma_abs(0)*corr_taumin;

    int backtest_sum=20;
    for (int i=0; i<ndettests; i++) {
        Nm[i]=getDetectedPhotons(0);
        Nm05[i]=getDetectedPhotons(0.5*IperParticle);
        Nm1[i]=getDetectedPhotons(1.0*IperParticle);
        Nm2[i]=getDetectedPhotons(5.0*IperParticle);
        Nm3[i]=getDetectedPhotons(10.0*IperParticle);
        if (i%backtest_sum==0 && i>0) {
            fprintf(f, "%15.10lf, %15.10lf, %15.10lf, %15.10lf, %15.10lf, %15.10lf\n",double(i)*corr_taumin , statisticsAverage(&(Nm[i-backtest_sum]), backtest_sum) , statisticsAverage(&(Nm05[i-backtest_sum]), backtest_sum), statisticsAverage(&(Nm1[i-backtest_sum]), backtest_sum), statisticsAverage(&(Nm2[i-backtest_sum]), backtest_sum), statisticsAverage(&(Nm3[i-backtest_sum]), backtest_sum));
        }
    }
    fprintf(f, "\n\n");
    for (  double N=0; N<30; N++) {
        int32_t* Nm=(int32_t*)calloc(ndettests, sizeof(int32_t));
        for (int i=0; i<ndettests; i++) {
            Nm[i]=getDetectedPhotons(N*IperParticle);
        }

        fprintf(f, "%lf, %15.10lf, %15.10lf,  %15.10lf\n", N, statisticsAverage(Nm, ndettests), statisticsVariance(Nm, ndettests), N*IperParticle);
        free(Nm);
    }

    double backMean=statisticsAverage(Nm, ndettests);
    double backVar=statisticsVariance(Nm, ndettests);
    double backMean05=statisticsAverage(Nm05, ndettests);
    double backVar05=statisticsVariance(Nm05, ndettests);
    double backMean1=statisticsAverage(Nm1, ndettests);
    double backVar1=statisticsVariance(Nm1, ndettests);
    double backMean2=statisticsAverage(Nm2, ndettests);
    double backVar2=statisticsVariance(Nm2, ndettests);
    double backMean3=statisticsAverage(Nm3, ndettests);
    double backVar3=statisticsVariance(Nm3, ndettests);
    free(Nm);
    free(Nm1);
    free(Nm2);
    free(Nm3);

    fclose(f);
    std::cout<<" done!\n";
    std::string detbacktestfn=fn;



    sprintf(fn, "%s%sdetbackgroundtest.plt", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    f=fopen(fn, "w");
    fprintf(f, "reset\n");
    fprintf(f, "bmean=%15.10lf\n", backMean);
    fprintf(f, "bvar=%15.10lf\n", backVar);
    fprintf(f, "bmean05=%15.10lf\n", backMean05);
    fprintf(f, "bvar05=%15.10lf\n", backVar05);
    fprintf(f, "bmean1=%15.10lf\n", backMean1);
    fprintf(f, "bvar1=%15.10lf\n", backVar1);
    fprintf(f, "bmean2=%15.10lf\n", backMean2);
    fprintf(f, "bvar2=%15.10lf\n", backVar2);
    fprintf(f, "bmean3=%15.10lf\n", backMean3);
    fprintf(f, "bvar3=%15.10lf\n", backVar3);
    fprintf(f, "dt=%15.10lf\n", corr_taumin);
    for (int plt=0; plt<2; plt++) {
        if (plt==0) {
            fprintf(f, "set terminal pdfcairo color solid font \"%s, 8\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
            fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"detbackgroundtest.pdf").c_str());
        } else if (plt==1) {
            fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
            fprintf(f, "set output\n");
        }

        fprintf(f, "set size noratio\n");

        fprintf(f, "set xlabel \"time t [seconds]\"\n");
        fprintf(f, "set ylabel \"background signal\"\n");
        fprintf(f, "set title \"background/1-photon signal: %s\"\n", description.c_str());
        fprintf(f, "plot \"%s\" using 1:2 index 0 title sprintf('background, (%%f +/- %%f)', bmean, sqrt(bvar)) with lines lc rgb \"red\""
        ", \"%s\" using 1:4 index 0 title sprintf('N=1, (%%f +/- %%f)', bmean1, sqrt(bvar1)) with lines lc rgb \"blue\""
        ", bmean notitle with lines lc rgb \"dark-red\""
        //", bmean-sqrt(bvar) notitle with lines lc rgb \"dark-red\""
        //", bmean+sqrt(bvar) notitle with lines lc rgb \"dark-red\""
        ", bmean1 notitle with lines lc rgb \"dark-blue\""
        //", bmean1-sqrt(bvar1) notitle with lines lc rgb \"dark-blue\""
        //", bmean1+sqrt(bvar1) notitle with lines lc rgb \"dark-blue\""
        "\n", extract_file_name(detbacktestfn).c_str(), extract_file_name(detbacktestfn).c_str());
        if (plt==1) fprintf(f, "pause -1\n");

        fprintf(f, "set title \"N-photon signal: %s\"\n", description.c_str());
        fprintf(f, "plot \"%s\" using 1:2 index 0 title sprintf('background, (%%f +/- %%f)', bmean, sqrt(bvar)) with lines lc rgb \"red\""
        ", \"%s\" using 1:3 index 0 title sprintf('N=0.5, (%%f +/- %%f)', bmean05, sqrt(bvar05)) with lines lc rgb \"green\""
        ", \"%s\" using 1:4 index 0 title sprintf('N=1, (%%f +/- %%f)', bmean1, sqrt(bvar1)) with lines lc rgb \"blue\""
        ", \"%s\" using 1:5 index 0 title sprintf('N=5, (%%f +/- %%f)', bmean2, sqrt(bvar2)) with lines lc rgb \"magenta\""
        ", \"%s\" using 1:6 index 0 title sprintf('N=10, (%%f +/- %%f)', bmean3, sqrt(bvar3)) with lines lc rgb \"yellow\""
        ", bmean notitle with lines lc rgb \"dark-red\""
        //", bmean-sqrt(bvar) notitle with lines lc rgb \"dark-red\""
        //", bmean+sqrt(bvar) notitle with lines lc rgb \"dark-red\""
        ", bmean05 notitle with lines lc rgb \"dark-green\""
        //", bmean05-sqrt(bvar05) notitle with lines lc rgb \"dark-green\""
        //", bmean05+sqrt(bvar05) notitle with lines lc rgb \"dark-green\""
        ", bmean1 notitle with lines lc rgb \"dark-blue\""
        //", bmean1-sqrt(bvar1) notitle with lines lc rgb \"dark-red\""
        //", bmean1+sqrt(bvar1) notitle with lines lc rgb \"dark-red\""
        ", bmean2 notitle with lines lc rgb \"dark-magenta\""
        //", bmean2-sqrt(bvar2) notitle with lines lc rgb \"dark-magenta\""
        //", bmean2+sqrt(bvar2) notitle with lines lc rgb \"dark-magenta\""
        ", bmean3 notitle with lines lc rgb \"dark-yellow\""
        //", bmean3-sqrt(bvar3) notitle with lines lc rgb \"dark-yellow\""
        //", bmean3+sqrt(bvar3) notitle with lines lc rgb \"dark-yellow\""
        "\n", extract_file_name(detbacktestfn).c_str(), extract_file_name(detbacktestfn).c_str(), extract_file_name(detbacktestfn).c_str(), extract_file_name(detbacktestfn).c_str(), extract_file_name(detbacktestfn).c_str());
        if (plt==1) fprintf(f, "pause -1\n");


        fprintf(f, "set xlabel \"average particle number\"\n");
        fprintf(f, "set ylabel \"signal\"\n");
        fprintf(f, "set title \"particle number vs. signal: %s\"\n", description.c_str());
        fprintf(f, "plot [0:10] \"%s\" using 1:2:3 index 1 notitle with lines lc rgb \"red\""
        ", bmean notitle with lines lc rgb \"dark-red\""
        "\n", extract_file_name(detbacktestfn).c_str());
        if (plt==1) fprintf(f, "pause -1\n");

        fprintf(f, "set xlabel \"average particle number\"\n");
        fprintf(f, "set ylabel \"signal\"\n");
        fprintf(f, "set title \"particle number vs. signal: %s\"\n", description.c_str());
        fprintf(f, "plot \"%s\" using 1:2:3 index 1 notitle with lines lc rgb \"red\""
        ", bmean notitle with lines lc rgb \"dark-red\""
        "\n", extract_file_name(detbacktestfn).c_str());
        if (plt==1) fprintf(f, "pause -1\n");

        fprintf(f, "set xlabel \"average particle number\"\n");
        fprintf(f, "set ylabel \"signal\"\n");
        fprintf(f, "set title \"particle number vs. signal: %s\"\n", description.c_str());
        fprintf(f, "plot \"%s\" using 1:2:3 index 1 notitle with yerrorbars lc rgb \"red\""
        ", bmean notitle with lines lc rgb \"dark-red\""
        "\n", extract_file_name(detbacktestfn).c_str());
        if (plt==1) fprintf(f, "pause -1\n");

    }

    fclose(f);
    std::cout<<" done!\n";

    sprintf(fn, "%s%spsf.dat", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    f=fopen(fn, "w");
    long npoints=1000;
    for (  long i=-npoints; i<=npoints; i++) {
        double x=psfplot_xmax*double(i)/double(npoints);
        double y=psfplot_ymax*double(i)/double(npoints);
        double z=psfplot_zmax*double(i)/double(npoints);

        if (partner) {
            fprintf(f, "%15.10lf, %15.10lf, %15.10lf,  %15.10lf, %15.10lf, %15.10lf,  %15.10lf, %15.10lf, %15.10lf\n", x, illuminationEfficiency(x,0,0,expsf_r0,expsf_z0), detectionEfficiency(x,0,0), y, illuminationEfficiency(0,y,0,expsf_r0,expsf_z0), detectionEfficiency(0,y,0), z, illuminationEfficiency(0,0,z,expsf_r0,expsf_z0), detectionEfficiency(0,0,z),       illuminationEfficiency(x,0,0,expsf_r02,expsf_z02), illuminationEfficiency(0,y,0,expsf_r02,expsf_z02), illuminationEfficiency(0,0,z,expsf_r02,expsf_z02));
        } else {
            fprintf(f, "%15.10lf, %15.10lf, %15.10lf,  %15.10lf, %15.10lf, %15.10lf,  %15.10lf, %15.10lf, %15.10lf\n", x, illuminationEfficiency(x,0,0,expsf_r0,expsf_z0), detectionEfficiency(x,0,0), y, illuminationEfficiency(0,y,0,expsf_r0,expsf_z0), detectionEfficiency(0,y,0), z, illuminationEfficiency(0,0,z,expsf_r0,expsf_z0), detectionEfficiency(0,0,z));
        }
    }
    fclose(f);
    std::cout<<" done!\n";
    std::string psffn=fn;


    sprintf(fn, "%s%spsfxy.dat", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    f=fopen(fn, "w");
    long npointsi=100;
    for (  long i=-npointsi; i<=npointsi; i++) {
        for (  long j=-npointsi; j<=npointsi; j++) {
            double x=psfplot_xmax*double(j)/double(npointsi);
            double y=psfplot_ymax*double(i)/double(npointsi);

            fprintf(f, "%15.10lf ", illuminationEfficiency(x,y,0,expsf_r0,expsf_z0));
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n\n");
    std::cout<<" .";
    for (  long i=-npointsi; i<=npointsi; i++) {
        for (  long j=-npointsi; j<=npointsi; j++) {
            double x=psfplot_xmax*double(j)/double(npointsi);
            double z=psfplot_zmax*double(i)/double(npointsi);

            fprintf(f, "%15.10lf ", illuminationEfficiency(x,0,z,expsf_r0,expsf_z0));
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n\n");
    std::cout<<" .";
    for (  long i=-npointsi; i<=npointsi; i++) {
        for (  long j=-npointsi; j<=npointsi; j++) {
            double y=psfplot_ymax*double(j)/double(npointsi);
            double z=psfplot_zmax*double(i)/double(npointsi);

            fprintf(f, "%15.10lf ", illuminationEfficiency(0,y,z,expsf_r0,expsf_z0));
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n\n");
    std::cout<<" .";
    for (  long i=-npointsi; i<=npointsi; i++) {
        for (  long j=-npointsi; j<=npointsi; j++) {
            double x=psfplot_xmax*double(j)/double(npointsi);
            double y=psfplot_ymax*double(i)/double(npointsi);

            fprintf(f, "%15.10lf ", detectionEfficiency(x,y,0));
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n\n");
    std::cout<<" .";
    for (  long i=-npointsi; i<=npointsi; i++) {
        for (  long j=-npointsi; j<=npointsi; j++) {
            double x=psfplot_xmax*double(j)/double(npointsi);
            double z=psfplot_zmax*double(i)/double(npointsi);

            fprintf(f, "%15.10lf ", detectionEfficiency(x,0,z));
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n\n");
    std::cout<<" .";
    for (  long i=-npointsi; i<=npointsi; i++) {
        for (  long j=-npointsi; j<=npointsi; j++) {
            double y=psfplot_ymax*double(j)/double(npointsi);
            double z=psfplot_zmax*double(i)/double(npointsi);

            fprintf(f, "%15.10lf ", detectionEfficiency(0,y,z));
        }
        fprintf(f, "\n");
    }
    std::cout<<" .";
    fprintf(f, "\n\n");

    for (  long i=-npointsi; i<=npointsi; i++) {
        for (  long j=-npointsi; j<=npointsi; j++) {
            double x=psfplot_xmax*double(j)/double(npointsi);
            double y=psfplot_ymax*double(i)/double(npointsi);

            fprintf(f, "%15.10lf ", illuminationEfficiency(x,y,0,expsf_r0,expsf_z0)*detectionEfficiency(x,y,0));
        }
        fprintf(f, "\n");
    }
    std::cout<<" .";
    fprintf(f, "\n\n");
    for (  long i=-npointsi; i<=npointsi; i++) {
        for (  long j=-npointsi; j<=npointsi; j++) {
            double x=psfplot_xmax*double(j)/double(npointsi);
            double z=psfplot_zmax*double(i)/double(npointsi);

            fprintf(f, "%15.10lf ", illuminationEfficiency(x,0,z,expsf_r0,expsf_z0)*detectionEfficiency(x,0,z));
        }
        fprintf(f, "\n");
    }
    std::cout<<" .";
    fprintf(f, "\n\n");
    for (  long i=-npointsi; i<=npointsi; i++) {
        for (  long j=-npointsi; j<=npointsi; j++) {
            double y=psfplot_ymax*double(j)/double(npointsi);
            double z=psfplot_zmax*double(i)/double(npointsi);

            fprintf(f, "%15.10lf ", illuminationEfficiency(0,y,z,expsf_r0,expsf_z0)*detectionEfficiency(0,y,z));
        }
        fprintf(f, "\n");
    }
    std::cout<<" .";
    fprintf(f, "\n\n");
    fclose(f);
    std::cout<<" done!\n";
    std::string psfxyfn=fn;


    sprintf(fn, "%s%spsf.plt", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    f=fopen(fn, "w");
    fprintf(f, "npoints=%s\n", inttostr(npointsi).c_str());
    fprintf(f, "deltax=%lf\n", psfplot_xmax/double(npoints));
    fprintf(f, "deltay=%lf\n", psfplot_ymax/double(npoints));
    fprintf(f, "deltaz=%lf\n", psfplot_zmax/double(npoints));
    fprintf(f, "deltaxi=%lf\n", psfplot_xmax/double(npointsi));
    fprintf(f, "deltayi=%lf\n", psfplot_ymax/double(npointsi));
    fprintf(f, "deltazi=%lf\n", psfplot_zmax/double(npointsi));
    fprintf(f, "g(x,w)=exp(-2.0*x*x/w/w)\n");
    for (int plt=0; plt<2; plt++) {
        if (plt==0) {
            fprintf(f, "set terminal pdfcairo color solid font \"%s, 5\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
            fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"psf.pdf").c_str());
        } else if (plt==1) {
            fprintf(f, "set terminal wxt font \"%s, 5\"\n", GNUPLOT_FONT);
            fprintf(f, "set output\n");
        }

        fprintf(f, "set size noratio\n");
        for (int mp=0; mp<2; mp++) {
            if (mp==0) fprintf(f, "set multiplot layout 1,3\n");
            fprintf(f, "set title \"PSF: x-direction: %s\"\n", description.c_str());
            fprintf(f, "set xlabel \"position x [micron]\"\n");
            fprintf(f, "set ylabel \"intensity / detection probability [0..1]\"\n");
            fprintf(f, "set samples 1000\n");
            fprintf(f, "wx=%lf\n", detpsf_r0);
            fprintf(f, "wy=%lf\n", detpsf_z0);
            fprintf(f, "wz=%lf\n", detpsf_z0);
            fprintf(f, "fit Ax*g(x, wx) \"%s\" using 1:(($2)*($3)) via wx,Ax\n", extract_file_name(psffn).c_str());
            fprintf(f, "fit Ay*g(x, wy) \"%s\" using 4:(($5)*($6)) via wy,Ay\n", extract_file_name(psffn).c_str());
            fprintf(f, "fit Az*g(x, wz) \"%s\" using 7:(($8)*($9)) via wz,Az\n", extract_file_name(psffn).c_str());
            fprintf(f, "plot \"%s\" using 1:2 title \"illumination\" with lines, "
            "\"%s\" using 1:3 title \"detection\" with lines, "
            "\"%s\" using 1:(($2)*($3)) title \"ill * det\" with lines, "
            "Ax*g(x, wx) title sprintf('gaussian fit, w=%%f micron', wx) with lines"
            "\n", extract_file_name(psffn).c_str(), extract_file_name(psffn).c_str(), extract_file_name(psffn).c_str());
            if (mp==1 && plt==1) fprintf(f, "pause -1\n");
            fprintf(f, "set title \"PSF: y-direction: %s\"\n", description.c_str());
            fprintf(f, "set xlabel \"position y [micron]\"\n");
            fprintf(f, "set ylabel \"intensity / detection probability [0..1]\"\n");
            fprintf(f, "plot \"%s\" using 4:5 title \"illumination\" with lines, "
            "\"%s\" using 4:6 title \"detection\" with lines, "
            "\"%s\" using 4:(($5)*($6)) title \"ill * det\" with lines, "
            "Ay*g(x, wy) title sprintf('gaussian fit, w=%%f micron', wy) with lines"
            "\n", extract_file_name(psffn).c_str(), extract_file_name(psffn).c_str(), extract_file_name(psffn).c_str());
            if (mp==1 && plt==1) fprintf(f, "pause -1\n");
            fprintf(f, "set title \"PSF: z-direction: %s\"\n", description.c_str());
            fprintf(f, "set xlabel \"position z [micron]\"\n");
            fprintf(f, "set ylabel \"intensity / detection probability [0..1]\"\n");
            fprintf(f, "plot \"%s\" using 7:8 title \"illumination\" with lines, "
            "\"%s\" using 7:9 title \"detection\" with lines, "
            "\"%s\" using 7:(($8)*($9)) title \"ill * det\" with lines, "
            "Az*g(x, wz) title sprintf('gaussian fit, w=%%f micron', wz) with lines"
            "\n", extract_file_name(psffn).c_str(), extract_file_name(psffn).c_str(), extract_file_name(psffn).c_str());
            if (mp==0) fprintf(f, "unset multiplot\n");
            if (plt==1) fprintf(f, "pause -1\n");
        }
        for (int transc=0; transc<2; transc++) {
            std::string trans="3";
            if (transc==1) trans="(log($3)/log(10))";
            std::string title="";
            if (transc==1) title=" (log10-scale)";
            fprintf(f, "set multiplot layout 1,3\n");
            fprintf(f, "set xlabel \"position x [micron]\"\n");
            fprintf(f, "set ylabel \"position y [micron]\"\n");
            fprintf(f, "set size ratio deltayi/deltaxi\n");
            fprintf(f, "set title \"xy-cut: illumination %s%s\"\n", description.c_str(), title.c_str());
            fprintf(f, "plot [0:2*npoints*deltaxi] [0:2*npoints*deltayi] \"%s\" using (($1)*deltaxi):(($2)*deltayi):%s matrix index 0 notitle with image"
            "\n", extract_file_name(psfxyfn).c_str(), trans.c_str());
            fprintf(f, "set title \"xy-cut: detection %s%s\"\n", description.c_str(), title.c_str());
            fprintf(f, "plot [0:2*npoints*deltaxi] [0:2*npoints*deltayi] \"%s\" using (($1)*deltaxi):(($2)*deltayi):%s matrix index 3 notitle with image"
            "\n", extract_file_name(psfxyfn).c_str(), trans.c_str());
            fprintf(f, "set title \"xy-cut: ill*det %s%s\"\n", description.c_str(), title.c_str());
            fprintf(f, "plot [0:2*npoints*deltaxi] [0:2*npoints*deltayi] \"%s\" using (($1)*deltaxi):(($2)*deltayi):%s matrix index 6 notitle with image"
            "\n", extract_file_name(psfxyfn).c_str(), trans.c_str());
            fprintf(f, "unset multiplot\n");
            if (plt==1) fprintf(f, "pause -1\n");

            fprintf(f, "set multiplot layout 1,3\n");
            fprintf(f, "set xlabel \"position x [micron]\"\n");
            fprintf(f, "set ylabel \"position z [micron]\"\n");
            fprintf(f, "set size ratio deltazi/deltaxi\n");
            fprintf(f, "set title \"xz-cut: illumination %s%s\"\n", description.c_str(), title.c_str());
            fprintf(f, "plot [0:2*npoints*deltaxi] [0:2*npoints*deltazi] \"%s\" using (($1)*deltaxi):(($2)*deltazi):%s matrix index 1 notitle with image"
            "\n", extract_file_name(psfxyfn).c_str(), trans.c_str());
            fprintf(f, "set title \"xz-cut: detection %s%s\"\n", description.c_str(), title.c_str());
            fprintf(f, "plot [0:2*npoints*deltaxi] [0:2*npoints*deltazi] \"%s\" using (($1)*deltaxi):(($2)*deltazi):%s matrix index 4 notitle with image"
            "\n", extract_file_name(psfxyfn).c_str(), trans.c_str());
            fprintf(f, "set title \"xz-cut: ill*det %s%s\"\n", description.c_str(), title.c_str());
            fprintf(f, "plot [0:2*npoints*deltaxi] [0:2*npoints*deltazi] \"%s\" using (($1)*deltaxi):(($2)*deltazi):%s matrix index 7 notitle with image"
            "\n", extract_file_name(psfxyfn).c_str(), trans.c_str());
            fprintf(f, "unset multiplot\n");
            if (plt==1) fprintf(f, "pause -1\n");

            fprintf(f, "set multiplot layout 1,3\n");
            fprintf(f, "set xlabel \"position y [micron]\"\n");
            fprintf(f, "set ylabel \"position z [micron]\"\n");
            fprintf(f, "set size ratio deltazi/deltayi\n");
            fprintf(f, "set title \"yz-cut: illumination %s%s\"\n", description.c_str(), title.c_str());
            fprintf(f, "plot [0:2*npoints*deltayi] [0:2*npoints*deltazi] \"%s\" using (($1)*deltayi):(($2)*deltazi):%s matrix index 2 notitle with image"
            "\n", extract_file_name(psfxyfn).c_str(), trans.c_str());
            fprintf(f, "set title \"yz-cut: detection %s%s\"\n", description.c_str(), title.c_str());
            fprintf(f, "plot [0:2*npoints*deltayi] [0:2*npoints*deltazi] [0:npoints*deltazi] \"%s\" using (($1)*deltayi):(($2)*deltazi):%s matrix index 5 notitle with image"
            "\n", extract_file_name(psfxyfn).c_str(), trans.c_str());
            fprintf(f, "set title \"yz-cut: ill*det %s%s\"\n", description.c_str(), title.c_str());
            fprintf(f, "plot [0:2*npoints*deltayi] [0:2*npoints*deltazi] [0:npoints*deltazi] \"%s\" using (($1)*deltayi):(($2)*deltazi):%s matrix index 8 notitle with image"
            "\n", extract_file_name(psfxyfn).c_str(), trans.c_str());
            fprintf(f, "unset multiplot\n");
            if (plt==1) fprintf(f, "pause -1\n");
        }
    }

    fclose(f);
    std::cout<<" done!\n";

    if (!online_correlation) {
        if (save_timeseries) {
            sprintf(fn, "%s%sts.dat", basename.c_str(), object_name.c_str());
            std::cout<<"writing '"<<fn<<"' ...";
            f=fopen(fn, "w");
            double t=0;
            for (unsigned long long i=0; i<timesteps; i++) {
                fprintf(f, "%15.10lf, %d\n", t, timeseries[i]);
                t=t+corr_taumin;
            }
            fclose(f);
            std::string tsfn=fn;
            std::cout<<" done!\n";

            sprintf(fn, "%s%stsplot.plt", basename.c_str(), object_name.c_str());
            std::cout<<"writing '"<<fn<<"' ...";
            f=fopen(fn, "w");
            for (int plt=0; plt<2; plt++) {
                if (plt==0) {
                    fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                    fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"tsplot.pdf").c_str());
                } else if (plt==1) {
                    fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
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
            std::cout<<" done!\n";
        }
    }

    double* bts_time=NULL;
    bool free_btstime=false;
    double* bts_1=NULL;
    bool free_bts1=false;
    double* bts_2=NULL;
    bool free_bts2=false;
    int bts_N=0;

    if (save_binning) {
        if (online_correlation) {
            // in online-correlation-mode no complete timeseries is built up, only a binned one, so we save that
            sprintf(fn, "%s%sbts.dat", basename.c_str(), object_name.c_str());
            std::cout<<"writing '"<<fn<<"' ...";
            f=fopen(fn, "w");
            double t=0;
            int b=round(save_binning_time/corr_taumin);
            bts_N=timesteps/b-1;
            bts_time=(double*)calloc(bts_N+20, sizeof(double));
            free_btstime=true;
            bts_1=(double*)calloc(bts_N+20, sizeof(double));
            free_bts1=true;
            for (unsigned long long i=0; i<timesteps/b-1; i++) {
                fprintf(f, "%15.10lf, %d\n", t, binned_timeseries[i]);
                if (bts_time) bts_time[i]=t;
                if (bts_1) bts_1[i]=binned_timeseries[i];
                t=t+save_binning_time;
            }
            fclose(f);
            std::string tsfn=fn;
            std::string tsotherfn="";
            if (partner) {
                char fn1[255];
                sprintf(fn1,"%s%sbts.dat", basename.c_str(), partner->get_object_name().c_str());
                tsotherfn=fn1;
                if (bts_N>0 && partner->save_binning && partner->online_correlation && save_binning_time== partner->save_binning_time && corr_taumin==partner->corr_taumin && timesteps==partner->timesteps) {
                    bts_2=(double*)calloc(bts_N+20, sizeof(double));
                    free_bts2=true;
                    if (bts_2) for (unsigned long long i=0; i<timesteps/b-1; i++) {
                        bts_2[i]=partner->binned_timeseries[i];
                    }
                }
            }
            std::cout<<" done!\n";

            sprintf(fn, "%s%sbtsplot.plt", basename.c_str(), object_name.c_str());
            std::cout<<"writing '"<<fn<<"' ...";
            f=fopen(fn, "w");
            for (int plt=0; plt<2; plt++) {
                if (plt==0) {
                    fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                    fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"btsplot.pdf").c_str());
                } else if (plt==1) {
                    fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
                    fprintf(f, "set output\n");
                }
                fprintf(f, "set xlabel \"time [seconds]\"\n");
                fprintf(f, "set ylabel \"photon count [photons/%lfsec]\"\n", corr_taumin*b);
                fprintf(f, "set title \"object description: %s\"\n", description.c_str());
                if (partner) {
                    fprintf(f, "plot \"%s\" title \"CR1: %s\" with steps lc rgb \"dark-green\", "
                                    "\"%s\" title \"CR2: %s\" with steps lc rgb \"dark-red\"\n", extract_file_name(tsfn).c_str(), object_name.c_str(), extract_file_name(tsotherfn).c_str(), partner->get_object_name().c_str());
                } else {
                    fprintf(f, "plot \"%s\" with steps\n", extract_file_name(tsfn).c_str());
                }
                if (plt==1) fprintf(f, "pause -1\n");
                fprintf(f, "set xlabel \"time [seconds]\"\n");
                fprintf(f, "set ylabel \"photon count [kcps]\"\n");
                fprintf(f, "set title \"object description: %s\"\n", description.c_str());
                if (partner) {
                    fprintf(f, "plot \"%s\" using 1:(($2)/%lf/1000.0) title \"CR1: %s\" with steps lc rgb \"dark-green\", "
                                    "\"%s\" using 1:(($2)/%lf/1000.0) title \"CR2: %s\" with steps lc rgb \"dark-red\"\n", extract_file_name(tsfn).c_str(), corr_taumin*double(b), object_name.c_str(), extract_file_name(tsotherfn).c_str(), corr_taumin*double(b), partner->get_object_name().c_str());
                } else {
                    fprintf(f, "plot \"%s\" using 1:(($2)/%lf/1000.0) with steps\n", extract_file_name(tsfn).c_str(), corr_taumin*b);
                }
                if (plt==1) fprintf(f, "pause -1\n");



            }
            fclose(f);
            std::cout<<" done!\n";
        } else {
            // if we are not in online-correlation mode we may simply calculate the binned timeseries from the
            // complete timeseries
            sprintf(fn, "%s%sbts.dat", basename.c_str(), object_name.c_str());
            std::cout<<"writing '"<<fn<<"' ...";
            f=fopen(fn, "w");
            double t=0;
            unsigned long long b=round(save_binning_time/corr_taumin);
            if (timesteps>b) {
                bts_N=timesteps/b-1;
                if (bts_N>0) {
                    bts_time=(double*)calloc(bts_N+20, sizeof(double));
                    free_btstime=true;
                    bts_1=(double*)calloc(bts_N+20, sizeof(double));
                    free_bts1=true;
                }
                int ii=0;
                for (unsigned long long i=0; i<timesteps-b; i=i+b) {
                    register long int  ts=0;
                    for (unsigned long long j=0; j<b; j++) {
                        ts=ts+timeseries[i+j];
                    }
                    fprintf(f, "%15.10lf, %ld\n", t, ts);
                    if (bts_time) bts_time[ii]=t;
                    if (bts_1) bts_1[ii]=ts;
                    //fprintf(stdout, "%15.10lf, %lu\n", t, ts);
                    t=t+corr_taumin*b;
                    ii++;
                }

                if (bts_N>0 && partner && partner->save_binning && !partner->online_correlation && save_binning_time== partner->save_binning_time && corr_taumin==partner->corr_taumin && timesteps==partner->timesteps) {
                    bts_2=(double*)calloc(bts_N+20, sizeof(double));
                    free_bts2=true;
                    ii=0;
                    for (unsigned long long i=0; i<timesteps-b; i=i+b) {
                        register long int  ts=0;
                        for (unsigned long long j=0; j<b; j++) {
                            ts=ts+partner->timeseries[i+j];
                        }
                        if (bts_2) bts_2[ii]=ts;
                        ii++;
                    }
                }
            }
            fclose(f);
            std::string tsfn=fn;
            std::string tsotherfn="";
            if (partner) {
                char fn1[255];
                sprintf(fn1,"%s%sbts.dat", basename.c_str(), partner->get_object_name().c_str());
                tsotherfn=fn1;
            }
            std::cout<<" done!\n";
            sprintf(fn, "%s%sbtsplot.plt", basename.c_str(), object_name.c_str());
            std::cout<<"writing '"<<fn<<"' ...";
            f=fopen(fn, "w");
            fprintf(f, "dt_corr=%15.10lf\n", corr_taumin);
            fprintf(f, "dt=%15.10lf\n", save_binning_time);
            for (int plt=0; plt<2; plt++) {
                //std::cout<<plt<<" 1\n";
                if (plt==0) {
                    fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                    fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"btsplot.pdf").c_str());
                } else if (plt==1) {
                    fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
                    fprintf(f, "set output\n");
                }
                //std::cout<<plt<<" 2\n";
                fprintf(f, "set xlabel \"time [seconds]\"\n");
                fprintf(f, "set ylabel \"photon count [photons/%lfus]\"\n", corr_taumin*1e6);
                fprintf(f, "set title \"counts per corrtaumin, object description: %s\"\n", description.c_str());
                if (partner) {
                    fprintf(f, "plot \"%s\" using 1:($2*dt_corr/dt) title \"CR1: %s\" with steps lc rgb \"dark-green\", "
                                    "\"%s\" using 1:($2*dt_corr/dt) title \"CR2: %s\" with steps lc rgb \"dark-red\"\n", extract_file_name(tsfn).c_str(), object_name.c_str(), extract_file_name(tsotherfn).c_str(), partner->get_object_name().c_str());
                } else {
                    fprintf(f, "plot \"%s\" using 1:($2*dt_corr/dt) with steps\n", extract_file_name(tsfn).c_str());
                }
                //std::cout<<plt<<" 3\n";

                if (plt==1) fprintf(f, "pause -1\n");
                fprintf(f, "set xlabel \"time [seconds]\"\n");
                fprintf(f, "set ylabel \"photon count [photons/%lfms]\"\n", corr_taumin*b*1e3);
                fprintf(f, "set title \"counts per bin time, object description: %s\"\n", description.c_str());
                if (partner) {
                    fprintf(f, "plot \"%s\" title \"CR1: %s\" with steps lc rgb \"dark-green\", "
                                    "\"%s\" title \"CR2: %s\" with steps lc rgb \"dark-red\"\n", extract_file_name(tsfn).c_str(), object_name.c_str(), extract_file_name(tsotherfn).c_str(), partner->get_object_name().c_str());
                } else {
                    fprintf(f, "plot \"%s\" with steps\n", extract_file_name(tsfn).c_str());
                }
                //std::cout<<plt<<" 4\n";

                if (plt==1) fprintf(f, "pause -1\n");
                fprintf(f, "set xlabel \"time [seconds]\"\n");
                fprintf(f, "set ylabel \"photon count [kcps]\"\n");
                fprintf(f, "set title \"counts, object description: %s\"\n", description.c_str());
                if (partner) {
                    fprintf(f, "plot \"%s\" using 1:(($2)/%lf/1000.0) title \"CR1: %s\" with steps lc rgb \"dark-green\", "
                                    "\"%s\" using 1:(($2)/%lf/1000.0) title \"CR2: %s\" with steps lc rgb \"dark-red\"\n", extract_file_name(tsfn).c_str(), corr_taumin*double(b), object_name.c_str(), extract_file_name(tsotherfn).c_str(), corr_taumin*double(b), partner->get_object_name().c_str());
                } else {
                    fprintf(f, "plot \"%s\" using 1:(($2)/%lf/1000.0) with steps\n", extract_file_name(tsfn).c_str(), corr_taumin*double(b));
                }
                //std::cout<<plt<<" 5\n";

                if (plt==1) fprintf(f, "pause -1\n");
            }
            fclose(f);
            std::cout<<" done!\n";
        }
    }











    sprintf(fn, "%s%salvacf.asc", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    f=fopen(fn, "w");
    std::string alvmode="SINGLE AUTO CH0";
    double* corr2=NULL;
    if (partner && timesteps==partner->timesteps && corr_taumin==partner->corr_taumin && S==partner->S && P==partner->P && m==partner->m) {
        alvmode="DUAL AUTO CH0";
        corr2=partner->corr;
    }
    alv5000WriteHeader(f, description, duration, alvmode.c_str(), lambda_ex, 0, 0);
    alv5000WriteCorrelation(f, istart, slots, corr_tau, corr, corr2);
    alv5000WriteCountrate(f, bts_N, bts_time, bts_1, bts_2, true);
    fclose(f);
    std::cout<<" done!\n";

    if (partner && timesteps==partner->timesteps && corr_taumin==partner->corr_taumin && S==partner->S && P==partner->P && m==partner->m) {
        sprintf(fn, "%s%salvccf.asc", basename.c_str(), object_name.c_str());
        std::cout<<"writing '"<<fn<<"' ...";
        f=fopen(fn, "w");
        alvmode="SINGLE CROSS CH0";
        alv5000WriteHeader(f, description, duration, alvmode.c_str(), lambda_ex, 0, 0);
        alv5000WriteCorrelation(f, istart, slots, corr_tau, corr_fccs);
        alv5000WriteCountrate(f, bts_N, bts_time, bts_1, bts_2, true);
        fclose(f);
        std::cout<<" done!\n";
    }
    if (free_btstime && bts_time) {
        free(bts_time);
        bts_time=NULL;
    }
    if (free_bts1 && bts_1) {
        free(bts_1);
        bts_1=NULL;
    }
    if (free_bts2 && bts_2) {
        free(bts_2);
        bts_2=NULL;
    }
    bts_N=0;



   if (save_arrivaltimes) {

        // if we are not in online-correlation mode we may simply calculate the binned timeseries from the
        // complete timeseries
        sprintf(fn, "%s%sarrivals.dat", basename.c_str(), object_name.c_str());
        std::cout<<"writing '"<<fn<<"' ...";
        f=fopen(fn, "w");
        for (size_t i=0; i<arrival_times.size(); i++) {
            fprintf(f, "%15.10lf\n", arrival_times[i]);
        }
        fclose(f);
        std::string atfn=fn;
        std::cout<<" done!\n";
        sprintf(fn, "%s%sarrivalsplot.plt", basename.c_str(), object_name.c_str());
        std::cout<<"writing '"<<fn<<"' ...";
        f=fopen(fn, "w");
        for (int plt=0; plt<2; plt++) {
            if (plt==0) {
                fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"sarrivalsplot.pdf").c_str());
            } else if (plt==1) {
                fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
                fprintf(f, "set output\n");
            }
            fprintf(f, "set xlabel \"time [seconds]\"\n");
            fprintf(f, "set ylabel \"photon count [photons/%lfus]\"\n", corr_taumin*1e6);
            fprintf(f, "set title \"photons per corrtaumin, object description: %s\"\n", description.c_str());
            fprintf(f, "plot \"%s\" using 1:(1.0) with impulses\n", extract_file_name(atfn).c_str());
            if (plt==1) fprintf(f, "pause -1\n");
        }
        fclose(f);
        std::cout<<" done!\n";
    }

    sprintf(fn, "%s%sdetpsf.tif", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    psf_rz_image_detection.save_tinytifffloat(fn);
    std::cout<<" done!\n";




    sprintf(fn, "%s%sdetill.tif", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    psf_rz_image_illumination.save_tinytifffloat(fn);
    std::cout<<" done!\n";


    sprintf(fn, "%s%sdetpsfi.tif", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    psf_rz_image_detection.save_tinytiffuint16scaled(fn);
    std::cout<<" done!\n";


    sprintf(fn, "%s%sdetilli.tif", basename.c_str(), object_name.c_str());
    std::cout<<"writing '"<<fn<<"' ...";
    psf_rz_image_illumination.save_tinytiffuint16scaled(fn);
    std::cout<<" done!\n";


    std::vector<std::vector<std::string> > plot_withs;
    if (plot_with.size()>0) {
        plot_withs.push_back(plot_with);
        std::map<int, std::vector<std::string> >::iterator it;
        for (it=plot_with_more.begin(); it!=plot_with_more.end(); ++it) {
            if (it->second.size()>0) {
                plot_withs.push_back(it->second);
            }
        }
    }

    for (int iip=0; iip<plot_withs.size(); iip++) {
        std::vector<std::string> plot_with=plot_withs[iip];
        if (plot_with.size()>0) {

            sprintf(fn, "%s%scorrplotwith%d_simple.plt", basename.c_str(), object_name.c_str(), iip);
            std::cout<<"writing '"<<fn<<"' ...";
            f=fopen(fn, "w");
            for (int plt=0; plt<2; plt++) {
                if (plt==0) {
                    fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                    fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"corrplotwith"+inttostr(iip)+"_simple.pdf").c_str());
                } else if (plt==1) {
                    fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
                    fprintf(f, "set output\n");
                }

                fprintf(f, "set logscale x\n");
                fprintf(f, "set title \"object description: %s\"\n", description.c_str());
                fprintf(f, "plot \"%s\" title \"simulation data: %s (%s)\" with points", extract_file_name(corrfn).c_str(), description.c_str(), object_name.c_str());
                for (size_t i=0; i<plot_with.size(); i++) {
                    std::string pw=strstrip(plot_with[i]);
                    if (measmap.count(pw)>0) {
                        FCSMeasurement* fcspw=dynamic_cast<FCSMeasurement*>(measmap[pw]);
                        char fn2[1024];
                        sprintf(fn2, "%s%scorr.dat", basename.c_str(), pw.c_str());
                        if (fcspw) {
                            fprintf(f, ", \\\n    \"%s\" title \"simulation data: %s (%s)\" with points", extract_file_name(std::string(fn2)).c_str(),fcspw->get_description().c_str(), pw.c_str());
                        }
                    }
                }
                fprintf(f, "\n");
            }
            fprintf(f, "pause -1\n");
            fclose(f);
            std::cout<<" done!\n";



            sprintf(fn, "%s%spsfplotwith%d.plt", basename.c_str(), object_name.c_str(), iip);
            std::cout<<"writing '"<<fn<<"' ...";
            f=fopen(fn, "w");
            fprintf(f, "npoints=%s\n", inttostr(npointsi).c_str());
            fprintf(f, "deltax=%lf\n", psfplot_xmax/double(npoints));
            fprintf(f, "deltay=%lf\n", psfplot_ymax/double(npoints));
            fprintf(f, "deltaz=%lf\n", psfplot_zmax/double(npoints));
            fprintf(f, "deltaxi=%lf\n", psfplot_xmax/double(npointsi));
            fprintf(f, "deltayi=%lf\n", psfplot_ymax/double(npointsi));
            fprintf(f, "deltazi=%lf\n", psfplot_zmax/double(npointsi));
            fprintf(f, "g(x,w)=exp(-2.0*x*x/w/w)\n");
            std::string plotx, ploty, plotz;
            for (int i=0; i<int(plot_with.size()); i++) {
                std::string pw=strstrip(plot_with[i]);
                if (measmap.count(pw)>0) {
                    FCSMeasurement* fcspw=dynamic_cast<FCSMeasurement*>(measmap[pw]);
                    char fn2[1024];
                    char plt[8192];
                    sprintf(fn2, "%s%spsf.dat", basename.c_str(), pw.c_str());
                    if (fcspw) {
                        fprintf(f, "fit Ax%d*g(x, wx%d) \"%s\" using 1:(($2)*($3)) via wx%d, Ax%d\n", i, i, extract_file_name(std::string(fn2)).c_str(), i, i);
                        fprintf(f, "fit Ay%d*g(x, wy%d) \"%s\" using 4:(($5)*($6)) via wy%d, Ay%d\n", i, i, extract_file_name(std::string(fn2)).c_str(), i, i);
                        fprintf(f, "fit Az%d*g(x, wz%d) \"%s\" using 7:(($8)*($9)) via wz%d, Az%d\n", i, i, extract_file_name(std::string(fn2)).c_str(), i, i);

                        sprintf(plt, ",\\\n  \"%s\" using 1:2 title \"%s: illumination\" with lines, "
                                     "  \"%s\" using 1:3 title \"%s: detection\" with lines, "
                                     "  \"%s\" using 1:(($2)*($3)) title \"%s: ill * det\" with lines, "
                                     "  Ax%d*g(x, wx%d) title sprintf('%s: gaussian fit, w=%%f micron', wx%d) with lines"
                                     "", extract_file_name(std::string(fn2)).c_str(), fcspw->get_description().c_str(), extract_file_name(std::string(fn2)).c_str(), fcspw->get_description().c_str(), extract_file_name(std::string(fn2)).c_str(), fcspw->get_description().c_str(), i, i, fcspw->get_description().c_str(), i);
                        plotx=plotx+std::string(plt);
                        sprintf(plt, ",\\\n  \"%s\" using 4:5 title \"%s: illumination\" with lines, "
                                     "  \"%s\" using 4:6 title \"%s: detection\" with lines, "
                                     "  \"%s\" using 4:(($5)*($6)) title \"%s: ill * det\" with lines, "
                                     "  Ay%d*g(x, wy%d) title sprintf('%s: gaussian fit, w=%%f micron', wy%d) with lines"
                                     "", extract_file_name(std::string(fn2)).c_str(), fcspw->get_description().c_str(), extract_file_name(std::string(fn2)).c_str(), fcspw->get_description().c_str(), extract_file_name(std::string(fn2)).c_str(), fcspw->get_description().c_str(), i, i, fcspw->get_description().c_str(), i);
                        ploty=ploty+std::string(plt);
                        sprintf(plt, ",\\\n  \"%s\" using 7:8 title \"%s: illumination\" with lines, "
                                     "  \"%s\" using 7:9 title \"%s: detection\" with lines, "
                                     "  \"%s\" using 7:(($8)*($9)) title \"%s: ill * det\" with lines, "
                                     "  Az%d*g(x, wz%d) title sprintf('%s: gaussian fit, w=%%f micron', wz%d) with lines"
                                     "", extract_file_name(std::string(fn2)).c_str(), fcspw->get_description().c_str(), extract_file_name(std::string(fn2)).c_str(), fcspw->get_description().c_str(), extract_file_name(std::string(fn2)).c_str(), fcspw->get_description().c_str(), i, i, fcspw->get_description().c_str(), i);
                        plotz=plotz+std::string(plt);
                    }
                }
            }


            for (int plt=0; plt<2; plt++) {
                if (plt==0) {
                    fprintf(f, "set terminal pdfcairo color solid font \"%s, 5\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                    fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"psfplotwith"+inttostr(iip)+".pdf").c_str());
                } else if (plt==1) {
                    fprintf(f, "set terminal wxt font \"%s, 5\"\n", GNUPLOT_FONT);
                    fprintf(f, "set output\n");
                }

                fprintf(f, "set size noratio\n");
                for (int mp=0; mp<2; mp++) {
                    if (mp==0) fprintf(f, "set multiplot layout 1,3\n");
                    fprintf(f, "set title \"PSF: x-direction: %s\"\n", description.c_str());
                    fprintf(f, "set xlabel \"position x [micron]\"\n");
                    fprintf(f, "set ylabel \"intensity / detection probability [0..1]\"\n");
                    fprintf(f, "set samples 1000\n");
                    fprintf(f, "wx=%lf\n", detpsf_r0);
                    fprintf(f, "wy=%lf\n", detpsf_z0);
                    fprintf(f, "wz=%lf\n", detpsf_z0);



                    fprintf(f, "fit g(x, wx) \"%s\" using 1:(($2)*($3)) via wx\n", extract_file_name(psffn).c_str());
                    fprintf(f, "fit g(x, wy) \"%s\" using 4:(($5)*($6)) via wy\n", extract_file_name(psffn).c_str());
                    fprintf(f, "fit g(x, wz) \"%s\" using 7:(($8)*($9)) via wz\n", extract_file_name(psffn).c_str());
                    fprintf(f, "plot \"%s\" using 1:2 title \"illumination\" with lines, "
                    "\"%s\" using 1:3 title \"detection\" with lines, "
                    "\"%s\" using 1:(($2)*($3)) title \"ill * det\" with lines, "
                    "g(x, wx) title sprintf('gaussian fit, w=%%f micron', wx) with lines"
                    , extract_file_name(psffn).c_str(), extract_file_name(psffn).c_str(), extract_file_name(psffn).c_str());
                    fprintf(f, "%s\n\n", plotx.c_str());
                    if (mp==1 && plt==1) fprintf(f, "pause -1\n");
                    fprintf(f, "set title \"PSF: y-direction: %s\"\n", description.c_str());
                    fprintf(f, "set xlabel \"position y [micron]\"\n");
                    fprintf(f, "set ylabel \"intensity / detection probability [0..1]\"\n");
                    fprintf(f, "plot \"%s\" using 4:5 title \"illumination\" with lines, "
                    "\"%s\" using 4:6 title \"detection\" with lines, "
                    "\"%s\" using 4:(($5)*($6)) title \"ill * det\" with lines, "
                    "g(x, wy) title sprintf('gaussian fit, w=%%f micron', wy) with lines"
                    , extract_file_name(psffn).c_str(), extract_file_name(psffn).c_str(), extract_file_name(psffn).c_str());
                    fprintf(f, "%s\n\n", plotx.c_str());
                    if (mp==1 && plt==1) fprintf(f, "pause -1\n");
                    fprintf(f, "set title \"PSF: z-direction: %s\"\n", description.c_str());
                    fprintf(f, "set xlabel \"position z [micron]\"\n");
                    fprintf(f, "set ylabel \"intensity / detection probability [0..1]\"\n");
                    fprintf(f, "plot \"%s\" using 7:8 title \"illumination\" with lines, "
                    "\"%s\" using 7:9 title \"detection\" with lines, "
                    "\"%s\" using 7:(($8)*($9)) title \"ill * det\" with lines, "
                    "g(x, wz) title sprintf('gaussian fit, w=%%f micron', wz) with lines"
                    , extract_file_name(psffn).c_str(), extract_file_name(psffn).c_str(), extract_file_name(psffn).c_str());
                    fprintf(f, "%s\n\n", plotx.c_str());
                    if (mp==0) fprintf(f, "unset multiplot\n");
                    if (plt==1) fprintf(f, "pause -1\n");
                }
            }

            fclose(f);
            std::cout<<" done!\n";




               if (save_arrivaltimes) {
                    sprintf(fn, "%s%sarrivals.dat", basename.c_str(), object_name.c_str());
                    std::string atfn=fn;
                    std::cout<<" done!\n";
                    sprintf(fn, "%s%sarrivalsplotwith%d.plt", basename.c_str(), object_name.c_str(), iip);
                    std::cout<<"writing '"<<fn<<"' ...";
                    f=fopen(fn, "w");
                    for (int plt=0; plt<2; plt++) {
                        if (plt==0) {
                            fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                            fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"sarrivalsplotwith"+inttostr(iip)+".pdf").c_str());
                        } else if (plt==1) {
                            fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
                            fprintf(f, "set output\n");
                        }
                        fprintf(f, "set xlabel \"time [seconds]\"\n");
                        fprintf(f, "set ylabel \"photon count [photons/%lfus]\"\n", corr_taumin*1e6);
                        fprintf(f, "set title \"photons per corrtaumin, object description: %s\"\n", description.c_str());
                        float impheight=1;
                        fprintf(f, "plot \"%s\" using 1:(1.0) title \"%s\" with impulses", extract_file_name(atfn).c_str(), get_object_name().c_str());
                        for (int i=0; i<int(plot_with.size()); i++) {
                            std::string pw=strstrip(plot_with[i]);
                            if (measmap.count(pw)>0) {
                                FCSMeasurement* fcspw=dynamic_cast<FCSMeasurement*>(measmap[pw]);
                                char fn2[1024];
                                char plt[8192];
                                sprintf(fn2, "%s%sarrivals.dat", basename.c_str(), pw.c_str());
                                std::string atfn2=fn2;
                                if (fcspw) {
                                    impheight=-1.0*impheight;
                                    fprintf(f, ",\\\n \"%s\" using 1:(%f) title \"%s\" with impulses", extract_file_name(atfn2).c_str(), impheight, fcspw->get_object_name().c_str());
                                }
                            }
                        }
                       fprintf(f, "\n\n");
                       if (plt==1) fprintf(f, "pause -1\n");
                    }
                    fclose(f);
                    std::cout<<" done!\n";
                }








                if (save_binning) {
                    unsigned long long b=round(save_binning_time/corr_taumin);

                    // if we are not in online-correlation mode we may simply calculate the binned timeseries from the
                    // complete timeseries
                    sprintf(fn, "%s%sbts.dat", basename.c_str(), object_name.c_str());
                    std::cout<<"writing '"<<fn<<"' ...";
                    std::string tsfn=fn;
                    std::cout<<" done!\n";
                    sprintf(fn, "%s%sbtsplotwith%d.plt", basename.c_str(), object_name.c_str(), iip);
                    std::cout<<"writing '"<<fn<<"' ...";
                    f=fopen(fn, "w");
                    fprintf(f, "dt_corr=%15.10lf\n", corr_taumin);
                    fprintf(f, "dt=%15.10lf\n", save_binning_time);
                    for (int plt=0; plt<2; plt++) {
                        if (plt==0) {
                            fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                            fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"btsplotwith"+inttostr(iip)+".pdf").c_str());
                        } else if (plt==1) {
                            fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
                            fprintf(f, "set output\n");
                        }
                        fprintf(f, "set xlabel \"time [seconds]\"\n");
                        fprintf(f, "set ylabel \"photon count [photons/%lfus]\"\n", corr_taumin*1e6);
                        fprintf(f, "set title \"counts per corrtaumin, object description: %s\"\n", description.c_str());
                        fprintf(f, "plot \"%s\" using 1:($2*dt_corr/dt) title \"%s\" with steps", extract_file_name(tsfn).c_str(), get_object_name().c_str());
                        for (int i=0; i<int(plot_with.size()); i++) {
                            std::string pw=strstrip(plot_with[i]);
                            if (measmap.count(pw)>0) {
                                FCSMeasurement* fcspw=dynamic_cast<FCSMeasurement*>(measmap[pw]);
                                char fn2[1024];
                                char plt[8192];
                                sprintf(fn2, "%s%sbts.dat", basename.c_str(), pw.c_str());
                                std::string atfn2=fn2;
                                if (fcspw) {
                                    fprintf(f, ",\\\n \"%s\" using 1:($2*dt_corr/dt) title \"%s\" with steps", extract_file_name(atfn2).c_str(), fcspw->get_object_name().c_str());
                                }
                            }
                        }
                        fprintf(f, "\n\n");
                        if (plt==1) fprintf(f, "pause -1\n");
                        fprintf(f, "set xlabel \"time [seconds]\"\n");
                        fprintf(f, "set ylabel \"photon count [photons/%lfms]\"\n", corr_taumin*b*1e3);
                        fprintf(f, "set title \"counts per bin time, object description: %s\"\n", description.c_str());
                        fprintf(f, "plot \"%s\" title \"%s\" with steps", extract_file_name(tsfn).c_str(), get_object_name().c_str());
                        for (int i=0; i<int(plot_with.size()); i++) {
                            std::string pw=strstrip(plot_with[i]);
                            if (measmap.count(pw)>0) {
                                FCSMeasurement* fcspw=dynamic_cast<FCSMeasurement*>(measmap[pw]);
                                char fn2[1024];
                                char plt[8192];
                                sprintf(fn2, "%s%sbts.dat", basename.c_str(), pw.c_str());
                                std::string atfn2=fn2;
                                if (fcspw) {
                                    fprintf(f, ",\\\n \"%s\" title \"%s\" with steps", extract_file_name(atfn2).c_str(), fcspw->get_object_name().c_str());
                                }
                            }
                        }
                        fprintf(f, "\n\n");
                        if (plt==1) fprintf(f, "pause -1\n");
                        fprintf(f, "set xlabel \"time [seconds]\"\n");
                        fprintf(f, "set ylabel \"photon count [kcps]\"\n");
                        fprintf(f, "set title \"counts, object description: %s\"\n", description.c_str(), get_object_name().c_str());
                        fprintf(f, "plot \"%s\" using 1:(($2)/%lf/1000.0) title \"%s\" with steps", extract_file_name(tsfn).c_str(), corr_taumin*b, get_object_name().c_str());
                                                for (int i=0; i<int(plot_with.size()); i++) {
                            std::string pw=strstrip(plot_with[i]);
                            if (measmap.count(pw)>0) {
                                FCSMeasurement* fcspw=dynamic_cast<FCSMeasurement*>(measmap[pw]);
                                char fn2[1024];
                                char plt[8192];
                                sprintf(fn2, "%s%sbts.dat", basename.c_str(), pw.c_str());
                                std::string atfn2=fn2;
                                if (fcspw) {
                                    fprintf(f, ",\\\n \"%s\" using 1:(($2)/%lf/1000.0) title \"%s\" with steps", extract_file_name(atfn2).c_str(), corr_taumin*b, fcspw->get_object_name().c_str());
                                }
                            }
                        }
                        fprintf(f, "\n\n");

                        if (plt==1) fprintf(f, "pause -1\n");
                    }
                    fclose(f);
                    std::cout<<" done!\n";
                }



        }
    }

}

std::string FCSMeasurement::report(){
    std::string s=FluorescenceMeasurement::report();
    s+="pos_laser     = ["+floattostr(ex_x0)+", "+floattostr(ex_y0)+", "+floattostr(ex_z0)+"] um\n";
    s+="distance_laser_detector = ["+floattostr(img_x0-ex_x0)+", "+floattostr(img_y0-ex_y0)+", "+floattostr(img_z0-ex_z0)+"] um\n";
    s+="                        = "+floattostr(sqrt(gsl_pow_2(img_x0-ex_x0)+gsl_pow_2(img_y0-ex_y0)+gsl_pow_2(img_z0-ex_z0)))+" um\n";
    s+="expsf_r0 = "+floattostr(expsf_r0)+" microns\n";
    s+="expsf_z0 = "+floattostr(expsf_z0)+" microns\n";
    if (lambda_ex2>0) {
        s+="pos_laser2    = ["+floattostr(ex_x02)+", "+floattostr(ex_y02)+", "+floattostr(ex_z02)+"] um\n";
        s+="distance_laser2_detector = ["+floattostr(img_x0-ex_x02)+", "+floattostr(img_y0-ex_y02)+", "+floattostr(img_z0-ex_z02)+"] um\n";
        s+="                        = "+floattostr(sqrt(gsl_pow_2(img_x0-ex_x02)+gsl_pow_2(img_y0-ex_y02)+gsl_pow_2(img_z0-ex_z02)))+" um\n";
        s+="expsf_r02 = "+floattostr(expsf_r02)+" microns\n";
        s+="expsf_z02 = "+floattostr(expsf_z02)+" microns\n";
    }

    s+="pos_detection = ["+floattostr(img_x0)+", "+floattostr(img_y0)+", "+floattostr(img_z0)+"] um\n";
    s+="detpsf_r0 = "+floattostr(detpsf_r0)+" microns\n";
    s+="detpsf_z0 = "+floattostr(detpsf_z0)+" microns\n";
    s+="pixel_size = "+floattostr(pixel_size)+" microns\n";
    s+="pixel_size_integrationdelta = "+floattostr(pixel_size_integrationdelta)+" microns\n";
    double psf_r0=1.0/sqrt(1.0/detpsf_r0/detpsf_r0+1.0/expsf_r0/expsf_r0);
    s+="confocal_psf_system = sqrt(1/expsf_r0^2 + 1/detpsf_r0^2) = "+floattostr(psf_r0)+" microns\n";
    double Veff=pow(M_PI, 1.5)*(psf_r0*psf_r0*psf_r0*detpsf_z0/detpsf_r0);
    s+="confocal_focus_volume (Veff) = pi^(3/2) * psf_system^3 * detpsf_z0/detpsf_r0 = "+floattostr(Veff)+" femto litre\n";

    double psf_r02=0;
    double Veff2=0;
    if (lambda_ex2>0) {
        psf_r02=1.0/sqrt(1.0/detpsf_r0/detpsf_r0+1.0/expsf_r02/expsf_r02);
        s+="confocal_psf_system2 = sqrt(1/expsf_r02^2 + 1/detpsf_r0^2) = "+floattostr(psf_r02)+" microns\n";
        Veff2=pow(M_PI, 1.5)*(psf_r02*psf_r02*psf_r02*detpsf_z0/detpsf_r0);
        s+="confocal_focus_volume2 (Veff) = pi^(3/2) * psf_system^3 * detpsf_z0/detpsf_r0 = "+floattostr(Veff)+" femto litre\n";
    }

    s+="illumination_distribution = "+ill_distribution_to_str(ill_distribution)+"\n";
    s+="detection_distribution = "+det_distribution_to_str(det_distribution)+"\n";
    s+="psf_rz_image_rwidth = "+inttostr(psf_rz_image_rwidth)+" pixels\n";
    s+="psf_rz_image_zwidth = "+inttostr(psf_rz_image_zwidth)+" pixels\n";
    s+="psf_rz_image_rresolution = "+floattostr(psf_rz_image_rresolution*1000.0)+" nm\n";
    s+="psf_rz_image_zresolution = "+floattostr(psf_rz_image_zresolution*1000.0)+" nm\n";


    double sum=0;
    if (dyn.size()>0) {
        for (size_t i=0; i<dyn.size(); i++) {
            s+="<"+dyn[i]->get_object_name()+" particles in Veff> = "+floattostr(Veff*dyn[i]->get_c_fluor()*6.022e23*1e-9/1e15)+"\n";
            sum=sum+Veff*dyn[i]->get_c_fluor()*6.022e23*1e-9/1e15;
        }
        s+="  SUM <particles in Veff> = "+floattostr(sum)+"\n";
    }
    if (lambda_ex2>0 && dyn.size()>0) {
        for (size_t i=0; i<dyn.size(); i++) {
            s+="<"+dyn[i]->get_object_name()+" particles in Veff2> = "+floattostr(Veff2*dyn[i]->get_c_fluor()*6.022e23*1e-9/1e15)+"\n";
            sum=sum+Veff2*dyn[i]->get_c_fluor()*6.022e23*1e-9/1e15;
        }
        s+="  SUM <particles in Veff2> = "+floattostr(sum)+"\n";
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
    s+="min_photons = "+inttostr(min_photons)+"\n";

    if (detector_type==1) {
        s+="detector_type = linear\n";
        s+="detector_resolution = "+inttostr(lindet_bits)+"bits  =>  range = 0..."+inttostr(pow(2,lindet_bits))+"\n";
        s+="detector_gain = "+floattostr(lindet_gain)+"\n";
        s+="lindet_readnoise = "+floattostr(lindet_readnoise)+"\n";
        s+="detector_variance_factor = "+floattostr(lindet_var_factor)+"\n";
    } else {
        s+="detector_type = photon_counting\n";
    }
    s+="detection_efficiency = "+floattostr(q_det*100.0)+"%\n";

    s+="lambda_ex = "+floattostr(lambda_ex)+" nm\n";
    if (lambda_ex2>0) s+="lambda_ex2 = "+floattostr(lambda_ex2)+" nm\n";
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
    s+="EPhoton_ex2 = "+floattostr(Ephoton2/1.602176487e-19/1e-3)+" meV\n";
    s+="I0 = "+floattostr(I0)+" uW/m^2  =  "+floattostr(I0/1e12)+" uW/micron^2  =  "+floattostr(I0/1e6)+" uW/mm^2\n";
    s+="P0 [on focus, i.e. on A=pi*(2*expsf_r0)^2] = "+floattostr(I0*(M_PI*gsl_pow_2(2.0*expsf_r0*1e-6)))+" uW\n";
    if (lambda_ex2>0) {
        s+="I02 = "+floattostr(I02)+" uW/m^2  =  "+floattostr(I02/1e12)+" uW/micron^2  =  "+floattostr(I02/1e6)+" uW/mm^2\n";
        s+="P02 [on focus, i.e. on A=pi*(2*expsf_r02)^2] = "+floattostr(I02*(M_PI*gsl_pow_2(2.0*expsf_r02*1e-6)))+" uW\n";
    }
    for (size_t i=0; i<dyn.size(); i++) {
        s+="  using "+dyn[i]->get_object_name()+" data:\n";
        s+="    init_spectrum = "+floattostr(dyn[i]->get_init_spectrum())+"\n";
        s+="    sigma_abs = "+format("%g", dyn[i]->get_init_sigma_abs(0),-1, false, 1e-30)+" m^2\n";
        s+="    q_fluor = "+floattostr(dyn[i]->get_init_q_fluor(0)*100)+" %\n";
        s+="    abs_photons/(molecule*s) @ lambda_ex  = "+floattostr(I0*1e-6/Ephoton*dyn[i]->get_init_sigma_abs(0))+"\n";
        s+="    fluor_photons/(molecule*s) @ lambda_ex = "+floattostr(dyn[i]->get_init_q_fluor(0)*I0*1e-6/Ephoton*dyn[i]->get_init_sigma_abs(0))+"\n";
        s+="    det_photons/(molecule*s) @ lambda_ex = "+floattostr(q_det*dyn[i]->get_init_q_fluor(0)*I0*1e-6/Ephoton*dyn[i]->get_init_sigma_abs(0))+"\n";
        s+="    absorbtion @ lambda_ex = "+floattostr(fluorophors->get_spectral_absorbance(dyn[i]->get_init_spectrum(), lambda_ex)*100.0)+" %\n";
        if (lambda_ex2>0) {
            s+="    abs_photons/(molecule*s) @ lambda_ex2  = "+floattostr(I02*1e-6/Ephoton2*dyn[i]->get_init_sigma_abs(0))+"\n";
            s+="    fluor_photons/(molecule*s) @ lambda_ex2 = "+floattostr(dyn[i]->get_init_q_fluor(0)*I02*1e-6/Ephoton2*dyn[i]->get_init_sigma_abs(0))+"\n";
            s+="    det_photons/(molecule*s) @ lambda_ex2 = "+floattostr(q_det*dyn[i]->get_init_q_fluor(0)*I02*1e-6/Ephoton2*dyn[i]->get_init_sigma_abs(0))+"\n";
            s+="    absorbtion @ lambda_ex2 = "+floattostr(fluorophors->get_spectral_absorbance(dyn[i]->get_init_spectrum(), lambda_ex2)*100.0)+" %\n";
        }
        if (det_wavelength_min>0 && det_wavelength_max>0) {
            s+="    deteted fraction of emission spectrum = "+floattostr(fluorophors->get_spectral_fluorescence(dyn[i]->get_init_spectrum(),det_wavelength_min, det_wavelength_max)*100.0)+" %\n";
        }
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
    s+="fccs_partner = "+fccs_partner+"\n";
    s+="save_arrivaltimes = "+booltostr(save_arrivaltimes)+"\n";
    if (save_arrivaltimes) {
        s+="arrivaltimes_onlyonce = "+booltostr(arrivaltimes_onlyonce)+"\n";
        s+="registered_arrival_times = "+inttostr(arrival_times.size())+"\n";
    }
    s+="detected_photons = "+inttostr(photoncounter)+"\n";
    s+="estimated_memory_consumtion_all = "+floattostr((mem+sum)/1024.0/1024.0)+" MBytes\n";
    s+="runtime correlation = "+floattostr(correlation_runtime)+" secs\n";

    return s;
}

 std::string FCSMeasurement::ill_distribution_to_str(int i) const {
    if (i==0) return "gaussian";
    if (i==1) return "gaussian_spim";
    if (i==2) return "slit_spim";
    return "unknown";
 }

 std::string FCSMeasurement::det_distribution_to_str(int i) const {
    if (i==0) return "gaussian";
    if (i==1) return "square_pixel";
    if (i==2) return "gaussian_beam";
    if (i==3) return "gaussian_beam_pixel";
    return "unknown";
 }


int FCSMeasurement::str_to_ill_distribution(std::string i) const {
    std::string s=tolower(i);
    if (s=="gaussian") return 0;
    if (s=="gaussian_spim") return 1;
    if (s=="slit_spim") return 2;
    return strtol(s.c_str(), NULL, 10);
 }

int FCSMeasurement::str_to_det_distribution(std::string i) const {
    std::string s=tolower(i);
    if (s=="gaussian") return 0;
    if (s=="square_pixel") return 1;
    if (s=="gaussian_beam") return 2;
    if (s=="gaussian_beam_pixel") return 3;
    return strtol(s.c_str(), NULL, 10);
 }

bool FCSMeasurement::getTimeSeries(int32_t** timeseries, int64_t& timeseries_size) {
    if (online_correlation) return false;
    *timeseries=this->timeseries;
    timeseries_size=this->timeseries_size;
    return timeseries_ended;
}

bool FCSMeasurement::getLastNPhotons(int64_t& N) {
    if (online_correlation) {
            N=lastNPhotons;
        return true;
    }
    return false;
}


bool FCSMeasurement::depends_on(const FluorescenceMeasurement* other) const {
    return get_fccs_partner_object()==other;
}

FCSMeasurement* FCSMeasurement::get_fccs_partner_object() const {
    if (measmap.count(fccs_partner)>0) {
        return dynamic_cast<FCSMeasurement*>(measmap[fccs_partner]);
    } else if (fccs_partner.size()>0) {
        throw MeasurementException("Did not find FCCS partner: "+fccs_partner);
    }
    return NULL;
}
