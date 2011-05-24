#include "fluorescenceimaging.h"

#include <unistd.h>

FluorescenceImaging::FluorescenceImaging(FluorophorManager* fluorophors, std::string objectname):
    FluorescenceMeasurement(fluorophors, objectname)
{
    psf_r0=0.1;
    psf_z0=0.7;
    exposure_time=50e-6;
    imaging_delay=0.001e-3;
    I0=200/0.25e-12;

    // Rho 6G
    lambda_ex=527;
    sigma_abs=2.2e-20;
    q_fluor=0.98;

    q_det=0.1;

    img_z0=2;
    img_x0=2;
    img_y0=2;
    img_dx=10;
    img_dy=10;
    img_xpixel=30;
    img_ypixel=30;

    avgImage=NULL;
    cntImage=NULL;
    current_image=0;
    basename="";
    save_binary_images=false;
}

FluorescenceImaging::~FluorescenceImaging()
{
    if (avgImage!=NULL) gsl_matrix_free (avgImage);
    if (cntImage!=NULL) gsl_matrix_free (cntImage);
    avgImage=NULL;
    cntImage=NULL;
}

void FluorescenceImaging::read_config_internal(jkINIParser2& parser) {
    FluorescenceMeasurement::read_config_internal(parser);
    psf_r0=parser.getSetAsDouble("psf_r0", psf_r0);
    psf_z0=parser.getSetAsDouble("psf_z0", psf_z0);
    exposure_time=parser.getSetAsDouble("exposure_time", exposure_time);
    imaging_delay=parser.getSetAsDouble("imaging_delay", imaging_delay);
    I0=parser.getSetAsDouble("I0", I0);

    // Rho 6G
    lambda_ex=parser.getSetAsDouble("lambda_ex", lambda_ex);
    sigma_abs=parser.getSetAsDouble("sigma_abs", sigma_abs);
    q_fluor=parser.getSetAsDouble("q_fluor", q_fluor);

    q_det=parser.getSetAsDouble("q_det", q_det);

    img_z0=parser.getSetAsDouble("img_z0", img_z0);
    img_x0=parser.getSetAsDouble("img_x0", img_x0);
    img_y0=parser.getSetAsDouble("img_y0", img_y0);
    img_dx=parser.getSetAsDouble("img_dx", img_dx);
    img_dy=parser.getSetAsDouble("img_dy", img_dy);
    img_xpixel=parser.getSetAsDouble("img_xpixel", img_xpixel);
    img_ypixel=parser.getSetAsDouble("img_ypixel", img_ypixel);

    save_binary_images=parser.getSetAsBool("save_binary", save_binary_images);
    init(exposure_time, imaging_delay);
}


void FluorescenceImaging::init(double exposure_time, double imaging_delay){
    FluorescenceMeasurement::init();
    this->current_image=0;
    this->exposure_time=exposure_time;
    this->imaging_delay=imaging_delay;
    if (avgImage!=NULL) gsl_matrix_free (avgImage);
    if (cntImage!=NULL) gsl_matrix_free (cntImage);
    //images=(gsl_matrix**)calloc(n_images, sizeof(gsl_matrix*));
    avgImage=gsl_matrix_calloc(img_xpixel, img_ypixel);
    cntImage=gsl_matrix_calloc(img_xpixel, img_ypixel);

    // we start by acquiring an image
    acquiringImage=true;
    timeStepComplete=sim_time+this->exposure_time;
    //for (unsigned int i=0; i<this->n_images; i++) {
    Ephoton=6.626e-34*2.99e8/(lambda_ex*1e-9);
}

void FluorescenceImaging::propagate(){
    FluorescenceMeasurement::propagate();
    if (!acquiringImage) {
        if (sim_time>timeStepComplete) {
            // if delay is over: start to acquire the next image
            timeStepComplete=sim_time+exposure_time;
            acquiringImage=true;
        }
    }

    if (acquiringImage) {
        if (sim_time>timeStepComplete) {
            // change mode to delay calculate count image (poissonian distribution!) and save the lately acquired image
            timeStepComplete=sim_time+imaging_delay;
            acquiringImage=false;

            for (unsigned int ix=0; ix<img_xpixel; ix++) {
                for (unsigned int iy=0; iy<img_ypixel; iy++) {
                    gsl_matrix_set(cntImage, ix, iy, gsl_ran_poisson(rng, gsl_matrix_get(avgImage, ix, iy)));
                }
            }
            save_current_images();

        } else { // add intensity to the current images
            double n0=I0*1e-6/Ephoton*sim_timestep;//*q_fluor*sigma_abs;
            for (size_t d=0; d<dyn.size(); d++) { // go through all dynamics objects that provide data for this measurement object
                FluorophorDynamics::walkerState* dy=dyn[d]->get_walker_state();
                unsigned long wc=dyn[d]->get_walker_count();
                if (!dyn[d]->end_of_trajectory_reached()) for (unsigned long i=0; i<wc; i++) { // iterate through all walkers in the d-th dynamics object
                    if (dy[i].exists) {
                        double x0,y0,z0;
                        x0=dy[i].x;
                        y0=dy[i].y;
                        z0=dy[i].z;
                        //dy->get_walker_position(i, &x0, &y0, &z0);
                        if (fabs(z0-img_z0)<5.0*psf_z0)  {
                            unsigned int ix0, iy0, ix1, iy1;
                            double pdx=(img_dx/(double)img_xpixel);
                            double pdy=(img_dy/(double)img_ypixel);
                            ix0=GSL_MIN(GSL_MAX(0, (int)floor(((x0-img_x0)-5.0*psf_r0)/pdx)), img_xpixel);
                            ix1=GSL_MIN(GSL_MAX(0, (int)floor(((x0-img_x0)+5.0*psf_r0)/pdx)), img_xpixel);
                            iy0=GSL_MIN(GSL_MAX(0, (int)floor(((y0-img_y0)-5.0*psf_r0)/pdy)), img_ypixel);
                            iy1=GSL_MIN(GSL_MAX(0, (int)floor(((y0-img_y0)+5.0*psf_r0)/pdy)), img_ypixel);
                            //std::cout<<(x0-img_x0)<<"  "<<(y0-img_y0)<<"     "<<ix0<<".."<<ix1<<"   "<<iy0<<".."<<iy1<<std::endl;
                            for (unsigned int ix=ix0; ix<ix1; ix++) {
                                double x=img_x0+(double)ix*pdx;
                                for (unsigned int iy=iy0; iy<iy1; iy++) {
                                    double y=img_y0+(double)iy*pdy;
                                    double nphot=n0;
                                    double dxs=gsl_pow_2(x0-x);
                                    double dys=gsl_pow_2(y0-y);
                                    nphot=nphot*dyn[d]->get_walker_sigma_times_qfl(i);
                                    nphot=nphot*exp(-2.0*(dxs+dys)/gsl_pow_2(psf_r0)-2.0*gsl_pow_2(z0-img_z0)/gsl_pow_2(psf_z0));
                                    if (nphot>0) {
                                        gsl_matrix_set(avgImage, ix, iy, gsl_matrix_get(avgImage, ix, iy)+nphot);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void FluorescenceImaging::save_current_images(){
    //if (i<n_images-1) std::cout<<"now generating image "<<i+2<<"/"<<n_images<<" ... ";
    for (int j=0; j<2; j++) {
        char fn[1024];
        if (j==0) sprintf(fn, "%s%simg_avg%05d.dat", basename.c_str(), object_name.c_str(), current_image);
        else sprintf(fn, "%s%simg_cnt%05d.dat", basename.c_str(), object_name.c_str(), current_image);
        std::cout<<"writing '"<<fn<<"' ... ";
        if (!save_binary_images) {
            FILE* f=fopen(fn, "w");
            for (unsigned int iy=0; iy<img_ypixel; iy++) {
                for (unsigned int ix=0; ix<img_xpixel; ix++) {
                    if (ix>0) fprintf(f, ", ");
                    if (j==0) fprintf(f, "%10.10lf", gsl_matrix_get(avgImage, ix, iy));
                    else fprintf(f, "%10.0lf", gsl_matrix_get(cntImage, ix, iy));
                }
                fprintf(f, "\n");
            }
            fclose(f);
        } else {
            FILE* f=fopen(fn, "wb");
            for (unsigned int iy=0; iy<img_ypixel; iy++) {
                for (unsigned int ix=0; ix<img_xpixel; ix++) {
                    double d=0;
                    if (j==0) {
                        d=gsl_matrix_get(avgImage, ix, iy);
                        fwrite(&d, sizeof(d), 1, f);
                    } else {
                        uint32_t in=gsl_matrix_get(cntImage, ix, iy);
                        fwrite(&in, sizeof(in), 1, f);
                    }
                }
                fprintf(f, "\n");
            }
            fclose(f);
        }
        printf(" done!\n");
    }

    gsl_matrix_set_zero(avgImage);
    gsl_matrix_set_zero(cntImage);

    current_image++;
}

void FluorescenceImaging::save() {
}

std::string FluorescenceImaging::report() {
    std::string s=FluorescenceMeasurement::report();
    s+="image_pos0 = ["+floattostr(img_x0)+", "+floattostr(img_y0)+", "+floattostr(img_z0)+"] um\n";
    s+="image_size = "+floattostr(img_dx)+" * "+floattostr(img_dy)+" um^2;      "+inttostr(img_xpixel)+" * "+inttostr(img_ypixel)+" pixels\n";
    s+="pixel_size = "+floattostr(img_dx/(double)img_xpixel*1e3)+" * "+floattostr(img_dy/(double)img_ypixel*1e3)+" nm^2\n";
    double Veff=4.0/3.0*M_PI*psf_r0*psf_r0*psf_z0;
    s+="focus_volume (Veff) = "+floattostr(Veff)+" femto litre\n";
    s+="lambda_ex = "+floattostr(lambda_ex)+" nm\n";
    s+="EPhoton_ex = "+floattostr(Ephoton/1.602176487e-19/1e-3)+" meV\n";
    s+="I0 = "+floattostr(I0)+" uW/m^2  =  "+floattostr(I0/1e12)+" uW/micron^2\n";
//    s+="sigma_abs = "+floattostr(sigma_abs)+" m^2  =  "+floattostr(sigma_abs/1e-20)+" Angstrom^2\n";
//    s+="q_fluor = "+floattostr(q_fluor*100.0)+" %\n";
    s+="q_det = "+floattostr(q_det*100.0)+" %\n";
    s+="abs_photons/(molecule*s) = "+floattostr(I0*1e-6/Ephoton*sigma_abs)+"\n";
    s+="fluor_photons/(molecule*s) = "+floattostr(q_fluor*I0*1e-6/Ephoton*sigma_abs)+"\n";
    s+="det_photons/(molecule*s) = "+floattostr(q_det*q_fluor*I0*1e-6/Ephoton*sigma_abs)+"\n";
    s+="exposure_time = "+floattostr(exposure_time*1e3)+" msecs\n";
    s+="imaging_delay = "+floattostr(imaging_delay*1e3)+" msecs\n";
    s+="save_binary_images = "+booltostr(save_binary_images)+"\n";
    return s;
}
