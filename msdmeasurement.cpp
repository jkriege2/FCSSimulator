#include "msdmeasurement.h"

MSDMeasurement::MSDMeasurement(FluorophorManager* fluorophors, std::string objectname):
    FluorescenceMeasurement(fluorophors, objectname)
{
    msd_for_trajectories=10;
    msd_s=10;
    msd_p=16;
    msd_m=2;
}

MSDMeasurement::~MSDMeasurement()
{
    for (int i=0; i<msd_for_trajectories; i++) {
        delete msds[i];
    }
    msds.clear();
    trajectoryinfo.clear();
}


void MSDMeasurement::init() {
    FluorescenceMeasurement::init();
    clear();
    for (int i=0; i<msd_for_trajectories; i++) {
        msds.push_back(new MultiTauMSD<double>(msd_s, msd_m, msd_p, sim_timestep));
        MSDMeasurement::trajectory_info ti;
        ti.sum2_x=ti.sum_x=0;
        ti.sum2_y=ti.sum_y=0;
        ti.sum2_z=ti.sum_z=0;
        ti.cnt=0;
        trajectoryinfo.push_back(ti);
    }
}

void MSDMeasurement::propagate() {
    tick();
    FluorescenceMeasurement::propagate();
    int t=0;
    for (size_t d=0; d<dyn.size(); d++) { // go through all dynamics objects that provide data for this measurement object
        FluorophorDynamics::walkerState* dynn=dyn[d]->get_walker_state();
        unsigned long wc=dyn[d]->get_walker_count();
        if (!dyn[d]->end_of_trajectory_reached()) for (unsigned long i=0; i<wc; i++) { // iterate through all walkers in the d-th dynamics object
            if (dynn[i].exists) {
                msds[t]->contribute_step(dynn[i].x, dynn[i].y, dynn[i].z);
                trajectoryinfo[t].sum_x=trajectoryinfo[t].sum_x+dynn[i].x;
                trajectoryinfo[t].sum2_x=trajectoryinfo[t].sum2_x+dynn[i].x*dynn[i].x;
                trajectoryinfo[t].sum_y=trajectoryinfo[t].sum_y+dynn[i].y;
                trajectoryinfo[t].sum2_y=trajectoryinfo[t].sum2_y+dynn[i].y*dynn[i].y;
                trajectoryinfo[t].sum_z=trajectoryinfo[t].sum_z+dynn[i].z;
                trajectoryinfo[t].sum2_z=trajectoryinfo[t].sum2_z+dynn[i].z*dynn[i].z;
                trajectoryinfo[t].cnt=trajectoryinfo[t].cnt+1;
            }
            t++;
            if (t>=msd_for_trajectories) break;
        }
        if (t>=msd_for_trajectories) break;
    }


    tock();
    runtime=runtime+get_duration();
}

std::string MSDMeasurement::report() {
    std::string s=FluorescenceMeasurement::report();
    s+="msd_for_trajectories = "+inttostr(msd_for_trajectories)+"\n";
    s+="msd_s                = "+inttostr(msd_s)+"\n";
    s+="msd_p                = "+inttostr(msd_p)+"\n";
    s+="msd_m                = "+inttostr(msd_m)+"\n";

    return s;
}

void MSDMeasurement::save() {
    FILE* f;
    char fn[255];

    if (msd_for_trajectories<=0) return;

    sprintf(fn, "%s%strajectorystatistics.dat", basename.c_str(), object_name.c_str());
    printf("writing '%s' ...", fn);
    f=fopen(fn, "w");
    fprintf(f, "# trajectory no, steps, mean_x, std_x, mean_y, std_y, mean_z, std_z\n\n");
    for (int t=0; t<msd_for_trajectories; t++) {
        double s=trajectoryinfo[t].sum_x;
        double s2=trajectoryinfo[t].sum2_x;
        double cnt=trajectoryinfo[t].cnt;
        fprintf(f, "%d, %15.0lf, %15.10lf, %15.10lf", t, cnt, s/double(cnt), sqrt(1.0/(double(cnt)-1.0)*(s2-s*s/double(cnt))));
        s=trajectoryinfo[t].sum_y;
        s2=trajectoryinfo[t].sum2_y;
        fprintf(f, ", %15.0lf, %15.10lf", s/double(cnt), sqrt(1.0/(double(cnt)-1.0)*(s2-s*s/double(cnt))));
        s=trajectoryinfo[t].sum_z;
        s2=trajectoryinfo[t].sum2_z;
        fprintf(f, ", %15.0lf, %15.10lf\n", s/double(cnt), sqrt(1.0/(double(cnt)-1.0)*(s2-s*s/double(cnt))));
    }
    fclose(f);
    printf(" done!\n");


    double* taus=msds[0]->getMSDTau();
    int slots=msds[0]->getSlots();
    for (int t=0; t<msd_for_trajectories; t++) {
        MultiTauMSD<double>* msd=msds[t];
        msd->normalize();
    }

    sprintf(fn, "%s%smsd.dat", basename.c_str(), object_name.c_str());
    printf("writing '%s' ...", fn);
    f=fopen(fn, "w");
    for (int i=0; i<slots; i++) {
        fprintf(f, "%15.10lf", taus[i]);
        double s=0, s2=0;;
        for (int t=0; t<msd_for_trajectories; t++) {
            double* m=msds[t]->getMSD();
            fprintf(f, ", %15.10lf", m[i]);
            s+=m[i];
            s2+=m[i]*m[i];
        }
        fprintf(f, ", %15.10lf, %15.10lf\n", s/double(msd_for_trajectories), sqrt(1.0/(double(msd_for_trajectories)-1.0)*(s2-s*s/double(msd_for_trajectories))));
    }
    fclose(f);
    printf(" done!\n");
    std::string corrfn=fn;

    sprintf(fn, "%s%smsdplot.plt", basename.c_str(), object_name.c_str());
    printf("writing '%s' ...", fn);
    f=fopen(fn, "w");
    fprintf(f, "f(tau, D, alpha)=6*D*(tau**alpha)\n");
    fprintf(f, "D=10\n");
    fprintf(f, "fit f(x, D, 1) \"%s\" using 1:%d:%d via D\n", extract_file_name(corrfn).c_str(), msd_for_trajectories+2, msd_for_trajectories+3);
    fprintf(f, "Da=10\n");
    fprintf(f, "alpha=1\n");
    fprintf(f, "fit f(x, Da, alpha) \"%s\" using 1:%d:%d via Da, alpha\n", extract_file_name(corrfn).c_str(), msd_for_trajectories+2, msd_for_trajectories+3);
    fprintf(f, "set xlabel \"lag time [seconds]\"\n");
    fprintf(f, "set ylabel \"msd [micron^2]\"\n");

    for (int plt=0; plt<2; plt++) {
        if (plt==0) {
            fprintf(f, "set terminal pdfcairo color solid font \"Arial, 7\" linewidth 2 size 20cm,15cm\n");
            fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"msdplot.pdf").c_str());
        } else if (plt==1) {
            fprintf(f, "set terminal wxt font \"Arial, 8\"\n");
            fprintf(f, "set output\n");
        }
        fprintf(f, "set logscale x\n");
        fprintf(f, "set logscale y\n");

        fprintf(f, "set title \"MSD: %s\"\n", object_name.c_str());
        fprintf(f, "plot  \"%s\" using 1:%d:%d title 'average&std. dev.' with yerrorbars lc rgbcolor \"blue\"\n", extract_file_name(corrfn).c_str(), msd_for_trajectories+2, msd_for_trajectories+3);

        if (plt==1) fprintf(f, "pause -1\n");

        fprintf(f, "set title \"MSD: %s\"\n", object_name.c_str());
        fprintf(f, "plot \\\n");
        for (int p=0; p<msd_for_trajectories; p++) {
            fprintf(f, "   \"%s\" using 1:%d notitle with lines lc rgbcolor \"gray30\" lw 1, \\\n", extract_file_name(corrfn).c_str(), p+2);
        }
        fprintf(f, "   \"%s\" using 1:%d:%d title 'average&std. dev.' with yerrorbars lc rgbcolor \"blue\"\n", extract_file_name(corrfn).c_str(), msd_for_trajectories+2, msd_for_trajectories+3);

        if (plt==1) fprintf(f, "pause -1\n");

        fprintf(f, "set title \"MSD: %s\"\n", object_name.c_str());
        fprintf(f, "plot \\\n");
        for (int p=0; p<msd_for_trajectories; p++) {
            fprintf(f, "   \"%s\" using 1:%d notitle with lines lc rgbcolor \"gray30\" lw 1, \\\n", extract_file_name(corrfn).c_str(), p+2);
        }
        fprintf(f, "   \"%s\" using 1:%d:%d title 'average&std. dev.' with yerrorbars lc rgbcolor \"blue\", \\\n", extract_file_name(corrfn).c_str(), msd_for_trajectories+2, msd_for_trajectories+3);
        fprintf(f, "   f(x, D, 1) title sprintf('fit D=%%f, alpha=1', D) with lines lc rgbcolor \"green\" lw 2, \\\n");
        fprintf(f, "   f(x, Da, alpha) title sprintf('fit D=%%f, alpha=%%f', Da, alpha) with lines lc rgbcolor \"red\" lw 2\n");

        if (plt==1) fprintf(f, "pause -1\n");
    }

    fclose(f);
    printf(" done!\n");
}

void MSDMeasurement::read_config_internal(jkINIParser2& parser) {
    FluorescenceMeasurement::read_config_internal(parser);

    msd_for_trajectories=parser.getAsInt("msd_for_trajectories", msd_for_trajectories);
    msd_s=parser.getAsInt("msd_s", msd_s);
    msd_p=parser.getAsInt("msd_p", msd_p);
    msd_m=parser.getAsInt("msd_m", msd_m);

    //init();
}

void MSDMeasurement::clear() {
    for (size_t i=0; i<msds.size(); i++) {
        delete msds[i];
    }
    msds.clear();
}
