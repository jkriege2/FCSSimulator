#include "alvtools.h"
#include <math.h>
#include <iostream>

void alv5000WriteHeader(FILE* f, std::string object_name, double duration, const char* mode, double wavelength, double cr0, double cr1) {
    time_t rawtime;
    time (&rawtime);
    struct tm * timeinfo = localtime (&rawtime);;
    char timebuffer[256];

    fprintf(f, "ALV-5000/E-WIN Data\n");
    strftime (timebuffer,256,"%e-%b-%C",timeinfo);
    fprintf(f, "Date :	\"%s\"\n", timebuffer);
    strftime (timebuffer,256,"%T",timeinfo);
    fprintf(f, "Time :	\"%s\"\n", timebuffer);
    fprintf(f, "Samplename : 	\"%s\"\n", object_name.c_str());
    for (int i=0; i<=9; i++) fprintf(f, "SampMemo(%d) : 	\"\"\n",i);
    fprintf(f, "Temperature [K] :	     298.16000\n");
    fprintf(f, "Viscosity [cp]  :	       0.89000\n");
    fprintf(f, "Refractive Index:	       1.33200\n");
    fprintf(f, "Wavelength [nm] :	     %lf\n", wavelength);
    fprintf(f, "Angle [�]       :	      90.00000\n");
    fprintf(f, "Duration [s]    :	     %lf\n", duration);
    fprintf(f, "Runs            :	         1\n");
    fprintf(f, "Mode            :	\"%s\"\n", mode);
    fprintf(f, "MeanCR0 [kHz]   :	       %lf\n", cr0);
    fprintf(f, "MeanCR1 [kHz]   :	       %lf\n", cr1);
}


void alv5000WriteCorrelation(FILE* f, int istart, int slots, double* corr_tau, double* corr_1, double* corr_2) {
    if (!f || slots<=0 || !corr_tau || !corr_1) return;
    fprintf(f, "\n");
    fprintf(f, "\"Correlation\"\n");
    for (int i=istart; i<slots; i++) {
        if (corr_2) {
            if (isfinite(corr_tau[i])&&isfinite(corr_1[i])&&isfinite(corr_2[i])) fprintf(f, "%14.5lg\t%14.5lg\t%14.5lg\n", corr_tau[i]*1e3, corr_1[i], corr_2[i]);
        } else {
            if (isfinite(corr_tau[i])&&isfinite(corr_1[i])) fprintf(f, "%14.5lg\t%14.5lg\n", corr_tau[i]*1e3, corr_1[i]);
        }
    }
}

void alv5000WriteCountrate(FILE* f, int bts_N, double* bts_time, double* bts_1, double* bts_2, bool bts_is_photoncounts) {
    //std::cout<<"alv5000WriteCountrate("<<f<<", "<<bts_N<<", "<<bts_time<<", "<<bts_1<<", "<<bts_2<<", "<<bts_is_photoncounts<<")\n";
    if (f && bts_time && bts_1 && bts_N>1) {
        fprintf(f, "\n");
        fprintf(f, "\"Count Rate\"\n");
        double factor=1;
        if (bts_is_photoncounts) factor=1.0/(bts_time[1]-bts_time[0])/1000.0;
        for (int i=0; i<bts_N; i++) {
            double t=bts_time[i];
            double c1=bts_1[i]*factor;
            double c2=0;
            if (bts_2) c2=bts_2[i]*factor;
            if (!isfinite(c1)) c1=0;
            if (!isfinite(c2)) c2=0;
            if (bts_2) fprintf(f, "%14.5lf\t%14.5lf\t%14.5lf\n", t, c1, c2);
            else fprintf(f, "%14.5lf\t%14.5lf\n", t, c1);
        }
    }
}
