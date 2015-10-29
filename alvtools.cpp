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


#include "alvtools.h"
#include <math.h>
#include <iostream>
#include <time.h>

void alv5000WriteHeader(FILE* f, std::string object_name, double duration, const char* mode, double wavelength, double cr0, double cr1) {
    time_t rawtime;
    time (&rawtime);
    //struct tm * timeinfo = localtime (&rawtime);;
    //char timebuffer[256];

    fprintf(f, "ALV-5000/E-WIN Data\n");

    /*strftime (timebuffer,256,"%e-%b-%C",timeinfo);

    fprintf(f, "Date :	\"%s\"\n", timebuffer);

    strftime (timebuffer,256,"%T",timeinfo);

    fprintf(f, "Time :	\"%s\"\n", timebuffer);*/
    fprintf(f, "Date :	\"\"\n");
    fprintf(f, "Time :	\"\"\n");

    fprintf(f, "Samplename : 	\"%s\"\n", object_name.c_str());

    for (int i=0; i<=9; i++) fprintf(f, "SampMemo(%d) : 	\"\"\n",i);

    fprintf(f, "Temperature [K] :	     298.16000\n");

    fprintf(f, "Viscosity [cp]  :	       0.89000\n");

    fprintf(f, "Refractive Index:	       1.33200\n");

    fprintf(f, "Wavelength [nm] :	     %lf\n", wavelength);

    fprintf(f, "Angle [°]       :	      90.00000\n");

    fprintf(f, "Duration [s]    :	     %lf\n", duration);

    fprintf(f, "Runs            :	         1\n");

    fprintf(f, "Mode            :	\"%s\"\n", mode);

    fprintf(f, "MeanCR0 [kHz]   :	       %lf\n", cr0);

    fprintf(f, "MeanCR1 [kHz]   :	       %lf\n", cr1);

}


void alv5000WriteCorrelation(FILE* f, int istart, int slots, double* corr_tau, double* corr_1, double* corr_2, double subtractFromCF) {
    if (!f || slots<=0 || !corr_tau || !corr_1) return;
    fprintf(f, "\n");
    fprintf(f, "\"Correlation\"\n");
    for (int i=istart; i<slots; i++) {
        if (corr_2) {
            if (isfinite(corr_tau[i])&&isfinite(corr_1[i])&&isfinite(corr_2[i])) fprintf(f, "%14.5lg\t%14.5lg\t%14.5lg\n", corr_tau[i]*1e3, corr_1[i]-subtractFromCF, corr_2[i]-subtractFromCF);
        } else {
            if ( isfinite(corr_tau[i]) && isfinite(corr_1[i]) ) fprintf(f, "%14.5lg\t%14.5lg\n", corr_tau[i]*1e3, corr_1[i]-subtractFromCF);
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
        int NUP=bts_N;
        for (int i=bts_N-1; i>=0; i--) {
            if ((!bts_1 || (bts_1[i]==0)) && (!bts_2 || (bts_2[i]==0))) {
                NUP--;
            } else {
                break;
            }
        }
        if (NUP>0) {
            for (int i=0; i<NUP; i++) {
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
}


void qf3acorrWriteHeader(FILE *f, std::string object_name, double duration, double wavelength, int correlations, int channels, bool isFCCS)
{
    fprintf(f, "[Properties]\ncodec = ISO-8859-1\nfileformat_name = QF3ASCIICorrelationData\nfileformat_version = \"1.0\"\n");
    fprintf(f, "Samplename=\"%s\"\n", object_name.c_str());
    fprintf(f, "Wavelength_nm=%lf\n", wavelength);
    fprintf(f, "duration_seconds=%lf\n", duration);
    fprintf(f, "duration_ms=%lf\n", duration*1e3);
    fprintf(f, "runs=1\n");
    fprintf(f, "rateRuns=1\n");
    fprintf(f, "correlations=%d\n", correlations);
    fprintf(f, "channels=%d\n", channels);

    if (isFCCS) {
        if (correlations==1) {
            fprintf(f, "role=FCCS\n");
        } else {
            fprintf(f, "role0=ACF0\n");
            fprintf(f, "role1=ACF1\n");
            fprintf(f, "role2=FCCS\n");
            fprintf(f, "preferred_channel0=0\n");
            fprintf(f, "preferred_channel1=1\n");
            fprintf(f, "preferred_channel2=0\n");
        }
    } else {
        for (int i=0; i<correlations; i++) {
            fprintf(f, "role%d=ACF%d\n", i,i);
        }
        for (int i=0; i<correlations; i++) {
            fprintf(f, "preferred_channel%d=%d\n", i,i);
        }
    }


}


void qf3acorrWriteCorrelation(FILE *f, int istart, int slots, double *corr_tau, double *corr_1, double *corr_2, double *corr_12, double subtractFromCF)
{
    if (!f || slots<=0 || !corr_tau || !corr_1) return;
    fprintf(f, "\n");
    fprintf(f, "[CorrelationData]\n");
    for (int i=istart; i<slots; i++) {
        if (corr_2 && corr_12) {
            if (isfinite(corr_tau[i])&&isfinite(corr_1[i])&&isfinite(corr_2[i])&&isfinite(corr_12[i])) fprintf(f, "%14.5lg, %14.5lg, %14.5lg, %14.5lg\n", corr_tau[i], corr_1[i]-subtractFromCF, corr_2[i]-subtractFromCF, corr_12[i]-subtractFromCF);
        } else if (corr_2) {
            if (isfinite(corr_tau[i])&&isfinite(corr_1[i])&&isfinite(corr_2[i])) fprintf(f, "%14.5lg, %14.5lg, %14.5lg\n", corr_tau[i], corr_1[i]-subtractFromCF, corr_2[i]-subtractFromCF);
        } else {
            if ( isfinite(corr_tau[i]) && isfinite(corr_1[i]) ) fprintf(f, "%14.5lg, %14.5lg\n", corr_tau[i], corr_1[i]-subtractFromCF);
        }
    }

}


void qf3acorrWriteCountrate(FILE *f, int bts_N, double *bts_time, double *bts_1, double *bts_2, bool bts_is_photoncounts)
{
    if (f && bts_time && bts_1 && bts_N>1) {
        fprintf(f, "\n");
        fprintf(f, "[RateData]\n");
        double factor=1;
        if (bts_is_photoncounts) factor=1.0/(bts_time[1]-bts_time[0])/1000.0;
        int NUP=bts_N;
        for (int i=bts_N-1; i>=0; i--) {
            if ((!bts_1 || (bts_1[i]==0)) && (!bts_2 || (bts_2[i]==0))) {
                NUP--;
            } else {
                break;
            }
        }
        if (NUP>0) {
            for (int i=0; i<NUP; i++) {
                double t=bts_time[i];
                double c1=bts_1[i]*factor;
                double c2=0;
                if (bts_2) c2=bts_2[i]*factor;
                if (!isfinite(c1)) c1=0;
                if (!isfinite(c2)) c2=0;
                if (bts_2) fprintf(f, "%14.5lf, %14.5lf, %14.5lf\n", t, c1, c2);
                else fprintf(f, "%14.5lf, %14.5lf\n", t, c1);
            }
        }
    }
}
