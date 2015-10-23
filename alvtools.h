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


#ifndef ALVTOOLS_H
#define ALVTOOLS_H

#include <ctime>
#include <cstdio>
#include <cstdio>
#include <string>

void alv5000WriteHeader(FILE* f, std::string object_name, double duration_seconds, const char* mode, double wavelength_nanometers, double cr0=0, double cr1=0);

void alv5000WriteCorrelation(FILE* f, int istart, int slots, double* corr_tau_seconds, double* corr_1, double* corr_2=NULL, double subtractFromCF=0.0);
void alv5000WriteCountrate(FILE* f, int bts_N, double* bts_time_seconds, double* bts_1, double* bts_2, bool bts_is_photoncounts);

void qf3acorrWriteHeader(FILE* f, std::string object_name, double duration, double wavelength, int ncorrelations, int channels, bool isFCCS=false);

void qf3acorrWriteCorrelation(FILE* f, int istart, int slots, double* corr_tau_seconds, double* corr_1, double* corr_2=NULL, double* corr_12=NULL, double subtractFromCF=0.0);
void qf3acorrWriteCountrate(FILE* f, int bts_N, double* bts_time_seconds, double* bts_1, double* bts_2, bool bts_is_photoncounts);


#endif // ALVTOOLS_H
