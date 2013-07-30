#ifndef ALVTOOLS_H
#define ALVTOOLS_H

#include <ctime>
#include <cstdio>
#include <cstdio>
#include <string>

void alv5000WriteHeader(FILE* f, std::string object_name, double duration, const char* mode, double wavelength, double cr0=0, double cr1=0);

void alv5000WriteCorrelation(FILE* f, int istart, int slots, double* corr_tau, double* corr_1, double* corr_2=NULL);
void alv5000WriteCountrate(FILE* f, int bts_N, double* bts_time, double* bts_1, double* bts_2, bool bts_is_photoncounts);
#endif // ALVTOOLS_H
