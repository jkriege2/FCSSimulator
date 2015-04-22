#include "diffusiontools.h"

std::vector<FluorophorDynamics*> dyn;
std::vector<FluorescenceMeasurement*> meas;
std::map<std::string, FluorophorDynamics*> dynmap;
std::map<std::string, FluorescenceMeasurement*> measmap;
gsl_rng * global_rng=NULL;

std::map<std::string, std::string> preset_ini_params;

FluorophorManager::FluorophorManager(std::string database_path, bool test_spectra){
    this->database_path=database_path;
    this->test_spectra=test_spectra;
    spectral_interpolation_type=gsl_interp_cspline;
    init_fluorophor_database();
    init_spectra();
}

FluorophorManager::~FluorophorManager() {
    for (size_t i=0; i<spectra.size(); i++) {
        if (spectra[i].loaded) {
            gsl_spline_free(spectra[i].spline_abs);
            gsl_interp_accel_free(spectra[i].accel_abs);
            free(spectra[i].lambda);
            free(spectra[i].eff_abs);
        }
    }
    spectra.clear();
    spectra_map.clear();
}

void FluorophorManager::init_fluorophor_database() {
    std::string dirname=include_trailing_backslash(database_path)+"spectra"+std::string(PATHSEPARATOR_STRING);
    jkINIParser2 ini(dirname+"fluorophors.ini");
    fluorophor_database.clear();
    for (int i=0; i<ini.getGroupCount(); i++) {
        std::string g=tolower(ini.getGroupName(i));
        std::string gg=ini.getGroupName(i);
        fluorophor_database[g].fluorescence_efficiency=ini.getAsDouble(gg+".q_fluor", 1);
        fluorophor_database[g].fluorescence_lifetime=ini.getAsDouble(gg+".tau_fluor", 1e-9);
        if (ini.exists(gg+".sigma_abs")) {
            fluorophor_database[g].sigma_abs=ini.getAsDouble(gg+".sigma_abs", 2e-20);
        } else {
            fluorophor_database[g].sigma_abs=ini.getAsDouble(gg+".molar_extinction", 2e4)/10.0/6.022e23;
        }
        fluorophor_database[g].bleaching_propability=ini.getAsDouble(gg+".p_bleach", 0);
        fluorophor_database[g].triplet_lifetime=ini.getAsDouble(gg+".tau_trip", 0);
        fluorophor_database[g].triplet_propability=ini.getAsDouble(gg+".p_trip", 0);
        fluorophor_database[g].spectrum=-1;
        std::cout<<"added fluorophor '"<<g<<"' to fluorophor database \n";
    }
}


void FluorophorManager::add_spectrum(std::string filename) {
    Spectrum s;
    s.filename_abs=tolower(filename);
    s.loaded=false;
    spectra.push_back(s);
    std::string specname=change_file_ext(extract_file_name(tolower(filename)), "");
    spectra_map[specname]=spectra.size()-1;
    if (fluorophor_database.find(specname)!=fluorophor_database.end()) {
        fluorophor_database[specname].spectrum=spectra.size()-1;
    }
}

void FluorophorManager::load_spectrum(int ID) {
    if (ID<0) return;
    //std::cout<<"waiting for load_spectrum ID="<<ID<<" ...\n";
    //boost::mutex::scoped_lock scoped_lock(spectrum_load_lock);
    //std::cout<<"in load_spectrum ...\n";
    if (spectra[ID].loaded) {
        //std::cout<<"finished load_spectrum ...\n";
        return;
    }
    //std::cout<<"run load_spectrum ...\n";
    datatable tab;
    std::string filename=spectra[ID].filename_abs;
    //std::cout<<"loading spectrum: '"<<filename<<"'\n";
    tab.load_csv(filename);
    spectra[ID].val_count=tab.get_line_count();
    spectra[ID].lambda=(double*)calloc(spectra[ID].val_count, sizeof(double));
    spectra[ID].eff_abs=(double*)calloc(spectra[ID].val_count, sizeof(double));
    spectra[ID].eff_fl=(double*)calloc(spectra[ID].val_count, sizeof(double));
    double max_abs=tab.column_max(1);
    double sum_fl=1;
    if (tab.get_column_count()>2) sum_fl=tab.column_sum(2);
    double max_fl=1;
    if (tab.get_column_count()>2) max_fl=tab.column_max(2);
    //std::cout<<"max_abs="<<max_abs<<"   sum_fl="<<sum_fl<<std::endl;
    for (size_t i=0; i<spectra[ID].val_count; i++) {
        spectra[ID].lambda[i]=tab.get(0, i);
        spectra[ID].eff_abs[i]=tab.get(1, i)/max_abs;
        if (tab.get_column_count()>2) spectra[ID].eff_fl[i]=tab.get(2, i)/sum_fl;
    }
    spectra[ID].accel_abs= gsl_interp_accel_alloc ();
    spectra[ID].spline_abs=gsl_spline_alloc (spectral_interpolation_type, spectra[ID].val_count);
    gsl_spline_init(spectra[ID].spline_abs, spectra[ID].lambda, spectra[ID].eff_abs, spectra[ID].val_count);
    spectra[ID].accel_fl= gsl_interp_accel_alloc ();
    spectra[ID].spline_fl=gsl_spline_alloc (spectral_interpolation_type, spectra[ID].val_count);
    gsl_spline_init(spectra[ID].spline_fl, spectra[ID].lambda, spectra[ID].eff_fl, spectra[ID].val_count);


    if (tab.get_column_count()>2) {
        double dl=spectra[ID].lambda[spectra[ID].val_count-1]-spectra[ID].lambda[0];
        double int_fl=gsl_spline_eval_integ(spectra[ID].spline_fl, spectra[ID].lambda[0], spectra[ID].lambda[0]+dl, spectra[ID].accel_fl);
        for (size_t i=0; i<spectra[ID].val_count; i++) {
            spectra[ID].eff_fl[i]=spectra[ID].eff_fl[i]/int_fl;
        }
        /*gsl_spline_free(spectra[ID.spline_abs);
        gsl_interp_accel_free(spectra[ID].accel_abs);
        spectra[ID].accel_fl= gsl_interp_accel_alloc ();
        spectra[ID].spline_fl=gsl_spline_alloc (spectral_interpolation_type, spectra[ID].val_count);*/
        gsl_spline_init(spectra[ID].spline_fl, spectra[ID].lambda, spectra[ID].eff_fl, spectra[ID].val_count);
    }


    spectra[ID].loaded=true;
    //std::cout<<"finished load_spectrum ...\n";
}


void FluorophorManager::test_spectrum(int ID) {
    if (ID<0) return;
    if (!spectra[ID].loaded) {
        load_spectrum(ID)    ;
    }
    //std::cout<<"run load_spectrum ...\n";
    datatable tab;
    std::string filename=spectra[ID].filename_abs;
    tab.load_csv(filename);
    double max_abs=tab.column_max(1);
    double sum_fl=1;
    if (tab.get_column_count()>2) sum_fl=tab.column_sum(2);
    double max_fl=1;
    if (tab.get_column_count()>2) max_fl=tab.column_max(2);

    double dl=spectra[ID].lambda[spectra[ID].val_count-1]-spectra[ID].lambda[0];
    double int_fl=0;
    gsl_spline* spline=spectra[ID].spline_fl;
    gsl_interp_accel* acc=spectra[ID].accel_fl;
    int_fl=gsl_spline_eval_integ(spline, spectra[ID].lambda[0], spectra[ID].lambda[0]+dl, acc);
    int N=500;
    std::string fn1=change_file_ext(filename, ".intspec");
    std::cout<<"writing spectrum interpolation: "<<fn1<<std::endl;
    FILE* o=fopen(fn1.c_str(), "w");
    for (int i=0; i<N; i++) {
        double l=spectra[ID].lambda[0]+dl/(double)(N-1)*i;
        double e=get_spectral_absorbance(ID, l);
        double f=get_spectral_fluorescence(ID, l);
        fprintf(o, "%10.5lf %10.5lf %10.5lf\n", l, e, f);
    }
    fclose(o);
    o=fopen((fn1+".plt").c_str(), "w");
    fprintf(o, "plot \"%s\" using 1:2 title \"absorption, interpolation\" with lines, \"%s\" using 1:(($2)/%lf) title \"absorption,input data\" with points", extract_file_name(fn1).c_str(), extract_file_name(filename).c_str(), max_abs);
    fprintf(o, "\n\n\n");
    fprintf(o, "pause -1\n");
    if (tab.get_column_count()>2) {
        fprintf(o, "plot \"%s\" using 1:3 title \"emission, interpolation, integral=%lf\" with lines, \"%s\" using 1:(($3)/%lf) title \"emission,input data, sum=%lf\" with points\n", extract_file_name(fn1).c_str(), int_fl, extract_file_name(filename).c_str(), sum_fl, sum_fl);
        fprintf(o, "\n\n\n");
        fprintf(o, "pause -1\n");
    }
    fclose(o);
}

void FluorophorManager::init_spectra() {
    DIR *dp;
    struct dirent *ep;


    //test_spectra=false;
    for (size_t i=0; i<spectra.size(); i++) {
        gsl_spline_free(spectra[i].spline_abs);
        gsl_interp_accel_free(spectra[i].accel_abs);
        free(spectra[i].lambda);
        free(spectra[i].eff_abs);
    }
    spectra.clear();
    spectra_map.clear();
    std::string dirname=include_trailing_backslash("."+std::string(PATHSEPARATOR_STRING)+"spectra");
    dp = opendir (dirname.c_str());
    if (dp != NULL) {
        while ((ep = readdir(dp))) {
            //std::cout<<ep->d_name<<"   "<<extract_file_ext(ep->d_name)<<std::endl;
            if (extract_file_ext(ep->d_name)=="spec") {
                std::string fn=dirname;
                fn=fn+ep->d_name;
                std::cout<<"loading "<<fn<<std::endl;
                add_spectrum(fn);
            }
        }
        (void) closedir (dp);
    }

    if (test_spectra) {
        std::cout<<"creating test plotfiles for spectra:\n";
        for (size_t i=0; i<spectra.size(); i++) {
            std::cout<<" * loading spectrum i="<<i<<": "<<spectra[i].filename_abs<<std::endl;
            test_spectrum(i);
        }
    }
    std::cout<<"loading ... finished\n";

}


double FluorophorManager::get_spectral_fluorescence(int spectrum, double wavelength_start, double wavelength_end) {
    if (spectrum==-1) return 1.0;
    if (!spectra[spectrum].loaded) load_spectrum(spectrum);
    gsl_spline* spline=spectra[spectrum].spline_fl;
    gsl_interp_accel* acc=spectra[spectrum].accel_fl;
    return gsl_spline_eval_integ(spline, wavelength_start, wavelength_end, acc);
}


double gaussbeam_w(double z, double z0, double w0) {
    return w0*sqrt(1.0+z*z/z0/z0);
}

double gaussbeam_R(double z, double z0) {
    return z*(1.0+(z0*z0/z/z));
}

double FluorophorManager::get_spectral_absorbance(int spectrum, double wavelength) {
    if (spectrum==-1) return 1.0;
    if (!spectra[spectrum].loaded) load_spectrum(spectrum);
    gsl_spline* spline=spectra[spectrum].spline_abs;
    gsl_interp_accel* acc=spectra[spectrum].accel_abs;
    if (wavelength<spectra[spectrum].lambda[0]) return 0;
    if (wavelength>=spectra[spectrum].lambda[spectra[spectrum].val_count-1]) return 0;
    return GSL_MIN(1.0, GSL_MAX(0.0, gsl_spline_eval(spline, wavelength, acc)));
}

double FluorophorManager::get_spectral_fluorescence(int spectrum, double wavelength) {
    if (spectrum==-1) return 1.0;
    if (!spectra[spectrum].loaded) load_spectrum(spectrum);
    gsl_spline* spline=spectra[spectrum].spline_fl;
    gsl_interp_accel* acc=spectra[spectrum].accel_fl;
    if (wavelength<spectra[spectrum].lambda[0]) return 0;
    if (wavelength>=spectra[spectrum].lambda[spectra[spectrum].val_count-1]) return 0;
    return GSL_MIN(1.0, GSL_MAX(0.0, gsl_spline_eval(spline, wavelength, acc)));
}
