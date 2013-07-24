#include <iostream>
#include "browniandynamics.h"
#include "fluorescenceimaging.h"
#include "jkiniparser2.h"
#include "fcsmeasurement.h"
#include "fluorophordynamics.h"
#include "datatable.h"
#include "teebuf.h"
#include "dynamicsfromfiles2.h"
#include "gridrandomwalkdynamics.h"
#include "msdmeasurement.h"
#include "childdynamics.h"
#include "trajectoryplot.h"
#include "fretdynamics.h"
#include <gsl/gsl_matrix.h>


#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

FluorophorManager* fluorophors;

struct DynamicsSortFunctor {
  bool operator() (const FluorophorDynamics* i, const FluorophorDynamics* j) {
      return j->depends_on(i);
  }
} myDynamicsSortFunctor;

struct MeasurementSortFunctor {
  bool operator() (const FluorescenceMeasurement* i, const FluorescenceMeasurement* j) {
      return j->depends_on(i);
  }
} myMeasurementSortFunctor;


void do_sim(std::string inifilename) {

    std::ofstream logfilestream;
    ScopedDelayedStreamDoubler2 stdredirect;

    try {
        jkINIParser2 ini;
        ini.readFile(inifilename);
        std::string basename=ini.getSetAsString("simulation.basename", "");
        logfilestream.open((basename+"log.txt").c_str());
        stdredirect.redirect(std::cout, logfilestream);
    } catch (std::exception& E) {
        std::cout<<"Error: COULD NOT CREATE LOG-FILE \n   "<<E.what()<<std::endl;
    }

    clock_t start, endtime;
    double cpu_time_used;

    start = clock();

    try {
        dyn.clear();
        meas.clear();
        dynmap.clear();
        measmap.clear();
        cout<<"reading config '"<<inifilename<<"' ... ";
        jkINIParser2 ini;
        ini.readFile(inifilename);
        cout<<"done!\n\n";

        std::string basename=ini.getSetAsString("simulation.basename", "");
        //bool multithread=ini.getSetAsBool("simulation.multithread", false);
        double duration=ini.getSetAsDouble("simulation.duration", -1);

        cout<<"creating directory for output '"<<extract_file_path(basename+"config.ini").c_str()<<"' ... ";
        mk_all_dir(extract_file_path(basename+"config.ini").c_str());
        cout<<"done!\n\n";


        ini.print();


        cout<<"writing '"<<basename<<"config.ini' ... ";
        ini.writeFile(basename+"config.ini");
        cout<<"done!\n\n";

        // first we read all groups in the ini file into a list, as the inifile may be changed by the
        // following instructions ... so also groups may be added and then the iteration is skrewed.
        std::vector<std::string> groups_in_ini;
        for (unsigned long i=0; i<ini.getGroupCount(); i++) {
            groups_in_ini.push_back(ini.getGroupName(i));
        }

        double estimated_max_runtime=0;

        // create all dynamics objects
        for (unsigned long i=0; i<groups_in_ini.size(); i++) {
            std::string gname=groups_in_ini[i];
            std::string lgname=tolower(gname);
            std::string oname=ini.getAsString(gname+".object_name", lgname);
            std::string supergroup="";
            FluorophorDynamics* d=NULL;
            //std::cout<<"** i="<<i<<"/"<<groups_in_ini.size()<<" gname="<<gname<<"  lgname="<<lgname<<"  oname="<<oname<<std::endl;
            if (lgname.find("brownian")==0 && lgname.size()>8) {
                supergroup="brownian";
                d=new BrownianDynamics(fluorophors, oname);
                if (duration<=0) throw std::runtime_error("if using brownian dynamics: duration has to be >0!");
            } else if (lgname.find("grid")==0 && lgname.size()>4) {
                supergroup="grid";
                d=new GridRandomWalkDynamics(fluorophors, oname);
                if (duration<=0) throw std::runtime_error("if using grid dynamics: duration has to be >0!");
            } else if (lgname.find("dynfile")==0 && lgname.size()>7) {
                supergroup="dynfile";
                d=new DynamicsFromFiles2(fluorophors, oname);
            } else if (lgname.find("child")==0 && lgname.size()>5) {
                supergroup="child";
                d=new ChildDynamics(fluorophors, oname);
            } else if (lgname.find("fret")==0 && lgname.size()>4) {
                supergroup="fret";
                d=new FRETDynamics(fluorophors, oname);
            }
            if (d!=NULL) {
                dyn.push_back(d);
                dynmap[oname]=d;
                dynmap[lgname]=d;
                d->set_basename(basename);
                if (ini.groupExists(gname+".test_dynamics")) {
                    d->read_config(ini, gname+".test_dynamics", supergroup);
                    d->test(ini.getSetAsInt(gname+".test_dynamics.sim_steps", 10000), ini.getSetAsInt(gname+".test_dynamics.walkers", 200));
                }
                d->read_config(ini, gname, supergroup);
                std::cout<<"created "<<supergroup<<" dynamics object '"<<oname<<"' ("<<lgname<<")\n";
            }
        }

        std::cout<<"sorting dynamics objects ...";
        std::sort(dyn.begin(), dyn.end(), myDynamicsSortFunctor);
        std::cout<<" done!\n";

        std::cout<<"initializing dynamics objects:\n";
        for (size_t i=0; i<dyn.size(); i++) {
            dyn[i]->init();
            double rt=dyn[i]->estimate_runtime();
            if (rt>estimated_max_runtime) estimated_max_runtime=rt;
            std::cout<<"    initialized "<< dyn[i]->get_supergroup() <<" dynamics object '"<<dyn[i]->get_object_name()<<"' ("<<dyn[i]->get_group()<<")     est. runtime="<<rt<<"s,   est. max. runtime="<<estimated_max_runtime<<"s\n";
        }
        std::cout<<"    DONE!\n";

        if (duration<=0 && estimated_max_runtime>0)  {
            duration=estimated_max_runtime;
            ini.setProperty("simulation.duration", duration);
            std::cout<<"\nno fixed simulation runtime given!!!\n";
            std::cout<<"  simulation runtime was estimated to be "<<duration<<" seconds\n\n";
        }
        // create all measurement objects and connect them to the dynamics
        for (unsigned long i=0; i<groups_in_ini.size(); i++) {
            std::string gname=groups_in_ini[i];
            std::string lgname=tolower(gname);
            std::string oname=ini.getAsString(gname+".object_name", lgname);
            std::vector<std::string> sources=tokenize_string(ini.getAsString(gname+".sources", ""), ",");
            if (sources.size()<=0) sources=tokenize_string(ini.getAsString(gname+".source", ""), ",");
            std::string supergroup="";
            int object_number=extract_right_int(gname);
            FluorescenceMeasurement* m=NULL;
            //std::cout<<"## i="<<i<<"/"<<groups_in_ini.size()<<"  gname="<<gname<<"  lgname="<<lgname<<"  oname="<<oname<<"  sources.size="<<sources.size()<<std::endl;
            if (sources.size()>0) {
                if (lgname.find("fcs")==0 && lgname.size()>3) {
                    supergroup="fcs";
                    m=new FCSMeasurement(fluorophors, oname);
                } else if (lgname.find("imaging")==0 && lgname.size()>7) {
                    supergroup="imaging";
                    m=new FluorescenceImaging(fluorophors, oname);
                } else if (lgname.find("msd")==0 && lgname.size()>3) {
                    supergroup="msd";
                    m=new MSDMeasurement(fluorophors, oname);
                } else if (lgname.find("trajectoryplot")==0 && lgname.size()>14) {
                    supergroup="trajectoryplot";
                    m=new TrajectoryPlotMeasurement(fluorophors, oname);
                }
                if (m!=NULL) {
                    m->set_object_number(object_number);
                    meas.push_back(m);
                    measmap[oname]=m;
                    measmap[lgname]=m;
                    m->read_config(ini, gname, supergroup);
                    m->init();
                    std::cout<<"created & initialized "<<supergroup<<" measurement object '"<<oname<<"' connected to: ";
                    m->clearDynamics();
                    for (size_t s=0; s<sources.size(); s++) {
                        std::cout<<"'"<<sources[s]<<"'/";
                        m->addDynamics(dynmap[sources[s]]);
                        if (s>0) std::cout<<", ";
                        std::cout<<dynmap[tolower(strstrip(sources[s]))]->get_object_name();
                    }
                    m->init();
                    std::cout<<"created "<<supergroup<<" measurement object '"<<oname<<"' connected to: ";
                    std::cout<<std::endl;
                }
            }
        }

        std::cout<<"sorting measurement objects ...";
        std::sort(meas.begin(), meas.end(), myMeasurementSortFunctor);
        std::cout<<" done!\n";



        if (dyn.size()<=0)  throw std::runtime_error("you need at least one dynamics object!");
        if (meas.size()<=0)  throw std::runtime_error("you need at least one measurement object!");

        ofstream filestr;
        filestr.open ((basename+"config.txt").c_str(), fstream::out);
        for (size_t i=0; i<dyn.size(); i++) {
            dyn[i]->save_results();
            filestr<<"\n\n\n---------------------------------------------------------------------------------------\n";
            filestr<<"- dynamics: "<< dyn[i]->get_object_name() <<"\n";
            filestr<<"---------------------------------------------------------------------------------------\n";
            filestr<<dyn[i]->report()<<endl;
        }

        for (size_t i=0; i<meas.size(); i++) {
            meas[i]->save_results();
            filestr<<"\n\n\n---------------------------------------------------------------------------------------\n";
            filestr<<"- measurement: "<< meas[i]->get_object_name() <<"\n";
            filestr<<"---------------------------------------------------------------------------------------\n";
            filestr<<meas[i]->report()<<endl;
        }
        filestr.close();


        if (duration>0) {
            // run simulation while any of the sim_time properties is <duration
            while (dyn[0]->get_sim_time()<duration) {
                for (size_t i=0; i<dyn.size(); i++) {
                    dyn[i]->propagate();
                    dyn[i]->propagate_additional_walkers();
                }
                for (size_t i=0; i<meas.size(); i++) {
                    meas[i]->propagate();
                }
            }
            for (size_t i=0; i<dyn.size(); i++) {
                dyn[i]->finalize_sim();
            }
            for (size_t i=0; i<meas.size(); i++) {
                meas[i]->finalize_sim();
            }
        } else {
            // run the simulation until all trajectories have ended
            bool done=false;
            while (!done) {
                done=true;
                for (size_t i=0; i<dyn.size(); i++) {
                    dyn[i]->propagate();
                    dyn[i]->propagate_additional_walkers();
                    done=done && dyn[i]->end_of_trajectory_reached();
                }
                for (size_t i=0; i<meas.size(); i++) {
                    meas[i]->propagate();
                }
            }
            for (size_t i=0; i<dyn.size(); i++) {
                dyn[i]->finalize_sim();
            }
            for (size_t i=0; i<meas.size(); i++) {
                meas[i]->finalize_sim();
            }
        }
        filestr.open ((basename+"config.txt").c_str(), fstream::out);
        for (size_t i=0; i<dyn.size(); i++) {
            dyn[i]->save_results();
            filestr<<"\n\n\n---------------------------------------------------------------------------------------\n";
            filestr<<"- dynamics: "<< dyn[i]->get_object_name() <<"\n";
            filestr<<"---------------------------------------------------------------------------------------\n";
            filestr<<dyn[i]->report()<<endl;
        }

        for (size_t i=0; i<meas.size(); i++) {
            meas[i]->save_results();
            filestr<<"\n\n\n---------------------------------------------------------------------------------------\n";
            filestr<<"- measurement: "<< meas[i]->get_object_name() <<"\n";
            filestr<<"---------------------------------------------------------------------------------------\n";
            filestr<<meas[i]->report()<<endl;
        }
        filestr.close();

    } catch (datatable_exception& E) {
        std::cout<<"Error:\n   "<<E.get_message()<<std::endl;
    } catch (std::exception& E) {
        std::cout<<"Error:\n   "<<E.what()<<std::endl;
    }
    for (size_t i=0; i<dyn.size(); i++) {
        delete dyn[i];
    }
    for (size_t i=0; i<meas.size(); i++) {
        delete meas[i];
    }
    dyn.clear();
    meas.clear();
    dynmap.clear();
    measmap.clear();
    endtime = clock();
    cpu_time_used = ((double) (endtime - start)) / CLOCKS_PER_SEC;
    std::cout<<"\n\n=======================================================\n= simulation duration: "<<cpu_time_used<<" secs\n=======================================================\n";
}

int main(int argc, char* argv[])
{
    bool test_spectra=false;
    bool do_execute=true;
    if (argc>1) {
        for (int i=1; i<argc; i++) {
            std::string fn=argv[i];
            if (tolower(fn)=="--testspectra")  {
                std::cout<<"spectrum test mode activated!\n";
                std::cout<<"\n";
                test_spectra=true;
                do_execute=false;
            }
            if (tolower(fn)=="--help")  {
                std::cout<<"diffusion 4 FCS simulator\n";
                std::cout<<"   (c)2008-2013 bz J.W.Krieger <j.krieger@dkfz.de>\n";
                std::cout<<"\nusage:\n";
                std::cout<<"    diffusion4 [options] file1 [file2 [file3 ...] ] ]\n";
                std::cout<<"\noptions:\n";
                std::cout<<"    --help: this online help message\n";
                std::cout<<"    --testspectra: output test infor for fl. and abs. spectra in database\n";
                std::cout<<"\nfiles:\n";
                std::cout<<"    give a list of files to process. It is possible to use wildcards in\n";
                std::cout<<"    filenames: * matches at least one character and ? matches exactlz one\n";
                std::cout<<"    character.\n\n\n";
                do_execute=false;
            }
        }
    }
    //std::cout<<"test spectra: "<<booltostr(test_spectra)<<std::endl;
    fluorophors=new FluorophorManager(extract_file_path(argv[0]), test_spectra);

    if (do_execute) {
        if (argc>1) {
            std::vector<std::string> files;
            for (int i=1; i<argc; i++) {
                std::string fn=argv[i];
                if (fn.size()<1 || fn[0]!='-') {
                    std::vector<std::string> files1=listfiles_wildcard(argv[i]);
                    std::string fn=argv[i];
                    if (fn.find('*')!=std::string::npos || fn.find('?')!=std::string::npos) {
                        files1=listfiles_wildcard(argv[i]);
                    } else {
                        files1.clear();
                        files1.push_back(fn);
                    }

                     for (size_t j=0; j<files1.size(); j++) {
                        files.push_back(files1[j]);
                        std::cout<<"   will simluate '"<<files1[j]<<"' ...\n";
                    }
                }
            }
            for (unsigned int i=0; i<files.size(); i++) {
                std::cout<<"---------------------------------------------------------------------------------------------------------"<<std::endl<<std::endl;
                std::cout<<"--  simulating for "<<files[i]<<"   "<<i+1<<"/"<<files.size()<<std::endl;
                std::cout<<"---------------------------------------------------------------------------------------------------------"<<std::endl<<std::endl;
                do_sim(files[i]);
            }
        } else {
            do_sim("diffusion4.ini");
        }
    }

    //dyn.save_trajectories();
    return 0;
}
