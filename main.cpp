#include <iostream>
#include "browniandynamics.h"
#include "fluorescenceimaging.h"
#include "../../../LIB/trunk/jkINIParser2.h"
#include "fcsmeasurement.h"
#include "fluorophordynamics.h"
#include "../../../LIB/trunk/datatable.h"
#include "dynamicsfromfiles2.h"

#include <gsl/gsl_matrix.h>

#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <map>

using namespace std;

FluorophorManager* fluorophors;

std::vector<FluorophorDynamics*> dyn;
std::vector<FluorescenceMeasurement*> meas;
std::map<std::string, FluorophorDynamics*> dynmap;
std::map<std::string, FluorescenceMeasurement*> measmap;

void do_sim(std::string inifilename) {
    try {
        dyn.clear();
        meas.clear();
        dynmap.clear();
        measmap.clear();
        jkINIParser2 ini;
        ini.readFile(inifilename);
        ini.print();
        std::string basename=ini.getSetAsString("simulation.basename", "");
        //bool multithread=ini.getSetAsBool("simulation.multithread", false);
        double duration=ini.getSetAsDouble("simulation.duration", -1);

        mk_all_dir(extract_file_path(basename+"config.ini").c_str());
        cout<<"writing '"<<basename<<"config.ini' ... ";
        ini.writeFile(basename+"config.ini");
        cout<<"done!\n\n";

        // create all dynamics objects
        for (unsigned long i=0; i<ini.getGroupCount(); i++) {
            std::string gname=ini.getGroupName(i);
            std::string lgname=tolower(ini.getGroupName(i));
            std::string oname=ini.getAsString(gname+".object_name", lgname);
            std::string supergroup="";
            FluorophorDynamics* d=NULL;
            if (lgname.find("brownian")==0 && lgname.size()>8) {
                supergroup="brownian";
                d=new BrownianDynamics(fluorophors, oname);
                if (duration<=0) throw std::runtime_error("if using brownian dynamics: duration has to be >0!");
            } else if (lgname.find("dynfile")==0 && lgname.size()>7) {
                supergroup="dynfile";
                d=new DynamicsFromFiles2(fluorophors, oname);
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

        // create all measurement objects and connect them to the dynamics
        for (unsigned long i=0; i<ini.getGroupCount(); i++) {
            std::string gname=ini.getGroupName(i);
            std::string lgname=tolower(ini.getGroupName(i));
            std::string oname=ini.getAsString(gname+".object_name", lgname);
            std::vector<std::string> sources=tokenize_string(ini.getAsString(gname+".sources", ""), ",");
            std::string supergroup="";
            FluorescenceMeasurement* m=NULL;
            if (sources.size()>0) {
                if (lgname.find("fcs")==0 && lgname.size()>3) {
                    supergroup="fcs";
                    m=new FCSMeasurement(fluorophors, oname);
                } else if (lgname.find("imaging")==0 && lgname.size()>7) {
                    supergroup="imaging";
                    m=new FluorescenceImaging(fluorophors, oname);
                }
                if (m!=NULL) {
                    meas.push_back(m);
                    measmap[oname]=m;
                    measmap[lgname]=m;
                    m->read_config(ini, gname, supergroup);
                    std::cout<<"created "<<supergroup<<" measurement object '"<<oname<<"' connected to: ";
                    m->clearDynamics();
                    for (size_t s=0; s<sources.size(); s++) {
                        std::cout<<"'"<<sources[s]<<"'/";
                        m->addDynamics(dynmap[sources[s]]);
                        if (s>0) std::cout<<", ";
                        std::cout<<dynmap[tolower(strstrip(sources[s]))]->get_object_name();
                    }
                    std::cout<<std::endl;
                }
            }
        }

        if (dyn.size()<=0)  throw std::runtime_error("you need at least one dynamics object!");
        if (meas.size()<=0)  throw std::runtime_error("you need at least one measurement object!");

        ofstream filestr;
        filestr.open ((basename+"config.txt").c_str(), fstream::out);
        for (size_t i=0; i<dyn.size(); i++) {
            filestr<<"\n\n\n---------------------------------------------------------------------------------------\n";
            filestr<<"- dynamics: "<< dyn[i]->get_object_name() <<"\n";
            filestr<<"---------------------------------------------------------------------------------------\n";
            filestr<<dyn[i]->report()<<endl;
        }

        for (size_t i=0; i<meas.size(); i++) {
            meas[i]->save();
            filestr<<"\n\n\n---------------------------------------------------------------------------------------\n";
            filestr<<"- measurement: "<< meas[i]->get_object_name() <<"\n";
            filestr<<"---------------------------------------------------------------------------------------\n";
            filestr<<meas[i]->report()<<endl;
        }
        filestr.close();


        if (duration>0) {
            // run simulation while any of the sim_time properties is <duration
            // we do not have to call init(), as this is done by the read_config() methods
            while (dyn[0]->get_sim_time()<duration) {
                for (size_t i=0; i<dyn.size(); i++) {
                    dyn[i]->propagate();
                }
                for (size_t i=0; i<meas.size(); i++) {
                    meas[i]->propagate();
                }
            }
        } else {
            // run the simulation until all trajectories have ended
            // we do not have to call init(), as this is done by the read_config() methods
            bool done=false;
            while (!done) {
                done=true;
                for (size_t i=0; i<dyn.size(); i++) {
                    dyn[i]->propagate();
                    done=done && dyn[i]->end_of_trajectory_reached();
                }
                for (size_t i=0; i<meas.size(); i++) {
                    meas[i]->propagate();
                }
            }
        }
        filestr.open ((basename+"config.txt").c_str(), fstream::out);
        for (size_t i=0; i<dyn.size(); i++) {
            filestr<<"\n\n\n---------------------------------------------------------------------------------------\n";
            filestr<<"- dynamics: "<< dyn[i]->get_object_name() <<"\n";
            filestr<<"---------------------------------------------------------------------------------------\n";
            filestr<<dyn[i]->report()<<endl;
        }

        for (size_t i=0; i<meas.size(); i++) {
            meas[i]->save();
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

}

int main(int argc, char* argv[])
{
    fluorophors=new FluorophorManager;

    if (argc>1) {
        std::vector<std::string> files;
        for (int i=1; i<argc; i++) {
            std::vector<std::string> files1=listfiles_wildcard(argv[i]);
            for (size_t j=0; j<files1.size(); j++) {
                files.push_back(files1[j]);
                std::cout<<"   will simluate '"<<files1[j]<<"' ...\n";
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


    //dyn.save_trajectories();
    return 0;
}
