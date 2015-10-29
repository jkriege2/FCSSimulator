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


#include "trajectoryplot.h"

TrajectoryPlotMeasurement::TrajectoryPlotMeasurement(FluorophorManager* fluorophors, std::string objectname):
    FluorescenceMeasurement(fluorophors, objectname)
{
    trajectories_to_plot=5;
    max_steps=-1;
    avg_steps=100;
    avgCount=0;
}

TrajectoryPlotMeasurement::~TrajectoryPlotMeasurement()
{

    clear();
}


void TrajectoryPlotMeasurement::init() {
    FluorescenceMeasurement::init();
    clear();
    avgCount=0;
}

void TrajectoryPlotMeasurement::propagate() {
    tick();
    FluorescenceMeasurement::propagate();

    bool runThis=true;

    if (trajectories.size()>0)  {
        if (max_steps>0 && (int64_t)trajectories[0].size()>=max_steps) runThis=false;
    }

    if (runThis) {
        int t=0;
        for (size_t d=0; d<dyn.size(); d++) { // go through all dynamics objects that provide data for this measurement object
            const FluorophorDynamics::walkerState* dynn=dyn[d]->get_walker_state();
            unsigned long wc=dyn[d]->get_walker_count();
            if (!dyn[d]->end_of_trajectory_reached()) {
                for (unsigned long i=0; i<wc; i++) { // iterate through all walkers in the d-th dynamics object
                    if (dynn[i].exists) {
                        const float x=dynn[i].x;
                        const float y=dynn[i].y;
                        const float z=dynn[i].z;
                        TrajectoryPlotMeasurement::trajectory_info ti;
                        ti.t=0;
                        ti.x=0;
                        ti.y=0;
                        ti.z=0;
                        if (avgCount>0) ti=currentT[t];
                        ti.t=ti.t+get_sim_time();
                        ti.x=ti.x+x;
                        ti.y=ti.y+y;
                        ti.z=ti.z+z;
                        currentT[t]=ti;
                    }
                    t++;
                    if (t>=trajectories_to_plot) break;
                }
                if (t>=trajectories_to_plot) break;
            }
        }
        avgCount++;
        if (avgCount>=avg_steps) {
            for (int i=0; i<trajectories_to_plot; i++) {
                TrajectoryPlotMeasurement::trajectory_info ti=currentT[i];
                ti.t=ti.t/float(avgCount);
                ti.x=ti.x/float(avgCount);
                ti.y=ti.y/float(avgCount);
                ti.z=ti.z/float(avgCount);
                if (i>=(int64_t)trajectories.size())  {
                    std::vector<trajectory_info> inf;
                    inf.push_back(ti);
                    trajectories.push_back(inf);
                } else {
                    trajectories[i].push_back(ti);
                }
                ti.t=0;
                ti.x=0;
                ti.y=0;
                ti.z=0;
                currentT[i]=ti;
            }
            avgCount=0;
        }
    }


    tock();
    runtime=runtime+get_duration();
}

std::string TrajectoryPlotMeasurement::report() {
    std::string s=FluorescenceMeasurement::report();
    s+="trajectories = "+inttostr(trajectories_to_plot)+"\n";
    s+="avg_steps    = "+inttostr(avg_steps)+"\n";
    s+="max_steps    = "+inttostr(max_steps)+"\n";

    return s;
}

std::string TrajectoryPlotMeasurement::dot_get_properties() {
    std::string s=FluorescenceMeasurement::dot_get_properties();
    s+="trajectories = "+inttostr(trajectories_to_plot)+"<BR/>";
    s+="avg_steps    = "+inttostr(avg_steps)+"<BR/>";
    s+="max_steps    = "+inttostr(max_steps)+"<BR/>";

    return s;
}
void TrajectoryPlotMeasurement::save() {
    FILE* f;
    char fn[255];

    if (trajectories_to_plot<=0) return;

    sprintf(fn, "%s%straj.dat", basename.c_str(), object_name.c_str());
    size_t steps=0;
    int rowItems=0;
    float maxT=0;
    if (trajectories.size()>0) steps=trajectories[0].size();
    if (steps>0) {
        std::cout<<"writing '"<<fn<<"' ...";
        f=fopen(fn, "wb");
        rowItems=3*trajectories.size()+1;
        float* d=(float*)malloc(rowItems*sizeof(float));
        for (size_t t=0; t<steps; t++) {
            for (size_t i=0; i<trajectories.size(); i++) {
                TrajectoryPlotMeasurement::trajectory_info ti=trajectories[i].at(t);
                d[0]=ti.t;
                d[1+i*3]=ti.x;
                d[1+i*3+1]=ti.y;
                d[1+i*3+2]=ti.z;
                if (ti.t>maxT) maxT=ti.t;
            }
            fwrite(d,rowItems*sizeof(float),1,f);
        }
        free(d);
        fclose(f);
        std::cout<<" done!\n";
    }
    std::string trajfn=fn;

    if (rowItems>0 && steps>0) {
        sprintf(fn, "%s%strajplot.plt", basename.c_str(), object_name.c_str());
        std::cout<<"writing '"<<fn<<"' ...";
        f=fopen(fn, "w");
        fprintf(f, "set xlabel \"x [micron]\"\n");
        fprintf(f, "set ylabel \"y [micron]\"\n");
        fprintf(f, "set zlabel \"z [micron]\"\n");
        fprintf(f, "set size square\n");
        fprintf(f, "h(x,y) = (x*x+y*y<=r)?(1/0):(gamma*sqrt(r*r - (x-x0)**2 - (y-y0)**2 ) + z0)\n");
        fprintf(f, "x0 = 0\n");
        fprintf(f, "y0 = 0\n");
        fprintf(f, "z0 = 0\n");
        fprintf(f, "r=1\n");
        fprintf(f, "gamma=0\n");
        fprintf(f, "set isosamples 100\n");

        for (int plt=0; plt<2; plt++) {
            if (plt==0) {
                fprintf(f, "set terminal pdfcairo color solid font \"%s, 7\" linewidth 2 size 20cm,15cm\n", GNUPLOT_FONT);
                fprintf(f, "set output \"%s\"\n", extract_file_name(basename+object_name+"trajplot.pdf").c_str());
            } else if (plt==1) {
                fprintf(f, "set terminal wxt font \"%s, 8\"\n", GNUPLOT_FONT);
                fprintf(f, "set output\n");
            }
            fprintf(f, "set title \"%s:    all trajectories\"\n", object_name.c_str());
            for (size_t i=0; i<trajectories.size(); i++) {
                if (i>0) fprintf(f, ", \\\n");
                else fprintf(f, "splot ");
                fprintf(f, " \"%s\" binary format=\"%%%dfloat\" using %d:%d:%d title\"trajectory %d\" with lines", extract_file_name(trajfn).c_str(), rowItems, 1+1+int(i)*3, 1+1+int(i)*3+1, 1+1+int(i)*3+2, int(i));
            }
            fprintf(f, "\n\n");
            if (plt==1) fprintf(f, "pause -1\n");

            for (size_t i=0; i<trajectories.size(); i++) {
                fprintf(f, "set title \"%s:    trajectory %d\"\n", object_name.c_str(), int(i));
                fprintf(f, "splot \"%s\" binary format=\"%%%dfloat\" using %d:%d:%d:1 with lines palette \n#,h(x,y) title \"r=1micron sphere around (0,0,0)\" with lines ls 1, -h(x,y) notitle with lines ls 1\n", extract_file_name(trajfn).c_str(), rowItems, 1+1+int(i)*3, 1+1+int(i)*3+1, 1+1+int(i)*3+2);
                if (plt==1) fprintf(f, "pause -1\n");
            }
        }

        fclose(f);
        std::cout<<" done!\n";
    }
}

void TrajectoryPlotMeasurement::read_config_internal(jkINIParser2& parser) {
    FluorescenceMeasurement::read_config_internal(parser);

    trajectories_to_plot=parser.getAsInt("trajectories", trajectories_to_plot);
    avg_steps=parser.getAsInt("avg_steps", avg_steps);
    max_steps=parser.getAsInt("max_steps", max_steps);

    //init();
}

void TrajectoryPlotMeasurement::clear() {
    trajectories.clear();
    currentT.clear();
    avgCount=0;
}
