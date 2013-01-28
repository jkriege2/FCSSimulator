#include "dynamicsfromfiles2.h"
#include "../../../LIB/trunk/ticktock.h"
#include "../../../LIB/trunk/tools.h"


DynamicsFromFiles2::DynamicsFromFiles2(FluorophorManager* fluorophors, std::string object_name):
    FluorophorDynamics(fluorophors, object_name)
{
    shift_x=NULL;
    shift_y=NULL;
    shift_z=NULL;
    max_lines=-1;
    max_files=-1;
    trajectory_count=0;
    tmode=Sequential;
    shift_trajectories=true;
    separator_char=',';
    comment_char='#';
    randomdisplace_x_min=-1;
    randomdisplace_x_max=1;
    randomdisplace_y_min=-1;
    randomdisplace_y_max=1;
    randomdisplace_z_min=-1;
    randomdisplace_z_max=1;
    //time_factor=1;
    position_factor=1;
    abs_factor=1;
    max_columns=10;
    qfluor_factor=1;
    //time_between_trajectories=0;
    shiftmode=HalfTime;
    col_time=0;
    col_posx=1;
    col_posy=2;
    col_posz=3;
    col_px=4;
    col_py=5;
    col_pz=6;
    col_qfluor=7;
    col_abs=8;
    col_qmstate=9;
    col_px1=10;
    col_py1=11;
    col_pz1=12;
    p_fraction=1.0;
    num_start=1;
    num_stop=1;
    filenames="traj%.3d.dat";
    shiftmode=HalfTime;
    tmode=Sequential;
    parallelTrajectories=10;
}

DynamicsFromFiles2::~DynamicsFromFiles2()
{
    clear();
}



void DynamicsFromFiles2::read_config_internal(jkINIParser2& parser) {
    clear();
    FluorophorDynamics::read_config_internal(parser);
    int i=0;
    trajectory_count=0;

    shift_trajectories=parser.getSetAsBool("shift_trajectories", shift_trajectories);

    max_files=parser.getSetAsInt("max_files", max_files);
    col_time=parser.getSetAsInt("col_time", col_time);
    col_posx=parser.getSetAsInt("col_posx", col_posx);
    col_posy=parser.getSetAsInt("col_posy", col_posy);
    col_posz=parser.getSetAsInt("col_posz", col_posz);
    col_px=parser.getSetAsInt("col_px", col_px);
    col_py=parser.getSetAsInt("col_py", col_py);
    col_pz=parser.getSetAsInt("col_pz", col_pz);
    col_px1=parser.getSetAsInt("col_px1", col_px1);
    col_py1=parser.getSetAsInt("col_py1", col_py1);
    col_pz1=parser.getSetAsInt("col_pz1", col_pz1);
    p_fraction=parser.getSetAsDouble("p_fraction", p_fraction);
    col_qfluor=parser.getSetAsInt("col_qfluor", col_qfluor);
    col_abs=parser.getSetAsInt("col_abs", col_abs);
    col_qmstate=parser.getSetAsInt("col_qmstate", col_qmstate);
    parallelTrajectories=parser.getSetAsInt("parallel_trajectories", parallelTrajectories);
    randomdisplace_x_min=parser.getSetAsDouble("randomdisplace_x_min", parser.getSetAsDouble("randomdisplace_min", randomdisplace_x_min));
    randomdisplace_x_max=parser.getSetAsDouble("randomdisplace_x_max", parser.getSetAsDouble("randomdisplace_max", randomdisplace_x_max));
    randomdisplace_y_min=parser.getSetAsDouble("randomdisplace_y_min", parser.getSetAsDouble("randomdisplace_min", randomdisplace_y_min));
    randomdisplace_y_max=parser.getSetAsDouble("randomdisplace_y_max", parser.getSetAsDouble("randomdisplace_max", randomdisplace_y_max));
    randomdisplace_z_min=parser.getSetAsDouble("randomdisplace_z_min", parser.getSetAsDouble("randomdisplace_min", randomdisplace_z_min));
    randomdisplace_z_max=parser.getSetAsDouble("randomdisplace_z_max", parser.getSetAsDouble("randomdisplace_max", randomdisplace_z_max));

    // read trajectories, specified as file0=, file1=, file2=, ...
    trajectory_files.clear();
    while (parser.exists("file"+inttostr(i))) {
        trajectory_files.push_back(parser.getAsString("file"+inttostr(i)));
        trajectory_count++;
        i++;
    }
    if (i==0) { // read trajectories, specified as filenames=, num_start=, num_stop=
        num_start=parser.getSetAsInt("num_start", num_start);
        num_stop=parser.getSetAsInt("num_stop", num_stop);
        filenames=parser.getSetAsString("filename", filenames);
        std::cout<<filenames<<std::endl;
        if (filenames.find('*')!=std::string::npos || filenames.find('?')!=std::string::npos) {
            std::vector<std::string> l=listfiles_wildcard(filenames);
            int fc=l.size();
            if (max_files>0 && fc>max_files) fc=max_files;
            for (size_t j=0; j<fc; j++) {
                std::string fn=l[j];
                std::cout<<"testing "<<fn<<std::endl;
                if (file_exists(fn)) {
                    trajectory_files.push_back(fn);
                    std::cout<<"adding "<<fn<<std::endl;
                    trajectory_count++;
                }
            }
        } else {
            for (i=num_start; i<=num_stop; i++) {
                std::string fn=format(filenames, i);
                std::cout<<"testing "<<fn<<std::endl;
                if (file_exists(fn) && ((max_files<=0) || ((max_files>0) && (trajectory_count<max_files)))) {
                    trajectory_files.push_back(fn);
                    std::cout<<"adding "<<fn<<std::endl;
                    trajectory_count++;
                }
            }
        }
    }

    std::string d=chartostr(separator_char);
    if (separator_char==' ') d="space";
    if (separator_char=='\t') d="tab";
    std::string sc=strstrip(parser.getSetAsString("separator_char", d));
    if ((tolower(sc)=="space") || (sc.size()==0)) separator_char=' ';
    else if (tolower(sc)=="tab") separator_char='\t';
    else separator_char=sc[0];

    chartostr(comment_char);
    if (comment_char=='#') d="sharp";
    std::string cc=strstrip(parser.getSetAsString("comment_char", d));
    if ((tolower(cc)=="sharp")) comment_char='#';
    else comment_char=cc[0];

    //time_factor=parser.getSetAsDouble("time_factor", time_factor);
    position_factor=parser.getSetAsDouble("position_factor", position_factor);
    abs_factor=parser.getSetAsDouble("abs_factor", abs_factor);
    qfluor_factor=parser.getSetAsDouble("qfluor_factor", qfluor_factor);
    max_columns=parser.getSetAsInt("max_columns", max_columns);
    max_lines=parser.getSetAsInt("max_lines", max_lines);

    
    std::string t="sequential";
    if (tmode==Parallel) t="parallel";
    if (tmode==SequentialParallel) t="sequentialparallel";
    t=tolower(parser.getSetAsString("play_mode", t));
    if (t=="parallel") {
        tmode=Parallel;
    } else if (t=="sequentialparallel") {
        tmode=SequentialParallel;
    }
    
    
    if (shiftmode==HalfTime) t="halftime";
    else if (shiftmode==RandomDisplacedMean) t="mean_random";
    else if (shiftmode==EndEndDistanceCenter) t="end_end_center";
    else if (shiftmode==EndEndDistanceCenterRandom) t="end_end_center_random";
    else t="mean";
    std::string sm=tolower(parser.getSetAsString("shift_mode", t));
    if (sm=="halftime") {
        shiftmode=HalfTime;
    } else if (sm=="mean_random" || sm=="meanrandom" || sm=="random_mean" || sm=="randommean") {
        shiftmode=RandomDisplacedMean;
    } else if (sm=="end_end_center" || sm=="endendcenter" || sm=="endend") {
        shiftmode=EndEndDistanceCenter;
    } else if (sm=="end_end_center_random" || sm=="endendcenterrandom" || sm=="endendrandom") {
        shiftmode=EndEndDistanceCenterRandom;
    } else {
        shiftmode=Mean;
    }

    //time_between_trajectories=parser.getAsDouble("time_between_trajectories", time_between_trajectories);

    //init();
}

void DynamicsFromFiles2::init() {
    shift_x=(double*)calloc(trajectory_count, sizeof(double));
    shift_y=(double*)calloc(trajectory_count, sizeof(double));
    shift_z=(double*)calloc(trajectory_count, sizeof(double));
    file.clear();
    linecount.clear();
    columncount.clear();

    //double comx, comy, comz;
    //double minx, maxx, miny, maxy, minz, maxz;

    PublicTickTock tim, tim1;
    tim.tick();
    timing_load1=0;
    for (int i=0; i<trajectory_count; i++) {
        tim1.tick();
        std::cout<<"checking file "<<i+1<<"/"<<trajectory_count<<": "<<trajectory_files[i]<<" [sc="+chartoprintablestr(separator_char)+", cc="+chartoprintablestr(comment_char)+", ";
        //data[i].load_csv(trajectory_files[i], separator_char, comment_char);
        unsigned long long lcount=count_lines(trajectory_files[i], comment_char);
        if (max_lines>0 && lcount>max_lines) lcount=max_lines;
        if (lcount<1) throw FluorophorException(format("file '%s' does not contain data or does not exist", trajectory_files[i].c_str()));
        std::cout<<"lines="<<lcount<<"] ... "<<std::endl;

        FILE* f=fopen(trajectory_files[i].c_str(), "r");

        int j=0;
        std::vector<double> r=csv_readline(f, separator_char, comment_char);
        int ccount=r.size();
        double lastPos[3]={0,0,0};
        //double t0;
        while (r.size()>0 && ((max_lines<=0) || (j<max_lines))) {
            //if (j==0 && col_time>=0 && col_time<r.size()) t0=r[col_time];
            //if (j==1) set_sim_timestep(fabs(r[0]-t0)*time_factor);
            j++;
            if (shift_trajectories) {
                // calculate center of mass (COM) of trajectory and then shift it to 0
                if (shiftmode==Mean || shiftmode==RandomDisplacedMean) {
                    shift_x[i]=shift_x[i]+r[col_posx]*position_factor;
                    shift_y[i]=shift_y[i]+r[col_posy]*position_factor;
                    shift_z[i]=shift_z[i]+r[col_posz]*position_factor;
                } else if (shiftmode==HalfTime) {
                    if (j==lcount/2) {
                        shift_x[i]=r[col_posx]*position_factor;
                        shift_y[i]=r[col_posy]*position_factor;
                        shift_z[i]=r[col_posz]*position_factor;
                        break;
                    }
                } else if (shiftmode==EndEndDistanceCenter || shiftmode==EndEndDistanceCenterRandom) {
                    if (j<=1) {
                        shift_x[i]=r[col_posx]*position_factor;
                        shift_y[i]=r[col_posy]*position_factor;
                        shift_z[i]=r[col_posz]*position_factor;
                    }
                    lastPos[0]=r[col_posx]*position_factor;
                    lastPos[1]=r[col_posy]*position_factor;
                    lastPos[2]=r[col_posz]*position_factor;
                }
            } else {
                shift_x[i]=0;
                shift_y[i]=0;
                shift_z[i]=0;
                if (j>1) break;
            }
            /*std::cout<<j;
            for (size_t kk=0; kk<r.size(); kk++) {
                std::cout<<", "<<r[kk];
            }
            std::cout<<std::endl;*/
            r=csv_readline(f, separator_char, comment_char);
        }
        if (shift_trajectories && shiftmode==Mean) {
            shift_x[i]=shift_x[i]/double(j);
            shift_y[i]=shift_y[i]/double(j);
            shift_z[i]=shift_z[i]/double(j);
        } else if (shift_trajectories && shiftmode==RandomDisplacedMean) {
            //std::cout<<"shift_x: "<<randomdisplace_x_min<<" ... "<<randomdisplace_x_max<<std::endl;
            shift_x[i]=shift_x[i]/double(j)+gsl_ran_flat(rng, randomdisplace_x_min, randomdisplace_x_max);
            shift_y[i]=shift_y[i]/double(j)+gsl_ran_flat(rng, randomdisplace_y_min, randomdisplace_y_max);
            shift_z[i]=shift_z[i]/double(j)+gsl_ran_flat(rng, randomdisplace_z_min, randomdisplace_z_max);
        } else if (shift_trajectories && shiftmode==EndEndDistanceCenter) {
            //std::cout<<"** p1x="<<shift_x[i]<<"   p2x="<<lastPos[0]<<std::endl;
            shift_x[i]=(shift_x[i]+lastPos[0])/2.0;
            shift_y[i]=(shift_y[i]+lastPos[1])/2.0;
            shift_z[i]=(shift_z[i]+lastPos[2])/2.0;
        } else if (shift_trajectories && shiftmode==EndEndDistanceCenterRandom) {
            //std::cout<<"## p1x="<<shift_x[i]<<"   p2x="<<lastPos[0]<<std::endl;
            shift_x[i]=(shift_x[i]+lastPos[0])/2.0+gsl_ran_flat(rng, randomdisplace_x_min, randomdisplace_x_max);
            shift_y[i]=(shift_y[i]+lastPos[1])/2.0+gsl_ran_flat(rng, randomdisplace_y_min, randomdisplace_y_max);
            shift_z[i]=(shift_z[i]+lastPos[2])/2.0+gsl_ran_flat(rng, randomdisplace_z_min, randomdisplace_z_max);
        }
        fclose(f);
        f=fopen(trajectory_files[i].c_str(), "r");
        file.push_back(f);
        columncount.push_back(ccount);
        linecount.push_back(lcount);
        std::cout<<"ready!\n"<<std::endl;
        tim1.tock();
        timing_load1=timing_load1+tim1.get_duration();
    }
    tim.tock();
    timing_loadall=tim.get_duration();
    timing_load1=timing_load1/(double)trajectory_count;

    file_counter=0;
    line_counter=0;

    FluorophorDynamics::init();

    if (tmode==Sequential) {
        change_walker_count(trajectory_count);
    } else if (tmode==SequentialParallel) {
        change_walker_count(trajectory_count);
    } else { // tmode==Parallel
        change_walker_count(trajectory_count);
    }

    propagate();
}

void DynamicsFromFiles2::setAllDone() {
    for (int i=0; i<trajectory_count; i++) {
	walker_state[i].exists=false;
    }
    endoftrajectory=true;
}


void DynamicsFromFiles2::propagate(bool boundary_check) {
    FluorophorDynamics::propagate(boundary_check);
    if (file.size()<=0) return;
    endoftrajectory=trajectory_count<=0;
    if (endoftrajectory) return;
    if (tmode==Sequential) {
        for (int i=0; i<trajectory_count; i++) {
	    walker_state[i].exists= (i==file_counter);
	}
      
        int lc=linecount[file_counter]; //data[file_counter].get_line_count();
        int cc=columncount[file_counter]; //data[file_counter].get_column_count();
        //std::cout<<"fc="<<file_counter<<"   lc="<<lc<<"   cc="<<cc<<"   f="<<file[file_counter]<<"   f.size()="<<file.size()<<std::endl;
        std::vector<double> data =csv_readline(file[file_counter], separator_char, comment_char);
        walker_state[file_counter].time++;
        walker_state[file_counter].x=data[col_posx]*position_factor-shift_x[file_counter];
        walker_state[file_counter].y=data[col_posy]*position_factor-shift_y[file_counter];
        walker_state[file_counter].z=data[col_posz]*position_factor-shift_z[file_counter];
        walker_state[file_counter].p_x=init_p_x;
        walker_state[file_counter].p_y=init_p_y;
        walker_state[file_counter].p_z=init_p_z;
        walker_state[file_counter].qm_state=init_qm_state;
        walker_state[file_counter].type=init_type;
        walker_state[file_counter].spectrum=init_spectrum;
        //std::cout<<walker_state[0].x<<", "<<walker_state[0].y<<", "<<walker_state[0].z<<std::endl;
        if ( (cc>col_px) && (cc>col_py) && (cc>col_pz) &&
             (col_px>=0) && (col_py>=0) && (col_pz>=0) &&
             (max_columns>col_px) && (max_columns>col_py) && (max_columns>col_pz) ){
            walker_state[file_counter].p_x=data[col_px];
            walker_state[file_counter].p_y=data[col_py];
            walker_state[file_counter].p_z=data[col_pz];
        }
        if ( (p_fraction<1.0) &&
             (cc>col_px1) && (cc>col_py1) && (cc>col_pz1) &&
             (col_px1>=0) && (col_py1>=0) && (col_pz1>=0) &&
             (max_columns>col_px) && (max_columns>col_py) && (max_columns>col_pz) ){
            walker_state[file_counter].p_x=p_fraction*walker_state[file_counter].p_x+(1.0-p_fraction)*data[col_px1];
            walker_state[file_counter].p_y=p_fraction*walker_state[file_counter].p_y+(1.0-p_fraction)*data[col_py1];
            walker_state[file_counter].p_z=p_fraction*walker_state[file_counter].p_z+(1.0-p_fraction)*data[col_pz1];
            double sqrsum=sqrt(walker_state[file_counter].p_x*walker_state[file_counter].p_x+walker_state[file_counter].p_y*walker_state[file_counter].p_y+walker_state[file_counter].p_z*walker_state[file_counter].p_z);
            if (sqrsum>0) {
                walker_state[file_counter].p_x=walker_state[file_counter].p_x/sqrsum;
                walker_state[file_counter].p_y=walker_state[file_counter].p_y/sqrsum;
                walker_state[file_counter].p_z=walker_state[file_counter].p_z/sqrsum;
            }
        }
        // normalize p

        if ((cc>col_qfluor) && (col_qfluor>=0) && (max_columns>col_qfluor)) {
            walker_state[file_counter].q_fluor[file_counter]=data[col_qfluor]*qfluor_factor;
        }
        if ((cc>col_abs) && (col_abs>=0) && (max_columns>col_abs)) {
            walker_state[file_counter].sigma_abs[file_counter]=data[col_abs]*abs_factor;
        }
        if ((cc>col_qmstate) && (col_qmstate>=0) && (max_columns>col_qmstate)) {
            walker_state[file_counter].qm_state=(int)round(data[col_qmstate]);
        }
        //std::cout<<line_counter<<": x="<<walker_state[0].x<<" y="<<walker_state[0].y<<" z="<<walker_state[0].z<<" p_x="<<walker_state[0].p_x<<" p_y="<<walker_state[0].p_y<<" p_z="<<walker_state[0].p_z<<" q_fluor="<<walker_state[0].q_fluor<<" sigma_abs="<<walker_state[0].sigma_abs<<std::endl;
        line_counter++;
        if (line_counter>=lc) {
            if (file_counter<trajectory_count-1) {
                file_counter++;
                line_counter=0;
            } else {
                 line_counter--;
		 setAllDone();
            };
        }
    } else if (tmode==SequentialParallel) {
        
        if (file_counter<0 || file_counter+parallelTrajectories>trajectory_count) {
	    setAllDone();	  
	} else {
	    for (int i=0; i<trajectory_count; i++) {
                walker_state[i].exists= ((i>=file_counter) && (i<file_counter+parallelTrajectories));
	    }
	    
	    bool anyRunning=false;
	    
	    for (int f=file_counter; f<file_counter+parallelTrajectories; f++) {
		int lc=linecount[f]; //data[file_counter].get_line_count();
		int cc=columncount[f]; //data[file_counter].get_column_count();
		//std::cout<<"fc="<<file_counter<<"   lc="<<lc<<"   cc="<<cc<<"   f="<<file[file_counter]<<"   f.size()="<<file.size()<<std::endl;
		if (line_counter<lc) {
                    std::vector<double> data =csv_readline(file[f], separator_char, comment_char);
                    walker_state[f].time++;
                    walker_state[f].x=data[col_posx]*position_factor-shift_x[f];
                    walker_state[f].y=data[col_posy]*position_factor-shift_y[f];
                    walker_state[f].z=data[col_posz]*position_factor-shift_z[f];
                    walker_state[f].p_x=init_p_x;
                    walker_state[f].p_y=init_p_y;
                    walker_state[f].p_z=init_p_z;
                    walker_state[f].qm_state=init_qm_state;
                    walker_state[f].type=init_type;
                    walker_state[f].spectrum=init_spectrum;
                    //std::cout<<walker_state[0].x<<", "<<walker_state[0].y<<", "<<walker_state[0].z<<std::endl;
                    if ( (cc>col_px) && (cc>col_py) && (cc>col_pz) &&
                        (col_px>=0) && (col_py>=0) && (col_pz>=0) &&
                        (max_columns>col_px) && (max_columns>col_py) && (max_columns>col_pz) ){
                    walker_state[f].p_x=data[col_px];
                    walker_state[f].p_y=data[col_py];
                    walker_state[f].p_z=data[col_pz];
                    }
                    if ( (p_fraction<1.0) &&
                        (cc>col_px1) && (cc>col_py1) && (cc>col_pz1) &&
                        (col_px1>=0) && (col_py1>=0) && (col_pz1>=0) &&
                        (max_columns>col_px) && (max_columns>col_py) && (max_columns>col_pz) ){
                    walker_state[f].p_x=p_fraction*walker_state[f].p_x+(1.0-p_fraction)*data[col_px1];
                    walker_state[f].p_y=p_fraction*walker_state[f].p_y+(1.0-p_fraction)*data[col_py1];
                    walker_state[f].p_z=p_fraction*walker_state[f].p_z+(1.0-p_fraction)*data[col_pz1];
                    double sqrsum=sqrt(walker_state[f].p_x*walker_state[f].p_x+walker_state[f].p_y*walker_state[f].p_y+walker_state[f].p_z*walker_state[f].p_z);
                    if (sqrsum>0) {
                        walker_state[f].p_x=walker_state[f].p_x/sqrsum;
                        walker_state[f].p_y=walker_state[f].p_y/sqrsum;
                        walker_state[f].p_z=walker_state[f].p_z/sqrsum;
                    }
                    }
                    // normalize p

                    if ((cc>col_qfluor) && (col_qfluor>=0) && (max_columns>col_qfluor)) {
                        walker_state[f].q_fluor[f]=data[col_qfluor]*qfluor_factor;
                    }
                    if ((cc>col_abs) && (col_abs>=0) && (max_columns>col_abs)) {
                        walker_state[f].sigma_abs[f]=data[col_abs]*abs_factor;
                    }
                    if ((cc>col_qmstate) && (col_qmstate>=0) && (max_columns>col_qmstate)) {
                        walker_state[f].qm_state=(int)round(data[col_qmstate]);
                    }
                    //std::cout<<line_counter<<": x="<<walker_state[0].x<<" y="<<walker_state[0].y<<" z="<<walker_state[0].z<<" p_x="<<walker_state[0].p_x<<" p_y="<<walker_state[0].p_y<<" p_z="<<walker_state[0].p_z<<" q_fluor="<<walker_state[0].q_fluor<<" sigma_abs="<<walker_state[0].sigma_abs<<std::endl;
                    
                    anyRunning=true;
                } 
	    }
	    line_counter++;
            if (!anyRunning) {
                file_counter=file_counter+parallelTrajectories;
                line_counter=0;
                if (file_counter+parallelTrajectories>=trajectory_count) {
                    setAllDone();
                }
            }
	}
    } else { // tmode==Parallel
        bool anyon=false;
        for (file_counter=0; file_counter<trajectory_count; file_counter++) {
            int lc=linecount[file_counter]; //data[file_counter].get_line_count();
            int cc=columncount[file_counter]; //data[file_counter].get_column_count();
            std::vector<double> data =csv_readline(file[file_counter], separator_char, comment_char);
            if (line_counter<lc) {
                walker_state[file_counter].time++;
                walker_state[file_counter].x=data[col_posx]*position_factor-shift_x[file_counter];
                walker_state[file_counter].y=data[col_posy]*position_factor-shift_y[file_counter];
                walker_state[file_counter].z=data[col_posz]*position_factor-shift_z[file_counter];
                walker_state[file_counter].p_x=init_p_x;
                walker_state[file_counter].p_y=init_p_y;
                walker_state[file_counter].p_z=init_p_z;
                walker_state[file_counter].qm_state=init_qm_state;
                walker_state[file_counter].type=init_type;
                walker_state[file_counter].spectrum=init_spectrum;
                //std::cout<<line_counter<<", "<<lc<<", x="<<walker_state[file_counter].x<<", y="<<walker_state[file_counter].y<<", z="<<walker_state[file_counter].z<<std::endl;
                if ( (cc>col_px) && (cc>col_py) && (cc>col_pz) &&
                     (col_px>=0) && (col_py>=0) && (col_pz>=0) &&
                     (max_columns>col_px) && (max_columns>col_py) && (max_columns>col_pz) ){
                    walker_state[file_counter].p_x=data[col_px];
                    walker_state[file_counter].p_y=data[col_py];
                    walker_state[file_counter].p_z=data[col_pz];
                }
                if ( (p_fraction<1.0) &&
                     (cc>col_px1) && (cc>col_py1) && (cc>col_pz1) &&
                     (col_px1>=0) && (col_py1>=0) && (col_pz1>=0) &&
                     (max_columns>col_px) && (max_columns>col_py) && (max_columns>col_pz) ){
                    walker_state[file_counter].p_x=p_fraction*walker_state[file_counter].p_x+(1.0-p_fraction)*data[col_px1];
                    walker_state[file_counter].p_y=p_fraction*walker_state[file_counter].p_y+(1.0-p_fraction)*data[col_py1];
                    walker_state[file_counter].p_z=p_fraction*walker_state[file_counter].p_z+(1.0-p_fraction)*data[col_pz1];
                    double sqrsum=sqrt(walker_state[0].p_x*walker_state[0].p_x+walker_state[0].p_y*walker_state[0].p_y+walker_state[0].p_z*walker_state[0].p_z);
                    if (sqrsum>0) {
                        walker_state[0].p_x=walker_state[0].p_x/sqrsum;
                        walker_state[0].p_y=walker_state[0].p_y/sqrsum;
                        walker_state[0].p_z=walker_state[0].p_z/sqrsum;
                    }
                }
                // normalize p
                if ((cc>col_qfluor) && (col_qfluor>=0) && (max_columns>col_qfluor)) {
                    walker_state[file_counter].q_fluor[0]=data[col_qfluor]*qfluor_factor;
                }
                if ((cc>col_abs) && (col_abs>=0) && (max_columns>col_abs)) {
                    walker_state[file_counter].sigma_abs[0]=data[col_abs]*abs_factor;
                }
                if ((cc>col_qmstate) && (col_qmstate>=0) && (max_columns>col_qmstate)) {
                    walker_state[file_counter].qm_state=(int)round(data[col_qmstate]);
                }
            } else {
                walker_state[file_counter].exists=false; // deactivate fluorophor, if no more data exists
            }
            anyon = anyon || walker_state[file_counter].exists;
        }
        if (anyon) line_counter++; // if there is still some data available, increase linecount
        else endoftrajectory=true; // if there is no more data, signal this with     endoftrajectory=true!
    }
}

std::string DynamicsFromFiles2::report() {
    std::string s=FluorophorDynamics::report();
    s+="trajectory_count = "+inttostr(trajectory_count)+"\n";
    s+="trajectory_files = ";
    for (int i=0; i<trajectory_count; i++) {
        if (i>0) s+=",\n                   ";
        s+=trajectory_files[i]+"  [ lines="+inttostr(linecount[i])+" columns="+inttostr(columncount[i])+"shift_x="+floattostr(shift_x[i])+" shift_y="+floattostr(shift_y[i])+" shift_z="+floattostr(shift_z[i])+" ] microns";
    }
    s+="\n";
    s+="separator_char = "+chartoprintablestr(separator_char)+"\n";
    s+="comment_char = "+chartoprintablestr(comment_char)+"\n";
    if (shiftmode==Mean) s+="shift_mode = mean (center-of-mass)\n";
    if (shiftmode==RandomDisplacedMean) {
        s+="shift_mode = mean_random (randomly displaced center-of-mass)\n";
        s+="mean_shift_range_x = "+floattostr(randomdisplace_x_min)+" ... "+floattostr(randomdisplace_x_max)+" micron\n";
        s+="mean_shift_range_y = "+floattostr(randomdisplace_y_min)+" ... "+floattostr(randomdisplace_y_max)+" micron\n";
        s+="mean_shift_range_z = "+floattostr(randomdisplace_z_min)+" ... "+floattostr(randomdisplace_z_max)+" micron\n";
    }
    if (shiftmode==EndEndDistanceCenter) {
        s+="shift_mode = end_end_center (shift center of end-to-end line to 0)\n";
    }
    if (shiftmode==EndEndDistanceCenterRandom) {
        s+="shift_mode = end_end_center_random (shift center of end-to-end line to random position)\n";
        s+="mean_shift_range_x = "+floattostr(randomdisplace_x_min)+" ... "+floattostr(randomdisplace_x_max)+" micron\n";
        s+="mean_shift_range_y = "+floattostr(randomdisplace_y_min)+" ... "+floattostr(randomdisplace_y_max)+" micron\n";
        s+="mean_shift_range_z = "+floattostr(randomdisplace_z_min)+" ... "+floattostr(randomdisplace_z_max)+" micron\n";
    }
    if (shiftmode==HalfTime) s+="shift_mode = halftime (to position at file_duration/2)\n";
    s+="timing_loadall = "+floattostr(timing_loadall)+" secs\n";
    s+="timing_load1 = "+floattostr(timing_load1)+" secs\n";
    s+="max_files = "+inttostr(max_files)+"\n";
    s+="max_columns = "+inttostr(max_columns)+"\n";
    if (max_lines<=0) {
        s+="max_lines = [all lines]\n";
    } else {
        s+="max_lines = "+inttostr(max_lines)+"\n";
    }
    s+="col_time = "+inttostr(col_time)+"\n";
    s+="col_posx, col_posy, col_posz = "+inttostr(col_posx)+", "+inttostr(col_posy)+", "+inttostr(col_posz)+"\n";
    s+="col_px, col_py, col_pz = "+inttostr(col_px)+", "+inttostr(col_py)+", "+inttostr(col_pz)+"\n";
    s+="col_px1, col_py1, col_pz1 = "+inttostr(col_px1)+", "+inttostr(col_py1)+", "+inttostr(col_pz1)+"\n";
    s+="p_fraction = "+floattostr(p_fraction)+"\n";
    s+="col_abs = "+inttostr(col_abs)+"\n";
    s+="col_qmstate = "+inttostr(col_qmstate)+"\n";
    s+="col_qfluor = "+inttostr(col_qfluor)+"\n";
    if (tmode==Sequential) s+="play_mode = sequential\n";
    if (tmode==SequentialParallel) s+="play_mode = sequential_parallel\nparallel_trajectories = "+inttostr(parallelTrajectories)+"\n";
    if (tmode==Parallel) s+="play_mode = parallel\n";
    //s+="max_columns = "+inttostr(max_columns)+"\n";
    //s+="time_factor = "+floattostr(time_factor)+"\n";
    s+="position_factor = "+floattostr(position_factor)+"\n";
    s+="abs_factor = "+floattostr(abs_factor)+"\n";
    s+="qfluor_factor = "+floattostr(qfluor_factor)+"\n";
    s+="estimated_runtime = "+floattostr(estimate_runtime())+" secs\n";
    return s;
}

void DynamicsFromFiles2::clear() {
   if (file.size()>0) {
        for (size_t i=0; i<file.size(); i++) {
            fclose(file[i]);
        }
        file.clear();
    }
    if (shift_x!=NULL) free(shift_x);
    if (shift_y!=NULL) free(shift_y);
    if (shift_z!=NULL) free(shift_z);
    shift_x=NULL;
    shift_y=NULL;
    shift_z=NULL;
    columncount.clear();
    linecount.clear();
    trajectory_files.clear();
    trajectory_count=0;
}

double DynamicsFromFiles2::estimate_runtime() {
    double rt=0;
    if (tmode==Sequential) {
        std::cout<<"sequential runtime estimation";
        for (int fc=0; fc<trajectory_count; fc++) {
            int lc=linecount[fc];
            rt=rt+double(lc)*sim_timestep;
        }
    } else if (tmode==SequentialParallel) {
        std::cout<<"sequential parallel runtime estimation (parallelTrajectories="<<parallelTrajectories<<"  trajectory_count="<<trajectory_count<<")";
        for (int fc=0; fc<trajectory_count; fc=fc+parallelTrajectories) {
            double llc=0;
            if (fc+parallelTrajectories<=trajectory_count) {
                for (int kk=fc; kk<fc+parallelTrajectories; kk++) {
                    int lc=linecount[fc];
                    if (lc>llc) llc=lc;
                }
            }
            rt=rt+double(llc)*sim_timestep;
        }
    } else {
        std::cout<<"parallel runtime estimation";
        for (int fc=0; fc<trajectory_count; fc++) {
            int lc=linecount[fc];
            double d=double(lc)*sim_timestep;
            if (d>rt) rt=d;
        }
    }
    std::cout<<"   => runtime="<<rt<<std::endl;
    return rt;
}
