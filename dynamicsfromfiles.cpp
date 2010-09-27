#include "dynamicsfromfiles.h"
#include "..\lib\ticktock.h"


DynamicsFromFiles::DynamicsFromFiles(FluorophorManager* fluorophors):
    FluorophorDynamics(fluorophors)
{
    data=NULL;
    shift_x=NULL;
    shift_y=NULL;
    shift_z=NULL;
    trajectory_count=0;
    tmode=Sequential;
    shift_trajectories=true;
    separator_char=',';
    comment_char='#';
    time_factor=1;
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

}

DynamicsFromFiles::~DynamicsFromFiles()
{
    if (data!=NULL) {
        if (trajectory_count>1) delete[] data;
        else delete data;
    }
    if (shift_x!=NULL) free(shift_x);
    if (shift_y!=NULL) free(shift_y);
    if (shift_z!=NULL) free(shift_z);
    data=NULL;
    shift_x=NULL;
    shift_y=NULL;
    shift_z=NULL;
    trajectory_files.clear();
}


void DynamicsFromFiles::read_config_internal(jkINIParser2& parser) {
    clear();
    FluorophorDynamics::read_config(parser);
    int i=0;
    trajectory_count=0;

    shift_trajectories=parser.getSetAsBool("shift_trajectories", shift_trajectories);

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

    // read trajectories, specified as file0=, file1=, file2=, ...
    while (parser.exists("file"+inttostr(i))) {
        trajectory_files.push_back(parser.getAsString("file"+inttostr(i)));
        trajectory_count++;
        i++;
    }
    if (i==0) { // read trajectories, specified as filenames=, num_start=, num_stop=
        num_start=parser.getSetAsInt("num_start", 1);
        num_stop=parser.getSetAsInt("num_stop", 1);
        std::string filenames=parser.getSetAsString("filenames", "traj%.3d.dat");
        for (i=num_start; i<=num_stop; i++) {
            std::string fn=format(filenames, i);
            std::cout<<"testing "<<fn<<std::endl;
            if (file_exists(fn)) {
                trajectory_files.push_back(fn);
                std::cout<<"adding "<<fn<<std::endl;
                trajectory_count++;
            }
        }
    }

    std::string d=chartostr(comment_char);
    if (comment_char==' ') d="space";
    if (comment_char=='\t') d="tab";
    std::string sc=strstrip(parser.getSetAsString("separator_char", d));
    if ((tolower(sc)=="space") || (sc.size()==0)) separator_char=' ';
    else if (tolower(sc)=="tab") separator_char='\t';
    else separator_char=sc[0];

    chartostr(comment_char);
    if (comment_char=='#') d="sharp";
    std::string cc=strstrip(parser.getSetAsString("comment_char", d));
    if ((tolower(cc)=="sharp")) comment_char='#';
    else comment_char=cc[0];

    time_factor=parser.getSetAsDouble("time_factor", time_factor);
    position_factor=parser.getSetAsDouble("position_factor", position_factor);
    abs_factor=parser.getSetAsDouble("abs_factor", abs_factor);
    qfluor_factor=parser.getSetAsDouble("qfluor_factor", qfluor_factor);
    max_columns=parser.getSetAsInt("max_columns", max_columns);

    tmode=Sequential;
    if (tolower(parser.getSetAsString("play_mode", "sequential"))=="parallel") {
        tmode=Parallel;
    }
    shiftmode=HalfTime;
    if (tolower(parser.getSetAsString("shift_mode", "halftime"))=="mean") {
        shiftmode=Mean;
    }

    //time_between_trajectories=parser.getAsDouble("time_between_trajectories", time_between_trajectories);

    init();
}

void DynamicsFromFiles::init() {
    data=new datatable[trajectory_count];
    shift_x=(double*)calloc(trajectory_count, sizeof(double));
    shift_y=(double*)calloc(trajectory_count, sizeof(double));
    shift_z=(double*)calloc(trajectory_count, sizeof(double));

    double comx, comy, comz;
    double minx, maxx, miny, maxy, minz, maxz;

    PublicTickTock tim, tim1;
    tim.tick();
    timing_load1=0;
    for (int i=0; i<trajectory_count; i++) {
        tim1.tick();
        std::cout<<"loading file "<<trajectory_files[i]<<" [sc="+chartoprintablestr(separator_char)+", cc="+chartoprintablestr(comment_char)+"]...";
        data[i].load_csv(trajectory_files[i], separator_char, comment_char);
        std::cout<<"ready!    [lines="<<data[i].get_line_count()<<"  columns="<<data[i].get_column_count()<<"]"<<std::endl;
        if (i==0) set_sim_timestep(fabs(data[i].get(col_time,1)-data[i].get(col_time,0))*time_factor);
        if (shift_trajectories) {
            // calculate center of mass (COM) of trajectory and then shift it to 0
            if (shiftmode==Mean) {
                shift_x[i]=data[i].column_average(col_posx)*position_factor;
                shift_y[i]=data[i].column_average(col_posy)*position_factor;
                shift_z[i]=data[i].column_average(col_posz)*position_factor;
            } else if (shiftmode==HalfTime) {
                shift_x[i]=data[i].get(col_posx, data[i].get_line_count()/2)*position_factor;
                shift_y[i]=data[i].get(col_posy, data[i].get_line_count()/2)*position_factor;
                shift_z[i]=data[i].get(col_posz, data[i].get_line_count()/2)*position_factor;
            }

        } else {
            shift_x[i]=0;
            shift_y[i]=0;
            shift_z[i]=0;
        }
        std::cout<<" .\n";
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
        change_walker_count(1);
    } else { // tmode==Parallel
        change_walker_count(trajectory_count);
    }

    propagate();
}

void DynamicsFromFiles::propagate(bool boundary_check) {
    if (data==NULL) return;
    if (tmode==Sequential) {
        int lc=data[file_counter].get_line_count();
        int cc=data[file_counter].get_column_count();
        walker_state[0].time++;
        walker_state[0].x=data[file_counter].get(col_posx, line_counter)*position_factor-shift_x[file_counter];
        walker_state[0].y=data[file_counter].get(col_posy, line_counter)*position_factor-shift_y[file_counter];
        walker_state[0].z=data[file_counter].get(col_posz, line_counter)*position_factor-shift_z[file_counter];
        walker_state[0].p_x=init_p_x;
        walker_state[0].p_y=init_p_y;
        walker_state[0].p_z=init_p_z;
        walker_state[0].q_fluor=init_q_fluor;
        walker_state[0].tau_fl=init_tau_fl;
        walker_state[0].qm_state=init_qm_state;
        walker_state[0].sigma_abs=init_sigma_abs;
        walker_state[0].type=init_type;
        walker_state[0].spectrum=init_spectrum;
        //std::cout<<walker_state[0].x<<", "<<walker_state[0].y<<", "<<walker_state[0].z<<std::endl;
        if ( (cc>col_px) && (cc>col_py) && (cc>col_pz) &&
             (col_px>=0) && (col_py>=0) && (col_pz>=0) &&
             (max_columns>col_px) && (max_columns>col_py) && (max_columns>col_pz) ){
            walker_state[0].p_x=data[file_counter].get(col_px, line_counter);
            walker_state[0].p_y=data[file_counter].get(col_py, line_counter);
            walker_state[0].p_z=data[file_counter].get(col_pz, line_counter);
        }
        if ( (p_fraction<1.0) &&
             (cc>col_px1) && (cc>col_py1) && (cc>col_pz1) &&
             (col_px1>=0) && (col_py1>=0) && (col_pz1>=0) &&
             (max_columns>col_px) && (max_columns>col_py) && (max_columns>col_pz) ){
            walker_state[0].p_x=p_fraction*walker_state[0].p_x+(1.0-p_fraction)*data[file_counter].get(col_px1, line_counter);
            walker_state[0].p_y=p_fraction*walker_state[0].p_y+(1.0-p_fraction)*data[file_counter].get(col_py1, line_counter);
            walker_state[0].p_z=p_fraction*walker_state[0].p_z+(1.0-p_fraction)*data[file_counter].get(col_pz1, line_counter);
            double sqrsum=sqrt(walker_state[0].p_x*walker_state[0].p_x+walker_state[0].p_y*walker_state[0].p_y+walker_state[0].p_z*walker_state[0].p_z);
            if (sqrsum>0) {
                walker_state[0].p_x=walker_state[0].p_x/sqrsum;
                walker_state[0].p_y=walker_state[0].p_y/sqrsum;
                walker_state[0].p_z=walker_state[0].p_z/sqrsum;
            }
        }
        // normalize p

        if ((cc>col_qfluor) && (col_qfluor>=0) && (max_columns>col_qfluor)) {
            walker_state[0].q_fluor=data[file_counter].get(col_qfluor, line_counter)*qfluor_factor;
        }
        if ((cc>col_abs) && (col_abs>=0) && (max_columns>col_abs)) {
            walker_state[0].sigma_abs=data[file_counter].get(col_abs, line_counter)*abs_factor;
        }
        if ((cc>col_qmstate) && (col_qmstate>=0) && (max_columns>col_qmstate)) {
            walker_state[0].qm_state=(int)round(data[file_counter].get(col_qmstate, line_counter));
        }
        //std::cout<<line_counter<<": x="<<walker_state[0].x<<" y="<<walker_state[0].y<<" z="<<walker_state[0].z<<" p_x="<<walker_state[0].p_x<<" p_y="<<walker_state[0].p_y<<" p_z="<<walker_state[0].p_z<<" q_fluor="<<walker_state[0].q_fluor<<" sigma_abs="<<walker_state[0].sigma_abs<<std::endl;
        line_counter++;
        if (line_counter>=lc) {
            if (file_counter<trajectory_count-1) {
                file_counter++;
                line_counter=0;
            } else line_counter--;
        }
    } else { // tmode==Parallel
        for (file_counter=0; file_counter<trajectory_count; file_counter++) {
            int lc=data[file_counter].get_line_count();
            int cc=data[file_counter].get_column_count();
            if (line_counter<lc) {
                walker_state[file_counter].time++;
                walker_state[file_counter].x=data[file_counter].get(col_posx, line_counter)*position_factor-shift_x[file_counter];
                walker_state[file_counter].y=data[file_counter].get(col_posy, line_counter)*position_factor-shift_y[file_counter];
                walker_state[file_counter].z=data[file_counter].get(col_posz, line_counter)*position_factor-shift_z[file_counter];
                walker_state[file_counter].p_x=init_p_x;
                walker_state[file_counter].p_y=init_p_y;
                walker_state[file_counter].p_z=init_p_z;
                walker_state[file_counter].q_fluor=init_q_fluor;
                walker_state[file_counter].tau_fl=init_tau_fl;
                walker_state[file_counter].qm_state=init_qm_state;
                walker_state[file_counter].sigma_abs=init_sigma_abs;
                walker_state[file_counter].type=init_type;
                walker_state[file_counter].spectrum=init_spectrum;
                //std::cout<<walker_state[file_counter].x<<", "<<walker_state[file_counter].y<<", "<<walker_state[file_counter].z<<std::endl;
                if ( (cc>col_px) && (cc>col_py) && (cc>col_pz) &&
                     (col_px>=0) && (col_py>=0) && (col_pz>=0) &&
                     (max_columns>col_px) && (max_columns>col_py) && (max_columns>col_pz) ){
                    walker_state[file_counter].p_x=data[file_counter].get(col_px, line_counter);
                    walker_state[file_counter].p_y=data[file_counter].get(col_py, line_counter);
                    walker_state[file_counter].p_z=data[file_counter].get(col_pz, line_counter);
                }
                if ( (p_fraction<1.0) &&
                     (cc>col_px1) && (cc>col_py1) && (cc>col_pz1) &&
                     (col_px1>=0) && (col_py1>=0) && (col_pz1>=0) &&
                     (max_columns>col_px) && (max_columns>col_py) && (max_columns>col_pz) ){
                    walker_state[file_counter].p_x=p_fraction*walker_state[file_counter].p_x+(1.0-p_fraction)*data[file_counter].get(col_px1, line_counter);
                    walker_state[file_counter].p_y=p_fraction*walker_state[file_counter].p_y+(1.0-p_fraction)*data[file_counter].get(col_py1, line_counter);
                    walker_state[file_counter].p_z=p_fraction*walker_state[file_counter].p_z+(1.0-p_fraction)*data[file_counter].get(col_pz1, line_counter);
                    double sqrsum=sqrt(walker_state[0].p_x*walker_state[0].p_x+walker_state[0].p_y*walker_state[0].p_y+walker_state[0].p_z*walker_state[0].p_z);
                    if (sqrsum>0) {
                        walker_state[0].p_x=walker_state[0].p_x/sqrsum;
                        walker_state[0].p_y=walker_state[0].p_y/sqrsum;
                        walker_state[0].p_z=walker_state[0].p_z/sqrsum;
                    }
                }
                // normalize p
                if ((cc>col_qfluor) && (col_qfluor>=0) && (max_columns>col_qfluor)) {
                    walker_state[file_counter].q_fluor=data[file_counter].get(col_qfluor, line_counter)*qfluor_factor;
                }
                if ((cc>col_abs) && (col_abs>=0) && (max_columns>col_abs)) {
                    walker_state[file_counter].sigma_abs=data[file_counter].get(col_abs, line_counter)*abs_factor;
                }
                if ((cc>col_qmstate) && (col_qmstate>=0) && (max_columns>col_qmstate)) {
                    walker_state[file_counter].qm_state=(int)round(data[file_counter].get(col_qmstate, line_counter));
                }
            }
        }
        if (line_counter<data[0].get_line_count()-1) line_counter++;
    }
}

std::string DynamicsFromFiles::report() {
    std::string s=FluorophorDynamics::report();
    s+="trajectory_count = "+inttostr(trajectory_count)+"\n";
    s+="trajectory_files = ";
    for (int i=0; i<trajectory_count; i++) {
        if (i>0) s+=",\n                   ";
        s+=trajectory_files[i]+"  [ shift_x="+floattostr(shift_x[i])+" shift_y="+floattostr(shift_y[i])+" shift_z="+floattostr(shift_z[i])+" ] microns";
    }
    s+="\n";
    if (shiftmode==Mean) s+="shift_mode = mean (to center-of-mass)\n";
    if (shiftmode==HalfTime) s+="shift_mode = halftime (to position at file_duration/2)\n";
    s+="timing_loadall = "+floattostr(timing_loadall)+" secs\n";
    s+="timing_load1 = "+floattostr(timing_load1)+" secs\n";
    s+="max_columns = "+inttostr(max_columns)+"\n";
    s+="col_time = "+inttostr(col_time)+"\n";
    s+="col_posx, col_posy, col_posz = "+inttostr(col_posx)+", "+inttostr(col_posy)+", "+inttostr(col_posz)+"\n";
    s+="col_px, col_py, col_pz = "+inttostr(col_px)+", "+inttostr(col_py)+", "+inttostr(col_pz)+"\n";
    s+="col_px1, col_py1, col_pz1 = "+inttostr(col_px1)+", "+inttostr(col_py1)+", "+inttostr(col_pz1)+"\n";
    s+="p_fraction = "+floattostr(p_fraction)+"\n";
    s+="col_abs = "+inttostr(col_abs)+"\n";
    s+="col_qmstate = "+inttostr(col_qmstate)+"\n";
    s+="col_qfluor = "+inttostr(col_qfluor)+"\n";
    if (tmode==Sequential) s+="play_mode = sequential\n";
    if (tmode==Parallel) s+="play_mode = parallel\n";
    //s+="max_columns = "+inttostr(max_columns)+"\n";
    s+="time_factor = "+floattostr(time_factor)+"\n";
    s+="position_factor = "+floattostr(position_factor)+"\n";
    s+="abs_factor = "+floattostr(abs_factor)+"\n";
    s+="qfluor_factor = "+floattostr(qfluor_factor)+"\n";
    return s;
}

void DynamicsFromFiles::clear() {
    if (data!=NULL) {
        delete[] data;
    }
    data=NULL;
    if (shift_x!=NULL) free(shift_x);
    if (shift_y!=NULL) free(shift_y);
    if (shift_z!=NULL) free(shift_z);
    shift_x=NULL;
    shift_y=NULL;
    shift_z=NULL;
    trajectory_files.clear();
    trajectory_count=0;
}
