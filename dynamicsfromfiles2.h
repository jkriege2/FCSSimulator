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


#ifndef DYNAMICSFROMFILES_H
#define DYNAMICSFROMFILES_H

#include "fluorophordynamics.h"
#include<string>
#include<map>
#include "tools.h"

/** \brief this class loads several datafiles, containing trajectories and uses these trajectories
 *         for the dynamics simulation
 *  \ingroup diff4_dynamic
 *
 * Note that most of the parameters you may set for a simulation (sim_timestep ...) are not
 * usable with this class, as they are determined by the contents of the files.
 *
 * The trajectories may be played in parallel or one after the other. It is also possible to
 * center the trajectories around a certain point in space, i.e. the center of mass of the
 * trajectory may be shifted.
 *
 * All loaded files must have at least 3 columns (the simulation time with
 * \b equidistant steps does not have to be in the file but is defined by sim_timestep, and
 * then the position of the fluorophor. In parallel mode all files
 * have to have the same number of lines!!!
 *
 * Take care that it is importand that you don't simulate more steps than there are steps in the files!
 *
 * The file may have additional columns:
 *   - col 1: time
 *   - col 2-4: x,y,z position
 *   - col 5-7: polarisation vector in global coordinate system (cartesian unit vector!)
 *   - col 8: quantum efficiency [0..1]
 *   - col 9: absorbtion cross section in [m^2]
 *   - col 10: quantum state [integer]
 * .
 * 
 * 
 * Input files may also be in a simple binary format:
 * \verbatim
      contents                                size [bytes]  format
    +--------------------------------------+
    | intsize=sizeof(int)                  |  2             unsigned integer
    +--------------------------------------+
    | number_of_records                    |  8             unsigned integer
    +--------------------------------------+
    | x[0]                                 |  intsize       signed integer
    | y[0]                                 |  intsize       signed integer
    | z[0]                                 |  intsize       signed integer
    +--------------------------------------+
    | x[1]                                 |  intsize       signed integer
    | y[1]                                 |  intsize       signed integer
    | z[1]                                 |  intsize       signed integer
    +--------------------------------------+
    | ...                                  |  
    | ...                                  |  
    +--------------------------------------+
    | x[number_of_records-1]               |  intsize       signed integer
    | y[number_of_records-1]               |  intsize       signed integer
    | z[number_of_records-1]               |  intsize       signed integer
    +--------------------------------------+
    \endverbatim
 * This method is used to read file file \c dynfile.file_format=binary is set in the configuration file!
 */
class DynamicsFromFiles2 : public FluorophorDynamics
{
    public:
        /** Default constructor */
        DynamicsFromFiles2(FluorophorManager* fluorophors, std::string object_name);
        /** Default destructor */
        virtual ~DynamicsFromFiles2();

        /** \brief initialize the simulation environment (random walker positions ... */
        virtual void init();

        /** \brief output the next step from the trajectory
         *
         * \param boundary_check when \c boundary_check==false the walkers won't be reset when they reach
         *                       the border of the simulational box
         */
        virtual void propagate(bool boundary_check=true);

        /** \brief report the state of the object, as human-readable text */
        virtual std::string report();

        virtual double estimate_runtime();
        
        enum FileFormats {
            ffCSV=0,
            ffBinary=1
        };
        
        static std::string fileFormatToString(FileFormats format) {
            if (format==ffBinary) return "binary";
            return "csv";
        }
        
        static FileFormats stringToFileFormat(std::string format) {
            std::string f=tolower(format);
            if (f=="b" || f=="bin" || f=="binary" || f=="1") return ffBinary;
            return ffCSV;
        }
    protected:
        /** \brief read configuration from INI file */
        virtual void read_config_internal(jkINIParser2& parser);
	
	/** \brief set all walkers to non-existent and set end of trajectory true */
	void setAllDone();
        
        /** \brief read a single line or record from a trajectory file */
        virtual std::vector<double> dataFileReadLine(FILE* f);
        /** \brief open a trajectory file for reading */
        virtual FILE* dataFileOpen(const std::string& filename);
        /** \brief close an opened trajectory file */
        virtual void dataFileClose(FILE* file);
        /** \brief count the number of lines or records in a trajectory file */
        virtual unsigned long long dataFileCountLines(const std::string& filename);



        /** \brief these file handlers are used to read the data from the files */
        std::vector<FILE*> file;
        /** \brief conatins the linecount of each file */
        std::vector<unsigned long long> linecount;
        /** \brief contains the column count of each file  */
        std::vector<int> columncount;
        /** \brief coordinate shift (x-axis) of the trajectories */
        double* shift_x;
        /** \brief coordinate shift (y-axis) of the trajectories */
        double* shift_y;
        /** \brief coordinate shift (z-axis) of the trajectories */
        double* shift_z;
        /** \brief number of loaded trajectories */
        int trajectory_count;
        /** \brief release all memory allocated by this class and closes all files */
        virtual void clear();
        /** \brief the files to be loaded */
        std::vector<std::string> trajectory_files;
        /** \brief the possible simulation modes */
        enum TrajectoryMode {
            Parallel,
            Sequential,
	    SequentialParallel
        };
        /** \brief the possible shift modes */
        enum ShiftMode {
            Mean,
            HalfTime,
            RandomDisplacedMean,
            EndEndDistanceCenter,
            EndEndDistanceCenterRandom,
            NthTime
        };
        /** \brief range minimum for random mean displacement in x-direction */
        double randomdisplace_x_min;
        /** \brief range maximum for random mean displacement in x-direction */
        double randomdisplace_x_max;
        /** \brief range minimum for random mean displacement in y-direction */
        double randomdisplace_y_min;
        /** \brief range maximum for random mean displacement in y-direction */
        double randomdisplace_y_max;
        /** \brief range minimum for random mean displacement in z-direction */
        double randomdisplace_z_min;
        /** \brief range maximum for random mean displacement in z-direction */
        double randomdisplace_z_max;
        /** \brief how to shift the trajectories */
        ShiftMode shiftmode;
        /** \brief the mode in which to use several trajectory files */
        TrajectoryMode tmode;
        /** \brief indicates whether to shift the trajectories' COM to 0 or not */
        bool shift_trajectories;
        /** \brief counts the current file to output */
        int file_counter;
        /** \brief counts the line in the current file to output */
        int line_counter;
        /** \brief separator char in CSV files */
        char separator_char;
        /** \brief comment char in CSV files */
        char comment_char;
        /** \brief scaling factor of time column to seconds x[units of time_column]*factor=x[seconds] */
        //double time_factor;
        /** \brief scaling factor of position column to microns x[units of position_column]*factor=x[microns] */
        double position_factor;
        /** \brief scaling factor of absorbtion_corsssection column to meters^2 x[units of abs_column]*factor=x[m^2] */
        double abs_factor;
        /**  \brief scaling factor of fluorescence quantum efficiency column to [0..1] x[units of abs_column]*factor=x[m^2] */
        double qfluor_factor;
        /** \brief how many columns to read */
        int max_columns;
        /** \brief how many lines to read at most (or -1 if all)*/
        int max_lines;
        /** \brief maximum number of files to use, even if wildcard search finds more (-1 for all) )*/
        int max_files;
        /** \brief start numbered input files here */
        int num_start;
        /** \brief stop numbered input files here */
        int num_stop;
        /** \brief filename template */
        std::string filenames;
        /** \brief format of the input files */
        FileFormats fileformat;

        /** \brief column index of time column in file (default: 0)*/
        int col_time;
        /** \brief column index of x-position column in file (default: 1)*/
        int col_posx;
        /** \brief column index of y-position column in file (default: 2)*/
        int col_posy;
        /** \brief column index of z-position column in file (default: 3)*/
        int col_posz;
        /** \brief column index of absorption cross section column in file (default: 8)*/
        int col_abs;
        /** \brief column index of fluorescence quantum efficiency column in file (default: 7)*/
        int col_qfluor;
        /** \brief column index of qm_state column in file (default: 9)*/
        int col_qmstate;
        /** \brief column index of x-component of polarisation column in file (default: 4)*/
        int col_px;
        /** \brief column index of y-component of polarisation column in file (default: 5)*/
        int col_py;
        /** \brief column index of z-component of polarisation column in file (default: 6)*/
        int col_pz;
        /** \brief column index of x-component of polarisation column 1 in file (default: 4)*/
        int col_px1;
        /** \brief column index of y-component of polarisation column 1 in file (default: 5)*/
        int col_py1;
        /** \brief column index of z-component of polarisation column 1 in file (default: 6)*/
        int col_pz1;
        /** \brief if you use a linear combination of polarisation vectors, this compination is \f$ f\cdot\vec{p}+(1-f)\cdot\vec{p}_1 \f$
                   where \f$ f \f$ is this parameter. */
        double p_fraction;
        
        bool rotate_trajectories;
        std::vector<double> alpha1,alpha2,alpha3;
        int repeat_files;
	
	/** \brief number of trajectories to run in parallel in SequentialParallel mode */
	int parallelTrajectories;


        /** \brief how long did it take to load all files */
        double timing_loadall;
        /** \brief how long did it take to load 1 file (in average) */
        double timing_load1;
        
        struct BinaryFileInfo {
            int intSize;
            uint64_t records;
            uint16_t entries;
            
            BinaryFileInfo() {
                intSize=sizeof(int); 
                records=0;
                entries=3;
            }
        };
        std::map<FILE*, BinaryFileInfo> binFileInfo;
        void rotateTrajectory(int i, double alpha1, double alpha2, double alpha3);

};

#endif // DYNAMICSFROMFILES_H
