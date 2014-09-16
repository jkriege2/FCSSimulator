#! /bin/bash

function print_result {
	if [ $2 == -1 ]; then
		echo -e "building '$1': \033[0;33m NOT SELECTED \033[0m" 
	fi
	if [ $2 == -2 ]; then
		echo -e "building '$1': \033[0;31m ERROR IN CONFIGURE \033[0m" 
	fi
	if [ $2 == -3 ]; then
		echo -e "building '$1': \033[0;31m ERROR IN MAKE STEP \033[0m" 
	fi
	if [ $2 == -4 ]; then
		echo -e "building '$1': \033[0;31m ERROR IN INSTALL STEP \033[0m" 
	fi
	if [ $2 == -5 ]; then
		echo -e "building '$1': \033[0;33m INSTALL SCRIPT DOES NOT SUPPORT THIS ON YOUR PLATFORM \033[0m" 
	fi
	if [ $2 == 0 ]; then
		echo -e "building '$1': \033[1;32m OK \033[0m" 
	fi
}



echo -e "========================================================================\n"\
"= QuickFit 3 dependency build script                                   =\n"\
"========================================================================\n\n"\
"This script will try to build all the dependency libraries \n"\
"needed for QuickFit. All libraries will be built in static linking mode\n"\
"and stay inside their subdirectory of this extlibs directory, so they\n"\
"don't interfere with your system-libraries. The .pri-files in this di-\n"\
"rectory include the libraries in a way that preferences these libraries \n"\
"here!\n\nFirst we need to set some basics for the builds:"\

read -p "Do you want to build keep the build directories (y/n)? " -n 1 KEEP_BUILD_DIR
echo -e  "\n"
read -p "How many parallel builds do you want to use in make (1/2/3/...)? " -n 1  MAKE_PARALLEL_BUILDS
echo -e  "\n"
read -p "Do you want to use more agressive optimizations for the built libraries (y/n)? " -n 1  MAKE_AGRESSIVEOPTIMIZATIONS
echo -e  "\n"
read -p "Do you want to optimize libraries for your local machine? (y/n)? " -n 1  MAKE_COMPILEFORLOCAL
echo -e  "\n"

#sh ../output/get_bit_depth.sh
MORECFLAGS=" -mtune=generic -msse -msse2 -mmmx -m3dnow -mfpmath=sse "
if [ $MAKE_COMPILEFORLOCAL == "y" ] ; then
	if [ $MAKE_AGRESSIVEOPTIMIZATIONS == "y" ] ; then
		MORECFLAGS=" -mtune=native -msse -msse2 -mmmx -m3dnow -mfpmath=sse -ftree-vectorize -ftree-vectorizer-verbose=1"
	else
		MORECFLAGS=" -mtune=native -msse -msse2 -mmmx -m3dnow -mfpmath=sse "
	fi
else
	if [ $MAKE_AGRESSIVEOPTIMIZATIONS == "y" ] ; then
		MORECFLAGS=" -mtune=generic -msse -msse2 -mmmx -m3dnow -mfpmath=sse -ftree-vectorize -ftree-vectorizer-verbose=1"
	fi
fi

CURRENTDIR=${PWD}

ISMSYS=`uname -o`
echo $ISMSYS
if [ "$ISMSYS" != "${string/Msys/}" ] ; then
	echo -e "building in MSys environment on Windows!\n\n"
fi




libgslOK=-1
read -p "Do you want to build 'gsl' (y/n)? " -n 1 INSTALL_ANSWER
echo -e  "\n"
if [ $INSTALL_ANSWER == "y" ] ; then
	echo -e  "\n------------------------------------------------------------------------\n"\
	"-- BUILDING: gsl                                                       --\n"\
	"------------------------------------------------------------------------\n\n"\

	cd gsl
	mkdir build
	tar xvf gsl-1.16.tar.gz -C ./build/
	cd build/gsl-1.16
	./configure --enable-static --disable-shared --prefix=${CURRENTDIR}/gsl   CFLAGS="-fPIC ${MORECFLAGS}" CPPFLAGS="-fPIC ${MORECFLAGS}"
	libOK=$?
	if [ $libOK -eq 0 ] ; then
		make -j${MAKE_PARALLEL_BUILDS}
		
		libOK=$?
		if [ $libOK -eq 0 ] ; then		
			make -j${MAKE_PARALLEL_BUILDS} install
			libOK=$?
			if [ $libOK -ne 0 ] ; then		
				libOK=-4
			fi
		else
			libOK=-3
		fi
	else
	    libOK=-2
	fi
	

	cd ../../
	if [ $KEEP_BUILD_DIR == "n" ] ; then
		rm -rf build
	fi
	cd ${CURRENTDIR}
	
	libgslOK=$libOK

fi



echo -e  "\n------------------------------------------------------------------------\n"\
"-- BUILD RESULTS                                                       --\n"\
"------------------------------------------------------------------------\n\n"\

print_result "gsl" $libgslOK


