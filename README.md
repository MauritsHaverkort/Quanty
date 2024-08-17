Reporistry of Quanty. If you intend to work on part of the code please discuss this first with Maurits, he coordinates the different parts and tasks. Part of the codes are experimental and used to test ideas of students. The branch master contains the released versions. The branch develop contains the main development cases and is transformed to Master before a new release.

Quanty is intended to be compilable on all linux systems found in the wild. This includes 10 year old clusters at some insitute of a univeristy that runs fine and never has been updated. We thus can not assume that there are libraries available. We link to blas and lapack but should not assume other external libraries. (i.e. no boost, hdf5 or fft library linking. For hdf5 we should think of inlcuding things in the source code.

Quanty was originally written in C. Based upon many requests I transfered the code to C++. After a few years of experience I might redo this. When students implement reasonable but inefficient code in C I can convert and optimize this. C++ allows the creation of very inefficient code and I have a hard time optimizing this. So unless you work on a really seperated part of the code that does not relate to the core functions I would highly recomment to program in C. (even when the file ends as .cpp)

The code compiles on windows, linux and mac using either intel, gnu, Clang, ... compilers. We follow the C++11 standard to make sure we stay compatible with different systems.

In order to make the io functions work there is a specific part of the code for windows, linux or mac. You specify this by using
-D LUA_USE_LINUX 
-D LUA_USE_MACOSX

for linux or mac operating systems on the compile command line.

You need to compile the file quanty/LuaInterface/main.cpp this will link all headers and other files needed

Possible compile commands for the intel compiler on linux are:

icc quanty/LuaInterface/main.cpp -D LUA_USE_LINUX -O3 -ipo -static -qopenmp -I"${MKLROOT}/include" -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl -o Quanty


i.e. we link to openmp, blas, lapack and the standard libraries 

