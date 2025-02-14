================== Some general remarks                            =======================

Repository of Quanty. If you intend to work on part of the code please discuss this first with Maurits, he coordinates the different parts and tasks. 

The Master branch is for stable releases

The Develop branch is for the current version with updates. The develop branch should be stable and bug free, i.e. only finished projects should be pushed to develop. Historically not always true. 

Sub projects, new functions and algoritm updates are developed in their own branch and then merged to develop.

Part of the codes are experimental and used to test ideas of students. In order to facilitate sharing and testing we allow these to be merged into Develop and Master in directories with the name of the student. The functions in Quanty then are called Name.**** (you can find these for Yi, Martin, ...). Note that the use of these functions is "unsafe". We do not garante anything and discourage the use of such functions. They are intended to use only in collaboration with the person who wrote them. If they are useful, tested and complete we will include them into the main code base.

================== Compilation                                     =======================

Quanty is intended to be compilable on all linux windows and Macintosh systems found in the wild. This includes 10 year old clusters at some insitute of a univeristy that runs fine and never has been updated. We can and do not assume that there are libraries available. Quanty links to blas and lapack but should not assume other external libraries. (i.e. no boost, hdf5 or fft library linking.)

The code compiles on windows, linux and mac using either intel, gnu, Clang, ... compilers. We follow the C++11 standard to make sure we stay compatible with different systems.

In order to make the io functions work there is a specific part of the code for windows, linux or mac. You specify this by using
-D LUA_USE_LINUX 
-D LUA_USE_MACOSX
for linux or mac operating systems on the compile command line.

You need to compile the file quanty/LuaInterface/main.cpp this will link all headers and other files needed

Possible compile commands for the intel compiler on linux are:

icc quanty/LuaInterface/main.cpp -D LUA_USE_LINUX -O3 -ipo -static -qopenmp -I"${MKLROOT}/include" -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl -o Quanty

A possible command for GNU g++ is:

g++ LuaInterface/main.cpp -D LUA_USE_LINUX -DBLASunderscore -O3 -fopenmp -flto -o Quanty -lblas -llapack

Note the use of the BLASunderscore macro here, that Quanty exposes to use different conventions of blas and lapack libraries.



i.e. we link to openmp, blas, lapack and the standard libraries 


================== C v.s. C++                                      =======================

Quanty was originally written in C. Based upon many requests I transfered the code to C++. After a few years of experience I seriously tought about changing the code back to pure C. In the end I decided againts this, but this does require some restraints on the people who develop code.

First the problem: When students implement reasonable but inefficient code in C I can convert and optimize this. C++ allows the creation of very inefficient code and I have a hard time optimizing this. 

So unless you work on a really seperated part of the code that does not relate to the core functions the code should use C and not C++, even when the file ends with cpp. 

In practice:
-- use arrays, malloc and calloc instead of vectors, maps and strings.
-- use structures instead of classes

I am aware that the use of pointers can lead to memory loss and bugs because of accessing wrong memory. It also requires some basic understanding of how computers work on a hardware level. But if I want to explain numerical optimisations this is a good thing. How can one discuss the trade-offs off reusing memory v.s. less efficient use of memory or the difference between storing variables and recalculating values if one does not have controle over what is allocated in the first place. 

Furthermore even in high level languages like python one needs to understand the concept of pointers. How else do you explain the difference between 
a=1
b=a
a=3
print(b) -- output 1
 and
a={1,2,3}
b=a
a[2] = 5
print(b) -- output {1,5,3}

Last but not least, out of bound read and writes can give rise to security issues on computer servers. This should not be a big concern for our code. Quanty is not intended to run the databases of a bank or hospital. On your personal desktop Quanty is safe. On a central server you should give Quanty access to the full node.

==========================================================================================
================== Code style                                      =======================
==========================================================================================


================== Indentation                                     =======================

We prefer to have a common style for indentation throughout the code. Opening and closing brackets have their own line and we indent by two spaces.


================== Small special functions                         =======================

We do prefer the use of commonly known C and lua functions over the use of special functions defined for things that could have been done in one or two lines. The standard functions are known, the special written functions not.


================== object files                                   =======================

Although Quanty is getting big there is currently no reason to split the code into parts and separately compile them. First of all, if two functions have the same name, they probably should be the same function so there is no need to have two functions with the same name in different name spaces. I know a code with two implementations of spherical harmonics with different phase definitions in different name spaces. Not having namespaces prevents this and that is a good thing. 

Furthermore the compiler sees every function and can optimize. As long as compilation times are reasonable and memory usage ok we role with the current setup. (actually compilation of Quanty is quite fast compared to other modular codes I know)

==========================================================================================
================== Code setup                                      =======================
==========================================================================================

The main functions of Quanty can be found in the folder Core. Here you also find the fileStructDef.cpp. This file defines the structures we use for data storage in Quanty.

The C functions are unsafe and do no type checking (i.e. junk in junk out).

You find the lua based interface in the folder LuaInterface. The file main.cpp links the code together with includes (yes we have one big file we compile).

The main Lua interface based code can be found in the file ssp.cpp. The botom of this file registers the Quanty functions in the lua registry.


==========================================================================================

All in all, I hope we made something useful, have fun calculating and coding, be respectful and we appreciate feedback and bug-reports.
Best wishes,
Maurits 

