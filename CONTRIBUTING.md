# Contributing to Quanty

Thank you for wanting to work on Quanty. If you intend to work on part of the code, please discuss this with Maurits on the forum at [www.quanty.org](https://www.quanty.org) or by [e-mail](mailto:m.w.haverkort@thphys.uni-heidelberg.de) first; he coordinates the different parts and tasks.

## Branches

- The `master` branch is for stable releases and can be found at [git.quanty.org/Quanty/release](https://git.quanty.org/Quanty/release).
- The `develop` branch holds the current version with updates. It should be stable and free of bugs, i.e. only finished projects should be merged into `develop`. Historically this has not always been true. The `develop` branch can be found at [git.quanty.org/Quanty/develop](https://git.quanty.org/Quanty/develop).
- Subprojects, new functions and algorithm updates are developed in their own branch and then merged into `develop`.

## Experimental code

Part of the code is experimental and used to test new ideas. To make sharing and testing easier, we allow these parts to be merged into `develop` and `master`, in directories named after the main developer. The corresponding Quanty functions are then called `Name.****` (you can find these for Yi, Martin, …). Note that these functions are "unsafe". We do not guarantee anything and discourage their use. They are intended to be used only in collaboration with the person who wrote them. Once they are useful, tested and complete, we will include them in the main code base.

## Backward compatibility

Scripts that users wrote years ago must keep running. Please do not change the names of Lua functions or options, their default values or the format of their output. Add a new option or a new function instead.

## Tests

The test scripts are in the folder `TestScripts`. Every push is tested by the GitLab CI, which compiles Quanty with the Intel compiler on Linux and runs every test script on its own. A test passes when Quanty exits with exit code 0.

When you add a test script, register it in the matching `TEST_NAME` list in `.gitlab-ci.yml`. Otherwise the CI never runs it.

Tests run on Linux, while you may develop on another system. Do not use `math.random` in a test: its sequence differs between operating systems. Derive thresholds from the data the script builds, not from numbers read off a run on your own computer.

## C vs. C++

Quanty was originally written in C. Following many requests, I moved the code to C++. After a few years of experience I seriously thought about changing the code back to pure C. In the end I decided against this, but this does require some restraint from the people who develop code.

First the problem: when students implement reasonable but inefficient code in C, I can convert and optimize it. C++ allows the creation of very inefficient code, and I have a hard time optimizing that.

So unless you work on a truly separate part of the code that does not relate to the core functions, the code should use C and not C++, even when the file ends in `.cpp`. In practice:

- use arrays, `malloc` and `calloc` instead of vectors, maps and strings;
- use structures instead of classes.

I am aware that the use of pointers can lead to memory leaks and to bugs from accessing the wrong memory. It also requires some basic understanding of how computers work at the hardware level. But if I want to explain numerical optimizations, this is a good thing. How can one discuss the trade-offs of reusing memory versus a less efficient use of memory, or the difference between storing variables and recalculating values, if one does not have control over what is allocated in the first place?

Furthermore, even in high-level languages like Python one needs to understand the concept of pointers. How else do you explain the difference between

```python
a = 1
b = a
a = 3
print(b) # output 1
```

and

```python
a = [1, 2, 3]
b = a
a[1] = 5
print(b) # output [1, 5, 3]
```

Out-of-bounds reads and writes are bugs, but not security problems: a Quanty script can already do anything the user can; see [SECURITY.md](SECURITY.md).

## Code style

### Indentation

We prefer a common indentation style throughout the code. Opening and closing braces have their own line, and we indent by two spaces:

```c
if(n>0)
{
  for(i=0;i<n;i++)
  {
    a[i]=0.0;
  }
}
```

### Small special functions

We prefer commonly known C and Lua functions over special functions written for things that could have been done in one or two lines. Everybody knows the standard functions; nobody knows the specially written ones.

### Object files

Although Quanty is getting big, there is currently no reason to split the code into parts and compile them separately. First of all, if two functions have the same name, they should probably be the same function, so there is no need for two functions with the same name in different namespaces. I know of a code with two implementations of spherical harmonics, with different phase conventions, in different namespaces. Not having namespaces prevents this, and that is a good thing.

Furthermore, the compiler sees every function and can optimize. As long as compilation times are reasonable and memory usage is OK, we will stick with the current setup. (Actually, Quanty compiles quite fast compared to other modular codes I know.)

## Code layout

- The main functions of Quanty are in the folder `Core`. Here you also find `StructDef.cpp`, which defines the structures we use for data storage in Quanty. The C functions are unsafe and do no type checking (i.e. junk in, junk out).
- The Lua interface is in the folder `LuaInterface`. The file `main.cpp` links the code together with includes (yes, we compile one big file).
- Most of the Lua interface code is in the file `ssp.cpp`. The bottom of this file registers the Quanty functions in the Lua registry.
- The test scripts are in the folder `TestScripts`.

## Building

See [BUILDING.md](BUILDING.md).
