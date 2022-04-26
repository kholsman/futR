# comilation Error Notes


# Install Matrix if not yet installed

# Compilation error may related to RTools version: https://cran.r-project.org/bin/windows/Rtools/
    # To use rtools, download the installer from CRAN:
        # On Windows 64-bit: rtools40-x86_64.exe (https://cran.r-project.org/bin/windows/Rtools/rtools40-x86_64.exe includes both i386 and x64 compilers). Permanent url: rtools40-x86_64.exe.
        # On Windows 32-bit: rtools40-i686.exe (i386 compilers only). Permanent url: rtools40-i686.exe.

# write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)

#Restart R, and verify that make can be found, which should show the path to your Rtools installation.

Sys.which("make")
## "C:\\rtools40\\usr\\bin\\make.exe"


#Now try to install an R package from source:
  
  install.packages("jsonlite", type = "source")






https://stackoverflow.com/questions/42756640/what-is-this-compilation-error-in-the-tmb-tutorial

I have solved my problem now, which came in two parts. The first was not having Rtools installed in Windows. This is mentioned here:
  
  Why can I not build and run models in Windows, after installing TMB from CRAN? TMB on Windows requires Rtools. The PATH environment variable should point to the Rtools 'make' and 'gcc', and no other instances of 'make' or 'gcc'.

The second problem was that after installing Rtools I still could not find the compile, despite verifying that I could use it from the command line following these instructions:
  
  In some cases this PATH might be C:\RBuildTools... Further, we can check if g++ can be really called from R. For example, we can see the version of gcc in R as follows.

> system('g++ -v')
Using built-in specs.
COLLECT_GCC=c:\Rtools\GCC-46~1.3\bin\G__~1.EXE
COLLECT_LTO_WRAPPER=c:/rtools/gcc-46~1.3/bin/../libexec/gcc/i686-w64-mingw32/4.6.3/lto-wrapper.exe
Target: i686-w64-mingw32
Configured with: /data/gannet/ripley/Sources/mingw-test3/src/gcc/configure --host=i686-w64-mingw32 --build=x86_64-linux-gnu --target=i686-w64-mingw32 --with-sysroot=/data/gannet/ripley/Sources/mingw-test3/mingw32mingw32/mingw32 --prefix=/data/gannet/ripley/Sources/mingw-test3/mingw32mingw32/mingw32 --with-gmp=/data/gannet/ripley/Sources/mingw-test3/mingw32mingw32/prereq_install --with-mpfr=/data/gannet/ripley/Sources/mingw-test3/mingw32mingw32/prereq_install --with-mpc=/data/gannet/ripley/Sources/mingw-test3/mingw32mingw32/prereq_install --disable-shared --enable-static --enable-targets=all --enable-languages=c,c++,fortran --enable-libgomp --enable-sjlj-exceptions --enable-fully-dynamic-string --disable-nls --disable-werror --enable-checking=release --disable-win32-registry --disable-rpath --disable-werror CFLAGS='-O2 -mtune=core2 -fomit-frame-pointer' LDFLAGS=
  Thread model: win32
gcc version 4.6.3 20111208 (prerelease) (GCC)

> system('where make')
c:\Rtools\bin\make.exe
The error I received was:
  
  > compile('tutorial.cpp')
c:/Rtools/mingw_64/bin/g++  
  ...
c:/Rtools/mingw_64/bin/g++: not found
This suggested that the path was being hard coded somewhere. Following the instructions for using alternative gcc compilers, I discovered that the path was hard coded in my Makeconf file located at: "C:\Program Files\R\R-3.3.2\etc\x64\Makeconf"

I then commented out the old line and replaced it with the right path:
  
  # BINPREF ?= c:/Rtools/mingw_64/bin/
  BINPREF ?= C:/RBuildTools/3.4/mingw_64/bin/
  After this the compilation works without a restart.

Share
Follow