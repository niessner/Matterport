GAPS Users -

This directory contains all code for the GAPS software library.
There are several subdirectories:

    pkgs - source and include files for all packages (software libraries).
    apps - source files for several application and example programs. 
    makefiles - unix-style make file definitions
    lib - archive library (.lib) files (created during compilation).
    bin - executable files (created during compilation).

If you are using linux or cygwin and have gcc and OpenGL development
libraries installed, or if you are using MAC OS X with the xcode
development environment, you should be able to compile all the code by
typing "make clean; make" in this directory.  If you are using Windows
Visual Studio, then you will have to build your own solution files
(the code should compile and run in Visual Studio). For other
development platforms, you should edit the shared compilation settings
in the makefiles/Makefiles.std to meet your needs.

To write a program that uses the GAPS pkgs, then you should include
"-I XXX/gaps/pkgs" in your compile flags (CFLAGS) and "-L
XXX/gaps/lib" in your link flags (LDFLAGS), where XXX is the directory
where you installed the gaps software.

The software is distributed under the MIT license (see LICENSE.txt)
and thus can be used for any purpose without warranty, any liability,
or any suport of any kind.

- Tom Funkhouser



