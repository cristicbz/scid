#!/usr/local/bin/rdmd --shebang

/** Build system for SciD

    This is an experimental build system for SciD, using a regular D
    script instead of a makefile.

    Usage:
    To build the library and generate header (.di) files, run
    ---
    rdmd build
    ---
    To build only the library file, run
    ---
    rdmd build lib
    ---
    To only generate header files, run
    ---
    rdmd build headers
    ---
    To make the documentation, run
    ---
    rdmd build html
    ---
*/
import std.algorithm, std.array, std.exception, std.file, std.path, std.process,
    std.stdio, std.string, std.zip, std.getopt;



/** Various build directories.

    NOTE:
    Running "build clean" will recursively and mercilessly
    delete these directories.  Make sure you don't use the current
    directory, the root directory or any such thing.  Preferably,
    just leave them the way they are.
*/
version (Posix)
{
    immutable libDir    = "generated";          // Location of lib file
    immutable headerDir = "generated/headers";  // Location of .di files
    immutable htmlDir   = "generated/html";     // Location of .html files
}
version (Windows)
{
    immutable libDir    = r"generated";
    immutable headerDir = r"generated\headers";
    immutable htmlDir   = r"generated\html";
}


/** The name of the library. */
immutable libName   = "scid";

/** The top-level directory of the source files. */
immutable srcDir    = "scid";


int main(string[] args)
in { assert (args.length > 0); }
body
{
    string libBLAS, libLAPACK;
    string compilerName = "dmd";
    getopt(args, config.passThrough,
        "blas", &libBLAS,
        "lapack", &libLAPACK,
        "compiler", &compilerName
    );

    if(args.length < 2) {
        printUsage();
        return -1;
    }

    try
    {
        immutable flags = (args.length > 2) ? (join(args[2..$], " ")) : "";

        switch(args[1]) {
            case "lib":
            case "libs":
                buildLib(compilerName, flags, libBLAS, libLAPACK);
                break;
            case "demo":
                buildDemo(compilerName, flags, libBLAS, libLAPACK);
                break;
            case "headers":
                buildHeaders();
                break;
            case "html":
                buildHTML();
                break;
            case "clean":
                buildClean();
                break;
            default:
                printUsage();
        }

        return 0;
    }
    catch (Exception e)
    {
        stderr.writeln(e.msg);
        return 1;
    }
}



/** Build the library file. */
void buildLib(string compiler, string flags, string libBLAS, string libLAPACK)
{
    ensureDir(libDir);
    auto sources = getSources();

    version (Posix)     immutable libFile = "lib"~libName~".a";
    version (Windows)   immutable libFile = libName~".lib";

    auto buildCmd = compiler ~ " -lib " ~ flags;
    if(libBLAS.length) buildCmd ~= " " ~ libBLAS;
    if(libLAPACK.length) buildCmd ~= " " ~ libLAPACK;
    buildCmd ~= " " ~ std.string.join(sources, " ")
        ~" -od"~libDir~" -of"~libFile;
    writeln(buildCmd);
    enforce(system(buildCmd) == 0, "Error building library");
}

/** Build the demo file. */
void buildDemo(string compiler, string flags, string libBLAS, string libLAPACK)
{
    ensureDir(libDir);
    auto sources = getSources();

    auto buildCmd = compiler ~ " -unittest -version=demo " ~ flags;
    if(libBLAS.length) buildCmd ~= " " ~ libBLAS;
    if(libLAPACK.length) buildCmd ~= " " ~ libLAPACK;
    buildCmd ~= " " ~ std.string.join(sources, " ")
        ~" -ofdemo";
    writeln(buildCmd);
    enforce(system(buildCmd) == 0, "Error building demo");
}

/** Generate header files. */
void buildHeaders()
{
    ensureDir(headerDir);
    auto sources = getSources();
    foreach (s; sources)
    {
        immutable d = headerDir; // std.path.join(headerDir, dirname(s));
        ensureDir(d);

        immutable diName = baseName(s, ".d")~".di";
        immutable cmd = "dmd "~s~" -c -o- -H -Hf"~d~diName;
        writeln(cmd);
        enforce(system(cmd) == 0, "Error making header file: "~diName);
    }
}



/** Build documentation. */
void buildHTML()
{
    auto sources = getSources();
    sort(sources);

    ensureDir(htmlDir);
    unzip("candydoc.zip", htmlDir);

    string[] htmlFiles;
    string moduleList = "MODULES =\n";
    foreach (s; sources)
    {
        version (Posix)     auto slash = "/";
        version (Windows)   auto slash = "\\";
        htmlFiles ~= baseName(replace(s, slash, "_"),".d") ~ ".html";

        // Do not list the scid.internal.* modules in the
        // doc browser tree.
        if (std.string.indexOf(s, "internal") == -1)
            moduleList ~=
                "\t$(MODULE "
                ~baseName(replace(s, slash, "."), ".d")
                ~")\n";
    }

    immutable modulesDdoc = std.path.buildPath(htmlDir, "candydoc", "modules.ddoc");
    writeln("Writing "~modulesDdoc);
    std.file.write(modulesDdoc, moduleList);

    immutable candyDdoc = std.path.buildPath(htmlDir, "candydoc", "candy.ddoc");
    foreach (i; 0 .. sources.length)
    {
        immutable cmd =
            "dmd "~sources[i]~" "~candyDdoc~" "~modulesDdoc
            ~" -c -o- -D -Dd"~htmlDir~" -Df"~htmlFiles[i];
        writeln(cmd);
        enforce(system(cmd) == 0, "Error making HTML file: "~htmlFiles[i]);
    }
}



/** Remove build directories. */
void buildClean()
{
    void rm(string path)
    {
        if (!exists(path)) return;
        writeln("Removing: ", path);
        if (isDir(path)) rmdirRecurse(path);
        else std.file.remove(path);
    }

    rm(libDir);
    rm(headerDir);
    rm(htmlDir);
    rm(__FILE__~".deps");   // Clean up after rdmd as well
}



// Various utility functions


string[] getSources()
{
    static string[] sources;
    if (sources == null)
    {
        foreach (string f; dirEntries(srcDir, SpanMode.depth))
            if (isFile(f) && extension(f) == ".d") sources ~= f;
    }
    return sources;
}


void ensureDir(string dir)
{
    if (exists(dir)) enforce(isDir(dir), "Not a directory: "~dir);
    else mkdirRecurse(dir);
}


void unzip(string zipFile, string toDir)
{
    writeln("Unzipping "~zipFile~" to "~toDir);
    auto zip = new ZipArchive(std.file.read(zipFile));
    foreach (member; zip.directory)
    {
        if (member.name[$-1] == '/') continue;  // Skip directory names

        immutable f = std.path.buildPath(toDir, member.name);
        ensureDir(dirName(f));
        std.file.write(f, zip.expand(member));
    }
}

void printUsage() {
    enum usage =
"Build tool for SciD.
Usage:  build target [options...] [compiler flags...]

Targets:
    lib     = The SciD library
    demo    = The included demo/test program.
    headers = .di files for SciD.
    html    = The SciD documentation.
    clean   = Remove all files from previous build attempts.

Options:
    --compiler = The name of the compiler command (default = dmd).  If
                 another compiler (e.g. GDC or LDC) is used, then the
                 adapter script (gdmd or ldmd2) should be provided to this
                 command so that a DMD-style argument syntax is accepted.

    --blas     = A BLAS library file to be included in the generated library.
                 If none is provided, then a BLAS library will need to be
                 explicitly provided when building a SciD application.

    --lapack   = A LAPACK library file to be included in the generated library.
                 If none is provided, then a LAPACK library will need to be
                 explicitly provided when building a SciD application.

All other command line arguments are interpreted as compiler flags.
Example:

./build lib --compiler gdmd --blas /usr/lib/libblas.a
  --lapack /usr/lib/liblapack.a -m64 -O -inline -release";

    writeln(usage);
}
