# clas12-tof-efficiency
Efficiency determination tools for the CLAS12 Time-of-Flight system, analysis macros in clas12root

CTOF_eff.C - Raw track efficiency determination

CTOF_effUnified.C - Attempt at unifying all missing particle cases, plus single track, into one code, for ease of processing

All versions run in the same way, with two available methods;

CTOF_eff() runs with a hard-coded input file name, as previous test/development versions

CTOF_eff(TString infile) takes input file name as an argument (requires full or relative path to file)

Plotting Macros also provided
