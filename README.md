# clas12-tof-efficiency
Efficiency determination tools for the CLAS12 Time-of-Flight system, analysis macros in clas12root

CTOF_eff.C - Raw track efficiency determination

CTOF_effUnified.C - Unified version with all missing particle cases, and single track, in one code, for ease of processing

All versions run in the same way, with two available methods;

CTOF_eff() runs with hard-coded input file names, for testing and development

CTOF_eff(TString infile, TString outfile, const std::string databasefile) needs input and output file names as an argument, and a prepared root file with rcdb information for the corresponding run(s) studied (all require full or relative path to file)

Plotting Macros and JLab farm submission scripts also provided
