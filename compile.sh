#mkdir 
#g++ src/LcpSkim.cc -o LcpSkim -L/group/clas12/packages/hipo/1.3/lib -lhipo4 -llz4 -I/group/clas12/packages/hipo/1.3/hipo4 -I/group/clas12/users/rafopar/clas12AnaTools/include -L/group/clas12/users/rafopar/clas12AnaTools/lib -lclas12AnaTools `root-config --cflags --libs`
g++ src/LcpSkim.cc -o LcpSkim -L/group/clas12/packages/hipo/1.3/lib -lhipo4 -llz4 -I/group/clas12/packages/hipo/1.3/hipo4 
