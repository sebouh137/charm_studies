import ROOT
f = ROOT.TFile("/work/clas12/spaul/charm_hists/hadd.root")
h = f.Get("Kp_pim_2d_zoom")

xaxis = h.GetXaxis()
D0=1.869
dM=0.05
ixmin = xaxis.FindBin(D0-dM/2)
ixmax = xaxis.FindBin(D0+dM/2)

yaxis = h.GetYaxis()
lambdaC=2.286
dM=0.100
iymin = yaxis.FindBin(lambdaC-dM/2)
iymax = yaxis.FindBin(lambdaC+dM/2)

print(yaxis.GetBinUpEdge(iymax),yaxis.GetBinLowEdge(iymin),xaxis.GetBinUpEdge(ixmax),xaxis.GetBinLowEdge(ixmin))

N=h.Integral(ixmin,ixmax,iymin, iymax)
print("n events = ", N)
window = (yaxis.GetBinUpEdge(iymax)-yaxis.GetBinLowEdge(iymin))*(xaxis.GetBinUpEdge(ixmax)-xaxis.GetBinLowEdge(ixmin))
window*=1000*1000 #convert to MeV
print("events / (inv mass window) / (miss mass window) =", N/window, "/MeV^2")

#now for neighborhood
xaxis = h.GetXaxis()
D0=1.869
dM=0.05
ixmin = xaxis.FindBin(D0-dM*3/2)
ixmax = xaxis.FindBin(D0+dM*3/2)

yaxis = h.GetYaxis()
lambdaC=2.286
dM=0.100
iymin = yaxis.FindBin(lambdaC-dM*3/2)
iymax = yaxis.FindBin(lambdaC+dM*3/2)

window2 = (yaxis.GetBinUpEdge(iymax)-yaxis.GetBinLowEdge(iymin))*(xaxis.GetBinUpEdge(ixmax)-xaxis.GetBinLowEdge(ixmin))
window2*=1000*1000 #convert to MeV                                                                                                          
N2=h.Integral(ixmin,ixmax,iymin, iymax)

print("events / (inv mass window) / (miss mass window) in neighborhood =", (N2-N)/(window2-window), "/MeV^2")
