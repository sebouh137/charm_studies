import ROOT
f = ROOT.TFile("/work/clas12/spaul/charm_hists/hadd.root")
h = f.Get("Kp_pim_2d")

xaxis = h.GetXaxis()
D0=1.869
dM=0.02
ixmin = xaxis.FindBin(D0-dM)
ixmax = xaxis.FindBin(D0+dM)

yaxis = h.GetYaxis()
lambdaC=2.286
dM=0.02
iymin = xaxis.FindBin(lambdaC-dM)
iymax = xaxis.FindBin(lambdaC+dM)



print("n events = ", h.Integral(ixmin,ixmax,iymin, iymax))
window = (yaxis.GetBinUpEdge(iymax)-yaxis.GetBinLowEdge(iymin))*(xaxis.GetBinUpEdge(ixmax)-xaxis.GetBinLowEdge(ixmin))
window*=1000*1000 #convert to MeV
print("events / (inv mass window) / (miss mass window) =", h.Integral(ixmin,ixmax,iymin, iymax)/window, "/MeV^2")
