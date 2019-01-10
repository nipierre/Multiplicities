using DelimitedFiles
using Plots
using LaTeXStrings

z = [.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85]
xexp = [.007,.015,.025,.035,.05,.08,.12,.16,.29]
zFFred = [3,8,13,18,23,28,33,38,43,48,53,61]

PDF = readdlm("data/MSTW2008lo68cl.txt")
FFPion = readdlm("data/FFPiondata0.txt")
FFKaon = readdlm("data/FFKaondata0.txt")
ExpMultPion = readdlm("data/Mult2016Pion.txt")
ExpMultPionD = readdlm("data/Mult2006Pion.txt")

x = PDF[:,1]
Q2 = PDF[:,2]
u = PDF[:,3]
ub = PDF[:,4]
d = PDF[:,5]
db = PDF[:,6]
s = PDF[:,7]
sb = PDF[:,8]

zFF = zeros((15,66))
Dfavpi = zeros((15,66))
Dunfpi = zeros((15,66))
DfavK = zeros((15,66))
DunfK = zeros((15,66))
DstrK = zeros((15,66))
Multpip = zeros((9,12))
Multpim = zeros((9,12))
MultpipD = zeros((9,12))
MultpimD = zeros((9,12))

for i in 1:15
    for j in 1:66
        global zFF
        global Dfavpi
        global Dunfpi
        global DfavK
        global DunfK
        global DstrK
        zFF[i,j] = FFPion[j+(i-1)*66,1]
        Dfavpi[i,j] = FFPion[j+(i-1)*66,2]/zFF[i,j]
        Dunfpi[i,j] = FFPion[j+(i-1)*66,4]/zFF[i,j]
        DfavK[i,j] = FFKaon[j+(i-1)*66,2]/zFF[i,j]
        DunfK[i,j] = FFKaon[j+(i-1)*66,4]/zFF[i,j]
        DstrK[i,j] = FFKaon[j+(i-1)*66,4]/zFF[i,j]
    end
end

for i in 1:9
    for j in 1:12
        global Multpim
        Multpim[i,j] = ExpMultPion[j+24*(i-1),4]
        MultpimD[i,j] = ExpMultPionD[j+24*(i-1),4]
    end
    for j in 1:12
        global Multpip
        Multpip[i,j] = ExpMultPion[j+24*(i-1)+12,4]
        MultpipD[i,j] = ExpMultPionD[j+24*(i-1)+12,4]
    end
end

zmid = zeros(12)
for i in 1:12
    zmid[i] = (z[i+1]+z[i])/2
end

Mpisp = zeros(9)
Mpisd = zeros(9)
Mpispexp = zeros(9)
Mpisdexp = zeros(9)
Mpirp = zeros(9)
Mpird = zeros(9)
Mpirpexp = zeros(9)
Mpirdexp = zeros(9)
Mpirpp = zeros(9)
Mpirdp = zeros(9)
Mpirppexp = zeros(9)
Mpirdpexp = zeros(9)
Mpirpm = zeros(9)
Mpirdm = zeros(9)
Mpirpmexp = zeros(9)
Mpirdmexp = zeros(9)
MKsp = zeros(9)
MKsd = zeros(9)
MKrp = zeros(9)
MKrd = zeros(9)
MKrpp = zeros(9)
MKrdp = zeros(9)
MKrpm = zeros(9)
MKrdm = zeros(9)
Mpip = zeros(12)
Mpim = zeros(12)
Mpexp = zeros(12)
Mmexp = zeros(12)

for i in 1:9
    global Mpirp, Mpird, Mpirpexp, Mpirdexp
    global MKrp, MKrd
    global Mpip, Mpim
    global Mpexp, Mmexp
    for j in 1:12
        global Mpisp, Mpisd, Mpispexp, Mpisdexp, Mpirpp, Mpirdp, Mpirppexp, Mpirdpexp, Mpirpm, Mpirdm, Mpirpmexp, Mpirdmexp
        global MKsp, MKsd, MKrpp, MKrdp, MKrpm, MKrdm
        local Mppiplus, Mppiminus, Mdpiplus, Mdpiminus
        local MpKplus, MpKminus, MdKplus, MdKminus
        Mppiplus = ((4*u[i]+db[i])*Dfavpi[round(Int,Q2[i]),zFFred[j]]+(4*ub[i]+d[i]+s[i]+sb[i])*Dunfpi[round(Int,Q2[i]),zFFred[j]])/(4*(u[i]+ub[i])+(d[i]+db[i])+(s[i]+sb[i]))
        Mppiminus = ((4*ub[i]+d[i])*Dfavpi[round(Int,Q2[i]),zFFred[j]]+(4*u[i]+db[i]+s[i]+sb[i])*Dunfpi[round(Int,Q2[i]),zFFred[j]])/(4*(u[i]+ub[i])+(d[i]+db[i])+(s[i]+sb[i]))
        Mdpiplus = ((4*(u[i]+d[i])+(ub[i]+db[i]))*Dfavpi[round(Int,Q2[i]),zFFred[j]]+((u[i]+d[i])+4*(ub[i]+db[i])+2*(s[i]+sb[i]))*Dunfpi[round(Int,Q2[i]),zFFred[j]])/(5*(u[i]+ub[i]+d[i]+db[i])+2*(s[i]+sb[i]))
        Mdpiminus = ((4*(ub[i]+db[i])+(u[i]+d[i]))*Dfavpi[round(Int,Q2[i]),zFFred[j]]+((ub[i]+db[i])+4*(u[i]+d[i])+2*(s[i]+sb[i]))*Dunfpi[round(Int,Q2[i]),zFFred[j]])/(5*(u[i]+ub[i]+d[i]+db[i])+2*(s[i]+sb[i]))
        Mpisp[i] += (Mppiplus+Mppiminus)*(z[j+1]-z[j])
        Mpisd[i] += (Mdpiplus+Mdpiminus)*(z[j+1]-z[j])
        Mpispexp[i] += (Multpip[i,j]+Multpim[i,j])*(z[j+1]-z[j])
        Mpisdexp[i] += (MultpipD[i,j]+MultpimD[i,j])*(z[j+1]-z[j])
        Mpirpp[i] += Mppiplus*(z[j+1]-z[j])
        Mpirdp[i] += Mdpiplus*(z[j+1]-z[j])
        Mpirppexp[i] += Multpip[i,j]*(z[j+1]-z[j])
        Mpirdpexp[i] += MultpipD[i,j]*(z[j+1]-z[j])
        Mpirpm[i] += Mppiminus*(z[j+1]-z[j])
        Mpirdm[i] += Mdpiminus*(z[j+1]-z[j])
        Mpirpmexp[i] += Multpim[i,j]*(z[j+1]-z[j])
        Mpirdmexp[i] += MultpimD[i,j]*(z[j+1]-z[j])
        Mpip[j] = Mppiplus
        Mpim[j] = Mppiminus
        Mpexp[j] = Multpip[i,j]
        Mmexp[j] = Multpim[i,j]
        MpKplus = ((4*u[i])*DfavK[round(Int,Q2[i]),zFFred[j]]+(4*ub[i]+d[i]+s[i]+db[i])*DunfK[round(Int,Q2[i]),zFFred[j]]+sb[i]*DstrK[round(Int,Q2[i]),zFFred[j]])/(4*(u[i]+ub[i])+(d[i]+db[i])+(s[i]+sb[i]))
        MpKminus = ((4*ub[i])*Dfavpi[round(Int,Q2[i]),zFFred[j]]+(4*u[i]+db[i]+sb[i]+db[i])*DunfK[round(Int,Q2[i]),zFFred[j]]+s[i]*DstrK[round(Int,Q2[i]),zFFred[j]])/(4*(u[i]+ub[i])+(d[i]+db[i])+(s[i]+sb[i]))
        MdKplus = ((4*(u[i]+d[i]))*DfavK[round(Int,Q2[i]),zFFred[j]]+((u[i]+d[i])+5*(ub[i]+db[i])+2*(s[i]))*DunfK[round(Int,Q2[i]),zFFred[j]]+2*sb[i]*DstrK[round(Int,Q2[i]),zFFred[j]])/(5*(u[i]+ub[i]+d[i]+db[i])+2*(s[i]+sb[i]))
        MdKminus = ((4*(ub[i]+db[i]))*DfavK[round(Int,Q2[i]),zFFred[j]]+((ub[i]+db[i])+5*(u[i]+d[i])+2*(sb[i]))*DunfK[round(Int,Q2[i]),zFFred[j]]+2*s[i]*DstrK[round(Int,Q2[i]),zFFred[j]])/(5*(u[i]+ub[i]+d[i]+db[i])+2*(s[i]+sb[i]))
        MKsp[i] += (MpKplus+MpKminus)*(z[j+1]-z[j])
        MKsd[i] += (MdKplus+MdKminus)*(z[j+1]-z[j])
        MKrpp[i] += MpKplus*(z[j+1]-z[j])
        MKrdp[i] += MdKplus*(z[j+1]-z[j])
        MKrpm[i] += MpKminus*(z[j+1]-z[j])
        MKrdm[i] += MdKminus*(z[j+1]-z[j])
    end
    t = string("Pion Multiplicities at x =", xexp[i])
    plot(zmid, Mpip, lw=3, xlims = (0.,1),
                           ylims = (-0.1,3),
                           linecolor = :red,
                           xlabel = "z",
                           ylabel = L"M^{\pi}",
                           label = L"\pi^+_{Pred}")
    plot!(zmid, Mpim, lw=3, xlims = (0.,1),
                            ylims = (-0.1,3),
                            linecolor = :blue,
                            xlabel = "z",
                            ylabel = L"M^{\pi}",
                            label = L"\pi^-_{Pred}")
    plot!(zmid, Mpexp, lw=3, xlims = (0.,1),
                             ylims = (-0.1,3),
                             linecolor = :red,
                             linestyle = :dash,
                             xlabel = "z",
                             ylabel = L"M^{\pi}",
                             label = L"\pi^+_{Exp}")
    plot!(zmid, Mmexp, lw=3, xlims = (0.,1),
                             ylims = (-0.1,3),
                             linecolor = :blue,
                             linestyle = :dash,
                             xlabel = "z",
                             ylabel = L"M^{\pi}",
                             label = L"\pi^-_{Exp}",
                             title=t)
    savefig(string("plots/MultiplicitiesProtPi",i,".png"))
    Mpirp[i] = (Mpirpp[i]/Mpirpm[i])
    Mpirpexp[i] = (Mpirppexp[i]/Mpirpmexp[i])
    Mpirdexp[i] = (Mpirdpexp[i]/Mpirdmexp[i])
    Mpird[i] = (Mpirdp[i]/Mpirdm[i])
    MKrp[i] = (MKrpp[i]/MKrpm[i])
    MKrd[i] = (MKrdp[i]/MKrdm[i])
end

plot(x,Mpisp, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           ylims = (0.4,0.9),
           xlabel = "x",
           ylabel = L"\int M^{\pi^+}+M^{\pi^-} dz",
           label = "Proton Pred.")
plot!(x,Mpisd, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (0.4,0.9),
            xlabel = "x",
            ylabel = L"\int M^{\pi^+}+M^{\pi^-} dz",
            label = "Deuteron Pred.")
plot!(xexp,Mpispexp, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (0.4,0.9),
            xlabel = "x",
            ylabel = L"\int M^{\pi^+}+M^{\pi^-} dz",
            label = "Proton Exp.")
plot!(xexp,Mpisdexp, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (0.4,0.9),
            xlabel = "x",
            ylabel = L"\int M^{\pi^+}+M^{\pi^-} dz",
            label = "Deuteron Exp.")
savefig("plots/MultiplicitiesProtDeutSumPi.png")

plot(x,Mpirp, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           xlabel = "x",
           ylabel = L"\frac{\int M^{\pi^+} dz}{\int M^{\pi^-} dz}",
           label = "Proton Pred.")
plot!(x,Mpird, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            xlabel = "x",
            ylabel = L"\frac{\int M^{\pi^+} dz}{\int M^{\pi^-} dz}",
            label = "Deuteron Pred.")
plot!(xexp,Mpirpexp, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            xlabel = "x",
            ylabel = L"\frac{\int M^{\pi^+} dz}{\int M^{\pi^-} dz}",
            label = "Proton Exp.")
plot!(xexp,Mpirdexp, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            xlabel = "x",
            ylabel = L"\frac{\int M^{\pi^+} dz}{\int M^{\pi^-} dz}",
            label = "Deuteron Exp.")
savefig("plots/MultiplicitiesProtDeutRatioPi.png")

plot(x,MKsp, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           xlabel = "x",
           ylabel = L"\int M^{K^+}+M^{K^-} dz",
           label = "Proton")
plot!(x,MKsd, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            xlabel = "x",
            ylabel = L"\int M^{K^+}+M^{K^-} dz",
            label = "Deuteron")
savefig("plots/MultiplicitiesProtDeutSumK.png")

plot(x,MKrp, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           xlabel = "x",
           ylabel = L"\frac{\int M^{K^+} dz}{\int M^{K^-} dz}",
           label = "Proton")
plot!(x,MKrd, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            xlabel = "x",
            ylabel = L"\frac{\int M^{K^+} dz}{\int M^{K^-} dz}",
            label = "Deuteron")
savefig("plots/MultiplicitiesProtDeutRatioK.png")
