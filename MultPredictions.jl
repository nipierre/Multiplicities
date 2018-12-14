using DelimitedFiles
using Plots
using LaTeXStrings

z = [.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85]
zFFred = [3,8,13,18,23,28,33,38,43,48,53,61]

PDF = readdlm("data/MSTW2008lo68cl.txt")
FFPion = readdlm("data/FFPiondata0.txt")
FFKaon = readdlm("data/FFKaondata0.txt")

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

Mpisp = zeros(9)
Mpisd = zeros(9)
Mpirp = zeros(9)
Mpird = zeros(9)
MKsp = zeros(9)
MKsd = zeros(9)
MKrp = zeros(9)
MKrd = zeros(9)

for i in 1:9
    for j in 1:12
        global Mpisp, Mpisd, Mpirp, Mpird
        global MKsp, MKsd, MKrp, MKrd
        local Mppiplus, Mppiminus, Mdpiplus, Mdpiminus
        local MpKplus, MpKminus, MdKplus, MdKminus
        Mppiplus = ((4*u[i]+db[i])*Dfavpi[round(Int,Q2[i]),zFFred[j]]+(4*ub[i]+d[i]+s[i]+sb[i])*Dunfpi[round(Int,Q2[i]),zFFred[j]])/(4*(u[i]+ub[i])+(d[i]+db[i])+(s[i]+sb[i]))
        Mppiminus = ((4*ub[i]+d[i])*Dfavpi[round(Int,Q2[i]),zFFred[j]]+(4*u[i]+db[i]+s[i]+sb[i])*Dunfpi[round(Int,Q2[i]),zFFred[j]])/(4*(u[i]+ub[i])+(d[i]+db[i])+(s[i]+sb[i]))
        Mdpiplus = ((4*(u[i]+d[i])+(ub[i]+db[i]))*Dfavpi[round(Int,Q2[i]),zFFred[j]]+((u[i]+d[i])+4*(ub[i]+db[i])+2*(s[i]+sb[i]))*Dunfpi[round(Int,Q2[i]),zFFred[j]])/(5*(u[i]+ub[i]+d[i]+db[i])+2*(s[i]+sb[i]))
        Mdpiminus = ((4*(ub[i]+db[i])+(u[i]+d[i]))*Dfavpi[round(Int,Q2[i]),zFFred[j]]+((ub[i]+db[i])+4*(u[i]+d[i])+2*(s[i]+sb[i]))*Dunfpi[round(Int,Q2[i]),zFFred[j]])/(5*(u[i]+ub[i]+d[i]+db[i])+2*(s[i]+sb[i]))
        Mpisp[i] += (Mppiplus+Mppiminus)*(z[j+1]-z[j])
        Mpisd[i] += (Mdpiplus+Mdpiminus)*(z[j+1]-z[j])
        Mpirp[i] += (Mppiplus/Mppiminus)*(z[j+1]-z[j])
        Mpird[i] += (Mdpiplus/Mdpiminus)*(z[j+1]-z[j])
        MpKplus = ((4*u[i])*DfavK[round(Int,Q2[i]),zFFred[j]]+(4*ub[i]+d[i]+s[i]+db[i])*DunfK[round(Int,Q2[i]),zFFred[j]]+sb[i]*DstrK[round(Int,Q2[i]),zFFred[j]])/(4*(u[i]+ub[i])+(d[i]+db[i])+(s[i]+sb[i]))
        MpKminus = ((4*ub[i])*Dfavpi[round(Int,Q2[i]),zFFred[j]]+(4*u[i]+db[i]+sb[i]+db[i])*DunfK[round(Int,Q2[i]),zFFred[j]]+s[i]*DstrK[round(Int,Q2[i]),zFFred[j]])/(4*(u[i]+ub[i])+(d[i]+db[i])+(s[i]+sb[i]))
        MdKplus = ((4*(u[i]+d[i]))*DfavK[round(Int,Q2[i]),zFFred[j]]+((u[i]+d[i])+5*(ub[i]+db[i])+2*(s[i]))*DunfK[round(Int,Q2[i]),zFFred[j]]+2*sb[i]*DstrK[round(Int,Q2[i]),zFFred[j]])/(5*(u[i]+ub[i]+d[i]+db[i])+2*(s[i]+sb[i]))
        MdKminus = ((4*(ub[i]+db[i]))*DfavK[round(Int,Q2[i]),zFFred[j]]+((ub[i]+db[i])+5*(u[i]+d[i])+2*(sb[i]))*DunfK[round(Int,Q2[i]),zFFred[j]]+2*s[i]*DstrK[round(Int,Q2[i]),zFFred[j]])/(5*(u[i]+ub[i]+d[i]+db[i])+2*(s[i]+sb[i]))
        MKsp[i] += (MpKplus+MpKminus)*(z[j+1]-z[j])
        MKsd[i] += (MdKplus+MdKminus)*(z[j+1]-z[j])
        MKrp[i] += (MpKplus/MpKminus)*(z[j+1]-z[j])
        MKrd[i] += (MdKplus/MdKminus)*(z[j+1]-z[j])
    end
end

plot(x,Mpisp, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           ylims = (0.5,0.95),
           xlabel = "x",
           ylabel = L"\int M^{\pi^+}+M^{\pi^-} dz",
           label = "Proton")
plot!(x,Mpisd, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (0.5,0.95),
            xlabel = "x",
            ylabel = L"\int M^{\pi^+}+M^{\pi^-} dz",
            label = "Deuteron")
savefig("MultiplicitiesProtDeutSumPi.png")

plot(x,Mpirp, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           xlabel = "x",
           ylabel = L"\int M^{\pi^+}/M^{\pi^-} dz",
           label = "Proton")
plot!(x,Mpird, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            xlabel = "x",
            ylabel = L"\int M^{\pi^+}/M^{\pi^-} dz",
            label = "Deuteron")
savefig("MultiplicitiesProtDeutRatioPi.png")

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
savefig("MultiplicitiesProtDeutSumK.png")

plot(x,MKrp, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           xlabel = "x",
           ylabel = L"\int M^{K^+}/M^{K^-} dz",
           label = "Proton")
plot!(x,MKrd, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            xlabel = "x",
            ylabel = L"\int M^{K^+}/M^{K^-} dz",
            label = "Deuteron")
savefig("MultiplicitiesProtDeutRatioK.png")
