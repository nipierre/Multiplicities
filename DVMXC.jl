using DelimitedFiles
using Plots
using LaTeXStrings

z = [.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85]
xexp = [.007,.015,.025,.035,.05,.08,.12,.16,.29]

DIS = readdlm("data/DISFrac.dat")
Had = readdlm("data/HadFrac.dat")
DISa = readdlm("data/DISfractions.txt")
Hada = readdlm("data/HADfractions.txt")

DIS_SIDIS = zeros((9,6))
DIS_rho = zeros((9,6))
DIS_phi = zeros((9,6))
DIS_SIDISa = zeros((9,6))
DIS_rhoa = zeros((9,6))
DIS_phia = zeros((9,6))
Had_SIDIS = zeros((9,6,12))
Had_rho = zeros((9,6,12))
Had_phi = zeros((9,6,12))
Had_SIDISa = zeros((9,6,12))
Had_rhoa = zeros((9,6,12))
Had_phia = zeros((9,6,12))

for i in 1:54
    DIS_SIDIS[Int.(DIS[i,1]),Int.(DIS[i,2])] = DIS[i,3]
    DIS_rho[Int.(DIS[i,1]),Int.(DIS[i,2])] = DIS[i,4]
    DIS_phi[Int.(DIS[i,1]),Int.(DIS[i,2])] = DIS[i,5]
end

for i in 1:45
    DIS_SIDISa[Int.(DISa[i,1]+1),Int.(DISa[i,2]+1)] = DISa[i,3]
    DIS_rhoa[Int.(DISa[i,1]+1),Int.(DISa[i,2]+1)] = DISa[i,4]
    DIS_phia[Int.(DISa[i,1]+1),Int.(DISa[i,2]+1)] = DISa[i,5]
end

for i in 1:648
    Had_SIDIS[Int.(Had[i,1]),Int.(Had[i,2]),Int.(Had[i,3])] = Had[i,4]
    Had_rho[Int.(Had[i,1]),Int.(Had[i,2]),Int.(Had[i,3])] = Had[i,5]
    Had_phi[Int.(Had[i,1]),Int.(Had[i,2]),Int.(Had[i,3])] = Had[i,6]
end

for i in 1:540
    Had_SIDISa[Int.(Hada[i,1]+1),Int.(Hada[i,2]+1),Int.(Hada[i,3]+1)] = Hada[i,4]
    Had_rhoa[Int.(Hada[i,1]+1),Int.(Hada[i,2]+1),Int.(Hada[i,3]+1)] = Hada[i,5]
    Had_phia[Int.(Hada[i,1]+1),Int.(Hada[i,2]+1),Int.(Hada[i,3]+1)] = Hada[i,6]
end

zmid = zeros(12)
for i in 1:12
    zmid[i] = (z[i+1]+z[i])/2
end

DIS_s = zeros((9,6))
DIS_r = zeros((9,6))
DIS_p = zeros((9,6))
Had_s = zeros(12,45)
Had_r = zeros(12,45)
Had_p = zeros(12,45)

for i in 1:9
    for j in 1:5
        global DIS_s, DIS_r, DIS_p
        global Had_s, Had_r, Had_p
        DIS_s[i,j] = (DIS_SIDIS[i,j] > 0 && DIS_SIDISa[i,j] > 0) ? DIS_SIDIS[i,j]/DIS_SIDISa[i,j] : 0
        DIS_r[i,j] = (DIS_rho[i,j] > 0 && DIS_rhoa[i,j] > 0) ? DIS_rho[i,j]/DIS_rhoa[i,j] : 0
        DIS_p[i,j] = (DIS_phi[i,j] > 0 && DIS_phia[i,j] > 0) ? DIS_phi[i,j]/DIS_phia[i,j] : 0
        for k in 1:12
            Had_s[k,i+(j-1)*9] = (Had_SIDIS[i,j,k] > 0 && Had_SIDISa[i,j,k] > 0) ? Had_SIDIS[i,j,k]/Had_SIDISa[i,j,k] : 0
            Had_r[k,i+(j-1)*9] = (Had_rho[i,j,k] > 0 && Had_rhoa[i,j,k] > 0) ? Had_rho[i,j,k]/Had_rhoa[i,j,k] : 0
            Had_p[k,i+(j-1)*9] = (Had_phi[i,j,k] > 0 && Had_phia[i,j,k] > 0) ? Had_phi[i,j,k]/Had_phia[i,j,k] : 0
        end
    end
end

plot(xexp,DIS_s,layout=(2,3),seriestype=:scatter,
                                 xscale = :log10,
                                 xlims = (0.01,1),
                                 size = (1600,900),
                                 dpi = 150,
                                 ylims = (0.,2.),
                                 # ylims = (0.75,1.25),
                                 legend=false)
savefig("plots/DVMXC_DIS_SIDIS.png")

plot(xexp,DIS_r,layout=(2,3),seriestype=:scatter,
                                 xscale = :log10,
                                 xlims = (0.01,1),
                                 size = (1600,900),
                                 dpi = 150,
                                 ylims = (0.,2.),
                                 # ylims = (0.75,1.25),
                                 legend=false)
savefig("plots/DVMXC_DIS_rho.png")

plot(xexp,DIS_p,layout=(2,3),seriestype=:scatter,
                                 xscale = :log10,
                                 xlims = (0.01,1),
                                 size = (1600,900),
                                 dpi = 150,
                                 ylims = (0.,2.),
                                 # ylims = (0.75,1.25),
                                 legend=false)
savefig("plots/DVMXC_DIS_phi.png")

plot(zmid,Had_s,layout=(5,9),seriestype=:scatter,
                                size = (3200,1800),
                                dpi = 300,
                                ylims = (0.,2.),
                                # ylims = (0.75,1.25),
                                legend=false)
savefig("plots/DVMXC_Had_SIDIS.png")

plot(zmid,Had_r,layout=(5,9),seriestype=:scatter,
                                size = (3200,1800),
                                dpi = 300,
                                ylims = (0.,2.),
                                # ylims = (0.75,1.25),
                                legend=false)
savefig("plots/DVMXC_Had_rho.png")

plot(zmid,Had_p,layout=(5,9),seriestype=:scatter,
                                size = (3200,1800),
                                dpi = 300,
                                ylims = (0.,2.),
                                # ylims = (0.75,1.25),
                                legend=false)
savefig("plots/DVMXC_Had_phi.png")
