using DelimitedFiles
using Plots
using LaTeXStrings

z = [.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85]
xexp = [.007,.015,.025,.035,.05,.08,.12,.16,.29]
zFFred = [3,8,13,18,23,28,33,38,43,48,53,61]

# Mult1 = readdlm("data/Mult2016Hadron.txt")
# Mult2 = readdlm("data/Mult2016HadronNORICH.txt")
Mult1 = readdlm("data/Mult2016Hadron.txt")
Mult2 = readdlm("data/Mult2016HadronMTOT.txt")

Multp1 = zeros((9,12))
Multm1 = zeros((9,12))
Multp2 = zeros((9,12))
Multm2 = zeros((9,12))

for i in 1:9
    for j in 1:12
        global Multm1, Multm2
        Multm1[i,j] = Mult1[j+24*(i-1),4]
        Multm2[i,j] = Mult2[j+24*(i-1),4]
    end
    for j in 1:12
        global Multp1, Multp2
        Multp1[i,j] = Mult1[j+24*(i-1)+12,4]
        Multp2[i,j] = Mult2[j+24*(i-1)+12,4]
    end
end

zmid = zeros(12)
for i in 1:12
    zmid[i] = (z[i+1]+z[i])/2
end

Ms1 = zeros(9)
Ms2 = zeros(9)
Mp1 = zeros(12)
Mp2 = zeros(12)
Mr1 = zeros(9)
Mr2 = zeros(9)
Mm1 = zeros(12)
Mm2 = zeros(12)
Mrp1 = zeros(9)
Mrp2 = zeros(9)
Mrm1 = zeros(9)
Mrm2 = zeros(9)

for i in 1:9
    global Mr1, Mr2
    for j in 1:12
        global Ms1, Ms2, Mp1, Mp2, Mm1, Mm2, Mrp1, Mrp2, Mrm1, Mrm2
        Ms1[i] += (Multp1[i,j]+Multm1[i,j])*(z[j+1]-z[j])
        Ms2[i] += (Multp2[i,j]+Multm2[i,j])*(z[j+1]-z[j])
        Mrp1[i] += Multp1[i,j]*(z[j+1]-z[j])
        Mrp2[i] += Multp2[i,j]*(z[j+1]-z[j])
        Mrm1[i] += Multm1[i,j]*(z[j+1]-z[j])
        Mrm2[i] += Multm2[i,j]*(z[j+1]-z[j])
        Mp1[j] = Multp1[i,j]
        Mm1[j] = Multm1[i,j]
        Mp2[j] = Multp2[i,j]
        Mm2[j] = Multm2[i,j]
    end
    t = string("Pion Multiplicities at x =", xexp[i])
    plot(zmid, Mp1, lw=3, xlims = (0.,1),
                           ylims = (-0.1,3),
                           linecolor = :red,
                           xlabel = "z",
                           ylabel = L"M^{\pi}",
                           label = L"\pi^+_{Reference}")
    plot!(zmid, Mm1, lw=3, xlims = (0.,1),
                            ylims = (-0.1,3),
                            linecolor = :blue,
                            xlabel = "z",
                            ylabel = L"M^{\pi}",
                            label = L"\pi^-_{Reference}")
    plot!(zmid, Mp2, lw=3, xlims = (0.,1),
                             ylims = (-0.1,3),
                             linecolor = :red,
                             linestyle = :dash,
                             xlabel = "z",
                             ylabel = L"M^{\pi}",
                             label = L"\pi^+_{Test}")
    plot!(zmid, Mm2, lw=3, xlims = (0.,1),
                             ylims = (-0.1,3),
                             linecolor = :blue,
                             linestyle = :dash,
                             xlabel = "z",
                             ylabel = L"M^{\pi}",
                             label = L"\pi^-_{Test}",
                             title=t)
    savefig(string("plots/MultiplicitiesTest",i,".png"))
    Mr1[i] = (Mrp1[i]/Mrm1[i])
    Mr2[i] = (Mrp2[i]/Mrm2[i])
end

plot(xexp,Ms1, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           ylims = (0.7,1.2),
           xlabel = "x",
           ylabel = L"\int M^{\pi^+}+M^{\pi^-} dz",
           label = "Reference")
plot!(xexp,Ms2, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (0.7,1.2),
            xlabel = "x",
            ylabel = L"\int M^{\pi^+}+M^{\pi^-} dz",
            label = "Test")
savefig("plots/MultiplicitiesSumTest.png")

plot(xexp,Mr1, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           xlabel = "x",
           ylabel = L"\frac{\int M^{\pi^+} dz}{\int M^{\pi^-} dz}",
           label = "Reference")
plot!(xexp,Mr2, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            xlabel = "x",
            ylabel = L"\frac{\int M^{\pi^+} dz}{\int M^{\pi^-} dz}",
            label = "Test")
savefig("plots/MultiplicitiesRatioTest.png")
