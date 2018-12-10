using DelimitedFiles
using Plots
using LaTeXStrings

PDF = readdlm("data/MSTW2008lo68cl.txt")
FFPion = readdlm("data/evolvedCOMPASSFFPion.txt")

x = PDF[:,1]
Q2 = PDF[:,2]
u = PDF[:,3]
ub = PDF[:,4]
d = PDF[:,5]
db = PDF[:,6]
s = PDF[:,7]
sb = PDF[:,8]

z = [.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85]
Dfavpi = FFPion[:,]
Dfavpi = FFPion[:,]

for i in 1:9
    for i in 1:12
        Mppiplus = ((4*u[i]+db[i])*Dfavpi[i][k]+(4*ub[i]+d[i]+s[i]+sb[i])*Dunfpi[i][k])/(4*(u[i]+ub[i])+(d[i]+db[i])+(s[i]+sb[i]))
        Mppiminus = ((4*ub[i]+d[i])*Dfavpi[i][k]+(4*u[i]+db[i]+s[i]+sb[i])*Dunfpi[i][k])/(4*(u[i]+ub[i])+(d[i]+db[i])+(s[i]+sb[i]))
        Mdpiplus = ((4*u[i]*(u[i]+d[i])+db[i]*(ub[i]+db[i]))*Dfavpi[i][k]+(u[i]+d[i]+ub[i]+db[i]+2*(s[i]+sb[i]))*Dunfpi[i][k])/(5*(u[i]+ub[i]+d[i]+db[i])+2*(s[i]+sb[i]))
        Mdpiminus = ((4*ub[i]*(ub[i]+db[i])+d[i]*(u[i]+d[i]))*Dfavpi[i][k]+(ub[i]+db[i]+u[i]+d[i]+2*(s[i]+sb[i]))*Dunfpi[i][k])/(5*(u[i]+ub[i]+d[i]+db[i])+2*(s[i]+sb[i]))
        Mp[i] = (Mppiplus+Mppiminus)*(z[k+1]-z[k])
        Md[i] = (Mdpiplus+Mdpiminus)*(z[k+1]-z[k])
    end
end

plot(x,Mp, lw=3,
           xlabel = "x",
           ylabel = L"\int M^{\pi^+}+M^{\pi^-} dz",
           label = "Proton")
plot!(x,Md, lw=3,
            xlabel = "x",
            ylabel = L"\int M^{\pi^+}+M^{\pi^-} dz",
            label = "Deuteron")
savefig("MultiplicitiesProtDeut.png")
