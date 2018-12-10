using DelimitedFiles
using ImplicitEquations, Plots

target_data = readdlm("data/target-data.dat")
target_mc = readdlm("data/target-mc.dat")

x_data = target_data[:,8]
y_data = target_data[:,9]
z_data = target_data[:,1]
x_mc = target_mc[:,8]
y_mc = target_mc[:,9]
z_mc = target_mc[:,1]

Rd = 2
Rmc = 2
d = y_data - y_mc
sizeData = size(y_data)
intersection_Volume = 0
total_Volume = 0
slice_diff = zeros(25)
slice_cum = zeros(25)
zred = zeros(25)
y12 = fill(1.2,sizeData)


for i in 1:25
    global intersection_Volume
    global total_Volume
    global slice_diff
    global slice_cum
    global zred
    local iV, tV, hY
    tV = Rd^2*pi*abs(z_data[i]-z_data[i+1])
    if abs((d[i]^2+Rd^2-Rmc^2)/(2*abs(d[i])*Rd))<=1 && abs((d[i]^2+Rmc^2-Rd^2)/(2*abs(d[i])*Rmc))<=1 && (-abs(d[i])+Rd+Rmc)*(abs(d[i])+Rd-Rmc)*(abs(d[i])-Rd+Rmc)*(abs(d[i])+Rd+Rmc)>0
        iV = (Rd^2*acos((d[i]^2+Rd^2-Rmc^2)/(2*abs(d[i])*Rd))+Rmc^2*acos((d[i]^2+Rmc^2-Rd^2)/(2*abs(d[i])*Rmc))-(1/2)*sqrt((-abs(d[i])+Rd+Rmc)*(abs(d[i])+Rd-Rmc)*(abs(d[i])-Rd+Rmc)*(abs(d[i])+Rd+Rmc)))*abs(z_data[i]-z_data[i+1])
    else
        iV = tV
    end
    if (total_Volume-intersection_Volume != 0)
        hD = Rd/2+abs(y_data[i]-1.2)
        hMC = Rmc/2+abs(y_mc[i]-1.2)
        iV -= Rmc*acos((Rmc-hMC)/Rmc)-(Rmc-hMC)*sqrt(2*Rmc*hMC-hMC^2)
        tV -= Rd*acos((Rd-hD)/Rd)-(Rd-hD)*sqrt(2*Rd*hD-hD^2)
    end
    intersection_Volume += iV
    total_Volume += tV
    slice_diff[i] = (tV-iV)/tV
    slice_cum[i] = (total_Volume-intersection_Volume)/total_Volume
    zred[i] = z_data[i]
end

# p1 = plot(slice_diff)
# p2 = plot(slice_cum)
# plot(p1,p2,lw=3,title="Lost volume in intersection (%)")
plot(zred,slice_diff, lw=3,
                      xlabel = "z",
                      ylabel = "%",
                      label="Residual volume in intersection (%)")
plot!(zred,slice_cum, lw=3,
                      xlabel = "z",
                      ylabel = "%",
                      label="Cumulated residual volume in intersection (%)")
savefig("myplot.png")
plot(z_data,y_data, ribbon = (Rd,Rd),
                    fillalpha = 0.3,
                    xlabel = "z",
                    ylabel = "y",
                    label="Data Target")
plot!(z_mc,y_mc, ribbon = (Rmc,Rmc),
                 fillalpha = 0.3,
                 xlabel = "z",
                 ylabel = "y",
                 label="MC Target")
plot!(z_mc,y12, fillalpha = 0.3,
                 xlabel = "z",
                 ylabel = "y",
                 label="Cut y=1.2")
savefig("myplot2.png")

ygif = fill(1.2,2)
xgif = [-3,3]

anim = @animate for i=1:26
    f(x,y) = x^2 + (y-y_data[i])^2 - Rd^2
    g(x,y) = x^2 + (y-y_mc[i])^2 - Rmc^2
    r = ((Le(f,0)) & (Ge(g,0)))
    t = string("Residual volume in intersection at z =", z_data[i])
    plot(r, xlabel = "x",
            xlims = (-3,3),
            xticks = -3:0.5:3,
            ylabel = "y",
            ylims = (-3,1.2),
            yticks = -3:0.5:3,
            title=t)
    plot!(xgif, ygif, label="Cut y=1.2")
end
gif(anim, "mygif.gif", fps = 3)

println((total_Volume-intersection_Volume)/total_Volume)
