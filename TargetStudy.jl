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

R = 2
d = y_data - y_mc
sizeData = size(y_data)
intersection_Volume = 0
total_Volume = 0
slice_diff = zeros(25)
slice_cum = zeros(25)
zred = zeros(25)


for i in 1:25
    global intersection_Volume
    global total_Volume
    global slice_diff
    global slice_cum
    global zred
    local iV, tV
    iV = (2*2^2*acos(abs(d[i])/(2*R))-(abs(d[i])/2)*sqrt(4*R^2-d[i]^2))*abs(z_data[i]-z_data[i+1])
    tV = 4*pi*abs(z_data[i]-z_data[i+1])
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
plot(z_data,y_data, ribbon = (1,1),
                    fillalpha = 0.3,
                    xlabel = "z",
                    ylabel = "y",
                    label="Data Target")
plot!(z_mc,y_mc, ribbon = (1,1),
                 fillalpha = 0.3,
                 xlabel = "z",
                 ylabel = "y",
                 label="MC Target")
savefig("myplot2.png")

anim = @animate for i=1:26
    f(x,y) = x^2 + (y-y_data[i])^2 - 4
    g(x,y) = x^2 + (y-y_mc[i])^2 - 4
    r = (Le(f,0)) & (Ge(g,0))
    t = string("Residual volume in intersection at z =", z_data[i])
    plot(r, xlabel = "x",
            xlims = (-3,3),
            xticks = -3:0.5:3,
            ylabel = "y",
            ylims = (-3,3),
            yticks = -3:0.5:3,
            title=t)
end
gif(anim, "mygif.gif", fps = 3)

println((total_Volume-intersection_Volume)/total_Volume)
