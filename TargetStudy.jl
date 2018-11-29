using DelimitedFiles

target_data = readdlm("data/target-274508-274901.dat")
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
Intersection_Volume = 0
Total_Volume = 0

for i = 1:29
    Intersection_Volume += (2*2^2*acos(d[i]/(2*R))-(d[i]/2)*sqrt(4*R^2-d[i]^2))*abs(z_data[i]-z_data[i+1])
    Total_Volume += 4*pi*abs(z_data[i]-z_data[i+1])
end

show(Intersection_Volume/Total_Volume)
