function plot_various_bits()

close all

Nb=16;

[t, x, y, z] = read_from_file(Nb,'');
displayname=[num2str(Nb) ' bits'];
plot_read_data(t,x,y,z,displayname);

max(t)

[t, x, y, z] = read_from_file(Nb,'_1');
displayname=[num2str(Nb) ' bits'];
plot_read_data(t,x,y,z,displayname);
max(t)


[t, x, y, z] = read_from_file(Nb,'_2');
displayname=[num2str(Nb) ' bits'];
plot_read_data(t,x,y,z,displayname);
max(t)


[t, x, y, z] = read_from_file(Nb,'_3');
displayname=[num2str(Nb) ' bits'];
plot_read_data(t,x,y,z,displayname);
max(t)

end

