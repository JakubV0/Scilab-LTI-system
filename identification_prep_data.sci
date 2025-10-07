clc;
clear;

// User set values
start = 1;              // Index of first sample for identification
stop = -1;              // Index of last sample for identification, stop = -1 for full data
file_name = 'data.csv'; // Name of a csv table.

// Load data from csv file a transform them into vectors
data = csvRead(file_name);
t = data(2:$,1);
u = data(2:$,2);
y = data(2:$,3);

// Cut out valid data for identification
if stop < 0 then
    t = t(start:stop) - t(start);
    u = u(start:stop) - u(start);
    y = y(start:stop) - y(start);
    save('data.sav');
end

N = length(t);
n = (0:N-1)';

scf(1);
subplot(2,1,1);
plot(n,u,n,y);
subplot(2,1,2);
plot(t,u,t,y);

