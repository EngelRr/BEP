clear
close all
clc

profile on;

memory = 2;
states = zeros(8,memory);

for  i = 1:1e8
    poep = size(states,2);
end

for i = 1:1e8
    poep = memory;
end

profile viewer