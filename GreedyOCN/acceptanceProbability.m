function prob = acceptanceProbability(E,E_neighbour,Temp)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

if E_neighbour < E
    prob = 1;
else
    prob = exp(-(E_neighbour-E)/Temp);
end
