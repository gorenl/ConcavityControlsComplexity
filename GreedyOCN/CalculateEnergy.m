function E = CalculateEnergy(drca,theta)
%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

E = sum(drca(:,4).^(1-theta).*drca(:,5));