function [E_vec,Prob_vec,drca,n,i,lower_energy,median_delL_vec,median_delChi_vec]...
    = simulatedAnnealingLandscape(len,wid,res,stopT,theta,T0,alpha,DEM,drca,mapping)
%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

global t_area movie reality_check
%Calculate energy, del_L and del_chi for the initial conditions = input
%landscape
E = CalculateEnergy(drca,theta);
[L,chi] = CalcFlowLenAndChi(drca,theta);
[medianDelL,~,~,medianDelChi,~,~] = ...
    DeltaLDeltaChi(DEM,drca,mapping,wid,L,chi);

E_vec = E;
median_delL_vec = medianDelL;
median_delChi_vec = medianDelChi;


no_outlet_vec_ind = drca(:,1)~=drca(:,2);
no_outlet_vec = drca(no_outlet_vec_ind,1);

T_vec = T0;

Prob_vec = 1;
i = 1;

lower_energy = 0;

if movie
   max_chi = max(chi);
   max_L = max(L);
   mv_name1 = strcat('mv',num2str(theta),'_chi_highres.mp4');
   vidfile1 = VideoWriter(mv_name1,'MPEG-4');
   open(vidfile1);
   im = plotLandscape(DEM,drca,mapping,t_area,11);
   hold on
   [xd,yd] = ind2coord(DEM,drca(:,1));
   scatter(xd,yd,5,chi);
   colorbar;
   frame = getframe(im);
   writeVideo(vidfile1, frame);
   close 11
    
   mv_name2 = strcat('mv',num2str(theta),'_L_highres.mp4');
   vidfile2 = VideoWriter(mv_name2,'MPEG-4');
   open(vidfile2);
   im = plotLandscape(DEM,drca,mapping,t_area,12);
   hold on
   [xd,yd] = ind2coord(DEM,drca(:,1));
   scatter(xd,yd,5,L);
   colorbar;
   frame = getframe(im);
   writeVideo(vidfile2, frame);
   close 12 

   mv_name3 = strcat('mv',num2str(theta),'_White_highres.mp4');
   vidfile3 = VideoWriter(mv_name3,'MPEG-4');
   open(vidfile3);
   im = plotLandscape(DEM,drca,mapping,t_area,13);
   frame = getframe(im);
   writeVideo(vidfile2, frame);
   close 13 

end

%commented out critetion for number of iterations based on temperature 
%while T0*alpha^i > stopT

while i <= stopT % number of global iterations
    i = i+1 % just to follow code progression
    for j = 1:ceil(len*wid/res) % internal iterations. More nodes == more iterations
        %attempt a new network with a single flip
        drca_neighbour = generateNeighbour(drca,mapping,no_outlet_vec,wid,res);
        %calculate the energy of the new network
        E_neighbour = CalculateEnergy(drca_neighbour,theta);
        
        % Acceptance probability for greedy approach. Returns 1 if energy
        % is lower
        accept_prob = acceptanceProbabilityDecreaseOnly(E,E_neighbour,T0*alpha^i);
        
        % Acceptance probability for true SimulatedAnnealing
        %accept_prob = acceptanceProbability(E,E_neighbour,T0*alpha^i)
       
        if E_neighbour < E
            lower_energy = lower_energy+1;
        end
        
        % Could be used for post-processing
        Prob_vec = [Prob_vec accept_prob];
        T_vec = [T_vec T0*alpha^i];

        %If accepting the new network, updating parameters and calculating
        %its del_L and del_chi
        if accept_prob >= rand
            drca = drca_neighbour;
            E = E_neighbour;
            [L,chi] = CalcFlowLenAndChi(drca,theta);
            [medianDelL,~,~,medianDelChi,~,~] = ...
                 DeltaLDeltaChi(DEM,drca,mapping,wid,L,chi);
            %E_vec(i) = E;
            
            if reality_check
                plotLandscape(DEM,drca,mapping,t_area);
            end
        end
        E_vec = [E_vec E];
        median_delL_vec = [median_delL_vec medianDelL];
        median_delChi_vec = [median_delChi_vec medianDelChi];

        if movie && mod(j,len*wid/res/5) == 1
            [L,chi] = CalcFlowLenAndChi(drca,theta);
            [xd,yd] = ind2coord(DEM,drca(:,1));

            %chi background
            im = plotLandscape(DEM,drca,mapping,t_area,11);
            hold on
            scatter(xd,yd,5,chi);
            colorbar
            frame = getframe(im);
            writeVideo(vidfile1, frame);
            close 11

            %L background
            im = plotLandscape(DEM,drca,mapping,t_area,12);
            hold on
            scatter(xd,yd,5,L);
            colorbar
            frame = getframe(im);
            writeVideo(vidfile2, frame);
            close 12

            %No background
            im = plotLandscape(DEM,drca,mapping,t_area,13);
            frame = getframe(im);
            writeVideo(vidfile3, frame);
            close 13

        end
    end
    
end
n = length(E_vec);
if movie
    close(vidfile1)
    close(vidfile2)
    close(vidfile3)
end
