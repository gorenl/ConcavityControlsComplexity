function drca_neighbour = generateNeighbour(drca,mapping,no_outlet_vec,wid,res)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

% function to generate a new network following a single valid (no loops) 
% random edge flip

%tic

drca_neighbour = drca;

%Defining a random donor that switches its receiver.
don_to_switch_ind = randi(length(no_outlet_vec),1);
don_to_switch = no_outlet_vec(don_to_switch_ind);
row_of_don = mapping(don_to_switch);


%remove drainage area of donor from its former receiver and all downstream
%nodes
old_rec = drca(row_of_don,2);
row_of_old_rec = mapping(old_rec);
drca_neighbour(row_of_old_rec,4) = drca_neighbour(row_of_old_rec,4)-...
    drca(row_of_don,4);
while drca(row_of_old_rec,1) ~= drca(row_of_old_rec,2) %while not an outlet
    %temp_don = old_rec;
    old_rec = drca(row_of_old_rec,2);
    row_of_old_rec = mapping(old_rec);
    if drca_neighbour(row_of_old_rec,4)<=drca(row_of_don,4)
        'stop'
    end
    drca_neighbour(row_of_old_rec,4) = drca_neighbour(row_of_old_rec,4)-...
        drca(row_of_don,4);
end

%the ids of the neighbours of donor i are:
%  i-1-wid   i-1   i-1+wid
%  i-wid      i    i+wid
%  i+1-wid   i+1   1+1+wid
%place them in an array

potential_rec=[don_to_switch-1-wid, don_to_switch-1,don_to_switch-1+wid,...
    don_to_switch-wid,don_to_switch+wid,...
    don_to_switch+1-wid, don_to_switch+1,don_to_switch+1+wid];

%remove current receiver from the list above
ind = potential_rec==drca(row_of_don,2);
potential_rec(ind) = [];

%Remove donors of donors from list above
ind_of_pot_rec_to_remove = ones(1,length(potential_rec));
for i = 1:length(potential_rec)
    row_of_pot_rec = mapping(potential_rec(i));
    while drca(row_of_pot_rec,1) ~= drca(row_of_pot_rec,2) %not outlet
        if drca(row_of_pot_rec,2) == don_to_switch % the donor is a rec of this neighbour
            ind_of_pot_rec_to_remove(i) = 0;
            break
        end
        row_of_pot_rec = mapping(drca(row_of_pot_rec,2));
    end
end

potential_rec = potential_rec(logical(ind_of_pot_rec_to_remove));



%randomly choose a new receiver
if ~isempty(potential_rec)
    rec_to_switch_ind = randi(length(potential_rec),1);
    rec_to_switch = potential_rec(rec_to_switch_ind);

    % find the row in which the receiver is a donor
    row_of_rec = mapping(rec_to_switch);

    %update the donor
    drca_neighbour(row_of_don,2) = rec_to_switch;
    %update donor catchment based on new receiver catchment
    drca_neighbour(row_of_don,3) = drca_neighbour(row_of_rec,3);
    %update length of link from don to rec
    if rec_to_switch == don_to_switch - 1 ||...
            rec_to_switch == don_to_switch + 1 || ...
            rec_to_switch == don_to_switch - wid ||...
            rec_to_switch == don_to_switch + wid
        drca_neighbour(row_of_don,5) = res;
    else
        drca_neighbour(row_of_don,5) = sqrt(2*res);
    end

    %add drainage area of donor from its new receiver and all downstream
    %nodes

    drca_neighbour(row_of_rec,4) = drca_neighbour(row_of_rec,4)+ ...
        drca_neighbour(row_of_don,4);
    row_of_temp_rec = row_of_rec;
    %temp_rec = rec_to_switch;
    while drca(row_of_temp_rec,1) ~= drca(row_of_temp_rec,2) %while not an outlet
        %temp_don = temp_rec;

        temp_rec = drca(row_of_temp_rec,2);
        row_of_temp_rec = mapping(temp_rec);
        drca_neighbour(row_of_temp_rec,4) = drca_neighbour(row_of_temp_rec,4)+...
            drca(row_of_don,4);
    end
else
    drca_neighbour=drca;
end

%toc




