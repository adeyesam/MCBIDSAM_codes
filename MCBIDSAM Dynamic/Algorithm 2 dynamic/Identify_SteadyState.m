function sel_indx = Identify_SteadyState(data,u)
%% identify steady state points 
% Here, steadiness of all input and output variables is used as basis for checking for steady state
szp = 10;                                           % Size of each piece
spn = 0.001;                                        % Min rel change to be considered off-steady state
np = floor(size(data,2)/szp);                       % Number of pieces
sel_indx = [];                                      % Indices of selected datapoints
yss = [];
uss = [];
N = size(data,1);
lengthu = size(u,1);

for i = 1:np
    indxp = (i-1)*szp+1:i*szp;                      % Indices of elts in current block
    
    up = u(:,indxp);
    yp = data(:,indxp);
    for iu = 1:lengthu
        mean_up(iu,i) = mean(up(iu,:));             % Mean of current piece
        range_up(iu,i) = range(up(iu,:));           % Range of current piece
    end
    for iy = 1:N
        mean_yp(iy,i) = mean(yp(iy,:));             % Mean of current piece
        range_yp(iy) = range(yp(iy,:));             % Range of current piece
    end

    if i==1
        dmean = true;
        if all(range_yp>spn) || all(dmean)          % Selection based on either large range or diffence in mean of successive pieces
            sel_indx = [sel_indx indxp(end)];
            uss = [uss up(:,1)];
            yss = [yss mean_yp(:,i)];
        end
    else
        dmean = abs((mean_yp(:,i)-mean_yp(:,i-1))./mean_yp(:,i))>spn;           % Is difference in mean of consecutive pieces above threshold?
        if all(range_up(:,i)==0) && ~all(range_up(:,i)==range_up(:,i-1))        % No change in any input var and not same as previous piece
            if all(range_yp>spn) || all(dmean)                                  % Selection based on either large range or diffence in mean of successive pieces
                sel_indx = [sel_indx indxp(1)];
                uss = [uss up(:,1)];
                yss = [yss mean_yp(:,i)];
            end
        end
    end
end


% data_ss = data(:,sel_indx);          
% u_ss = u(:,sel_indx);                

% figure;plot(data(1,:))
% xline(sel_indx)
end