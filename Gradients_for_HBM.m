

addpath(genpath('...\BrainSpace-0.1.2'))    %add the path of BrainSpace toolbox
lable = load('..\...mat');  %matrix that maps brain regions to Yeo's seven networks
sublist = textread('...\sublist.txt', '%s');  %subject list

%% group-averaged connectivity matrix
group_matrix = zeros(length(lable));
for iss = 1:length(sublist) 
    
    load([sublist{iss}, '\..mat'], 'correlation_matrix')
    temp = isnan(correlation_matrix);   %set nan to 0
    correlation_matrix(temp==1) = 0;
    
    correlation_matrix = zscore(correlation_matrix);   %normalize connectivity matrix
    group_matrix = group_matrix + correlation_matrix;
    
end
group_matrix = group_matrix/length(sublist);
group_matrix(1:length(lable)+1:end) = 0;


%% group template, Brainspace
group_matrix = real(0.5*log((1+group_matrix)./(1-group_matrix)));   %Fisher-Z tranformation
group_matrix(group_matrix<0)=0;

temp = zeros(size(group_matrix,1));    %matrix needs to be connected
temp(group_matrix == 0) = 1;
temp = sum(temp);
row_zero = find(temp == size(group_matrix,1));
group_matrix(row_zero,:) = rand(length(row_zero),size(group_matrix,1))/10000;
group_matrix(:,row_zero) = rand(size(group_matrix,1), length(row_zero))/10000;

for ir = 1 : size(group_matrix,1)
    [Y,~] = sort(group_matrix(:,ir));
    group_matrix(group_matrix(:,ir) < Y(round(size(group_matrix,1)*0.9)), ir) = 0;     %remain top 10% connections
end
gm_temp = GradientMaps('kernel', 'cs', 'approach', 'dm','n_components',10);
gm_temp = gm_temp.fit(group_matrix');


gm_p1_system_mean = zeros(length(sublist), 7);
gm_p2_system_mean = zeros(length(sublist), 7);
gm_p1_ratio = zeros(length(sublist), 1);
gm_p2_ratio = zeros(length(sublist), 1);
gm_within_system_dispersion = zeros(length(sublist), 7);
gm_between_system_dispersion = zeros(length(sublist), 7);
for ija = 1 : length(sublist)

    
    %% indivudal gradients
    load([sublist{ija}, '\...mat'], 'correlation_matrix')
    correlation_matrix = zscore(correlation_matrix);
    s_corrcoef = correlation_matrix;
    temp = isnan(s_corrcoef);
    s_corrcoef(temp==1) = 0;
    
    %PA
    s_corrcoef(1:size(s_corrcoef,1)+1:end) = 0;
    s_corrcoef = real(0.5*log( (1+s_corrcoef)./(1-s_corrcoef) ));    
    s_corrcoef(s_corrcoef<0)=0;
    
    temp = zeros(size(s_corrcoef,1));
    temp(s_corrcoef == 0) = 1;
    temp = sum(temp);
    row_zero = find(temp == size(s_corrcoef,1));
    s_corrcoef(row_zero,:) = rand(length(row_zero),size(s_corrcoef,1))/10000;
    s_corrcoef(:,row_zero) = rand(size(s_corrcoef,1), length(row_zero))/10000;
    
    for ir = 1 : size(s_corrcoef,1)
        [Y,~] = sort(s_corrcoef(:,ir));
        s_corrcoef(s_corrcoef(:,ir) < Y(round(size(s_corrcoef,1)*0.9)), ir) = 0;
    end
    
    
    gm = GradientMaps('kernel', 'cs', 'approach', 'dm', 'alignment', 'pa');  %two alignment techniques: Procrustes analysis (pa) and joint alignment (ja)
    gm = gm.fit(s_corrcoef', 'reference', gm_temp.gradients{1});

    gm_p1_ratio(ija) = gm.lambda{1,1}(1)/sum(gm.lambda{1,1});
    gm_p2_ratio(ija) = gm.lambda{1,1}(2)/sum(gm.lambda{1,1});
    for ibs = 1 : 7
        
        col = find(lable == ibs);
        gm_p1_system_mean(ija, ibs) = mean(gm.aligned{1,1}(col,1));
        gm_p2_system_mean(ija, ibs) = mean(gm.aligned{1,1}(col,2));
    
        
        centroid = median(gm.aligned{1,1}(col,1:2));
        centroid_temp = median(gm_temp.gradients{1,1}(col,1:2));
        distance = zeros(length(col),1);
        distance_temp = zeros(length(col),1);
        for idc = 1 : length(col)
            distance(idc) = ((gm.aligned{1,1}(col(idc),1)-centroid(1))^2 + (gm.aligned{1,1}(col(idc),2)-centroid(2))^2 )^0.5;
            distance_temp(idc) = ((gm.aligned{1,1}(col(idc),1)-centroid_temp(1))^2 + (gm.aligned{1,1}(col(idc),2)-centroid_temp(2))^2 )^0.5;
        end
        gm_within_system_dispersion(ija, ibs) = mean(distance);
        clear centroid_temp distance distance_temp
        
        distance = zeros(7,1);
        for ibss = 1 : 7
            col_temp = find(lable == ibss);
            centroid_ind = median(gm.aligned{1,1}(col_temp,1:2));
            distance(ibss) = ( (centroid(1)-centroid_ind(1))^2 + (centroid(2)-centroid_ind(2))^2 )^0.5;
     
        end
        gm_between_system_dispersion(ija, ibs) = sum(distance)/6;
        clear distance centroid_temp
    end
    
end