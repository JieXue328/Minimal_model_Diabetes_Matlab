%% mean glucose over doses
load('Gly_data_input_d125.mat')
[r,c] = size(Glu_input)

figure
plot(t_gly,Glu_input(r,:)','r-','linewidth',1.5)
clear
hold on 

load('Gly_data_input_d250.mat')
[r,c] = size(Glu_input)
plot(t_gly,Glu_input(r,:)','b-','linewidth',1.5)
clear

load('Gly_data_input_d500.mat')
[r,c] = size(Glu_input)
plot(t_gly,Glu_input(r,:)','g-','linewidth',1.5)
clear

load('Gly_data_input_d750.mat')
[r,c] = size(Glu_input)
plot(t_gly,Glu_input(r,:)','y-','linewidth',1.5)
clear

load('Gly_data_input_d1000.mat')
[r,c] = size(Glu_input)
plot(t_gly,Glu_input(r,:)','c-','linewidth',1.5)
clear

xlabel('Time');ylabel('Meanvalue');
title('Comparison of MeanGlu at different doses')
legend ('dose=1.25','dose=2.5','dose=5','dose=7.5','dose=10',2 )

%% mean insulin over doses
load('Gly_data_input_d125.mat')
[r,c] = size(Glu_input)

figure
plot(t_gly,Ins_input(r,:)','r-','linewidth',1.5)
clear
hold on 

load('Gly_data_input_d250.mat')
[r,c] = size(Glu_input)
plot(t_gly,Ins_input(r,:)','b-','linewidth',1.5)
clear

load('Gly_data_input_d500.mat')
[r,c] = size(Glu_input)
plot(t_gly,Ins_input(r,:)','g-','linewidth',1.5)
clear

load('Gly_data_input_d750.mat')
[r,c] = size(Glu_input)
plot(t_gly,Ins_input(r,:)','y-','linewidth',1.5)
clear

load('Gly_data_input_d1000.mat')
[r,c] = size(Glu_input)
plot(t_gly,Ins_input(r,:)','c-','linewidth',1.5)
clear

xlabel('Time');ylabel('Meanvalue');
title('Comparison of MeanIns at different doses')
legend ('dose=1.25','dose=2.5','dose=5','dose=7.5','dose=10',2 )

%% mean Cpeptide over doses
load('Gly_data_input_d125.mat')
[r,c] = size(Glu_input)

figure
plot(t_gly,Cpep_input(r,:)','r-','linewidth',1.5)
clear
hold on 

load('Gly_data_input_d250.mat')
[r,c] = size(Glu_input)
plot(t_gly,Cpep_input(r,:)','b-','linewidth',1.5)
clear

load('Gly_data_input_d500.mat')
[r,c] = size(Glu_input)
plot(t_gly,Cpep_input(r,:)','g-','linewidth',1.5)
clear

load('Gly_data_input_d750.mat')
[r,c] = size(Glu_input)
plot(t_gly,Cpep_input(r,:)','y-','linewidth',1.5)
clear

load('Gly_data_input_d1000.mat')
[r,c] = size(Glu_input)
plot(t_gly,Cpep_input(r,:)','c-','linewidth',1.5)
clear

xlabel('Time');ylabel('Meanvalue');
title('Comparison of MeanCpep at different doses')
legend ('dose=1.25','dose=2.5','dose=5','dose=7.5','dose=10',2 )

%% mean GLY over doses
load('Gly_data_input_d125.mat')
[r,c] = size(Glu_input)

figure
plot(t_gly,Gly_input(r,:)','r-','linewidth',1.5)
clear
hold on 

load('Gly_data_input_d250.mat')
[r,c] = size(Glu_input)
plot(t_gly,Gly_input(r,:)','b-','linewidth',1.5)
clear

load('Gly_data_input_d500.mat')
[r,c] = size(Glu_input)
plot(t_gly,Gly_input(r,:)','g-','linewidth',1.5)
clear

load('Gly_data_input_d750.mat')
[r,c] = size(Glu_input)
plot(t_gly,Gly_input(r,:)','y-','linewidth',1.5)
clear

load('Gly_data_input_d1000.mat')
[r,c] = size(Glu_input)
plot(t_gly,Gly_input(r,:)','c-','linewidth',1.5)
clear

xlabel('Time');ylabel('Meanvalue');
title('Comparison of MeanGly at different doses')
legend ('dose=1.25','dose=2.5','dose=5','dose=7.5','dose=10',2)

