plot(t_gly,Mean_Glu,'-o',t_gly,Mean_Ins,'-o',t_gly,Mean_Cpep,'-o',t_gly,Mean_Gly,'-o')
plot(t_gly,Glu_input,'-o')
plot(t_gly,Ins_input,'-o')
plot(t_gly,Cpep_input,'-o')
plot(t_gly,Gly_input,'-o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plotregression(Mean_Gly,Mean_Glu)

figure
plotregression( Mean_Gly, Mean_Ins)

figure
plotregression(Mean_Gly, Mean_Cpep)

