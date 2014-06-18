plot(t_gly,Mean_Glu,'-o',t_gly,Mean_Insulin,'-o',t_gly,Mean_Cpep,'-o',t_gly,Mean_Gly,'-o')
plot(t_gly,Glu_input,'-o')
plot(t_gly,Insulin_input,'-o')
plot(t_gly,Cpept_input,'-o')
plot(t_gly,Gly_input,'-o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotregression(Mean_Glu, Mean_Insulin)
plotregression(Mean_Glu, Mean_Cpep)
plotregression(Mean_Insulin, Mean_Cpep)