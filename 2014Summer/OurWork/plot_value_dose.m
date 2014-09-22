[r,c] = size(Glu_input)

figure
plot(t_gly,Glu_input','-o');
hold on
plot(t_gly,Glu_input(r,:)','linewidth',2)
xlabel('Time');ylabel('Glu');
title('Glu-Time course profile at dose = ')
hold off

figure
plot(t_gly,Ins_input,'-o')
hold on
plot(t_gly,Ins_input(r,:)','linewidth',2)
xlabel('Time');ylabel('Ins');
title('Ins-Time course profile at dose = ')
hold off

figure
plot(t_gly,Cpep_input,'-o')
hold on
plot(t_gly,Cpep_input(r,:)','linewidth',2)
xlabel('Time');ylabel('Cpep');
title('Cpep-Time course profile at dose = ')
hold off

figure
plot(t_gly,Gly_input,'-o')
hold on
plot(t_gly,Gly_input(r,:)','linewidth',2)
xlabel('Time');ylabel('Gly');
title('Gly-Time course profile at dose = ')
hold off

figure
plot(t_gly,Glu_input(r,:)',t_gly,Ins_input(r,:)',t_gly,Cpep_input(r,:)',t_gly,Gly_input(r,:)')
xlabel('Time');ylabel('Mean value');
legend('Mean Glu','Mean Ins', 'Mean Cpep', 'Mean Gly',2)
title('Mean value-Time course profile at dose = ')
