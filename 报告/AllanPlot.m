function AllanPlot(tau, data)
loglog(tau,data,'-r','LineWidth',2,...
    'MarkerFaceColor','k','MarkerSize',8);  %设置线型，颜色，宽度
xlabel('\tau (s)','fontsize',15);
ylabel('\sigma(\tau)','fontsize',15,'Rotation',90);
title('Allan deviation','fontsize',15);
grid on
axis equal 
end