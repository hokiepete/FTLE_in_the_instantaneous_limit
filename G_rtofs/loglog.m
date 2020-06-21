close all
clear all
clc
width = 1%#5+3/8
height = 4.5/1.61803399
set(0,'defaulttextinterpreter','latex')
figSize = [width, height]
load wrf_plot_rmse_data.mat
%t = log(abs(time(2:end)));
%rmse3 = log(rmse3(2:end)*3600);
%rmse2 = log(rmse2(2:end)*3600);
%rmse1 = log(rmse1(2:end)*3600);

t = abs(time(2:end));
rmse3 = rmse3(2:end)*3600;
rmse2 = rmse2(2:end)*3600;
rmse1 = rmse1(2:end)*3600;


figure('units','inch','position',[1,1,width,height]);
hold on
plot(t,rmse1,'b-','LineWidth',2)
plot(t,rmse2,'m-','LineWidth',2)
plot(t,rmse3,'k-','LineWidth',2)
ylabel('FTLE field root mean-squared error','FontSize',14)
xlabel('$|T|$','FontSize',14)
axis('tight')
grid
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

%ticklabel_format(style='sci',axis='y', scilimits=(0,0))

%plt.savefig('rtofs_rmse_log.png', transparent=False, bbox_inches='tight',pad_inches=0.03)
%plt.savefig('rtofs_rmse_log.eps', transparent=False, bbox_inches='tight',pad_inches=0.03)