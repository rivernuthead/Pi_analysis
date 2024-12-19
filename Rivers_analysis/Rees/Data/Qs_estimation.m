%Author: Marco Redolfi
%Created on: 07-Oct-2013
	
Qs_Q_fname='Qs-Q_30.txt'
%'Qs-Q.txt' is computed on the basis of dem 00 by using 20 mm as a representative sediment size
%'Qs-Q_30.txt' refers to a sediment size od 30 mm
%% Reading of Hydro_Data
data=importdata([pwd,'/../Hydro_Data/Hydrograph.txt']);
t=data.data(:,1);
Q=data.data(:,3);

	
	data=importdata([pwd,'/../Hydro_Data/Survey_times.txt']);
	
	%Points of the scan day (midnight) in the Q vector
	start_point=data.data(:,2)+1;
	end_point=data.data(:,3)+1;
	
	Delta_t=15; %Sampling time of the hydrograph [min]
	t_start=data.data(:,2)*Delta_t;
	t_end=data.data(:,3)*Delta_t;


	
	
	mid_point=fix((2*start_point+0*end_point)/2);
	
	for j=1:length(t_start)-1
		
		Qbetw=Q(mid_point(j):mid_point(j+1))
		
		Qbetw_50=Qbetw-10;
		Qbetw_50(Qbetw_50<0)=0;
		vol_W(j)=sum(Qbetw)*(3600)/1E6;
		Qmax(j)=max(Qbetw);
		
		vol_Qs1(j)=sum(Qbetw.^1.5)
		vol_Qs2(j)=sum(Qbetw.^1.7)
		
		vol_Qs3(j)=sum(Qbetw.^1.9)
		vol_Qs4(j)=sum(Qbetw.^2.1)
		
	end
	
	vol_Qs1=vol_Qs1/max(vol_Qs1);
	vol_Qs2=vol_Qs2/max(vol_Qs2);
	vol_Qs3=vol_Qs3/max(vol_Qs3);
	vol_Qs4=vol_Qs4/max(vol_Qs4);

	t_start=t_start/(60*24);
	t_end=t_end/(60*24);
	t=t/(60*24);
	
%% Plotting of the hydrograph and the survey time
	
fig=figure()
	plot(t,Q,'-r')
	grid on 
	xlabel('t [days]')
	ylabel('Q [m^3/s]')
	
	hold on
		for j=1:length(t_start)		
			plot([t_start(j),t_end(j)],[Q(start_point(j)),Q(start_point(j))],'d-b','LineWidth',2)
		end
	hold on
	xlim([0 300])

	saveas(fig,[pwd,'/Hydrograph.png'],'png')
	
%We can read a Qs(Q) file (created by Anal/Avg_sec.m), compute the Qs for each time step and integrate.



%% Estimation an plotting of solid trasport

%Qs(Q) curve caclulated with the method of the Average Section
data=importdata([pwd,'/../Hydro_Data/',Qs_Q_fname])
Q_curve=data.data(:,1);
Qs_curve=data.data(:,2);

Qs_curve(Q_curve==0)=[]; %Removal of initial null values
Q_curve(Q_curve==0)=[];

%Calculation of Qs for each curve by considering the linearly interpolated
%curve
Qs=interp1(Q_curve,Qs_curve,Q,'linear','extrap'); 

% Plotting
fig=figure
	[pl,h1,h2]=plotyy(t,Q,t,cumsum(Qs)*900/1E6)
	set(get(pl(1),'Ylabel'),'String','Q [m^3/s]') 
    set(get(pl(2),'Ylabel'),'String','Vs [10^3 m^3]') 
	set(pl(2),'ycolor','r')
	
	set(h2,'Color','r')
	set(h2,'LineWidth',1.)
	set(pl(1),'ylim',[0 450])
	set(pl(2),'ylim',[0 200])
	%set(pl(2),'YLabelColor','r')
	xlabel('t [days]')
	grid on
	set(pl(1),'xlim',[0 300])
	set(pl(2),'xlim',[0 300])
	saveas(fig,[pwd,'/Hydrograph-transport.png'],'png')
	