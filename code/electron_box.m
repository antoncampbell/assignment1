clc
%close all
clear all

%General Program Parameters
global G;
G.Mode=6;
% %Mode=1: 1c)i)
% %Mode=2: 1c)ii)
% 
% %Mode=3: 2a)
% %Mode=4: 2b)
% %Mode=5: 2c)
% 
% %Mode=6: 3a)
% %Mode=7: 3b)
% %Mode=8: 3c)

num_elec=10000;
num_elec_plot=10;%64;
pause_time=0.001;
multiplier=5;
old_plot_elec=zeros(num_elec,2);

%Assignment Parameters
m0=9.10938356e-31; % [kg]
mn=0.26*m0; % [kg]
kB=1.38064852e-23; % [J*K^-1]
Temp=300; % [K]
Tau_mn=0.2e-12; % [s]
global P;
P.box_size_x=200e-9; % [m]
P.box_size_y=100e-9; % [m]
P.thermal_speed=sqrt(kB*Temp/mn); %FIXME: is for 1d -- may be right
P.time_step=(P.box_size_y/(multiplier*100))/P.thermal_speed;
P.distancemax=6*(P.box_size_y/(multiplier*100));
P.Pscat=1-exp(-P.time_step/Tau_mn);
P.MFP=sqrt(2)*P.thermal_speed*Tau_mn;

% Matrix Names
global S;
S.px=1; S.py=2; S.vx=3; S.vy=4; S.num_dt=5;


% Matrix to Track path lengths
scat_info=zeros(1,2);



%% Generate Boxes
number_boxes=1;
% boxes=zeros(1,2,4);
% boxes(1,1,:)=P.box_size_x.*[.2 .4 .4 .5];%
% boxes(1,2,:)=P.box_size_y.*[0 .2 .2 0];
boxes=zeros(1,5);
boxes(1,:)=P.box_size_y.*[-0.5 1 3 1 0]; %top boundary
boxes(2,:)=P.box_size_y.*[-0.5 -1 3 1 0]; %bottom boundary
if(G.Mode==6||G.Mode==7||G.Mode==8)
    boxes(3,:)=P.box_size_y.*[0.8 -0.1 0.4 0.5 1/P.box_size_y];% x0 y0 width height is_diffuse --- (x0 y0) is bottom left
    boxes(4,:)=P.box_size_y.*[0.8 0.6 0.4 0.5 0];% x0 y0 width height is_diffuse --- (x0 y0) is bottom left
end


%% Create matrix of electrons
electrons=zeros(num_elec,5);
for i1=1:num_elec
    valid=false;
    while(~valid)
        electrons(i1,:)=create_electon;
        valid=true;%assume valid after creating electron
        for ii=1:size(boxes,1)
            e_y=electrons(i1,S.py);
            e_x=electrons(i1,S.px);
            x0=boxes(ii,1);
            y0=boxes(ii,2);
            x1=boxes(ii,1)+boxes(ii,3);
            y1=boxes(ii,2)+boxes(ii,4);
            inside_box=(e_y>y0)&(e_y<y1)&(e_x>x0)&(e_x<x1);
            valid=~inside_box&valid;
        end
    end
    
end
%TEMP
%electrons(1,:)=[0.2*P.box_size_y 0.3*P.box_size_y P.thermal_speed 0];





%% Run for over time and plot
figure(G.Mode)

%hold on when plotting particle trajectories
if(G.Mode==1||G.Mode==4||G.Mode==6)
    hold on;
end

%plot boxes
if(G.Mode==6)
    for ii=1:size(boxes,1)
        x0=boxes(ii,1);
        y0=boxes(ii,2);
        x1=boxes(ii,1)+boxes(ii,3);
        y1=boxes(ii,2)+boxes(ii,4);
        plot([x0 x1 x1 x0 x0],[y0 y0 y1 y1 y0],'Color','black');
    end
end

%initialize matrix that keeps track of limited number of electrons
for eii=1:num_elec_plot
    old_plot_elec(eii,1)=electrons(eii,S.px);
    old_plot_elec(eii,2)=electrons(eii,S.py);
end
color_map=hsv(num_elec_plot);

%loop over time
for t=1:1000*multiplier
   time=t*P.time_step;
   [electrons,scat_info]=move_electrons(electrons,boxes,scat_info);

   %Calculate Temperature
   v_speed=sqrt(electrons(:,S.vx).^2+electrons(:,S.vy).^2);
   vx_MS=sum(electrons(:,S.vx).^2)/size(electrons,1);
   vy_MS=sum(electrons(:,S.vy).^2)/size(electrons,1);   
   vt_MS=vx_MS+vy_MS;
   Temp_Calc=mn*vt_MS/(2*kB);
   
   Time_Temp(t)=Temp_Calc;
   
   %Plot '2-D plot of particle trajectories'
   if(G.Mode==1||G.Mode==4||G.Mode==6)
       %plot
       if(mod(t,multiplier)==0)
           for eii=1:num_elec_plot
    %            subplot(2,1,1);
    %            hold on;
    %            figure(1)
    %             hold on;
                loop_around=abs(old_plot_elec(eii,1)-electrons(eii,S.px))>(P.box_size_x/2);
                old_plot_elec(eii,1)=~loop_around.*old_plot_elec(eii,1)+loop_around.*electrons(eii,S.px);
    %             plot([electrons(eii,S.ox) electrons(eii,S.px)],...
    %                  [electrons(eii,S.oy) electrons(eii,S.py)],...
    %                  'LineStyle','-','Color',color_map(eii,:)) 
                plot([old_plot_elec(eii,1) electrons(eii,S.px)],...
                     [old_plot_elec(eii,2) electrons(eii,S.py)],...
                     'LineStyle','-','Color',color_map(eii,:))

                old_plot_elec(eii,1)=electrons(eii,S.px);
                old_plot_elec(eii,2)=electrons(eii,S.py);
    %             hold off;
           end
       end
       %old_plot_elec=(
    %    plot(electrons(:,S.px),electrons(:,S.py),'o')

       xlim([0 P.box_size_x])
       ylim([0 P.box_size_y])
    %    xlim([-0.1*P.box_size_x 1.1*P.box_size_x])
    %    ylim([-0.1*P.box_size_y 1.1*P.box_size_y])
   elseif(G.Mode==2||G.Mode==5)
        %plot(Time_Temp)
        plot(P.time_step:P.time_step:t*P.time_step,Time_Temp)
        title('Tempeture vs. Time')
   elseif(G.Mode==3)
%        xlim([0 6*P.thermal_speed])
       histogram(v_speed,100)
       xlim([0 6*P.thermal_speed])
       ylim([0 400])   
       
       
   elseif(G.Mode==7)
       hist3([electrons(:,S.px),electrons(:,S.py)],[50 50],'CdataMode','auto')
       view(0,90)
       
   elseif(G.Mode==8)
%            vx_MS=sum(electrons(:,S.vx).^2)/size(electrons,1);
%            vy_MS=sum(electrons(:,S.vy).^2)/size(electrons,1);   
%            vt_MS=vx_MS+vy_MS;
%            Temp_Calc=mn*vt_MS/(2*kB);
       
       
       size_div=P.box_size_y/50;
       numb_x=round(P.box_size_x/size_div)+1; numb_y=round(P.box_size_y/size_div)+1;
       temp_map=zeros(numb_x,numb_y,2);
       for ee=1:size(electrons,1)
           x_round=ceil(electrons(ee,S.px)./size_div)+1;
           y_round=ceil(electrons(ee,S.py)./size_div)+1;
           temp_map(x_round,y_round,1)=temp_map(x_round,y_round,1)+(electrons(ee,S.vx).^2)+(electrons(ee,S.vy).^2);
           temp_map(x_round,y_round,2)=temp_map(x_round,y_round,2)+1;
       end
       temp_map_val=temp_map(:,:,1)./(temp_map(:,:,2)+0.01);
       surf(temp_map_val');
       view(0,90)
   end
   

   pause(pause_time)

   
end

%% Perform Final Evaluation
if(G.Mode==1)
    fprintf('The thermal velocity is: %4.3f m/s\n',P.thermal_speed);
    fprintf('The mean free path is: %4.3s m\n',P.MFP);
elseif(G.Mode==5)
    tmp_scat_info=scat_info(2:end,:);
    Calc_Tau_mn=mean(P.time_step.*tmp_scat_info(:,1))
    Calc_MFP=mean(P.time_step.*tmp_scat_info(:,1).*tmp_scat_info(:,2))
end

%FIXME: fill in once modes are changed
% if(G.Mode==1)
%     
% elseif(G.Mode==2)
%     
% elseif(G.Mode==3)
% 
% end

%%
function [electron]=create_electon
    global G;
    global S;
    global P;
    

    electron(1,S.px)=rand()*P.box_size_x;
    electron(1,S.py)=rand()*P.box_size_y;

    if((G.Mode==1)||(G.Mode==2))
        angle=2*pi*rand();
        speed=sqrt(2)*P.thermal_speed;
        electron(1,S.vx)=speed*cos(angle);
        electron(1,S.vy)=speed*sin(angle);
    elseif(G.Mode>=3)
        electron(1,S.vx)=P.thermal_speed*randn();
        electron(1,S.vy)=P.thermal_speed*randn();
    end
    electron(1,S.num_dt)=0;

end

%%
function [electrons,scat_info]=move_electrons(electrons,boxes,scat_info)
    global G;
    global S;
    global P;

    
    
     new_py=electrons(:,S.py)+electrons(:,S.vy)*P.time_step;
     new_px=electrons(:,S.px)+electrons(:,S.vx)*P.time_step;
     
     one=1;
     distancemax=P.distancemax;
     inside_box=zeros(size(electrons,1),1);
     bottom=inside_box; top=inside_box; left=inside_box; right=inside_box; diffusive=inside_box;
     for ii=1:size(boxes,1)
        x0=boxes(ii,1);
        y0=boxes(ii,2);
        x1=boxes(ii,1)+boxes(ii,3);
        y1=boxes(ii,2)+boxes(ii,4);
        inside_box=ii.*((new_py>y0)&(new_py<y1)&(new_px>x0)&(new_px<x1)&(inside_box==0))+one.*inside_box;

%FIXME: remove; was used for debugging
%         tmp=0;
%         for jj=1:length(inside_box)
%             tmp=tmp+inside_box(jj);
%         end
%         if(tmp>=1)
%             
%            xxx=1;
%         end
        
       
        bottom=(inside_box==ii)&(abs(new_py-y0)<distancemax)|bottom;
        top=(inside_box==ii)&(abs(new_py-y1)<distancemax)|top;
        left=(inside_box==ii)&(abs(new_px-x0)<distancemax)|left;
        right=(inside_box==ii)&(abs(new_px-x1)<distancemax)|right;
%         left=(inside_box==ii)&(abs(new_px-x0)<distancemax)&~bottom&~top|left;
%         right=(inside_box==ii)&(abs(new_px-x1)<distancemax)&~bottom&~top|right;
        %diffusive=(inside_box==ii).*boxes(ii,5)|diffusive;
     end
     
     diffusive=(inside_box>=1).*boxes(inside_box+double((inside_box==0)),5);
     
     
    if (G.Mode>=6)      
        b_rand_abs=abs(randn(size(electrons,1),1));
        b_rand_nor=(randn(size(electrons,1),1));
        
        bounce_vertical=(inside_box>=1)&(top|bottom);
        bounce_lateral=(inside_box>=1)&(left|right);       
        
        %x dim
        electrons(:,S.vx)=~diffusive.*(~bounce_lateral.*electrons(:,S.vx)-bounce_lateral.*electrons(:,S.vx))...
                           +diffusive.*((inside_box==0).*electrons(:,S.vx)...
                                       -bounce_lateral.*b_rand_abs.*P.thermal_speed.*(sign(electrons(:,S.vx)))...
                                       +bounce_vertical.*b_rand_nor.*P.thermal_speed);
        electrons(:,S.px)=mod(electrons(:,S.px)+electrons(:,S.vx)*P.time_step, P.box_size_x);


        %y dim
        electrons(:,S.vy)=~diffusive.*(~bounce_vertical.*electrons(:,S.vy)-bounce_vertical.*electrons(:,S.vy))...
                          +diffusive.*((inside_box==0).*electrons(:,S.vy)...
                                       -bounce_vertical.*b_rand_abs.*P.thermal_speed.*(sign(electrons(:,S.vy)))...
                                       +bounce_lateral.*b_rand_nor.*P.thermal_speed);
        electrons(:,S.py)=electrons(:,S.py)+electrons(:,S.vy)*P.time_step;
        
    else
        %x dim
        bounce_lateral=(inside_box>=1)&(left|right);
        electrons(:,S.vx)=~bounce_lateral.*electrons(:,S.vx)-bounce_lateral.*electrons(:,S.vx);
        electrons(:,S.px)=mod(electrons(:,S.px)+electrons(:,S.vx)*P.time_step, P.box_size_x);


        %y dim
        bounce_vertical=(inside_box>=1)&(top|bottom);
        electrons(:,S.vy)=~bounce_vertical.*electrons(:,S.vy)-bounce_vertical.*electrons(:,S.vy);
        electrons(:,S.py)=electrons(:,S.py)+electrons(:,S.vy)*P.time_step;
    end

    %% Scattering Functionality
    if (G.Mode>=3)
        scattering=(P.Pscat>rand(size(electrons,1),1));
        rand_x=randn(size(electrons,1),1);
        rand_y=randn(size(electrons,1),1);
        
        if (G.Mode==5)
            for ss=1:length(scattering)
                if(scattering(ss)==1)
                    num_dt=electrons(ss,S.num_dt);
                    speed_1=sqrt(electrons(ss,S.vx).^2+electrons(ss,S.vy).^2);
                    scat_info(end+1,:)=[num_dt speed_1];
                end
            end
        end
        
        electrons(:,S.num_dt)=~scattering.*(electrons(:,S.num_dt)+1)+scattering.*0;
        electrons(:,S.vx)=scattering.*P.thermal_speed.*rand_x+~scattering.*electrons(:,S.vx);
        electrons(:,S.vy)=scattering.*P.thermal_speed.*rand_y+~scattering.*electrons(:,S.vy);

    end

    

end
