clear;
close all;
fignum=1;
vis=1;

%path:
folder='/Users/';  titletag = ' '; in_vel = 10; %in mm/yr or km/Myr;

%field='nodemat';outname='matplot_';collabel='material type';       tclr='r';         %
%field='nodtemp';outname='tempplot_';collabel='temperature (°C)';    tclr='w';         % 
%field='nodevel';outname='vxplot_';collabel='horizontal velocity (mmyr^{-1})'; tclr='k'; dir=1;  % horizontal velocity
field='nodevel';outname='vyplot_';collabel='vertical velocity (mmyr^{-1})'; tclr='k'; dir=0;  % vertical velocity
% field='nodtempin';outname='tempplot_';collabel='temperature (�C)';    tclr='w';         % 


%starttime   %endtime
tmin=200;   tmax=200;
%output frequency:
step = 200; 

% plot boundaries (in km)
z = 0; % zoom or not
if z==1
    xmin =  600; xmax = 1800;
    ymin = -160; ymax = 20;    
else; xmin = -0; xmax = 3080; ymin = -660; ymax = 20;
end

%temperature flag
itemp=1;
geotherms=[100 300 500 700 900 1100 1200 1300]; % (degree C)
%velocity flag
ivelo=0;  vrate=20;  S=1; %vrate sets resampling rate, S scales the arrow lengths
%string flag
istring =1;


%file output
iout=1;
%filename tag (optional)
tag='';
%additional automatic tags
if z==1;Ztag='zoom';else;Ztag='';end
if itemp==1;Ttag='T';else;Ttag='_';end
if ivelo==1;vtag='v';else;vtag='_';end

%model size
n=161;
m=587; 

%material colormap (Cramerie 2018)
mymap = ([          0.78039 0.69804 0.60000;   % continental upper crust [light tan]
    0.00000 0.57255 0.27059;   % lith [green]
          0.77647 0.61176 0.42745;   % continental lower crust [dark tan]=
          1.00000 1.00000 1.00000;   % slm [white
          0.98431 0.69020 0.23137;   % continental sediment [dark-mustard]
          0.98431 0.69020 0.23137;   % oceanic sediments [dark-mustard]    
          0.00000 0.57255 0.27059;   % cml [green]
          0.96863 0.57647 0.11765;   % terrain upper crust [orange]
          0.45882 0.29804 0.14118;   % terrain mid-crust [mid-brown]
          0.25882 0.12941 0.04314]); % terrain lower crust [mid-brown]

%main loop
for i=tmin:step:tmax
    % building standard time format
    if     i<10;    timeform=['00000' num2str(i)];
    elseif i<100;   timeform=['0000'  num2str(i)];
    elseif i<1000;  timeform=['000'   num2str(i)]; 
    elseif i<10000; timeform=['00'    num2str(i)];  
    elseif i<100000;timeform=['0'     num2str(i)];          
    end

    % material and nodal read
    file=strcat(folder,field,timeform,'.out');
    temp=importdata(file," ",1);
    data=getfield(temp,'data');
    % getting timestep info out of header
    blabla=getfield(temp,'textdata');bla=char(blabla);
    bl=split(bla,'   ');
    dtstep=str2num(bl{8});dtstep=round(dtstep/60/60/24/365/1e6,2);
    clear temp
    x=reshape(data(:,1),[n,m]);
    y=reshape(data(:,2),[n,m]);
    if isequal(field,'nodemat');plotdata=reshape(data(:,3),[n,m]);end
    if isequal(field,'nodtemp');plotdata=reshape(data(:,3),[n,m]);end
    if (isequal(field,'nodevel') && dir==1);plotdata=reshape(data(:,3),[n,m]);plotdata=plotdata*60*60*24*365*1000;end
    if (isequal(field,'nodevel') && dir==0);plotdata=reshape(data(:,4),[n,m]);plotdata=plotdata*60*60*24*365*1000;end
    %clear data

    % temperature read
    if (itemp==1 && isequal(field,'nodtemp'))
        temperature=plotdata;
    elseif itemp==1
        file=strcat(folder,'nodtemp',timeform,'.out');
        temp=importdata(file," ",1);
        data=getfield(temp,'data');
        clear temp
        temperature=reshape(data(:,3),[n,m]);
        clear data
    end
    
    % nodal velocity read
    if ivelo==1
        file=strcat(folder,'nodevel',timeform,'.out');
        temp=importdata(file," ",1);
        data=getfield(temp,'data');
        clear temp

        vx=reshape(data(:,3),[n,m]);vx=vx*60*60*24*365*1000;
        vy=reshape(data(:,4),[n,m]);vy=vy*60*60*24*365*1000;

        velo_x=x(1:vrate:n,1:vrate:m);
        velo_y=y(1:vrate:n,1:vrate:m);
        velo_vy=vy(1:vrate:n,1:vrate:m);
        velo_vx=vx(1:vrate:n,1:vrate:m);
        clear data
    end
    
    if istring==1
        file=strcat(folder,'strings',timeform,'.out');
        temp=importdata(file," ",1);
        string=getfield(temp,'data');
        clear temp
    end
    
    % extracting the top of the model for plotting
    surface(1,:)=x(1,:)/1000;
    surface(2,:)=y(1,:)/1000;
    
    
    % Plotting
    fhg = figure(fignum);

    pcolor(x/1000,y/1000,plotdata); shading flat; axis equal tight;
    if isequal(field,'nodemat');colormap(mymap);caxis(caxis);end
    if isequal(field,'nodtemp');load('buda.mat');colormap(buda);caxis(caxis);colorbar;ylabel(colorbar,collabel);end
    if isequal(field,'nodevel');colormap('redblue');caxis([-250 250]);colorbar;ylabel(colorbar,collabel);end
    hold on;
    plot(surface(1,:),surface(2,:),'k');
    if istring ==1; hold on; plot(string(:,1)/1000, string(:,2)/1000, '--b'); end 
    if itemp==1;[C,h]=contour(x/1000,y/1000,temperature,geotherms,tclr,'LineWidth',0.75);clabel(C,h,'LabelSpacing',1700,'color',tclr);end
    if ivelo==1;quiver(velo_x/1000,velo_y/1000,velo_vx,velo_vy,S,'k');end
    xlim([xmin xmax]);ylim([ymin ymax]);
    step=num2str(dtstep); conv=num2str(dtstep*in_vel);
    xlabel('Along profile distance (km)');ylabel('Depth (km)');
    title([titletag step ' Myr; ']);

    if iout==1;nom=strcat(strcat(folder,outname,Ttag,vtag,Ztag,tag,'_',timeform));print('-dpng','-r600','-painters',nom);end
end