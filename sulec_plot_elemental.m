close;
clear;
vis=1;

% path:
 folder='/Users/';  titletag = 'Initial collision test: '; in_vel = 50; fignum=2; %in mm/yr or km/Myr

% which field to plot
field='elvisco';outname='visco_';     collabel='log_1_0(\mu_e_f_f)';          tclr='r';         % log10 viscosity
%field='elpstrn';outname='pstrn_';     collabel='cummulative plastic strain';  tclr='r';         % plastic strain
%field='elstrrt';outname='strrt_';     collabel='log_1_0 (strain rate)';       tclr='w';         % strainrate tensor (second invariant)
%field='eldensi';outname='density';    collabel='\rho (kgm^{-3})';             tclr='w';         % density
% field='elsigma';outname='sigma_';     collabel='deviatoric stress';           tclr='r';         % deviatoric stress (second invariant)
% field='eltempe';outname='eltemp_';    collabel='temperature (°C)';            tclr='w';         % temperature
% field='elpress';outname='press_';     collabel='pressure (MPa)';              tclr='r';         % pressure
% field='elprdif';outname='diffpress_'; collabel='differential pressure (MPa)'; tclr='r';         % differential pressure
%field='eltresc';outname='stress_';    collabel='differential stress';         tclr='r';         % differential stress
 %field='eldissp';outname='dissip_';    collabel='dissipation';                tclr='r';          %dissipation rate
%file
% starttime % endtime
tmin=2000;   tmax=2000;
% output frequencny:
step = 200;

% plot boundaries (in km)
z = 0; % zoom or not
if z==1
    xmin =  800; xmax = 1300;
    ymin = -160; ymax = 20;    
else; xmin = 0; xmax = 3080; ymin = -660; ymax = 20; 
end

%temperature flag
itemp=1;
geotherms=[100 500 900 1300]; % (degree C)
%velocity flag
ivelo=1;     vrate=30; %vrate sets resampling rate
iparallelCreep=0;
istring =1;

%filename tag (optional)
tag="";
%additional automatic tags
if z==1;Ztag='zoom';else;Ztag='';end
if itemp==1;Ttag='T';else;Ttag='_';end
if ivelo==1;vtag='v';else;vtag='_';end

%minmum-maximum value plotted in a textbox onto the figure
annote=1;

%output saving flag
iout=0;

%domain size
n=161; nel=n-1;
m=587; mel=m-1;

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
    elx=reshape(data(:,1),[mel,nel]);elx=elx';
    ely=reshape(data(:,2),[mel,nel]);ely=ely';
    elfield=reshape(data(:,3),[mel,nel]);elfield=elfield';
    clear data
    if isequal(field,'elstrrt');elfield=log10(elfield);end
    
    % temperature read
    if itemp==1
        file=strcat(folder,'nodtemp',timeform,'.out');
        temp=importdata(file," ",1);
        data=getfield(temp,'data');
        clear temp
        %data=table2array(nodtemp);
        x=reshape(data(:,1),[n,m]);
        y=reshape(data(:,2),[n,m]);
        temperature=reshape(data(:,3),[n,m]);
        clear data
    end
    
    % nodal velocity read
    if ivelo==1
        file=strcat(folder,'nodevel',timeform,'.out');
        temp=importdata(file," ",1);
        data=getfield(temp,'data');
        clear temp
        %data=table2array(nodevel);
        x=reshape(data(:,1),[n,m]);
        y=reshape(data(:,2),[n,m]);
        vx=reshape(data(:,3),[n,m]);
        vy=reshape(data(:,4),[n,m]);
        velo_x=x(1:vrate:n,1:vrate:m);
        velo_y=y(1:vrate:n,1:vrate:m);
        velo_vy=vy(1:vrate:n,1:vrate:m);
        velo_vx=vx(1:vrate:n,1:vrate:m);
       clear data
    end
    
    % extracting the top of the model for plotting
    surface(1,:)=elx(1,:)/1000;
    surface(2,:)=ely(1,:)/1000;
    
    if istring==1
        file=strcat(folder,'strings',timeform,'.out');
        temp=importdata(file," ",1);
        string=getfield(temp,'data');
        clear temp
    end

    fhg = figure(fignum);
    if (vis); set(fignum,'position',[1 41 1920 963]); else; set(fignum,'position',[1 41 1920 963],'Visible','off'); end    
    pcolor(elx/1000,ely/1000,elfield); shading flat; axis equal tight;
    xlabel('Along profile distance (km)');ylabel('Depth (km)');
    if isequal(field,'eldensi');load('bamako.mat');colormap(bamako);caxis(caxis);caxis([3000 3300]);end
    if isequal(field,'elvisco');load('oslo.mat');colormap(flipud(oslo));caxis(caxis);caxis([19 30]);end
    if isequal(field,'elstrrt');load('buda.mat');colormap(buda);caxis(caxis);caxis([-21 -12]);end
    if isequal(field,'elpstrn');colormap(flipud(gray));caxis(caxis);caxis([0 1.5]);end
    if isequal(field,'eltresc');load('oslo.mat');colormap(flipud(oslo));caxis(caxis);caxis([0 2e9]);end
    if isequal(field,'elsigma');load('bamako.mat');colormap(flipud(bamako));caxis(caxis);caxis([0 4.5e8]);end
    colorbar;ylabel(colorbar,collabel);
    hold on;
    plot(surface(1,:),surface(2,:),'k');
    if istring ==1; hold on; plot(string(:,1)/1000, string(:,2)/1000, 'w'); end 
    if itemp==1;[C,h]=contour(x/1000,y/1000,temperature,geotherms,tclr,'LineWidth',0.75);clabel(C,h,'LabelSpacing',1700,'color',tclr);end
    if ivelo==1;quiver(velo_x/1000,velo_y/1000,velo_vx,velo_vy,'k');end
    xlim([xmin xmax]);ylim([ymin ymax]);
    step=num2str(dtstep); conv=num2str(dtstep*in_vel);
    xlabel('Along profile distance (km)');ylabel('Depth (km)');
    title([titletag step ' Myr']);

    % Puts minimum and maximum values of the plotted field on the figure in
    % a textbox on the bottom left corner
    if (annote)
        if (isequal(field,'elvisco') || isequal(field,'eldensi') || isequal(field,'elstrrt')|| isequal(field,'eldissp'))
            minel = num2str(min(min(elfield)),'%.6f');
            maxel = num2str(max(max(elfield)),'%.6f');
        elseif (isequal(field,'elsigma') || isequal(field,'eltresc'))
            minel = num2str(min(min(elfield)),'%.2e');
            maxel = num2str(max(max(elfield)),'%.2e');
        end
        annotation('textbox',[0.75 0.25 0.1 0.1],'String',['max = ' maxel; 'min = ' minel],'BackGroundColor',[1 1 1],'HorizontalAlignment','right');
    end
    
    if (iout==1);nom=strcat(strcat(folder,outname,Ttag,vtag,Ztag,tag,'_',timeform));print('-dpdf','-r600','-painters',nom);end
end
