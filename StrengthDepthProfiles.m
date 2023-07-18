%Plot theoretical strength vs depth for our models to compare with modeled
%stresses.

%Solve equation: 2nd inv of dev stress tensor = 1/(2^((n-1)/n)) * strrt^(1/n) * 1/(A^(1/n) * d^(m/n))*e^(E+VP/nRT). see Taras Gerya pg 74,75

field='eltresc';outname='stress_';    collabel='differential stress';         tclr='r';         % differential stress
folder='/Users/nataliehummel/Desktop/matlab/diffDisl5_outs/';  titletag = 'Initial collision test: '; in_vel = 50; fignum=2; %in mm/yr or km/Myr

%domain size
n=161; nel=n-1;
m=587; mel=m-1;

t=2400;
   if     t<10;    timeform=['00000' num2str(t)];
    elseif t<100;   timeform=['0000'  num2str(t)];
    elseif t<1000;  timeform=['000'   num2str(t)]; 
    elseif t<10000; timeform=['00'    num2str(t)];  
    elseif t<100000;timeform=['0'     num2str(t)];          
   end
    
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
 
 file=strcat(folder,'nodtemp',timeform,'.out');
        temp=importdata(file," ",1);
        data=getfield(temp,'data');
        clear temp
        x=reshape(data(:,1),[n,m]);
        y=reshape(data(:,2),[n,m]);
        temperature=reshape(data(:,3),[n,m]);
        clear data
 
%load strain rate
field='elstrrt';outname='strrt_';     collabel='log_1_0 (strain rate)';       tclr='w';
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
    elstrrt=reshape(data(:,3),[mel,nel]);elstrrt=elstrrt';
    clear data
    
%strength vs depth plots
%
col = 400; %302 is trench
C=2000000;
phicrust = 2;
pcrust=tan(phicrust*pi/180);
Cl=20000000;
philith= 10;
philithstr = 20;
nl=tan(philith*pi/180);
nlstr=tan(philithstr*pi/180);
HO=1000;
R=8.314;
colX=num2str(elx(1,col));

% plot strength in subducted portion of the slab
%
if col < 300
    file=strcat(folder,'strings',timeform,'.out');
        temp=importdata(file," ",1);
        stringPoints=getfield(temp,'data');
        clear temp
     gs=[];
     fs=[];
     is=[];
     elfield2 = elfield;
     stringMark=440;
     for k=stringMark %set location of cross section
                point = [stringPoints(k,1),stringPoints(k,2)]; %[x,y]
                nextPoint = [stringPoints(k+1,1),stringPoints(k+1,2)];
                prevPoint = [stringPoints(k-1,1),stringPoints(k-1,2)];
                dip = atan((prevPoint(1,2)-nextPoint(1,2))/(prevPoint(1,1)-nextPoint(1,1)));
                crossLine = point;
                crossMoments = zeros(1,77);
                for i=1:77    %string is 3 km down, line continues 77 km down through plate.
                    if ((prevPoint(1,1)-nextPoint(1,1))<0)
                        newPoint = point+[i*1000*sin(dip), -i*1000*cos(dip)]; 
                    else
                        newPoint = point+[i*1000*sin(-1*dip), i*1000*cos(-1*dip)];
                    end
                    xcoord=elx(1); 
                    ycoord=ely(1);
                    f=1;
                    g=1;
                    for j=1:(m-1) 
                        if (abs(xcoord-newPoint(1))>abs(newPoint(1)-elx(1,j))); xcoord=elx(1,j); f = f+1; end
                    end
                    for j=1:(n-1)
                          if abs(ycoord-newPoint(2))>abs(ely(j,1)-newPoint(2)); ycoord=ely(j,1);g = g+1; end              
                    end
                    elfield2(g,f) = 10000000000; %color point (g,f) on elfield map so that we can plot the area we are calculating over 
                    gs = [gs,g];
                    fs=[fs,f];
                    is = [is,i];
                end
     end
     elfieldSample = [];
     T=[];
     P=[];
     strrt=[];
     for i =1:77
        elfieldSample = [elfieldSample, elfield(gs(i),fs(i))];
        thisP = 9.8 .* 3200 .* -1.*ely(gs(i), fs(i));
        T=[T,273+temperature(gs(i), fs(i))];
        strrt = [strrt, elstrrt(gs(i), fs(i))];
        P = [P,thisP];
     end
     
    

    plot(0.5.*elfieldSample, is.*-1); %differential stress/2 vs depth
    hold on;
    plot(((Cl+nl.*P)./(sqrt(nl^2+1)-nl)), is.*-1);%eq for stress vs depth plastic failure envelope, weakened
    hold on;
    plot(((Cl+nlstr.*P)./(sqrt(nlstr^2+1)-nlstr)), is.*-1);%eq for stress vs depth plastic failure envelope, unweakened
    hold on;
    %dislocation creep law:
    plot(min(10^10,(0.5*(strrt.*(5.33.*(10^-19).*(HO^1.2).*exp((-480000-(P.*11.*(10^-6)))./(R.*T))).^-1).^(1/3.5))), is.*-1); %the minimum is so thatthe bounds of the plot are not dragged out to 10^23 on the x axis
    hold on;
    %diffusion creep law:  
    plot(min(10^10,(0.5*(strrt.*(1.5.*(10^-18).*(0.005^-3).*(HO).*exp((-335000-(P.*4.*(10^-6)))./(R.*T))).^-1))), is.*-1);
    xlabel('Strength (differential stress/2)');ylabel('Depth (km)');
    %hold on;
    annotation('textbox',[0.75 0.25 0.1 0.1],'String',['String point = ' stringMark], 'BackGroundColor',[1 1 1],'HorizontalAlignment','right');
    annotation('textbox',[0.75 0.53 0.1 0.1],'String',['model time step = ' num2str(t)], 'BackGroundColor',[1 1 1],'HorizontalAlignment','right');
    legend('true strength','dislocation creep', 'diffusion creep');
    ylim([-120 0])
end
