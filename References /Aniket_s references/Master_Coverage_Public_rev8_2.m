%Jan 4 2018 - Adam Schroeder - a.schroeder.me@gmail.com
%This is code for simulating how robots move in a discretized area using
%different types of biased or unbiased random walks. The walks are biased
%by a virtual 'pheromone' that robots emit at a constant rate as they move
%through the domain. Other agents can detect and are repulsed by this 
%pheromone, which diffuses and evaporates over time.
% Different metrics related to area coverage are tracked.

%Some of the options included have only been used in older versions of the
%code and may not directly work if engaged, but the code should run as
%written

%Plots are produced of the agents paths, percent of area covered over time
%and the pheromone field. Parameters displayed at the end include the
%percent of area covered, the entropy of the number of visits to each
%location, and the number of 'simple' and 'smart' pop-up threats detected

%This is the code used for "Efficient Spatial Coverage by a Robot
%Swarm based on an Ant Foraging Model and the Lévy Distribution" in
%Springer's Journal of Swarm Intelligence, 2017

%The main parameters that a user would control are:
%'mode' to choose gradient following with a constant path length, gradient
%following with Levy flight, or unbiased Levy flight
%'D' to set the diffusion rate of pheromone
%'gamma' to set the evaporation rate of pheromone
%'noise' absolute magnitude of random noise added when new path is
%calculated
%'levyalpha' parameter used to find the path length for levy flight
%'time' simulated time
%'agents' number of agents
%'bc' sets boundary condition for pheromones (dirichlet, von Neumann, or
%periodic)
%'magictorus' sets if agents are allowed to move across boundaries or not

%An FTCS numerical scheme is used for diffusion. This version of the code
%does not use a GPU for the diffusion calculation, but subsequent work has
%shown this to give a significant speed increase for this type of
%calculation, which is the most computationally intense portion of this
%code.

%Subsequent work with the code has identified discretization effects which
%are masked if the noise is high enough. If the noise is low, the agents
%can be observed traveling in nearly vertical or nearly horizontal lines.
%This is caused by how pheromone is deposited and the gradient calculated
%in any grid. 

%This cleaner Github version was based on my internal file rev8, v2

tic

clc
clear
clear plot

%display settings
update1=1; %if one, update agent position plot in real time, which of course slows the code
updateevery=1; %update agent position on plot every so many time steps (used with update1)
displayoutputs=1; %if one, output plots will be displayed
exportcount=1; %if one, cells visited at each time step and frequency will be exported to excel
exportfrequency=1; %sets how often data will be sampled (must be limited due to excel column limit)
filename = 'Periodic_MagicTorus_PointOne.xlsx'; %only needed if exportcount=1

%mode settings
mode=1; %0=gradient following with constant path length, 1=gradient following with Levy flight, 2=pure levy flight
constantlength=1; %only needed for mode=0
levyalpha=1.0; %only needed for mode=1 or mode=2 Levy flight parameter
maxpath=141; %if levy flight selects a path longer than

%general settings
deltat=.1; %time step
gamma=0.0001; %evaporation rate=percent pheromone/second
deposit=1; %deposition rate = units pheromone/second
travel=1; %travel speed = units distance/second
noise=0.1; %magnitude of noise, if gradient strong, noise will be neglible, if weak, noise will be significant
gamma_effective=gamma*deltat; %adjusted to time step
deposit_effective=deposit*deltat; %adjusted to time step
travel_effective=travel*deltat; %adjusted to time step
coverage=1.0; %program stops when this percent has been visited at least once.
time=10; %maximum time steps allowed
agents=10; %number of randomly initialized agents
magictorus=0; %if zero agents can't pass through boundaries

%diffusion settings
D=0.001; %diffusion rate
resolution=4; %units per unit length
deltax=1/resolution;
deltay=1/resolution;
bc=0; %0 means constant value - set 'edge' value in loop, 1 means constant flux of zero, 2 means periodic (agents+pher)

%check if FTCS diffusion will be stable
if D>deltax*deltax/(4*deltat);
    disp('diffusion is unstable due to deltax,deltat,diffusion constant')
    pause
end

%search area boundaries
xmin=0;
xmax=100;
ymin=0;
ymax=100;
m=2:(xmax/deltax-1); %used as indices for later diffusion calcs
n=2:(ymax/deltay-1);

%Pop Up Threat settings
popduration=10; %duration of popup threat
popdetection=5; %detection range of agent for pop-up threat
popquantity=round(time/popduration); %number of pop-up occurences

%create grid for pheromone, use pre-assigned resolution
xgv=linspace(xmin,xmax,(xmax-xmin)/deltax);%assume grid cell size of 1 unit
ygv=linspace(ymin,ymax,(ymax-ymin)/deltay);%assume grid cell size of 1 unit

%create grid for area coverage, resolution of one unit
xvv=linspace(xmin,xmax,(xmax-xmin));%assume grid cell size of 1 unit
yvv=linspace(ymin,ymax,(ymax-ymin));%assume grid cell size of 1 unit

% The following q, r, and s while loops are used if you want to loop
% through different parameter combinations, e.g. diffusion rate, evap rate,
% noise, using 'q' and 'r', or just want to perform multiple statistical runs,
% using 's'. By default, these are all set to one
q=1; %multiple set index
qmax=1; %number of sets
while q<=qmax
    
    r=1; %multiple run index
    rmax=1; %number of runs
    while r<=rmax
        
        s=1; %statistical run, no variation
        smax=1; %number of statistical runs
        while s<=smax;
            %tracks pop-up threats detected and initializes pop-up threats
            smartpopdetected=zeros(1,popquantity); %will track if detected, if one, means detected
            smartpoplocationx=xmax.*rand(1,popquantity);
            smartpoplocationy=ymax.*rand(1,popquantity);
            simplepopdetected=zeros(1,popquantity); %will track if detected, if one, means detected
            simplepoplocationx=xmax.*rand(1,popquantity);
            simplepoplocationy=ymax.*rand(1,popquantity);
            
            %these diffusion constants inside of loops need if D value
            %being varied
            alpha=D*deltat/deltax^2;
            beta=D*deltat/deltay^2;
            
            %initialize agents
            Nx=xmax.*rand(1,agents); %initialize population randomly
            Ny=xmax.*rand(1,agents); %initialize population randomly
            
            pathremaining=zeros(1,agents); %will give remaining path length
            agentheading=zeros(1,agents); %will track agent headings
            pheromone=zeros(numel(xgv),numel(ygv)); %initial matrix of known/unknown cells
            deltapheromone=zeros(numel(xgv),numel(ygv)); %will be used for changes to pheromone in each time step
            visited=zeros(numel(xvv),numel(yvv)); %once a cell has been visited once, will turn to one
            frequency=zeros(numel(xvv),numel(yvv)); %will track how many times each cell has been visited
            
            %monitors cells that have not been visited
            for k=1:numel(Nx)
                currentvisitcellx=floor(Nx(1,k))+1; %finds current cell x
                currentvisitcelly=floor(Ny(1,k))+1; %finds current cell y
                visited(currentvisitcelly,currentvisitcellx)=1;
                frequency(currentvisitcelly,currentvisitcellx)=frequency(currentvisitcelly,currentvisitcellx)+1;
            end
            t=0; %time index
            
            %agent position graphing
            if displayoutputs==1;
                c = linspace(1,10,numel(Nx));
                clf(figure(2))
                clf(figure(3))
                clf(figure(4))
                figure(2); %display position
                scatter(Nx(1,:),Ny(1,:),[],c);
                xlim([xmin xmax])
                ylim([ymin ymax])
                hold on;
            end
            
            u=1; %update index
            loop=1; %loop index
            stray=0; %stray counter
            export=1; %will track when to export data
            
            count=ones(1,time/exportfrequency); %tracks percent of locations visited at least once
            
            while t<=time
                
                %% this section of code calculates the gradient if needed, and updates agent positions
                
                %now find gradient
                Nx_prev=Nx; %will store previous positions
                Ny_prev=Ny; %will store previous positions
                for k=1:numel(Nx); %now update agent positions
                    currentcellx=floor(Nx(1,k)/deltax)+1; %finds current cell x %gives column number
                    currentcelly=floor(Ny(1,k)/deltay)+1; %finds current cell y %gives row number
                    %should never be outside of grid
                    if Nx(1,k)<=xmin||Nx(1,k)>xmax||Ny(1,k)<=ymin||Ny(1,k)>ymax %if outside of grid, move back toward center of grid
                        agentheading(k)=atan2(-(Ny(1,k)-50),-(Nx(1,k)-50)); %heading is summation of noise and reverse gradient
                        pathremaining(k)=0;
                        stray=stray+1;
                    else %if not outside of grid
                        if pathremaining(k)<=0 %if a new path needs calculated
                            if mode~=2 %if not pure levy flight, calculate gradient
                                up=currentcelly+1; %gives row number
                                down=currentcelly-1; %gives row number
                                left=currentcellx-1; %gives column number
                                right=currentcellx+1; %gives column number
                                %all set to one just to track for special circumstances
                                phero1=1;
                                phero3=1;
                                phero5=1;
                                phero7=1;
                                %look for special conditions at edge of grid
                                if up==ymax/deltay+1;
                                    if magictorus==0
                                        phero1=0;
                                    else
                                        up=1;
                                    end
                                end
                                if down<=ymin/deltay+1;
                                    if magictorus==0
                                        phero5=0;
                                    else
                                        down=ymax/deltay;
                                    end
                                end
                                if right==xmax/deltax+1;
                                    if magictorus==0
                                        phero3=0;
                                    else
                                        right=1;
                                    end
                                end
                                if left<=xmin/deltax+1;
                                    if magictorus==0
                                        phero7=0;
                                    else
                                        left=xmax/deltax;
                                    end
                                end
                                
                                %now set pheromone levels if no special conditions apply
                                if phero1~=0;
                                    phero1=pheromone(up,currentcellx);
                                end
                                
                                if phero5~=0;
                                    phero5=pheromone(down,currentcellx);
                                end
                                
                                if phero3~=0;
                                    phero3=pheromone(currentcelly,right);
                                end
                                
                                if phero7~=0;
                                    phero7=pheromone(currentcelly,left);
                                end
                                
                                %these are reverse gradient values (pheromone repulsive)
                                gradientx=(phero7-phero3)/(2*deltax);
                                gradienty=(phero5-phero1)/(2*deltay);
                                
                            else %for pure levy flight
                                gradientx=0;
                                gradienty=0;
                            end
                            
                            %these are the noise values
                            noiseangle=2*pi*rand;
                            noisex=noise*cos(noiseangle);
                            noisey=noise*sin(noiseangle);
                            agentheading(k)=atan2(gradienty+noisey,gradientx+noisex); %heading is summation of noise and reverse gradient
                            if mode==0 %non-levy flight
                                pathremaining(k)=constantlength-travel_effective;
                            else %levy flight
                                pathremaining(k)=rand^(-1/levyalpha)-travel_effective;
                                if pathremaining(k)>maxpath
                                    pathremaining(k)=maxpath;
                                end
                            end
                            
                            %if choices cause path to stray outside of grid,
                            %recalculate
                            if magictorus==0
                                while Nx(1,k)+travel_effective*cos(agentheading(k))<xmin||Nx(1,k)+travel_effective*cos(agentheading(k))>xmax||Ny(1,k)+travel_effective*sin(agentheading(k))<ymin||Ny(1,k)+travel_effective*sin(agentheading(k))>ymax
                                    noiseangle=2*pi*rand;
                                    noisex=noise*cos(noiseangle);
                                    noisey=noise*sin(noiseangle);
                                    agentheading(k)=atan2(noisey,noisex); %heading is summation of noise and reverse gradient
                                    if mode==0 %non-Levy flight
                                        pathremaining(k)=constantlength-travel_effective;
                                    else %Levy flight
                                        pathremaining(k)=rand^(-1/levyalpha)-travel_effective;
                                        if pathremaining(k)>maxpath
                                            pathremaining(k)=maxpath;
                                        end
                                    end
                                end %end while new path generated moves outside grid
                            end
                            %end %end if a new path needs calculated
                            %end %end if outside grid
                        else %if a new path doesn't need calculated
                            if Nx(1,k)+travel_effective*cos(agentheading(k))<xmin||Nx(1,k)+travel_effective*cos(agentheading(k))>xmax||Ny(1,k)+travel_effective*sin(agentheading(k))<ymin||Ny(1,k)+travel_effective*sin(agentheading(k))>ymax
                                if magictorus==0
                                    while Nx(1,k)+travel_effective*cos(agentheading(k))<xmin||Nx(1,k)+travel_effective*cos(agentheading(k))>xmax||Ny(1,k)+travel_effective*sin(agentheading(k))<ymin||Ny(1,k)+travel_effective*sin(agentheading(k))>ymax
                                        noiseangle=2*pi*rand;
                                        noisex=noise*cos(noiseangle);
                                        noisey=noise*sin(noiseangle);
                                        agentheading(k)=atan2(noisey,noisex); %heading is summation of noise and reverse gradient
                                        if mode==0% non-Levy flight
                                            pathremaining(k)=constantlength-travel_effective;
                                        else %Levy flight
                                            pathremaining(k)=rand^(-1/levyalpha)-travel_effective;
                                            if pathremaining(k)>maxpath
                                                pathremaining(k)=maxpath;
                                            end
                                        end %end mode determination
                                    end %end while
                                    Nx(1,k)=Nx(1,k)+travel_effective*cos(agentheading(k)); 
                                    Ny(1,k)=Ny(1,k)+travel_effective*sin(agentheading(k)); 
                                else %if agents can cross boundaries
                                    Nx(1,k)=Nx(1,k)+travel_effective*cos(agentheading(k)); 
                                    Ny(1,k)=Ny(1,k)+travel_effective*sin(agentheading(k)); 
                                end
                            else %continue on path
                                Nx(1,k)=Nx(1,k)+travel_effective*cos(agentheading(k));
                                Ny(1,k)=Ny(1,k)+travel_effective*sin(agentheading(k));
                                pathremaining(k)=pathremaining(k)-travel_effective; %decrease length of path remaining
                            end %if current path will move out of grid
                        end %if a new path needs created
                    end %if outside of grid (shouldn't be needed now)
                end %end for agent loop
                
                deltapheromone=zeros(numel(xgv),numel(ygv));
                
                %% monitors cells that have not been visited and will also add               
                for k=1:numel(Nx)  
                    
                    %% this deposits pheromone in cells that have been traveled through
                    endcellx=floor(Nx(1,k)/deltax)+1; %finds end cell x %gives column
                    endcelly=floor(Ny(1,k)/deltay)+1; %finds end cell y %gives row
                    startcellx=floor(Nx_prev(1,k)/deltax)+1; %finds start cell x %gives column
                    startcelly=floor(Ny_prev(1,k)/deltay)+1; %finds start cell y %gives row
                    currentcellx=startcellx; %will increment through entire path %gives column
                    currentcelly=startcelly; %will increment through entire path %gives row
                    currentx=Nx_prev(1,k); %will increment through entire path
                    currenty=Ny_prev(1,k); %will increment through entire path
                    visitedcells=abs(endcellx-startcellx)+abs(endcelly-startcelly);
                    temp0=0; %will track how many visited cells have been found
                    slope=(Ny(1,k)-Ny_prev(1,k))/(Nx(1,k)-Nx_prev(1,k));
                    heading=atan2(Ny(1,k)-Ny_prev(1,k),Nx(1,k)-Nx_prev(1,k));
                    
                    while temp0<visitedcells&&mode~=2;
                        if heading>=0&&heading<=pi/2 %first quadrant, check right boundary, slope +
                            BR=(currentcelly-1)*deltay; %bottom right corner
                            TR=BR+deltay; %top right corner
                            IR=(currentcellx*deltax-currentx)*slope+currenty; %y intercept = run*slope (which equals rise) + current y
                            if IR>=BR&&IR<=TR %advance right
                                currentcellx=currentcellx+1; %x changes but y stays same
                                if currentcellx>xmax/deltax
                                    currentcellx=currentcellx-xmax/deltax;
                                end
                                currentx=(currentcellx-1)*deltax;
                                currenty=IR;
                            else %advance up
                                currentcelly=currentcelly+1; %y changes but x stays same
                                if currentcelly>ymax/deltay
                                    currentcelly=currentcelly-ymax/deltay;
                                end
                            end
                            
                            
                        else if heading>pi/2 %quadrant 2,check left boundary, slope -
                                BL=(currentcelly-1)*deltay; %bottom right corner
                                TL=BL+deltay; %top right corner
                                IL=((currentcellx-1)*deltax-currentx)*slope+currenty; %y intercept = run*slope (which equals rise) + current y
                                if IL>=BL&&IL<=TL %advance left
                                    if currentcellx<1
                                        currentcellx=currentcellx+xmax/deltax;
                                    end
                                    currentcellx=currentcellx-1; %x changes but y stays same
                                    currentx=(currentcellx)*deltax;
                                    currenty=IL;
                                else %advance up
                                    currentcelly=currentcelly+1; %y changes but x stays same
                                    if currentcelly>ymax/deltay
                                        currentcelly=currentcelly-ymax/deltay;
                                    end
                                end
                                
                            else if heading<=-pi/2 %check quadrant 3, left boundary, slope +
                                    BL=(currentcelly-1)*deltay; %bottom right corner
                                    TL=BL+deltay; %top right corner
                                    IL=((currentcellx-1)*deltax-currentx)*slope+currenty; %y intercept = run*slope (which equals rise) + current y
                                    if IL>=BL&&IL<=TL %advance left
                                        currentcellx=currentcellx-1; %x changes but y stays same
                                        if currentcellx<1
                                            currentcellx=currentcellx+xmax/deltax;
                                        end
                                        currentx=(currentcellx)*deltax;
                                        currenty=IL;
                                    else %advance down
                                        currentcelly=currentcelly-1; %y changes but x stays same
                                        if currentcelly<1
                                            currentcelly=currentcelly+ymax/deltay;
                                        end
                                    end
                                    
                                else %must be quadrant 4, right boundary, slope -
                                    BR=(currentcelly-1)*deltay; %bottom right corner
                                    TR=BR+deltay; %top right corner
                                    IR=(currentcellx*deltax-currentx)*slope+currenty; %y intercept = run*slope (which equals rise) + current y
                                    if IR>=BR&&IR<=TR %advance right
                                        currentcellx=currentcellx+1; %x changes but y stays same
                                        if currentcellx>xmax/deltax
                                            currentcellx=currentcellx-xmax/deltax;
                                        end
                                        currentx=(currentcellx-1)*deltax;
                                        currenty=IR;
                                    else %advance down
                                        currentcelly=currentcelly-1; %y changes but x stays same
                                        if currentcelly<1
                                            currentcelly=currentcelly+ymax/deltay;
                                        end
                                    end
                                end% end quadrant 3 and 4
                            end %end quadrant 2
                        end %end quadrant 1
                        
                        if currentcellx*deltax<=xmin||currentcellx*deltax>xmax||currentcelly*deltay<=ymin||currentcelly*deltay>ymax
                        else
                            deltapheromone(currentcelly,currentcellx)=deltapheromone(currentcelly,currentcellx)+deposit_effective/(visitedcells+1)/deltax; %add pheromone to visited cell
                        end
                        
                        temp0=temp0+1;
                    end %end while
                    deltapheromone(startcelly,startcellx)=deltapheromone(startcelly,startcellx)+deposit_effective/(visitedcells+1)/deltax;
                    
                    %% this tracks visits in the passed-through cells
                    endvisitcellx=floor(Nx(1,k))+1; %finds end cell x %gives column
                    endvisitcelly=floor(Ny(1,k))+1; %finds end cell y %gives row
                    startvisitcellx=floor(Nx_prev(1,k))+1; %finds start cell x %gives column
                    startvisitcelly=floor(Ny_prev(1,k))+1; %finds start cell y %gives row
                    currentvisitcellx=startvisitcellx; %will increment through entire path %gives column
                    currentvisitcelly=startvisitcelly; %will increment through entire path %gives row
                    currentvisitx=Nx_prev(1,k); %will increment through entire path
                    currentvisity=Ny_prev(1,k); %will increment through entire path
                    visitedvisitcells=abs(endvisitcellx-startvisitcellx)+abs(endvisitcelly-startvisitcelly);
                    temp1=0; %will track how many visited cells have been found
                    
                    
                    while temp1<visitedvisitcells;
                        if heading>=0&&heading<=pi/2 %first quadrant, check right boundary, slope +
                            BR=(currentvisitcelly-1); %bottom right corner
                            TR=BR+1; %top right corner
                            IR=(currentvisitcellx-currentvisitx)*slope+currentvisity; %y intercept = run*slope (which equals rise) + current y
                            if IR>=BR&&IR<=TR %advance right
                                currentvisitcellx=currentvisitcellx+1; %x changes but y stays same
                                if currentvisitcellx>xmax
                                    currentvisitcellx=currentvisitcellx-xmax;
                                end
                                currentvisitx=(currentvisitcellx-1);
                                currentvisity=IR;
                            else %advance up
                                currentvisitcelly=currentvisitcelly+1; %y changes but x stays same
                                if currentvisitcelly>ymax
                                    currentvisitcelly=currentvisitcelly-ymax;
                                end
                            end
                            
                            
                        else if heading>pi/2 %quadrant 2,check left boundary, slope -
                                BL=(currentvisitcelly-1); %bottom right corner
                                TL=BL+1; %top right corner
                                IL=((currentvisitcellx-1)-currentvisitx)*slope+currentvisity; %y intercept = run*slope (which equals rise) + current y
                                if IL>=BL&&IL<=TL %advance left
                                    currentvisitcellx=currentvisitcellx-1; %x changes but y stays same
                                    if currentvisitcellx<1
                                        currentvisitcellx=currentvisitcellx+xmax;
                                    end
                                    currentvisitx=(currentvisitcellx);
                                    currentvisity=IL;
                                else %advance up
                                    currentvisitcelly=currentvisitcelly+1; %y changes but x stays same
                                    if currentvisitcelly>ymax
                                        currentvisitcelly=currentvisitcelly-ymax;
                                    end
                                end
                                
                            else if heading<=-pi/2 %check quadrant 3, left boundary, slope +
                                    BL=(currentvisitcelly-1); %bottom right corner
                                    TL=BL+1; %top right corner
                                    IL=((currentvisitcellx-1)-currentvisitx)*slope+currentvisity; %y intercept = run*slope (which equals rise) + current y
                                    if IL>=BL&&IL<=TL %advance left
                                        currentvisitcellx=currentvisitcellx-1; %x changes but y stays same
                                        if currentvisitcellx<1
                                            currentvisitcellx=currentvisitcellx+xmax;
                                        end
                                        currentvisitx=(currentvisitcellx);
                                        currentvisity=IL;
                                    else %advance down
                                        currentvisitcelly=currentvisitcelly-1; %y changes but x stays same
                                        if currentvisitcelly<1
                                            currentvisitcelly=currentvisitcelly+ymax;
                                        end
                                    end
                                    
                                else %must be quadrant 4, right boundary, slope -
                                    BR=(currentvisitcelly-1); %bottom right corner
                                    TR=BR+1; %top right corner
                                    IR=(currentvisitcellx-currentvisitx)*slope+currentvisity; %y intercept = run*slope (which equals rise) + current y
                                    if IR>=BR&&IR<=TR %advance right
                                        currentvisitcellx=currentvisitcellx+1; %x changes but y stays same
                                        if currentvisitcellx>xmax
                                            currentvisitcellx=currentvisitcellx-xmax;
                                        end
                                        currentvisitx=(currentvisitcellx-1);
                                        currentvisity=IR;
                                    else %advance down
                                        currentvisitcelly=currentvisitcelly-1; %y changes but x stays same
                                        if currentvisitcelly<1
                                            currentvisitcelly=currentvisitcelly+ymax;
                                        end
                                    end
                                end% end quadrant 3 and 4
                            end %end quadrant 2
                        end %end quadrant 1
                        
                        if currentvisitcellx<=xmin||currentvisitcellx>xmax||currentvisitcelly<=ymin||currentvisitcelly>ymax
                        else
                            visited(currentvisitcelly,currentvisitcellx)=1;
                            if currentvisitcellx~=startvisitcellx||currentvisitcelly~=startvisitcelly %only adds to frequency if not already in cell
                                frequency(currentvisitcelly,currentvisitcellx)=frequency(currentvisitcelly,currentvisitcellx)+1;
                            end
                        end
                        
                        temp1=temp1+1;
                    end %end while
                    
                    %% if robots allowed to move through boundaries(magictorus option)
                    % they may end up outside of domain limits and their position
                    %needs reset to the opposite boundary
                    if Nx(1,k)>xmax
                        Nx(1,k)=Nx(1,k)-xmax;
                    else if Nx(1,k)<0
                            Nx(1,k)=Nx(1,k)+xmax;
                        end
                    end
                    if Ny(1,k)>ymax
                        Ny(1,k)=Ny(1,k)-ymax;
                    else if Ny(1,k)<0
                            Ny(1,k)=Ny(1,k)+ymax;
                        end
                    end
                    
                    %% this section of code checks if pop-up threats have been detected
                    temp5=min(floor(t/popduration+1),time/popduration); %finds which pop-up threat is active
                    %%if new pop-up threat detected (when remainder of t/popduration=0) just because agent
                    %%there initially, regenerate point
                    while rem(t,popduration)==0&&((Nx(k)-smartpoplocationx(temp5))^2+(Ny(k)-smartpoplocationy(temp5))^2)^0.5<=popdetection
                        smartpoplocationx(temp5)=xmax*rand(1);
                        smartpoplocationy(temp5)=ymax*rand(1);
                    end
                    if smartpopdetected(temp5)==0&&((Nx(k)-smartpoplocationx(temp5))^2+(Ny(k)-smartpoplocationy(temp5))^2)^0.5<=popdetection
                        smartpopdetected(temp5)=1;
                    end
                    if simplepopdetected(temp5)==0&&((Nx(k)-simplepoplocationx(temp5))^2+(Ny(k)-simplepoplocationy(temp5))^2)^0.5<=popdetection
                        simplepopdetected(temp5)=1;
                    end
                end %end for loop
                
                %% this section of code tracks the portion of cells that have been visited, and associated entropy and standard deviation of visitation freqency
                t=t+deltat; %moved from below
                if export==exportfrequency/deltat
                    count(1,loop)=sum(sum(visited));
                    e=numel(unique(frequency));
                    clear freqvalues;
                    freqvalues=unique(frequency);
                    temp3=0;
                    temp4=0;
                    avgfreq=1/(numel(xvv)*numel(yvv))*numel(Nx);
                    for i=1:e
                        temp1=sum(sum(frequency==freqvalues(i,1))); %number of cells with that many visits
                        temp2=temp1/(numel(xvv)*numel(yvv)); %this is probability of this frequency
                        temp3=temp3+temp2*log2(temp2); %this finds the entropy component
                        temp4=temp4+temp1*(freqvalues(i,1)/t-avgfreq)^2; %this finds the sum of all deviations
                    end
                    entropy(1,loop)=-temp3;
                    deviation(1,loop)=temp4/(numel(xvv)*numel(yvv));
                    export=0; %reset export timer
                    
                    loop=loop+1;
                end
                export=export+1; %increment export timer
                
                %% this section of code applies pheromone diffision (applied differently depending on boundary condition)
                if mode~=2; %not needed for pure levy
                    %deltapheromone=deltapheromone+pheromone.*(-gamma_effective); %should be negative because evaporation, trying to increase speed
                    pheromonenew=zeros(numel(xgv),numel(ygv));
                    
                    pheromonenew(m,n)=alpha*(pheromone(m+1,n)+pheromone(m-1,n))+beta*(pheromone(m,n+1)+pheromone(m,n-1))+(1-2*alpha-2*beta)*pheromone(m,n);
                    if bc==0 %use constant value of zero
                        %edge=mean(mean(pheromone));
                        edge=0;
                        pheromonenew(1,n)=alpha*(pheromone(1+1,n)+edge)+beta*(pheromone(1,n+1)+pheromone(1,n-1))+(1-2*alpha-2*beta)*pheromone(1,n);
                        pheromonenew(xmax/deltax,n)=alpha*(edge+pheromone(xmax/deltax-1,n))+beta*(pheromone(xmax/deltax,n+1)+pheromone(xmax/deltax,n-1))+(1-2*alpha-2*beta)*pheromone(xmax/deltax,n);
                        pheromonenew(m,1)=alpha*(pheromone(m+1,1)+pheromone(m-1,1))+beta*(pheromone(m,1+1)+edge)+(1-2*alpha-2*beta)*pheromone(m,1);
                        pheromonenew(m,ymax/deltay)=alpha*(pheromone(m+1,ymax/deltay)+pheromone(m-1,ymax/deltay))+beta*(edge+pheromone(m,ymax/deltay-1))+(1-2*alpha-2*beta)*pheromone(m,ymax/deltay);
                        
                        %solve for special case of corners
                        pheromonenew(1,1)=alpha*(pheromone(2,1)+edge)+beta*(pheromone(1,2)+edge)+(1-2*alpha-2*beta)*pheromone(1,1);
                        pheromonenew(1,ymax/deltay)=alpha*(pheromone(2,ymax/deltay)+edge)+beta*(edge+pheromone(1,ymax/deltay-1))+(1-2*alpha-2*beta)*pheromone(1,ymax/deltay);
                        pheromonenew(xmax/deltax,ymax/deltay)=alpha*(edge+pheromone(xmax/deltax-1,ymax/deltay))+beta*(edge+pheromone(xmax/deltax,ymax/deltay-1))+(1-2*alpha-2*beta)*pheromone(xmax/deltax,ymax/deltay);
                        pheromonenew(xmax/deltax,1)=alpha*(edge+pheromone(xmax/deltax-1,1))+beta*(pheromone(xmax/deltax,1+1)+edge)+(1-2*alpha-2*beta)*pheromone(xmax/deltax,1);
                    else if bc==1 %use constant flux of zero
                            pheromonenew(1,n)=pheromone(1,n)+alpha*((pheromone(1+1,n)-pheromone(1,n))-(0))+beta*((pheromone(1,n+1)-pheromone(1,n))-(pheromone(1,n)-pheromone(1,n-1)));
                            pheromonenew(xmax/deltax,n)=pheromone(xmax/deltax,n)+alpha*((0)-(pheromone(xmax/deltax,n)-pheromone(xmax/deltax-1,n)))+beta*((pheromone(xmax/deltax,n+1)-pheromone(xmax/deltax,n))-(pheromone(xmax/deltax,n)-pheromone(xmax/deltax,n-1)));
                            pheromonenew(m,1)=pheromone(m,1)+alpha*((pheromone(m+1,1)-pheromone(m,1))-(pheromone(m,1)-pheromone(m-1,1)))+beta*((pheromone(m,1+1)-pheromone(m,1))-(0));
                            pheromonenew(m,ymax/deltay)=pheromone(m,ymax/deltay)+alpha*((pheromone(m+1,ymax/deltay)-pheromone(m,ymax/deltay))-(pheromone(m,ymax/deltay)-pheromone(m-1,ymax/deltay)))+beta*((0)-(pheromone(m,ymax/deltay)-pheromone(m,ymax/deltay-1)));
                            
                            %solve for special case of corners
                            pheromonenew(1,1)=pheromone(1,1)+alpha*((pheromone(1+1,1)-pheromone(1,1))-(0))+beta*((pheromone(1,1+1)-pheromone(1,1))-(0));
                            pheromonenew(1,ymax/deltay)=pheromone(1,ymax/deltay)+alpha*((pheromone(1+1,ymax/deltay)-pheromone(1,ymax/deltay))-(0))+beta*((0)-(pheromone(1,ymax/deltay)-pheromone(1,ymax/deltay-1)));
                            pheromonenew(xmax/deltax,ymax/deltay)=pheromone(xmax/deltax,ymax/deltay)+alpha*((0)-(pheromone(xmax/deltax,ymax/deltay)-pheromone(xmax/deltax-1,ymax/deltay)))+beta*((0)-(pheromone(xmax/deltax,ymax/deltay)-pheromone(xmax/deltax,ymax/deltay-1)));
                            pheromonenew(xmax/deltax,1)=pheromone(xmax/deltax,1)+alpha*((0)-(pheromone(xmax/deltax,1)-pheromone(xmax/deltax-1,1)))+beta*((pheromone(xmax/deltax,1+1)-pheromone(xmax/deltax,1))-(0));
                            
                        else %use perioidic boundary conditions
                            pheromonenew(1,n)=alpha*(pheromone(1+1,n)+pheromone(xmax/deltax,n))+beta*(pheromone(1,n+1)+pheromone(1,n-1))+(1-2*alpha-2*beta)*pheromone(1,n);
                            pheromonenew(xmax/deltax,n)=alpha*(pheromone(1,n)+pheromone(xmax/deltax-1,n))+beta*(pheromone(xmax/deltax,n+1)+pheromone(xmax/deltax,n-1))+(1-2*alpha-2*beta)*pheromone(xmax/deltax,n);
                            pheromonenew(m,1)=alpha*(pheromone(m+1,1)+pheromone(m-1,1))+beta*(pheromone(m,1+1)+pheromone(m,ymax/deltay))+(1-2*alpha-2*beta)*pheromone(m,1);
                            pheromonenew(m,ymax/deltay)=alpha*(pheromone(m+1,ymax/deltay)+pheromone(m-1,ymax/deltay))+beta*(pheromone(m,1)+pheromone(m,ymax/deltay-1))+(1-2*alpha-2*beta)*pheromone(m,ymax/deltay);
                            
                            %solve for special case of corners
                            pheromonenew(1,1)=alpha*(pheromone(2,1)+pheromone(xmax/deltax,1))+beta*(pheromone(1,2)+pheromone(1,ymax/deltay))+(1-2*alpha-2*beta)*pheromone(1,1);
                            pheromonenew(1,ymax/deltay)=alpha*(pheromone(2,ymax/deltay)+pheromone(xmax/deltax,ymax/deltay))+beta*(pheromone(1,1)+pheromone(1,ymax/deltay-1))+(1-2*alpha-2*beta)*pheromone(1,ymax/deltay);
                            pheromonenew(xmax/deltax,ymax/deltay)=alpha*(pheromone(1,ymax/deltay)+pheromone(xmax/deltax-1,ymax/deltay))+beta*(pheromone(xmax/deltax,1)+pheromone(xmax/deltax,ymax/deltay-1))+(1-2*alpha-2*beta)*pheromone(xmax/deltax,ymax/deltay);
                            pheromonenew(xmax/deltax,1)=alpha*(pheromone(1,1)+pheromone(xmax/deltax-1,1))+beta*(pheromone(xmax/deltax,1+1)+pheromone(xmax/deltax,ymax/deltay))+(1-2*alpha-2*beta)*pheromone(xmax/deltax,1);
                            
                        end
                    end
                    
                    pheromone=pheromonenew+(deltapheromone+pheromone.*(-gamma_effective));
                end %end if
                
                
                %% update plots while the simulation is running
                if displayoutputs==1 %plot is updated at the end of the run
                    figure(2); %display position
                    scatter(Nx(1,:),Ny(1,:),[],c);
                    title('agent position history')
                    
                    if update1==1||u*deltat==updateevery; %with these options, plot is updated in real time
                        h(1)=scatter(Nx(1,:),Ny(1,:),[],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'SizeData',75);
                        
                        figure(4) %display pheromone distribution, replace figure 1
                        image(pheromone,'CDataMapping','scaled')
                        set(gca,'Ydir','Normal')
                        title('pheromone field')
                        colorbar;
                        axis off
                        
                        delete(h(1));
                        drawnow;
                        u=0;
                    end
                end
                
                u=u+1;
            end
            
            %% export data and display plots
            %plot naming for export
            rowname=strcat('A',num2str(s));
            sheetname1=strcat('cover',',q-',num2str(q),',r-',num2str(r));
            sheetname2=strcat('entropy',',q-',num2str(q),',r-',num2str(r));
            sheetname3=strcat('dev.',',q-',num2str(q),',r-',num2str(r));
            sheetname4=strcat('sipop',',q-',num2str(q),',r-',num2str(r));
            sheetname5=strcat('smpop',',q-',num2str(q),',r-',num2str(r));
            
            if displayoutputs==1 %plots are updated or generated at the end of the un
                figure(2); %display position
                scatter(Nx(1,:),Ny(1,:),[],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
                title('agent position history')
                xlim([xmin xmax])
                ylim([ymin ymax])
                hold on;
                
                figure(3); %display percent of area explored
                plot(1:loop-1,(count(:)/(numel(xvv)*numel(yvv))),'b-');
                title('area exploration rate')
                xlabel('time,t')
                ylabel('% of area explored')
                
                figure(4) %display pheromone distribution
                image(pheromone,'CDataMapping','scaled')
                title('pheromone field')
                set(gca,'Ydir','Normal')
                colorbar;
                axis off
                
                hold on;
            end
            
            if exportcount==1;
                warning('off','MATLAB:xlswrite:AddSheet')
                xlswrite(filename,count,sheetname1,rowname)
                xlswrite(filename,entropy,sheetname2,rowname)
                xlswrite(filename,deviation,sheetname3,rowname)
                xlswrite(filename,simplepopdetected,sheetname4,rowname)
                xlswrite(filename,smartpopdetected,sheetname5,rowname)
            end
            
            s=s+1;
        end %while statistical run
        
        r=r+1;
    end %while multiple run
    
    q=q+1;
end %while set

toc

sum(sum(count))/10000
sum(entropy)
sum(smartpopdetected)
sum(simplepopdetected)







