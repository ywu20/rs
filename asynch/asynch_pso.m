%% Particle Swarm Optimization Simulation

%% Initialization
function [success time] = pso(initial_nest)

%initialize pso.
iterations =100;
inertia = 3600; %speed

%how birds respond
correction_factor = 2.0;
swarms = size(initial_nest,1);

% ---- initial swarm position -----
swarm=zeros(swarms,9);
% 1,2 calculated positions.
% 3,4 actual position
% 5,6 random factors that determine next position
% 7 concentration
% 10 distance
for i = 1 : swarms
swarm(i, 1:2) =initial_nest(i,:) ; %swarm size
swarm(i, 3:4) =initial_nest(i,:);
swarm(i,8:9) = initial_nest(i,:);
swarm(i,7)=fobj(swarm(i,1:2));
swarm(i,10) = 0; %total distance travelled
end

[temp gbest] = min(swarm(:,7));
swarm(:, 5) = 0;          % initial velocity
swarm(:, 6) = 0;          % initial velocity
%% Iterations
time = 0; maxd = 0; T = 0;
gbest = min(i,7);
tStart=tic;
tEnd=0;
%while (tEnd<0.15) % experiment based on cpu time
%while T < 500 * size(initial_nest) % experiment based on generation
%while sum(abs(swarm(gbest,3:4)-[10 10])<= [2 2]) ~=2 %experiment based on whether it's found or not
 a= sqrt((swarm(gbest,2)-10) * (swarm(gbest,2)-10) + (swarm(gbest,3)-10)* (swarm(gbest,3)-10));
while a<=1
    T=T+1;
    pso_first = 1;
    min_d = swarm(1,10);
    for i=2:swarms
        if min_d > swarm (i,10)
            min_d = swarm (i,10);
            pso_first = i;
        end
    end

    % generating random factors
 swarm(pso_first, 5) = rand*inertia*swarm(pso_first, 5) + correction_factor*rand*(swarm(pso_first, 3)... 
    - swarm(pso_first, 1)) + correction_factor*rand*(swarm(gbest, 3) - swarm(pso_first, 1));   % u velocity parameters, updated best velocity of swarm. 
    swarm(pso_first, 6) = rand*inertia*swarm(pso_first, 6) + correction_factor*rand*(swarm(pso_first, 4)...
     - swarm(pso_first, 2)) + correction_factor*rand*(swarm(gbest, 4) - swarm(pso_first, 2));   % v velocity parameters
 
    %-- updating position of Swarms ---
    swarm(pso_first, 1) = swarm(pso_first, 3) + swarm(pso_first, 5)/1.2  ;  %update u position column 1&2 stores previous position.
    swarm(pso_first, 2) = swarm(pso_first, 4) + swarm(pso_first, 6)/1.2 ;    %update v position
    

    % if out of bound go to a random place
        if swarm (pso_first,1) < 1 
            %swarm (pso_first,1) = 1;
            swarm (pso_first,1) = randi(50);
        end
        if swarm(pso_first,1) > 50 
            %swarm(pso_first,1)=50;
            swarm(pso_first,1)=randi(50);
        end
        if swarm (pso_first,2) <1 
            %swarm (pso_first,2) = 1;
            swarm (pso_first,2) = randi(50);
        end
        if swarm(pso_first,2) > 50
            %swarm(pso_first,2)=50;
            swarm(pso_first,2)=randi(50);
        end    
        if (swarm(pso_first, 5) == 0 && swarm(pso_first,6) ==0)
            swarm(pso_first,1)=randi(50);
            swarm(pso_first,2)=randi(50);
        end

        value = fobj(swarm(pso_first,1:2));          %Objective function has infinite solution. By taking 500 swarms we have consider the 5000 best solutions.  
        d = sqrt((swarm(pso_first,3)-swarm(pso_first,1)) * (swarm(pso_first,3)-swarm(pso_first,1)) + (swarm(pso_first,4)-swarm(pso_first,2)) * (swarm(pso_first,4)-swarm(pso_first,2)));
        %swarm(pso_first,8) = swarm(pso_first,1); swarm (pso_first,9) = swarm(pso_first,2);
        swarm (pso_first,10) =  swarm (pso_first,10) + d;
   
        if value < swarm(pso_first, 7)           % Always True
            swarm(pso_first, 3) = swarm(pso_first, 1);    % update best pospso_firsttpso_firston of u, updated position.
            swarm(pso_first, 4) = swarm(pso_first, 2);    % update best postions of v,
            swarm(pso_first, 7) = value;          % best updated minimum value Pbest (the best in a iteration)   
        else
            swarm (pso_first,10) =  swarm (pso_first,10) + d;
        end
    [temp, gbest] = min(swarm(:, 7));        % gbest position, global best
    time=time+maxd/inertia;
    %--- updating velocity vectors
   
 a= sqrt((gbest(1,1)-10) * (gbest(1,1)-10)) + ((gbest(1,2)-10)* (gbest(1,2)-10));
 %if a(1,1)<= 2  && a(1,2)<= 2 
 %   success = 1;
 %   break
 %else 
 %   success = 0;
 %end
     %% Plotting the swarm
%figure(1)
 %    clf    
  %  plot(swarm(:, 1), swarm(:, 2), 'x')   % drawing swarm movements
   % axis([0 50 0 50])
 %pause(.1)

 tEnd=toc(tStart);
end
success = 1;
time = max(swarm(:,10))/inertia;
end