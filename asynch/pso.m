%% Particle Swarm Optimization Simulation

%next steps:
    %early stop when found
    %time things
    %draw diagram
% Find minimum of  the objective function
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
for i = 1 : swarms
swarm(i, 1:2) =initial_nest(i,:) ; %swarm size
swarm(i, 3:4) =initial_nest(i,:);
swarm(i,8:9) = initial_nest(i,:);
swarm(i,7)=fobj(swarm(i,1:2));
end
[temp gbest] = min(swarm(:,7));
swarm(:, 5) = 0;          % initial velocity
swarm(:, 6) = 0;          % initial velocity
%% Iterations
time = 0; maxd = 0; T = 0;
  tStart=tic;
        tEnd=0;
        
%while tEnd<0.15
%while T<500
while sum(abs(swarm(gbest,3:4)-[10 10])<= [2 2]) ~=2
    T=T+1;
    %-- position of Swarms ---
    for i = 1 : swarms
        
        swarm(i, 1) = swarm(i, 1) + swarm(i, 5)/1.2  ;  %update u position column 1&2 stores previous position.
        swarm(i, 2) = swarm(i, 2) + swarm(i, 6)/1.2 ;    %update v position
        if swarm (i,1) < 1 
            swarm (i,1) = 1;
        end
        if swarm(i,1) > 50 
            swarm(i,1)=50;
        end
        if swarm (i,2) <1 
            swarm (i,2) = 1;
        end
        if swarm(i,2) > 50
            swarm(i,2)= 50;
        end    
       % u = swarm(i, 1);
        %v = swarm(i, 2);
        value = fobj(swarm(i,1:2));          %Objective function has infinite solution. By taking 500 swarms we have consider the 5000 best solutions.  
        d = sqrt((swarm(i,8)-swarm(i,1)) * (swarm(i,8)-swarm(i,1)) + (swarm(i,9)-swarm(i,2)) * (swarm(i,9)-swarm(i,2)));
        swarm(i,8) = swarm(i,1); swarm (i,9) = swarm(i,2);
        if d>maxd
            maxd = d
        end 
        if value < swarm(i, 7)           % Always True
            swarm(i, 3) = swarm(i, 1);    % update best position of u, updated position.
            swarm(i, 4) = swarm(i, 2);    % update best postions of v,
            swarm(i, 7) = value;          % best updated minimum value Pbest (the best in a iteration)
        end
    end
    [temp, gbest] = min(swarm(:, 7));        % gbest position, global best
    time=time+maxd/inertia;
    %--- updating velocity vectors
    for i = 1 : swarms
        swarm(i, 5) = rand*inertia*swarm(i, 5) + correction_factor*rand*(swarm(i, 3)... 
            - swarm(i, 1)) + correction_factor*rand*(swarm(gbest, 3) - swarm(i, 1));   % u velocity parameters, updated best velocity of swarm. 
        swarm(i, 6) = rand*inertia*swarm(i, 6) + correction_factor*rand*(swarm(i, 4)...
            - swarm(i, 2)) + correction_factor*rand*(swarm(gbest, 4) - swarm(i, 2));   % v velocity parameters
         tEnd=toc(tStart);
    end
    
    %% Plotting the swarm
   % clf    
  %  plot(swarm(:, 1), swarm(:, 2), 'x')   % drawing swarm movements
   % axis([0 50 0 50])
%pause(1)
end
 a= abs(gbest-[10 10]);
 if a(1,1)<= 2  && a(1,2)<= 2 
    success = 1;
 else 
    success = 0;
 end

end