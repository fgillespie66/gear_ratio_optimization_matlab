clear
close all

%% Dependencies
restoredefaultpath               % "clean slate" for your matlab path
addpath(genpath('../casadi_mac')) % make sure you have added your OS-specific casadi folder to MATLAB-Optimization
import casadi.*

%% Output Options
save_data = false;
show_individual_plots = true;
save_surface_plot = false;

%% Derive dynamics

%actuator base parameters
gear1_ratios = [9.5]; %HIP
gear2_ratios = [9.5]; %KNEE

max_height_data = zeros(length(gear1_ratios), length(gear2_ratios));

for i = 1:length(gear1_ratios)
    i;
    gear_ratio1 = gear1_ratios(i);
    for j = 1:length(gear2_ratios)
        gear_ratio2 = gear2_ratios(j);

        motor_base_torque = 2.75; %Nm based on mini cheetah motor, saturation torque
        motor_base_free_speed = 190; %rads per second mini cheetah motor
        motor_torque_intercept = 3.677777777777778; %y-intercept of the power line Nm
        
        %mini cheetah parameters
        body_width = 0.45;
        body_mass = 3.5; %mass of body in kg 
        %need for a constraint later
        l = 0.209; %make sure this matches derive_leg_2DOF
        
        %actuator torque-speed approx
        %first actuator
        p1 = [0,motor_torque_intercept*gear_ratio1];
        p2 = [motor_base_free_speed/gear_ratio1, 0];
        m_ts = (p2(2)-p1(2)) / (p2(1)-p1(1));
        b_ts = motor_torque_intercept*gear_ratio1;
        %second actuator
        p1k = [0,motor_torque_intercept*gear_ratio2];
        p2k = [motor_base_free_speed/gear_ratio2, 0];
        m_tsk = (p2k(2)-p1k(2)) / (p2k(1)-p1k(1));
        b_tsk = motor_torque_intercept*gear_ratio2;
        
        %derive dynamics
        [kinematics,dynamics] = derive_leg_2DOF(gear_ratio1, gear_ratio2, body_mass/2); 
        
        %% Formulate Optimization
        % via trapezoidal Collocation
        
        opti = casadi.Opti(); % Optimization problem
        
        step_scaling = 5;
        N  = 50*step_scaling;   % number of control intervals
        dt = 0.025/step_scaling; % dynamics dt
        T  = N*dt; % duration of stance phase
        
        % ---- decision variables ---------
        X = opti.variable(6,N+1); % [y, theta, dy, dtheta] 
        U = opti.variable(2,N);   % NOW WE HAVE TWO JOINT TORQUES
        F = opti.variable(1,N);   % vertical reaction force
        
        % ---- objective          ---------
        g = -9.81;
        terminal_COM          = kinematics.COM(X(:,end)); 
        terminal_com_y_height = terminal_COM(2);
        terminal_com_y_vel    = terminal_COM(4);
        
        time_to_peak   = -terminal_com_y_vel / g;           % find time when
        projectile_motion = @(t,y,v,a) y + v*t + 0.5*a*t^2; % anonymous function for projectile motion
        max_com_height = projectile_motion(time_to_peak,terminal_com_y_height,terminal_com_y_vel,g);
       
        opti.minimize(-max_com_height); % maximize peak com height
        %opti.minimize(-max_com_height-0.1*terminal_com_y_vel); % maximize peak com height
        
        % ---- dynamic constraints --------
        for k=1:N % loop over control intervals
            Xk  = X(:,k); 
            Xk1 = X(:,k+1);
            Uk  = U(:,k);
            Fk  = F(:,k);
            Ak1 = dynamics.A(Xk1);
            bk1 = dynamics.b(Xk1, Uk, Fk);
            opti.subject_to( Xk1(1:3) - Xk(1:3) == dt*Xk1(4:6) ) % Euler integration - position
            opti.subject_to( Ak1*(Xk1(4:6)-Xk(4:6))  == dt*bk1 ) % Euler integration - velocity
        end
        
        % ---- path constraints -----------
        for k=1:N % loop over control intervals
            Xk  = X(:,k); 
            Xk1 = X(:,k+1);
            Uk  = U(:,k);
            Uk1  = Uk(1);
            Uk2  = Uk(2);
            Fk  = F(:,k);
            Vm1 = Xk(5);
            Vm2 = Xk(6);
            opti.subject_to( Xk1(1) == 0 )      % Foot on ground
            opti.subject_to( Fk >= 0 )          % Unilateral force constraint (no pulling the ground)
            %torques
            opti.subject_to( -motor_base_torque*gear_ratio1 <= Uk1 <= motor_base_torque*gear_ratio1)   % Control limits INCLUDE GEAR RATIO
            opti.subject_to( -motor_base_torque*gear_ratio2 <= Uk2 <= motor_base_torque*gear_ratio2)
            %velocities constraint
            opti.subject_to( -motor_base_free_speed/gear_ratio1 <= Vm1 <= motor_base_free_speed/gear_ratio1)
            opti.subject_to( -motor_base_free_speed/gear_ratio2 <= Vm2 <= motor_base_free_speed/gear_ratio2)
        
            %torque speed curve line torque <= m * velocity + b
            opti.subject_to( (Uk1 - m_ts * Vm1) <= (b_ts) )
            opti.subject_to( (Uk2 - m_tsk * Vm2) >= (-b_tsk) ) 
        
            %add in constraint that whole leg must stay above the ground 
            opti.subject_to( 0 <= l*sin(Xk1(2)) + l*sin(Xk1(3)) + Xk1(1))
            opti.subject_to( 0 <= l*sin(Xk1(3)) + Xk1(1))
        
            %add in some joint limits to reduce search space
            opti.subject_to( 0 <= Xk1(2) <= pi/2 )
            opti.subject_to( pi/2 <= Xk1(3) <= 3*pi/2 )% Knee above ground, avoid singularity
        end
        
        % ---- boundary conditions --------
        opti.subject_to( X(:,1) == [0;deg2rad(60);deg2rad(120);0;0;0] );
        %opti.subject_to( X(:,1) == [0;deg2rad(45);deg2rad(135);0;0;0] );
        
        % ---- initial values for solver ---
        %opti.set_initial(X, );
        %opti.set_initial(U, );
        %opti.set_initial(F, );
        
        % ---- solve NLP              ------
        p_opts = struct('expand',true); % expand to casadi variables to SX (10x speedup)
        opti.solver('ipopt',p_opts);    % set numerical backend
        sol = opti.solve();             % actual solve
        
        %% Simulate Forward the Dynamics
        
        % ---- post-processing        ------
        t = 0:dt:(N*dt);
        z = sol.value(X);
        Fs = sol.value(F);
        
        %---- simulating forward the dyanmics ----
        terminal_COM_sol = kinematics.COM(z(:,end));
        yi = 0;
        vi    = terminal_COM_sol(4);
        t_peak   = -vi / g;  
        zf = z(:,end);
        theta1f = zf(2);
        theta2f = zf(3);
        
        %simulate forward projectile motion 
        MAX_HEIGHT = projectile_motion(full(t_peak),yi,vi,g)
        max_height_data(i, j) = full(MAX_HEIGHT);
        

        if show_individual_plots
            t2 = dt:dt:2*full(t_peak); %simulate to right before impact
            ys = [];
            t1s = [];
            t2s = [];
            gs = [];
            zeros = [];
            for time = t2
                yc = projectile_motion(time,yi,vi,g); %use this with yi = 0 but vi is the same! 
                ys = [ys, full(yc)];
                t1s = [t1s, theta1f];
                t2s = [t2s, theta2f];
                gs = [gs, g];
                zeros = [zeros, 0.0];
            end
            
            %% Animate the Solution
            
            %prep the arrays for plotting
            t2 = t2+N*dt;
            zs = [ys; t1s; t2s; gs; zeros; zeros];
            t_fs = [t, t2];
            z_fs = [z, zs];
            
            %animate the solution
            figure;
            speed = 1;
            animate_simple(t_fs,z_fs,kinematics,speed, gear_ratio1, body_width);
            figure;
            animate_simple(t,z,kinematics,speed, gear_ratio1, body_width);
            
            %% Plot Actuation Efforts + True Torque Speed Curve
            
            %populate the torque curve and the speed curve 
            motor_torques = sol.value(U);
            V1s = [];
            V2s = [];
            U1s = [];
            U2s = [];
            for x = z
                V1s = [V1s, x(5)];
                V2s = [V2s, -x(6)];
            end
            for U12 = motor_torques
                U1s = [U1s, U12(1)];
                U2s = [U2s, -U12(2)];
            end
            
            %remove the initial configuration of the leg from state vector
            V1s = V1s(1:end-1); %DONT REMOVE THE FIRST ONE
            V2s = V2s(1:end-1);
            
            %formatting the figure
            sz = 25;
            %plot
            figure;
            hold on
            %create color map
            len = length(V1s);
            red = [216, 17, 89]/255;
            pink = [143, 45, 86]/255;
            colors_tau1 = [linspace(red(1),pink(1),len)', linspace(red(2),pink(2),len)', linspace(red(3),pink(3),len)'];
            blue = [0, 0, 0]/255;
            lightblue = [20, 20, 20]/255;
            colors_tau2 = [linspace(blue(1),lightblue(1),len)', linspace(blue(2),lightblue(2),len)', linspace(blue(3),lightblue(3),len)'];
            
            line1 = [m_ts,b_ts];
            line2 = [0, motor_base_torque*gear_ratio1, motor_base_free_speed/gear_ratio1, motor_base_torque*gear_ratio1];
            [x_int,y_int] = line_intersection(line1,line2);
            ar2 = area([-5,ceil(motor_base_free_speed/gear_ratio1/5)*5],[ceil(motor_base_torque*gear_ratio1/5)*5,ceil(motor_base_torque*gear_ratio1/5)*5]);
            ar = area([-5,x_int,motor_base_free_speed/gear_ratio1],[motor_base_torque*gear_ratio1,motor_base_torque*gear_ratio1, 0]);
            ar.EdgeColor = '#808080';
            ar2.EdgeColor = '#808080';
            ar.FaceColor = '#ffdda6';
            ar2.FaceColor = '#bfecf2';

            xgrid = -5:0.5:ceil(motor_base_free_speed/gear_ratio1/5)*5; 
            ygrid = 0:0.5:ceil(motor_base_torque*gear_ratio1/5)*5; 


           xl = arrayfun(@(x)xline(x,'Color','#808080','LineWidth',0.2),xgrid);
           yl = arrayfun(@(y)yline(y,'Color','#808080','LineWidth',0.2),ygrid);


            %plot the tau-s constraints
            %plot([0,motor_base_free_speed/gear_ratio1],[motor_torque_intercept*gear_ratio1,0]);
            %plot([0,motor_base_free_speed/gear_ratio1],[motor_base_torque*gear_ratio1,motor_base_torque*gear_ratio1]);
            %plot([motor_base_free_speed/gear_ratio1,motor_base_free_speed/gear_ratio1],[0,motor_base_torque*gear_ratio1]);
            %plot([0,motor_base_free_speed/gear_ratio2],[motor_torque_intercept*gear_ratio2,0]);
            %plot([0,motor_base_free_speed/gear_ratio2],[motor_base_torque*gear_ratio2,motor_base_torque*gear_ratio2]);
            %plot([motor_base_free_speed/gear_ratio2,motor_base_free_speed/gear_ratio2],[0,motor_base_torque*gear_ratio2]);
            
            %plot the trajectory
            scatter(V1s, U1s, sz, colors_tau1, 'filled');
            scatter(V2s, U2s, sz, colors_tau2, 'filled');
            %plot(V1s, U1s);
            %plot(V2s, U2s);
            
            %fplot(@(x) m_ts * x + b_ts); %draw the power limit line
            xlabel("Motor Velocity (rads/s)");
            ylabel("Motor Torques (Nm)");
            title("Actuation Command Overlayed on TS Curve during Stance");
            legend('Unattainable Region', 'Motor Operation Region', 'Hip Torques', 'Knee Torques');
            grid on
            
            %plot the "gem" of strange looping action we found
            figure
            p = plot(V1s, U1s, '-o', 'LineWidth', 2);
            p.MarkerFaceColor = [1 0.5 0];
            p.MarkerSize = 5;
            p.MarkerEdgeColor = [1 0.5 0];
            xlabel("Motor Velocity (rads/s)");
            ylabel("Motor Torques (Nm)");
            title("Hip Joint Path Through T-S Plane Hip Gear: " + gear_ratio1 + "; Knee Gear: " + gear_ratio2);
            grid on
            
            %plot trajectories as a function of time
            LW = 5
            figure
            subplot(3,1,1);
            hold on
            plot(t(:,1:end-1), U1s, 'LineWidth',LW);
            plot(t(:,1:end-1), U2s, 'LineWidth',LW);
            xlabel("Simulation Time - Stance Phase (s)");
            ylabel("Motor Torques (Nm)");
            title("Torque Trajectories of a Single Jump; Hip Gear: " + gear_ratio1 + "; Knee Gear: " + gear_ratio2);
            legend('Hip Joint', 'Knee Joint');
            grid on
            hold off
            subplot(3,1,2);
            hold on
            plot(t(:,1:end-1), V1s, 'LineWidth',LW);
            plot(t(:,1:end-1), V2s, 'LineWidth',LW);
            xlabel("Simulation Time - Stance Phase (s)");
            ylabel("Motor Velocities (rads/s)");
            title("Joint Velocities of a Single Jump; Hip Gear: " + gear_ratio1 + "; Knee Gear: " + gear_ratio2);
            legend('Hip Joint', 'Knee Joint');
            grid on
            hold off
            subplot(3,1,3, 'LineWidth',LW);
            plot(t(:,1:end-1), Fs, 'LineWidth',LW);
            xlabel("Simulation Time - Stance Phase (s)");
            ylabel("Ground Reaction Force (N)");
            title("Ground Reaction Force vs Simulation Time; Hip Gear: " + gear_ratio1 + "; Knee Gear: " + gear_ratio2);
            grid on
        end
    end
end

%% 

if save_data
    % write height data to csv
    curr_time = datestr(datetime('now'));
    filename = pwd + "/Data/" + curr_time + "-heightData.xls";

    % gear ratios, then blank line, then height data
    data = [gear1_ratios; zeros(1, length(gear1_ratios)); max_height_data]
    writematrix(data, filename)
end

%%

% create gear ratio vs max height plot
fig = figure;
hold on
surface(gear1_ratios, gear2_ratios, max_height_data);


title("Peak Jump Height vs Hip and Knee Gear Ratio for 2D Hopper")
xlabel("Hip Gear Ratio");
ylabel("Knee Gear Ratio")
zlabel("Peak Height (m)");
grid on

% to do: fix filename
if save_surface_plot
    filename = pwd + "/Graphics/GearPlot-"+ curr_time + ".png";
    saveas(fig, pwd + filename);
end



%% PLOTS WE SHOULD MAKE

%1DOF CASE
% - max height vs. gear ratio (line plot)
% - trajectory torque as a function of time for the best run
% - graph of the constraints and the torque vs. speed to show how it rides
% the curve 
% - add in mass and do a mass vs optimal gear ratio plot (if we search over
% mass what's the optimal gear ratio as a function of body mass / payload
% mass)

%2DOF CASE
% - max height vs gear ratio 1 vs gear ratio 2 (surface plot) 
% - torque trajectories as a function of time for both joints 
% - graph of constraints 
% - plot of this looping action it does






%% TESTING CONSTRAINTS GRAPHICALLY
    %figure
    %hold on
    %plot([0,motor_base_free_speed/gear_ratio],[motor_torque_intercept*gear_ratio,0]); %draw the power limit line
    %fplot(@(x) m_ts * x + b_ts);
    %legend('one', 'two');
%CORRECTION FACTOR TELLS US IT'S JUST THE y-INTERCEPT THAT'S WRONG

%CITE MATT AND https://github.com/tamaskis/line_intersection-MATLAB











%% NOTES FOR OTHER THINGS


%PLAN
%first expand this to simulate the whole leg forward 
%second expand this to a more complex leg model
%third expand this to a mini cheetah 
%CITE the interias paper, yanran's paper, RUS's work, MATT's work, anything
%else we referenced including the Matt Kelly Optimization tutorial, matt
%and charles' PDF tutorial

%remember, feedback control is REQUIRED just to take out errors in your
%modeling!!! 

%can we use a MATLAB function as a constraint? if we write three
%constraints 



%SOME OTHER NOTES
%add the line as a COST? 
%we learned that the constraint doesn't actually work in this case 
% might need a quper quadratic cost ?? 


%print the constraint violations at each time step 

%TRY DECREASING THE TIMESTEP??? 
%problem gets worse @ higher gear ratios, 
%decreasing timestep significantly helped 

%% OK HERE's HOW WE FIXED IT 

%look @ our integration scheme 
% V_(k+1) = Vk + Uk*dt

%so the arrays look like this !!! 
% [v1, v2, v2, v4 ... vN]
% [u1, u2, u3, ... uN-1] 

%the torques go with the FIRST time step because they get APPLIED to
%calculate the next time step 

%which also asks the question should we constrain Uk to Vk or Vk+1 ??? 


%GEMS
%flinging itself like a trebuchet through the floor is a GEM 
