function relay_server()
    % Clean up
    clear; clc; close all;

    % Create separate servers for gut model and pacemaker
    gut_server = tcpserver("0.0.0.0", 8081, "ByteOrder", "little-endian");
    pmk_server = tcpserver("0.0.0.0", 8082, "ByteOrder", "little-endian");

    fprintf("Waiting for gut model connection on port 8081...\n");
    while ~gut_server.Connected
        pause(0.1);
    end
    fprintf("Gut model connected.\n");

    fprintf("Waiting for pacemaker connection on port 8082...\n");
    while ~pmk_server.Connected
        pause(0.1);
    end
    fprintf("Pacemaker connected.\n");

    % Setup optional plot
    figure;
    h = animatedline('Color','b');
    ylim([-100000 -75000]);
    grid on;
    tick = 0;

    while true
        % Check disconnections
        if ~gut_server.Connected
            fprintf("Gut model disconnected. Waiting for reconnection...\n");
            while ~gut_server.Connected
                pause(0.1);
            end
            fprintf("Gut model reconnected.\n");
        end
        if ~pmk_server.Connected
            fprintf("Pacemaker disconnected. Waiting for reconnection...\n");
            while ~pmk_server.Connected
                pause(0.1);
            end
            fprintf("Pacemaker reconnected.\n");
        end

        % Relay if both are connected
        if gut_server.NumBytesAvailable >= 8
            gut_signal = read(gut_server, 1, "double");
            write(pmk_server, gut_signal, "double");
            addpoints(h, tick, gut_signal);
            center = tick - 800;  % Scroll so latest sample is at 3/4 of 1000
            xlim([center, center + 1000]);
            tick = tick + 1;
        end

        if pmk_server.NumBytesAvailable >= 1
            pacing = read(pmk_server, 1, "uint8");
            write(gut_server, pacing, "uint8");

    % If pacing occurred, draw vertical line at current tick
    if pacing == 1
    % if pacing == 1 && tick > 0
        disp('pace recieved');
        xline(tick, 'r', 'Alpha', 0.5);  % Red vertical line
    end
        end

        drawnow limitrate;
        pause(0.03);
    end
end
