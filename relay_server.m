% === relay_server.m ===
function relay_server(mode, filename, sampling_f, channel, signal)
    REALTIME_DATASET_MODE = 1;
    GUT_MODEL_MODE = 2;

    MODE = GUT_MODEL_MODE;

    figure;
    h = animatedline('Color','b');
    ylim([-100000 -75000]);
    grid on;
    tick = 0;
    last_pace = 0;

    

    if (exist("mode", "var"))
        MODE = mode;
    end

    while true              
        if MODE == GUT_MODEL_MODE
            fprintf("Waiting for model connection on port 8081...\n");
            while ~signal_server.Connected, pause(0.1); end
            while signal_server.NumBytesAvailable < 1, pause(0.1); end
            mode = read(signal_server, 1, "uint8");
            fprintf("Received mode from gut model: %d\n", mode);
            fprintf("Gut model running in GUT_MODEL_MODE.\n");
            break;
        elseif MODE == REALTIME_DATASET_MODE
            fprintf("Realtime Dataset Mode Start\n")
            pmk_server = tcpserver("0.0.0.0", 8082, "ByteOrder", "little-endian");

            fprintf("Waiting for pacemaker connection on port 8082...\n");
            while ~pmk_server.Connected, pause(0.1); end
            fprintf("Pacemaker connected.\n");

            write(pmk_server, MODE, "uint8");
            write(pmk_server, uint8(filename), "uint8");
            write(pmk_server, sampling_f, "int32");
            write(pmk_server, channel, "int32");
            fprintf("Signal information sent.\n");
            Ts = 1/sampling_f;
            
            pause(0.1);

            % Send data at ~32 Hz
            for i = 1:length(signal)
                if ~pmk_server.Connected
                    fprintf("Pacemaker disconnected. Waiting for reconnection...\n");
                    while ~pmk_server.Connected
                        pause(0.1);
                    end
                    fprintf("Pacemaker reconnected.\n");
                end
                
                addpoints(h, tick, signal(i));
                center = tick - 800;
                xlim([center, center + 1000]);
                tick = tick + 1;

                tic;
                write(pmk_server, signal(i), "double");
                % disp(signal(1));              

                while toc < Ts
                    % active wait (spins CPU but tighter timing)
                end
        
                % fprintf("%d\n",i);
            end
            break;
        end
    end
end
