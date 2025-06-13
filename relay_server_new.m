function relay_server_new(mode, filename, sampling_f, channel, signal, signal_server, pmk_server)
    REALTIME_DATASET_MODE = 1;
    GUT_MODEL_MODE = 2;

    MODE = GUT_MODEL_MODE;

    figure;
    h = animatedline('Color','b');
    h = animatedline('Color','b');
    % ylim([-100000 -75000]);
    ylim([-10 10]);
    grid on;
    tick = 0;
    last_pace = 0;    

    if (exist("mode", "var"))
        MODE = mode;
    end

    if MODE == GUT_MODEL_MODE
        signal_server = tcpserver("0.0.0.0", 8081, "ByteOrder", "little-endian", "ConnectionChangedFcn", @on_connection);
        pmk_server = tcpserver("0.0.0.0", 8082, "ByteOrder", "little-endian", "ConnectionChangedFcn", @on_connection);
        
        fprintf("Waiting for gut model connection on port 8081...\n");
        while ~signal_server.Connected, pause(0.1); end
        fprintf("Gut model connected on port 8081!\n");
        fprintf("Waiting for pacemaker connection on port 8082...\n");
        while ~pmk_server.Connected, pause(0.1); end
        fprintf("Pacemaker connected on port 8082!\n");
        write(pmk_server, MODE, "uint8");
        % 
        % % while signal_server.NumBytesAvailable < 1, pause(0.1); end
        % % mode = read(signal_server, 1, "uint8");
        % fprintf("Received mode from gut model: %d\n", mode);
        % fprintf("Gut model running in GUT_MODEL_MODE.\n");
        % while signal_server.NumBytesAvailable < 1, pause(0.1); end
        % filename = read(signal_server, 1, "uint8");
        % fprintf("filename : %s",filename)
        % while signal_server.NumBytesAvailable < 8, pause(0.1); end
        % sampling_f = read(signal_server, 1, "uint32");
        % while signal_server.NumBytesAvailable < 8, pause(0.1); end
        % channel = read(signal_server, 1, "uint32");

        while true
            % if ~signal_server.Connected
            %     fprintf("Gut model disconnected. Waiting for reconnection...\n");
            %     while ~signal_server.Connected
            %         pause(0.1);
            %     end
            %     fprintf("Gut model reconnected.\n");
            % end
            % if signal_server.NumBytesAvailable >= 1
            %     raw = read(signal_server, signal_server.NumBytesAvailable, "uint8");
            %     fprintf("Bytes received: %d\n", length(raw));
            % 
            %     disp(raw);  % Look at how many bytes and what they are
            % end
            % 
            if signal_server.NumBytesAvailable >= 8
                signal = read(signal_server, 1, "double");
                write(pmk_server, signal, "double");
                fprintf("Received: %.4f\n", signal);
                addpoints(h, tick, signal);
                center = tick - 800;
                xlim([center, center + 1000]);
                tick = tick + 1;
            end
            if pmk_server.NumBytesAvailable >= 8
                pacing = read(pmk_server, 1, "uint8");
                write(signal_server, pacing, "uint8");
                last_pace = pacing;
                addpoints(h, tick, signal);
                center = tick - 800;
                xlim([center, center + 1000]);
                tick = tick + 1;
            end
            drawnow limitrate;
        end

        % pmk_server = tcpserver("0.0.0.0", 8082, "ByteOrder", "little-endian");
        % fprintf("Waiting for pacemaker connection on port 8082...\n");
        % while ~pmk_server.Connected, pause(0.1); end
        % write(pmk_server, MODE, "uint8");
        % % write(pmk_server, uint8(filename), "uint8");
        % % write(pmk_server, sampling_f, "int32");
        % % write(pmk_server, channel, "int32");
        % % fprintf("Signal information sent.\n");
        % Ts = 1/sampling_f;
        
        pause(0.1);

        % while true
        %     if ~signal_server.Connected
        %         fprintf("Gut model disconnected. Waiting for reconnection...\n");
        %         while ~signal_server.Connected
        %             pause(0.1);
        %         end
        %         fprintf("Gut model reconnected.\n");
        %     end
        % 
        %     if ~pmk_server.Connected
        %         fprintf("Pacemaker disconnected. Waiting for reconnection...\n");
        %         while ~pmk_server.Connected
        %             pause(0.1);
        %         end
        %         fprintf("Pacemaker reconnected.\n");
        %     end
        % 
        %     if signal_server.NumBytesAvailable >= 8
        %         signal = read(signal_server, 1, "double");
        %         write(pmk_server, signal, "double");
        %         addpoints(h, tick, signal);
        %         center = tick - 800;
        %         xlim([center, center + 1000]);
        %         tick = tick + 1;
        %     end
        % 
        %     if pmk_server.NumBytesAvailable >= 1
        %         pacing = read(pmk_server, 1, "uint8");
        %         write(signal_server, pacing, "uint8");
        %         last_pace = pacing;
        %         if pacing == 1
        %             disp('Pacing signal relayed.');
        %             xline(tick, 'r', 'Alpha', 0.5);
        %         end
        %     end
        % 
        %     drawnow limitrate;
        % end
    elseif MODE == REALTIME_DATASET_MODE
        pmk_server = tcpserver("0.0.0.0", 8082, "ByteOrder", "little-endian");
        fprintf("Realtime Dataset Mode Start\n")

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
    end
end

function on_connection(server, ~)
    if numel(server.Clients) > 1
        % Close the newest (last) client
        disp("Too many clients. Closing extra connection.");
        flush(server.Clients(end));
        clear(server.Clients(end));  % Disconnect it
    else
        disp("Client accepted.");
    end
end



% function signal_thread(h)
%     fprintf("Start signal thread...");
%     try
%         fprintf("[Signal] Waiting for gut model on port 8081...\n");
%         s_server = tcpserver("0.0.0.0", 8081, "ByteOrder", "little-endian");
% 
%         fprintf("[Signal] Connected to gut model.\n");
% 
%         fprintf("[Signal] Waiting for pacemaker on port 9091...\n");
%         p_client = tcpclient("127.0.0.1", 9091, "ByteOrder", "little-endian");
% 
%         fprintf("[Signal] Connected to pacemaker.\n");
% 
%         tick = 0;
%         while true
%             if s_server.NumBytesAvailable >= 8
%                 signal = read(s_server, 1, "double");
%                 write(p_client, signal, "double");
% 
%                 % Plot
%                 addpoints(h, tick, signal);
%                 xlim([tick - 800, tick + 200]);
%                 drawnow limitrate;
%                 tick = tick + 1;
%             end
%             pause(0.01);
%         end
%     catch ME
%         fprintf("[Signal Error] %s\n", ME.message);
%     end
% end   
% 
% 
% 
% function relay_signal_thread(h)
% 
% 
%     signal_server = tcpserver("0.0.0.0", 8081, "ByteOrder", "little-endian");
%     pmk_server = tcpserver("0.0.0.0", 8082, "ByteOrder", "little-endian");
%     tick = 0;
% 
%     fprintf("Waiting for model connection on port 8081...\n");
%     while ~signal_server.Connected, pause(0.1); end % wait for connect
%     fprintf("Gut model connected.\n");
%     % while signal_server.NumBytesAvailable < 1, pause(0.1); end % wait for mode
%     % mode = read(signal_server, 1, "uint8");
%     % fprintf("Received mode from gut model: %d\n", mode);
%     % fprintf("Gut model running in GUT_MODEL_MODE.\n");
%     try
%         while true
%         % Wait for data from gut model
%         if signal_server.Connected && pmk_server.Connected && signal_server.NumBytesAvailable >= 8
%             signal = read(signal_server, 1, "double");
%             write(pmk_server, signal, "double");
% 
%             % Plot signal
%             addpoints(h, tick, signal);
%             center = tick - 800;
%             xlim([center, center + 1000]);
%             tick = tick + 1;
% 
%             drawnow limitrate;
%         end        
%         pause(0.01); % prevent CPU hog
%     end
% 
%     catch ME
%         fprintf("[Signal Thread] Error: %s\n", ME.message);
%     end
% 
% end
% 
% function relay_pacing_thread(pmk_server, signal_server)
%     fprintf("Waiting for pacemaker connection on port 8082...\n");
%     while ~pmk_server.Connected, pause(0.1); end
%     fprintf("Pacemaker connected.\n");
% 
%     write(pmk_server, MODE, "uint8");
%     write(pmk_server, uint8(filename), "uint8");
%     write(pmk_server, sampling_f, "int32");
%     write(pmk_server, channel, "int32");
%     fprintf("Signal information sent.\n");
% 
%     pause(0.1);
% 
%     Ts = 1/sampling_f;
%     while true
%         try
%             % Wait for pacing data from pacemaker
%             if signal_server.Connected && pmk_server.Connected && pmk_server.NumBytesAvailable >= 1
%                 pacing = read(pmk_server, 1, "uint8");
%                 write(signal_server, pacing, "uint8");
% 
%                 % Plot pacing event
%                 if pacing == 1
%                     fprintf("Pacing signal relayed.\n");
%                     xline(gca, 'r', 'Alpha', 0.5);  % Add red vertical line
%                     drawnow limitrate;
%                 end
%             end
%         catch ME
%             fprintf("[Pacing Thread] Error: %s\n", ME.message);
%         end
% 
%         pause(0.01); % prevent CPU hog
%     end
% end
