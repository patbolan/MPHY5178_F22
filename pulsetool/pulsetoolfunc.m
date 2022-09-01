function pulsetoolfunc(action)
% pulsetoolfunc - the message handling switchboard for pulsetool.m

% PJB 06.10.00
% Revision comments in pulsetool.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Specify all the globals that will be updated
global filename;
global b1 phase Nrf;
global b1max Tp pulsephase;

global T1 T2;
global MOX MOY MOZ;
global simbw offsetsteps offset;
global offsets Mx My Mz;
global B1 b1tx b1ty b1null;
global pauseval;

switch(action)

    case 'loadfile'
        [fname,pname] = uigetfile('*.RF','Varian Pulse', 'pulseshapes\');
        pulsetoolfunc('gui_to_data');
        if fname~=0
            filename = strcat(pname,fname);
            [b1,phase,Nrf] = rdpls(filename);
        end
        pulsetoolfunc('data_to_gui');
        pulsetoolfunc('plotpulse');

    case 'calcflipangle'
        % Calculate the integral of the amp. Can be used to
        %	calculate flip angles: B1max = flipangle/(integral*Tp)
        % Note that this takes a few seconds to display
        pulsetoolfunc('gui_to_data');
        pulsetoolfunc('gui_to_data');
        integral = 0;
        for idx=1:Nrf
            integral = integral + b1(idx)/Nrf * cos(phase(idx));
        end
        disp(' ');
        disp('**** Flip Angle ****');
        disp('INTEGRAL is the integral of the x-component of the pulse amplitude.');
        s = sprintf('     INTEGRAL = %f', integral);
        disp(s);
        disp('It can be used to calculate the flip angle: ');
        disp('     flip angle = B1max(Hz) * INTEGRAL * Tp * 360');
        disp(sprintf('For B1max = %f, flip angle (deg) is %f', b1max, b1max*integral*Tp*360));
        disp(sprintf('On-resonance 90deg at gamB1 = %.1f', b1max * 90 /(b1max*integral*Tp*360) ));
        disp(' ');

    case 'plotpulse'
        % Need to clear all the local variables
        clear time dt phi_dot alpha Weff K;
        pulsetoolfunc('gui_to_data')

        time = b1;
        for idx = 1:Nrf
            time(idx) = Tp * idx * 1000000 / Nrf;
        end

        hDisplay = figure(101);
        set(gcf, 'Name', 'Input Pulse');

        % Plot the B1
        subplot(2,1,1);
        plot(time, b1);
        ylabel('Amplitude');

        % Plot phase
        subplot(2,1,2)
        plot(time, phase);
        ylabel('\phi (rad)');

        %movegui(hDisplay, 'northwest')



    case 'data_to_gui'
        ptrhdl = findobj('Tag', 'EditB1max');
        set(ptrhdl, 'String', b1max);

        ptrhdl = findobj('Tag', 'EditTp');
        set(ptrhdl, 'String', Tp);

        ptrhdl = findobj('Tag', 'EditPulsephase');
        set(ptrhdl, 'String', pulsephase);

        ptrhdl = findobj('Tag', 'EditFilename');
        set(ptrhdl, 'String', filename);


        ptrhdl = findobj('Tag', 'EditT1');
        set(ptrhdl, 'String', T1);

        ptrhdl = findobj('Tag', 'EditT2');
        set(ptrhdl, 'String', T2);

        ptrhdl = findobj('Tag', 'EditSimBW');
        set(ptrhdl, 'String', simbw*2);

        ptrhdl = findobj('Tag', 'EditOffsetsteps');
        set(ptrhdl, 'String', offsetsteps);


        ptrhdl = findobj('Tag', 'EditMOX');
        set(ptrhdl, 'String', MOX);

        ptrhdl = findobj('Tag', 'EditMOY');
        set(ptrhdl, 'String', MOY);

        ptrhdl = findobj('Tag', 'EditMOZ');
        set(ptrhdl, 'String', MOZ);


        ptrhdl = findobj('Tag', 'EditOffset');
        set(ptrhdl, 'String', offset);

        ptrhdl = findobj('Tag', 'EditPause');
        set(ptrhdl, 'String', pauseval);



    case 'gui_to_data'
        ptrhdl = findobj('Tag', 'EditB1max');
        b1max = eval(get(ptrhdl, 'String'));

        ptrhdl = findobj('Tag', 'EditTp');
        Tp = eval(get(ptrhdl, 'String'));

        ptrhdl = findobj('Tag', 'EditPulsephase');
        pulsephase = eval(get(ptrhdl, 'String'));
        % Convert to radians
        pulsephase = pulsephase * pi / 180;

        ptrhdl = findobj('Tag', 'EditFilename');
        filename = get(ptrhdl, 'String');



        ptrhdl = findobj('Tag', 'EditT1');
        T1 = eval(get(ptrhdl, 'String'));

        ptrhdl = findobj('Tag', 'EditT2');
        T2 = eval(get(ptrhdl, 'String'));

        ptrhdl = findobj('Tag', 'EditSimBW');
        simbw = eval(get(ptrhdl, 'String'))./2;

        ptrhdl = findobj('Tag', 'EditOffsetsteps');
        offsetsteps = eval(get(ptrhdl, 'String'));



        ptrhdl = findobj('Tag', 'EditMOX');
        MOX = eval(get(ptrhdl, 'String'));

        ptrhdl = findobj('Tag', 'EditMOY');
        MOY = eval(get(ptrhdl, 'String'));

        ptrhdl = findobj('Tag', 'EditMOZ');
        MOZ = eval(get(ptrhdl, 'String'));


        ptrhdl = findobj('Tag', 'EditOffset');
        offset = eval(get(ptrhdl, 'String'));

        ptrhdl = findobj('Tag', 'EditPause');
        pauseval = eval(get(ptrhdl, 'String'));

        % Caclulate complex components of the pulse
        b1tx = b1;
        b1ty = b1;
        b1null = b1;
        for idx = 1:Nrf
            b1tx(idx) = 2 * pi * b1max * ( b1(idx) * cos(pulsephase + phase(idx)) );
            b1ty(idx) = 2 * pi * b1max * ( b1(idx) * sin(pulsephase+ phase(idx)) );
            b1null(idx)=0;
        end


    case 'fft'
        hFFT = figure(103);
        Nf = 10000;
        SW =(Nrf-1)/Tp;
        freq = -(SW/2):SW/(Nf-1):SW/2;
        plot(freq ./ 1000, abs(fftshift(fft(b1, Nf)))./(Nrf))        
        %set(gca, 'Xlim', 1000*[min(offsets(:)) max(offsets(:))]);

        set(gca, 'Xlim', [-simbw simbw]./ 1000);
        set(gca, 'Box', 'off');
        legend('FFT');
        xlabel('offset freq (kHz)');


    case 'simulate'
        hSim = figure(102);
        set(gcf, 'Name', 'Simulation Results');

        % Check Booleans
        bInitM0 = get(findobj('Tag', 'CheckboxInitM0'), 'Value');
        bRefocus = get(findobj('Tag', 'CheckboxRefocus'), 'Value');

        % Reset the offset
        offsets = zeros(offsetsteps,1);

        dt = Tp/Nrf;
        MO = [MOX MOY MOZ];

        for idx = 1:offsetsteps+1

            offsets(idx) = simbw * 2 * pi * (2*(idx-1)/offsetsteps-1);

            % Use old magnetization if requested
            if ~bInitM0
                MO = [Mx(idx) My(idx) Mz(idx)];
            end

            M = blochRK4_loop(b1tx, b1ty, Tp, offsets(idx), 1/T1, 1/T2, MO);

            
            
            % Refocusing
            if(bRefocus)
                M2 = blochRK4_loop(b1null, b1null, Tp/2, -offsets(idx), 1/T1, 1/T2, M);
                M=M2;
            end;

            Mx(idx) = M(1);
            My(idx) = M(2);
            Mz(idx) = M(3);
            Mxy(idx) = sqrt(M(1)^2 + M(2)^2);

        end

        % Check plotting options
        bShowMx = get(findobj('Tag', 'CheckboxMx'), 'Value');
        bShowMy = get(findobj('Tag', 'CheckboxMy'), 'Value');
        bShowMz = get(findobj('Tag', 'CheckboxMz'), 'Value');
        bShowMxy = get(findobj('Tag', 'CheckboxMxy'), 'Value');
        bShowPhs = get(findobj('Tag', 'CheckboxPhase'), 'Value');

        offsets = offsets / (1000 * 2 * pi);
        legendstr='';
        cla;
        hold on;
        if(bShowMx)
            plot(offsets, Mx, '-r');
            legendstr = strvcat(legendstr, 'Mx');
        end;
        if (bShowMy)
            hold on;
            plot(offsets, My, '-b');
            legendstr = strvcat(legendstr, 'My');
        end;
        if (bShowMz)
            plot(offsets, Mz, '-k');
            legendstr = strvcat(legendstr, 'Mz');
        end;
        if (bShowMxy)
            plot(offsets, Mxy, '-m');
            legendstr = strvcat(legendstr, 'Mxy');
        end;
        if (bShowPhs)
            plot(offsets, atan2(Mx,My)/pi, '-g');
            legendstr = strvcat(legendstr, 'Phs');
        end;

        legend(legendstr);
        %plot(offsets, Mx, offsets, My, offsets, Mz, offsets, Mxy);
        %legend('Mx','My', 'Mz', 'Mxy');
        xlabel('offset freq (kHz)');
        hold off;

    case 'rotate'
        figure(103);
        eye = get(gca, 'CameraPosition');
        theta = pi/4;
        rotate = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];

        %eye(1) = eye(1)*cos(theta) + eye(2)*sin(theta);
        %eye(2) = -eye(1)*sin(theta) + eye(2)*cos(theta);
        eye = eye*rotate;
        set(gca, 'CameraPosition', eye);
        set(gca, 'DataAspectRatio', [1 1 1.25]);

    case 'play'
        % This will show an animation of the rotation of M
        %	over a single pulse.
        pulsetoolfunc('gui_to_data');
        figure(103);
        rotaxis = gca;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set up axes
        cla;
        as = 1; %axis size
        textpos = 1.1; % offset from line to text
        axis([-as as -as as -as as]);
        hx = line([0 1], [0 0], [0 0], 'Color', 'k');
        hy = line([0 0], [0 1], [0 0], 'Color', 'k');
        hz = line([0 0], [0 0], [0 1], 'Color', 'k');
        text(textpos, 0, 0, 'X');
        text(0, textpos, 0, 'Y');
        text(0, 0, textpos, 'Z');set(gca, 'Visible', 'off');
        viewpos = 2;
        %set(gca, 'CameraPosition', [viewpos viewpos viewpos/2]);
        %set(gca, 'CameraTarget', [0 0 0]);
        set(rotaxis, 'Visible', 'off');
        set(rotaxis, 'DataAspectRatio', [1 1 1]);
        set(rotaxis, 'Projection', 'orthographic');
        set(rotaxis, 'CameraPosition', [9 12 9]);

        M = [MOX MOY MOZ];
        offsetrad = offset * 2 * pi;
        B1effmax = max(offsetrad, b1max*2*pi);
        B1eff = [0 0 offsetrad]/B1effmax;
        Mline = line(...
            'XData', [0 M(1)], ...
            'YData', [0 M(2)], ...
            'ZData', [0 M(3)], ...
            'Color', 'r', 'LineWidth', 2);

        % Get Checkbox values
        ptrhdl = findobj('Tag', 'CheckboxB1eff');
        bShowBeff = get(ptrhdl, 'Value');
        ptrhdl = findobj('Tag', 'CheckboxTrace');
        bShowTrace = get(ptrhdl, 'Value');
        ptrhdl = findobj('Tag', 'CheckboxFrame');
        bB1effFrame = get(ptrhdl, 'Value');
        ptrhdl = findobj('Tag', 'CheckboxShadow');
        bShadow = get(ptrhdl, 'Value');

        if bShowBeff == 1
            B1line = line(...
                'XData', [0 B1eff(1)], ...
                'YData', [0 B1eff(2)], ...
                'ZData', [0 B1eff(3)], ...
                'Color', 'b', 'LineWidth', 2);
        end

        % Draw shadow
        if bShadow == 1
            shadow = line(...
                'XData', [0 B1eff(1)], ...
                'YData', [0 B1eff(2)], ...
                'ZData', [0 0], ...
                'Color', [.5 0 0], 'LineWidth', 2);
        end


        timecounter = text(0, 0, -textpos, 't=0');
        for idx = 1:Nrf

            % Simulate Bloch
            M0 = M;
            %M = blochRK4_loop(b1tx(idx), b1ty(idx), Tp/Nrf, offsetrad, 1/T1, 1/T2, MO);
            M = blochRK4(M0, b1tx(idx), b1ty(idx), offsetrad, 1/T1, 1/T2, Tp/Nrf)
            
            
            set(Mline,'XData',[0 M(1)]);
            set(Mline,'YData',[0 M(2)]);
            set(Mline,'ZData',[0 M(3)]);

            % Update Shadow
            if bShadow
                set(shadow ,'XData',[0 M(1)]);
                set(shadow ,'YData',[0 M(2)]);
            end
            if bShowTrace ==1
                % This makes a trace, but in doing so adds huge numbers
                %	of graphics primitives and slows the drawing down lots.

                %epsilon = 0.001;
                %line([M(1) M(1)-epsilon], [M(2) M(2)-epsilon], [M(3) M(3)-epsilon], 'Color', 'r');
                line([0 M(1)], [0 M(2)], [0 M(3)], 'Color', 'r');
            end

            if bShowBeff == 1
                W(1) = b1tx(idx);
                W(2) = b1ty(idx);
                % Recaculate Weffective.
                % There's a component due to the changing phi!
                if idx>1
                    W(3) = offsetrad + (phase(idx)-phase(idx-1))/(Tp/Nrf);
                else
                    W(3) = offsetrad;
                end
                W = W/B1effmax;

                set(B1line,'XData',[0 W(1)]);
                set(B1line,'YData',[0 W(2)]);
                set(B1line,'ZData',[0 W(3)]);

                if bShowTrace ==1
                    % This makes a trace, but in doing so adds huge numbers
                    %	of graphics primitives and slows the drawing down lots.
                    line([W(1) W(1)], [W(2) W(2)], [W(3) W(3)], 'Color', 'b');
                end

                % This will set the camera to look along the B1eff vector.
                %	The magetization should always rotate clockwise, and never
                % 	change length in the image.
                if bB1effFrame
                    if max(abs(W)) ~= 0
                        W = W/max(abs(W));
                        set(rotaxis, 'CameraPosition', [W(1) W(2) W(3)]);
                        set(rotaxis, 'DataAspectRatio', [1 1 1]);
                    end
                end

            end

            % Show time of pulse at bottom
            s = sprintf('t = %.0f \\mus', Tp*idx/Nrf*1000000);
            set(timecounter, 'String', s);


            % The built-in function pause only works with integer
            %	number of seconds. These values will vary based on
            %	the machine speed.
            pause(pauseval)
%             pause(0.0000001);
%             for idx = 1:pauseval
%                 cos(idx);
%             end
        end


end