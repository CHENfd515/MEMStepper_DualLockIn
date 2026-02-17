classdef MEMStepper_DualLockIn < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        MEMSStepperControlUIFigure      matlab.ui.Figure
        GridLayout                      matlab.ui.container.GridLayout
        CapacitiveMotionDetectionPanel  matlab.ui.container.Panel
        GridLayout17                    matlab.ui.container.GridLayout
        PlotCapLossButton               matlab.ui.control.Button
        StartStopCapButton              matlab.ui.control.Button
        CapValue                        matlab.ui.control.Label
        ClearCapButton                  matlab.ui.control.Button
        AxesCapReadout                  matlab.ui.control.UIAxes
        AxesCapFitting                  matlab.ui.control.UIAxes
        GridLayout6                     matlab.ui.container.GridLayout
        RecordMovieButton               matlab.ui.control.Button
        StartPreviewButton              matlab.ui.control.Button
        ContinuousButton                matlab.ui.control.Button
        SnapshotButton                  matlab.ui.control.Button
        StatusPanel                     matlab.ui.container.Panel
        GridLayout16                    matlab.ui.container.GridLayout
        VzGauge                         matlab.ui.control.SemicircularGauge
        VzGaugeLabel                    matlab.ui.control.Label
        V1Gauge                         matlab.ui.control.SemicircularGauge
        V1GaugeLabel                    matlab.ui.control.Label
        V3Gauge                         matlab.ui.control.SemicircularGauge
        V3GaugeLabel                    matlab.ui.control.Label
        V2Gauge                         matlab.ui.control.SemicircularGauge
        V2GaugeLabel                    matlab.ui.control.Label
        MEMScontrolreadoutPanel         matlab.ui.container.Panel
        GridLayout15                    matlab.ui.container.GridLayout
        SR844ConnectButton              matlab.ui.control.Button
        SR844DropDown                   matlab.ui.control.DropDown
        SR844Label                      matlab.ui.control.Label
        ZeroVzButton                    matlab.ui.control.Button
        FlashVzButton                   matlab.ui.control.Button
        SetVzButton                     matlab.ui.control.Button
        VzField                         matlab.ui.control.NumericEditField
        VdField                         matlab.ui.control.NumericEditField
        WaitTimesEditField              matlab.ui.control.NumericEditField
        WaitTimesEditFieldLabel         matlab.ui.control.Label
        SetVdButton                     matlab.ui.control.Button
        SR830Label                      matlab.ui.control.Label
        OUTPUTOFFButton                 matlab.ui.control.Button
        TabGroup3                       matlab.ui.container.TabGroup
        DirectTab                       matlab.ui.container.Tab
        V1Label                         matlab.ui.control.NumericEditField
        V2Label                         matlab.ui.control.NumericEditField
        V3Label                         matlab.ui.control.NumericEditField
        V3KnobLabel                     matlab.ui.control.Label
        V3Knob                          matlab.ui.control.Knob
        V2KnobLabel                     matlab.ui.control.Label
        V2Knob                          matlab.ui.control.Knob
        V1KnobLabel                     matlab.ui.control.Label
        V1Knob                          matlab.ui.control.Knob
        SetButton                       matlab.ui.control.Button
        SteppingTab                     matlab.ui.container.Tab
        SyncCapFittingCheckBox          matlab.ui.control.CheckBox
        MinusQuarterButton              matlab.ui.control.Button
        PlusQuarterButton               matlab.ui.control.Button
        MinusHalfButton                 matlab.ui.control.Button
        PlusHalfButton                  matlab.ui.control.Button
        Minus18Button                   matlab.ui.control.Button
        ZeroStepsButton                 matlab.ui.control.Button
        StepsEditField                  matlab.ui.control.NumericEditField
        StepsEditFieldLabel             matlab.ui.control.Label
        PhaseEditField                  matlab.ui.control.NumericEditField
        PhaseEditFieldLabel             matlab.ui.control.Label
        DriveVoltageLabel               matlab.ui.control.Label
        DriveVoltageEditField           matlab.ui.control.NumericEditField
        Plus18Button                    matlab.ui.control.Button
        InitButton                      matlab.ui.control.Button
        PositioningTab                  matlab.ui.container.Tab
        ZeroCapButton                   matlab.ui.control.Button
        StepsCapEditField               matlab.ui.control.NumericEditField
        StepsCapEditFieldLabel          matlab.ui.control.Label
        CloseLoopSwitch                 matlab.ui.control.Switch
        CloseLoopSwitchLabel            matlab.ui.control.Label
        LogCapDataButton                matlab.ui.control.Button
        ReadChannelDropDown             matlab.ui.control.DropDown
        ReadChannelDropDownLabel        matlab.ui.control.Label
        MicroStepDropDown               matlab.ui.control.DropDown
        MicroStepDropDownLabel          matlab.ui.control.Label
        GOButton                        matlab.ui.control.Button
        SpeedstepsEditFieldLabel        matlab.ui.control.Label
        SpeedstepsEditField             matlab.ui.control.NumericEditField
        StepsPhyEditField               matlab.ui.control.NumericEditField
        StepsPhyEditFieldLabel          matlab.ui.control.Label
        TargetStepEditFieldLabel        matlab.ui.control.Label
        TargetStepEditField             matlab.ui.control.NumericEditField
        InitResetButton                 matlab.ui.control.Button
        DriveVoltageEditField_2Label    matlab.ui.control.Label
        DriveVoltageEditField_2         matlab.ui.control.NumericEditField
        PositioningFeedbackTab          matlab.ui.container.Tab
        PIDwaitsEditFieldLabel          matlab.ui.control.Label
        PIDwaitsEditField               matlab.ui.control.NumericEditField
        LogPlotSwitch                   matlab.ui.control.Switch
        LogPlotSwitchLabel              matlab.ui.control.Label
        KdEditField                     matlab.ui.control.NumericEditField
        KdEditFieldLabel                matlab.ui.control.Label
        KiEditField                     matlab.ui.control.NumericEditField
        KiEditFieldLabel                matlab.ui.control.Label
        KpEditField                     matlab.ui.control.NumericEditField
        KpEditFieldLabel                matlab.ui.control.Label
        SCANButton                      matlab.ui.control.Button
        ScanStepEditFieldLabel_2        matlab.ui.control.Label
        ScanStepEnd                     matlab.ui.control.NumericEditField
        ScanMicrostepDropDown           matlab.ui.control.DropDown
        ScanMicrostepLabel              matlab.ui.control.Label
        FilterWinEditField              matlab.ui.control.NumericEditField
        FilterWinEditFieldLabel         matlab.ui.control.Label
        CANCELButton                    matlab.ui.control.Button
        PIDMicrostepDropDown            matlab.ui.control.DropDown
        FinalMicrostepLabel             matlab.ui.control.Label
        StepMinLossEditField            matlab.ui.control.NumericEditField
        StepMinLossEditFieldLabel       matlab.ui.control.Label
        FinalStepEditFieldLabel         matlab.ui.control.Label
        FinalStepEditField              matlab.ui.control.NumericEditField
        GOButton_2                      matlab.ui.control.Button
        CurrentStepEditField            matlab.ui.control.NumericEditField
        CurrentStepLabel                matlab.ui.control.Label
        SpeedstepsEditField_2Label      matlab.ui.control.Label
        SpeedstepsEditField_2           matlab.ui.control.NumericEditField
        ScanStepEditFieldLabel          matlab.ui.control.Label
        ScanStepSrt                     matlab.ui.control.NumericEditField
        InitResetButton_2               matlab.ui.control.Button
        DriveVoltageEditField_3Label    matlab.ui.control.Label
        DriveVoltageEditField_3         matlab.ui.control.NumericEditField
        SR830ConnectButton              matlab.ui.control.Button
        SR830DropDown                   matlab.ui.control.DropDown
        ImageMotionDetectionPanel       matlab.ui.container.Panel
        GridLayout2                     matlab.ui.container.GridLayout
        ClearGraphsButton               matlab.ui.control.Button
        UIAxes6                         matlab.ui.control.UIAxes
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxes3                         matlab.ui.control.UIAxes
        UIAxes4                         matlab.ui.control.UIAxes
        UIAxes5                         matlab.ui.control.UIAxes
        TabGroup2                       matlab.ui.container.TabGroup
        CameraControlTab                matlab.ui.container.Tab
        GridLayout3                     matlab.ui.container.GridLayout
        MicroscopeObjectiveDropDown     matlab.ui.control.DropDown
        MicroscopeObjectiveDropDownLabel  matlab.ui.control.Label
        GridLayout8                     matlab.ui.container.GridLayout
        ResolutionLabel                 matlab.ui.control.Label
        Resolution                      matlab.ui.control.DropDown
        GridLayout7                     matlab.ui.container.GridLayout
        DeviceLabel                     matlab.ui.control.Label
        DeviceName                      matlab.ui.control.DropDown
        GainSlider                      matlab.ui.control.Slider
        GainSliderLabel                 matlab.ui.control.Label
        ExposureSliderLabel             matlab.ui.control.Label
        ExposureSlider                  matlab.ui.control.Slider
        MotionDetectionTab              matlab.ui.container.Tab
        GridLayout4                     matlab.ui.container.GridLayout
        GridLayout10                    matlab.ui.container.GridLayout
        radialpitchumEditField          matlab.ui.control.NumericEditField
        radialpitchumEditFieldLabel     matlab.ui.control.Label
        AlignpeakCheckBox               matlab.ui.control.CheckBox
        ChooseradialrangeButton         matlab.ui.control.Button
        ChoosethetarangeButton          matlab.ui.control.Button
        PickcircumferenceButton         matlab.ui.control.Button
        MiscTab                         matlab.ui.container.Tab
        GridLayout9                     matlab.ui.container.GridLayout
        GridLayout13                    matlab.ui.container.GridLayout
        MinScaleEditField               matlab.ui.control.NumericEditField
        MinScaleEditFieldLabel          matlab.ui.control.Label
        UseGPUCheckBox                  matlab.ui.control.CheckBox
        GridLayout12                    matlab.ui.container.GridLayout
        VariableNameEditField           matlab.ui.control.EditField
        VarNameLabel                    matlab.ui.control.Label
        PreviewRefDesignButton          matlab.ui.control.Button
        AlignsnapshottoRefDesignExperimentalButton  matlab.ui.control.Button
        GridLayout11                    matlab.ui.container.GridLayout
        ReferenceDesignEditFieldLabel   matlab.ui.control.Label
        ReferenceDesignEditField        matlab.ui.control.EditField
        PickButton                      matlab.ui.control.Button
        SavesnapshottofileButton        matlab.ui.control.Button
        SyncSnapshotinworkspaceCheckBox  matlab.ui.control.CheckBox
        UIAxes                          matlab.ui.control.UIAxes
    end

    properties (Access = private)
        % --- Video Acquisition & Imaging ---
        video           % Videoinput object for image capture
        hImage          % Image handle for GUI visualization
        hImage2         % Secondary image handle (e.g., for processed frames)
        res             % Video resolution settings [Width, Height]
        C               % Calculated center coordinates for motion reference
        snapshot        % Temporary buffer for a single captured frame
        rect_th         % User-selected ROI for theta (angular) displacement detection
        rect_r          % User-selected ROI for radial (linear) displacement detection
        
        % --- Hardware Communication (VISA) ---
        LockIn          % VISA object for SR830 (typically for driving/control)
        LockIn_Read     % VISA object for secondary Lock-In (SR844) for sensing
        prev_closed = false % Safety flag to verify previous VISA session closure
        positioningSTOP = true    % Global interrupt flag to stop feedback and movement
        plotCapFBloss = false     % Global plot flag for detecting feedback process
        default_wait_time = 0.003 % 3ms, Ensures the Lock-in amplifier's integration (Time Constant) has settled.
        
        % --- Fitting & Coordinate Data ---
        k_fit           % Slope of the linear fit: d(Cap)/d(Step)
        b_fit           % Intercept of the linear fit
        X               % Current X-coordinate in Cartesian space
        Y               % Current Y-coordinate in Cartesian space
        THETA           % Angular position in polar coordinates
        R               % Radial position in polar coordinates
        Xi              % Initial X-coordinate reference
        Yi              % Initial Y-coordinate reference
        
        % --- Threshold & Targeted Coordinates ---
        THETAth         % Threshold or target theta value
        Rth             % Threshold or target radial value
        Xth             % Threshold or target X position
        Yth             % Threshold or target Y position
        
        % --- Calculated Displacements ---
        disp_th         % Measured theta displacement, in degrees
        disp_r          % Measured radial displacement, in micrometers (um)
        
        % --- Phase Data ---
        phase_th        % Control phase corresponding to theta displacement
        phase_r         % Control phase corresponding to radial displacement
        
        % --- Recording & Media ---
        refImg          % Reference image for digital image correlation or comparison
        vw              % VideoWriter object for saving video files to disk
        volt_amp_ratio  % Voltage amplification ratio (Input-to-Output gain)
        
        % --- Capacitance Feedback & Plotting ---
        tim                 % Timer object for periodic capacitance data acquisition
        capData             % Buffer storing raw capacitance readout time-series
        capPlot             % Graphics handle for the capacitance readout plot
        capFitData          % Matrix storing [Step, Cap] pairs for fitting and FB analysis
        StepsCapZero = -inf % Cap equal point

        % --- Log File ---
        LogFileName = 'mems.log';

        % --- Timer Update: Parameters for app.zeroingCap() Scan Process ---
        scanTimer           % Dedicated timer controlling non-blocking scan steps
        scanIndex = 0       % Current scan iteration index
        totalSteps = 0      % Total number of scan steps to execute
        scanData            % Buffer storing scan results as [Step, Capacitance]
        stepInc             % Increment per scan step in physical step units
        zeroVolt            % Voltage amplitude used during zeroing scan
        zeroMicrostep       % Microstep resolution used during scan
        zeroSpeed           % Motion speed used during scan
        isScanning = false  % Flag indicating whether zeroing scan is active

    end

    
    methods (Static)
        function v=vphase(phase, shift)
            phase = mod(phase-shift, 2*pi);
            v = zeros(size(phase));
            v(phase<=4*pi/3) = -sin(phase(phase<=4*pi/3)*1.5);
            v(phase<=pi) = 1;
            v(phase<=pi/3) = sin(phase(phase<=pi/3)*1.5);
        end
    end


    methods (Access = private)
    
        function updateAvailableDevices(app)
            info = imaqhwinfo('winvideo');
            app.DeviceName.Items = arrayfun(@(x) x.DeviceName, info.DeviceInfo,'UniformOutput',false);
            app.DeviceName.ItemsData = arrayfun(@(x) x, info.DeviceInfo, 'UniformOutput', false);
        end
        
        function updateAvailableResolution(app)
            if ~isempty(app.DeviceName.Value)
                app.Resolution.Items = app.DeviceName.Value.SupportedFormats;
                app.Resolution.Value = app.DeviceName.Value.DefaultFormat;
                app.Resolution.UserData = app.DeviceName.Value;
            end
        end
        
        function cameraConfig(app, exp, gain)
            
            prev=imaqfind();
            while ~isempty(prev)
                stop(prev(1));
                stoppreview(prev(1));
                delete(prev(1));
                prev=prev(2:end);
            end
            app.Resolution.Value
            app.video = videoinput('winvideo', app.Resolution.UserData.DeviceID, app.Resolution.Value);
            src = getselectedsource(app.video);
            
            app.video.FramesPerTrigger = 1;
            app.video.TriggerRepeat = Inf;
            
            try
                if isprop(src, 'ExposureMode')
                    src.ExposureMode = 'manual';
                end
                if isprop(src, 'Exposure')
                    src.Exposure = exp;
                end
                if isprop(src, 'GainMode')
                    src.GainMode = 'manual';
                end
                if isprop(src, 'Gain')
                    src.Gain = gain;
                end
                if isprop(src, 'GammaMode')
                    src.GammaMode = 'manual';
                end
                if isprop(src, 'Gamma')
                    src.Gamma = 50;
                end
            catch
            end
        end
    
        function cameraStartPreview(app)
            app.res = app.video.VideoResolution;
            if isempty(app.hImage)
                app.hImage=image(app.UIAxes, zeros(fliplr(app.res)));
            end
            axis(app.UIAxes, 'equal');
            axis(app.UIAxes, [1, app.res(1), 1, app.res(2)])
            preview(app.video, app.hImage);
        end
        
        function cameraStopPreview(app)
            if isa(app.video, 'videoinput')
                stoppreview(app.video);
            end
            app.video=[];
        end
     
        function scanUSBPort(app)
            % SCANUSBPORT Scans system for VISA-compliant instruments and filters for Lock-in Amplifiers.
            
            % Retrieve a table of all connected VISA devices (USB, GPIB, TCP/IP, etc.)
            visa_table = visadevlist;
            
            % Handle case where no VISA devices are detected by the OS
            if isempty(visa_table)
                app.SR830DropDown.Items = {'No Device Found'};
                app.SR844DropDown.Items = {'No Device Found'};
                return;
            end
            
            % Extract Resource Names (e.g., 'GPIB0::8::INSTR') and Model identifiers from the table
            addresses = visa_table.ResourceName;
            models = visa_table.Model;
            
            % Initialize a cell array to store formatted strings for the UI dropdown
            filteredItems = {};
            
            % Iterate through all detected models to find supported Stanford Research (SR) Lock-ins
            for i = 1:length(models)
                currentModel = char(models(i));
                
                % Filter for specific supported models: SR860, SR830, or SR844
                if contains(currentModel, 'SR860') || ...
                   contains(currentModel, 'SR830') || ...
                   contains(currentModel, 'SR844')
                    
                    % Create a user-friendly display string: "Address (Model)"
                    % Example: "GPIB0::8::INSTR (SR830)"
                    displayName = sprintf('%s (%s)', addresses{i}, currentModel);
                    filteredItems{end+1} = displayName; %#ok<AGROW>
                end
            end
            
            % Update UI components if supported instruments are found
            if ~isempty(filteredItems)
                % Populate the visible text in the DropDown menus
                app.SR830DropDown.Items = filteredItems;
                app.SR844DropDown.Items = filteredItems;
                
                % Map the 'ItemsData' to raw VISA addresses for backend communication
                app.SR830DropDown.ItemsData = addresses(ismember(addresses, extractBefore(filteredItems, " (")));
                app.SR844DropDown.ItemsData = addresses(ismember(addresses, extractBefore(filteredItems, " (")));
            else
                % Handle case where VISA devices exist but none match the supported SR models
                app.SR830DropDown.Items = {'No Supported Lock-in Found'};
                app.SR844DropDown.Items = {'No Supported Lock-in Found'};
            end
        end

        function prepare(app)
            [app.X,app.Y]=ndgrid(1:app.res(1), 1:app.res(2)); 
            [t,r]=cart2pol(app.X-app.C(1),app.Y-app.C(2));
            [app.THETA, app.R]=ndgrid(linspace(min(t(:)), max(t(:)), 501), linspace(min(r(:)), max(r(:)), 601));
            [app.Xi, app.Yi]=pol2cart(app.THETA, app.R);
            
            % Prepare Axes2
            cla(app.UIAxes2);
            app.hImage2=imagesc(app.UIAxes2, app.THETA(:,1), app.R(1,:), zeros(size(app.THETA))');
            set(app.UIAxes2, 'YDir', 'normal');
            axis(app.UIAxes2, 'tight');
        end
        
        
        function process(app,f)
            if length(app.C)~=2
                return;
            end
            
            if size(f,3)~=1
                f=rgb2gray(f);
            end
            
            finterp = interp2(app.X'-app.C(1), app.Y'-app.C(2), double(f), app.Xi, app.Yi);
%             h=pcolor(app.UIAxes2, app.THETA, app.R, finterp);
%             set(h, 'EdgeColor', 'none');
            app.hImage2.CData = finterp';
            
%             if length(app.rsel)==2 && all(app.rsel>=app.R(1,1)) && all(app.rsel<=app.R(1,end))
            if isa(app.rect_th, 'images.roi.Rectangle') && ~isempty(app.rect_th.Position)
                % Detect theta movement
                rect = app.rect_th.Position;
                rect(3:4)=rect(1:2)+rect(3:4); % [xmin,ymin,xmax,ymax]
            
                [~,ind1]=min(abs(app.THETA(:,1)-rect(1)));
                [~,ind2]=min(abs(app.THETA(:,1)-rect(3)));
                [~,ind3]=min(abs(app.R(1,:)-rect(2)));
                [~,ind4]=min(abs(app.R(1,:)-rect(4)));
                
                f2 = finterp(ind1:ind2, ind3:ind4);
                Iavg = mean(f2, 2);
                
                th = app.THETA(~isnan(Iavg), 1);
                iavg = Iavg(~isnan(Iavg));
                
                if length(iavg)>4
                
                    plot(app.UIAxes3, th*180/pi, iavg, '.-');
                    axis(app.UIAxes3, 'tight');
                    
                    ift = fftshift(fft(iavg));
                    
                    [~,ind]=sort(abs(ift),'descend');
                    phase = angle(ift(ind(2)));
                    freq = abs(ind(3)-ind(2))/2*1/(th(2)-th(1))/length(th);
                    [~,ind]=max(iavg);
                    peak=th(ind);
                    
                                    
                    app.phase_th(end+1) = phase;
                    
                    % Convert to actual angle
                    if app.AlignpeakCheckBox.Value
                        if isempty(app.disp_th)
                            % First point, define as same as peak
                            theta = peak;
                        else
                            % peak + (phase-phase0)/2/pi/freq
                            % Select the same branch as the first point
                            tt = phase + app.disp_th(1)/180*pi*2*pi*freq - app.phase_th(1); % Shifted phase
                            theta = peak + (mod(tt - peak*2*pi*freq + pi, 2*pi) - pi) / (2*pi*freq); % Make sure it's in same winding as peak
                        end
                    else
                        theta = (mod(phase+pi,2*pi)-pi)/(2*pi*freq);
                    end
                    app.disp_th(end+1) = theta/pi*180; % Convert to degrees
                    
                    % Plot
                    plot(app.UIAxes5, app.disp_th);
                end
            end
            
            % Radial extraction
            if isa(app.rect_r, 'images.roi.Rectangle') && ~isempty(app.rect_r.Position)
                rect = app.rect_r.Position;
                rect(3:4)=rect(1:2)+rect(3:4); % [xmin,ymin,xmax,ymax]
                
                [~,ind1]=min(abs(app.THETA(:,1)-rect(1)));
                [~,ind2]=min(abs(app.THETA(:,1)-rect(3)));
                [~,ind3]=min(abs(app.R(1,:)-rect(2)));
                [~,ind4]=min(abs(app.R(1,:)-rect(4)));
                
                f2 = finterp(ind1:ind2, ind3:ind4);
                Iavg = mean(f2, 1);
                
                r = app.R(1, ~isnan(Iavg));
                iavg = Iavg(~isnan(Iavg));
                
                if length(iavg)>4
                
                    plot(app.UIAxes6, r, iavg, '.-');
                    axis(app.UIAxes6, 'tight');
                    
                    ift = fftshift(fft(iavg));
                    
                    [~,ind]=sort(abs(ift),'descend');
                    phase = angle(ift(ind(2)));
                    freq = abs(ind(3)-ind(2))/2*1/(r(2)-r(1))/length(r);
                    
                    app.phase_r(end+1) = phase;
                    
                    app.disp_r(end+1) = (mod(phase-app.phase_r(1)+pi, 2*pi) - pi) / (2*pi) * app.radialpitchumEditField.Value;
                    
                    % Plot
                    plot(app.UIAxes4, app.disp_r);
                end
            end
        end
        
        function updateVoltages(app)
            while(app.LockIn.BytesAvailable>0)
                fread(app.LockIn, app.LockIn.BytesAvailable);
            end
            app.V1Gauge.Value=str2double(query(app.LockIn, 'auxv? 1', '%s\n', '%s\n')) * app.volt_amp_ratio;
            app.V2Gauge.Value=str2double(query(app.LockIn, 'auxv? 2', '%s\n', '%s\n')) * app.volt_amp_ratio;
            app.V3Gauge.Value=str2double(query(app.LockIn, 'auxv? 3', '%s\n', '%s\n')) * app.volt_amp_ratio;
            app.VzGauge.Value=str2double(query(app.LockIn, 'auxv? 4', '%s\n', '%s\n')) * app.volt_amp_ratio;
        end
    
        function setSteppingEnable(app, en)
            app.ZeroStepsButton.Enable=en;
            app.Minus18Button.Enable=en;
            app.MinusHalfButton.Enable=en;
            app.MinusQuarterButton.Enable=en;
            app.Plus18Button.Enable=en;
            app.PlusHalfButton.Enable=en;
            app.PlusQuarterButton.Enable=en;
        end
        
        function setPositioningEnable(app, en)
            app.GOButton.Enable = en;
        end

        function setPositioningEnable_2(app, en)
            app.SCANButton.Enable = en;
        end

        function setPositioningEnable_3(app, en)
            app.GOButton_2.Enable = en;
        end
        
        function setVoltages(app, volt, phase)
            if ~isa(app.LockIn, 'visa') || ~strcmp(app.LockIn.Status,'open')
                errordlg('USB connection not established');
                return;
            end
            fprintf(app.LockIn, sprintf('auxv 1,%f\n', volt/app.volt_amp_ratio*MEMStepper.vphase(phase/180*pi, 0)));
            fprintf(app.LockIn, sprintf('auxv 2,%f\n', volt/app.volt_amp_ratio*MEMStepper.vphase(phase/180*pi, -2*pi/3)));
            fprintf(app.LockIn, sprintf('auxv 3,%f\n', volt/app.volt_amp_ratio*MEMStepper.vphase(phase/180*pi, 2*pi/3)));
            
            app.updateVoltages();
        end
        
        function setVz(app, vz)
            fprintf(app.LockIn, sprintf('auxv 4,%f\n', vz/app.volt_amp_ratio));
        end

        function alignReference(app, snap, ref, mag, minscale)
            % Scale snapshot to be on same scale as reference
            snap = imresize(snap, mag);
%             ref = imcomplement(ref);
                        
            % Multiscale rotation search to save time
            scale = [12, 10, 8, 6, 5, 4];
            rotrange = [-180, 180]; % Initial rotation range
            rotpts = 13; % Search points within the range for each scale
            
            figure(77907);
            ax=subplot(1,2,2);
            imshowpair(imgradient(snap), imgradient(ref));
            
            sc_final=[];
            first=true;
            
            for sc = scale
                fprintf('Scale=/%d, resolution=%f deg, ', sc, diff(rotrange)/(rotpts-1));
                % Downsample
                snap_sc = imresize(snap, 1/sc);
                ref_sc = imresize(ref, 1/sc);
                
                if sc>1 && (min(size(ref_sc))<20 || min(size(snap_sc))<20)
                    % Scale too coarse, skip the scale
                    fprintf('skipped\n');
                    rotpts = rotpts * 2;
                    continue;
                end
                
                if (max(size(snap_sc))>500)
                    break;
                end
                
                if sc < minscale
                    fprintf('scale too small\n');
                    break;
                end
                
                % Take gradient and binarize, then map to -1/1
%                 snap_sc = snap_sc;
%                 ref_sc =  ref_sc;
                
%                 snap_sc = int8(imbinarize(uint8(imgradient(snap_sc))))*127;
%                 ref_sc = int8(imbinarize(uint8(imgradient(ref_sc))))*127;

%                   snap_sc = imbinarize(snap_sc);
%                   ref_sc = imbinarize(ref_sc);
%                 snap_sc = snap_sc - mean(snap_sc(:));
%                 ref_sc = ref_sc - mean(ref_sc(:));
%                 snap_sc = snap_sc - int8(mean(snap_sc(:)));
%                 ref_sc = ref_sc - int8(mean(ref_sc(:)));
                snap_sc = imgradient(snap_sc);
                ref_sc = imgradient(ref_sc);
                
%                 PSF = fspecial('gaussian',5,1); % create PSF
%                 snap_sc = imfilter(snap_sc,PSF,'replicate','same','conv');
%                 ref_sc = imfilter(ref_sc,PSF,'replicate','same','conv');
                
%                 snap_sc = imbinarize(uint8(snap_sc));
%                 ref_sc = imbinarize(uint8(ref_sc));
                
                bestrot=[];
                bestcorr=0;
                bestsnaprot=[];
                if first
                    % Double points for first range
                    range = linspace(rotrange(1), rotrange(2), 2*rotpts);
                else
                    range = linspace(rotrange(1), rotrange(2), rotpts);
                end
                
                first=false;
                for rot = range
                    rot
                    rot_mat = [cosd(rot) -sind(rot);sind(rot) cosd(rot)];
                    tform = rigid2d(rot_mat, [0 0]);
                    snap_sc_rot = imwarp(snap_sc, tform, 'FillValues', 0);
                    
%                     snap_sc_rot = imrotate(snap_sc, rot, "bilinear", "crop");
                    if app.UseGPUCheckBox.Value
                       corr = gather(xcorr2(gpuArray(single(ref_sc)), gpuArray(single(snap_sc_rot))));
                    else
                       corr = xcorr2(single(ref_sc), single(snap_sc_rot));
                    end
                    
                    ax=subplot(1,2,1);
                    surf(ax, corr, 'EdgeColor', 'none');
                    if max(corr(:)) > max(bestcorr(:))
                        bestcorr=corr;
                        bestrot=rot;
                        bestsnaprot = snap_sc_rot;
                    end
                    drawnow;
                end
                
                fprintf('best=%f deg\n', bestrot);
                
                % Extract translation vector
                [ssr,snd] = max(bestcorr(:));
                [ypeak,xpeak] = ind2sub(size(bestcorr),snd);
                offset = [(xpeak-size(bestsnaprot,2)) (ypeak-size(bestsnaprot,1))];
                
                
                % Create coordinate system
                Rsnap = imref2d(size(bestsnaprot),sc,sc);
                Rref = imref2d(size(ref_sc),sc,sc);
                
                [bestsnaprot_t, Rsnap_t]=imtranslate(bestsnaprot, Rsnap, offset*sc, 'FillValues', 0, 'OutputView', 'full');
                
                
                ax=subplot(1,2,2);
                imshowpair(bestsnaprot_t, Rsnap_t, ref_sc, Rref);
                
                % Refile rotation angle search range
                % Shrink range by 4
                range = 0.25*(rotrange(2)-rotrange(1));
                rotrange = [bestrot-range, bestrot+range];

                sc_final = sc;
                offset_final = offset*sc;
                
%                 [~,ind]=min(abs(range-bestrot));
%                 rotrange = [range(max(ind-2,1)), range(min(ind+2,end))];
            end
            fprintf('Done\n');
            
            if isempty(sc_final)
                % Nothing has been done!
                return;
            end
            
        
            % Now correct the full image
            
            
            % world offset
            % Note that imwarp rotates wrt top left corner
            % whereas imrotate rotates around center
            
            rot_w = [cosd(bestrot) -sind(bestrot);sind(bestrot) cosd(bestrot)];
            tform_w = rigid2d(rot_w, [0 0]);
            snap_aligned =imwarp(snap, tform_w, 'FillValues', 0);
            
            
%             offset_w = [size(snap,2)/2  size(snap,1)/2]*(eye(2)-rot_w) + offset_final;
            
            
            
            Rsnap = imref2d(size(snap_aligned),1,1);
            Rref = imref2d(size(ref),1,1);
            
            [snap_aligned, Rsnap] = imtranslate(snap_aligned, Rsnap, offset_final, 'FillValues', 0, 'OutputView', 'full');
            
            
            
%             snap_aligned = imrotate(snap, bestrot, "bilinear", "crop");
            
            ax=subplot(1,2,2);
            imshowpair(snap_aligned, Rsnap, ref, Rref);
        end

        % function v = readLockin(app)
        %     if isa(app.LockIn_Read, 'visa') && strcmp(app.LockIn_Read.Status,'open')
        %         try
        %             % 1=X, 2=Y, 3=R, 4=Theta
        %             channelIdx = app.ReadChannelDropDown.Value;
        %             switch channelIdx
        %             case 'X'
        %                 channelIdx = 1;
        %             case 'Y'
        %                 channelIdx = 2;
        %             case 'R'
        %                 channelIdx = 3;
        %             case 'Theta'
        %                 channelIdx = 4;
        %             otherwise
        %                 channelIdx = 3; % default
        %             end
        %             cmd = sprintf('OUTP? %d', channelIdx);
        % 
        %             raw_str = query(app.LockIn_Read, cmd); 
        % 
        %             if isempty(raw_str) || strcmp(raw_str, "")
        %                 fprintf('Warning: Received empty string, retrying Channel %d...\n', channelIdx);
        %                 pause(0.1);
        %                 raw_str = query(app.LockIn_Read, cmd);
        %             end
        % 
        %             clean_str = regexp(raw_str, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match');
        %             if ~isempty(clean_str)
        %                 v = str2double(clean_str{1});
        % 
        %                 if channelIdx == 4
        %                     app.CapValue.Text = sprintf("%+9.3f deg", v);
        %                 else
        %                     app.CapValue.Text = sprintf("%+9.3f uV", v*1e6);
        %                 end
        % 
        %                 app.capData(end+1) = v;
        %                 if strcmp(app.StartStopCapButton.Text{1},"Stop")
        %                     plot(app.AxesCapReadout, app.capData, '-');
        %                 end
        %             else
        %                 v = 0;
        %                 app.CapValue.Text = "No Data";
        %             end
        % 
        %         catch ME
        %             fprintf('Error during Channel %d query: %s\n', channelIdx, ME.message);
        %             v = 0;
        %         end
        %     else
        %         v = 0;
        %     end
        % end

        % Timer Update: stable version
        function v = readLockin(app)

            v = 0;
        
            if ~isvalid(app) 
                return; end
            if isempty(app.LockIn_Read)
                return; end
            if ~isa(app.LockIn_Read, 'visa')
                return; end
            if ~strcmp(app.LockIn_Read.Status,'open')
                return;
            end
        
            try
                channelStr = app.ReadChannelDropDown.Value;
                switch channelStr
                    case 'X'
                        channelIdx = 1;
                    case 'Y'
                        channelIdx = 2;
                    case 'R'
                        channelIdx = 3;
                    case 'Theta'
                        channelIdx = 4;
                    otherwise
                        channelIdx = 3; % default: R
                end
        
                cmd = sprintf('OUTP? %d', channelIdx);
                raw_str = query(app.LockIn_Read, cmd);
                if isempty(raw_str)
                    return; end
        
                clean_str = regexp(raw_str, ...
                    '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match');
                if isempty(clean_str)
                    return; end
        
                v = str2double(clean_str{1});
        
                % ---- Check ----
                if ~isvalid(app)
                    return; end
        
                if channelIdx == 4
                    app.CapValue.Text = sprintf("%+9.3f deg", v);
                else
                    app.CapValue.Text = sprintf("%+9.3f uV", v*1e6);
                end
        
                app.capData(end+1) = v;
        
                if strcmp(app.StartStopCapButton.Text{1},"Stop")
                    plot(app.AxesCapReadout, app.capData, '-');
                end
        
            catch
                fprintf('Error during Channel %d query: %s\n', channelIdx, ME.message)
                v = 0;
            end
        end


        % function timerCallback(app, src, evt)
        %     app.readLockin();
        % end

        % Timer Update: stable version
        function timerCallback(app)
        
            if ~isvalid(app) || strcmp(app.tim.Running,'off')
                return;
            end
        
            try
                app.readLockin();
            catch ME
                disp(ME.message)
            end
        
        end


        function addCapFittingPoint(app)
            if app.TabGroup3.SelectedTab == app.SteppingTab
                pause(app.WaitTimesEditField.Value);
                % Prevent data racing
                app.tim.stop();
                app.readLockin();
                if strcmp(app.StartStopCapButton.Text{1},"Stop")
                    app.tim.start();
                end
                x = app.StepsEditField.Value;
                y = app.capData(end);
                app.capFitData(end+1,:) = [x,y];
                plot(app.AxesCapFitting, app.capFitData(:,1), app.capFitData(:,2), '.', 'MarkerSize',5);
            end
        end

        function addCapFittingPointwithPolyfit(app, wait_time)
            if app.TabGroup3.SelectedTab == app.PositioningFeedbackTab
                % pause(app.WaitTimescaponlyEditField.Value); % note: setup value in org code
                pause(wait_time); % Prevent data racing

                app.tim.stop();
                app.readLockin();
                if strcmp(app.StartStopCapButton.Text{1},"Stop")
                    app.tim.start();
                end
                if ~isempty(app.capData)
                    x = app.CurrentStepEditField.Value;
                    y = app.capData(end);
                    app.capFitData(end+1,:) = [x,y];

                    % --- Old Plot with Fitting Line ---
                    % plot(app.AxesCapFitting, app.capFitData(:,1), app.capFitData(:,2), '.', 'MarkerSize',5);
    
                    % --- Updated Plot with Fitting Line ---
                    x_data = app.capFitData(:,1);
                    y_data = app.capFitData(:,2);
                    
                    % Perform linear fitting: y = p(1)*x + p(2)
                    if length(x_data) >= 2
                        p = polyfit(x_data, y_data, 1);
                        app.k_fit = p(1);
                        app.b_fit = p(2);
                        
                        % Generate y values for the fitting line
                        y_fit = polyval(p, x_data);
                        
                        % 1. Plot Raw Data: Red open circles ('or')
                        plot(app.AxesCapFitting, x_data, y_data, 'o', ...
                            'Color', 'r', ...
                            'MarkerSize', 7, ...
                            'MarkerFaceColor', 'none', ...
                            'LineWidth', 1); 
                        
                        hold(app.AxesCapFitting, 'on'); 
                        
                        % 2. Plot Fitting Line: Black solid line ('-k')
                        plot(app.AxesCapFitting, x_data, y_fit, '-k', 'LineWidth', 1.5);
                        
                        hold(app.AxesCapFitting, 'off'); 
                        
                        % 3. Update Legend with k and b values
                        legend_str = sprintf('%.2e*x+%.2e', app.k_fit, app.b_fit);
                        legend(app.AxesCapFitting, 'Raw Data', legend_str, 'Location', 'northwest');
                    else
                        % If not enough data points, just plot the red hollow circles
                        plot(app.AxesCapFitting, x_data, y_data, 'o', ...
                            'Color', 'r', ...
                            'MarkerSize', 7, ...
                            'MarkerFaceColor', 'none');
                    end
                end
            end
        end

        function goTo(app, volt, target, microstep, speed)
            % NOTE: microstep should not be '1'!
            curr_step = app.StepsPhyEditField.Value;
            total_steps = round(abs(target-curr_step)*(microstep));
            wait_time = 1/(speed * microstep);
            step_inc = 1/microstep*sign(target-curr_step);
            
            for i=1:total_steps
                curr_step = curr_step + step_inc;
                phase = curr_step * 120;
                
                app.setVoltages(volt, phase);
                % fprintf("goTo setVoltages(%.3e, %.3e) @%s [app.capData(end)=%.3e]\n", volt, phase, datetime('now','Format','mm:ss'), app.capData(end));
                pause(wait_time);
                if isa(app.video, 'videoinput') && isrunning(app.video)
                    trigger(app.video);
                end
                app.StepsPhyEditField.Value = curr_step;
                if app.StepsCapZero > -100
                    app.StepsCapEditField.Value = curr_step - app.StepsCapZero;
                end
            end
        end

        function goTo_fit(app, volt, start_pos, target_pos, microstep, speed, update_k_b)
            dist = target_pos - start_pos;
            total_steps = round(abs(dist) * microstep);
            
            if total_steps == 0
                return;
            end
            
            % Calculate the signed increment per microstep and the sampling interval
            % wait_time ensures the Lock-in amplifier's integration (Time Constant) has settled.
            step_inc = (1 / microstep) * sign(dist);
            % wait_time = max(0.5, 1/(speed * microstep));
            wait_time = max(app.default_wait_time, 1/speed);
            fprintf("[LOG] goTo_fit (scan version):\twait_time=%.4f s\n", wait_time);
            
            % --- Execution Loop ---
            for i = 0:total_steps
                if app.positioningSTOP
                    break
                elseif app.SCANButton.Enable % Check if the system is in scanning mode

                    curr_step = start_pos + i * step_inc;
                    
                    phase = curr_step * 120;
                    app.setVoltages(volt, phase);
                    app.CurrentStepEditField.Value = curr_step;
                    
                    % Updates the Sensitivity (k_fit) and Offset (b_fit) for subsequent feedback.
                    if update_k_b
                        app.addCapFittingPointwithPolyfit(wait_time);
                    else
                        pause(wait_time);
                    end
              
                    % Note: Video triggering can be enabled here for synchronized visual inspection
                    % if isa(app.video, 'videoinput') && isrunning(app.video)
                    %     trigger(app.video);
                    % end
                end
            end
        end


        % function goTo_fb(app, volt, target_pos, microstep2)
        %     % This function is designed for scan best PID sets
        % 
        %     % --- 1. Control & Target Initialization ---                
        %     % --- Start overall timer ---
        %     total_timer = tic; 
        % 
        %     navg = max(1, app.FilterWinEditField.Value);
        %     step_tolerance = abs(app.StepMinLossEditField.Value);
        %     demanded_x = target_pos;
        %     y_theoretical = app.k_fit * demanded_x + app.b_fit;
        % 
        %     % speed = app.SpeedstepsEditField_2.Value;
        %     wait_time = max(app.default_wait_time, app.PIDwaitsEditField.Value);
        %     if wait_time ~= (1/app.SpeedstepsEditField_2.Value)
        %         fprintf("[WARING] goTo_fb (scan version):\tmodeling sampling time is not equal to feedback sampling time!\n")
        %     end
        %     fprintf("[LOG] goTo_fb (scan version):\twait_time=%.4f s. (>= default %.4f s)\n", wait_time, app.default_wait_time);
        % 
        %     % --- 2. Define Parameter Sweep Combinations ---
        %     % Kp_list = [0.005, 0.15, 0.3, 0.45];
        %     % Ki_list = [0.001, 0.002, 0.008, 0.01];
        %     % Kd_list = [0.01, 0.02, 0.08, 0.1];
        %     Kp_list = [0.3, 0.35, 0.4, 0.45, 0.5]; 
        %     Ki_list = [0.002, 0.0005, 0.0015]; 
        %     Kd_list = [0.05, 0.1, 0.15, 0.2];
        %     pts_per_config = 300; 
        % 
        %     % Generate combination matrix (N x 3)
        %     [KP, KI, KD] = meshgrid(Kp_list, Ki_list, Kd_list);
        %     configs = [KP(:), KI(:), KD(:)]; 
        %     zero_config = [0, 0, 0];
        %     configs = [configs; zero_config];
        % 
        %     % % Define combination matrix
        %     % specific_configs = [
        %     %     0.3,   0.002, 0.02;    % best (G3)
        %     %     0.3,   0.005, 0.02;    % +I
        %     %     0.3,   0.01,  0.02;    % ++I
        %     %     0.3,   0.002, 0.05;    % +D
        %     %     0.3,   0.002, 0.09;    % ++D
        %     %     0.005, 0.01,  0.01;     % vibrate
        %     % ];
        %     % 
        %     % configs = [specific_configs; zero_config];
        % 
        %     num_configs = size(configs, 1);
        % 
        %     % --- 3. Memory & Log Initialization ---
        %     app.capFitData = zeros(num_configs * pts_per_config, 2); 
        %     pt_count = 0;
        % 
        %     % Open log file
        %     fid = fopen(app.LogFileName, 'a');
        %     if fid == -1, error('[LOG] Could not create log file.'); end
        % 
        %     % Force initial position for the very first group
        %     force_org_step = 5;
        %     fprintf("[LOG] Forcing initial x=%d before starting sweeps.\n", force_org_step);
        %     app.setVoltages(volt, force_org_step * 120);
        %     app.CurrentStepEditField.Value = force_org_step;
        % 
        %     % --- 4. Main Sweep Loop ---
        %     for c = 1:num_configs
        %         % Apply current PID parameters
        %         Kp = configs(c, 1);
        %         Ki = configs(c, 2);
        %         Kd = configs(c, 3);
        %         app.KpEditField.Value = Kp;
        %         app.KiEditField.Value = Ki;
        %         app.KdEditField.Value = Kd;
        % 
        %         % Reset PID state when switching parameters
        %         integral_error = 0;
        %         last_error_step = 0;
        % 
        %         % MODIFICATION: Initialize PID pause counter
        %         pid_pause_counter = 0;
        % 
        %         dt = wait_time * navg;
        % 
        %         % --- 5. Data Collection (Inner Loop) ---
        %         for s = 1:pts_per_config
        %             if app.positioningSTOP, break; end 
        % 
        %             % --- A. Data Acquisition ---
        %             app.tim.stop(); 
        %             temp_data = zeros(1, navg);
        %             for n = 1:navg
        %                 pause(wait_time); 
        %                 app.readLockin();
        %                 temp_data(n) = app.capData(end);
        %             end
        %             if strcmp(app.StartStopCapButton.Text{1}, "Stop") && strcmp(app.tim.Running, 'off')
        %                 app.tim.start(); 
        %             end
        % 
        %             % --- B. Error Calculation ---
        %             y_actual_filtered = mean(temp_data);
        %             current_x = app.CurrentStepEditField.Value;
        %             cap_loss = y_theoretical - y_actual_filtered; 
        %             error_step = cap_loss / app.k_fit;
        % 
        %             % --- C. PID Control / Force Step Logic ---    
        %             pt_count = pt_count + 1; % Total counter
        %             fprintf("[DEBUG] pt_count=%d, pid_pause_counter=%d [mod(pt_count, 100)=%d]\n", pt_count, pid_pause_counter, mod(pt_count, 100));
        %             % Every pts_per_config points, force step and pause PID
        %             if (mod(pt_count, pts_per_config) == 1) && (Kp ~=0)
        %                 % fprintf("[DEBUG] pid_pause_counter=%d\n", pid_pause_counter);
        %                 new_x = force_org_step;
        %                 status_str = 'FORCE STEP';
        % 
        %                 % Reset PID state
        %                 integral_error = 0;
        %                 last_error_step = 0;
        % 
        %                 % Pause PID calculations for the next 10 points
        %                 pid_pause_counter = 10;
        % 
        %             elseif pid_pause_counter > 0
        %                 % --- PID Paused Phase (P=I=D=0) ---
        %                 P_out = 0; I_out = 0; D_out = 0; dynamic_step = 0;
        %                 new_x = current_x; % Hold position
        %                 status_str = 'PID PAUSED';
        %                 pid_pause_counter = pid_pause_counter - 1;
        % 
        %             elseif microstep2 > 256
        %                 % --- Normal PID Adjustment ---
        %                 P_out = Kp * error_step;
        %                 if abs(error_step) > step_tolerance
        %                     integral_error = integral_error + (error_step * dt);
        %                 end
        %                 I_out = Ki * integral_error;
        %                 derivative_error = (error_step - last_error_step) / dt;
        %                 D_out = Kd * derivative_error;
        %                 dynamic_step = P_out + I_out + D_out;
        %                 last_error_step = error_step;
        % 
        %                 max_safe_step = 10; 
        %                 if abs(dynamic_step) > max_safe_step
        %                     dynamic_step = sign(dynamic_step) * max_safe_step;
        %                 end
        % 
        %                 new_x = current_x + dynamic_step;
        %                 status_str = 'PID ADJ';
        %             else
        %                 % --- Non-PID Adjustment ---
        %                 P_out = 0; I_out = 0; D_out = 0; dynamic_step = 0;
        %                 if abs(error_step) <= step_tolerance
        %                     status_str = 'STABLE';
        %                     new_x = current_x; 
        %                 else
        %                     status_str = 'ADJUSTING';
        %                     new_x = current_x + sign(error_step) * (1 / microstep2);
        %                 end
        %             end
        % 
        %             % --- D. Execute Hardware Move ---
        %             if ~strcmp(status_str, 'STABLE')
        %                 app.setVoltages(volt, new_x * 120);
        %                 app.CurrentStepEditField.Value = new_x;
        %                 current_x = new_x; % Update local variable for recording
        %             end
        % 
        %             % --- E. Recording & UI Update ---
        %             app.capFitData(pt_count, :) = [current_x, y_actual_filtered];
        % 
        %             % [ ... UI Plotting logic remains the same ... ]
        %             if isprop(app, 'plotCapFBloss') && app.plotCapFBloss
        %                 ax_fb = app.AxesCapReadout;
        %                 hFBTrace = findobj(ax_fb, 'Tag', 'fb_loss_trace');
        %                 if isempty(hFBTrace) || ~isgraphics(hFBTrace)
        %                     cla(ax_fb); hold(ax_fb, 'on');
        %                     yline(ax_fb, y_theoretical, '--r', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        %                     hFBTrace = plot(ax_fb, NaN, NaN, '-o', 'Color', [0 0.5 0], ...
        %                         'LineWidth', 1.2, 'MarkerSize', 5, ...
        %                         'MarkerFaceColor', [0 0.8 0], 'Tag', 'fb_loss_trace');
        %                     hold(ax_fb, 'off'); grid(ax_fb, 'on');
        %                     ylabel(ax_fb, 'Capacitance'); xlabel(ax_fb, 'Sample Index');
        %                 end
        %                 currentX = hFBTrace.XData; currentY = hFBTrace.YData;
        %                 newXData = [currentX(~isnan(currentX)), pt_count];
        %                 newYData = [currentY(~isnan(currentY)), y_actual_filtered];
        %                 if length(newXData) > 20
        %                     newXData = newXData(end-19:end); newYData = newYData(end-19:end);
        %                 end
        %                 set(hFBTrace, 'XData', newXData, 'YData', newYData);
        %                 if abs(app.k_fit) > 0
        %                     equiv_step_loss = -cap_loss / app.k_fit;
        %                 else
        %                     equiv_step_loss = 0;
        %                 end
        %                 legStr = sprintf('C_{act}\nΔC:%.2e\nΔStep:%.3f', cap_loss, equiv_step_loss);
        %                 legend(ax_fb, legStr, 'Location', 'northoutside', 'FontSize', 9, 'FontWeight', 'bold');
        %             end
        % 
        %             % --- F. Fixed-Width Formatted Logging ---
        %             step_loss_from_cap = cap_loss / app.k_fit;
        %             fmt = ['Cap(theor-act): [%11.3e - %11.3e = %11.3e] | ', ...
        %                    'Step Loss (=Cap/k_fit): [%7.3f] | ', ...
        %                    'Gains(P,I,D): [%.2e, %.2e, %.2e] | ', ...
        %                    'PID_Out(P+I+D): [%7.4f + %7.4f + %7.4f = %+10.4f] @ %s [%-15s] @ Step [%.4e/%.4e]\n'];
        %             % Update logging arguments based on status
        %             if strcmp(status_str, 'FORCE STEP') || strcmp(status_str, 'PID PAUSED')
        %                 log_args = {y_theoretical, y_actual_filtered, cap_loss, ...
        %                             step_loss_from_cap, Kp, Ki, Kd, ...
        %                             0, 0, 0, 0, ... % Outputs are effectively zero or forced
        %                             datetime('now','Format','HH:mm:ss'), status_str, current_x, demanded_x};
        %             else
        %                 log_args = {y_theoretical, y_actual_filtered, cap_loss, ...
        %                             step_loss_from_cap, Kp, Ki, Kd, ...
        %                             P_out, I_out, D_out, dynamic_step, ...
        %                             datetime('now','Format','HH:mm:ss'), status_str, current_x, demanded_x};
        %             end
        % 
        %             fprintf(fmt, log_args{:});      % Command Window
        %             fprintf(fid, fmt, log_args{:}); % Log File
        %             drawnow limitrate; 
        %         end
        %         if app.positioningSTOP, break; end
        %     end
        % 
        %     % --- 6. Cleanup ---
        %     fclose(fid);
        %     app.capFitData = app.capFitData(1:pt_count, :);
        %     total_time = toc(total_timer);
        %     fprintf("\n[DONE] Auto-sweep completed. Total points: %d. Total Time: %.2f seconds\n", pt_count, total_time);
        % end
        % 

        function goTo_fb(app, volt, target_pos, microstep2)
            % GOTO_FB (PID): High-precision PID feedback with optimized memory and formatted logging.

            % --- 1. Control & Target Initialization ---
            navg = max(1, app.FilterWinEditField.Value);
            step_tolerance = abs(app.StepMinLossEditField.Value);
            demanded_x = target_pos; % Theoretical target position
            y_theoretical = app.k_fit * demanded_x + app.b_fit; % Theoretical target capacitance

            % Sampling and system settling interval
            % wait_time = max(0.01, 1/(speed * microstep2)); 
            % wait_time = max(app.default_wait_time, 1/speed);
            wait_time = max(app.default_wait_time, app.PIDwaitsEditField.Value);
            if wait_time ~= (1/app.SpeedstepsEditField_2.Value)
                fprintf("[WARING] goTo_fb:\tmodeling sampling time is not equal to feedback sampling time!\n")
            end
            fprintf("[LOG] goTo_fb:\twait_time=%.4f s. (>= default %.4f s)\n", wait_time, app.default_wait_time);

            % --- 2. PID Coefficients ---
            max_safe_step = 10;
            % Kp = app.KpEditField.Value;      % Proportional Gain
            % Ki = app.KiEditField.Value;      % Integral Gain
            % Kd = app.KdEditField.Value;      % Derivative Gain

            % PID State Variables
            integral_error = 0;
            last_error_step = 0;
            dt = wait_time * navg; % Effective loop time step

            % --- 3. Performance Optimization: Memory Pre-allocation ---
            % Prevents the O(N) overhead of dynamic array resizing
            max_pts = 50000; 
            app.capFitData = zeros(max_pts, 2); 
            pt_count = 0;

            % --- 4. UI/Graphics Setup (Full Recovery with Optimized Handles) ---
            % Log file
            fid = fopen(app.LogFileName, 'a');
            if fid == -1
                error('[LOG] Could not create log file.');
            end

            % Prepare Plot
            if strcmp(app.LogPlotSwitch.Value, "On")
                ax = app.AxesCapFitting;
                % A. Cleanup: Clear redundant data points from previous sessions
                hOldPoints = findobj(ax, 'Type', 'line', 'Color', 'r', 'Marker', 'o');
                if ~isempty(hOldPoints), delete(hOldPoints); end
                % B. Handle Management: Seek existing specialized plot handles
                hTrace = findobj(ax, 'Tag', 'trace_line');
                hPoint = findobj(ax, 'Tag', 'current_pt');
                hTheo  = findobj(ax, 'Tag', 'target_pt');
                % C. Initialization: Create plot objects if they don't exist
                if isempty(hTrace)
                    hold(ax, 'on');
                    % Target crosshair (Horizontal & Vertical references)
                    yline(ax, y_theoretical, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');
                    xline(ax, demanded_x, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');
                    % Target point (The static goal)
                    hTheo = plot(ax, demanded_x, y_theoretical, 'o', 'Color', 'r', 'MarkerFaceColor', 'b', 'MarkerSize', 8, ...
                        'DisplayName', 'Target Point', 'Tag', 'target_pt');
                    % Trace line (Recent historical path)
                    hTrace = plot(ax, NaN, NaN, 'o-', 'Color', [0 0 0.5], 'MarkerSize', 4, 'DisplayName', 'Trace Path', 'Tag', 'trace_line');
                    % Current position (Real-time marker)
                    hPoint = plot(ax, NaN, NaN, '+r', 'LineWidth', 2, 'MarkerSize', 12, 'DisplayName', 'Current Pos', 'Tag', 'current_pt');
                    hold(ax, 'off');

                    % Restore the legend for user clarity
                    legend(ax, [hTheo, hTrace, hPoint], 'Location', 'northwest', 'FontSize', 7);
                else
                    % If handles exist but target has changed, update Target Point and X/Y lines
                    set(hTheo, 'XData', demanded_x, 'YData', y_theoretical);
                    % (Note: To update yline/xline dynamically you'd need their handles, 
                    % but usually 'cla' or deleting old ones at start is simpler).
                end
                max_dist_x = 0; max_dist_y = 0;
            end

            last_status = "";

            % --- 5. Main PID Feedback Loop ---
            while (~app.positioningSTOP)
                % --- A. Data Acquisition ---
                Kp = app.KpEditField.Value;      % Proportional Gain
                Ki = app.KiEditField.Value;      % Integral Gain
                Kd = app.KdEditField.Value;      % Derivative Gain

                app.tim.stop(); % Pause background timer to avoid VISA conflicts
                temp_data = zeros(1, navg);
                for n = 1:navg
                    pause(wait_time); 
                    app.readLockin();
                    temp_data(n) = app.capData(end);
                end
                if strcmp(app.StartStopCapButton.Text{1}, "Stop")
                    if strcmp(app.tim.Running, 'off'), app.tim.start(); end
                else
                    % fprintf('Timer is already active, skipping start.\n');
                end


                % --- B. Error Calculation (Capacitance Domain) ---
                y_actual_filtered = mean(temp_data);
                current_x = app.CurrentStepEditField.Value;
                cap_loss = y_theoretical - y_actual_filtered; 
                error_step = cap_loss / app.k_fit;

                % --- C. PID Control (Step Domain) ---    
                if microstep2 > 256 % <CASE A> PID Continuous Feedback Control (High-Resolution)
                    % Proportional term
                    P_out = Kp * error_step;

                    % Integral term (with Anti-windup deadzone)
                    % Even if this if-condition isn't met, integral_error already equals 0
                    if abs(error_step) > step_tolerance
                        integral_error = integral_error + (error_step * dt);
                    end
                    I_out = Ki * integral_error;

                    % Derivative term
                    derivative_error = (error_step - last_error_step) / dt;
                    D_out = Kd * derivative_error;

                    % Combined PID Output
                    dynamic_step = P_out + I_out + D_out;
                    last_error_step = error_step;

                    % Safety Clamping 
                    if abs(dynamic_step) > max_safe_step
                        dynamic_step = sign(dynamic_step) * max_safe_step;
                    end

                    new_x = current_x + dynamic_step;
                    status_str = 'PID ADJ';

                else % <CASE B> Fixed Step Logic (Standard Resolution)
                    % We reset PID components to zero for clean logging in this mode
                    P_out = 0; I_out = 0; D_out = 0; dynamic_step = 0;

                    if abs(error_step) <= step_tolerance
                        status_str = 'STABLE';
                        new_x = current_x; 
                    else
                        status_str = 'ADJUSTING';
                        % Direction is the sign of the step error
                        direction = sign(error_step);
                        new_x = current_x + direction * (1 / microstep2);
                    end
                end

                % --- Execute Hardware Move ---
                if ~strcmp(status_str, 'STABLE')
                    app.setVoltages(volt, new_x * 120);
                    app.CurrentStepEditField.Value = new_x;
                    current_x = new_x;
                end

                % --- E. Fast Data Recording ---
                pt_count = pt_count + 1;
                if pt_count <= max_pts
                    app.capFitData(pt_count, :) = [current_x, y_actual_filtered];
                else
                    app.capFitData(end+1, :) = [current_x, y_actual_filtered];
                end

                % --- F. High-Efficiency Visualization ---
                % --- app.AxesCapReadout ---
                if isprop(app, 'plotCapFBloss') && app.plotCapFBloss
                    ax_fb = app.AxesCapReadout; % Target axes for real-time tracking

                    % Find the handle using Tag to avoid redundant plot objects
                    hFBTrace = findobj(ax_fb, 'Tag', 'fb_loss_trace');

                    % If handle is missing (first run or cleared), initialize the plot
                    if isempty(hFBTrace) || ~isgraphics(hFBTrace)
                        cla(ax_fb);
                        hold(ax_fb, 'on');

                        % Draw horizontal dashed line as the setpoint reference
                        yline(ax_fb, y_theoretical, '--r', 'LineWidth', 1.2, 'HandleVisibility', 'off');

                        % Initialize the dynamic trace line with specific markers
                        hFBTrace = plot(ax_fb, NaN, NaN, '-o', 'Color', [0 0.5 0], ...
                            'LineWidth', 1.2, 'MarkerSize', 5, ...
                            'MarkerFaceColor', [0 0.8 0], 'Tag', 'fb_loss_trace');

                        hold(ax_fb, 'off');
                        grid(ax_fb, 'on');
                        % Ensure labels are set for clarity
                        ylabel(ax_fb, 'Capacitance');
                        xlabel(ax_fb, 'Sample Index');
                    end

                    % Keep latest 20 points
                    currentX = hFBTrace.XData;
                    currentY = hFBTrace.YData;

                    % Concatenate new data (filtering out NaNs from initialization)
                    newXData = [currentX(~isnan(currentX)), pt_count];
                    newYData = [currentY(~isnan(currentY)), y_actual_filtered];

                    % Slice the array to keep only the last 20 elements
                    if length(newXData) > 20
                        newXData = newXData(end-19:end);
                        newYData = newYData(end-19:end);
                    end

                    % Update the graphics object directly for high performance
                    set(hFBTrace, 'XData', newXData, 'YData', newYData);

                    % Calculate equivalent positioning error (Step Loss)
                    % Step Loss = Capacitance Loss / Sensitivity (k_fit)
                    if abs(app.k_fit) > 0
                        equiv_step_loss = -cap_loss / app.k_fit;
                    else
                        equiv_step_loss = 0; % Prevent division by zero
                    end

                    % Construct legend string with both Cap and Step error metrics
                    % Using Consolas or Monospace font in UI is recommended for stable digits
                    legStr = sprintf('C_{act}\nΔC:%.2e\nΔStep:%.3f', cap_loss, equiv_step_loss);
                    legend(ax_fb, legStr, 'Location', 'northoutside', 'FontSize', 9, 'FontWeight', 'bold');
                end

                % --- app.AxesCapFitting ---
                if strcmp(app.LogPlotSwitch.Value, "On")
                    trace_idx = max(1, pt_count-14):pt_count;
                    set(hTrace, 'XData', app.capFitData(trace_idx, 1), 'YData', app.capFitData(trace_idx, 2));
                    set(hPoint, 'XData', current_x, 'YData', y_actual_filtered);

                    % Adaptive Zooming
                    max_dist_x = max(max_dist_x, abs(current_x - demanded_x));
                    max_dist_y = max(max_dist_y, abs(y_actual_filtered - y_theoretical));
                    xlim(ax, [demanded_x - 1.5*max_dist_x - 1e-5, demanded_x + 1.5*max_dist_x + 1e-5]);
                    ylim(ax, [y_theoretical - 1.5*max_dist_y - 1e-11, y_theoretical + 1.5*max_dist_y + 1e-11]);

                    if ~strcmp(status_str, last_status)
                        title(ax, sprintf('Status: %s', status_str));
                        last_status = status_str;
                    end
                    drawnow limitrate; 

                    % --- G. Fixed-Width Formatted Logging (Full PID Breakdown) ---
                    % Calculate the equivalent step loss (physical error based on capacitance)
                    % Equivalent Step Loss = Capacitance Loss / Sensitivity (k_fit)
                    step_loss_from_cap = cap_loss / app.k_fit;

                    fmt = ['Cap(theor-act): [%11.3e - %11.3e = %11.3e] | ', ...
                           'Step Loss (=Cap/k_fit): [%7.3f] | ', ...
                           'Gains(P,I,D): [%.2e, %.2e, %.2e] | ', ...
                           'PID_Out(P+I+D): [%7.4f + %7.4f + %7.4f = %+10.4f] @ %s [%-9s] @ Step [%.4e/%.4e]\n'];


                    fprintf(fmt, ...
                            y_theoretical, ...     % Theoretical Cap
                            y_actual_filtered, ... % Actual Cap
                            cap_loss, ...          % Cap Loss
                            step_loss_from_cap, ...% Step Loss (Calculated from Cap error)
                            Kp, Ki, Kd, ...        % CURRENT PID PARAMETERS (Gains)
                            P_out, ...             % Proportional term
                            I_out, ...             % Integral term
                            D_out, ...             % Derivative term
                            dynamic_step, ...      % Final clamped output
                            datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS'), ...
                            status_str, current_x, demanded_x);


                end
            end

            % Clean up: Trim unused pre-allocated rows and close log file
            fclose(fid);
            app.capFitData = app.capFitData(1:pt_count, :);
        end


        
        % Button pushed function (LogCapDataButton) calls
        function logCapData(app)
            % Persistent variable to track the last logged data point index
            persistent lastSavedIdx
            
            % Initialize index to 0 for the first time the button is pushed
            if isempty(lastSavedIdx)
                lastSavedIdx = 0;
                tic;
            end

            t_stamp = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
            t_str = char(t_stamp); % Convert to character array for fprintf

            % 1. Check for new data points in the app.capData array
            currentDataLen = length(app.capData);
            if currentDataLen <= lastSavedIdx
                fprintf('[LOG] No new data available since the last log @ %s\n', t_str);
                return;
            end

            % 2. Slice the array to get only the new data points
            newData = app.capData(lastSavedIdx+1 : currentDataLen);
            
            % 3. Open 'capdata.log' in Append mode ('a')
            fid = fopen('capdata.log', 'a');
            if fid == -1
                uialert(app.UIFigure, 'Failed to open capdata.log', 'File Error');
                return;
            end
            
            try
                % Write a header for this specific logging batch
                fprintf(fid, '\n>> LOG BATCH: %s | Count: %d <<\n', t_str, length(newData));
                fprintf(fid, '--------------------------------------------------\n');
                
                % 4. Write data points with the timestamp
                % Determine if the current channel is degrees or voltage
                isTheta = strcmp(app.ReadChannelDropDown.Value, 'Theta');
                
                for i = 1:length(newData)
                    if isTheta
                        % Format for Angular data
                        fprintf(fid, '[%s] Pt %06d: %+.5e deg\n', t_str, lastSavedIdx + i, newData(i));
                    else
                        % Format for Voltage data (converted to uV)
                        capVolt = newData(i);
                        if ~(isempty(app.k_fit) || app.k_fit == 0)
                            step_theoretical = (capVolt - app.b_fit) / app.k_fit;
                        else
                            step_theoretical = 0;
                        end
                        fprintf(fid, '[%s] Pt %06d: %+.5e V [step=%7.3f]\n', t_str, lastSavedIdx + i, capVolt, step_theoretical);
                    end
                end
                
                fprintf(fid, '--------------------------------------------------\n');
                fclose(fid);
                
                % 5. Update the persistent index to the end of the current array
                lastSavedIdx = currentDataLen;
                fprintf('[LOG] Appended %d points to capdata.log @ %s\n', length(newData), t_str);
                
            catch ME
                % Ensure the file handle is closed in case of an error
                if fid ~= -1, fclose(fid); end
                rethrow(ME);
            end
            fprintf('[LOG] Loged successfully @ %s\n', t_str);
        end

            
        function goTo_fb_scan(app, volt, target, microstep, speed)
            % parameters for scan loop
            scan_loop_step = app.StepsPhyEditField.Value;
            total_steps = round(abs(target-scan_loop_step)*(microstep));
            step_inc = 1/microstep*sign(target-scan_loop_step);
        
            % wait_time_for_loop = 1/(speed * microstep);
            wait_time_for_loop = 1/speed;
            wait_time = max(app.default_wait_time, app.PIDwaitsEditField.Value);
            if wait_time ~= (1/app.SpeedstepsEditField_2.Value)
                fprintf("[WARING] goTo_fb_scan:\tmodeling sampling time is not equal to feedback sampling time!\n")
            end
            fprintf("[LOG] goTo_fb_scan:\tscan loop wait_time=%.4f s. (>= default %.4f s)\n", wait_time, app.default_wait_time);
            
        
            % parameters for PID feedback loop
            navg = max(1, app.FilterWinEditField.Value);
            step_tolerance = abs(app.StepMinLossEditField.Value);
            integral_error = 0;
            last_error_step = 0;
            dt = wait_time * navg;
            fb_loop_step = scan_loop_step;
        
            % feedback log
            fid = fopen(app.LogFileName, 'a');
            if fid == -1, error('[LOG] Could not create log file.'); end

            fprintf("[NOTE] app.logCapData() in use, suggest to push 'StartStopCapButton' first.\n")
            if wait_time_for_loop/wait_time < 1
                fprintf("[WARNING] SCAN wait time (%.3e) should larger than feedback sampling wait time (%.3e)\n", wait_time_for_loop, wait_time);
                total_cnt = 1;
            else
                total_cnt = wait_time_for_loop/wait_time;
            end
            for i=1:total_steps
                scan_loop_step = scan_loop_step + step_inc;
        
                cnt = 0;
                
                fprintf("[LOG] Starting feedback loop with total_cnt=%d for step %d/%d (scan_loop_step=%.3e).\n", total_cnt, i, total_steps, scan_loop_step);
                while cnt < total_cnt
                    if strcmp(app.CloseLoopSwitch.Value, "Off")
                        fprintf("[LOG] Feedback loop terminated.\n\n");
                        break
                    end
                    y_theoretical = app.k_fit * scan_loop_step + app.b_fit; % Theoretical target capacitance
                    Kp = app.KpEditField.Value;
                    Ki = app.KiEditField.Value;
                    Kd = app.KdEditField.Value;
        
                    app.tim.stop(); 
                    temp_data = zeros(1, navg);
                    for n = 1:navg
                        pause(wait_time); 
                        app.readLockin();
                        temp_data(n) = app.capData(end);
                    end
                    if strcmp(app.StartStopCapButton.Text{1}, "Stop")
                        if strcmp(app.tim.Running, 'off'), app.tim.start(); end
                    else
                        % fprintf('Timer is already active, skipping start.\n');
                    end
        
                    y_actual_filtered = mean(temp_data);
                    cap_loss = y_theoretical - y_actual_filtered; 
                    error_step = cap_loss / app.k_fit;
        
                    % PID Calculations
                    P_out = Kp * error_step;
                    if abs(error_step) > step_tolerance
                        integral_error = integral_error + (error_step * dt);
                    end
                    I_out = Ki * integral_error;
                    derivative_error = (error_step - last_error_step) / dt;
                    D_out = Kd * derivative_error;
                    dynamic_step = P_out + I_out + D_out;
                    last_error_step = error_step;
        
                    % Safety Clamping
                    max_safe_step = 1; 
                    if abs(dynamic_step) > max_safe_step
                        dynamic_step = sign(dynamic_step) * max_safe_step;
                    end
        
                    % Update Position
                    new_x = fb_loop_step + dynamic_step;
                    app.setVoltages(volt, new_x * 120);
                    % app.StepsPhyEditField.Value = new_x;
                    fb_loop_step = new_x;
        
                    % Log
                    step_loss_from_cap = cap_loss / app.k_fit;
                    fmt = ['Cap(theor-act): [%11.3e - %11.3e = %11.3e] | ', ...
                           'Step Loss (=Cap/k_fit): [%7.3f] | ', ...
                           'Gains(P,I,D): [%.2e, %.2e, %.2e] | ', ...
                           'PID_Out(P+I+D): [%7.4f + %7.4f + %7.4f = %+10.4f] @ %s @ Step [%.4e/%.4e]\n'];

                    log_args = {y_theoretical, y_actual_filtered, cap_loss, ...
                                step_loss_from_cap, Kp, Ki, Kd, ...
                                P_out, I_out, D_out, dynamic_step, ...
                                datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS'), fb_loop_step, scan_loop_step};

                    fprintf(fmt, log_args{:});      % Command Window
                    fprintf(fid, fmt, log_args{:}); % Log file

                    cnt = cnt + 1;
                end
                app.logCapData();
                fprintf("[LOG] Feedback loop completed for step %d/%d. Moving to next step.\n\n", i, total_steps);
                app.StepsPhyEditField.Value = scan_loop_step;
            end

            fclose(fid);


        end



        function filteredData = myMedfilt(app, data, windowSize)
            % CUSTOMMEDFILT1 Applies a moving median filter to data.
            %
            % Inputs:
            %   data: 1D array of data points
            %   windowSize: Size of the sliding window (should be odd for best results)
            
            dataLength = length(data);
            filteredData = zeros(size(data));
            
            % Ensure windowSize is odd
            if mod(windowSize, 2) == 0
                windowSize = windowSize + 1;
            end
            
            padSize = floor(windowSize / 2);
            
            % Pad the data to handle edges (replicating edge values)
            paddedData = [repmat(data(1), padSize, 1); data; repmat(data(end), padSize, 1)];
            
            % Apply sliding window
            for i = 1:dataLength
                % Extract the current window
                window = paddedData(i : i + windowSize - 1);
                
                % Calculate the median of the window
                filteredData(i) = median(window);
            end
        end
        
        function zeroingCap(app, volt, microstep, speed)
            % zeroingCap: Performs a scan to find the inflection point of capacitance
            % and moves the apparatus to that position.
            app.StartStopCapButton.Enable = false; % prevent data racing

            searchStep_srt = 10; 
            searchStep_end = -searchStep_srt;
            app.StepsCapZero = -inf;
            app.StepsCapEditField.Value = -inf;

            % --- Move to starting position ---
            app.goTo(volt, searchStep_srt, microstep, speed);
            fprintf('[LOG] Starting scan to find inflection point...\n');

            start_step = app.StepsPhyEditField.Value;
            total_steps = round(abs(searchStep_end - start_step) * (microstep));
            if total_steps == 0, return; end
            step_inc = (1 / microstep) * sign(searchStep_end - start_step);
            app.scanData = zeros(total_steps, 2);
            wait_time = 1 / (speed * microstep);

            % 1. Data Collection Loop
            for i = 1:total_steps
                curr_step = app.StepsPhyEditField.Value + step_inc;
                phase = curr_step * 120;
                app.setVoltages(volt, phase);

                pause(wait_time); 
                app.tim.stop(); 
                app.readLockin();
                if ~isempty(app.capData)
                    app.scanData(i, :) = [curr_step, app.capData(end)];
                    plot(app.AxesCapReadout, app.capData, '-');
                else
                    app.scanData(i, :) = [curr_step, 0]; % Fallback if no data
                end
                if strcmp(app.StartStopCapButton.Text{1}, "Stop")
                    if strcmp(app.tim.Running, 'off'), app.tim.start(); end
                else
                    % fprintf('Timer is already active, skipping start.\n');
                end

                % Update the physical position display in the UI
                app.StepsPhyEditField.Value = curr_step;
                drawnow; 
            end

            % 2. Data Processing (Post-Scan)
            fprintf('[LOG] Data collection finished, processing data...\n');
            validIdx = app.scanData(:,1) ~= 0; 
            dataToProcess = app.scanData(validIdx, :);

            if isempty(dataToProcess) || size(dataToProcess, 1) < 5
                warning('Not enough data to find inflection point.');
                return;
            end

            winLen = 5;
            cleanedValues = app.myMedfilt(dataToProcess(:,2), winLen);

            % 3. Find Inflection Point Logic
            totalPoints = length(cleanedValues);

            % Sliding Window for Inflection Point
            windowSize = min(10, round(totalPoints / 6)); 
            if windowSize < 4
                windowSize = 4;
            end

            if totalPoints < windowSize * 2
                warning('Not enough data points for inflection detection.');
                peakIdx = round(totalPoints/2);
            else
                minDiff = inf;
                peakIdx = round(totalPoints/2);

                for i = windowSize+1 : totalPoints-windowSize
                    headMean = mean(cleanedValues(i-windowSize : i-1));
                    tailMean = mean(cleanedValues(i : i+windowSize-1));

                    currentDiff = abs(tailMean - headMean);

                    if currentDiff < minDiff
                        minDiff = currentDiff;
                        peakIdx = i;
                    end
                end
            end

            % Extract the specific step position of the found extremum
            inflectionStep = dataToProcess(peakIdx, 1);
            inflectionValue = cleanedValues(peakIdx);

            fprintf('[LOG] Zeroing cap scan complete.\n      Inflection Point found at Step: %.2f, Value: %.2e\n', ...
                inflectionStep, inflectionValue);

            % 4. Move to Final Position
            app.StepsCapZero = inflectionStep;
            app.StepsCapEditField.Value = 0;

            % Move device to the identified inflection point
            app.goTo(volt, inflectionStep, microstep, speed);
            app.StartStopCapButton.Enable = true; % prevent data racing
        end
        
        % % Timer Update: stable version (tool function for app.zeroingCap())
        % function finishScan(app)
        %     fprintf('[LOG] Processing scan data...\n');
        %     validIdx = app.scanData(:,1) ~= 0;
        %     dataToProcess = app.scanData(validIdx,:);
        % 
        %     if size(dataToProcess,1) < 5
        %         warning('Not enough data.');
        %         app.isScanning = false;
        %         return;
        %     end
        % 
        %     winLen = 5;
        %     cleanedValues = app.myMedfilt(dataToProcess(:,2), winLen);
        % 
        %     totalPoints = length(cleanedValues);
        %     windowSize = min(10, round(totalPoints/6));
        %     windowSize = max(windowSize,4);
        % 
        %     minDiff = inf;
        %     peakIdx = round(totalPoints/2);
        % 
        %     for i = windowSize+1 : totalPoints-windowSize
        % 
        %         headMean = mean(cleanedValues(i-windowSize:i-1));
        %         tailMean = mean(cleanedValues(i:i+windowSize-1));
        %         currentDiff = abs(tailMean - headMean);
        % 
        %         if currentDiff < minDiff
        %             minDiff = currentDiff;
        %             peakIdx = i;
        %         end
        %     end
        % 
        %     inflectionStep = dataToProcess(peakIdx,1);
        %     inflectionValue = cleanedValues(peakIdx);
        % 
        %     fprintf('[LOG] Inflection at %.3f, value %.3e\n', ...
        %         inflectionStep, inflectionValue);
        % 
        %     app.StepsCapZero = inflectionStep;
        %     app.StepsCapEditField.Value = 0;
        % 
        %     app.goTo(app.zeroVolt, inflectionStep, ...
        %         app.zeroMicrostep, app.zeroSpeed);
        % 
        %     app.isScanning = false;
        % 
        % end
        % function scanStep(app)
        %     if ~app.isScanning
        %         stop(app.scanTimer);
        %         delete(app.scanTimer);
        %         return;
        %     end
        % 
        %     app.scanIndex = app.scanIndex + 1;
        % 
        %     if app.scanIndex > app.totalSteps
        %         stop(app.scanTimer);
        %         delete(app.scanTimer);
        %         app.finishScan();
        %         return;
        %     end
        % 
        %     curr_step = app.StepsPhyEditField.Value + app.stepInc;
        %     phase = curr_step * 120;
        % 
        %     app.setVoltages(app.zeroVolt, phase);
        % 
        %     v = app.readLockin();
        % 
        %     app.scanData(app.scanIndex,:) = [curr_step, v];
        % 
        %     app.StepsPhyEditField.Value = curr_step;
        % 
        % end
        % % Timer Update: stable version
        % function zeroingCap(app, volt, microstep, speed)
        % 
        %     if app.isScanning
        %         return;
        %     end
        % 
        %     app.isScanning = true;
        % 
        %     searchStep_srt = 10;
        %     searchStep_end = -searchStep_srt;
        % 
        %     app.goTo(volt, searchStep_srt, microstep, speed);
        % 
        %     start_step = app.StepsPhyEditField.Value;
        %     app.totalSteps = round(abs(searchStep_end - start_step) * microstep);
        % 
        %     if app.totalSteps == 0
        %         app.isScanning = false;
        %         return;
        %     end
        % 
        %     app.stepInc = (1 / microstep) * sign(searchStep_end - start_step);
        % 
        %     app.scanData = zeros(app.totalSteps, 2);
        %     app.scanIndex = 0;
        % 
        %     app.zeroVolt = volt;
        %     app.zeroMicrostep = microstep;
        %     app.zeroSpeed = speed;
        % 
        %     wait_time = 1 / (speed * microstep);
        % 
        %     app.scanTimer = timer( ...
        %         'ExecutionMode','fixedRate', ...
        %         'Period', wait_time, ...
        %         'BusyMode','drop', ...
        %         'TimerFcn', @(~,~) app.scanStep());
        %     start(app.scanTimer);
        % 
        % end
        % 




    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            image(app.UIAxes, zeros(1280,960));
            axis(app.UIAxes, 'tight');
            colormap(app.UIAxes,'gray');
            app.C=[];
            app.hImage=[];
            app.hImage2=[];
            app.snapshot=[];
            app.vw=[];
            app.scanUSBPort();
            app.LockIn=[];
            app.capData=[];
            app.capFitData=[];
            app.setSteppingEnable(false);
            app.setPositioningEnable(false);
            app.updateAvailableDevices();
%             app.DeviceName.Value = app.DeviceName.ItemsData{end};
            app.updateAvailableResolution();

            app.volt_amp_ratio = 25;
            % app.tim = timer('Period',0.25,...
            %     'ExecutionMode', 'fixedSpacing', ...
            %     'TasksToExecute', Inf);
            % app.tim.TimerFcn = @app.timerCallback;

            % Timer Update: stable version
            app.tim = timer( ...
                'Period', 0.25, ...
                'ExecutionMode', 'fixedRate', ...
                'BusyMode', 'drop', ...
                'TasksToExecute', Inf);
            app.tim.TimerFcn = @(~,~) app.timerCallback();

        end

        % Button pushed function: StartPreviewButton
        function StartPreviewButtonPushed(app, event)
            if strcmpi(app.StartPreviewButton.Text, "Start Preview")
                app.cameraConfig(app.ExposureSlider.Value, app.GainSlider.Value);
                app.cameraStartPreview();
                app.StartPreviewButton.Text = "Stop Preview";
                app.PickcircumferenceButton.Enable=true;
                app.SnapshotButton.Enable=true;
                app.ContinuousButton.Enable=true;
                app.RecordMovieButton.Enable=true;
                app.ChoosethetarangeButton.Enable=true;
                app.ChooseradialrangeButton.Enable=true;
            else
                app.cameraStopPreview();
                app.StartPreviewButton.Text = "Start Preview";
                app.PickcircumferenceButton.Enable=false;
                app.SnapshotButton.Enable=false;
                app.RecordMovieButton.Enable=false;
                app.ContinuousButton.Enable=false;
                app.ChoosethetarangeButton.Enable=false;
                app.ChooseradialrangeButton.Enable=false;
            end
        end

        % Button pushed function: PickcircumferenceButton
        function PickcircumferenceButtonPushed(app, event)
            if length(app.res)~=2
                errordlg('Start video first');
                return;   
            end
            
            for i = length(app.UIAxes.Children):-1:1
                if ~isa(app.UIAxes.Children(i), 'matlab.graphics.primitive.Image')
                    delete(app.UIAxes.Children(i));
                end
            end
            [x,y]=ginputuiax(app.UIAxes, 3);
            if length(x)~=3
                return;
            end
            
            x=round(x);
            y=round(y);
            
            % Get snapshot and find maximum in gradient
            F=getsnapshot(app.video);
            if size(F,3)~=1
                F=rgb2gray(F);
            end
            G = imgradient(double(F));
            w=10;
            for i=1:3
                Gcrop = G(max(1,y(i)-w):min(end,y(i)+w),max(1,x(i)-w):min(end,x(i)+w)); % Window centered at (x,y)
                [~,ind]=max(Gcrop(:));
                [yy,xx] = ind2sub(size(Gcrop), ind);
                y(i)=yy+max(1,y(i)-w)-1;
                x(i)=xx+max(1,x(i)-w)-1;
            end
            
            [r, app.C] = fit_circle_through_3_points([x(:),y(:)]);
            if any(isnan(app.C))
                return;
            end
            rectangle(app.UIAxes, 'Position', [app.C(1)-r, app.C(2)-r, 2*r, 2*r], 'Curvature', [1,1],'EdgeColor','r', 'LineWidth', 2);
            
            % Calculate necessary matrices to accelerate processing
            app.prepare();
            
            app.rect_r=[];
            app.rect_th=[];
        end

        % Value changing function: ExposureSlider
        function ExposureSliderValueChanging(app, event)
            value = event.Value;
            if isa(app.video, 'videoinput') && ~isrunning(app.video)
                src = getselectedsource(app.video);
                try
                    src.Exposure = value;
                catch
                end
            end
        end

        % Value changing function: GainSlider
        function GainSliderValueChanging(app, event)
            value = event.Value;
            if isa(app.video, 'videoinput') && ~isrunning(app.video)
                src = getselectedsource(app.video);
                try
                    src.Gain = value;
                catch
                end
            end
        end

        % Button pushed function: SnapshotButton
        function SnapshotButtonPushed(app, event)
            if ~isa(app.video, 'videoinput')
                errordlg('Start preview first');
                return;
            end
            
            app.snapshot = getsnapshot(app.video);
            
            if length(app.C)==2
                % Process one frame
                app.process(app.snapshot);
                return;
            else
                % Directly plot
                app.hImage2=image(app.UIAxes2, app.snapshot);
                axis(app.UIAxes2, 'tight');
            end
            
            if app.SyncSnapshotinworkspaceCheckBox.Value
                % Save snapshot to workspace
                assignin('base', app.VariableNameEditField.Value, app.snapshot);
            end
            
        end

        % Button pushed function: ContinuousButton
        function ContinuousButtonPushed(app, event)
            if ~isa(app.video, 'videoinput')
                errordlg('Start preview first');
                return;
            end
            
            if length(app.C)~=2
                errordlg('Pick circumference first');
                return;
            end
            
            if strcmpi(app.ContinuousButton.Text, 'Continuous')
                % Setup continuous acquisition+analysis
                app.RecordMovieButton.Enable = false;
                app.SnapshotButton.Enable = false;
                triggerconfig(app.video, 'immediate');
                
                % Set callback
                app.video.TriggerFcn = {@triggerfn, app};
                start(app.video);
                
                app.ContinuousButton.Text = {'Stop','Continuous'};
            else
                app.ContinuousButton.Text = 'Continuous';
                stop(app.video);
                app.video.TriggerFcn = {};
                app.RecordMovieButton.Enable = true;
                app.SnapshotButton.Enable = true;
            end
            
            function triggerfn(vid, event, app)
                F=getdata(vid,1);
                app.process(F);
                flushdata(vid,'all'); % Remove all frames stored in memory
                if ~isempty(app.hImage)
                    app.hImage.CData = F;
                end
                drawnow;
            end
        end

        % Button pushed function: ChoosethetarangeButton
        function ChoosethetarangeButtonPushed(app, event)
            if isempty(app.snapshot)
                errordlg('Take snapshot first');
                return;
            end
            for i = length(app.UIAxes2.Children):-1:1
                if isa(app.UIAxes2.Children(i), 'images.roi.Rectangle') && strcmp(app.UIAxes2.Children(i).Label, 'Theta')
                    delete(app.UIAxes2.Children(i));
                end
            end
%             [~,y]=ginputuiax(app.UIAxes2,2,'Y');
%             if length(y)~=2
%                 return;
%             end
%             line(app.UIAxes2,[app.THETA(1,1),app.THETA(end,1)],[1,1]*y(1),'Color','red');
%             line(app.UIAxes2,[app.THETA(1,1),app.THETA(end,1)],[1,1]*y(2),'Color','red');
%             app.rsel = sort(y);

            app.rect_th = drawrectangle(app.UIAxes2, 'Label','Theta', 'LineWidth', 1);
            addlistener(app.rect_th, 'DeletingROI', @deletecb);
            app.process(app.snapshot);
            
            % Called before deletion to make sure we don't use it again
            function deletecb(src,evt)
                a=ancestor(src, 'Figure');
                ap=a.RunningAppInstance;
                ap.rect_th=[];
            end
        end

        % Button pushed function: ChooseradialrangeButton
        function ChooseradialrangeButtonPushed(app, event)
            if isempty(app.snapshot)
                errordlg('Take snapshot first');
                return;
            end
            for i = length(app.UIAxes2.Children):-1:1
                if isa(app.UIAxes2.Children(i), 'images.roi.Rectangle') && strcmp(app.UIAxes2.Children(i).Label, 'Radial')
                    delete(app.UIAxes2.Children(i));
                end
            end
            app.rect_r = drawrectangle(app.UIAxes2, 'Label','Radial', 'LineWidth', 1);
            addlistener(app.rect_r, 'DeletingROI', @deletecb);
            app.process(app.snapshot);
            
            % Called before deletion to make sure we don't use it again
            function deletecb(src,evt)
                a=ancestor(src, 'Figure');
                ap=a.RunningAppInstance;
                ap.rect_r=[];
            end
        end

        % Button pushed function: ClearGraphsButton
        function ClearGraphsButtonPushed(app, event)
            app.disp_r=[];
            app.disp_th=[];
            app.phase_r=[];
            app.phase_th=[];
            cla(app.UIAxes5);
            cla(app.UIAxes4);
        end

        % Value changed function: AlignpeakCheckBox
        function AlignpeakCheckBoxValueChanged(app, event)
            value = app.AlignpeakCheckBox.Value;
            app.disp_th=[];
            app.phase_th=[];
            cla(app.UIAxes5);
        end

        % Drop down opening function: SR830DropDown
        function SR830DropDownOpening(app, event)
            app.scanUSBPort();
        end

        % Button pushed function: SR830ConnectButton
        function SR830ConnectButtonPushed(app, event)
            if isa(app.LockIn, 'visa')
                fclose(app.LockIn);
                delete(app.LockIn);
            end

            prev = instrfind(Type='visa-gpib');
            if ~isempty(prev) && ~app.prev_closed
                for p=prev
                    fclose(p);
                end
                app.prev_closed = true;
            end
            try
                app.LockIn=visa('NI', app.SR830DropDown.Value, 'EOSMode', 'read&write', 'EOSCharCode', 'CR', 'Timeout', 0.5);
                fopen(app.LockIn);
                app.updateVoltages();
                if isa(app.LockIn, 'visa') && strcmp(app.LockIn.Status, 'open')
                    % Log detailed success info
                    fprintf('[LOG] SR830 successfully connected.\n');
                    fprintf('      Resource: %s\n', app.SR830DropDown.Value);
                    fprintf('      Status:   OPEN\n');
                    % fprintf('      Purpose:  High-frequency Readout via OUTP? 1\n');
                    
                    % Optional: Query identification string to verify instrument model
                    idn = query(app.LockIn, '*IDN?');
                    fprintf('      Instrument ID: %s\n', strtrim(idn));
                else
                    % This part might not be reached if fopen throws an error, but kept for safety
                    error('VISA status is not open.');
                end
                app.V1Knob.Value = app.V1Gauge.Value;
                app.V2Knob.Value = app.V2Gauge.Value;
                app.V3Knob.Value = app.V3Gauge.Value;
                app.V1Label.Value = app.V1Gauge.Value;
                app.V2Label.Value = app.V2Gauge.Value;
                app.V3Label.Value = app.V3Gauge.Value;

                app.SetButton.Enable = true;
                app.InitButton.Enable = true;
                % app.StartStopCapButton.Enable = true;
                app.SetVdButton.Enable = true;
                app.SetVzButton.Enable = true;
                app.ZeroVzButton.Enable = true;
                app.FlashVzButton.Enable = true;

                app.VdField.Value = str2double(query(app.LockIn, "auxv? 4\n"));
                app.VzField.Value = app.VzGauge.Value;
            catch e
                % Detailed error logging
                fprintf('[ERROR] Failed to connect to SR830 at %s.\n', app.SR830DropDown.Value);
                fprintf('        Reason: %s\n', e.message);
                errordlg(sprintf('SR830 connection not established: %s', e.message));
                return;
            end
        end

        % Button pushed function: OUTPUTOFFButton
        function OUTPUTOFFButtonPushed(app, event)
            app.V1Knob.Value=0;
            app.V2Knob.Value=0;
            app.V3Knob.Value=0;
            app.VzField.Value=0;
            app.V1Label.Value = 0;
            app.V2Label.Value = 0;
            app.V3Label.Value = 0;
            app.setSteppingEnable(false); % Stepping needs to be reinitialized
            app.setPositioningEnable(false); % Positioning needs to be reinitialized
            if isa(app.LockIn, 'visa') && strcmp(app.LockIn.Status,'open')
                fprintf(app.LockIn, 'auxv 1,0\n');
                fprintf(app.LockIn, 'auxv 2,0\n');
                fprintf(app.LockIn, 'auxv 3,0\n');
                fprintf(app.LockIn, 'auxv 4,0\n');
                pause(0.1);
                app.updateVoltages();
            end
        end

        % Value changing function: V1Knob
        function V1KnobValueChanging(app, event)
            app.V1Label.Value = event.Value;
        end

        % Value changing function: V2Knob
        function V2KnobValueChanging(app, event)
            app.V2Label.Value = event.Value;
        end

        % Value changing function: V3Knob
        function V3KnobValueChanging(app, event)
            app.V3Label.Value = event.Value;
        end

        % Button pushed function: SetButton
        function SetButtonPushed(app, event)
            if isa(app.LockIn, 'visa') && strcmp(app.LockIn.Status,'open')
                fprintf(app.LockIn, sprintf('auxv 1,%f\n', app.V1Knob.Value/app.volt_amp_ratio));
                fprintf(app.LockIn, sprintf('auxv 2,%f\n', app.V2Knob.Value/app.volt_amp_ratio));
                fprintf(app.LockIn, sprintf('auxv 3,%f\n', app.V3Knob.Value/app.volt_amp_ratio));
                app.updateVoltages();
            end
            
            app.setSteppingEnable(false); % Stepping needs to be reinitialized
            app.setPositioningEnable(false); % Positioning needs to be reinitialized
        end

        % Value changed function: V1Label
        function V1LabelValueChanged(app, event)
            app.V1Knob.Value = app.V1Label.Value;
        end

        % Value changed function: V2Label
        function V2LabelValueChanged(app, event)
            app.V2Knob.Value = app.V2Label.Value;
        end

        % Value changed function: V3Label
        function V3LabelValueChanged(app, event)
            app.V3Knob.Value = app.V3Label.Value;
        end

        % Button pushed function: InitButton
        function InitButtonPushed(app, event)
            app.PhaseEditField.Value = 0;
            app.StepsEditField.Value = 0;
            app.setVoltages(app.DriveVoltageEditField.Value, app.PhaseEditField.Value);
            app.setSteppingEnable(true);
            if app.SyncCapFittingCheckBox.Value
                app.addCapFittingPoint();
            end
        end

        % Button pushed function: Minus18Button
        function Minus18ButtonPushed(app, event)
            app.PhaseEditField.Value = mod(app.PhaseEditField.Value - 15, 360);
            app.StepsEditField.Value = app.StepsEditField.Value - 0.125;
            app.setVoltages(app.DriveVoltageEditField.Value, app.PhaseEditField.Value);
            if app.SyncCapFittingCheckBox.Value
                app.addCapFittingPoint();
            end
        end

        % Button pushed function: Plus18Button
        function Plus18ButtonPushed(app, event)
            app.PhaseEditField.Value = mod(app.PhaseEditField.Value + 15, 360);
            app.StepsEditField.Value = app.StepsEditField.Value + 0.125;
            app.setVoltages(app.DriveVoltageEditField.Value, app.PhaseEditField.Value);
            if app.SyncCapFittingCheckBox.Value
                app.addCapFittingPoint();
            end
        end

        % Button pushed function: MinusHalfButton
        function MinusHalfButtonPushed(app, event)
            app.PhaseEditField.Value = mod(app.PhaseEditField.Value - 60, 360);
            app.StepsEditField.Value = app.StepsEditField.Value - 0.5;
            app.setVoltages(app.DriveVoltageEditField.Value, app.PhaseEditField.Value);
            if app.SyncCapFittingCheckBox.Value
                app.addCapFittingPoint();
            end
        end

        % Button pushed function: PlusHalfButton
        function PlusHalfButtonPushed(app, event)
            app.PhaseEditField.Value = mod(app.PhaseEditField.Value + 60, 360);
            app.StepsEditField.Value = app.StepsEditField.Value + 0.5;
            app.setVoltages(app.DriveVoltageEditField.Value, app.PhaseEditField.Value);
            if app.SyncCapFittingCheckBox.Value
                app.addCapFittingPoint();
            end
        end

        % Button pushed function: MinusQuarterButton
        function MinusQuarterButtonPushed(app, event)
            app.PhaseEditField.Value = mod(app.PhaseEditField.Value - 30, 360);
            app.StepsEditField.Value = app.StepsEditField.Value - 0.25;
            app.setVoltages(app.DriveVoltageEditField.Value, app.PhaseEditField.Value);
            if app.SyncCapFittingCheckBox.Value
                app.addCapFittingPoint();
            end
        end

        % Button pushed function: PlusQuarterButton
        function PlusQuarterButtonPushed(app, event)
            app.PhaseEditField.Value = mod(app.PhaseEditField.Value + 30, 360);
            app.StepsEditField.Value = app.StepsEditField.Value + 0.25;
            app.setVoltages(app.DriveVoltageEditField.Value, app.PhaseEditField.Value);
            if app.SyncCapFittingCheckBox.Value
                app.addCapFittingPoint();
            end
        end

        % Button pushed function: ZeroStepsButton
        function ZeroStepsButtonPushed(app, event)
            app.StepsEditField.Value = 0;
        end

        % Drop down opening function: Resolution
        function ResolutionDropDownOpening(app, event)
            app.updateAvailableResolution();
        end

        % Value changed function: DeviceName
        function DeviceNameValueChanged(app, event)
            app.updateAvailableResolution();
            cla(app.UIAxes);
            app.C=[];
            app.hImage=[];
        end

        % Value changed function: Resolution
        function ResolutionValueChanged(app, event)
            cla(app.UIAxes);
            app.C=[];
            app.hImage=[];
        end

        % Drop down opening function: DeviceName
        function DeviceNameDropDownOpening(app, event)
            app.updateAvailableDevices();
        end

        % Button pushed function: SavesnapshottofileButton
        function SavesnapshottofileButtonPushed(app, event)
            if isnumeric(app.snapshot) && ~isempty(app.snapshot)
                [file, path] = uiputfile('*.png', 'Save snapshot as', 'Fsnap.png');
                if file ~= 0
                    imwrite(app.snapshot, fullfile(path, file));
                end
            end
        end

        % Button pushed function: PickButton
        function PickButtonPushed(app, event)
            [file,path] = uigetfile('*.png', 'Pick reference image');
            app.ReferenceDesignEditField.Value = [path file];
        end

        % Button pushed function: PreviewRefDesignButton
        function PreviewRefDesignButtonPushed(app, event)
            if isfile(app.ReferenceDesignEditField.Value)
                G = imread(app.ReferenceDesignEditField.Value);
                G = imbinarize(G);
                G = imcomplement(G);
                h=figure;
                imshow(G);
            end
        end

        % Button pushed function: 
        % AlignsnapshottoRefDesignExperimentalButton
        function AlignsnapshottoRefDesignExperimentalButtonPushed(app, event)
            if isfile(app.ReferenceDesignEditField.Value) && isnumeric(app.snapshot) && ~isempty(app.snapshot)
                G = imread(app.ReferenceDesignEditField.Value);
                app.alignReference(app.snapshot, G, str2double(app.MicroscopeObjectiveDropDown.Value), app.MinScaleEditField.Value);
            end
        end

        % Button pushed function: InitResetButton
        function InitResetButtonPushed(app, event)
            if ~isa(app.LockIn, 'visa') || ~strcmp(app.LockIn.Status,'open')
                errordlg('Initializing failed: USB connection not established.');
                return;
            end
            app.StepsPhyEditField.Value = 0;
            app.StepsCapZero = -inf;
            app.StepsCapEditField.Value = -inf;
            app.TargetStepEditField.Value = 0;
            app.capFitData = []; 
            % app.k_fit = 0; % Don't Clean fit data
            % app.b_fit = 0; % Don't Clean fit data

            app.setVoltages(0, 0); % Turn off voltage to reset the motor
            pause(0.5);
            app.setVoltages(app.DriveVoltageEditField_2.Value, 0);

            app.setPositioningEnable(true);
            app.LogCapDataButton.Enable = true;
            app.CloseLoopSwitch.Value = "Off";
            app.CloseLoopSwitch.Enable = true;
            app.ZeroCapButton.Enable = true;
            app.ReadChannelDropDown.Editable = "on";
        end

        % Button pushed function: GOButton
        function GOButtonPushed(app, event)
            % "for loop" params
            speed = app.SpeedstepsEditField.Value;
            if speed <= 0
                errordlg(sprintf('Speed must larger than 0! [speed=%.3f]', speed));
                return;
            end
            microstep = str2double(app.MicroStepDropDown.Value);
            target = round(app.TargetStepEditField.Value*microstep)/microstep; % Round to nearest microstep
            
            % prevent channel switching while moving MEMS
            app.ReadChannelDropDown.Editable = "off";
            app.ZeroCapButton.Enable = false; % prevent data racing
            % swith mode
            if strcmp(app.CloseLoopSwitch.Value, "On")
                % check if model is calibrated
                if isempty(app.k_fit) || app.k_fit == 0
                    uialert(ancestor(app.SCANButton, 'figure'),'Model calibration failed. Use Positioning(Feedback)Tab to build scan model first.', 'Error');
                return;
                else
                    fprintf("Calibration Successful: Slope (k) = %.3e, Intercept (b) = %.3e\n", app.k_fit, app.b_fit);
                end
                
                app.goTo_fb_scan(app.DriveVoltageEditField_2.Value, target, microstep, speed)
            else
                app.goTo(app.DriveVoltageEditField_2.Value, target, microstep, speed);
            end
            app.ReadChannelDropDown.Editable = "on";
            app.ZeroCapButton.Enable = true; % prevent data racing
        end

        % Button pushed function: RecordMovieButton
        function RecordMovieButtonPushed(app, event)

             if strcmpi(app.RecordMovieButton.Text, "Record Movie")
                
                % Save video
                [file, path] = uiputfile({'*.mp4';'*.avi';'*.*'}, 'Save movie as', sprintf('movie-%s.mp4', datestr(datetime('now'), "yymmdd-HHMMSS")));
                if file==0
                    return;
                end
                filename = fullfile(path, file);
                [~,~,ext] = fileparts(filename);
                if strcmpi(ext, '.avi')
                    app.vw = VideoWriter(filename, "Uncompressed AVI");
                else
                    app.vw = VideoWriter(filename, "MPEG-4");
                    app.vw.Quality = 85;
                end
                app.vw.FrameRate = 30;
                open(app.vw);
                
                triggerconfig(app.video, 'manual');
                % Set callback
                app.video.TriggerFcn = {@triggerfn, app};
                start(app.video);
                    
                app.RecordMovieButton.Text = "Stop Record";
                app.SnapshotButton.Enable=false;
                app.ContinuousButton.Enable=false;
             else
                
                stop(app.video);
                app.video.TriggerFcn = {};
                
                close(app.vw);
                
                app.RecordMovieButton.Text = "Record Movie";
                app.SnapshotButton.Enable=true;
                app.ContinuousButton.Enable=true;
             end
             
             
            function triggerfn(vid, event, app)
                F=getdata(vid,1);
                writeVideo(app.vw, F);
                flushdata(vid,'all'); % Remove all frames stored in memory
                if ~isempty(app.hImage)
                    app.hImage.CData = F;
                    drawnow;
                end
            end
        end

        % Button pushed function: SetVdButton
        function SetVdButtonPushed(app, event)
            if isa(app.LockIn, 'visa') && strcmp(app.LockIn.Status,'open')
                fprintf(app.LockIn, "auxv 4,%f\n", app.VdField.Value);
            end
        end

        % Button pushed function: StartStopCapButton
        function StartStopCapButtonPushed(app, event)
            if strcmp(app.StartStopCapButton.Text{1}, "Start")
                % Start update
                app.StartStopCapButton.Text{1} = 'Stop';
                app.StartStopCapButton.BackgroundColor = 'red';
                app.StartStopCapButton.FontColor = [1, 1, 1];
                app.ClearCapButton.Enable = false;
                app.tim.start();
            else
                % % Timer Update: to stop app.zeroingCap()
                % app.isScanning = false;
                % if ~isempty(app.scanTimer) && isvalid(app.scanTimer)
                %     stop(app.scanTimer);
                %     delete(app.scanTimer);
                %     fprintf("[LOG] app.scanTimer is stopped.\n")
                % end

                app.StartStopCapButton.Text{1} = 'Start';
                app.StartStopCapButton.BackgroundColor = 'black';
                app.StartStopCapButton.FontColor = [1, 1, 1];
                app.ClearCapButton.Enable = true;
                app.tim.stop();
            end
        end

        % Button pushed function: ClearCapButton
        function ClearCapButtonPushed(app, event)
            app.capData = [];
            app.capFitData = [];
            cla(app.AxesCapReadout, 'reset');
            title(app.AxesCapReadout, {'Capacitance Readout'; ''})
            xlabel(app.AxesCapReadout, 'Time')
            ylabel(app.AxesCapReadout, 'Signal')
            % cla(app.AxesCapFitting);
        end

        % Button pushed function: SetVzButton
        function SetVzButtonPushed(app, event)
            if isa(app.LockIn, 'visa') && strcmp(app.LockIn.Status,'open')
                app.setVz(app.VzField.Value);
                app.updateVoltages();
            end
        end

        % Button pushed function: FlashVzButton
        function FlashVzButtonPushed(app, event)
            if isa(app.LockIn, 'visa') && strcmp(app.LockIn.Status,'open')
                for i=1:20
                    app.setVz(app.VzField.Value);
                    app.updateVoltages();
                    pause(0.1);
                    app.setVz(0);
                    app.updateVoltages();
                    pause(0.1);
                end
                app.setVz(0);
                app.updateVoltages();
            end
        end

        % Button pushed function: ZeroVzButton
        function ZeroVzButtonPushed(app, event)
            if isa(app.LockIn, 'visa') && strcmp(app.LockIn.Status,'open')
                app.setVz(0);
                app.updateVoltages();
            end
        end

        % Button pushed function: SR844ConnectButton
        function SR844ConnectButtonPushed(app, event)
            if isa(app.LockIn_Read, 'visa')
                fclose(app.LockIn_Read);
                delete(app.LockIn_Read);
            end
            prev = instrfind(Type='visa-gpib');
            if ~isempty(prev) && ~app.prev_closed
                for p=prev
                    fclose(p);
                end
                app.prev_closed = true;
            end
            try
                app.LockIn_Read=visa('NI', app.SR844DropDown.Value, 'EOSMode', 'read&write', 'EOSCharCode', 'LF', 'Timeout', 0.5);
                fopen(app.LockIn_Read);
                if isa(app.LockIn_Read, 'visa') && strcmp(app.LockIn_Read.Status, 'open')
                    % Log detailed success info
                    fprintf('[LOG] SR844 successfully connected.\n');
                    fprintf('      Resource: %s\n', app.SR844DropDown.Value);
                    fprintf('      Status:   OPEN\n');
                    fprintf('      Purpose:  High-frequency Readout via OUTP? 1\n');
                    
                    % Optional: Query identification string to verify instrument model
                    idn = query(app.LockIn_Read, '*IDN?');
                    fprintf('      Instrument ID: %s\n', strtrim(idn));

                    app.StartStopCapButton.Enable = true;
                else
                    % This part might not be reached if fopen throws an error, but kept for safety
                    error('VISA status is not open.');
                end

            catch e
                % Detailed error logging
                fprintf('[ERROR] Failed to connect to SR844 at %s.\n', app.SR844DropDown.Value);
                fprintf('        Reason: %s\n', e.message);
                errordlg(sprintf('SR844 connection not established: %s', e.message));
                return;
            end
        end

        % Button pushed function: InitResetButton_2
        function InitResetButton_2Pushed(app, event)
            % 1. Reset UI Fields
            app.CurrentStepEditField.Value = 0;
            app.ScanStepSrt.Value = 0;
            app.ScanStepEnd.Value = 7;
            app.FinalStepEditField.Value = 0;
            app.FilterWinEditField.Value = 1;
            % app.StepMinLossEditField.Value = 0.05; % close "I" when StepLoss < 5%
            % app.SpeedstepsEditField_2.Value = 0.2;
            
            % 2. Reset Data Properties
            % Clear the array so old points don't appear in the next 'goTo' call
            app.capFitData = []; 
            app.k_fit = 0; % Clean fit data
            app.b_fit = 0;
        
            % 3. Clear the Plot
            % 'reset' ensures legends and hold states are also cleared
            app.plotCapFBloss = false;
            cla(app.AxesCapReadout, 'reset');
            title(app.AxesCapReadout, {'Capacitance Readout'; ''})
            xlabel(app.AxesCapReadout, 'Time')
            ylabel(app.AxesCapReadout, 'Signal')

            cla(app.AxesCapFitting, 'reset');
            title(app.AxesCapFitting, 'Fitting');
            xlabel(app.AxesCapFitting, 'Steps');
            ylabel(app.AxesCapFitting, 'Signal');
            
            % 4. Reset Hardware State
            app.setVoltages(0, 0); % Turn off voltage to reset the motor
    
            pause(0.5);
            app.setVoltages(app.DriveVoltageEditField_3.Value, 0);
            
            % 5. Enable Controls
            app.positioningSTOP = false;
            app.CANCELButton.BackgroundColor = 'red';
            app.CANCELButton.Enable = false;
            app.LogPlotSwitch.Enable = true;
            app.setPositioningEnable_2(true);
            app.setPositioningEnable_3(false);
        end

        % Button pushed function: GOButton_2
        function GOButton_2Pushed(app, event)
            % Retrieve parameters from UI components
            drive_volt = app.DriveVoltageEditField_3.Value;
            final_pos  = app.FinalStepEditField.Value;
            if app.PIDMicrostepDropDown.Value == "inf"
                microstep2 = 512;
            else
                microstep2 = str2double(app.PIDMicrostepDropDown.Value);
            end
            if isempty(app.StepMinLossEditField.Value)
                app.StepMinLossEditField.Value = 0.0001; % Default step_tolerance
            end
            
            % Disable button after GO button pushed
            app.LogPlotSwitch.Enable = false;
            app.GOButton_2.Enable = false;
            app.CANCELButton.Enable = true;

            % Move from the end of the scan range to the user-defined final position
            if app.CurrentStepEditField.Value ~= final_pos
                microstep1 = str2double(app.ScanMicrostepDropDown.Value);
                speed_move = 1; % faster
                app.goTo_fit(drive_volt, app.CurrentStepEditField.Value, final_pos, microstep1, speed_move, false);
            end

            if app.positioningSTOP
                app.positioningSTOP = false;
                app.CANCELButton.BackgroundColor = 'red';
                app.goTo_fb(drive_volt, final_pos, microstep2);
            else
                app.goTo_fb(drive_volt, final_pos, microstep2);
            end
        end

        % Button pushed function: CANCELButton
        function CANCELButtonPushed(app, event)
            app.positioningSTOP = true;
            app.CANCELButton.BackgroundColor = '#C0C0C0';
            app.GOButton_2.Enable = true;
            app.CANCELButton.Enable = false;
        end

        % Button pushed function: SCANButton
        function SCANButtonPushed(app, event)
            % SCANBUTTONPUSHED Executes a two-stage motion sequence:
            % Phase 1: Model Calibration (Scanning to find Sensitivity k and Offset b)
            % Phase 2: Open-loop Targeted Positioning
            
            % --- 1. Parameter Retrieval ---
            % Extract drive settings from the User Interface
            microstep1 = str2double(app.ScanMicrostepDropDown.Value);
            speed = app.SpeedstepsEditField_2.Value;
            speed_move = 1; % faster
            drive_volt = app.DriveVoltageEditField_3.Value;
            
            % Define the trajectory coordinates for the motion sequence
            origin_pos   = 0;                               % Hardware zero reference
            scan_pos_srt = app.ScanStepSrt.Value;           % Start boundary for linear fitting
            scan_pos_end = app.ScanStepEnd.Value;           % End boundary for linear fitting
            final_pos    = app.FinalStepEditField.Value;    % Ultimate target position
            
            % --- 2. Calibration & Model Building ---

            % (Step A) Fast move from origin to the start of the scan range
            if origin_pos ~= scan_pos_srt
                app.goTo_fit(drive_volt, origin_pos, scan_pos_srt, microstep1, speed_move, false);
            end

            % (Step B Precision scan across the designated range
            % update_k_b = true: Perform real-time Polyfit to calculate k_fit and b_fit
            app.goTo_fit(drive_volt, scan_pos_srt, scan_pos_end, microstep1, speed, true);
            idx_end_part1 = size(app.capFitData, 1);
            app.goTo_fit(drive_volt, scan_pos_end, scan_pos_srt, microstep1, speed, true);
            
            % Split the data
            data_part1 = app.capFitData(1:idx_end_part1, :);
            data_part2 = app.capFitData(idx_end_part1+1:end, :);
            % Fits
            p1 = polyfit(data_part1(:,1), data_part1(:,2), 1);
            p2 = polyfit(data_part2(:,1), data_part2(:,2), 1);
            p_all = polyfit(app.capFitData(:,1), app.capFitData(:,2), 1);
            app.k_fit = p_all(1); app.b_fit = p_all(2);
            % (Step C) Point-by-Point Hysteresis Analysis
            data_part2_flipped = flipud(data_part2);
            n_points = min(size(data_part1, 1), size(data_part2_flipped, 1));
            dy_list = zeros(n_points, 1);
            for i = 1:n_points
                % Direct comparison: Y_backward - Y_forward at the exact same X
                y_fwd = data_part1(i, 2);
                y_bwd = data_part2_flipped(i, 2);
                dy_list(i) = abs(y_bwd - y_fwd);
            end
            % Calculate Mean Hysteresis Metrics
            mean_dy = mean(dy_list);
            mean_dx = mean_dy / abs(app.k_fit);

            % Log
            fid = fopen(app.LogFileName, 'w');
            if fid == -1, error('[LOG] Could not create log file.'); end
            logAndPrint = @(fmt, varargin) cellfun(@(target) fprintf(target, fmt, varargin{:}), {1, fid}, 'UniformOutput', false);
            logAndPrint("[LOG] SCAN_PART Channel               : %s\n", app.ReadChannelDropDown.Value);
            logAndPrint("[LOG] SCAN_PART_1 (Forward %d points) : k1 = %.3e,\tb1 = %.3e\n", size(data_part1, 1), p1(1), p1(2));
            logAndPrint("[LOG] SCAN_PART_2 (Backward %d points): k2 = %.3e,\tb2 = %.3e\n", size(data_part2_flipped, 1), p2(1), p2(2));
            logAndPrint("[LOG] SCAN_TOTAL  (Symmetric)         : k_total = %.3e,\tb_total = %.3e\n", app.k_fit, app.b_fit);
            logAndPrint("[LOG] HYSTERESIS_AVG: Mean Delta_Cap  (dy) = %.3e\n", mean_dy);
            logAndPrint("[LOG] HYSTERESIS_AVG: Mean Delta_Step (dx) = %.4f steps\n", mean_dx);
            fclose(fid);

            for i = 1:size(app.capFitData, 1)
                fprintf("\t(%.6f, %.6e)\n", app.capFitData(i, 1), app.capFitData(i, 2));
            end


            % --- 3. Targeted Positioning ---
            % Move from the end of the scan range to the user-defined final position
            if app.CurrentStepEditField.Value ~= final_pos
                % Model is already built; move to target without further updating k/b
                app.goTo_fit(drive_volt, app.CurrentStepEditField.Value, final_pos, microstep1, speed_move, false);
            end
            
            % Enable subsequent feedback controls
            app.setPositioningEnable_3(true);
            app.CANCELButton.Enable = true;
        end

        % Button pushed function: PlotCapLossButton
        function PlotCapLossButtonPushed(app, event)
            app.plotCapFBloss = ~(app.plotCapFBloss);
            if app.plotCapFBloss
                app.PlotCapLossButton.BackgroundColor = [0, 0.5, 0];
                app.PlotCapLossButton.FontColor = [1, 1, 1];
            else
                app.PlotCapLossButton.BackgroundColor = [1, 1, 1];
                app.PlotCapLossButton.FontColor = [0, 0, 0];
            end
        end

        % Button pushed function: LogCapDataButton
        function LogCapDataButtonPushed(app, event)
            app.logCapData();
        end

        % Button pushed function: ZeroCapButton
        function ZeroCapButtonPushed(app, event)
            % "for loop" params
            volt = app.DriveVoltageEditField_2.Value;
            speed = app.SpeedstepsEditField.Value;
            if speed <= 0
                errordlg(sprintf('Speed must larger than 0! [speed=%.3f]', speed));
                return;
            end
            microstep = str2double(app.MicroStepDropDown.Value);
            % disable GOButton
            app.GOButton.Enable = false;
            app.ReadChannelDropDown.Editable = "off";

            % search for the zero capacitance step 
            if ~strcmp(app.ReadChannelDropDown.Value, 'R')
                errordlg('Must use R channel to zero capacitance step.');
                return;
            else
                app.zeroingCap(volt, microstep, speed);
            end

            % enable button
            app.ReadChannelDropDown.Editable = "on";
            app.GOButton.Enable = true;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create MEMSStepperControlUIFigure and hide until all components are created
            app.MEMSStepperControlUIFigure = uifigure('Visible', 'off');
            app.MEMSStepperControlUIFigure.Position = [100 100 2036 779];
            app.MEMSStepperControlUIFigure.Name = 'MEMS Stepper Control';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.MEMSStepperControlUIFigure);
            app.GridLayout.ColumnWidth = {'6x', 171, '5x', '4x'};
            app.GridLayout.RowHeight = {170, '1x', '1x'};

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout);
            title(app.UIAxes, 'Camera View')
            app.UIAxes.XTick = [];
            app.UIAxes.YTick = [];
            app.UIAxes.Box = 'on';
            app.UIAxes.Layout.Row = [2 3];
            app.UIAxes.Layout.Column = [1 2];

            % Create TabGroup2
            app.TabGroup2 = uitabgroup(app.GridLayout);
            app.TabGroup2.Layout.Row = 1;
            app.TabGroup2.Layout.Column = 1;

            % Create CameraControlTab
            app.CameraControlTab = uitab(app.TabGroup2);
            app.CameraControlTab.Title = 'Camera Control';

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.CameraControlTab);
            app.GridLayout3.ColumnWidth = {'1x', 60, '1x', '1x'};

            % Create ExposureSlider
            app.ExposureSlider = uislider(app.GridLayout3);
            app.ExposureSlider.Limits = [-10 10];
            app.ExposureSlider.MajorTicks = [-10 -5 0 5 10];
            app.ExposureSlider.ValueChangingFcn = createCallbackFcn(app, @ExposureSliderValueChanging, true);
            app.ExposureSlider.Layout.Row = 1;
            app.ExposureSlider.Layout.Column = 3;
            app.ExposureSlider.Value = -5;

            % Create ExposureSliderLabel
            app.ExposureSliderLabel = uilabel(app.GridLayout3);
            app.ExposureSliderLabel.HorizontalAlignment = 'right';
            app.ExposureSliderLabel.Layout.Row = 1;
            app.ExposureSliderLabel.Layout.Column = 2;
            app.ExposureSliderLabel.Text = 'Exposure';

            % Create GainSliderLabel
            app.GainSliderLabel = uilabel(app.GridLayout3);
            app.GainSliderLabel.HorizontalAlignment = 'right';
            app.GainSliderLabel.Layout.Row = 2;
            app.GainSliderLabel.Layout.Column = 2;
            app.GainSliderLabel.Text = 'Gain';

            % Create GainSlider
            app.GainSlider = uislider(app.GridLayout3);
            app.GainSlider.Limits = [0 300];
            app.GainSlider.ValueChangingFcn = createCallbackFcn(app, @GainSliderValueChanging, true);
            app.GainSlider.Layout.Row = 2;
            app.GainSlider.Layout.Column = 3;

            % Create GridLayout7
            app.GridLayout7 = uigridlayout(app.GridLayout3);
            app.GridLayout7.ColumnWidth = {'1x'};
            app.GridLayout7.Layout.Row = 1;
            app.GridLayout7.Layout.Column = 1;

            % Create DeviceName
            app.DeviceName = uidropdown(app.GridLayout7);
            app.DeviceName.Items = {'1', '2', '3', '4', '5', ''};
            app.DeviceName.DropDownOpeningFcn = createCallbackFcn(app, @DeviceNameDropDownOpening, true);
            app.DeviceName.ValueChangedFcn = createCallbackFcn(app, @DeviceNameValueChanged, true);
            app.DeviceName.Layout.Row = 2;
            app.DeviceName.Layout.Column = 1;
            app.DeviceName.Value = '3';

            % Create DeviceLabel
            app.DeviceLabel = uilabel(app.GridLayout7);
            app.DeviceLabel.HorizontalAlignment = 'center';
            app.DeviceLabel.Layout.Row = 1;
            app.DeviceLabel.Layout.Column = 1;
            app.DeviceLabel.Text = 'Device #';

            % Create GridLayout8
            app.GridLayout8 = uigridlayout(app.GridLayout3);
            app.GridLayout8.ColumnWidth = {'1x'};
            app.GridLayout8.Layout.Row = 2;
            app.GridLayout8.Layout.Column = 1;

            % Create Resolution
            app.Resolution = uidropdown(app.GridLayout8);
            app.Resolution.DropDownOpeningFcn = createCallbackFcn(app, @ResolutionDropDownOpening, true);
            app.Resolution.ValueChangedFcn = createCallbackFcn(app, @ResolutionValueChanged, true);
            app.Resolution.Layout.Row = 2;
            app.Resolution.Layout.Column = 1;

            % Create ResolutionLabel
            app.ResolutionLabel = uilabel(app.GridLayout8);
            app.ResolutionLabel.HorizontalAlignment = 'center';
            app.ResolutionLabel.Layout.Row = 1;
            app.ResolutionLabel.Layout.Column = 1;
            app.ResolutionLabel.Text = 'Resolution';

            % Create MicroscopeObjectiveDropDownLabel
            app.MicroscopeObjectiveDropDownLabel = uilabel(app.GridLayout3);
            app.MicroscopeObjectiveDropDownLabel.HorizontalAlignment = 'center';
            app.MicroscopeObjectiveDropDownLabel.Layout.Row = 1;
            app.MicroscopeObjectiveDropDownLabel.Layout.Column = 4;
            app.MicroscopeObjectiveDropDownLabel.Text = {'Microscope'; 'Objective'};

            % Create MicroscopeObjectiveDropDown
            app.MicroscopeObjectiveDropDown = uidropdown(app.GridLayout3);
            app.MicroscopeObjectiveDropDown.Items = {'5x', '10x', '40x', '100x'};
            app.MicroscopeObjectiveDropDown.ItemsData = {'0.6047', '0.2402', '0.07269', '0.02939', '', ''};
            app.MicroscopeObjectiveDropDown.Layout.Row = 2;
            app.MicroscopeObjectiveDropDown.Layout.Column = 4;
            app.MicroscopeObjectiveDropDown.Value = '0.2402';

            % Create MotionDetectionTab
            app.MotionDetectionTab = uitab(app.TabGroup2);
            app.MotionDetectionTab.Title = 'Motion Detection';

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.MotionDetectionTab);
            app.GridLayout4.ColumnWidth = {'1x', '1x', '2.5x'};

            % Create PickcircumferenceButton
            app.PickcircumferenceButton = uibutton(app.GridLayout4, 'push');
            app.PickcircumferenceButton.ButtonPushedFcn = createCallbackFcn(app, @PickcircumferenceButtonPushed, true);
            app.PickcircumferenceButton.Enable = 'off';
            app.PickcircumferenceButton.Layout.Row = 1;
            app.PickcircumferenceButton.Layout.Column = 1;
            app.PickcircumferenceButton.Text = {'Pick '; 'circumference'};

            % Create ChoosethetarangeButton
            app.ChoosethetarangeButton = uibutton(app.GridLayout4, 'push');
            app.ChoosethetarangeButton.ButtonPushedFcn = createCallbackFcn(app, @ChoosethetarangeButtonPushed, true);
            app.ChoosethetarangeButton.Enable = 'off';
            app.ChoosethetarangeButton.Layout.Row = 1;
            app.ChoosethetarangeButton.Layout.Column = 2;
            app.ChoosethetarangeButton.Text = {'Choose'; 'theta range'};

            % Create ChooseradialrangeButton
            app.ChooseradialrangeButton = uibutton(app.GridLayout4, 'push');
            app.ChooseradialrangeButton.ButtonPushedFcn = createCallbackFcn(app, @ChooseradialrangeButtonPushed, true);
            app.ChooseradialrangeButton.Enable = 'off';
            app.ChooseradialrangeButton.Layout.Row = 2;
            app.ChooseradialrangeButton.Layout.Column = 2;
            app.ChooseradialrangeButton.Text = {'Choose'; 'radial range'};

            % Create AlignpeakCheckBox
            app.AlignpeakCheckBox = uicheckbox(app.GridLayout4);
            app.AlignpeakCheckBox.ValueChangedFcn = createCallbackFcn(app, @AlignpeakCheckBoxValueChanged, true);
            app.AlignpeakCheckBox.Text = {'Align'; 'peak'};
            app.AlignpeakCheckBox.Layout.Row = 2;
            app.AlignpeakCheckBox.Layout.Column = 1;

            % Create GridLayout10
            app.GridLayout10 = uigridlayout(app.GridLayout4);
            app.GridLayout10.ColumnWidth = {'2x', '1x'};
            app.GridLayout10.RowHeight = {'1x'};
            app.GridLayout10.Layout.Row = 1;
            app.GridLayout10.Layout.Column = 3;

            % Create radialpitchumEditFieldLabel
            app.radialpitchumEditFieldLabel = uilabel(app.GridLayout10);
            app.radialpitchumEditFieldLabel.HorizontalAlignment = 'center';
            app.radialpitchumEditFieldLabel.Layout.Row = 1;
            app.radialpitchumEditFieldLabel.Layout.Column = 1;
            app.radialpitchumEditFieldLabel.Text = {'radial pitch'; '(um)'};

            % Create radialpitchumEditField
            app.radialpitchumEditField = uieditfield(app.GridLayout10, 'numeric');
            app.radialpitchumEditField.Layout.Row = 1;
            app.radialpitchumEditField.Layout.Column = 2;
            app.radialpitchumEditField.Value = 8;

            % Create MiscTab
            app.MiscTab = uitab(app.TabGroup2);
            app.MiscTab.Title = 'Misc.';

            % Create GridLayout9
            app.GridLayout9 = uigridlayout(app.MiscTab);
            app.GridLayout9.ColumnWidth = {'1.5x', '1x', '1.5x', '1x', '1x'};

            % Create SyncSnapshotinworkspaceCheckBox
            app.SyncSnapshotinworkspaceCheckBox = uicheckbox(app.GridLayout9);
            app.SyncSnapshotinworkspaceCheckBox.Text = {'Sync Snapshot in '; 'workspace'};
            app.SyncSnapshotinworkspaceCheckBox.Layout.Row = 1;
            app.SyncSnapshotinworkspaceCheckBox.Layout.Column = [1 2];
            app.SyncSnapshotinworkspaceCheckBox.Value = true;

            % Create SavesnapshottofileButton
            app.SavesnapshottofileButton = uibutton(app.GridLayout9, 'push');
            app.SavesnapshottofileButton.ButtonPushedFcn = createCallbackFcn(app, @SavesnapshottofileButtonPushed, true);
            app.SavesnapshottofileButton.Layout.Row = [1 2];
            app.SavesnapshottofileButton.Layout.Column = 2;
            app.SavesnapshottofileButton.Text = {'Save'; 'snapshot'; 'to file'};

            % Create GridLayout11
            app.GridLayout11 = uigridlayout(app.GridLayout9);
            app.GridLayout11.ColumnWidth = {'1.3x', '1x', '1x'};
            app.GridLayout11.RowHeight = {'1x'};
            app.GridLayout11.Layout.Row = 1;
            app.GridLayout11.Layout.Column = [3 4];

            % Create PickButton
            app.PickButton = uibutton(app.GridLayout11, 'push');
            app.PickButton.ButtonPushedFcn = createCallbackFcn(app, @PickButtonPushed, true);
            app.PickButton.Layout.Row = 1;
            app.PickButton.Layout.Column = 3;
            app.PickButton.Text = 'Pick';

            % Create ReferenceDesignEditField
            app.ReferenceDesignEditField = uieditfield(app.GridLayout11, 'text');
            app.ReferenceDesignEditField.Layout.Row = 1;
            app.ReferenceDesignEditField.Layout.Column = 2;

            % Create ReferenceDesignEditFieldLabel
            app.ReferenceDesignEditFieldLabel = uilabel(app.GridLayout11);
            app.ReferenceDesignEditFieldLabel.HorizontalAlignment = 'right';
            app.ReferenceDesignEditFieldLabel.Layout.Row = 1;
            app.ReferenceDesignEditFieldLabel.Layout.Column = 1;
            app.ReferenceDesignEditFieldLabel.Text = {'Reference'; 'Design'};

            % Create AlignsnapshottoRefDesignExperimentalButton
            app.AlignsnapshottoRefDesignExperimentalButton = uibutton(app.GridLayout9, 'push');
            app.AlignsnapshottoRefDesignExperimentalButton.ButtonPushedFcn = createCallbackFcn(app, @AlignsnapshottoRefDesignExperimentalButtonPushed, true);
            app.AlignsnapshottoRefDesignExperimentalButton.Layout.Row = 2;
            app.AlignsnapshottoRefDesignExperimentalButton.Layout.Column = 4;
            app.AlignsnapshottoRefDesignExperimentalButton.Text = {'Align snapshot'; 'to Ref. Design'; '(Experimental)'};

            % Create PreviewRefDesignButton
            app.PreviewRefDesignButton = uibutton(app.GridLayout9, 'push');
            app.PreviewRefDesignButton.ButtonPushedFcn = createCallbackFcn(app, @PreviewRefDesignButtonPushed, true);
            app.PreviewRefDesignButton.Layout.Row = 1;
            app.PreviewRefDesignButton.Layout.Column = 5;
            app.PreviewRefDesignButton.Text = {'Preview'; 'Ref. Design'};

            % Create GridLayout12
            app.GridLayout12 = uigridlayout(app.GridLayout9);
            app.GridLayout12.RowHeight = {'1x'};
            app.GridLayout12.Layout.Row = 2;
            app.GridLayout12.Layout.Column = 1;

            % Create VarNameLabel
            app.VarNameLabel = uilabel(app.GridLayout12);
            app.VarNameLabel.HorizontalAlignment = 'right';
            app.VarNameLabel.Layout.Row = 1;
            app.VarNameLabel.Layout.Column = 1;
            app.VarNameLabel.Text = {'Var'; 'Name'};

            % Create VariableNameEditField
            app.VariableNameEditField = uieditfield(app.GridLayout12, 'text');
            app.VariableNameEditField.Layout.Row = 1;
            app.VariableNameEditField.Layout.Column = 2;
            app.VariableNameEditField.Value = 'Fsnap';

            % Create UseGPUCheckBox
            app.UseGPUCheckBox = uicheckbox(app.GridLayout9);
            app.UseGPUCheckBox.Text = 'Use GPU';
            app.UseGPUCheckBox.Layout.Row = 2;
            app.UseGPUCheckBox.Layout.Column = 5;

            % Create GridLayout13
            app.GridLayout13 = uigridlayout(app.GridLayout9);
            app.GridLayout13.RowHeight = {'1x'};
            app.GridLayout13.Layout.Row = 2;
            app.GridLayout13.Layout.Column = 3;

            % Create MinScaleEditFieldLabel
            app.MinScaleEditFieldLabel = uilabel(app.GridLayout13);
            app.MinScaleEditFieldLabel.HorizontalAlignment = 'right';
            app.MinScaleEditFieldLabel.Layout.Row = 1;
            app.MinScaleEditFieldLabel.Layout.Column = 1;
            app.MinScaleEditFieldLabel.Text = {'Min'; 'Scale'};

            % Create MinScaleEditField
            app.MinScaleEditField = uieditfield(app.GridLayout13, 'numeric');
            app.MinScaleEditField.Limits = [1 16];
            app.MinScaleEditField.RoundFractionalValues = 'on';
            app.MinScaleEditField.Layout.Row = 1;
            app.MinScaleEditField.Layout.Column = 2;
            app.MinScaleEditField.Value = 2;

            % Create ImageMotionDetectionPanel
            % app.ImageMotionDetectionPanel = uipanel(app.GridLayout);
            % app.ImageMotionDetectionPanel.TitlePosition = 'centertop';
            % app.ImageMotionDetectionPanel.Title = 'Image Motion Detection';
            % app.ImageMotionDetectionPanel.Layout.Row = 2;
            % app.ImageMotionDetectionPanel.Layout.Column = [3 4];

            % [FIX] The modified code (the second section with [FIX] annotations) follows the constraint of MATLAB R2024 and ensures compatibility with both versions: 
            % 1. First, create the uipanel and directly bind it to the parent app.GridLayout at the time of creation (specify the parent as the first input argument), while setting basic properties (Title, BorderType, TitlePosition) simultaneously.
            app.ImageMotionDetectionPanel = uipanel(app.GridLayout, 'Title', 'Image Motion Detection', 'BorderType', 'line', 'TitlePosition', 'centertop');
            % 2. After the uipanel is successfully created and associated with its parent GridLayout, set the Layout.Row and Layout.Column properties to define the component's position in the GridLayout.
            app.ImageMotionDetectionPanel.Layout.Row = 2;
            app.ImageMotionDetectionPanel.Layout.Column = [3 4];

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.ImageMotionDetectionPanel);
            app.GridLayout2.ColumnWidth = {'1x', '1x', '1x', 50};

            % Create UIAxes5
            app.UIAxes5 = uiaxes(app.GridLayout2);
            title(app.UIAxes5, '\theta Displacement')
            xlabel(app.UIAxes5, 'Frame')
            ylabel(app.UIAxes5, '\Delta\theta (^\circ)')
            zlabel(app.UIAxes5, 'Z')
            app.UIAxes5.Box = 'on';
            app.UIAxes5.Layout.Row = 1;
            app.UIAxes5.Layout.Column = 2;

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.GridLayout2);
            title(app.UIAxes4, 'Radial Displacement')
            xlabel(app.UIAxes4, 'Frame')
            ylabel(app.UIAxes4, '\DeltaR (um)')
            app.UIAxes4.Box = 'on';
            app.UIAxes4.Layout.Row = 2;
            app.UIAxes4.Layout.Column = 2;

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.GridLayout2);
            title(app.UIAxes3, '\theta Detection')
            xlabel(app.UIAxes3, '\theta')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.YTick = [];
            app.UIAxes3.Box = 'on';
            app.UIAxes3.Layout.Row = 1;
            app.UIAxes3.Layout.Column = 3;

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.GridLayout2);
            title(app.UIAxes2, 'De-polar view')
            app.UIAxes2.XTick = [];
            app.UIAxes2.YTick = [];
            app.UIAxes2.BoxStyle = 'full';
            app.UIAxes2.Box = 'on';
            app.UIAxes2.Layout.Row = [1 2];
            app.UIAxes2.Layout.Column = 1;

            % Create UIAxes6
            app.UIAxes6 = uiaxes(app.GridLayout2);
            title(app.UIAxes6, 'Radial Detection')
            xlabel(app.UIAxes6, '\theta')
            zlabel(app.UIAxes6, 'Z')
            app.UIAxes6.YTick = [];
            app.UIAxes6.Box = 'on';
            app.UIAxes6.Layout.Row = 2;
            app.UIAxes6.Layout.Column = 3;

            % Create ClearGraphsButton
            app.ClearGraphsButton = uibutton(app.GridLayout2, 'push');
            app.ClearGraphsButton.ButtonPushedFcn = createCallbackFcn(app, @ClearGraphsButtonPushed, true);
            app.ClearGraphsButton.Layout.Row = [1 2];
            app.ClearGraphsButton.Layout.Column = 4;
            app.ClearGraphsButton.Text = {'Clear'; 'Graphs'};

            % Create MEMScontrolreadoutPanel
            app.MEMScontrolreadoutPanel = uipanel(app.GridLayout);
            app.MEMScontrolreadoutPanel.TitlePosition = 'centertop';
            app.MEMScontrolreadoutPanel.Title = 'MEMS control & readout';
            app.MEMScontrolreadoutPanel.Layout.Row = 1;
            app.MEMScontrolreadoutPanel.Layout.Column = [3 4];

            % Create GridLayout15
            app.GridLayout15 = uigridlayout(app.MEMScontrolreadoutPanel);
            app.GridLayout15.ColumnWidth = {'0.5x', '0.5x', '3.01x', '0.5x', '0.5x'};
            app.GridLayout15.RowHeight = {'1x', '1x', '1x', '1x', '1x'};
            app.GridLayout15.RowSpacing = 5.5;
            app.GridLayout15.Padding = [10 5.5 10 5.5];

            % Create SR830DropDown
            app.SR830DropDown = uidropdown(app.GridLayout15);
            app.SR830DropDown.Items = {'a', 'b'};
            app.SR830DropDown.DropDownOpeningFcn = createCallbackFcn(app, @SR830DropDownOpening, true);
            app.SR830DropDown.Layout.Row = 2;
            app.SR830DropDown.Layout.Column = 1;
            app.SR830DropDown.Value = 'a';

            % Create SR830ConnectButton
            app.SR830ConnectButton = uibutton(app.GridLayout15, 'push');
            app.SR830ConnectButton.ButtonPushedFcn = createCallbackFcn(app, @SR830ConnectButtonPushed, true);
            app.SR830ConnectButton.Layout.Row = 3;
            app.SR830ConnectButton.Layout.Column = 1;
            app.SR830ConnectButton.Text = 'Connect';

            % Create TabGroup3
            app.TabGroup3 = uitabgroup(app.GridLayout15);
            app.TabGroup3.Layout.Row = [1 5];
            app.TabGroup3.Layout.Column = 3;

            % Create DirectTab
            app.DirectTab = uitab(app.TabGroup3);
            app.DirectTab.Title = 'Direct';

            % Create SetButton
            app.SetButton = uibutton(app.DirectTab, 'push');
            app.SetButton.ButtonPushedFcn = createCallbackFcn(app, @SetButtonPushed, true);
            app.SetButton.Enable = 'off';
            app.SetButton.Position = [471 22 92 58];
            app.SetButton.Text = 'Set';

            % Create V1Knob
            app.V1Knob = uiknob(app.DirectTab, 'continuous');
            app.V1Knob.Limits = [-250 250];
            app.V1Knob.ValueChangingFcn = createCallbackFcn(app, @V1KnobValueChanging, true);
            app.V1Knob.Position = [52 27 60 60];

            % Create V1KnobLabel
            app.V1KnobLabel = uilabel(app.DirectTab);
            app.V1KnobLabel.HorizontalAlignment = 'center';
            app.V1KnobLabel.Position = [69 7 25 22];
            app.V1KnobLabel.Text = 'V1';

            % Create V2Knob
            app.V2Knob = uiknob(app.DirectTab, 'continuous');
            app.V2Knob.Limits = [-250 250];
            app.V2Knob.ValueChangingFcn = createCallbackFcn(app, @V2KnobValueChanging, true);
            app.V2Knob.Position = [196 25 60 60];

            % Create V2KnobLabel
            app.V2KnobLabel = uilabel(app.DirectTab);
            app.V2KnobLabel.HorizontalAlignment = 'center';
            app.V2KnobLabel.Position = [210 3 25 22];
            app.V2KnobLabel.Text = 'V2';

            % Create V3Knob
            app.V3Knob = uiknob(app.DirectTab, 'continuous');
            app.V3Knob.Limits = [-250 250];
            app.V3Knob.ValueChangingFcn = createCallbackFcn(app, @V3KnobValueChanging, true);
            app.V3Knob.Position = [337 26 60 60];

            % Create V3KnobLabel
            app.V3KnobLabel = uilabel(app.DirectTab);
            app.V3KnobLabel.HorizontalAlignment = 'center';
            app.V3KnobLabel.Position = [354 6 25 22];
            app.V3KnobLabel.Text = 'V3';

            % Create V3Label
            app.V3Label = uieditfield(app.DirectTab, 'numeric');
            app.V3Label.ValueChangedFcn = createCallbackFcn(app, @V3LabelValueChanged, true);
            app.V3Label.HorizontalAlignment = 'center';
            app.V3Label.FontSize = 9;
            app.V3Label.Position = [350 48 33 16];

            % Create V2Label
            app.V2Label = uieditfield(app.DirectTab, 'numeric');
            app.V2Label.ValueChangedFcn = createCallbackFcn(app, @V2LabelValueChanged, true);
            app.V2Label.HorizontalAlignment = 'center';
            app.V2Label.FontSize = 9;
            app.V2Label.Position = [207 47 33 16];

            % Create V1Label
            app.V1Label = uieditfield(app.DirectTab, 'numeric');
            app.V1Label.ValueChangedFcn = createCallbackFcn(app, @V1LabelValueChanged, true);
            app.V1Label.HorizontalAlignment = 'center';
            app.V1Label.FontSize = 9;
            app.V1Label.Position = [65 46 33 16];

            % Create SteppingTab
            app.SteppingTab = uitab(app.TabGroup3);
            app.SteppingTab.Title = 'Stepping';

            % Create InitButton
            app.InitButton = uibutton(app.SteppingTab, 'push');
            app.InitButton.ButtonPushedFcn = createCallbackFcn(app, @InitButtonPushed, true);
            app.InitButton.Enable = 'off';
            app.InitButton.Position = [29 11 53 44];
            app.InitButton.Text = 'Init';

            % Create Plus18Button
            app.Plus18Button = uibutton(app.SteppingTab, 'push');
            app.Plus18Button.ButtonPushedFcn = createCallbackFcn(app, @Plus18ButtonPushed, true);
            app.Plus18Button.Enable = 'off';
            app.Plus18Button.Position = [211 13 53 22];
            app.Plus18Button.Text = {'+0.125'; ''};

            % Create DriveVoltageEditField
            app.DriveVoltageEditField = uieditfield(app.SteppingTab, 'numeric');
            app.DriveVoltageEditField.Position = [83 71 41 24];
            app.DriveVoltageEditField.Value = 20;

            % Create DriveVoltageLabel
            app.DriveVoltageLabel = uilabel(app.SteppingTab);
            app.DriveVoltageLabel.HorizontalAlignment = 'right';
            app.DriveVoltageLabel.Position = [29 63 46 40];
            app.DriveVoltageLabel.Text = {'Drive'; 'Voltage'};

            % Create PhaseEditFieldLabel
            app.PhaseEditFieldLabel = uilabel(app.SteppingTab);
            app.PhaseEditFieldLabel.HorizontalAlignment = 'right';
            app.PhaseEditFieldLabel.FontSize = 18;
            app.PhaseEditFieldLabel.Position = [260 65 79 26];
            app.PhaseEditFieldLabel.Text = 'Phase';

            % Create PhaseEditField
            app.PhaseEditField = uieditfield(app.SteppingTab, 'numeric');
            app.PhaseEditField.Editable = 'off';
            app.PhaseEditField.FontSize = 20;
            app.PhaseEditField.Position = [356 68 56 25];

            % Create StepsEditFieldLabel
            app.StepsEditFieldLabel = uilabel(app.SteppingTab);
            app.StepsEditFieldLabel.HorizontalAlignment = 'right';
            app.StepsEditFieldLabel.FontSize = 18;
            app.StepsEditFieldLabel.Position = [282 37 52 22];
            app.StepsEditFieldLabel.Text = 'Steps';

            % Create StepsEditField
            app.StepsEditField = uieditfield(app.SteppingTab, 'numeric');
            app.StepsEditField.Editable = 'off';
            app.StepsEditField.FontSize = 20;
            app.StepsEditField.Position = [356 36 56 25];

            % Create ZeroStepsButton
            app.ZeroStepsButton = uibutton(app.SteppingTab, 'push');
            app.ZeroStepsButton.ButtonPushedFcn = createCallbackFcn(app, @ZeroStepsButtonPushed, true);
            app.ZeroStepsButton.Enable = 'off';
            app.ZeroStepsButton.Position = [83 11 53 44];
            app.ZeroStepsButton.Text = {'Zero'; 'Steps'};

            % Create Minus18Button
            app.Minus18Button = uibutton(app.SteppingTab, 'push');
            app.Minus18Button.ButtonPushedFcn = createCallbackFcn(app, @Minus18ButtonPushed, true);
            app.Minus18Button.Enable = 'off';
            app.Minus18Button.Position = [152 13 53 22];
            app.Minus18Button.Text = {'-0.125'; ''};

            % Create PlusHalfButton
            app.PlusHalfButton = uibutton(app.SteppingTab, 'push');
            app.PlusHalfButton.ButtonPushedFcn = createCallbackFcn(app, @PlusHalfButtonPushed, true);
            app.PlusHalfButton.Enable = 'off';
            app.PlusHalfButton.Position = [211 71 53 22];
            app.PlusHalfButton.Text = '+0.5';

            % Create MinusHalfButton
            app.MinusHalfButton = uibutton(app.SteppingTab, 'push');
            app.MinusHalfButton.ButtonPushedFcn = createCallbackFcn(app, @MinusHalfButtonPushed, true);
            app.MinusHalfButton.Enable = 'off';
            app.MinusHalfButton.Position = [152 71 53 22];
            app.MinusHalfButton.Text = '-0.5';

            % Create PlusQuarterButton
            app.PlusQuarterButton = uibutton(app.SteppingTab, 'push');
            app.PlusQuarterButton.ButtonPushedFcn = createCallbackFcn(app, @PlusQuarterButtonPushed, true);
            app.PlusQuarterButton.Enable = 'off';
            app.PlusQuarterButton.Position = [211 43 53 22];
            app.PlusQuarterButton.Text = '+0.25';

            % Create MinusQuarterButton
            app.MinusQuarterButton = uibutton(app.SteppingTab, 'push');
            app.MinusQuarterButton.ButtonPushedFcn = createCallbackFcn(app, @MinusQuarterButtonPushed, true);
            app.MinusQuarterButton.Enable = 'off';
            app.MinusQuarterButton.Position = [152 43 53 22];
            app.MinusQuarterButton.Text = '-0.25';

            % Create SyncCapFittingCheckBox
            app.SyncCapFittingCheckBox = uicheckbox(app.SteppingTab);
            app.SyncCapFittingCheckBox.Text = 'Sync Cap Fitting';
            app.SyncCapFittingCheckBox.Position = [443 69 110 22];

            % Create PositioningTab
            app.PositioningTab = uitab(app.TabGroup3);
            app.PositioningTab.Title = 'Positioning';

            % Create DriveVoltageEditField_2
            app.DriveVoltageEditField_2 = uieditfield(app.PositioningTab, 'numeric');
            app.DriveVoltageEditField_2.Position = [62 75 41 24];
            app.DriveVoltageEditField_2.Value = 50;

            % Create DriveVoltageEditField_2Label
            app.DriveVoltageEditField_2Label = uilabel(app.PositioningTab);
            app.DriveVoltageEditField_2Label.HorizontalAlignment = 'right';
            app.DriveVoltageEditField_2Label.Position = [8 67 46 40];
            app.DriveVoltageEditField_2Label.Text = {'Drive'; 'Voltage'};

            % Create InitResetButton
            app.InitResetButton = uibutton(app.PositioningTab, 'push');
            app.InitResetButton.ButtonPushedFcn = createCallbackFcn(app, @InitResetButtonPushed, true);
            app.InitResetButton.BackgroundColor = [0.149 0.149 0.149];
            app.InitResetButton.FontColor = [1 1 1];
            app.InitResetButton.Position = [26 8 80 44];
            app.InitResetButton.Text = 'Init/Reset';

            % Create TargetStepEditField
            app.TargetStepEditField = uieditfield(app.PositioningTab, 'numeric');
            app.TargetStepEditField.Position = [162 75 41 24];

            % Create TargetStepEditFieldLabel
            app.TargetStepEditFieldLabel = uilabel(app.PositioningTab);
            app.TargetStepEditFieldLabel.HorizontalAlignment = 'right';
            app.TargetStepEditFieldLabel.Position = [108 68 46 40];
            app.TargetStepEditFieldLabel.Text = {'Target'; 'Step'};

            % Create StepsPhyEditFieldLabel
            app.StepsPhyEditFieldLabel = uilabel(app.PositioningTab);
            app.StepsPhyEditFieldLabel.HorizontalAlignment = 'right';
            app.StepsPhyEditFieldLabel.FontSize = 16;
            app.StepsPhyEditFieldLabel.Position = [427 83 89 22];
            app.StepsPhyEditFieldLabel.Text = 'Steps (Phy)';

            % Create StepsPhyEditField
            app.StepsPhyEditField = uieditfield(app.PositioningTab, 'numeric');
            app.StepsPhyEditField.Editable = 'off';
            app.StepsPhyEditField.HorizontalAlignment = 'left';
            app.StepsPhyEditField.FontSize = 20;
            app.StepsPhyEditField.Position = [527 80 74 25];

            % Create SpeedstepsEditField
            app.SpeedstepsEditField = uieditfield(app.PositioningTab, 'numeric');
            app.SpeedstepsEditField.Position = [162 43 41 24];
            app.SpeedstepsEditField.Value = 1;

            % Create SpeedstepsEditFieldLabel
            app.SpeedstepsEditFieldLabel = uilabel(app.PositioningTab);
            app.SpeedstepsEditFieldLabel.HorizontalAlignment = 'right';
            app.SpeedstepsEditFieldLabel.Position = [108 35 46 40];
            app.SpeedstepsEditFieldLabel.Text = {'Speed'; 'step/s'};

            % Create GOButton
            app.GOButton = uibutton(app.PositioningTab, 'push');
            app.GOButton.ButtonPushedFcn = createCallbackFcn(app, @GOButtonPushed, true);
            app.GOButton.Enable = 'off';
            app.GOButton.Position = [241 8 80 38];
            app.GOButton.Text = 'GO';

            % Create MicroStepDropDownLabel
            app.MicroStepDropDownLabel = uilabel(app.PositioningTab);
            app.MicroStepDropDownLabel.HorizontalAlignment = 'right';
            app.MicroStepDropDownLabel.Position = [94 2 58 30];
            app.MicroStepDropDownLabel.Text = {'Micro'; 'Step'};

            % Create MicroStepDropDown
            app.MicroStepDropDown = uidropdown(app.PositioningTab);
            app.MicroStepDropDown.Items = {'2', '4', '8', '16', '32', '64', '128', '256'};
            app.MicroStepDropDown.Position = [162 4 67 22];
            app.MicroStepDropDown.Value = '16';

            % Create ReadChannelDropDownLabel
            app.ReadChannelDropDownLabel = uilabel(app.PositioningTab);
            app.ReadChannelDropDownLabel.HorizontalAlignment = 'right';
            app.ReadChannelDropDownLabel.Position = [471 8 50 30];
            app.ReadChannelDropDownLabel.Text = {'Read'; 'Channel'};

            % Create ReadChannelDropDown
            app.ReadChannelDropDown = uidropdown(app.PositioningTab);
            app.ReadChannelDropDown.Items = {'X', 'Y', 'R', 'Theta'};
            app.ReadChannelDropDown.Position = [530 12 67 22];
            app.ReadChannelDropDown.Value = 'R';

            % Create LogCapDataButton
            app.LogCapDataButton = uibutton(app.PositioningTab, 'push');
            app.LogCapDataButton.ButtonPushedFcn = createCallbackFcn(app, @LogCapDataButtonPushed, true);
            app.LogCapDataButton.Enable = 'off';
            app.LogCapDataButton.Position = [241 60 82 38];
            app.LogCapDataButton.Text = 'LogCapData';

            % Create CloseLoopSwitchLabel
            app.CloseLoopSwitchLabel = uilabel(app.PositioningTab);
            app.CloseLoopSwitchLabel.HorizontalAlignment = 'center';
            app.CloseLoopSwitchLabel.Enable = 'off';
            app.CloseLoopSwitchLabel.Position = [362 24 66 22];
            app.CloseLoopSwitchLabel.Text = 'Close Loop';

            % Create CloseLoopSwitch
            app.CloseLoopSwitch = uiswitch(app.PositioningTab, 'slider');
            app.CloseLoopSwitch.Enable = 'off';
            app.CloseLoopSwitch.Position = [386 12 18 8];

            % Create StepsCapEditFieldLabel
            app.StepsCapEditFieldLabel = uilabel(app.PositioningTab);
            app.StepsCapEditFieldLabel.HorizontalAlignment = 'right';
            app.StepsCapEditFieldLabel.FontSize = 16;
            app.StepsCapEditFieldLabel.Position = [426 47 90 22];
            app.StepsCapEditFieldLabel.Text = 'Steps (Cap)';

            % Create StepsCapEditField
            app.StepsCapEditField = uieditfield(app.PositioningTab, 'numeric');
            app.StepsCapEditField.Editable = 'off';
            app.StepsCapEditField.HorizontalAlignment = 'left';
            app.StepsCapEditField.FontSize = 20;
            app.StepsCapEditField.Position = [527 45 74 25];
            app.StepsCapEditField.Value = -Inf;

            % Create ZeroCapButton
            app.ZeroCapButton = uibutton(app.PositioningTab, 'push');
            app.ZeroCapButton.ButtonPushedFcn = createCallbackFcn(app, @ZeroCapButtonPushed, true);
            app.ZeroCapButton.Enable = 'off';
            app.ZeroCapButton.Position = [334 59 82 38];
            app.ZeroCapButton.Text = 'Zero Cap';

            % Create PositioningFeedbackTab
            app.PositioningFeedbackTab = uitab(app.TabGroup3);
            app.PositioningFeedbackTab.Title = 'Positioning (Feedback)';

            % Create DriveVoltageEditField_3
            app.DriveVoltageEditField_3 = uieditfield(app.PositioningFeedbackTab, 'numeric');
            app.DriveVoltageEditField_3.HorizontalAlignment = 'center';
            app.DriveVoltageEditField_3.Position = [52 83 26 24];
            app.DriveVoltageEditField_3.Value = 50;

            % Create DriveVoltageEditField_3Label
            app.DriveVoltageEditField_3Label = uilabel(app.PositioningFeedbackTab);
            app.DriveVoltageEditField_3Label.HorizontalAlignment = 'center';
            app.DriveVoltageEditField_3Label.Position = [2 83 46 28];
            app.DriveVoltageEditField_3Label.Text = {'Drive'; 'Voltage'};

            % Create InitResetButton_2
            app.InitResetButton_2 = uibutton(app.PositioningFeedbackTab, 'push');
            app.InitResetButton_2.ButtonPushedFcn = createCallbackFcn(app, @InitResetButton_2Pushed, true);
            app.InitResetButton_2.BackgroundColor = [0 0 0];
            app.InitResetButton_2.FontColor = [1 1 1];
            app.InitResetButton_2.Position = [6 8 88 33];
            app.InitResetButton_2.Text = 'Init/Reset';

            % Create ScanStepSrt
            app.ScanStepSrt = uieditfield(app.PositioningFeedbackTab, 'numeric');
            app.ScanStepSrt.HorizontalAlignment = 'center';
            app.ScanStepSrt.Position = [125 83 26 23];

            % Create ScanStepEditFieldLabel
            app.ScanStepEditFieldLabel = uilabel(app.PositioningFeedbackTab);
            app.ScanStepEditFieldLabel.HorizontalAlignment = 'center';
            app.ScanStepEditFieldLabel.Position = [83 80 33 30];
            app.ScanStepEditFieldLabel.Text = {'Scan'; 'Step'};

            % Create SpeedstepsEditField_2
            app.SpeedstepsEditField_2 = uieditfield(app.PositioningFeedbackTab, 'numeric');
            app.SpeedstepsEditField_2.HorizontalAlignment = 'center';
            app.SpeedstepsEditField_2.Position = [52 49 26 23];
            app.SpeedstepsEditField_2.Value = 2;

            % Create SpeedstepsEditField_2Label
            app.SpeedstepsEditField_2Label = uilabel(app.PositioningFeedbackTab);
            app.SpeedstepsEditField_2Label.HorizontalAlignment = 'center';
            app.SpeedstepsEditField_2Label.Position = [2 41 46 40];
            app.SpeedstepsEditField_2Label.Text = {'Speed'; '(step/s)'};

            % Create CurrentStepLabel
            app.CurrentStepLabel = uilabel(app.PositioningFeedbackTab);
            app.CurrentStepLabel.HorizontalAlignment = 'center';
            app.CurrentStepLabel.FontWeight = 'bold';
            app.CurrentStepLabel.FontColor = [0 0 0];
            app.CurrentStepLabel.Position = [496 62 102 22];
            app.CurrentStepLabel.Text = 'Current Step';

            % Create CurrentStepEditField
            app.CurrentStepEditField = uieditfield(app.PositioningFeedbackTab, 'numeric');
            app.CurrentStepEditField.Editable = 'off';
            app.CurrentStepEditField.HorizontalAlignment = 'center';
            app.CurrentStepEditField.FontWeight = 'bold';
            app.CurrentStepEditField.FontColor = [0 0 0];
            app.CurrentStepEditField.Position = [508 81 79 23];

            % Create GOButton_2
            app.GOButton_2 = uibutton(app.PositioningFeedbackTab, 'push');
            app.GOButton_2.ButtonPushedFcn = createCallbackFcn(app, @GOButton_2Pushed, true);
            app.GOButton_2.BackgroundColor = [0 0 1];
            app.GOButton_2.FontColor = [1 1 1];
            app.GOButton_2.Enable = 'off';
            app.GOButton_2.Position = [515 6 55 36];
            app.GOButton_2.Text = 'GO';

            % Create FinalStepEditField
            app.FinalStepEditField = uieditfield(app.PositioningFeedbackTab, 'numeric');
            app.FinalStepEditField.HorizontalAlignment = 'center';
            app.FinalStepEditField.Position = [247 83 36 23];

            % Create FinalStepEditFieldLabel
            app.FinalStepEditFieldLabel = uilabel(app.PositioningFeedbackTab);
            app.FinalStepEditFieldLabel.HorizontalAlignment = 'center';
            app.FinalStepEditFieldLabel.Position = [212 75 32 40];
            app.FinalStepEditFieldLabel.Text = {'Final'; 'Step'};

            % Create StepMinLossEditFieldLabel
            app.StepMinLossEditFieldLabel = uilabel(app.PositioningFeedbackTab);
            app.StepMinLossEditFieldLabel.HorizontalAlignment = 'center';
            app.StepMinLossEditFieldLabel.Position = [224 39 60 40];
            app.StepMinLossEditFieldLabel.Text = {'Step'; 'Min.Loss'};

            % Create StepMinLossEditField
            app.StepMinLossEditField = uieditfield(app.PositioningFeedbackTab, 'numeric');
            app.StepMinLossEditField.HorizontalAlignment = 'center';
            app.StepMinLossEditField.Position = [282 45 69 25];
            app.StepMinLossEditField.Value = 0.0001;

            % Create FinalMicrostepLabel
            app.FinalMicrostepLabel = uilabel(app.PositioningFeedbackTab);
            app.FinalMicrostepLabel.HorizontalAlignment = 'center';
            app.FinalMicrostepLabel.Position = [224 7 57 30];
            app.FinalMicrostepLabel.Text = {'PID'; 'Microstep'};

            % Create PIDMicrostepDropDown
            app.PIDMicrostepDropDown = uidropdown(app.PositioningFeedbackTab);
            app.PIDMicrostepDropDown.Items = {'2', '4', '8', '16', '32', '64', '128', '256', 'inf'};
            app.PIDMicrostepDropDown.Position = [286 10 64 24];
            app.PIDMicrostepDropDown.Value = 'inf';

            % Create CANCELButton
            app.CANCELButton = uibutton(app.PositioningFeedbackTab, 'push');
            app.CANCELButton.ButtonPushedFcn = createCallbackFcn(app, @CANCELButtonPushed, true);
            app.CANCELButton.BackgroundColor = [1 0 0];
            app.CANCELButton.FontColor = [1 1 1];
            app.CANCELButton.Enable = 'off';
            app.CANCELButton.Position = [443 65 55 40];
            app.CANCELButton.Text = 'CANCEL';

            % Create FilterWinEditFieldLabel
            app.FilterWinEditFieldLabel = uilabel(app.PositioningFeedbackTab);
            app.FilterWinEditFieldLabel.HorizontalAlignment = 'center';
            app.FilterWinEditFieldLabel.Position = [286 74 48 40];
            app.FilterWinEditFieldLabel.Text = {'Filter'; 'Win.'};

            % Create FilterWinEditField
            app.FilterWinEditField = uieditfield(app.PositioningFeedbackTab, 'numeric');
            app.FilterWinEditField.HorizontalAlignment = 'center';
            app.FilterWinEditField.Position = [327 83 25 23];
            app.FilterWinEditField.Value = 1;

            % Create ScanMicrostepLabel
            app.ScanMicrostepLabel = uilabel(app.PositioningFeedbackTab);
            app.ScanMicrostepLabel.HorizontalAlignment = 'center';
            app.ScanMicrostepLabel.Position = [105 7 55 30];
            app.ScanMicrostepLabel.Text = {'Scan'; 'Microstep'};

            % Create ScanMicrostepDropDown
            app.ScanMicrostepDropDown = uidropdown(app.PositioningFeedbackTab);
            app.ScanMicrostepDropDown.Items = {'2', '4', '8', '16', '32', '64', '128', '256', 'inf'};
            app.ScanMicrostepDropDown.Position = [165 9 59 24];
            app.ScanMicrostepDropDown.Value = '2';

            % Create ScanStepEnd
            app.ScanStepEnd = uieditfield(app.PositioningFeedbackTab, 'numeric');
            app.ScanStepEnd.HorizontalAlignment = 'center';
            app.ScanStepEnd.Position = [165 83 26 23];
            app.ScanStepEnd.Value = 7;

            % Create ScanStepEditFieldLabel_2
            app.ScanStepEditFieldLabel_2 = uilabel(app.PositioningFeedbackTab);
            app.ScanStepEditFieldLabel_2.HorizontalAlignment = 'center';
            app.ScanStepEditFieldLabel_2.Position = [151 80 14 30];
            app.ScanStepEditFieldLabel_2.Text = '~';

            % Create SCANButton
            app.SCANButton = uibutton(app.PositioningFeedbackTab, 'push');
            app.SCANButton.ButtonPushedFcn = createCallbackFcn(app, @SCANButtonPushed, true);
            app.SCANButton.BackgroundColor = [0 0 1];
            app.SCANButton.FontColor = [1 1 1];
            app.SCANButton.Enable = 'off';
            app.SCANButton.Position = [443 7 55 36];
            app.SCANButton.Text = 'SCAN';

            % Create KpEditFieldLabel
            app.KpEditFieldLabel = uilabel(app.PositioningFeedbackTab);
            app.KpEditFieldLabel.HorizontalAlignment = 'center';
            app.KpEditFieldLabel.Position = [356 76 22 40];
            app.KpEditFieldLabel.Text = 'Kp';

            % Create KpEditField
            app.KpEditField = uieditfield(app.PositioningFeedbackTab, 'numeric');
            app.KpEditField.HorizontalAlignment = 'center';
            app.KpEditField.Position = [378 83 49 23];
            app.KpEditField.Value = 0.35;

            % Create KiEditFieldLabel
            app.KiEditFieldLabel = uilabel(app.PositioningFeedbackTab);
            app.KiEditFieldLabel.HorizontalAlignment = 'center';
            app.KiEditFieldLabel.Position = [356 39 25 40];
            app.KiEditFieldLabel.Text = 'Ki';

            % Create KiEditField
            app.KiEditField = uieditfield(app.PositioningFeedbackTab, 'numeric');
            app.KiEditField.HorizontalAlignment = 'center';
            app.KiEditField.Position = [378 46 50 23];
            app.KiEditField.Value = 0.0005;

            % Create KdEditFieldLabel
            app.KdEditFieldLabel = uilabel(app.PositioningFeedbackTab);
            app.KdEditFieldLabel.HorizontalAlignment = 'center';
            app.KdEditFieldLabel.Position = [356 2 25 40];
            app.KdEditFieldLabel.Text = 'Kd';

            % Create KdEditField
            app.KdEditField = uieditfield(app.PositioningFeedbackTab, 'numeric');
            app.KdEditField.HorizontalAlignment = 'center';
            app.KdEditField.Position = [377 9 50 23];
            app.KdEditField.Value = 0.05;

            % Create LogPlotSwitchLabel
            app.LogPlotSwitchLabel = uilabel(app.PositioningFeedbackTab);
            app.LogPlotSwitchLabel.HorizontalAlignment = 'center';
            app.LogPlotSwitchLabel.Position = [165 56 49 22];
            app.LogPlotSwitchLabel.Text = 'Log/Plot';

            % Create LogPlotSwitch
            app.LogPlotSwitch = uiswitch(app.PositioningFeedbackTab, 'slider');
            app.LogPlotSwitch.Position = [182 45 16 7];
            app.LogPlotSwitch.Value = 'On';

            % Create PIDwaitsEditField
            app.PIDwaitsEditField = uieditfield(app.PositioningFeedbackTab, 'numeric');
            app.PIDwaitsEditField.HorizontalAlignment = 'center';
            app.PIDwaitsEditField.Position = [125 47 26 23];
            app.PIDwaitsEditField.Value = 0.5;

            % Create PIDwaitsEditFieldLabel
            app.PIDwaitsEditFieldLabel = uilabel(app.PositioningFeedbackTab);
            app.PIDwaitsEditFieldLabel.HorizontalAlignment = 'center';
            app.PIDwaitsEditFieldLabel.Position = [74 38 52 40];
            app.PIDwaitsEditFieldLabel.Text = {'PID'; 'wait (s)'};

            % Create OUTPUTOFFButton
            app.OUTPUTOFFButton = uibutton(app.GridLayout15, 'push');
            app.OUTPUTOFFButton.ButtonPushedFcn = createCallbackFcn(app, @OUTPUTOFFButtonPushed, true);
            app.OUTPUTOFFButton.FontColor = [1 0 0];
            app.OUTPUTOFFButton.Layout.Row = [4 5];
            app.OUTPUTOFFButton.Layout.Column = [1 2];
            app.OUTPUTOFFButton.Text = {'OUTPUT'; 'OFF'};

            % Create SR830Label
            app.SR830Label = uilabel(app.GridLayout15);
            app.SR830Label.HorizontalAlignment = 'center';
            app.SR830Label.Layout.Row = 1;
            app.SR830Label.Layout.Column = 1;
            app.SR830Label.Text = 'SR830';

            % Create SetVdButton
            app.SetVdButton = uibutton(app.GridLayout15, 'push');
            app.SetVdButton.ButtonPushedFcn = createCallbackFcn(app, @SetVdButtonPushed, true);
            app.SetVdButton.Enable = 'off';
            app.SetVdButton.Layout.Row = 2;
            app.SetVdButton.Layout.Column = 4;
            app.SetVdButton.Text = 'Set Vd';

            % Create WaitTimesEditFieldLabel
            app.WaitTimesEditFieldLabel = uilabel(app.GridLayout15);
            app.WaitTimesEditFieldLabel.HorizontalAlignment = 'center';
            app.WaitTimesEditFieldLabel.Layout.Row = 1;
            app.WaitTimesEditFieldLabel.Layout.Column = 5;
            app.WaitTimesEditFieldLabel.Text = 'Wait Time (s)';

            % Create WaitTimesEditField
            app.WaitTimesEditField = uieditfield(app.GridLayout15, 'numeric');
            app.WaitTimesEditField.Limits = [0 5];
            app.WaitTimesEditField.ValueDisplayFormat = '%5.2g';
            app.WaitTimesEditField.HorizontalAlignment = 'center';
            app.WaitTimesEditField.Layout.Row = 2;
            app.WaitTimesEditField.Layout.Column = 5;
            app.WaitTimesEditField.Value = 0.5;

            % Create VdField
            app.VdField = uieditfield(app.GridLayout15, 'numeric');
            app.VdField.Limits = [0 4];
            app.VdField.ValueDisplayFormat = '%4.2f';
            app.VdField.HorizontalAlignment = 'center';
            app.VdField.Layout.Row = 1;
            app.VdField.Layout.Column = 4;
            app.VdField.Value = 3;

            % Create VzField
            app.VzField = uieditfield(app.GridLayout15, 'numeric');
            app.VzField.Limits = [-250 250];
            app.VzField.ValueDisplayFormat = '%4.2f';
            app.VzField.HorizontalAlignment = 'center';
            app.VzField.Layout.Row = [4 5];
            app.VzField.Layout.Column = 4;

            % Create SetVzButton
            app.SetVzButton = uibutton(app.GridLayout15, 'push');
            app.SetVzButton.ButtonPushedFcn = createCallbackFcn(app, @SetVzButtonPushed, true);
            app.SetVzButton.Enable = 'off';
            app.SetVzButton.Layout.Row = 3;
            app.SetVzButton.Layout.Column = 5;
            app.SetVzButton.Text = 'Set Vz';

            % Create FlashVzButton
            app.FlashVzButton = uibutton(app.GridLayout15, 'push');
            app.FlashVzButton.ButtonPushedFcn = createCallbackFcn(app, @FlashVzButtonPushed, true);
            app.FlashVzButton.Enable = 'off';
            app.FlashVzButton.Layout.Row = 5;
            app.FlashVzButton.Layout.Column = 5;
            app.FlashVzButton.Text = 'Flash Vz';

            % Create ZeroVzButton
            app.ZeroVzButton = uibutton(app.GridLayout15, 'push');
            app.ZeroVzButton.ButtonPushedFcn = createCallbackFcn(app, @ZeroVzButtonPushed, true);
            app.ZeroVzButton.Enable = 'off';
            app.ZeroVzButton.Layout.Row = 4;
            app.ZeroVzButton.Layout.Column = 5;
            app.ZeroVzButton.Text = 'Zero Vz';

            % Create SR844Label
            app.SR844Label = uilabel(app.GridLayout15);
            app.SR844Label.HorizontalAlignment = 'center';
            app.SR844Label.Layout.Row = 1;
            app.SR844Label.Layout.Column = 2;
            app.SR844Label.Text = 'SR844';

            % Create SR844DropDown
            app.SR844DropDown = uidropdown(app.GridLayout15);
            app.SR844DropDown.Items = {'a', 'b'};
            app.SR844DropDown.Layout.Row = 2;
            app.SR844DropDown.Layout.Column = 2;
            app.SR844DropDown.Value = 'a';

            % Create SR844ConnectButton
            app.SR844ConnectButton = uibutton(app.GridLayout15, 'push');
            app.SR844ConnectButton.ButtonPushedFcn = createCallbackFcn(app, @SR844ConnectButtonPushed, true);
            app.SR844ConnectButton.Layout.Row = 3;
            app.SR844ConnectButton.Layout.Column = 2;
            app.SR844ConnectButton.Text = 'Connect';

            % Create StatusPanel
            app.StatusPanel = uipanel(app.GridLayout);
            app.StatusPanel.TitlePosition = 'centertop';
            app.StatusPanel.Title = 'Status';
            app.StatusPanel.Layout.Row = 3;
            app.StatusPanel.Layout.Column = 3;

            % Create GridLayout16
            app.GridLayout16 = uigridlayout(app.StatusPanel);
            app.GridLayout16.ColumnWidth = {'1x', '1x', '1x', '1x'};
            app.GridLayout16.RowHeight = {'1x', 15};

            % Create V2GaugeLabel
            app.V2GaugeLabel = uilabel(app.GridLayout16);
            app.V2GaugeLabel.HorizontalAlignment = 'center';
            app.V2GaugeLabel.Layout.Row = 2;
            app.V2GaugeLabel.Layout.Column = 2;
            app.V2GaugeLabel.Text = 'V2';

            % Create V2Gauge
            app.V2Gauge = uigauge(app.GridLayout16, 'semicircular');
            app.V2Gauge.Limits = [-250 250];
            app.V2Gauge.Orientation = 'east';
            app.V2Gauge.ScaleDirection = 'counterclockwise';
            app.V2Gauge.Layout.Row = 1;
            app.V2Gauge.Layout.Column = 2;

            % Create V3GaugeLabel
            app.V3GaugeLabel = uilabel(app.GridLayout16);
            app.V3GaugeLabel.HorizontalAlignment = 'center';
            app.V3GaugeLabel.Layout.Row = 2;
            app.V3GaugeLabel.Layout.Column = 3;
            app.V3GaugeLabel.Text = 'V3';

            % Create V3Gauge
            app.V3Gauge = uigauge(app.GridLayout16, 'semicircular');
            app.V3Gauge.Limits = [-250 250];
            app.V3Gauge.Orientation = 'east';
            app.V3Gauge.ScaleDirection = 'counterclockwise';
            app.V3Gauge.Layout.Row = 1;
            app.V3Gauge.Layout.Column = 3;

            % Create V1GaugeLabel
            app.V1GaugeLabel = uilabel(app.GridLayout16);
            app.V1GaugeLabel.HorizontalAlignment = 'center';
            app.V1GaugeLabel.Layout.Row = 2;
            app.V1GaugeLabel.Layout.Column = 1;
            app.V1GaugeLabel.Text = 'V1';

            % Create V1Gauge
            app.V1Gauge = uigauge(app.GridLayout16, 'semicircular');
            app.V1Gauge.Limits = [-250 250];
            app.V1Gauge.Orientation = 'east';
            app.V1Gauge.ScaleDirection = 'counterclockwise';
            app.V1Gauge.Layout.Row = 1;
            app.V1Gauge.Layout.Column = 1;

            % Create VzGaugeLabel
            app.VzGaugeLabel = uilabel(app.GridLayout16);
            app.VzGaugeLabel.HorizontalAlignment = 'center';
            app.VzGaugeLabel.Layout.Row = 2;
            app.VzGaugeLabel.Layout.Column = 4;
            app.VzGaugeLabel.Text = 'Vz';

            % Create VzGauge
            app.VzGauge = uigauge(app.GridLayout16, 'semicircular');
            app.VzGauge.Limits = [-250 250];
            app.VzGauge.Orientation = 'east';
            app.VzGauge.ScaleDirection = 'counterclockwise';
            app.VzGauge.Layout.Row = 1;
            app.VzGauge.Layout.Column = 4;

            % Create GridLayout6
            app.GridLayout6 = uigridlayout(app.GridLayout);
            app.GridLayout6.Layout.Row = 1;
            app.GridLayout6.Layout.Column = 2;

            % Create SnapshotButton
            app.SnapshotButton = uibutton(app.GridLayout6, 'push');
            app.SnapshotButton.ButtonPushedFcn = createCallbackFcn(app, @SnapshotButtonPushed, true);
            app.SnapshotButton.Enable = 'off';
            app.SnapshotButton.Layout.Row = 1;
            app.SnapshotButton.Layout.Column = 2;
            app.SnapshotButton.Text = 'Snapshot';

            % Create ContinuousButton
            app.ContinuousButton = uibutton(app.GridLayout6, 'push');
            app.ContinuousButton.ButtonPushedFcn = createCallbackFcn(app, @ContinuousButtonPushed, true);
            app.ContinuousButton.Enable = 'off';
            app.ContinuousButton.Layout.Row = 2;
            app.ContinuousButton.Layout.Column = 1;
            app.ContinuousButton.Text = 'Continuous';

            % Create StartPreviewButton
            app.StartPreviewButton = uibutton(app.GridLayout6, 'push');
            app.StartPreviewButton.ButtonPushedFcn = createCallbackFcn(app, @StartPreviewButtonPushed, true);
            app.StartPreviewButton.WordWrap = 'on';
            app.StartPreviewButton.Layout.Row = 1;
            app.StartPreviewButton.Layout.Column = 1;
            app.StartPreviewButton.Text = 'Start Preview';

            % Create RecordMovieButton
            app.RecordMovieButton = uibutton(app.GridLayout6, 'push');
            app.RecordMovieButton.ButtonPushedFcn = createCallbackFcn(app, @RecordMovieButtonPushed, true);
            app.RecordMovieButton.WordWrap = 'on';
            app.RecordMovieButton.Enable = 'off';
            app.RecordMovieButton.Layout.Row = 2;
            app.RecordMovieButton.Layout.Column = 2;
            app.RecordMovieButton.Text = 'Record Movie';

            % Create CapacitiveMotionDetectionPanel
            app.CapacitiveMotionDetectionPanel = uipanel(app.GridLayout);
            app.CapacitiveMotionDetectionPanel.TitlePosition = 'centertop';
            app.CapacitiveMotionDetectionPanel.Title = 'Capacitive Motion Detection';
            app.CapacitiveMotionDetectionPanel.Layout.Row = 3;
            app.CapacitiveMotionDetectionPanel.Layout.Column = 4;

            % Create GridLayout17
            app.GridLayout17 = uigridlayout(app.CapacitiveMotionDetectionPanel);
            app.GridLayout17.ColumnWidth = {'1x', '1x', '1x', '1x'};
            app.GridLayout17.RowHeight = {'1x', '1x', '1x', '1x'};

            % Create AxesCapFitting
            app.AxesCapFitting = uiaxes(app.GridLayout17);
            title(app.AxesCapFitting, 'Fitting')
            xlabel(app.AxesCapFitting, 'Steps')
            ylabel(app.AxesCapFitting, 'Signal')
            zlabel(app.AxesCapFitting, 'Z')
            app.AxesCapFitting.Box = 'on';
            app.AxesCapFitting.Layout.Row = [1 3];
            app.AxesCapFitting.Layout.Column = [3 4];

            % Create AxesCapReadout
            app.AxesCapReadout = uiaxes(app.GridLayout17);
            title(app.AxesCapReadout, {'Capacitance Readout'; ''})
            xlabel(app.AxesCapReadout, 'Time')
            ylabel(app.AxesCapReadout, 'Signal')
            zlabel(app.AxesCapReadout, 'Z')
            app.AxesCapReadout.Box = 'on';
            app.AxesCapReadout.Layout.Row = [1 3];
            app.AxesCapReadout.Layout.Column = [1 2];

            % Create ClearCapButton
            app.ClearCapButton = uibutton(app.GridLayout17, 'push');
            app.ClearCapButton.ButtonPushedFcn = createCallbackFcn(app, @ClearCapButtonPushed, true);
            app.ClearCapButton.Layout.Row = 4;
            app.ClearCapButton.Layout.Column = 3;
            app.ClearCapButton.Text = {'Clear'; 'Cap Readout'};

            % Create CapValue
            app.CapValue = uilabel(app.GridLayout17);
            app.CapValue.HorizontalAlignment = 'center';
            app.CapValue.FontName = 'Consolas';
            app.CapValue.FontWeight = 'bold';
            app.CapValue.FontColor = [0.9608 0.4667 0.1608];
            app.CapValue.Layout.Row = 4;
            app.CapValue.Layout.Column = 1;
            app.CapValue.Text = '0.0 uV';

            % Create StartStopCapButton
            app.StartStopCapButton = uibutton(app.GridLayout17, 'push');
            app.StartStopCapButton.ButtonPushedFcn = createCallbackFcn(app, @StartStopCapButtonPushed, true);
            app.StartStopCapButton.BackgroundColor = [0 0 0];
            app.StartStopCapButton.FontColor = [1 1 1];
            app.StartStopCapButton.Enable = 'off';
            app.StartStopCapButton.Layout.Row = 4;
            app.StartStopCapButton.Layout.Column = 2;
            app.StartStopCapButton.Text = {'Start'; 'Update'};

            % Create PlotCapLossButton
            app.PlotCapLossButton = uibutton(app.GridLayout17, 'push');
            app.PlotCapLossButton.ButtonPushedFcn = createCallbackFcn(app, @PlotCapLossButtonPushed, true);
            app.PlotCapLossButton.Layout.Row = 4;
            app.PlotCapLossButton.Layout.Column = 4;
            app.PlotCapLossButton.Text = 'Trace Feedback';

            % Show the figure after all components are created
            app.MEMSStepperControlUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = MEMStepper_DualLockIn

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.MEMSStepperControlUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.MEMSStepperControlUIFigure)
        end
    end
end