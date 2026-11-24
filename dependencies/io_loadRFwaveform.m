% io_loadRFwaveform.m
% Jamie Near, McGill University 2014.
%
% USAGE:
% [RF_struct]=io_loadRFwaveform(filename,type,f0);
%
% DESCRIPTION:
% Initialize an RF pulse structure to contain an
% RF Pulse waveform as well as its accompanying header
% information.  This function finds the time-bandwidth
% product (tbw) and the time-w1 product (tw1) of the pulse,
% and stores this information in the header fields of the
% output RF structure. For phase-modulated pulses, the inversion plot
% suggests an amplitude and lets the user select a different point.
%
% INPUTS:
% filename  = filename of RF pulse waveform text file.  Can be in Siemens
%             format (.pta), Varian/Agilent format (.RF), Bruker format
%             (.inv, .ref or .exc) or a plain text file (.txt, with two,
%             three or four columns (amplitude, phase, time-step
%             (optional),and gradient strength (optional).  Filename can
%             also be the name of a two, three or four column matlab vector
%             specifing the phase, amplitude, time vectors (optionally) and
%             gradient (optionally) of an RF pulse waveform.  Note about
%             Units: RF amplitude waveform [arbitrary units], B1 intensity
%             [kHz],  phase waveform [degrees], time step vector [arbitrary
%             units],  gradient waveform [G/cm].
% type      = Excitation ('exc'), Refocusing ('ref') or Inversion ('inv').
%             Alternatively, 'type' can specify the exact flip angle [in
%             degrees] of the pulse.  This can be useful if you want to use
%             a pulse that does not have an exact flip angle of 90 or 180
%             degrees (Thanks to Eric Plitman for help with this feature).
% f0        = centre frequency of the rf pulse [Hz].  Optional. Default=0.
%
% OUTPUTS:
% RF_struct = RF pulse waveform in FID-A rf pulse structure format.

function RF_struct=io_loadRFwaveform(filename,type,f0)

if nargin<3
    f0=0;
end

if isnumeric(type)
    %If 'type' was specified to be exactly 90 degrees, set the type to 'exc';
    if type==90
        type='exc';
    end

    %If 'type' was specified to be exactly 180 degrees, set the type to 'inv'
    %(functionally there is no real difference between 'ref' and 'inv', so I
    %just chose 'Inv';
    if type==180
        type='inv';
    end
end

%Now read in the waveform:
if ischar(filename)
    if exist(filename,'file')
        [~,~,ext]=fileparts(filename);
        switch ext
            case '.pta'
                disp('Siemens format .pta RF pulse file detected!! Loading waveform now.');
                rf=io_readpta(filename);
            case '.RF'
                disp('Varian/Agilent format .RF RF pulse file detected!! Loading waveform now.');
                rf=io_readRF(filename);
            case '.inv'
                disp('Bruker format .inv RF pulse file detected!! Loading waveform now.');
                rf=io_readRFBruk(filename);
            case '.rfc'
                disp('Bruker format .pta RF pulse file detected!! Loading waveform now.');
                rf=io_readRFBruk(filename);
            case '.exc'
                disp('Bruker format .exc RF pulse file detected!! Loading waveform now.');
                rf=io_readRFBruk(filename);
            case '.txt'
                disp('Basic .txt format RF pulse file detected!! Loading waveform now.');
                rf=io_readRFtxt(filename);
            otherwise
                error('ERROR:  RF Pulse file not recognized.  Aborting!');
        end
    else
        error('ERROR:  File not found!  Aborting!');
    end
else
    if ismatrix(filename) && (size(filename,2) >= 2) && (size(filename,2) <= 4)
        disp('Input is an RF waveform already in the matlab workspace.  Loading waveform now.');
        if size(filename,2)==2
            rf=filename;
            rf(:,3)=ones(length(filename(:,1)),1);
        elseif size(filename,2)==3 || size(filename,2)==4
            rf=filename;
        end
    end
end


%Find out if the pulse is phase modulated.  If it is not, then we can
%determine the time-w1 product of the pulse quite simply.  If it is phase
%modulated (adiabatic, etc) then the determination of the time-w1 product
%will need to me more interactive.
a=(round(rf(:,1))==180)|(round(rf(:,1))==0);

if sum(a)<length(rf(:,1))
    isPhsMod=true;
else
    isPhsMod=false;
end

%If there are any phase discontinuities in the phase function that are
%equal to a 360 degree jump we can remove these.  This will make it easier
%for rf_resample to do it's job later on:
jumps=diff(rf(:,1));
jumpsAbs=(abs(jumps)>355 & abs(jumps)<365);  %Assume jumps within this range are exactly = 360 degrees.
jumpIndex=find(jumpsAbs);
for n=1:length(jumpIndex)
    rf(jumpIndex(n)+1:end,1)=rf(jumpIndex(n)+1:end,1)-(360*(jumps(jumpIndex(n))/abs(jumps(jumpIndex(n)))));
end

%scale amplitude function so that maximum value is 1:
rf(:,2)=rf(:,2)./max(rf(:,2));


Tp=0.005;  %assume a 5 ms rf pulse;

if ~isPhsMod
    %The pulse is not phase modulated, so we can calculate the w1max:
    %find the B1 max of the pulse in [kHz]:
    if ischar(type)
        switch type
            case 'exc'
                flipCyc=0.25; %90 degrees is 0.25 cycles;
            case 'ref'
                flipCyc=0.5; %180 degress is 0.5 cycles;
            case 'inv'
                flipCyc=0.5; %180 degrees is 0.5 cycles;
        end
    elseif isnumeric(type) %assume that a flip angle (in degrees) was given
        flipCyc=type/360;
    end
    intRF=sum(rf(:,2).*((-2*(rf(:,1)>179))+1))/length(rf(:,2));
    if intRF~=0
        w1max=flipCyc/(intRF*Tp); %w1max is in [Hz]
    else
        w1max=0;
    end
else
    %Suggest the first local minimum of Mz.
    [mv,sc]=bes(rf,Tp*1000,'b',f0/1000,0,5,40000);
    first_min=find(diff(mv(3,:))>0,1,'first');
    if isempty(first_min)
        [~,selected_ind]=min(mv(3,:));
    else
        selected_ind=first_min+1;
    end
    h=figure(101);
    set(h,'Visible','on','Name','Choose RF inversion point','NumberTitle','off');
    clf(h);
    plot(sc,mv(3,:));
    hold on;
    point_marker=plot(sc(selected_ind),mv(3,selected_ind),'ro', ...
        'MarkerFaceColor','r','MarkerSize',8);
    hold off;
    xlabel('\omega_{1} (kHz)');
    ylabel('{\it{M}}_{z}');
    title(sprintf('RF inversion calibration step (for a %.2f ms pulse)',Tp*1000));
    drawnow;

    choice=questdlg(sprintf(['Suggested point: %.4f kHz.\n' ...
        'Do you want to change the point where full inversion is achieved?'], ...
        sc(selected_ind)), ...
        'Choose RF inversion point','Yes','No','No');
    if strcmp(choice,'Yes')
        while true
            figure(h);
            title(sprintf('RF inversion calibration step (for a %.2f ms pulse): click a point',Tp*1000));
            [clicked_w1,~,button]=ginput(1);
            if isempty(clicked_w1) || button~=1
                error('RF inversion point selection was cancelled.');
            end
            [~,candidate_ind]=min(abs(sc-clicked_w1));
            set(point_marker,'XData',sc(candidate_ind),'YData',mv(3,candidate_ind));
            title(sprintf('RF inversion calibration step (for a %.2f ms pulse)',Tp*1000));
            drawnow;

            choice=questdlg(sprintf('Selected point: %.4f kHz.\nUse this point?', ...
                sc(candidate_ind)), ...
                'Confirm RF inversion point','Yes','No','Yes');
            if strcmp(choice,'Yes')
                selected_ind=candidate_ind;
                break;
            elseif ~strcmp(choice,'No')
                error('RF inversion point selection was cancelled.');
            end
        end
    elseif ~strcmp(choice,'No')
        error('RF inversion point selection was cancelled.');
    end
    fprintf('RF inversion point: %.4f kHz.\n',sc(selected_ind));
    w1max=sc(selected_ind)*1000; %convert kHz to Hz
end
tw1=Tp*w1max;
%now it's time to find out the time-bandwidth product:
%First make a high resolution plot the pulse profile over a wide bandwidth:
[mv,sc]=bes(rf,Tp*1000,'f',w1max/1000,-5+f0/1000,5+f0/1000,100000);
if ischar(type)
    switch type
        case 'exc'
            index=find(mv(3,:)<0.5);
            bw=sc(index(end))-sc(index(1));
            %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
        case 'ref'
            index=find(mv(3,:)<0);
            bw=sc(index(end))-sc(index(1));
            %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
        case 'inv'
            index=find(mv(3,:)<0);
            bw=sc(index(end))-sc(index(1));
            %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
    end
elseif isnumeric(type)
    mz=cos(type); %Find out the Mz value immediately following the pulse.
    thr=(1+mz)/2;  %Find out the Mz value mid-way between 1 and the mz (half-max):
    index=find(mv(3,:)<thr);  %Find the indices of the corresponding "full width"
    bw=sc(index(end))-sc(index(1));  %Now find the bandwidth at that point ("Full width at half max").
end


%Now make a very high resolution plot the pulse profile over a narrower bandwidth:
[mv,sc]=bes(rf,Tp*1000,'f',w1max/1000,-bw+f0/1000,bw+f0/1000,100000);
if ischar(type)
    switch type
        case 'exc'
            index=find(mv(3,:)<0.5);
            bw=sc(index(end))-sc(index(1));
            %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
        case 'ref'
            index=find(mv(3,:)<0);
            bw=sc(index(end))-sc(index(1));
            %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
        case 'inv'
            index=find(mv(3,:)<0);
            bw=sc(index(end))-sc(index(1));
            %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
    end
elseif isnumeric(type)
    mz=cos(type); %Find out the Mz value immediately following the pulse.
    thr=(1+mz)/2;  %Find out the Mz value mid-way between 1 and the mz (half-max):
    index=find(mv(3,:)<thr);  %Find the indices of the corresponding "full width"
    bw=sc(index(end))-sc(index(1));  %Now find the bandwidth at that point ("Full width at half max").
end

%Now store the output structure;
RF_struct.waveform=rf;
RF_struct.type=type;
RF_struct.f0=f0;
RF_struct.w1max=w1max;
if size(rf,2)>3 && any(rf(:,4))
    RF_struct.tbw='N/A - gradient modulated pulse';
    RF_struct.isGM=true;
    RF_struct.tthk=bw*Tp; %This is the time x sliceThickness product for
                          %gradient modulated pulses.  It is in units [cm.s]
else
    RF_struct.tbw=bw*Tp*1000;
    RF_struct.isGM=false;
    RF_struct.tthk='N/A - frequency selective pulse';
end
RF_struct.tw1=tw1;
