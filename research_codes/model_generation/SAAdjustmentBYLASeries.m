function varargout = SAAdjustmentBYLASeries_DBG(varargin) 
% SAADJUSTMENTBYLASERIES MATLAB code for SAAdjustmentBYLASeries.fig
%      SAADJUSTMENTBYLASERIES, by itself, creates a new SAADJUSTMENTBYLASERIES or raises the existing
%      singleton*.
%
%      H = SAADJUSTMENTBYLASERIES returns the handle to a new SAADJUSTMENTBYLASERIES or the handle to
%      the existing singleton*.
%
%      SAADJUSTMENTBYLASERIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SAADJUSTMENTBYLASERIES.M with the given input arguments.
% 
%      SAADJUSTMENTBYLASERIES('Property','Value',...) creates a new SAADJUSTMENTBYLASERIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SAAdjustmentBYLASeries_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SAAdjustmentBYLASeries_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SAAdjustmentBYLASeries

% Last Modified by GUIDE v2.5 30-Jul-2019 13:30:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SAAdjustmentBYLASeries_OpeningFcn, ...
                   'gui_OutputFcn',  @SAAdjustmentBYLASeries_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SAAdjustmentBYLASeries is made visible.
function SAAdjustmentBYLASeries_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SAAdjustmentBYLASeries (see VARARGIN)

% Choose default command line output for SAAdjustmentBYLASeries
clc;
handles.output = hObject;
LVWM_config;
% handles.dicomDir = dicomDir;
handles.resultDir = resultDir;
handles.workingDir = workingDir;

%%%load results
cd(resultDir);
load imDesired;
cd(workingDir);

%SXSliceSorted = imDesired.SXSlice;
%LVOTSliceSorted = imDesired.LVOTSlice;
% use the saved SXSliceSorted and LVOTSliceSorted

totalSXSliceLocation = size(SXSliceSorted, 2);
totalLVOTSliceLocation = size(LVOTSliceSorted,2);
%TimeEndOfSystole = imDesired.TimeEndOfSystole; 
%TimeEndOfDiastole = imDesired.TimeEndOfDiastole; 
%TimeEarlyOfDiastole =  imDesired.TimeEarlyOfDiastole;


% TimeInstanceSelected = patientConfigs.timeInstanceSelected;
% TimeInstanceSelected = TimeEndOfDiastole;
list_phase = {'early_diastole', 'end_diastole', 'end_systole'};
[idx, ~] = listdlg('ListString', list_phase);
%msgStr = sprintf('TimeInstanceSelected:  %s, please check', list_phase{idx});
%msgbox(msgStr, 'warning');

%% assume all short-axis image will have same cardiac phases
if strcmp(list_phase{idx}, 'early_diastole')
     timeInstanceSelected =  SXSliceSorted(1,1).TimeEarlyOfDiastole;
     timeInstanceSelected_LVOT = LVOTSliceSorted(1,1).TimeEarlyOfDiastole;
     timeInstanceSelected_4CH = LVOTSliceSorted(1,1).TimeEarlyOfDiastole;
     timeInstanceSelected_1CH = LVOTSliceSorted(1,1).TimeEarlyOfDiastole;
elseif  strcmp(list_phase{idx}, 'end_diastole')
     timeInstanceSelected = SXSliceSorted(1,1).TimeEndOfDiastole;
     timeInstanceSelected_LVOT = LVOTSliceSorted(1,1).TimeEndOfDiastole;
     timeInstanceSelected_4CH =  LVOTSliceSorted(1,1).TimeEndOfDiastole;
     timeInstanceSelected_1CH =  LVOTSliceSorted(1,1).TimeEndOfDiastole;
elseif strcmp(list_phase{idx}, 'end_systole')
     timeInstanceSelected = SXSliceSorted(1,1).TimeEndOfSystole;
     timeInstanceSelected_LVOT = LVOTSliceSorted(1,1).TimeEndOfSystole;
     timeInstanceSelected_4CH =  LVOTSliceSorted(1,1).TimeEndOfSystole;
     timeInstanceSelected_1CH =  LVOTSliceSorted(1,1).TimeEndOfSystole;
else
     disp('could not determine the cardiac phase, quite')
     return;
end


%%also need to figure out phase dir
cd(resultDir);
if ~exist(list_phase{idx},'dir')
    mkdir(list_phase{idx});
    cd(list_phase{idx});
    phase_resultDir = pwd();    
end
cd(list_phase{idx});
phase_resultDir = pwd();
cd(workingDir);
handles.phase_resultDir = phase_resultDir;

cd(phase_resultDir);
load DataSegSA;
load DataSegLA;
cd(workingDir);

totalTimeInstance = size(SXSliceSorted(1,1).SXSlice,1);


%%%now figure out what others 
handles.totalSXSliceLocation = totalSXSliceLocation;
handles.usuableSXSlice = totalSXSliceLocation;
handles.apexSliceIndex = totalSXSliceLocation-2;
handles.SASlicePositionApex = totalSXSliceLocation;
handles.totalLVOTSliceLocation = totalLVOTSliceLocation;

handles.totalTimeInstance = totalTimeInstance;
handles.timeInstanceSelected = timeInstanceSelected;
handles.timeInstanceSelected_LVOT = timeInstanceSelected_LVOT;
handles.timeInstanceSelected_4CH = timeInstanceSelected_4CH;
handles.timeInstanceSelected_1CH = timeInstanceSelected_1CH;

handles.sampleN = sampleN;
handles.SASliceDistance = 10; %%%will be updated later
handles.LASliceSelected = totalLVOTSliceLocation;

handles.SXSliceSorted = SXSliceSorted;
handles.DataSegSA = DataSegSA;
handles.DataSegSAOrig = DataSegSA;
handles.DataSegLA = DataSegLA;
handles.DataSegLAOrig = DataSegLA;
handles.LVOTSliceSorted = LVOTSliceSorted;
handles.sliceNoToBeCorrected = [];
handles.figureSA = [];
handles.figureLA3D = [];
handles.imdataSA = [];
handles.imDataLA3D = [];
handles.infoSA = [];
handles.infoLA = [];
handles.crossPx = [];
handles.crossPy = [];
handles.crossVec = [];




%%% figure out the right LA data if using LVOT, 1Ch or 4CH
% LVOTSelected = 1;
prompt = {'1 for LVOT, 2 for 4CH, 3 for OneCH'};
dlg_title = 'choose logitudial view';
num_lines = 1;
def = {'1'};
LVOTSelected = inputdlg(prompt,dlg_title,num_lines,def);
LVOTSelected = str2num(LVOTSelected{1});

LVOT_index   = 1;
FourCH_index = 2;
OneCH_index  = 3;

if LVOTSelected == 1
    handles.imDataLA  =  LVOTSliceSorted(1,LVOT_index).LXSlice(1, timeInstanceSelected_LVOT).imData;
    handles.imInfo1LA =  LVOTSliceSorted(1,LVOT_index).LXSlice(1, timeInstanceSelected_LVOT).imInfo;
elseif LVOTSelected == 2
    handles.imDataLA  =  LVOTSliceSorted(1,FourCH_index).LXSlice(1, timeInstanceSelected_4CH).imData;
    handles.imInfo1LA = LVOTSliceSorted(1,FourCH_index).LXSlice(1, timeInstanceSelected_4CH).imInfo;
elseif LVOTSelected == 3 
    handles.imDataLA  =  LVOTSliceSorted(1,OneCH_index).LXSlice(1, timeInstanceSelected_1CH).imData;
    handles.imInfo1LA = LVOTSliceSorted(1,OneCH_index).LXSlice(1, timeInstanceSelected_1CH).imInfo;
end
    

handles.LVOT_index   = LVOT_index;
handles.FourCH_index = FourCH_index;
handles. OneCH_index  = OneCH_index;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SAAdjustmentBYLASeries wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SAAdjustmentBYLASeries_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in move_plus.
function move_plus_Callback(hObject, eventdata, handles)
% hObject    handle to move_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sliceNoToBeCorrected = handles.sliceNoToBeCorrected;
DataSegSA = handles.DataSegSA;
DataSegSAOrig = handles.DataSegSAOrig;
endo_lvurrent = DataSegSA(sliceNoToBeCorrected).endo_lv;
endo_rvurrent = DataSegSA(sliceNoToBeCorrected).endo_rv;
endo_orig = DataSegSAOrig(sliceNoToBeCorrected).endo_lv;
endo_orig = DataSegSAOrig(sliceNoToBeCorrected).endo_rv;
epi_current = DataSegSA(sliceNoToBeCorrected).epi_c;

%%%move along the cutting plane
disVec = handles.crossVec;
endo_lvurrent(1,:) = endo_lvurrent(1,:) + disVec(1);
endo_rvurrent(1,:) = endo_rvurrent(1,:) + disVec(1);
endo_lvurrent(2,:) = endo_lvurrent(2,:) + disVec(2);
endo_rvurrent(2,:) = endo_rvurrent(2,:) + disVec(2);

epi_current(1,:) = epi_current(1,:) + disVec(1);
epi_current(2,:) = epi_current(2,:) + disVec(2);

imInfo = infoExtract(handles.infoSA);
endo_lvReal = TransformCurvesFromImToRealSpace(endo_lvurrent,imInfo);
endo_rvReal = TransformCurvesFromImToRealSpace(endo_rvurrent,imInfo);
epi_cReal = TransformCurvesFromImToRealSpace(epi_current,imInfo);
%%update the result
DataSegSA(sliceNoToBeCorrected).endo_lv = endo_lvurrent;
DataSegSA(sliceNoToBeCorrected).endo_rv = endo_rvurrent;
DataSegSA(sliceNoToBeCorrected).endo_lvReal = endo_lvReal;
DataSegSA(sliceNoToBeCorrected).endo_rvReal = endo_rvReal;
DataSegSA(sliceNoToBeCorrected).epi_c = epi_current;
DataSegSA(sliceNoToBeCorrected).epi_cReal = epi_cReal;
handles.DataSegSA = DataSegSA;

figure(handles.figureSA); 
% clf(handles.figureSA);
% hold on;
imshow(handles.imdataSA,[]); hold on;
plot(endo_lvurrent(1,:),endo_lvurrent(2,:),'r-', 'LineWidth',2);
plot(endo_rvurrent(1,:),endo_rvurrent(2,:),'r-', 'LineWidth',2);
plot(epi_current(1,:),epi_current(2,:),'b-', 'LineWidth',2);
plot(endo_orig(1,:), endo_orig(2,:),'w-');

%%%%3D figure
figure(handles.figureLA3D);
DicomeRealDisplay(handles.figureLA3D, handles.imDataLA3D); hold on;
delete(handles.curve3Dh_lv)
delete(handles.curve3Dh_rv)
curve3Dh_lv = plot3(endo_lvReal(1,:),endo_lvReal(2,:), endo_lvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
curve3Dh_rv = plot3(endo_rvReal(1,:),endo_rvReal(2,:), endo_rvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
hold off
handles.curve3Dh_lv = curve3Dh_lv;
handles.curve3Dh_rv = curve3Dh_rv;
%%%update the whole structure
guidata(hObject, handles);

% --- Executes on button press in move_minus.
function move_minus_Callback(hObject, eventdata, handles)
% hObject    handle to move_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sliceNoToBeCorrected = handles.sliceNoToBeCorrected;
DataSegSA = handles.DataSegSA;
DataSegSAOrig = handles.DataSegSA;
endo_lvurrent = DataSegSA(sliceNoToBeCorrected).endo_lv;
endo_rvurrent = DataSegSA(sliceNoToBeCorrected).endo_rv;
endo_orig = DataSegSAOrig(sliceNoToBeCorrected).endo_lv;
endo_orig = DataSegSAOrig(sliceNoToBeCorrected).endo_rv;
epi_current = DataSegSA(sliceNoToBeCorrected).epi_c;

%%%move along the cutting plane
disVec = handles.crossVec;
endo_lvurrent(1,:) = endo_lvurrent(1,:) - disVec(1);
endo_rvurrent(1,:) = endo_rvurrent(1,:) - disVec(1);
endo_lvurrent(2,:) = endo_lvurrent(2,:) - disVec(2);
endo_rvurrent(2,:) = endo_rvurrent(2,:) - disVec(2);

epi_current(1,:) = epi_current(1,:) - disVec(1);
epi_current(2,:) = epi_current(2,:) - disVec(2);

imInfo = infoExtract(handles.infoSA);
endo_lvReal = TransformCurvesFromImToRealSpace(endo_lvurrent,imInfo);
endo_rvReal = TransformCurvesFromImToRealSpace(endo_rvurrent,imInfo);
epi_cReal = TransformCurvesFromImToRealSpace(epi_current,imInfo);
%%update the result
DataSegSA(sliceNoToBeCorrected).endo_lv = endo_lvurrent;
DataSegSA(sliceNoToBeCorrected).endo_rv = endo_rvurrent;
DataSegSA(sliceNoToBeCorrected).endo_lvReal = endo_lvReal;
DataSegSA(sliceNoToBeCorrected).endo_rvReal = endo_rvReal;
DataSegSA(sliceNoToBeCorrected).epi_c = epi_current;
DataSegSA(sliceNoToBeCorrected).epi_cReal = epi_cReal;
handles.DataSegSA = DataSegSA;

figure(handles.figureSA); 
% clf(handles.figureSA);
% hold on;
imshow(handles.imdataSA,[]); hold on;
plot(endo_lvurrent(1,:),endo_lvurrent(2,:),'r-', 'LineWidth',2);
plot(endo_rvurrent(1,:),endo_rvurrent(2,:),'r-', 'LineWidth',2);
plot(epi_current(1,:),epi_current(2,:),'b-', 'LineWidth',2);
plot(endo_orig(1,:), endo_orig(2,:),'w-');

%%%%3D figure
figure(handles.figureLA3D); hold off;
DicomeRealDisplay(handles.figureLA3D, handles.imDataLA3D); hold on;
delete(handles.curve3Dh_lv)
delete(handles.curve3Dh_rv)
curve3Dh_lv = plot3(endo_lvReal(1,:),endo_lvReal(2,:), endo_lvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
curve3Dh_rv = plot3(endo_rvReal(1,:),endo_rvReal(2,:), endo_rvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
handles.curve3Dh_lv = curve3Dh_lv;
handles.curve3Dh_rv = curve3Dh_rv;


%%%update the whole structure
guidata(hObject, handles);



% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% move inward
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sliceNoToBeCorrected = handles.sliceNoToBeCorrected;
DataSegSA = handles.DataSegSA;
DataSegSAOrig = handles.DataSegSA;
endo_lvurrent = DataSegSA(sliceNoToBeCorrected).endo_lv;
endo_rvurrent = DataSegSA(sliceNoToBeCorrected).endo_rv;
endo_orig = DataSegSAOrig(sliceNoToBeCorrected).endo_lv;
endo_orig = DataSegSAOrig(sliceNoToBeCorrected).endo_rv;
epi_current = DataSegSA(sliceNoToBeCorrected).epi_c;

%%%move along the cutting plane
disVec = handles.crossVec;
disVec = [-disVec(2), disVec(1)]; %% perpendicular
endo_lvurrent(1,:) = endo_lvurrent(1,:) - disVec(1);
endo_rvurrent(1,:) = endo_rvurrent(1,:) - disVec(1);
endo_lvurrent(2,:) = endo_lvurrent(2,:) - disVec(2);
endo_rvurrent(2,:) = endo_rvurrent(2,:) - disVec(2);

epi_current(1,:) = epi_current(1,:) - disVec(1);
epi_current(2,:) = epi_current(2,:) - disVec(2);

imInfo = infoExtract(handles.infoSA);
endo_lvReal = TransformCurvesFromImToRealSpace(endo_lvurrent,imInfo);
endo_rvReal = TransformCurvesFromImToRealSpace(endo_rvurrent,imInfo);
epi_cReal = TransformCurvesFromImToRealSpace(epi_current,imInfo);
%%update the result
DataSegSA(sliceNoToBeCorrected).endo_lv = endo_lvurrent;
DataSegSA(sliceNoToBeCorrected).endo_rv = endo_rvurrent;
DataSegSA(sliceNoToBeCorrected).endo_lvReal = endo_lvReal;
DataSegSA(sliceNoToBeCorrected).endo_rvReal = endo_rvReal;
DataSegSA(sliceNoToBeCorrected).epi_c = epi_current;
DataSegSA(sliceNoToBeCorrected).epi_cReal = epi_cReal;
handles.DataSegSA = DataSegSA;

figure(handles.figureSA); 
% clf(handles.figureSA);
% hold on;
imshow(handles.imdataSA,[]); hold on;
plot(endo_lvurrent(1,:),endo_lvurrent(2,:),'r-', 'LineWidth',2);
plot(endo_rvurrent(1,:),endo_rvurrent(2,:),'r-', 'LineWidth',2);
plot(epi_current(1,:),epi_current(2,:),'b-', 'LineWidth',2);
plot(endo_orig(1,:), endo_orig(2,:),'w-');

%%%%3D figure
figure(handles.figureLA3D); hold off;
DicomeRealDisplay(handles.figureLA3D, handles.imDataLA3D); hold on;
delete(handles.curve3Dh_lv)
delete(handles.curve3Dh_rv)
curve3Dh_lv = plot3(endo_lvReal(1,:),endo_lvReal(2,:), endo_lvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
curve3Dh_rv = plot3(endo_rvReal(1,:),endo_rvReal(2,:), endo_rvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
handles.curve3Dh_lv = curve3Dh_lv;
handles.curve3Dh_rv = curve3Dh_rv;



%%%update the whole structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% move_outward
% move inward
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sliceNoToBeCorrected = handles.sliceNoToBeCorrected;
DataSegSA = handles.DataSegSA;
DataSegSAOrig = handles.DataSegSA;
endo_lvurrent = DataSegSA(sliceNoToBeCorrected).endo_lv;
endo_rvurrent = DataSegSA(sliceNoToBeCorrected).endo_rv;
endo_orig = DataSegSAOrig(sliceNoToBeCorrected).endo_lv;
endo_orig = DataSegSAOrig(sliceNoToBeCorrected).endo_rv;
epi_current = DataSegSA(sliceNoToBeCorrected).epi_c;

%%%move along the cutting plane
disVec = handles.crossVec;
disVec = -[-disVec(2), disVec(1)]; %% perpendicular
endo_lvurrent(1,:) = endo_lvurrent(1,:) - disVec(1);
endo_rvurrent(1,:) = endo_rvurrent(1,:) - disVec(1);
endo_lvurrent(2,:) = endo_lvurrent(2,:) - disVec(2);
endo_rvurrent(2,:) = endo_rvurrent(2,:) - disVec(2);

epi_current(1,:) = epi_current(1,:) - disVec(1);
epi_current(2,:) = epi_current(2,:) - disVec(2);

imInfo = infoExtract(handles.infoSA);
endo_lvReal = TransformCurvesFromImToRealSpace(endo_lvurrent,imInfo);
endo_rvReal = TransformCurvesFromImToRealSpace(endo_rvurrent,imInfo);
epi_cReal = TransformCurvesFromImToRealSpace(epi_current,imInfo);
%%update the result
DataSegSA(sliceNoToBeCorrected).endo_lv = endo_lvurrent;
DataSegSA(sliceNoToBeCorrected).endo_rv = endo_rvurrent;
DataSegSA(sliceNoToBeCorrected).endo_lvReal = endo_lvReal;
DataSegSA(sliceNoToBeCorrected).endo_rvReal = endo_rvReal;
DataSegSA(sliceNoToBeCorrected).epi_c = epi_current;
DataSegSA(sliceNoToBeCorrected).epi_cReal = epi_cReal;
handles.DataSegSA = DataSegSA;

figure(handles.figureSA); 
% clf(handles.figureSA);
% hold on;
imshow(handles.imdataSA,[]); hold on;
plot(endo_lvurrent(1,:),endo_lvurrent(2,:),'r-', 'LineWidth',2);
plot(endo_rvurrent(1,:),endo_rvurrent(2,:),'r-', 'LineWidth',2);
plot(epi_current(1,:),epi_current(2,:),'b-', 'LineWidth',2);
plot(endo_orig(1,:), endo_orig(2,:),'w-');

%%%%3D figure
figure(handles.figureLA3D); hold off;
DicomeRealDisplay(handles.figureLA3D, handles.imDataLA3D); hold on;
delete(handles.curve3Dh_lv)
delete(handles.curve3Dh_rv)
curve3Dh_lv = plot3(endo_lvReal(1,:),endo_lvReal(2,:), endo_lvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
curve3Dh_rv = plot3(endo_rvReal(1,:),endo_rvReal(2,:), endo_rvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
handles.curve3Dh_lv = curve3Dh_lv;
handles.curve3Dh_rv = curve3Dh_rv;

%%%update the whole structure
guidata(hObject, handles);

function input_slice_NO_Callback(hObject, eventdata, handles)
% hObject    handle to input_slice_NO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_slice_NO as text
%        str2double(get(hObject,'String')) returns contents of input_slice_NO as a double
handles.sliceNoToBeCorrected = str2double(get(hObject,'String'));
str = sprintf('slice No %d: is being selected',handles.sliceNoToBeCorrected);
disp(str);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function input_slice_NO_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_slice_NO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_displayImages.
function btn_displayImages_Callback(hObject, eventdata, handles)
% hObject    handle to btn_displayImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%display figure 1 and figure 2 with SA and LA views
sliceNoToBeCorrected = handles.sliceNoToBeCorrected;
timeInstanceSelected = handles.timeInstanceSelected;
SXSliceSorted = handles.SXSliceSorted;
DataSegSA = handles.DataSegSAOrig;
% dicomDir = handles.dicomDir;

imData = SXSliceSorted(1,sliceNoToBeCorrected).SXSlice(1,timeInstanceSelected).imData;
imInfo1 = SXSliceSorted(1,sliceNoToBeCorrected).SXSlice(1, timeInstanceSelected).imInfo;
imInfo = infoExtract(imInfo1);
sliceLocationStr = sprintf('%s',imInfo.SliceLocation);
endo_lvurrent = DataSegSA(sliceNoToBeCorrected).endo_lv;
endo_rvurrent = DataSegSA(sliceNoToBeCorrected).endo_rv;
epi_current = DataSegSA(sliceNoToBeCorrected).epi_c;

if isempty(handles.figureSA)
    handles.figureSA=figure();
else
    figure(handles.figureSA);
end
imshow(imData,[]);hold on;
plot(endo_lvurrent(1,:), endo_lvurrent(2,:),'r-');
plot(endo_rvurrent(1,:), endo_rvurrent(2,:),'r-');
plot(epi_current(1,:), epi_current(2,:),'r-');
title(sliceLocationStr);

%%%3D figure with LA view
imDataLA  = handles.imDataLA;
imInfo1LA = handles.imInfo1LA;
imInfoLA = infoExtract(imInfo1LA);
imDataT = MRIMapToRealWithImageAndHeadData(imDataLA, imInfoLA);

if isempty(handles.figureLA3D)
    handles.figureLA3D = figure(); hold on;
    DicomeRealDisplay(handles.figureLA3D, imDataT);
    endo_lv = DataSegSA(sliceNoToBeCorrected).endo_lvReal;
    endo_rv = DataSegSA(sliceNoToBeCorrected).endo_rvReal;
    epi_c = DataSegSA(sliceNoToBeCorrected).epi_cReal;
    %%%%plot 3D curves
    curve3Dh_lv = plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
    curve3Dh_rv = plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
    % plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
    handles.curve3Dh_lv = curve3Dh_lv;
    handles.curve3Dh_rv = curve3Dh_rv;
end

%%%find the cutting plane of the long axis view
[px,py] = findCrossLines(imData, imDataLA, imInfo, imInfoLA);
figure(handles.figureSA);hold on;
plot(px(1),py(1),'r+');
plot(px(2),py(2),'ro');
line(px, py);

handles.infoSA = imInfo1;
handles.imdataSA = imData;
handles.imDataLA3D = imDataT;

handles.crossPx = px;
handles.crossPy = py;
handles.crossVec = NormalizationVec([px(1)-px(2) py(1)-py(2)]);
guidata(hObject, handles);


% --- Executes on button press in BCTogether3DLA.
function BCTogether3DLA_Callback(hObject, eventdata, handles)
% hObject    handle to BCTogether3DLA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LVOTSliceSorted = handles.LVOTSliceSorted;
timeInstanceSelected = handles.timeInstanceSelected;
DataSegSA = handles.DataSegSA;
usuableSXSlice = handles.usuableSXSlice;

%%%3D figure with LA view
imDataLA  = handles.imDataLA;
imInfo1LA = handles.imInfo1LA;
imInfoLA = infoExtract(imInfo1LA);
imDataT = MRIMapToRealWithImageAndHeadData(imDataLA, imInfoLA);

h3D = figure(); hold on;
DicomeRealDisplay(h3D, imDataT);

for imIndex = 1 : usuableSXSlice
    endo_lv = DataSegSA(imIndex).endo_lvReal;
    endo_rv = DataSegSA(imIndex).endo_rvReal;
    epi_c = DataSegSA(imIndex).epi_cReal;
    if ~isempty(endo_lv)  
        plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
        plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
        plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
    end
end

% --- Executes on button press in SAVE.
function SAVE_Callback(hObject, eventdata, handles)
% hObject    handle to SAVE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
workingDir = handles.workingDir;
phase_resultDir = handles.phase_resultDir;
cd(phase_resultDir);
DataSegSA = handles.DataSegSA;
save DataSegSA DataSegSA;
cd(workingDir);
disp('saving updated boundaries: done');


% --- Executes on button press in load_LASA_BCs.
function load_LASA_BCs_Callback(hObject, eventdata, handles)
% hObject    handle to load_LASA_BCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DataSegSA = handles.DataSegSA;
DataSegLA = handles.DataSegLA;
usuableSXSlice = handles.usuableSXSlice;

h3D_LAAdjust = figure(); hold on;

for imIndex = 1 : usuableSXSlice
    endo_lv = DataSegSA(imIndex).endo_lvReal;
    endo_rv = DataSegSA(imIndex).endo_rvReal;
%     epi_c = DataSegSA(imIndex).epi_cReal;
    plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '--', 'Color', 'b', 'LineWidth',1);
    plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '--', 'Color', 'b', 'LineWidth',1);
%     plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
end

for imIndex = 1 : 1
    endo_lv = DataSegLA(imIndex).endo_lvReal;
    endo_rv = DataSegLA(imIndex).endo_rvReal;
    epi_c = DataSegLA(imIndex).epi_cReal;
%     plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
%     plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
    curve3Dh_lv = plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',1);
    curve3Dh_rv = plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',1);
end
figure(h3D_LAAdjust); axis equal;

%%now figure out for each LA the up-down vector and in-out-vector
SXSliceSorted = handles.SXSliceSorted;
LVOTSliceSorted = handles.LVOTSliceSorted;

totalSXSliceLocation = handles.totalSXSliceLocation;
SXSlice_mid = int8(totalSXSliceLocation/2);
SX_epi_cReal = DataSegSA(SXSlice_mid).epi_cReal;
LVOT_epi_cReal = DataSegLA(handles.LVOT_index).epi_cReal;
FourCH_epi_cReal = DataSegLA(handles.FourCH_index).epi_cReal;
OneCH_epi_cReal = DataSegLA(handles.OneCH_index).epi_cReal;
dis_vec_up_down =  normal_plane_from_BC(SX_epi_cReal);
dis_vec_in_out_LVOT = normal_plane_from_BC(LVOT_epi_cReal);
dis_vec_in_out_4CH = normal_plane_from_BC(FourCH_epi_cReal);
dis_vec_in_out_1CH = normal_plane_from_BC(OneCH_epi_cReal);

handles.h3D_LAAdjust = h3D_LAAdjust;
handles.dis_vec_up_down = dis_vec_up_down;
handles.dis_vec_in_out_LVOT = dis_vec_in_out_LVOT;
handles.dis_vec_in_out_4CH = dis_vec_in_out_4CH;
handles.dis_vec_in_out_1CH = dis_vec_in_out_1CH;
handles.curve3Dh_lv = curve3Dh_lv;
handles.curve3Dh_rv = curve3Dh_rv;

guidata(hObject, handles);


% --- Executes on button press in LA_move_up.
function LA_move_up_Callback(hObject, eventdata, handles)
% hObject    handle to LA_move_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h3D_LAAdjust = handles.h3D_LAAdjust;
curve3Dh_lv = handles.curve3Dh_lv;
curve3Dh_rv = handles.curve3Dh_rv;
figure(h3D_LAAdjust); hold on;

DataSegSA = handles.DataSegSA;
DataSegLA = handles.DataSegLA;
usuableSXSlice = handles.usuableSXSlice;

figure(h3D_LAAdjust); hold on; axis equal; %% replot
for imIndex = 1 : usuableSXSlice
    endo_lv = DataSegSA(imIndex).endo_lvReal;
    endo_rv = DataSegSA(imIndex).endo_rvReal;
    plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '--', 'Color', 'b', 'LineWidth',1);
    plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '--', 'Color', 'b', 'LineWidth',1);
end

%% need to show the LVOT bcs
if handles.radioBtn_LVOT.Value == 1
    endo_lvReal = DataSegLA(handles.LVOT_index).endo_lvReal;
    endo_rvReal = DataSegLA(handles.LVOT_index).endo_rvReal;
    epi_cReal = DataSegLA(handles.LVOT_index).epi_cReal;
elseif handles.radioBtn_4CH.Value == 1
    endo_lvReal = DataSegLA(handles.FourCH_index).endo_lvReal;
    endo_rvReal = DataSegLA(handles.FourCH_index).endo_rvReal;
    epi_cReal = DataSegLA(handles.FourCH_index).epi_cReal;
elseif handles.radioBtn_1CH.Value == 1
    endo_lvReal = DataSegLA(handles.OneCH_index).endo_lvReal;
    endo_rvReal = DataSegLA(handles.OneCH_index).endo_rvReal;
    epi_cReal = DataSegLA(handles.OneCH_index).epi_cReal;
end
%figure(h3D_LAAdjust);
%curve3Dh = plot3(endo_lvReal(1,:),endo_lvReal(2,:), endo_lvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
%curve3Dh = plot3(endo_rvReal(1,:),endo_rvReal(2,:), endo_rvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);

PixelSpacing = handles.SXSliceSorted(1).SXSlice(1).imInfo.PixelSpacing;
%%now we will move up and down
dis_vec_up_down = handles.dis_vec_up_down;
dis_to_move = PixelSpacing(1);

vec_to_move = dis_to_move*dis_vec_up_down;
endo_lvReal(1,:) = endo_lvReal(1,:) + vec_to_move(1);
endo_rvReal(1,:) = endo_rvReal(1,:) + vec_to_move(1);
endo_lvReal(2,:) = endo_lvReal(2,:) + vec_to_move(2);
endo_rvReal(2,:) = endo_rvReal(2,:) + vec_to_move(2);
endo_lvReal(3,:) = endo_lvReal(3,:) + vec_to_move(3);
endo_rvReal(3,:) = endo_rvReal(3,:) + vec_to_move(3);

epi_cReal(1,:) = epi_cReal(1,:) + vec_to_move(1);
epi_cReal(2,:) = epi_cReal(2,:) + vec_to_move(2);
epi_cReal(3,:) = epi_cReal(3,:) + vec_to_move(3);

figure(h3D_LAAdjust);
if exist('curve3Dh_lv', 'var')
    delete(curve3Dh_lv);
    delete(curve3Dh_rv);
end
 hold on; axis equal;
curve3Dh_lv = plot3(endo_lvReal(1,:),endo_lvReal(2,:), endo_lvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2); hold on;
curve3Dh_rv = plot3(endo_rvReal(1,:),endo_rvReal(2,:), endo_rvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2); hold on;

if handles.radioBtn_LVOT.Value == 1
     DataSegLA(handles.LVOT_index).endo_lvReal = endo_lvReal;
     DataSegLA(handles.LVOT_index).endo_rvReal = endo_rvReal;
     DataSegLA(handles.LVOT_index).epi_cReal = epi_cReal;
elseif handles.radioBtn_4CH.Value == 1
     DataSegLA(handles.FourCH_index).endo_lvReal = endo_lvReal;
     DataSegLA(handles.FourCH_index).endo_rvReal = endo_rvReal;
     DataSegLA(handles.FourCH_index).epi_cReal = epi_cReal;
elseif handles.radioBtn_1CH.Value == 1
     DataSegLA(handles.OneCH_index).endo_lvReal = endo_lvReal;
     DataSegLA(handles.OneCH_index).endo_rvReal = endo_rvReal;
     DataSegLA(handles.OneCH_index).epi_cReal = epi_cReal;
end

handles.curve3Dh_lv = curve3Dh_lv;
handles.curve3Dh_rv = curve3Dh_rv;
handles.DataSegLA = DataSegLA;
guidata(hObject, handles);

% --- Executes on button press in LA_move_down.
function LA_move_down_Callback(hObject, eventdata, handles)
% hObject    handle to LA_move_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h3D_LAAdjust = handles.h3D_LAAdjust;
curve3Dh_lv = handles.curve3Dh_lv;
curve3Dh_rv = handles.curve3Dh_rv;
figure(h3D_LAAdjust); hold on;

DataSegSA = handles.DataSegSA;
DataSegLA = handles.DataSegLA;
usuableSXSlice = handles.usuableSXSlice;

figure(h3D_LAAdjust); hold on; axis equal; %% replot
for imIndex = 1 : usuableSXSlice
    endo_lv = DataSegSA(imIndex).endo_lvReal;
    endo_rv = DataSegSA(imIndex).endo_rvReal;
    plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '--', 'Color', 'b', 'LineWidth',1);
    plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '--', 'Color', 'b', 'LineWidth',1);
end

%% need to show the LVOT bcs
if handles.radioBtn_LVOT.Value == 1
    endo_lvReal = DataSegLA(handles.LVOT_index).endo_lvReal;
    endo_rvReal = DataSegLA(handles.LVOT_index).endo_rvReal;
    epi_cReal = DataSegLA(handles.LVOT_index).epi_cReal;
elseif handles.radioBtn_4CH.Value == 1
    endo_lvReal = DataSegLA(handles.FourCH_index).endo_lvReal;
    endo_rvReal = DataSegLA(handles.FourCH_index).endo_rvReal;
    epi_cReal = DataSegLA(handles.FourCH_index).epi_cReal;
elseif handles.radioBtn_1CH.Value == 1
    endo_lvReal = DataSegLA(handles.OneCH_index).endo_lvReal;
    endo_rvReal = DataSegLA(handles.OneCH_index).endo_rvReal;
    epi_cReal = DataSegLA(handles.OneCH_index).epi_cReal;
end
% figure(h3D_LAAdjust);
% curve3Dh = plot3(endo_lvReal(1,:),endo_lvReal(2,:), endo_lvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
% curve3Dh = plot3(endo_rvReal(1,:),endo_rvReal(2,:), endo_rvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);

PixelSpacing = handles.SXSliceSorted(1).SXSlice(1).imInfo.PixelSpacing;
%%now we will move up and down
dis_vec_up_down = -handles.dis_vec_up_down;
dis_to_move = PixelSpacing(1);

vec_to_move = dis_to_move*dis_vec_up_down;
endo_lvReal(1,:) = endo_lvReal(1,:) + vec_to_move(1);
endo_rvReal(1,:) = endo_rvReal(1,:) + vec_to_move(1);
endo_lvReal(2,:) = endo_lvReal(2,:) + vec_to_move(2);
endo_rvReal(2,:) = endo_rvReal(2,:) + vec_to_move(2);
endo_lvReal(3,:) = endo_lvReal(3,:) + vec_to_move(3);
endo_rvReal(3,:) = endo_rvReal(3,:) + vec_to_move(3);

epi_cReal(1,:) = epi_cReal(1,:) + vec_to_move(1);
epi_cReal(2,:) = epi_cReal(2,:) + vec_to_move(2);
epi_cReal(3,:) = epi_cReal(3,:) + vec_to_move(3);

figure(h3D_LAAdjust);
if exist('curve3Dh_lv', 'var')
    delete(curve3Dh_lv);
    delete(curve3Dh_rv);
end
hold on; axis equal;
curve3Dh_lv = plot3(endo_lvReal(1,:),endo_lvReal(2,:), endo_lvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
curve3Dh_rv = plot3(endo_rvReal(1,:),endo_rvReal(2,:), endo_rvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);

if handles.radioBtn_LVOT.Value == 1
     DataSegLA(handles.LVOT_index).endo_lvReal = endo_lvReal;
     DataSegLA(handles.LVOT_index).endo_rvReal = endo_rvReal;
     DataSegLA(handles.LVOT_index).epi_cReal = epi_cReal;
elseif handles.radioBtn_4CH.Value == 1
     DataSegLA(handles.FourCH_index).endo_lvReal = endo_lvReal;
     DataSegLA(handles.FourCH_index).endo_rvReal = endo_rvReal;
     DataSegLA(handles.FourCH_index).epi_cReal = epi_cReal;
elseif handles.radioBtn_1CH.Value == 1
     DataSegLA(handles.OneCH_index).endo_lvReal = endo_lvReal;
     DataSegLA(handles.OneCH_index).endo_rvReal = endo_rvReal;
     DataSegLA(handles.OneCH_index).epi_cReal = epi_cReal;
end

handles.DataSegLA = DataSegLA;
handles.curve3Dh_lv = curve3Dh_lv;
handles.curve3Dh_rv = curve3Dh_rv;
guidata(hObject, handles);

% --- Executes on button press in LA_move_in.
function LA_move_in_Callback(hObject, eventdata, handles)
% hObject    handle to LA_move_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h3D_LAAdjust = handles.h3D_LAAdjust;
curve3Dh_lv = handles.curve3Dh_lv;
curve3Dh_rv = handles.curve3Dh_rv;
figure(h3D_LAAdjust); hold on;

DataSegSA = handles.DataSegSA;
DataSegLA = handles.DataSegLA;
usuableSXSlice = handles.usuableSXSlice;

figure(h3D_LAAdjust); hold on; axis equal; %% replot
for imIndex = 1 : usuableSXSlice
    endo_lv = DataSegSA(imIndex).endo_lvReal;
    endo_rv = DataSegSA(imIndex).endo_rvReal;
    plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '--', 'Color', 'b', 'LineWidth',1);
    plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '--', 'Color', 'b', 'LineWidth',1);
end

%% need to show the LVOT bcs
dis_vec = [0 0 0];
if handles.radioBtn_LVOT.Value == 1
    endo_lvReal = DataSegLA(handles.LVOT_index).endo_lvReal;
    endo_rvReal = DataSegLA(handles.LVOT_index).endo_rvReal;
    epi_cReal = DataSegLA(handles.LVOT_index).epi_cReal;
    dis_vec = handles.dis_vec_in_out_LVOT;
elseif handles.radioBtn_4CH.Value == 1
    endo_lvReal = DataSegLA(handles.FourCH_index).endo_lvReal;
    endo_rvReal = DataSegLA(handles.FourCH_index).endo_rvReal;
    epi_cReal = DataSegLA(handles.FourCH_index).epi_cReal;
    dis_vec = handles.dis_vec_in_out_4CH;
elseif handles.radioBtn_1CH.Value == 1
    endo_lvReal = DataSegLA(handles.OneCH_index).endo_lvReal;
    endo_rvReal = DataSegLA(handles.OneCH_index).endo_rvReal;
    epi_cReal = DataSegLA(handles.OneCH_index).epi_cReal;
    dis_vec = handles.dis_vec_in_out_1CH;
end
% figure(h3D_LAAdjust);
% curve3Dh = plot3(endo_lvReal(1,:),endo_lvReal(2,:), endo_lvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
% curve3Dh = plot3(endo_rvReal(1,:),endo_rvReal(2,:), endo_rvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);

PixelSpacing = handles.SXSliceSorted(1).SXSlice(1).imInfo.PixelSpacing;
%%now we will move up and down

dis_to_move = PixelSpacing(1);

vec_to_move = dis_to_move*dis_vec;
endo_lvReal(1,:) = endo_lvReal(1,:) + vec_to_move(1);
endo_rvReal(1,:) = endo_rvReal(1,:) + vec_to_move(1);
endo_lvReal(2,:) = endo_lvReal(2,:) + vec_to_move(2);
endo_rvReal(2,:) = endo_rvReal(2,:) + vec_to_move(2);
endo_lvReal(3,:) = endo_lvReal(3,:) + vec_to_move(3);
endo_rvReal(3,:) = endo_rvReal(3,:) + vec_to_move(3);

epi_cReal(1,:) = epi_cReal(1,:) + vec_to_move(1);
epi_cReal(2,:) = epi_cReal(2,:) + vec_to_move(2);
epi_cReal(3,:) = epi_cReal(3,:) + vec_to_move(3);

figure(h3D_LAAdjust);
if exist('curve3Dh_lv', 'var')
    delete(curve3Dh_lv);
    delete(curve3Dh_rv);
end
hold on; axis equal;
curve3Dh_lv = plot3(endo_lvReal(1,:),endo_lvReal(2,:), endo_lvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
curve3Dh_rv = plot3(endo_rvReal(1,:),endo_rvReal(2,:), endo_rvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);

if handles.radioBtn_LVOT.Value == 1
     DataSegLA(handles.LVOT_index).endo_lvReal = endo_lvReal;
     DataSegLA(handles.LVOT_index).endo_rvReal = endo_rvReal;
     DataSegLA(handles.LVOT_index).epi_cReal = epi_cReal;
elseif handles.radioBtn_4CH.Value == 1
     DataSegLA(handles.FourCH_index).endo_lvReal = endo_lvReal;
     DataSegLA(handles.FourCH_index).endo_rvReal = endo_rvReal;
     DataSegLA(handles.FourCH_index).epi_cReal = epi_cReal;
elseif handles.radioBtn_1CH.Value == 1
     DataSegLA(handles.OneCH_index).endo_lvReal = endo_lvReal;
     DataSegLA(handles.OneCH_index).endo_rvReal = endo_rvReal;
     DataSegLA(handles.OneCH_index).epi_cReal = epi_cReal;
end

handles.DataSegLA = DataSegLA;
handles.curve3Dh_lv = curve3Dh_lv;
handles.curve3Dh_rv = curve3Dh_rv;
guidata(hObject, handles);

% --- Executes on button press in LA_move_out.
function LA_move_out_Callback(hObject, eventdata, handles)
% hObject    handle to LA_move_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h3D_LAAdjust = handles.h3D_LAAdjust;
curve3Dh_lv = handles.curve3Dh_lv;
curve3Dh_rv = handles.curve3Dh_rv;
figure(h3D_LAAdjust); hold on; 

DataSegSA = handles.DataSegSA;
DataSegLA = handles.DataSegLA;
usuableSXSlice = handles.usuableSXSlice;

figure(h3D_LAAdjust); hold on; axis equal; %% replot
for imIndex = 1 : usuableSXSlice
    endo_lv = DataSegSA(imIndex).endo_lvReal;
    endo_rv = DataSegSA(imIndex).endo_rvReal;
    plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '--', 'Color', 'b', 'LineWidth',1);
    plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '--', 'Color', 'b', 'LineWidth',1);
end

%% need to show the LVOT bcs
dis_vec = [0 0 0];
if handles.radioBtn_LVOT.Value == 1
    endo_lvReal = DataSegLA(handles.LVOT_index).endo_lvReal;
    endo_rvReal = DataSegLA(handles.LVOT_index).endo_rvReal;
    epi_cReal = DataSegLA(handles.LVOT_index).epi_cReal;
    dis_vec = handles.dis_vec_in_out_LVOT;
elseif handles.radioBtn_4CH.Value == 1
    endo_lvReal = DataSegLA(handles.FourCH_index).endo_lvReal;
    endo_rvReal = DataSegLA(handles.FourCH_index).endo_rvReal;
    epi_cReal = DataSegLA(handles.FourCH_index).epi_cReal;
    dis_vec = handles.dis_vec_in_out_4CH;
elseif handles.radioBtn_1CH.Value == 1
    endo_lvReal = DataSegLA(handles.OneCH_index).endo_lvReal;
    endo_rvReal = DataSegLA(handles.OneCH_index).endo_rvReal;
    epi_cReal = DataSegLA(handles.OneCH_index).epi_cReal;
    dis_vec = handles.dis_vec_in_out_1CH;
end
% figure(h3D_LAAdjust);
% curve3Dh = plot3(endo_lvReal(1,:),endo_lvReal(2,:), endo_lvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
% curve3Dh = plot3(endo_rvReal(1,:),endo_rvReal(2,:), endo_rvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);

PixelSpacing = handles.SXSliceSorted(1).SXSlice(1).imInfo.PixelSpacing;
%%now we will move up and down

dis_to_move = PixelSpacing(1);

vec_to_move = dis_to_move*(-dis_vec);
endo_lvReal(1,:) = endo_lvReal(1,:) + vec_to_move(1);
endo_rvReal(1,:) = endo_rvReal(1,:) + vec_to_move(1);
endo_lvReal(2,:) = endo_lvReal(2,:) + vec_to_move(2);
endo_rvReal(2,:) = endo_rvReal(2,:) + vec_to_move(2);
endo_lvReal(3,:) = endo_lvReal(3,:) + vec_to_move(3);
endo_rvReal(3,:) = endo_rvReal(3,:) + vec_to_move(3);

epi_cReal(1,:) = epi_cReal(1,:) + vec_to_move(1);
epi_cReal(2,:) = epi_cReal(2,:) + vec_to_move(2);
epi_cReal(3,:) = epi_cReal(3,:) + vec_to_move(3);

figure(h3D_LAAdjust);
if exist('curve3Dh_lv', 'var')
    delete(curve3Dh_lv);
    delete(curve3Dh_rv);
end
hold on; axis equal;
curve3Dh_lv = plot3(endo_lvReal(1,:),endo_lvReal(2,:), endo_lvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
curve3Dh_rv = plot3(endo_rvReal(1,:),endo_rvReal(2,:), endo_rvReal(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);

if handles.radioBtn_LVOT.Value == 1
     DataSegLA(handles.LVOT_index).endo_lvReal = endo_lvReal;
     DataSegLA(handles.LVOT_index).endo_rvReal = endo_rvReal;
     DataSegLA(handles.LVOT_index).epi_cReal = epi_cReal;
elseif handles.radioBtn_4CH.Value == 1
     DataSegLA(handles.FourCH_index).endo_lvReal = endo_lvReal;
     DataSegLA(handles.FourCH_index).endo_rvReal = endo_rvReal;
     DataSegLA(handles.FourCH_index).epi_cReal = epi_cReal;
elseif handles.radioBtn_1CH.Value == 1
     DataSegLA(handles.OneCH_index).endo_lvReal = endo_lvReal;
     DataSegLA(handles.OneCH_index).endo_rvReal = endo_rvReal;
     DataSegLA(handles.OneCH_index).epi_cReal = epi_cReal;
end

handles.DataSegLA = DataSegLA;
handles.curve3Dh_lv = curve3Dh_lv;
handles.curve3Dh_rv = curve3Dh_rv;
guidata(hObject, handles);
