function [PF,A,S,ORT,Nbits]=eoe_lmd(varargin)
% EOE-LMD
% The Local Mean Decompostion Using Empirical Optimal Envelope
% Syntax:
% [PF,A,S]=eoe_lmd(x);
% [PF,A,S,ORT,Nbits] = eoe_lmd(X,...,'Option_name',Option_value,...)
% ----------------------------------------------------------
% INPUTS: 
% x: Signal, row or column vecter will be ok
% OUTPUTS:
% (1)PF: Product Functions, each row of the Matrix PF is a product function
% (2)A: Instantaneous Amplitude (IA), each row of the Matrix is the instantaneous amplitude
% corresponding to the product function in the Matrix PF
% (3)S: Pure Frequency Modulated Functions Signal(PFMS), each row of the Matrix is the 
% PFMS corresponding to the product function in the Matrix PF, Those PFMS
% can be used to extract Instantaneous Frequency using methods like Hilbert Transform
% or Direct Quadrature Method or CTEO, etc.
% (4)ORT: The othoronality of the PFs and the original signal x
% (5)Nbits: The iteration numbers for sifting each PF
%
% -------------------------------------------------------------
% Option Descriptions:
% (1) MAX_PFS:   maximum number of pfs extracted (default: 10)
% (2) MAX_ITERATIONS:   maximum number of sifting iterations for the computation of each pf (default: 30)
% (3) INTERP_METHOD_EOE:   the interpolation method for empirical optimal method (default {'spline', 'pchip'})
% (4) MIN_EXTREMA: minimum number of extrema of the residual signal for the next decomposition (default: 6)
% (5) ERR_S: The threshold for obtaining a PFMS, when the the s_i(t) satisfied: -1-ERR_S<s_i(t)<1+ERR_S, we can stop the sifting process (default: 0.01)
% (6) DISPLAY: when DISPLAY=0, the function will show the sifting process of PF, otherwise it will not display those details (default: 1)
% (7) MMP:Minimum peak prominence for the function 'findpeaks' (default: 1e-6)
% -----------------------------------------------
% For problems with the code, please contact the author:
% Author: Linshan Jia
% Contact: jialinshan123@126.com
% Xi'an Jiaotong University
% 2018-11-19
%==========================================================================

%% ********************* Main Function: eoe_lmd ******************
[x,MAX_PFS,MAX_ITERATIONS,INTERP_METHOD_EOE,MIN_EXTREMA,ERR_S,DISPLAY,MPP]= init_lmd(varargin{:}); % initialize the LMD
x=x(:)';  % Transform the signal into a row vector
L_x=length(x); % Save the length of the signal
E_x=sum(x.^2);% The energy of the original signal
xx=x; % Save a copy of the original signal
i=0;
Nbits=[];%Save the number of iterations needed for each PF

% pre-allocation
PF=zeros(MAX_PFS+1,L_x);
A=zeros(MAX_PFS,L_x);
S=zeros(MAX_PFS,L_x);

% main loop of lmd
while i<=MAX_PFS
    i=i+1;
    [pf,a,s,N_Gpfi]=lmd_pf(x,ERR_S,MAX_ITERATIONS,INTERP_METHOD_EOE,DISPLAY,MPP);
    % S_abnormal is a logical variable, when it were 1, it means that there
    % are too much iterations for obtaining the PF, there must be something
    % wrong with the sifting process.
    %---------
    % stop criterions are saved in the function: stop_lmd
    % When N_Gpfi==0, it means that the 'pf' is the signal x, and the
    % sifting process is not started.
    
    if N_Gpfi==0 %|| S_abnormal
        PF(i,:)=x;
        break;
    end
    % save the pf, a, s
    PF(i,:)=pf;
    A(i,:)=a;
    S(i,:)=s;
    x=x-pf;
    %-----------
    stop=stop_lmd(x,E_x,MIN_EXTREMA,MPP);% check new x
    if stop
%         i=i+1;
        PF(i,:)=x;
        break;
    end
    %----------
    % save the iteration numbers
    Nbits=[Nbits,N_Gpfi];
end % end of the main loop
% clear the extra space
PF(i+1:MAX_PFS+1,:)=[];
A(i:MAX_PFS,:)=[];
S(i:MAX_PFS,:)=[];
ORT=lmd_io(xx,PF);
%EOF
end
%% **************************************************************
% ************************************************************************
% *****************built-in functions*************************************
% ************************************************************************
%% The 1st function: init_lmd
% Modified based on EMD toolbox written by G. Rilling
function [x,MAX_PFS,MAX_ITERATIONS,INTERP_METHOD_EOE,MIN_EXTREMA,ERR_S,DISPLAY,MPP]= init_lmd(varargin)
x = varargin{1}; % The first arg should be always the orignal signal 
if nargin == 2   % If there are only 2 inputs, it means that the second arg should be a 'struct' type
    if isstruct(varargin{2})
        inopts = varargin{2};
    else
        error('when using 2 arguments the first one is the analyzed signal X and the second one is a struct object describing the options')
    end
elseif nargin > 2 % However, if there are more than two args, it means that you want set parameters using the traditional method: ('PropertyName','PropertyValue')
    try
        inopts = struct(varargin{2:end});
    catch
        error('bad argument syntax')
    end
end

opt_fields={'MAX_PFS','MAX_ITERATIONS','INTERP_METHOD_EOE','MIN_EXTREMA','ERR_S','DISPLAY','MPP'};
% default for stopping, there are 6 options totally.
defopts.MAX_PFS=10;
defopts.MAX_ITERATIONS=30;
defopts.INTERP_METHOD_EOE={'spline','pchip'};
defopts.MIN_EXTREMA=6;
defopts.ERR_S=0.01;
defopts.DISPLAY=1; %display the sifting process
defopts.MPP=1e-6;
%---------------------------
% If there are only one input, or no input
opts = defopts;
if nargin==1
    inopts = defopts;
elseif nargin == 0
  error('not enough arguments')
end
%-------------------------
%If uesers just want to set some of the options
names = fieldnames(inopts);
for nom = names'
  if ~any(strcmpi(char(nom), opt_fields)) % compare the PropertyNames
    error(['bad option field name: ',char(nom)])
  end
  if ~isempty(eval(['inopts.',char(nom)])) % empty values are discarded
    eval(['opts.',upper(char(nom)),' = inopts.',char(nom),';'])
  end
end
%---------------------------
MAX_PFS=opts.MAX_PFS;
MAX_ITERATIONS=opts.MAX_ITERATIONS;
INTERP_METHOD_EOE=opts.INTERP_METHOD_EOE;
MIN_EXTREMA=opts.MIN_EXTREMA;
ERR_S=opts.ERR_S;
DISPLAY=opts.DISPLAY;
MPP=opts.MPP;
%--------------------------
end

%% The 2nd function: stop_lmd
function stop=stop_lmd(x,E_x,MIN_EXTREMA,MPP)
% This function is used to stop the LMD
% the initial value of the variable 'stop' is zero.
stop=0;
% The first stop criterion: If there are not enough extrema, then stop it
[p_max,~]=findpeaks(x,'MinPeakProminence',MPP);
[p_min,~]=findpeaks(-x,'MinPeakProminence',MPP);
N_extrema=length(p_max)+length(p_min);
if N_extrema<=MIN_EXTREMA
    stop=1;
    return;
end
% The second stop criterion: If the residual energy is pretty lettle, the stop it
R_energy= sum(x.^2)/E_x;
if R_energy<0.01 ||  R_energy>1.1
    stop=1;
    return;
end
%EOF
end

%% The 3rd function: lmd_pf
% This function is designed for sifting a pf
function [pf,a,s,k]=lmd_pf(x,ERR_S,MAX_ITERATIONS,INTERP_METHOD_EOE,DISPLAY,MPP)
%ERR_S is the Threshold for PFMS
%MAX_ITERATIONS is the max iterations for a PF sifting process
L=length(x);
a11_temp=ones(1,L);
k=0;
xx=x;
ORT=[];
%---------------------------------------------------------------------------
% When there is not enough extrema, we only produce a one row vector which is equal to the signal x.   it means that the signal x has been a residue which has no extra PF to extract.2012-11-21
[~,locs1]=findpeaks(x,'MinPeakProminence',MPP);
[~,locs2]=findpeaks(-x,'MinPeakProminence',MPP);
index_combine=union(locs1,locs2);
% information 1: Terminating the loop because of no enough extrema
N_index=length(index_combine);
if N_index<=6   % when there are not enough extrema, we should terminate the function
    if DISPLAY
        info1=sprintf('Because of there are not enough extrema for computing the envelope, that is NORMAL ');
        disp(info1);
    end
    
    pf=x;
    a=[];
    s=[];
else    
%--------------------------------------------------------------------------
% Otherwise, we will need to extract the PF from the signal x, and now
% let's start the sifting process.
    while k<MAX_ITERATIONS
        %----------------
        k=k+1;
        % compute the local mean m11 and the amplitude a11
        [m11,a11,upp,low]=lmd_lma(x,k,INTERP_METHOD_EOE); % the outputs "upp" and "low" are just for debuging
        a11_temp=a11_temp.*a11;
        h11=x-m11;
        s11=h11./a11; % PFMS
        % orthogonality test
        %     pf_tep=s11.*a11_temp;
        %     ort=lmd_ort_xy(pf_tep,xx);
        %     ORT=[ORT,ort];
        %     % end of the orthogonality test
        %    -----------------
        x=s11;

        %---------------------------------------------------
        stop_pf=stop_lmd_pf(a11,k,ERR_S,DISPLAY);
        if stop_pf
            break;
        end
    end
    % save a and s
    a=a11_temp;
    s=s11;
    pf=a.*s;
    if DISPLAY
        fprintf('after %d iterations, we get a pf\n\n',k);
    end
end
% EOF
end
%% stop_lmd_pf
function stop_pf=stop_lmd_pf(a11,k,ERR_S,DISPLAY)
% The function is used to determine to stop the sifting process
L=length(a11);
stop_pf=0;
j=k;
%------------------------------------------
% Jonathan S. Smith
a11_new=a11(3:end-3);
H=find((a11_new<=1+ERR_S)&(a11_new>=1-ERR_S));
L_H=length(H);
% 修正这个停止条件，这个条件原来是在包络精度不够的情形下提出的，但是，现在完全可以直接进行到底而不出现错误，因此，我们可以加紧这个条件！！！！
La=length(find(a11_new>1+ERR_S))+length(find(a11_new<1-ERR_S));
if La==0                    %(L_H>0.99*L)
    if DISPLAY
        info1=sprintf('Because of getting a qualified PFMS, that is NORMAL');
        disp(info1);
    end
    stop_pf=1;
    return
end
end
%% The 4th function: lmd_lma
% This function is designed for calculating the local mean and signal
% amplitude, using the EOE algorithm
function [m11,a11,upper_enve,lower_enve]=lmd_lma(x,k,INTERP_METHOD_EOE)
%--------------------------------
Fst_Interp=INTERP_METHOD_EOE{1};
Scd_Interp=INTERP_METHOD_EOE{2};
if k==1
        [upper_enve,lower_enve]=eoe(x,Fst_Interp,Fst_Interp);
else
        [upper_enve,lower_enve]=eoe(x,Scd_Interp,Scd_Interp);
end

% calculate the local mean 'm11' and local amplitude 'a11'
m11=   (upper_enve+lower_enve)/2;
a11=abs(upper_enve-lower_enve)/2;

end

%% The 5th function: eoe
% Empirical Optimal Envelope (EOE) is a novel envelope construvtion method
% which is designed for getting the optimal envelope that is tangential to
% the original signal itself.
function [upper_envelope,lower_envelope,Z11,Z22]=eoe(x,INTERP1,INTERP2)
% -----------------------------------------------------
if nargin==1
    % the default interpolation method is 'spline'
    INTERP1='spline'; 
    INTERP2='spline'; 
end
if nargin>3
    error('Too many inputs, error！')
end

%-------------------------------------------
x=x(:);
x=x';
L=length(x);
t=1:L;
MPP=1e-6;
%-------------------------------------------
[vx_max,kx_max]=findpeaks(x,'MinPeakProminence',MPP);
[vx_min,kx_min]=findpeaks(-x,'MinPeakProminence',MPP);
vx_min=-vx_min;

[kx_max,kx_min,vx_max,vx_min]=boundary_lmd(x,kx_max,kx_min,vx_max,vx_min);
%------------------------------------------------------------
% initial envelope
upper_envelope=interp1(kx_max,vx_max,t,INTERP1);
lower_envelope=interp1(kx_min,vx_min,t,INTERP1);

ERR=0.015;
MAX_N=10;
% ---------------------------upper envelope--------------------------------
counter1=0;
Z1=[];
Z11=[];
while 1
    counter1=counter1+1;
    d_xe=x-upper_envelope;
    [vx_max_dxe,kx_max]=findpeaks(d_xe,'MinPeakProminence',MPP);
    kx_max1=[1,kx_max,L];

    v_upp=x(kx_max1);
    v_upp=[vx_max(1),v_upp(2:end-1),vx_max(end)];

    upper_envelope=interp1(kx_max1,v_upp,t,INTERP2);
    z1=length(find(x-upper_envelope>0));
    Z1=[Z1,z1];
    z11=length(find(vx_max_dxe<0));
    Z11=[Z11,z11];
    
    if (max(vx_max_dxe)<(ERR))&&(min(vx_max_dxe)>(-ERR))
        break;
    end
    
    if counter1>=MAX_N
        break;
    end
    v_upp=[];
end

%---------------------------------
here_upp=find(x-upper_envelope>0);
upper_envelope(here_upp)=x(here_upp);
%----------------------------------

%---------------------------lower envelope---------------------------------
counter2=0;
Z2=[];
Z22=[];
while 1
    counter2=counter2+1;
    d_xe=lower_envelope-x;
    
    [vx_max_dxe,kx_max]=findpeaks(d_xe,'MinPeakProminence',MPP);
    kx_max1=[1,kx_max,L];
    v_low=x(kx_max1);
    v_low=[vx_min(1),v_low(2:end-1),vx_min(end)];
    

    lower_envelope=interp1(kx_max1,v_low,t,INTERP2);
    z2=length(find(lower_envelope-x>0));
    Z2=[Z2,z2];
    
    z22=length(find(vx_max_dxe<0));
    Z22=[Z22,z22];
    
    if (max(vx_max_dxe)<(ERR))&&(min(vx_max_dxe)>(-ERR))
        break;
    end
    if counter2>=MAX_N
        break;
    end
    v_low=[];
end
%-----------------------------
here_low=find(lower_envelope-x>0);
lower_envelope(here_low)=x(here_low);
%-----------------------------

%----------------------------------Modify--------------------------------
de=upper_envelope-lower_envelope;
err_pts=find(de<=0);
if ~isempty(err_pts)
%     warning('there is one or more intersections between upper and lower envelope')
    %----------------------------------------------------------
    for i=1:length(err_pts)
        mdf_v=max([de(err_pts+1),de(err_pts-1)]);
        if mdf_v<0
            mdf_v=-mdv_v;
        elseif mdf_v==0
            mdf_v=1/max(x);
        end
        mdf_k=err_pts(i);
        upper_envelope(mdf_k)=upper_envelope(err_pts(i))+0.5*mdf_v;
        lower_envelope(mdf_k)=lower_envelope(err_pts(i))-0.5*mdf_v;
    end
end
% debug
% disp(counter1)
% disp(counter2)
end

%% The 6th function: boundary_lmd
% this function is used to deal with the end effect of lmd and eoe
% algorithms.
function [kx_max,kx_min,vx_max,vx_min]=boundary_lmd(x,kx_max,kx_min,vx_max,vx_min)
L=length(x);
%---------------------left side--------------------------
if kx_max(1)<kx_min(1)
    
    if x(1)>vx_min(1)
        vx_max=[vx_max(1),vx_max];
        vx_min=[vx_min(1),vx_min];
    
    else
        vx_max=[vx_max(1),vx_max];
        vx_min=[x(1),vx_min];
    end

else  
    
    if x(1)>vx_max(1)
        vx_max=[x(1),vx_max];
        vx_min=[vx_min(1),vx_min];
    
    else
        vx_max=[vx_max(1),vx_max];
        vx_min=[vx_min(1),vx_min];
    end

end

%---------------------right side--------------------------
if kx_max(end)>kx_min(end)
    if x(end)<vx_min(end)
        vx_min=[vx_min,x(end)];
        vx_max=[vx_max,vx_max(end)];
    else
        vx_min=[vx_min,vx_min(end)];
        vx_max=[vx_max,vx_max(end)];
    end

else  
    if x(end)>vx_max(end)
        vx_min=[vx_min,vx_min(end)];
        vx_max=[vx_max,x(end)];
    else
        vx_min=[vx_min,vx_min(end)];
        vx_max=[vx_max,vx_max(end)];
    end
end
%------------------------------
kx_min=[1,kx_min,L];
kx_max=[1,kx_max,L];
end

%% The 7th function: Compute the index of orthogonality
% ** Copied from emd toolbox by G.Rilling and P.Flandrin
% ** http://perso.ens-lyon.fr/patrick.flandrin/emd.html
function ort = lmd_io(x,pfs)
% ort = IO(x,pfs) computes the index of orthogonality
% inputs : - x   : analyzed signal
%          - pfs  : production function
n = size(pfs,1);
s = 0;
for i = 1:n
    for j =1:n
        if i~=j
            s = s + abs(sum(pfs(i,:).*conj(pfs(j,:)))/sum(x.^2));
        end
    end
end
ort = 0.5*s;
end
%% The 8th function: Compute the orthogonality between x and y
% according to the Paper:
% Huang N E, Shen Z, Long S R, et al. The Empirical Mode Decomposition and 
% the Hilbert Spectrum for Nonlinear and Non-Stationary Time Series Analysis[J].
% Proceedings Mathematical Physical & Engineering Sciences, 1998, 454(1971):903-995.
function OC=lmd_ort_xy(x,y)
x=x(:)';
y=y(:);
xy=x*y;
OC=xy/(sum(x.^2)+sum(y.^2));
end
