function H=Generate_Channel_frequency_selective_LTV(f,TX_pos,RX_pos,scenario,v_RX,v_TX,Yt,Zt,Yr,Zr,a,Ts,Tc)
% Function for the generation of a MIMO channel matrix at mmWaves
% frequencies in case of linear time invariant channel. This function
% refers to equation (6) in:
%
% S. Buzzi, C. D'Andrea , "A Clustered Statistical MIMO Millimeter Wave
% Channel Model",submitted to IEEE Wireless Communications Letters
%
% License: This code is licensed under the GPLv2 License.If you in any way 
% use this code for research that results in publications, please cite our
% original article listed above.

%% INPUT PARAMETERS
% f: Carrier Frequency (e.g., 73e9 Hz)

% TX_pos: 3-dimensional vector for the position of the transmitter (e.g.,
% [0 0 7] for a transmitter placed at the origin at 7 meters height)

% RX_pos: 3-dimensional vector for the position of the receiver (e.g.,
% [30 0 1] for a transmitter placed at distance 30m from the transmitter at 
% 1 meters height)

% scenario: variable that contains information about the use-case scenario,
% it assumes the values:
% - scenario==1  ==> 'Open square'
% - scenario==2  ==> 'Street Canyon'
% - scenario==3  ==> 'Indoor Office'
% - scenario==4  ==> 'Shopping mall'
% Yt and Zt: number of antennas on the y and x axes of the transmitter
% planar array (e.g., 10 and 5 for a planar array with 50 antennas)

% Yr and Zr: number of antennas on the y and x axes of the receiver
% planar array (e.g., 10 and 8 for a planar array with 80 antennas)

% a: is a vector contains the samples of the convolution between the
% transmit and receive shaping pulses sampling with step Ts (e.g., raised
% cosine pulse, for the convolution between two RRC shaping pulse)

% Ts: sampling time for the filter 

% Tc: sampling time for the channel( e.g, 1e-9 s)

% N.B: Tc should never be smaller than Ts to avoid bad precision results. An
% interpolation routine is indeed used in order to come with missing
% samples


%% OUTPUT Parameter
% H: 3-dimensional array of dimension N_R x N_T x P representing the
% samples of the matrix-valued channel impulse response taken with a 
% sampling frequency equal to 1/Tc (e.g., Tc=1e-9 s) 

%% Initial configurations

Nr=Yr*Zr; % number of antennas of the receiver array
Nt=Yt*Zt; % number of antennas of the transmitter array

lambda=3e08/f; % wavelength
d_array=lambda/2; % distance between elements of planar array( we consider half-wavelength arrays)
k=2*pi/lambda; % wavenumber
kd=k*d_array;

% standard deviations of Laplacian distribution ( we consider 5�) 
dev_standard_phit=5/180*pi;
dev_standard_thetat=5/180*pi;
dev_standard_phir=5/180*pi;
dev_standard_thetar=5/180*pi;

% variance of the complex gain of the paths
alpha_variance=1;

% heigths of receiver and transmitter
hr=RX_pos(3);
ht=TX_pos(3);

% coordinate in the 3D space with TX in the origin of reference system

x=RX_pos(1)-TX_pos(1);
y=RX_pos(2)-TX_pos(2);
z=RX_pos(3)-TX_pos(3);

% Polar Coordinate

d=sqrt(x^2+y^2+z^2); % distance between transmitter and receiver
theta_r=atan(z/d);
phi_r=atan(y/x);


%% Parameter depending to scenario
% These expressions are in the equations (6) and (7) in he article listed
% above
if (scenario==1 || scenario==2)
   P_LOS=min(20/d,1)*(1-exp(-d/39))+exp(-d/39);

elseif (scenario==3 || scenario==4)
    if d<1.2
        P_LOS=1;
    elseif d>2 && d<6.5
        P_LOS=exp(-(d-1.3)/4.7);
    elseif d>=6.5
        P_LOS=0.32*exp(-(d-6.5)/32.6);
    end
else
    disp('ERROR: INVALID SCENARIO');
    H=NaN;
    return
end

%% Calculation of Number of scattering cluster depending from distance between transmitter and receiver
lambda_cluster=1.9; % mean of number of cluster (Poisson random variable)
Ncl=max(1,poissrnd(lambda_cluster)); % number of cluster
Nray=randi(30,Ncl,1); % Number of rays for each cluster for each cluster
gamma=sqrt(Nr*Nt/(sum(Nray,1)));

for ii=1:Ncl,
    
    % Generation of position in azimuth and elevation of the i-th cluster
    % as random variables with uniform distribution in appropriate
    % intervals
    phit_pos_i=rand*pi-pi/2; % uniform in [-pi/2,pi/2]
    thetat_pos_i=rand*pi-pi/2; % uniform in[-pi/2,pi/2]
    phir_pos_i=rand*pi; % uniform in [0,2*pi]
    thetar_pos_i=rand*pi-pi/2; % uniform in[-pi/2,pi/2]
    % Definition of Rmax as the maximum r-position of the cluster to avoid
    % considering a cluster underground
    
    Rmax=sqrt(ht^2+(ht*tan(pi/2+thetat_pos_i))^2);
    
    % Definition of r-position of the cluster
    if thetat_pos_i>0
        r=rand*(7*d/4 -1)+1; % uniform random variable in (1,7/4*d)
    else
        r=rand*(min(7*d/4,Rmax)-1)+1; % uniform random variable in (1,(min(7*d/4,Rmax)) to avoid considering a cluster underground
    end
        
    tau_temp=zeros(Nray(ii,1),1);
    phit_temp=zeros(Nray(ii,1),1);
    phir_temp=zeros(Nray(ii,1),1);
    thetat_temp=zeros(Nray(ii,1),1);
    thetar_temp=zeros(Nray(ii,1),1);
    r_temp=zeros(Nray(ii,1),1);
    
    for ll=1:Nray(ii,1),

        % Generation of the position of l-th path in the i-th cluster as a
        % Laplace random variable with mean as the position of the i-th
        % cluster in azimuth and elevation for receiver and transmitter for
        % calculation of path delay
        thetat_temp(ll,1)=Laplace_distribution(dev_standard_thetat)+thetat_pos_i;
        phit_temp(ll,1)=Laplace_distribution(dev_standard_phit)+phit_pos_i;
        thetar_temp(ll,1)=Laplace_distribution(dev_standard_thetar)+thetar_pos_i;
        phir_temp(ll,1)=Laplace_distribution(dev_standard_phir)+phir_pos_i;        
        
        % Calculation of the length of the path ii-ll as in equation (2) in
        % the article listed above
        
        x_s=r*cos(thetat_temp(ll,1))*cos(phit_temp(ll,1));
        y_s=r*cos(thetat_temp(ll,1))*sin(phit_temp(ll,1));
        z_s=r*sin(thetat_temp(ll,1));
        
        r_temp(ll,1)=r + norm([(x_s-x), (y_s-y), (z_s-z)]);
        
        
        % calculation of delay of the path ii-ll 
        
        tau_temp(ll,1)=(r_temp(ll,1))/3e08;
    end
    Path_delays{ii}=tau_temp;
    Path_lengths{ii}=r_temp;
    Thetat{ii}=thetat_temp;
    Phit{ii}=phit_temp;
    Thetar{ii}=thetar_temp;
    Phir_temp{ii}=phir_temp;
end
%% Computation of channel matrix

% maximum delay
max_delay_vector=zeros(Ncl,1);
for ii=1:Ncl
    temp=Path_delays{ii};
    max_delay_vector(ii,1)=max(temp);
end
max_delay=max(max_delay_vector);

tt=0:Tc:max(max_delay);
dim_t=length(tt);
% Introduction of correlation in time with a correlation matrix

R=zeros(dim_t);
for hh=1:dim_t
    for pp=1:dim_t
        R(hh,pp)=0.995^(abs(hh-pp));
    end
end

for ii=1:Ncl
   % calculation of the temporary complex path gain
    alpha_temp=(randn(dim_t,Nray(ii,1))+1j*randn(dim_t,Nray(ii,1)))*sqrt(alpha_variance/2); 
    alpha_Rays=zeros(dim_t,Nray(ii,1));
    for ll=1:Nray(ii,1)
        alpha_Rays(:,ll)=sqrtm(R)*alpha_temp(:,ll);
    end 
    Path_gains{ii}=alpha_Rays;
end

PP=ceil(max_delay/(Tc) + length(a)*Ts/Tc*100); % definition of temporary length for the matrix H_temp
H_temp=zeros(Nr,Nt,PP,dim_t);

%% Consideration of LOS path as in equation (10) in the article listed above
delta=sqrtm(R)*rand(dim_t,1)*2*pi; % expression for delta in equation (7)
for t=1:dim_t
    u=rand;    
    if u<=P_LOS
        alpha_LOS=exp(1j*delta(t));
        attenuation_LOS=Evaluation_Path_loss(d,f,scenario,1);
        theta_t_LOS=theta_r;
        theta_r_LOS=-theta_r;
        phi_t_LOS=phi_r;
        phi_r_LOS=-phi_r;
        
        % express the speeds of receiver in m/s
        v_TX_ms=v_TX*1e3/3600;
        v_RX_ms=v_RX*1e3/3600;
        ni_LOS=-f/3e8*(v_TX_ms*cos(phi_t_LOS)*cos(theta_t_LOS)+v_RX_ms*cos(phi_r_LOS)*cos(theta_r_LOS));
        tau_LOS=d/3e8;
        
        % Array responses as in article listed above
        at=Array_response(Yt,Zt,phi_t_LOS,theta_t_LOS,kd);
        ar=Array_response(Yr,Zr,phi_r_LOS,theta_r_LOS,kd);
        
        % Calculation of LOS component as in equation (7) in the article
        % listed above
        Q_LOS=sqrt(Nr*Nt)*alpha_LOS*sqrt(attenuation_LOS)*exp(-1j*2*pi*ni_LOS*tt(t))*ar*at';
        
        % product between temporary channel matrix and the interpolated
        % filter a
        gg=floor(tau_LOS/Tc); % index of filter "approximately" corresponding to the delay tau_matrix(ii,ll)
        gg_up=floor((tau_LOS+(length(a)-1)*Ts)/Tc);% index of filter "approximately" corresponding to the delay tau_matrix(ii,ll)+(length(a)-1)*Ts
        
        % Interpolation of filter a
        Ts_step=0:Ts:(length(a)-1)*Ts;
        Tc_step=(gg+1)*Tc-tau_LOS:Tc:gg_up*Tc-tau_LOS;
        a_interp=interp1(Ts_step,a,Tc_step,'spline');
        % Construction of channel matrix corresponding to the path i-th
        % l-th
        for ee=1:length(a_interp)
            H_temp(:,:,ee+gg,t)= H_temp(:,:,gg+ee,t) + Q_LOS .*a_interp(ee);
        end

    end

    for ii=1:Ncl,
        
        alpha_ii=Path_gains{ii};
        tau_ii=Path_delays{ii};
        r_ii=Path_lengths{ii};
        theta_ii=Thetat{ii};
        phit_ii=Phit{ii};
        thetar_ii=Thetar{ii};
        phir_ii=Phir_temp{ii};
        
        for ll=1:Nray(ii,1),   
            
            phit_ll=phit_ii(ll,1);
            phir_ll=phir_ii(ll,1);
            thetat_ll=theta_ii(ll,1);
            thetar_ll=thetar_ii(ll,1);
            tau_ll=tau_ii(ll,1);
            r_ll=r_ii(ll,1);
            alpha_ll=alpha_ii(t,ll);

            % Array responses as in article listed above
            at=Array_response(Yt,Zt,phit_ll,thetat_ll,kd);
            ar=Array_response(Yr,Zr,phir_ll,thetar_ll,kd);

            % Calculation of the attenuation of each path as in equation (2)
            % in article listed above
            attenuation_il=Evaluation_Path_loss(r_ll,f,scenario,0);
            
            % expression of the speeds of receiver in m/s
            v_TX_ms=v_TX*1e3/3600; 
            v_RX_ms=v_RX*1e3/3600;
            % calculation of Doppler shift associated to the path (i,l)^th
            ni_il=-f/3e8*(v_RX_ms*cos(thetar_ll)*cos(phir_ll)+v_TX_ms*cos(thetat_ll)*cos(phit_ll)); % doppler shift
            % Temporary channel matrix
            QQ=gamma*alpha_ll*sqrt(attenuation_il)*exp(-1j*2*pi*ni_il)*ar*at';

            % product between temporary channel matrix and the interpolated
            % filter a
            gg=floor(tau_ll/Tc); % index of filter "approximately" corresponding to the delay tau_matrix(ii,ll)
            gg_up=floor((tau_ll+(length(a)-1)*Ts)/Tc);% index of filter "approximately" corresponding to the delay tau_matrix(ii,ll)+(length(a)-1)*Ts
            % Interpolation of filter a
            Ts_step=0:Ts:(length(a)-1)*Ts;
            Tc_step=(gg+1)*Tc-tau_ll:Tc:gg_up*Tc-tau_ll;
            a_interp=interp1(Ts_step,a,Tc_step,'spline');

            % Construction of channel matrix corresponding to the path i-th
            % l-th
            for ee=1:length(a_interp)
                H_temp(:,:,ee+gg,t)= H_temp(:,:,gg+ee,t) + QQ .*a_interp(ee); 
            end
        end
    end
end

% Selection of effective channel matrix with the identification of indmin
% and indmax that indicate the start and end of discrete channel response,
% so we consider the relative delay as explained in article listed above
indmin=1;
indmax=PP;

% Calculation of Frobenius norm of matrices in third dimension
Frobenius_norm=zeros(PP,1);
for pp=1:PP
    for t=1:dim_t
       Frobenius_norm(pp,1)=Frobenius_norm(pp,1)+norm(H_temp(:,:,pp,t),'fro'); 
    end
end

 % We choice as indmin the firts index of matrix H_temp in which we have a
 % non zero contribute 
for pp=1:PP
   if Frobenius_norm(pp,1)~=0
       indmin=pp;
       break
   end
end

 % We choice as indmaxn the last index of matrix H_temp in which we have a
 % non zero contribute 
for qq=0:PP-1
   if Frobenius_norm(PP-qq,1)~=0 
       indmax=PP-qq;
       break
   end
end

% Selection of effective channel matrix
H=H_temp(:,:,indmin:indmax,:);
end