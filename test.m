%% 基本定义
% phi是闭环模型
% G是开环模型
% 存在等式 phi = G/(1+G)
% C_v是速度控制器 phi_i是电流闭环模型 G_iv是机械动力学
% 存在等式 G = C_v*phi_i*G_iv

clear
close all

s = tf('s')

% 生成激励信号
F1 = ([1:0.5:22, 24:2:40, 50:10:120,200,250,333,500]);%先给出频率点
T1 = 1000./F1;
T = round(T1);%周期取整
F = 1000./T;       

% 导入输入-输出数据
data = load('data/vofa+_radar_motor_25_5Hzbase_1Hzapm_spdP0.05I0.002_cur300.csv');
figure
input = data(:,1);
output = data(:,2);
plot(input);hold on
plot(output);hold on
legend('input','output')

% 导入幅值-相位数据，转化为idfrd模型
% data = load('data/radar_motor_0.5A_cur300.log');
% W = data(:,1)*2*pi;
% % FRE = data(:,1);
% AMP = data(:,2);
% PHA = data(:,3);
% figure
% subplot(211)
% plot(AMP);
% 
% subplot(212)
% plot(PHA);
% response = AMP.*exp(1j*PHA*pi/180);
% cl = idfrd(response,W,0) % 连续
% figure
% bode(cl)

% 根据位置闭环模型反求速度闭环模型（连续）
% load sysid.mat
% % bode(pitch_0_pos_cl)

% P = 5.0;
% deg2rad = pi/180;
% Hz2Deg = 360/(1*4);
% Gain_all = deg2rad*P*Hz2Deg;
% Gpcl = pos_cl;
% figure
% bode(Gpcl)
% Gvcl = Gpcl/(Gain_all*1/s*(1-Gpcl))
% figure
% bode(Gvcl)

% % 根据位置闭环模型反求速度闭环模型（离散）
% % load sysid.mat
% % % bode(pitch_0_pos_cl)
% s = tf('s')
% P = 5.0;
% deg2rad = pi/180;
% Hz2Deg = 360/(1*4);
% Gain_all = deg2rad*P*Hz2Deg;
% Gpcl = pitch_0_pos_cl;
% figure
% bode(Gpcl)
% int = 1/s;
% intd = c2d(int,0.001);
% Gvcl = Gpcl/(Gain_all*intd*(1-Gpcl))
% figure
% bode(Gvcl)

% 查看工具箱输出的模型
% load("sysid.mat")
% [b,a] = ss2tf(ss4.A,ss4.B,ss4.C,ss4.D);
% Phi=tf(b,a)
% figure
% bode(Phi)
%% 
% G = Phi/(1-Phi)
% bode(G)
% margin(G)
% [Gm,Pm,Wcg,Wcp] = margin(G)
% 
% Kp_v = 0.05;
% Ki_v = 0.002;
% C_v = (Kp_v * s + Kp_v * Ki_v)/s
