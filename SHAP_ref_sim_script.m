
%% 读取网络结构、稳态优化参数等
fnameid=fopen('model_folder.txt');
fname = textscan(fnameid, '%s');
fclose(fnameid);
par.mfolder=fname{1}{1};
par=options_input(par);

% 读取仿真模型和优化算法的参数，不涉及气网的结构，全部返回到par中，在仿真中则不应用优化参数，仅仿真
par=options_input(par);

%load static model (nodes, pipes, comps, gnodes) from xls fliles，
% 读取稳态网络对象，把所有信息存到txt中，若txt存在若存在txt则从新的xls读取，没有xls则从csv，在保存成新的txt
if(exist([par.mfolder '\input_network.txt'])~=2 || par.out.update_from_xls==1)
    if(exist([par.mfolder '\input_network_nodes.xls'])==2 || exist([par.mfolder '\input_network_nodes.xlsx'])==2), [par.ss.n0]=gas_model_reader_xls(par); else
        [par.ss.n0]=gas_model_reader_csv(par); end
end

%first input check
par=check_input_1(par);
if(par.flag==1), disp(par.message), return; end

%% 稳态ss优化
if(par.out.doss==1 || par.out.dosim==1)
    % 单位转换 保存一个针对ss优化的气网对象数组par.ss.n0---划分网格之前的网格结构对象
    [par.ss.n0]=gas_model_reader_new(par.mfolder);   % load data from text file 从新保存的txt中
    % 再次check对象是否过大
    if(par.ss.n0.nv>par.out.maxnv || par.ss.n0.ne>par.out.maxne || par.ss.n0.nc>par.out.maxng || par.ss.n0.ng>par.out.maxng), return; end

    %if(par.ss.n0.nv>20), return; end;
    % 划分求解网格后重新定义nodes等
    [par.ss.n]=gas_model_reconstruct_new(par.ss.n0,par.tr.lmax,1);

    %optimization parameters
    par.ss.m.cdw=50;                      %comp. ratio derivative penalty scale
    par.ss.m.ddw=50;                      %flex demand derivative penalty scale
    par.ss.m.odw=100;                     %objective scale
    par.ss.m.maxiter=400;                 %
    par.ss.m.opt_tol=1e-6;                %

    %model specifications
    % 归一化、进一步的参数定义、定义模型的一些细节、约束等
    [par.ss]=model_spec(par.ss);

    %demand function specifications
    par.ss=econ_spec(par.ss,par.mfolder);

    % Solve optimization
    [par.ss_start]=static_opt_base_ends(par.ss,1); %起点
    par.ss=par.ss_start;...同上
        [par.ss_terminal]=static_opt_base_ends(par.ss,-1); %终点
    [par]=process_output_ss_nofd(par);
    [par]=process_output_ss_terminal_nofd(par);
    [par]=process_output_ss_start_nofd(par);

    %exit if not solved
    if(par.ss.ip_info.status~=0), disp('Steady state optimization not feasible'),
        fid=fopen([par.mfolder '\output_log.txt'],'w');
        fprintf(fid,['Steady-state solve status: ' num2str(par.ss.ip_info.status) '\n']);
        fclose(fid);
    end

    %process steady-state output 后处理稳态优化的结果
    % if(par.out.steadystateonly==1), [par]=process_output_ss(par); if(par.out.intervals_out>0), gas_out_plots_i(par); end; return; end
    % if(par.out.ss_check_exit==1), return; end
end

%% 从ss优化结果出发，以tr里给出的边界信息，开始仿真

% 定义tr.n n0 容纳tran_sim网络对象
[par.tr.n0]=gas_model_reader_new(par.mfolder);                   % load data from text file
if(par.ss.n0.nv>par.out.maxnv || par.ss.n0.ne>par.out.maxne || par.ss.n0.nc>par.out.maxng || par.ss.n0.ng>par.out.maxng), return; end

[par.tr.n]=gas_model_reconstruct_new(par.tr.n0,par.tr.lmax,1);

%model specifications
[par.tr]=model_spec(par.tr);

%demand function specifications
par.tr=econ_spec(par.tr,par.mfolder); % 实际以上准备工作和ss中的东西完全一样，只不过为了区分

% 准备：setup 0：边界条件：气源压力、出口流量、压缩机压比

par.sim=par.ss;
par.sim.rtol0=1e-2; par.sim.atol0=1e-1;
par.sim.rtol1=1e-3; par.sim.atol1=1e-2;
par.sim.rtol=1e-5; par.sim.atol=1e-3;  %error tolerances for simulation
%par.sim.startup=1/4.2;     %startup time (fraction of horizon)
par.sim.startup=1/8;        %startup time (fraction of horizon)
par.sim.nperiods=2;         %number of periods after startup
par.sim.solsteps=24*6*2;        %solution steps per period
par.sim.fromss=1;

% 压缩机动作
% load('E:\working\Matlab_prjs\Gas_Line_case\src\ccccc.mat');
% 1.60000000000000
% 1.18399927993301
% 1.16290675532091
% 1.10910660796442
% 1.00097383792682

% 1.60000000000000
% 1.22613672921252
% 1.20293731120827
% 1.20293731182393
% 1.00776213532559

cc0_start = par.ss_start.cc0(:,1);
[all_n_pr,~] =size(cc0_start);
cc0_terminal = par.ss_terminal.cc0(:,1);
par.id_v = find(abs(cc0_start - cc0_terminal)>=0.005);
par.id_uv = find(abs(cc0_start - cc0_terminal)<0.005);
cc0_start_iterp = cc0_start(par.id_v,:);
cc0_terminal_iterp = cc0_terminal(par.id_v,:);
% interval_pr = cc0_terminal_iterp-cc0_start_iterp;
num_cpoints=5;
upbound=zeros(4,5);
lowbound=zeros(4,5);
for i=1:length(cc0_start_iterp) % 四个压缩机 五个插值点
    centor = linspace(cc0_start_iterp(i),cc0_terminal_iterp(i),7); % 区间中心
    centor = centor(2:end-1);
    interval = centor(2)-centor(1); % 基本区间长度
    upbound(i,:) = centor+interval*0.49;
    lowbound(i,:) = centor-interval*0.49;
end
interval_pr=upbound-lowbound;


% id_1 = find(abs(interval_pr)>=0.005 & abs(interval_pr)<0.02);
% id_2 = find(abs(interval_pr)>=0.02 & abs(interval_pr)<0.1);
% id_3 = find(abs(interval_pr)>=0.1 & abs(interval_pr)<0.5);
%% 采样，并根据前后端ss优化结果，形成压缩机动作序列
n_control_tps = 5;
k = 10;
[n_pr,~] = size(par.id_v);
doe = lhsdesign(100,8);
doe = reshape(doe,100,4,2);
finalv_doe = zeros(100,5,12);
p1=12;
p2=20;
for i=1:100
    a = reshape(doe(i,:,:),4,2);
    len = (cc0_terminal-cc0_start)./2/0.2;

    t1 = [0,4,8,12,16,24];
    v1 = [1.6,1.6,1.6,1.6,1.6,1.6];
    t2_1 = a(1,1).*p1;
    t2_2 = a(1,2).*(p2-t2_1-len(2));

%     1
    t3_1 = a(2,1).*(p1);
    t3_2 = a(2,2).*(p2-t3_1-len(3));

    t4_1 = a(3,1).*(p1-t2_1)+t2_1;
    t4_2 = a(3,2).*(p2-t4_1-len(4));

    t5_1 = a(4,1).*(p1-t2_1)+t2_1;
    t5_2 = a(4,2).*(p2-t5_1-len(5));

    %2
%     t3_1 = a(2,1).*(p1);
%     t3_2 = a(2,2).*(p2-t3_1-len(3));
% 
%     t_ref=max(t3_1,t2_1);
% 
%     t4_1 = a(3,1).*(p1-t_ref)+t_ref;
%     t4_2 = a(3,2).*(p2-t4_1-len(4));
% 
%     t5_1 = a(4,1).*(p1-t_ref)+t_ref;
%     t5_2 = a(4,2).*(p2-t5_1-len(5));

    %3
%     t3_1 = a(2,1).*(p1);
%     t3_2 = a(2,2).*(p2-t3_1-len(3));
% 
%     t4_1 = a(3,1).*(p1);
%     t4_2 = a(3,2).*(p2-t4_1-len(4));
% 
%     t_ref1=max(t4_1,t3_1);
%     t_ref2=max(t2_1,t3_1);
%     t_ref=max(t_ref1,t_ref2);
% 
%     t5_1 = a(4,1).*(p1-t_ref)+t_ref;
%     t5_2 = a(4,2).*(p2-t5_1-len(5));


    t2=[0,t2_1,t2_1+len(2),t2_1+len(2)+t2_2,t2_1+len(2)+t2_2+len(2),24];
    v2=[cc0_start(2),cc0_start(2),cc0_start(2)+len(2).*0.2,cc0_start(2)+len(2).*0.2,...
        cc0_start(2)+len(2).*0.4,cc0_start(2)+len(2).*0.4];
    t3=[0,t3_1,t3_1+len(3),t3_1+len(3)+t3_2,t3_1+len(3)+t3_2+len(3),24];
    v3=[cc0_start(3),cc0_start(3),cc0_start(3)+len(3).*0.2,cc0_start(3)+len(3).*0.2,...
        cc0_start(3)+len(3).*0.4,cc0_start(3)+len(3).*0.4];
    t4=[0,t4_1,t4_1+len(4),t4_1+len(4)+t4_2,t4_1+len(4)+t4_2+len(4),24];
    v4=[cc0_start(4),cc0_start(4),cc0_start(4)+len(4).*0.2,cc0_start(4)+len(4).*0.2,...
        cc0_start(4)+len(4).*0.4,cc0_start(4)+len(4).*0.4];
    t5=[0,t5_1,t5_1+len(5),t5_1+len(5)+t5_2,t5_1+len(5)+t5_2+len(5),24];
    v5=[cc0_start(5),cc0_start(5),cc0_start(5)+len(5).*0.2,cc0_start(5)+len(5).*0.2,...
        cc0_start(5)+len(5).*0.4,cc0_start(5)+len(5).*0.4];

    t_ay=[t1;t2;t3;t4;t5];
    v_ay=[v1;v2;v3;v4;v5];
    finalv_doe(i,:,1:6) = t_ay;
    finalv_doe(i,:,7:12) = v_ay;
end

% 进行spline插值
% load('E:\working\Matlab_prjs\Gas_Line_case\src\SHAPtree_data.mat')
n_control_tps=25;
all_data.('id_v')=par.id_v;
all_data.('id_uv')=par.id_uv;
plotpos=[20,20,450,450];
for i=1:100
    cc00 = squeeze(finalv_doe(i,:,7:12));
    tt0=  squeeze(finalv_doe(i,:,1:6));

    % 插值成为25h长度
    xx_new = linspace(1,n_control_tps,25)'-1;
    yi = zeros(5,25);
    for iii=1:5
        yi(iii,:)=interp1(tt0(iii,:),cc00(iii,:),xx_new,'linear');
    end
    cc0=yi;
%     cc_uv=yi(1,:);
%     cc0(cc0<=1)=1;
%     cc0=x_all(i,:,:);
%     cc_v=squeeze(cc0);
%     cc0=zeros(5,25);
%     cc0(1,:)=cc_uv;
%     cc0(2:5,:)=cc_v';
%     f1=figure(1); clf
%     set(f1,'position',plotpos,'Color',[1 1 1]);
%     subaxis(1,1,1,'MarginLeft',0.1,'SpacingHoriz',0.05,'MarginRight',0.1),
%     for i_plot=2:5
%         plot(tt0(i_plot,:),cc00(i_plot,:),'LineWidth',3), axis('tight'), xlabel('Time(h)')
%         hold on
%     end
%     title('Compressor ratio','fontweight','bold')
%     legend('Comp_2#','Comp_3#','Comp_4#','Comp_5#');

    [par]=tran_sim_setup_0(par,cc0);
    % execute simulation
    [par.sim]=tran_sim_base_flat_noextd(par.sim);
    [par]=process_output_tr_nofd(par);
    all_data.(['x',num2str(i)])=cc0;
    all_data.(['y',num2str(i)])=par;
    all_data.(['m_cost',num2str(i)])=par.tr.m_cost(end-par.sim.solsteps:end,:);
    all_data.(['m_var',num2str(i)])=par.tr.m_var(end-par.sim.solsteps:end,:);
    all_data.(['m_mass',num2str(i)])=par.tr.m_mass(end-par.sim.solsteps:end,:);
    %   all_data.(['m_ts',num2str(i)])=par.tr.m_ts;
    disp(['simulating ' ' sample... ' num2str(i) '..................................' ])
end
%     shap_f = @(cc)tran_sim_shap(cc,par); % 输入cc 5*25 = 125 展开
%     Xall=zeros(k*n_pr*n_control_tps,25*all_n_pr); % 注意对应的顺序
%     for i=1:k*n_pr*n_control_tps
%         Xall(i,:)=reshape(all_data.(['x',num2str(i)]),1,125);
%     end
%     shap_explainer = shapley(shap_f,Xall);
%     ex = fit(shap_explainer,Xall(1,:));

% 序列树模型构建。MC测试。两种对象。

% 滚动仿真，3*5一个动作？，套用强化学习智能体（可以根据对象变动，类似这种管理），开写。R里边只有一项，约束；
% 转供对象四台压缩机，ADMM优化case搞定。开写。
% 如何调用构建Python或matlab的树


%% 后处理，作图
% 保存shap生成tree的数据集
par.id_v=all_data.('id_v');
par.id_uv=all_data.('id_uv');
shap_cost = zeros(100,25,4);
shap_var = zeros(100,25,4);
shap_mass = zeros(100,25,4);
shap_supp = zeros(100,25,4);
shap_error = zeros(100,4,4);
x_all = zeros(100,25,4);
cost = zeros(100,289);
mass = zeros(100,289);
var = zeros(100,289);
for i=1:100
    xi = all_data.(['x',num2str(i)])';
    xi = xi(:,par.id_v);
    x_all(i,:,:)=xi;
    yi = all_data.(['y',num2str(i)]);
    shap_cost(i,:,:)=yi.shap_cost(:,par.id_v);
    shap_var(i,:,:)=yi.shap_var(:,par.id_v);
    shap_mass(i,:,:)=yi.shap_mass(:,par.id_v);
    shap_supp(i,:,:)=yi.shap_supp(:,par.id_v); % 与以上不同
    shap_error(i,:,:)=yi.shap_error(:,par.id_v);
    cost(i,:)=all_data.(['m_cost',num2str(i)]);
    mass(i,:)=all_data.(['m_mass',num2str(i)]);
    var(i,:)=all_data.(['m_var',num2str(i)]);
end
SHAPtree_data_ref1.cost=cost;
SHAPtree_data_ref1.mass=mass;
SHAPtree_data_ref1.var=var;
SHAPtree_data_ref1.shap_cost=shap_cost;
SHAPtree_data_ref1.shap_mass=shap_mass;
SHAPtree_data_ref1.shap_var=shap_var;
SHAPtree_data_ref1.shap_error=shap_error;
SHAPtree_data_ref1.x_all=x_all;
cost_mean=mean(sum(cost,2));
mass_mean=mean(sum(mass,2));
[par]=process_output_tr_nofd(par);

if(par.out.intervals_out==0), gas_out_plots_nofd(par); end
if(par.out.intervals_out>0), gas_out_plots_i(par); end

%pause(inf)
if(par.out.closeafter==1), close all, exit, end