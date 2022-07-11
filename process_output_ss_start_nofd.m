function [par]=process_output_ss_start_nofd(par)
% Anatoly Zlotnik, February 2019

% mfolder=par.mfolder;
out=par.out;
% 
ss_start=par.ss_start; 

psi_to_pascal=ss_start.c.psi_to_pascal;
mpa_to_psi=1000000/psi_to_pascal;
ss_start.c.mpa_to_psi=mpa_to_psi;
mmscfd_to_kgps=ss_start.c.mmscfd_to_kgps;
hp_to_watt=745.7;
if(par.out.doZ==1), b1=ss_start.c.b1; b2=ss_start.c.b2; end

%process_start optimization output
%process_start dimensional solution (metric)
pp=zeros(length(ss_start.tt0),ss_start.n0.nv);
%s_nodal=[ss_start.m.pslack ss_start.m.pslack(:,1)];
s_nodal=ss_start.m.pslack;
if(par.ss_start.m.extension>0), s_nodal=ss_start.m.pslack; end
scf=find(ismember(ss_start.m.snodes,ss_start.m.comp_pos(ss_start.m.spos,1)));
s(scf,:)=ss_start.cc0(ss_start.m.spos,:).*s_nodal(scf,:);
pp(:,ss_start.m.snodes)=s';
pp(:,ss_start.m.dnodes)=ss_start.pp0';
qq=ss_start.qq0'*ss_start.m.Xs*ss_start.c.qsc; 
ff=ss_start.qq0'*ss_start.c.qsc; 
for j=1:length(ss_start.m.dpos)
    cpj=ss_start.n.comp_pos(ss_start.m.dpos(j),1);
    pp(:,cpj)=pp(:,cpj).*ss_start.cc0(ss_start.m.dpos(j),:)';
end
%nodal press_starture (before compress_startors)
pp_nodal=zeros(length(ss_start.tt0),ss_start.n.nv);
pp_nodal(:,ss_start.m.snodes)=s_nodal';
pp_nodal(:,ss_start.m.dnodes)=ss_start.pp0';

if(ss_start.m.doZ==1), pp=rho_to_p_nd(pp,ss_start.c.b1,ss_start.c.b2,ss_start.c.psc);
    pp_nodal=rho_to_p_nd(pp_nodal,ss_start.c.b1,ss_start.c.b2,ss_start.c.psc); end

%compute mass_start in pipes
if(ss_start.m.doZ==1), p_density=p_to_rho(pp',ss_start.c.b1,ss_start.c.b2,par.ss_start);
else p_density=ss_start.c.psc*pp'/(ss_start.c.gasR*ss_start.c.gasT); end
p_comp=par.ss_start.cc0;
out.ss_start.p_mass=pipe_mass(p_density,p_comp,par.ss_start.m);    %all pipes
out.ss_start.pipe_mass_0=(par.ss_start.n.disc_to_edge*out.ss_start.p_mass)';       %original pipes

ss_start.pnodin=pp_nodal*ss_start.c.psc;   %nodal press_startures (before compress_startion)
ss_start.pnodout=pp*ss_start.c.psc;     %nodal press_startures (after compress_startion)
ss_start.qqin=qq(:,ss_start.n.from_flows);
ss_start.qqout=qq(:,ss_start.n.to_flows);
%pipe inlet and outlet press_starture (compress_startors only at inlets)
for j=1:ss_start.n0.ne
    if(ss_start.n0.comp_bool(j)==1)
        %ss_start.ppin(:,j)=ss_start.pnodout(:,ss_start.n0.from_id(j));
        ss_start.ppin(:,j)=ss_start.pnodin(:,ss_start.n0.from_id(j));
        ss_start.ppout(:,j)=ss_start.pnodin(:,ss_start.n0.to_id(j));
    elseif(ss_start.n0.comp_bool(j)==0)
        ss_start.ppin(:,j)=ss_start.pnodin(:,ss_start.n0.from_id(j));
        ss_start.ppout(:,j)=ss_start.pnodin(:,ss_start.n0.to_id(j));
    end
end
%pipe inlet and outlet flux
ss_start.ffin=ff(:,ss_start.n.from_flows);
ss_start.ffout=ff(:,ss_start.n.to_flows);

out.ss_start.tt0=ss_start.tt0/3600;         %plotting time in hours
out.ss_start.qqinopt=ss_start.qqin;                  %flow boundary in
out.ss_start.qqoutopt=ss_start.qqout;                %flow boundary out
%if(par.out.plotnodal==1)
    out.ss_start.ppoptnodal=ss_start.pnodin(:,1:ss_start.n0.nv);
    out.ss_start.ppopt=ss_start.pnodout(:,1:ss_start.n0.nv);  %press_starture (nodal)
    out.ss_start.qqopt=[out.ss_start.qqinopt out.ss_start.qqoutopt];  %all boundary flows
%else
    %out.ppoptall=ss_start.pnodout(:,1:ss_start.n.nv);   %press_starture (all)
    %out.qqoptall=ss_start.qq0;                      %flows (all)
%end
    out.ss_start.ppinopt=ss_start.ppin; out.ss_start.ppoutopt=ss_start.ppout;
if(par.out.units==1), out.ss_start.ppopt=out.ss_start.ppopt/psi_to_pascal; out.ss_start.ppoptnodal=out.ss_start.ppoptnodal/psi_to_pascal; 
out.ss_start.ppinopt=ss_start.ppin/psi_to_pascal; out.ss_start.ppoutopt=ss_start.ppout/psi_to_pascal; 
out.ss_start.qqopt=out.ss_start.qqopt/mmscfd_to_kgps; out.ss_start.qqinopt=out.ss_start.qqinopt/mmscfd_to_kgps;  
out.ss_start.qqoutopt=out.ss_start.qqoutopt/mmscfd_to_kgps; out.ss_start.pipe_mass_0=out.ss_start.pipe_mass_0/mmscfd_to_kgps/86400; end

%market flow solution
ss_start.m.Yd=[ss_start.m.Yq1(1:ss_start.m.FN) ss_start.m.Yq1(1:ss_start.m.FN)];
% ss_start.m.Yd(ss_start.m.guniqueind,:)=ss_start.m.Yd(ss_start.m.guniqueind,:)+ss_start.m.gtod*ss_start.fd0;
%ss_start.m.Yd=interp1qr(ss_start.m.xd',ss_start.m.Yq(1:ss_start.m.FN,:)',ss_start.tt0)';
%ss_start.m.Yd(ss_start.m.guniqueind,:)=ss_start.m.Yd(ss_start.m.guniqueind,:)+ss_start.m.gtod*ss_start.fd0;
% ss_start.m.Ygd=ss_start.fd0(1:length(ss_start.m.gd),:);
% ss_start.m.Ygs=-ss_start.fd0(length(ss_start.m.gd)+1:length(ss_start.m.gall),:);

%compress_startor discharge press_startures
%if(par.out.dosim==1), out.csetsim=psim1(:,ss_start.m.comp_pos(:,1)); end
out.ss_start.csetopt=out.ss_start.ppopt(:,ss_start.m.comp_pos(:,1));

%process_start parameters (compress_startion ratios and demands)
out.ss_start.cc=ss_start.cc0'; out.ss_start.td=ss_start.m.xd/3600; 
out.ss_start.dbase=ss_start.m.Yq(1:length(ss_start.m.fn),:)'*ss_start.c.qsc;     %base flow "q"
% out.ss_start.gsub=ss_start.m.Yubs'*ss_start.c.qsc;   %upper bounds on sales
% out.ss_start.gslb=ss_start.m.Ylbs'*ss_start.c.qsc;   %lower bounds on sales
% out.ss_start.gdub=ss_start.m.Yubd'*ss_start.c.qsc;   %upper bounds on buys
% out.ss_start.gdlb=ss_start.m.Ylbd'*ss_start.c.qsc;   %lower bounds on buys
% out.ss_start.gdsol=ss_start.m.Ygd'*ss_start.c.qsc;   %gnode buyer solutions
% out.ss_start.gss_startol=-ss_start.m.Ygs'*ss_start.c.qsc;   %gnode seller solutions
% %gnode buyer and seller solutions for all original gnodes 
% GN0=length(ss_start.n0.phys_node);    %number of original gnodes
% out.ss_start.gdsol_all=zeros(2,GN0); out.ss_start.gdsol_all(:,ss_start.dmax_pos)=out.ss_start.gdsol;
% out.ss_start.gss_startol_all=zeros(2,GN0); out.ss_start.gss_startol_all(:,ss_start.smax_pos)=out.ss_start.gss_startol;
out.ss_start.dgflows=full(ss_start.m.Yd(ss_start.m.guniqueind,:))'*ss_start.c.qsc; %flow at nodes with gnodes
slinks=ss_start.m.comp_pos(ss_start.m.spos,2); 
out.ss_start.supp_flow=qq(:,slinks);    %supply flow 
out.ss_start.nonslack_flow=full(ss_start.m.Yd(1:length(ss_start.n0.nonslack_nodes),:))*ss_start.c.qsc;
out.ss_start.flows_all=zeros(ss_start.n0.nv,par.ss_start.m.N1+1);
out.ss_start.flows_all(ss_start.n0.slack_nodes,:)=-out.ss_start.supp_flow;
out.ss_start.flows_all(ss_start.n0.nonslack_nodes,:)=out.ss_start.nonslack_flow;
out.ss_start.dgflows_all=out.ss_start.flows_all(ss_start.n0.phys_node,:); %flow at all original gnodes
% out.ss_start.supp_flow_sim=sim.qq(:,slinks);    %supply flow 
% out.ss_start.flows_all_sim=zeros(sim.n0.nv,length(sim.tt));
% out.ss_start.flows_all_sim(sim.n0.slack_nodes,:)=-out.ss_start.supp_flow_sim;
% out.ss_start.flows_all_sim(sim.n0.nonslack_nodes,:)=full(sim.m.Yd(1:length(sim.n0.nonslack_nodes),:))*sim.c.qsc;
if(par.out.units==1), 
    out.ss_start.dbase=out.ss_start.dbase/mmscfd_to_kgps; 
%     out.ss_start.gsub=out.ss_start.gsub/mmscfd_to_kgps; 
%     out.ss_start.gslb=out.ss_start.gslb/mmscfd_to_kgps; out.ss_start.gdub=out.ss_start.gdub/mmscfd_to_kgps; out.ss_start.gdlb=out.ss_start.gdlb/mmscfd_to_kgps;
%     out.ss_start.gdsol=out.ss_start.gdsol/mmscfd_to_kgps; out.ss_start.gss_startol=out.ss_start.gss_startol/mmscfd_to_kgps;
%     out.ss_start.gdsol_all=out.ss_start.gdsol_all/mmscfd_to_kgps; out.ss_start.gss_startol_all=out.ss_start.gss_startol_all/mmscfd_to_kgps;
    out.ss_start.dgflows=out.ss_start.dgflows/mmscfd_to_kgps; out.ss_start.dgflows_all=out.ss_start.dgflows_all/mmscfd_to_kgps; 
    out.ss_start.dgflows_all=out.ss_start.dgflows_all/mmscfd_to_kgps;  out.ss_start.supp_flow=out.ss_start.supp_flow/mmscfd_to_kgps; 
    out.ss_start.flows_all=out.ss_start.flows_all/mmscfd_to_kgps;  
    %out.ss_start.supp_flow_sim/mmscfd_to_kgps; %out.ss_start.flows_all_sim/mmscfd_to_kgps;  
end

%check nodal flow balance
out.ss_start.flowbal=ss_start.n0.Amp*out.ss_start.qqoutopt'+ss_start.n0.Amm*out.ss_start.qqinopt'-out.ss_start.flows_all;
out.ss_start.flowbalrel=3*out.ss_start.flowbal./(abs(ss_start.n0.Amp*out.ss_start.qqoutopt')+abs(ss_start.n0.Amm*out.ss_start.qqinopt')+abs(out.ss_start.flows_all));
out.ss_start.flowbalrel(mean(out.ss_start.flowbal')./mean(out.ss_start.flowbalrel')<ss_start.m.opt_tol,:)=0; out.ss_start.flowbalrel=out.ss_start.flowbalrel';

%out.ss_start.flowbals=ss_start.n0.Amp*out.ss_start.qqoutsim'+ss_start.n0.Amm*out.ss_start.qqinsim'-out.ss_start.flows_all;
%out.ss_start.flowbalsrel=3*out.ss_start.flowbal./(abs(ss_start.n0.Amp*out.ss_start.qqoutsim')+abs(ss_start.n0.Amm*out.ss_start.qqinsim')+abs(out.ss_start.flows_all));
%out.ss_start.flowbalsrel(mean(out.ss_start.flowbal')./mean(out.ss_start.flowbalrel')<ss_start.m.opt_tol,:)=0;


%compress_startor power
% if(par.out.dosim==1)
%     cposs_startim=par.ss_start.m.comp_pos; m=ss_start.m.mpow;
%     out.ss_start.cccom=interp1qr(out.ss_start.tt0,ss_start.cc0',out.ss_start.ttcom);
%     qcompsim=interp1qr(out.ss_start.tt,sim.qq(:,cposs_startim(:,2)),ttcomsim); 
%     cpow_nd=(abs(qcompsim)).*((out.ss_start.cccom).^(2*m)-1);
%     out.ss_start.cpowsim=cpow_nd.*kron(ss_start.m.eff',ones(size(cpow_nd,1),1))*ss_start.c.mmscfd_to_hp/mmscfd_to_kgps;
% end
out.ss_start.ccopt=ss_start.cc0';
cposopt=par.ss_start.m.comp_pos; m=ss_start.m.mpow;
qcompopt=qq(:,cposopt(:,2)); %cpow_nd=(abs(qcompopt)).*((ss_start.cc0').^(m)-1);
%out.ss_start.cpowopt=cpow_nd.*kron(ss_start.m.eff',ones(size(cpow_nd,1),1))*ss_start.c.mmscfd_to_hp/mmscfd_to_kgps;
out.ss_start.cpowopt=(abs(qcompopt)).*((ss_start.cc0').^(m)-1)*ss_start.m.Wc;   %comp power in Watts

%process_start locational marginal price
% out.ss_start.lmptr=par.ss_start.lmp0(par.ss_start.m.flexnodes,:)'/2*par.ss_start.m.N/par.ss_start.m.odw*par.ss_start.c.Tsc*par.ss_start.c.Tsc/2;
% lmpss_start=par.ss_start.lmp0(par.ss_start.m.flexnodes,:)'/par.ss_start.m.odw*par.ss_start.c.Tsc*par.ss_start.c.Tsc/2;
% out.ss_start.lmptr=par.ss_start.lmp0'/2*par.ss_start.m.N/par.ss_start.m.odw*par.ss_start.c.Tsc*par.ss_start.c.Tsc/2;
% lmpss_start=par.ss_start.lmp0'/ss_start.m.odw*ss_start.c.Tsc*ss_start.c.Tsc/2;
if(ss_start.m.N==0), trmN=2; end
if(ss_start.m.N>1), trmN=ss_start.m.N; end
% out.ss_start.trlmp=ss_start.lmp0'/2*trmN;    %all lmps
% out.ss_start.trlmpnodal=out.ss_start.trlmp(:,1:length(ss_start.m.fn));
% out.ss_start.trlmpnodal_all=zeros(par.ss_start.m.N1+1,ss_start.n0.nv);
% out.ss_start.trlmpnodal_all(:,ss_start.n0.slack_nodes)=-ss_start.m.prslack';
% out.ss_start.trlmpnodal_all(:,ss_start.n0.nonslack_nodes)=out.ss_start.trlmpnodal;
% out.ss_start.gnodelmp=out.ss_start.trlmpnodal_all(:,ss_start.n0.phys_node);
%out.ss_start.gnodelmp=ss_start.lmp0(ss_start.n0.phys_node,:)'/2*trmN;     %lmps at all gnodes
% out.ss_start.dglmp=ss_start.lmp0(ss_start.m.guniqueind,:)'/2*trmN;     %lmps at market nodes
% out.ss_start.gdlmp=ss_start.lmp0(ss_start.m.gallind(1:length(ss_start.m.gd)),:)'/2*trmN;     %lmps at demand gnodes
% out.ss_start.gslmp=ss_start.lmp0(ss_start.m.gallind(length(ss_start.m.gd):length(ss_start.m.gall)),:)'/2*trmN;     %lmps at supply gnodes
% if(par.out.ss_start.dosim==1), out.ss_start.lmpss_start=par.ss_start.lmp0'; end
out.ss_start.Prd=ss_start.m.Prd; out.ss_start.Prs=ss_start.m.Prs; 
out.ss_start.Prslack=interp1qr(ss_start.m.xd',ss_start.m.Prslack',out.ss_start.tt0);  %bid and offer prices
out.ss_start.mult0_pmax=ss_start.mult0_pmax'/2*trmN*ss_start.c.psi_to_pascal/ss_start.c.psc*3600;    %output press_starture marginal prices ($/Psi/hr)
out.ss_start.mult0_cmax=ss_start.mult0_cmax'/2*trmN*3.6/0.75; %compress_startion marginal prices ($/hp)
if(par.out.units==1), 
%     out.ss_start.trlmp=out.ss_start.trlmp*mmscfd_to_kgps; 
%      out.ss_start.trlmpnodal=out.ss_start.trlmpnodal*mmscfd_to_kgps; 
%      out.ss_start.dglmp=out.ss_start.dglmp*mmscfd_to_kgps; out.ss_start.gnodelmp=out.ss_start.gnodelmp*mmscfd_to_kgps;
%      out.ss_start.gdlmp=out.ss_start.gdlmp*mmscfd_to_kgps; out.ss_start.gslmp=out.ss_start.gslmp*mmscfd_to_kgps; 
     out.ss_start.Prd=out.ss_start.Prd*mmscfd_to_kgps; out.ss_start.Prs=out.ss_start.Prs*mmscfd_to_kgps;
     out.ss_start.Prslack=out.ss_start.Prslack*mmscfd_to_kgps;
     out.ss_start.cpowopt=out.ss_start.cpowopt/hp_to_watt;
     %if(par.out.ss_start.dosim==1), out.ss_start.lmpss_start=out.ss_start.lmpss_start*mmscfd_to_kgps; end
end
% out.ss_start.lmpin=zeros(size(out.ss_start.ppinopt)); out.ss_start.lmpout=zeros(size(out.ss_start.ppoutopt));
% for j=1:ss_start.n0.ne
%     if(ismember(ss_start.n0.from_id(j),ss_start.m.pn))
%         out.ss_start.lmpin(:,j)=out.ss_start.Prslack(:,find(ss_start.m.pn==ss_start.n0.from_id(j))); else
%         out.ss_start.lmpin(:,j)=out.ss_start.trlmpnodal(:,find(ss_start.m.fn==ss_start.n0.from_id(j))); end
%     if(ismember(ss_start.n0.to_id(j),ss_start.m.pn))
%         out.ss_start.lmpout(:,j)=out.ss_start.Prslack(:,find(ss_start.m.pn==ss_start.n0.to_id(j))); else
%         out.ss_start.lmpout(:,j)=out.ss_start.trlmpnodal(:,find(ss_start.m.fn==ss_start.n0.to_id(j))); end
% end
%cmap=colormap;
%set(0,'DefaultAxesColorOrder',cmap(floor(rand(length(ss_start.m.gs),1)*64)+1,:))

out.ss_start.guniqueind=ss_start.m.guniqueind; out.ss_start.gunique=ss_start.m.gunique; out.ss_start.fn=ss_start.m.fn; out.ss_start.pn=ss_start.m.pn;
out.ss_start.n0=ss_start.n0; out.ss_start.n=ss_start.n; out.ss_start.gd=ss_start.m.gd; out.ss_start.gs=ss_start.m.gs; out.ss_start.FN=ss_start.m.FN; out.ss_start.PN=ss_start.m.PN;
out.ss_start.cn=ss_start.m.C; 
out.ss_start.mfolder=par.mfolder;

if(par.out.savecsvoutput==1)
        mfolder=par.mfolder;
%     if(par.out.intervals_out==0)
        pipe_cols=[1:out.ss_start.n0.ne-out.ss_start.n0.nc]; comp_cols=[out.ss_start.n0.ne-out.ss_start.n0.nc+1:out.ss_start.n0.ne];
        %dlmwrite([mfolder '\output_ss_start_tpts.csv'],double(out.ss_start.tt0),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_pipe-press_starture-in.csv'],double([pipe_cols;out.ss_start.ppinopt(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_pipe-press_starture-out.csv'],double([pipe_cols;out.ss_start.ppoutopt(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_comp-press_starture-in.csv'],double([1:out.ss_start.n0.nc;out.ss_start.ppinopt(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_comp-press_starture-out.csv'],double([1:out.ss_start.n0.nc;out.ss_start.ppoutopt(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_pipe-flow-in.csv'],double([pipe_cols;out.ss_start.qqinopt(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_pipe-flow-out.csv'],double([pipe_cols;out.ss_start.qqoutopt(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_comp-flow-in.csv'],double([1:out.ss_start.n0.nc;out.ss_start.qqinopt(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_comp-flow-out.csv'],double([1:out.ss_start.n0.nc;out.ss_start.qqoutopt(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_nodal-press_starture.csv'],double([[1:out.ss_start.n0.nv];out.ss_start.ppoptnodal(1,:)]),'precision',16,'delimiter',',');
        %dlmwrite([mfolder '\output_ss_start_gnode-physical-withdrawals.csv'],double([out.ss_start.gunique';out.ss_start.dgflows]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_nonslack-flows.csv'],double([out.ss_start.fn';out.ss_start.flows_all(ss_start.n0.nonslack_nodes,1)']),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_start_gnode-supply-flows.csv'],double([1:GN0;out.ss_start.gss_startol_all(1,:)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_start_gnode-demand-flows.csv'],double([1:GN0;out.ss_start.gdsol_all(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_slack-flows.csv'],double([out.ss_start.pn';out.ss_start.supp_flow(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_comp-ratios.csv'],double([[1:out.ss_start.cn];out.ss_start.ccopt(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_comp-discharge-press_starture.csv'],double([[1:out.ss_start.cn];out.ss_start.csetopt(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_comp-power.csv'],double([[1:out.ss_start.cn];out.ss_start.cpowopt(1,:)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_start_lmp-nodal-all.csv'],double([out.ss_start.fn';out.ss_start.trlmpnodal(1,:)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_start_lmp-gnodes.csv'],double([[1:GN0];out.ss_start.gnodelmp(1,:)]),'precision',16,'delimiter',',');
        %dlmwrite([mfolder '\output_ss_start_lmp-bidders.csv'],double([out.ss_start.gunique';out.ss_start.dglmp]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_start_pipe-lmp-in.csv'],double([pipe_cols;out.ss_start.lmpin(1,pipe_cols)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_start_pipe-lmp-out.csv'],double([pipe_cols;out.ss_start.lmpout(1,pipe_cols)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_start_comp-lmp-in.csv'],double([1:out.ss_start.n0.nc;out.ss_start.lmpin(1,comp_cols)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_start_comp-lmp-out.csv'],double([1:out.ss_start.n0.nc;out.ss_start.lmpout(1,comp_cols)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_start_comp-pmax-mp.csv'],double([[1:out.ss_start.n0.nc];out.ss_start.mult0_pmax(1,:)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_start_comp-hpmax-mp.csv'],double([[1:out.ss_start.n0.nc];out.ss_start.mult0_cmax(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_flowbalrel.csv'],double([[1:out.ss_start.n0.nv];out.ss_start.flowbalrel(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_start_pipe-mass.csv'],double([pipe_cols;out.ss_start.pipe_mass_0(1,pipe_cols)]),'precision',16,'delimiter',',');
end
    if(par.out.intervals_out>0)
        %inputs on intervals
        out.ss_start.int_qbar=ss_start.int_qbar; out.ss_start.int_dmax=ss_start.int_dmax; out.ss_start.int_dmin=ss_start.int_dmin;
        out.ss_start.int_smax=ss_start.int_smax; out.ss_start.int_smin=ss_start.int_smin; out.ss_start.int_cd=ss_start.int_cd;
        out.ss_start.int_cs=ss_start.int_cs; out.ss_start.int_cslack=ss_start.int_cslack; out.ss_start.int_pslack=ss_start.int_pslack;
        if(ss_start.intervals>0 && ss_start.units==1)
            out.ss_start.int_cd=out.ss_start.int_cd*mmscfd_to_kgps;
            out.ss_start.int_cs=out.ss_start.int_cs*mmscfd_to_kgps; 
            out.ss_start.int_cslack=out.ss_start.int_cslack*mmscfd_to_kgps;
        end
        %------------------
        %revert to full gnode index list
        out.ss_start.gdsol_all=zeros(2,ss_start.n0.ng);
        out.ss_start.gdsol_all(:,ss_start.dmax_pos)=out.ss_start.gdsol;
        out.ss_start.gss_startol_all=zeros(2,ss_start.n0.ng);
        out.ss_start.gss_startol_all(:,ss_start.smax_pos)=out.ss_start.gss_startol;
        
        In=par.out.intervals_out;  %24 intervals on optimization period (1 day)
        T=ss_start.c.T/3600;             %time in hours
        int_bounds=[0:T/In:T]; out.ss_start.int_bounds=int_bounds;
        [out.ss_start.dbase_int]=pts_to_int(ss_start.m.xd'/3600,out.ss_start.dbase,int_bounds');
        [out.ss_start.gsub_int]=pts_to_int(ss_start.m.xd'/3600,out.ss_start.gsub,int_bounds');
        [out.ss_start.gslb_int]=pts_to_int(ss_start.m.xd'/3600,out.ss_start.gslb,int_bounds');
        [out.ss_start.gdub_int]=pts_to_int(ss_start.m.xd'/3600,out.ss_start.gdub,int_bounds');
        [out.ss_start.gdlb_int]=pts_to_int(ss_start.m.xd'/3600,out.ss_start.gdlb,int_bounds');
        [out.ss_start.gdsol_int]=pts_to_int(out.ss_start.tt0,out.ss_start.gdsol_all,int_bounds');
        [out.ss_start.gss_startol_int]=pts_to_int(out.ss_start.tt0,out.ss_start.gss_startol_all,int_bounds');
        [out.ss_start.dgflows_int]=pts_to_int(out.ss_start.tt0,out.ss_start.dgflows_all,int_bounds');
        [out.ss_start.supp_flow_int]=pts_to_int(out.ss_start.tt0,out.ss_start.supp_flow,int_bounds');
        [out.ss_start.flows_all_int]=pts_to_int(out.ss_start.tt0,out.ss_start.flows_all',int_bounds');
        [out.ss_start.Prslack_int]=pts_to_int(ss_start.m.xd'/3600,ss_start.m.Prslack',int_bounds');
        [out.ss_start.Prs_int]=pts_to_int(ss_start.m.xd'/3600,ss_start.m.Prs',int_bounds');
        [out.ss_start.Prd_int]=pts_to_int(ss_start.m.xd'/3600,ss_start.m.Prd',int_bounds');
        %------------------
        [out.ss_start.ppinopt_int]=pts_to_int(out.ss_start.tt0,out.ss_start.ppinopt,int_bounds');
        [out.ss_start.ppoutopt_int]=pts_to_int(out.ss_start.tt0,out.ss_start.ppoutopt,int_bounds');
        [out.ss_start.qqinopt_int]=pts_to_int(out.ss_start.tt0,out.ss_start.qqinopt,int_bounds');
        [out.ss_start.qqoutopt_int]=pts_to_int(out.ss_start.tt0,out.ss_start.qqoutopt,int_bounds');
        [out.ss_start.ppoptnodal_int]=pts_to_int(out.ss_start.tt0,out.ss_start.ppoptnodal,int_bounds');
        [out.ss_start.dgflows_int]=pts_to_int(out.ss_start.tt0,out.ss_start.dgflows,int_bounds');
        [out.ss_start.supp_flow_int]=pts_to_int(out.ss_start.tt0,out.ss_start.supp_flow,int_bounds');
        [out.ss_start.ccopt_int]=pts_to_int(out.ss_start.tt0,out.ss_start.ccopt,int_bounds');
        [out.ss_start.csetopt_int]=pts_to_int(out.ss_start.tt0,out.ss_start.csetopt,int_bounds');
        [out.ss_start.cpowopt_int]=pts_to_int(out.ss_start.tt0,out.ss_start.cpowopt,int_bounds');
        [out.ss_start.trlmpnodal_int]=pts_to_int(out.ss_start.tt0,out.ss_start.trlmpnodal,int_bounds');
        [out.ss_start.dglmp_int]=pts_to_int(out.ss_start.tt0,out.ss_start.dglmp,int_bounds');
        [out.ss_start.lmpin_int]=pts_to_int(out.ss_start.tt0,out.ss_start.lmpin,int_bounds');
        [out.ss_start.lmpout_int]=pts_to_int(out.ss_start.tt0,out.ss_start.lmpout,int_bounds');
        [out.ss_start.mult0_pmax_int]=pts_to_int(out.ss_start.tt0,out.ss_start.mult0_pmax,int_bounds');
        [out.ss_start.mult0_cmax_int]=pts_to_int(out.ss_start.tt0,out.ss_start.mult0_cmax,int_bounds');
        [out.ss_start.flowbalrel_int]=pts_to_int(out.ss_start.tt0,out.ss_start.flowbalrel,int_bounds');
        [out.ss_start.pipe_mass_int]=pts_to_int(out.ss_start.tt0,out.ss_start.pipe_mass_0,int_bounds');
        pipe_cols=[1:out.ss_start.n0.ne-out.ss_start.n0.nc]; comp_cols=[out.ss_start.n0.ne-out.ss_start.n0.nc+1:out.ss_start.n0.ne];
        

        if(par.out.steadystateonly==0)
            %write files
            dlmwrite([mfolder '\output_int_pipe-press_starture-in.csv'],double([pipe_cols;out.ss_start.ppinopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-press_starture-out.ss_start.csv'],double([pipe_cols;out.ss_start.ppoutopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-press_starture-in.csv'],double([1:out.ss_start.n0.nc;out.ss_start.ppinopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-press_starture-out.ss_start.csv'],double([1:out.ss_start.n0.nc;out.ss_start.ppoutopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-flow-in.csv'],double([pipe_cols;out.ss_start.qqinopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-flow-out.ss_start.csv'],double([pipe_cols;out.ss_start.qqoutopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-flow-in.csv'],double([1:out.ss_start.n0.nc;out.ss_start.qqinopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-flow-out.ss_start.csv'],double([1:out.ss_start.n0.nc;out.ss_start.qqoutopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_nodal-press_starture.csv'],double([[1:out.ss_start.n0.nv];out.ss_start.ppoptnodal_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_gnode-physical-withdrawals.csv'],double([ss_start.m.guniqueind';out.ss_start.dgflows_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_gnode-supply-flows.csv'],double([ss_start.n0.phys_node';out.ss_start.gss_startol_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_gnode-demand-flows.csv'],double([ss_start.n0.phys_node';out.ss_start.gdsol_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_slack-flows.csv'],double([out.ss_start.pn';out.ss_start.supp_flow_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-ratios.csv'],double([[1:out.ss_start.cn];out.ss_start.ccopt_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-discharge-press_starture.csv'],double([[1:out.ss_start.cn];out.ss_start.csetopt_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-power.csv'],double([[1:out.ss_start.cn];out.ss_start.cpowopt_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_lmp-nodal-all.csv'],double([out.ss_start.fn';out.ss_start.trlmpnodal_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_lmp-bidders.csv'],douint_flow_constble([out.ss_start.gunique';out.ss_start.dglmp_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-lmp-in.csv'],double([pipe_cols;out.ss_start.lmpin_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-lmp-out.ss_start.csv'],double([pipe_cols;out.ss_start.lmpout_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-lmp-in.csv'],double([1:out.ss_start.n0.nc;out.ss_start.lmpin_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-lmp-out.ss_start.csv'],double([1:out.ss_start.n0.nc;out.ss_start.lmpout_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-pmax-mp.csv'],double([[1:out.ss_start.n0.nc];out.ss_start.mult0_pmax_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-hpmax-mp.csv'],double([[1:out.ss_start.n0.nc];out.ss_start.mult0_cmax_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_flowbalrel.csv'],double([[1:out.ss_start.n0.nv];out.ss_start.flowbalrel_int]),'precision',16,'delimiter',',');
        end
    end

if(ss_start.m.save_state==1)
dlmwrite([mfolder '\output_ss_start_state_save.csv'],double(full(ss_start.state_save)),'precision',16,'delimiter',',');
end

par.out.ss_start=out.ss_start;
%par.ss_start=ss_start;

%out.ss_start.mult0_pmax=ss_start.mult0_pmax/2*ss_start.m.N/(ss_start.c.psc/1000000)/mpa_to_psi;    %output press_starture marginal prices ($/
%out.ss_start.mult0_cmax=ss_start.mult0_cmax/2*ss_start.m.N*3.6/0.75; %compress_startion marginal prices ($/hp)

function [xints]=pts_to_int(tpts,xpts,ibnds)
    xbnds=interp1qr(tpts,xpts,ibnds); In=length(ibnds)-1;
    xints=(xbnds(1:In,:)+xbnds(2:In+1,:))/2;
return;