function [par]=process_output_ss_terminal_nofd(par)
% Anatoly Zlotnik, February 2019

% mfolder=par.mfolder;
out=par.out;
% 
ss_terminal=par.ss_terminal; 

psi_to_pascal=ss_terminal.c.psi_to_pascal;
mpa_to_psi=1000000/psi_to_pascal;
ss_terminal.c.mpa_to_psi=mpa_to_psi;
mmscfd_to_kgps=ss_terminal.c.mmscfd_to_kgps;
hp_to_watt=745.7;
if(par.out.doZ==1), b1=ss_terminal.c.b1; b2=ss_terminal.c.b2; end

%process_terminal optimization output
%process_terminal dimensional solution (metric)
pp=zeros(length(ss_terminal.tt0),ss_terminal.n0.nv);
%s_nodal=[ss_terminal.m.pslack ss_terminal.m.pslack(:,1)];
s_nodal=ss_terminal.m.pslack;
if(par.ss_terminal.m.extension>0), s_nodal=ss_terminal.m.pslack; end
scf=find(ismember(ss_terminal.m.snodes,ss_terminal.m.comp_pos(ss_terminal.m.spos,1)));
s(scf,:)=ss_terminal.cc0(ss_terminal.m.spos,:).*s_nodal(scf,:);
pp(:,ss_terminal.m.snodes)=s';
pp(:,ss_terminal.m.dnodes)=ss_terminal.pp0';
qq=ss_terminal.qq0'*ss_terminal.m.Xs*ss_terminal.c.qsc; 
ff=ss_terminal.qq0'*ss_terminal.c.qsc; 
for j=1:length(ss_terminal.m.dpos)
    cpj=ss_terminal.n.comp_pos(ss_terminal.m.dpos(j),1);
    pp(:,cpj)=pp(:,cpj).*ss_terminal.cc0(ss_terminal.m.dpos(j),:)';
end
%nodal press_terminalure (before compress_terminalors)
pp_nodal=zeros(length(ss_terminal.tt0),ss_terminal.n.nv);
pp_nodal(:,ss_terminal.m.snodes)=s_nodal';
pp_nodal(:,ss_terminal.m.dnodes)=ss_terminal.pp0';

if(ss_terminal.m.doZ==1), pp=rho_to_p_nd(pp,ss_terminal.c.b1,ss_terminal.c.b2,ss_terminal.c.psc);
    pp_nodal=rho_to_p_nd(pp_nodal,ss_terminal.c.b1,ss_terminal.c.b2,ss_terminal.c.psc); end

%compute mass_terminal in pipes
if(ss_terminal.m.doZ==1), p_density=p_to_rho(pp',ss_terminal.c.b1,ss_terminal.c.b2,par.ss_terminal);
else p_density=ss_terminal.c.psc*pp'/(ss_terminal.c.gasR*ss_terminal.c.gasT); end
p_comp=par.ss_terminal.cc0;
out.ss_terminal.p_mass=pipe_mass(p_density,p_comp,par.ss_terminal.m);    %all pipes
out.ss_terminal.pipe_mass_0=(par.ss_terminal.n.disc_to_edge*out.ss_terminal.p_mass)';       %original pipes

ss_terminal.pnodin=pp_nodal*ss_terminal.c.psc;   %nodal press_terminalures (before compress_terminalion)
ss_terminal.pnodout=pp*ss_terminal.c.psc;     %nodal press_terminalures (after compress_terminalion)
ss_terminal.qqin=qq(:,ss_terminal.n.from_flows);
ss_terminal.qqout=qq(:,ss_terminal.n.to_flows);
%pipe inlet and outlet press_terminalure (compress_terminalors only at inlets)
for j=1:ss_terminal.n0.ne
    if(ss_terminal.n0.comp_bool(j)==1)
        %ss_terminal.ppin(:,j)=ss_terminal.pnodout(:,ss_terminal.n0.from_id(j));
        ss_terminal.ppin(:,j)=ss_terminal.pnodin(:,ss_terminal.n0.from_id(j));
        ss_terminal.ppout(:,j)=ss_terminal.pnodin(:,ss_terminal.n0.to_id(j));
    elseif(ss_terminal.n0.comp_bool(j)==0)
        ss_terminal.ppin(:,j)=ss_terminal.pnodin(:,ss_terminal.n0.from_id(j));
        ss_terminal.ppout(:,j)=ss_terminal.pnodin(:,ss_terminal.n0.to_id(j));
    end
end
%pipe inlet and outlet flux
ss_terminal.ffin=ff(:,ss_terminal.n.from_flows);
ss_terminal.ffout=ff(:,ss_terminal.n.to_flows);

out.ss_terminal.tt0=ss_terminal.tt0/3600;         %plotting time in hours
out.ss_terminal.qqinopt=ss_terminal.qqin;                  %flow boundary in
out.ss_terminal.qqoutopt=ss_terminal.qqout;                %flow boundary out
%if(par.out.plotnodal==1)
    out.ss_terminal.ppoptnodal=ss_terminal.pnodin(:,1:ss_terminal.n0.nv);
    out.ss_terminal.ppopt=ss_terminal.pnodout(:,1:ss_terminal.n0.nv);  %press_terminalure (nodal)
    out.ss_terminal.qqopt=[out.ss_terminal.qqinopt out.ss_terminal.qqoutopt];  %all boundary flows
%else
    %out.ppoptall=ss_terminal.pnodout(:,1:ss_terminal.n.nv);   %press_terminalure (all)
    %out.qqoptall=ss_terminal.qq0;                      %flows (all)
%end
    out.ss_terminal.ppinopt=ss_terminal.ppin; out.ss_terminal.ppoutopt=ss_terminal.ppout;
if(par.out.units==1), out.ss_terminal.ppopt=out.ss_terminal.ppopt/psi_to_pascal; out.ss_terminal.ppoptnodal=out.ss_terminal.ppoptnodal/psi_to_pascal; 
out.ss_terminal.ppinopt=ss_terminal.ppin/psi_to_pascal; out.ss_terminal.ppoutopt=ss_terminal.ppout/psi_to_pascal; 
out.ss_terminal.qqopt=out.ss_terminal.qqopt/mmscfd_to_kgps; out.ss_terminal.qqinopt=out.ss_terminal.qqinopt/mmscfd_to_kgps;  
out.ss_terminal.qqoutopt=out.ss_terminal.qqoutopt/mmscfd_to_kgps; out.ss_terminal.pipe_mass_0=out.ss_terminal.pipe_mass_0/mmscfd_to_kgps/86400; end

%market flow solution
ss_terminal.m.Yd=[ss_terminal.m.Yq1(1:ss_terminal.m.FN) ss_terminal.m.Yq1(1:ss_terminal.m.FN)];
% ss_terminal.m.Yd(ss_terminal.m.guniqueind,:)=ss_terminal.m.Yd(ss_terminal.m.guniqueind,:)+ss_terminal.m.gtod*ss_terminal.fd0;
%ss_terminal.m.Yd=interp1qr(ss_terminal.m.xd',ss_terminal.m.Yq(1:ss_terminal.m.FN,:)',ss_terminal.tt0)';
%ss_terminal.m.Yd(ss_terminal.m.guniqueind,:)=ss_terminal.m.Yd(ss_terminal.m.guniqueind,:)+ss_terminal.m.gtod*ss_terminal.fd0;
% ss_terminal.m.Ygd=ss_terminal.fd0(1:length(ss_terminal.m.gd),:);
% ss_terminal.m.Ygs=-ss_terminal.fd0(length(ss_terminal.m.gd)+1:length(ss_terminal.m.gall),:);

%compress_terminalor discharge press_terminalures
%if(par.out.dosim==1), out.csetsim=psim1(:,ss_terminal.m.comp_pos(:,1)); end
out.ss_terminal.csetopt=out.ss_terminal.ppopt(:,ss_terminal.m.comp_pos(:,1));

%process_terminal parameters (compress_terminalion ratios and demands)
out.ss_terminal.cc=ss_terminal.cc0'; out.ss_terminal.td=ss_terminal.m.xd/3600; 
out.ss_terminal.dbase=ss_terminal.m.Yq(1:length(ss_terminal.m.fn),:)'*ss_terminal.c.qsc;     %base flow "q"
% out.ss_terminal.gsub=ss_terminal.m.Yubs'*ss_terminal.c.qsc;   %upper bounds on sales
% out.ss_terminal.gslb=ss_terminal.m.Ylbs'*ss_terminal.c.qsc;   %lower bounds on sales
% out.ss_terminal.gdub=ss_terminal.m.Yubd'*ss_terminal.c.qsc;   %upper bounds on buys
% out.ss_terminal.gdlb=ss_terminal.m.Ylbd'*ss_terminal.c.qsc;   %lower bounds on buys
% out.ss_terminal.gdsol=ss_terminal.m.Ygd'*ss_terminal.c.qsc;   %gnode buyer solutions
% out.ss_terminal.gss_terminalol=-ss_terminal.m.Ygs'*ss_terminal.c.qsc;   %gnode seller solutions
% %gnode buyer and seller solutions for all original gnodes 
% GN0=length(ss_terminal.n0.phys_node);    %number of original gnodes
% out.ss_terminal.gdsol_all=zeros(2,GN0); out.ss_terminal.gdsol_all(:,ss_terminal.dmax_pos)=out.ss_terminal.gdsol;
% out.ss_terminal.gss_terminalol_all=zeros(2,GN0); out.ss_terminal.gss_terminalol_all(:,ss_terminal.smax_pos)=out.ss_terminal.gss_terminalol;
out.ss_terminal.dgflows=full(ss_terminal.m.Yd(ss_terminal.m.guniqueind,:))'*ss_terminal.c.qsc; %flow at nodes with gnodes
slinks=ss_terminal.m.comp_pos(ss_terminal.m.spos,2); 
out.ss_terminal.supp_flow=qq(:,slinks);    %supply flow 
out.ss_terminal.nonslack_flow=full(ss_terminal.m.Yd(1:length(ss_terminal.n0.nonslack_nodes),:))*ss_terminal.c.qsc;
out.ss_terminal.flows_all=zeros(ss_terminal.n0.nv,par.ss_terminal.m.N1+1);
out.ss_terminal.flows_all(ss_terminal.n0.slack_nodes,:)=-out.ss_terminal.supp_flow;
out.ss_terminal.flows_all(ss_terminal.n0.nonslack_nodes,:)=out.ss_terminal.nonslack_flow;
out.ss_terminal.dgflows_all=out.ss_terminal.flows_all(ss_terminal.n0.phys_node,:); %flow at all original gnodes
% out.ss_terminal.supp_flow_sim=sim.qq(:,slinks);    %supply flow 
% out.ss_terminal.flows_all_sim=zeros(sim.n0.nv,length(sim.tt));
% out.ss_terminal.flows_all_sim(sim.n0.slack_nodes,:)=-out.ss_terminal.supp_flow_sim;
% out.ss_terminal.flows_all_sim(sim.n0.nonslack_nodes,:)=full(sim.m.Yd(1:length(sim.n0.nonslack_nodes),:))*sim.c.qsc;
if(par.out.units==1), 
    out.ss_terminal.dbase=out.ss_terminal.dbase/mmscfd_to_kgps; 
%     out.ss_terminal.gsub=out.ss_terminal.gsub/mmscfd_to_kgps; 
%     out.ss_terminal.gslb=out.ss_terminal.gslb/mmscfd_to_kgps; out.ss_terminal.gdub=out.ss_terminal.gdub/mmscfd_to_kgps; out.ss_terminal.gdlb=out.ss_terminal.gdlb/mmscfd_to_kgps;
%     out.ss_terminal.gdsol=out.ss_terminal.gdsol/mmscfd_to_kgps; out.ss_terminal.gss_terminalol=out.ss_terminal.gss_terminalol/mmscfd_to_kgps;
%     out.ss_terminal.gdsol_all=out.ss_terminal.gdsol_all/mmscfd_to_kgps; out.ss_terminal.gss_terminalol_all=out.ss_terminal.gss_terminalol_all/mmscfd_to_kgps;
    out.ss_terminal.dgflows=out.ss_terminal.dgflows/mmscfd_to_kgps; out.ss_terminal.dgflows_all=out.ss_terminal.dgflows_all/mmscfd_to_kgps; 
    out.ss_terminal.dgflows_all=out.ss_terminal.dgflows_all/mmscfd_to_kgps;  out.ss_terminal.supp_flow=out.ss_terminal.supp_flow/mmscfd_to_kgps; 
    out.ss_terminal.flows_all=out.ss_terminal.flows_all/mmscfd_to_kgps;  
    %out.ss_terminal.supp_flow_sim/mmscfd_to_kgps; %out.ss_terminal.flows_all_sim/mmscfd_to_kgps;  
end

%check nodal flow balance
out.ss_terminal.flowbal=ss_terminal.n0.Amp*out.ss_terminal.qqoutopt'+ss_terminal.n0.Amm*out.ss_terminal.qqinopt'-out.ss_terminal.flows_all;
out.ss_terminal.flowbalrel=3*out.ss_terminal.flowbal./(abs(ss_terminal.n0.Amp*out.ss_terminal.qqoutopt')+abs(ss_terminal.n0.Amm*out.ss_terminal.qqinopt')+abs(out.ss_terminal.flows_all));
out.ss_terminal.flowbalrel(mean(out.ss_terminal.flowbal')./mean(out.ss_terminal.flowbalrel')<ss_terminal.m.opt_tol,:)=0; out.ss_terminal.flowbalrel=out.ss_terminal.flowbalrel';

%out.ss_terminal.flowbals=ss_terminal.n0.Amp*out.ss_terminal.qqoutsim'+ss_terminal.n0.Amm*out.ss_terminal.qqinsim'-out.ss_terminal.flows_all;
%out.ss_terminal.flowbalsrel=3*out.ss_terminal.flowbal./(abs(ss_terminal.n0.Amp*out.ss_terminal.qqoutsim')+abs(ss_terminal.n0.Amm*out.ss_terminal.qqinsim')+abs(out.ss_terminal.flows_all));
%out.ss_terminal.flowbalsrel(mean(out.ss_terminal.flowbal')./mean(out.ss_terminal.flowbalrel')<ss_terminal.m.opt_tol,:)=0;


%compress_terminalor power
% if(par.out.dosim==1)
%     cposs_terminalim=par.ss_terminal.m.comp_pos; m=ss_terminal.m.mpow;
%     out.ss_terminal.cccom=interp1qr(out.ss_terminal.tt0,ss_terminal.cc0',out.ss_terminal.ttcom);
%     qcompsim=interp1qr(out.ss_terminal.tt,sim.qq(:,cposs_terminalim(:,2)),ttcomsim); 
%     cpow_nd=(abs(qcompsim)).*((out.ss_terminal.cccom).^(2*m)-1);
%     out.ss_terminal.cpowsim=cpow_nd.*kron(ss_terminal.m.eff',ones(size(cpow_nd,1),1))*ss_terminal.c.mmscfd_to_hp/mmscfd_to_kgps;
% end
out.ss_terminal.ccopt=ss_terminal.cc0';
cposopt=par.ss_terminal.m.comp_pos; m=ss_terminal.m.mpow;
qcompopt=qq(:,cposopt(:,2)); %cpow_nd=(abs(qcompopt)).*((ss_terminal.cc0').^(m)-1);
%out.ss_terminal.cpowopt=cpow_nd.*kron(ss_terminal.m.eff',ones(size(cpow_nd,1),1))*ss_terminal.c.mmscfd_to_hp/mmscfd_to_kgps;
out.ss_terminal.cpowopt=(abs(qcompopt)).*((ss_terminal.cc0').^(m)-1)*ss_terminal.m.Wc;   %comp power in Watts

%process_terminal locational marginal price
% out.ss_terminal.lmptr=par.ss_terminal.lmp0(par.ss_terminal.m.flexnodes,:)'/2*par.ss_terminal.m.N/par.ss_terminal.m.odw*par.ss_terminal.c.Tsc*par.ss_terminal.c.Tsc/2;
% lmpss_terminal=par.ss_terminal.lmp0(par.ss_terminal.m.flexnodes,:)'/par.ss_terminal.m.odw*par.ss_terminal.c.Tsc*par.ss_terminal.c.Tsc/2;
% out.ss_terminal.lmptr=par.ss_terminal.lmp0'/2*par.ss_terminal.m.N/par.ss_terminal.m.odw*par.ss_terminal.c.Tsc*par.ss_terminal.c.Tsc/2;
% lmpss_terminal=par.ss_terminal.lmp0'/ss_terminal.m.odw*ss_terminal.c.Tsc*ss_terminal.c.Tsc/2;
if(ss_terminal.m.N==0), trmN=2; end
if(ss_terminal.m.N>1), trmN=ss_terminal.m.N; end
% out.ss_terminal.trlmp=ss_terminal.lmp0'/2*trmN;    %all lmps
% out.ss_terminal.trlmpnodal=out.ss_terminal.trlmp(:,1:length(ss_terminal.m.fn));
% out.ss_terminal.trlmpnodal_all=zeros(par.ss_terminal.m.N1+1,ss_terminal.n0.nv);
% out.ss_terminal.trlmpnodal_all(:,ss_terminal.n0.slack_nodes)=-ss_terminal.m.prslack';
% out.ss_terminal.trlmpnodal_all(:,ss_terminal.n0.nonslack_nodes)=out.ss_terminal.trlmpnodal;
% out.ss_terminal.gnodelmp=out.ss_terminal.trlmpnodal_all(:,ss_terminal.n0.phys_node);
%out.ss_terminal.gnodelmp=ss_terminal.lmp0(ss_terminal.n0.phys_node,:)'/2*trmN;     %lmps at all gnodes
% out.ss_terminal.dglmp=ss_terminal.lmp0(ss_terminal.m.guniqueind,:)'/2*trmN;     %lmps at market nodes
% out.ss_terminal.gdlmp=ss_terminal.lmp0(ss_terminal.m.gallind(1:length(ss_terminal.m.gd)),:)'/2*trmN;     %lmps at demand gnodes
% out.ss_terminal.gslmp=ss_terminal.lmp0(ss_terminal.m.gallind(length(ss_terminal.m.gd):length(ss_terminal.m.gall)),:)'/2*trmN;     %lmps at supply gnodes
% if(par.out.ss_terminal.dosim==1), out.ss_terminal.lmpss_terminal=par.ss_terminal.lmp0'; end
out.ss_terminal.Prd=ss_terminal.m.Prd; out.ss_terminal.Prs=ss_terminal.m.Prs; 
out.ss_terminal.Prslack=interp1qr(ss_terminal.m.xd',ss_terminal.m.Prslack',out.ss_terminal.tt0);  %bid and offer prices
out.ss_terminal.mult0_pmax=ss_terminal.mult0_pmax'/2*trmN*ss_terminal.c.psi_to_pascal/ss_terminal.c.psc*3600;    %output press_terminalure marginal prices ($/Psi/hr)
out.ss_terminal.mult0_cmax=ss_terminal.mult0_cmax'/2*trmN*3.6/0.75; %compress_terminalion marginal prices ($/hp)
if(par.out.units==1), 
%     out.ss_terminal.trlmp=out.ss_terminal.trlmp*mmscfd_to_kgps; 
%      out.ss_terminal.trlmpnodal=out.ss_terminal.trlmpnodal*mmscfd_to_kgps; 
%      out.ss_terminal.dglmp=out.ss_terminal.dglmp*mmscfd_to_kgps; out.ss_terminal.gnodelmp=out.ss_terminal.gnodelmp*mmscfd_to_kgps;
%      out.ss_terminal.gdlmp=out.ss_terminal.gdlmp*mmscfd_to_kgps; out.ss_terminal.gslmp=out.ss_terminal.gslmp*mmscfd_to_kgps; 
     out.ss_terminal.Prd=out.ss_terminal.Prd*mmscfd_to_kgps; out.ss_terminal.Prs=out.ss_terminal.Prs*mmscfd_to_kgps;
     out.ss_terminal.Prslack=out.ss_terminal.Prslack*mmscfd_to_kgps;
     out.ss_terminal.cpowopt=out.ss_terminal.cpowopt/hp_to_watt;
     %if(par.out.ss_terminal.dosim==1), out.ss_terminal.lmpss_terminal=out.ss_terminal.lmpss_terminal*mmscfd_to_kgps; end
end
% out.ss_terminal.lmpin=zeros(size(out.ss_terminal.ppinopt)); out.ss_terminal.lmpout=zeros(size(out.ss_terminal.ppoutopt));
% for j=1:ss_terminal.n0.ne
%     if(ismember(ss_terminal.n0.from_id(j),ss_terminal.m.pn))
%         out.ss_terminal.lmpin(:,j)=out.ss_terminal.Prslack(:,find(ss_terminal.m.pn==ss_terminal.n0.from_id(j))); else
%         out.ss_terminal.lmpin(:,j)=out.ss_terminal.trlmpnodal(:,find(ss_terminal.m.fn==ss_terminal.n0.from_id(j))); end
%     if(ismember(ss_terminal.n0.to_id(j),ss_terminal.m.pn))
%         out.ss_terminal.lmpout(:,j)=out.ss_terminal.Prslack(:,find(ss_terminal.m.pn==ss_terminal.n0.to_id(j))); else
%         out.ss_terminal.lmpout(:,j)=out.ss_terminal.trlmpnodal(:,find(ss_terminal.m.fn==ss_terminal.n0.to_id(j))); end
% end
%cmap=colormap;
%set(0,'DefaultAxesColorOrder',cmap(floor(rand(length(ss_terminal.m.gs),1)*64)+1,:))

out.ss_terminal.guniqueind=ss_terminal.m.guniqueind; out.ss_terminal.gunique=ss_terminal.m.gunique; out.ss_terminal.fn=ss_terminal.m.fn; out.ss_terminal.pn=ss_terminal.m.pn;
out.ss_terminal.n0=ss_terminal.n0; out.ss_terminal.n=ss_terminal.n; out.ss_terminal.gd=ss_terminal.m.gd; out.ss_terminal.gs=ss_terminal.m.gs; out.ss_terminal.FN=ss_terminal.m.FN; out.ss_terminal.PN=ss_terminal.m.PN;
out.ss_terminal.cn=ss_terminal.m.C; 
out.ss_terminal.mfolder=par.mfolder;

if(par.out.savecsvoutput==1)
        mfolder=par.mfolder;
%     if(par.out.intervals_out==0)
        pipe_cols=[1:out.ss_terminal.n0.ne-out.ss_terminal.n0.nc]; comp_cols=[out.ss_terminal.n0.ne-out.ss_terminal.n0.nc+1:out.ss_terminal.n0.ne];
        %dlmwrite([mfolder '\output_ss_terminal_tpts.csv'],double(out.ss_terminal.tt0),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_pipe-press_terminalure-in.csv'],double([pipe_cols;out.ss_terminal.ppinopt(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_pipe-press_terminalure-out.csv'],double([pipe_cols;out.ss_terminal.ppoutopt(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_comp-press_terminalure-in.csv'],double([1:out.ss_terminal.n0.nc;out.ss_terminal.ppinopt(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_comp-press_terminalure-out.csv'],double([1:out.ss_terminal.n0.nc;out.ss_terminal.ppoutopt(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_pipe-flow-in.csv'],double([pipe_cols;out.ss_terminal.qqinopt(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_pipe-flow-out.csv'],double([pipe_cols;out.ss_terminal.qqoutopt(1,pipe_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_comp-flow-in.csv'],double([1:out.ss_terminal.n0.nc;out.ss_terminal.qqinopt(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_comp-flow-out.csv'],double([1:out.ss_terminal.n0.nc;out.ss_terminal.qqoutopt(1,comp_cols)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_nodal-press_terminalure.csv'],double([[1:out.ss_terminal.n0.nv];out.ss_terminal.ppoptnodal(1,:)]),'precision',16,'delimiter',',');
        %dlmwrite([mfolder '\output_ss_terminal_gnode-physical-withdrawals.csv'],double([out.ss_terminal.gunique';out.ss_terminal.dgflows]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_nonslack-flows.csv'],double([out.ss_terminal.fn';out.ss_terminal.flows_all(ss_terminal.n0.nonslack_nodes,1)']),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_terminal_gnode-supply-flows.csv'],double([1:GN0;out.ss_terminal.gss_terminalol_all(1,:)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_terminal_gnode-demand-flows.csv'],double([1:GN0;out.ss_terminal.gdsol_all(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_slack-flows.csv'],double([out.ss_terminal.pn';out.ss_terminal.supp_flow(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_comp-ratios.csv'],double([[1:out.ss_terminal.cn];out.ss_terminal.ccopt(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_comp-discharge-press_terminalure.csv'],double([[1:out.ss_terminal.cn];out.ss_terminal.csetopt(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_comp-power.csv'],double([[1:out.ss_terminal.cn];out.ss_terminal.cpowopt(1,:)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_terminal_lmp-nodal-all.csv'],double([out.ss_terminal.fn';out.ss_terminal.trlmpnodal(1,:)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_terminal_lmp-gnodes.csv'],double([[1:GN0];out.ss_terminal.gnodelmp(1,:)]),'precision',16,'delimiter',',');
        %dlmwrite([mfolder '\output_ss_terminal_lmp-bidders.csv'],double([out.ss_terminal.gunique';out.ss_terminal.dglmp]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_terminal_pipe-lmp-in.csv'],double([pipe_cols;out.ss_terminal.lmpin(1,pipe_cols)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_terminal_pipe-lmp-out.csv'],double([pipe_cols;out.ss_terminal.lmpout(1,pipe_cols)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_terminal_comp-lmp-in.csv'],double([1:out.ss_terminal.n0.nc;out.ss_terminal.lmpin(1,comp_cols)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_terminal_comp-lmp-out.csv'],double([1:out.ss_terminal.n0.nc;out.ss_terminal.lmpout(1,comp_cols)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_terminal_comp-pmax-mp.csv'],double([[1:out.ss_terminal.n0.nc];out.ss_terminal.mult0_pmax(1,:)]),'precision',16,'delimiter',',');
%         dlmwrite([mfolder '\output_ss_terminal_comp-hpmax-mp.csv'],double([[1:out.ss_terminal.n0.nc];out.ss_terminal.mult0_cmax(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_flowbalrel.csv'],double([[1:out.ss_terminal.n0.nv];out.ss_terminal.flowbalrel(1,:)]),'precision',16,'delimiter',',');
        dlmwrite([mfolder '\output_ss_terminal_pipe-mass.csv'],double([pipe_cols;out.ss_terminal.pipe_mass_0(1,pipe_cols)]),'precision',16,'delimiter',',');
end
    if(par.out.intervals_out>0)
        %inputs on intervals
        out.ss_terminal.int_qbar=ss_terminal.int_qbar; out.ss_terminal.int_dmax=ss_terminal.int_dmax; out.ss_terminal.int_dmin=ss_terminal.int_dmin;
        out.ss_terminal.int_smax=ss_terminal.int_smax; out.ss_terminal.int_smin=ss_terminal.int_smin; out.ss_terminal.int_cd=ss_terminal.int_cd;
        out.ss_terminal.int_cs=ss_terminal.int_cs; out.ss_terminal.int_cslack=ss_terminal.int_cslack; out.ss_terminal.int_pslack=ss_terminal.int_pslack;
        if(ss_terminal.intervals>0 && ss_terminal.units==1)
            out.ss_terminal.int_cd=out.ss_terminal.int_cd*mmscfd_to_kgps;
            out.ss_terminal.int_cs=out.ss_terminal.int_cs*mmscfd_to_kgps; 
            out.ss_terminal.int_cslack=out.ss_terminal.int_cslack*mmscfd_to_kgps;
        end
        %------------------
        %revert to full gnode index list
        out.ss_terminal.gdsol_all=zeros(2,ss_terminal.n0.ng);
        out.ss_terminal.gdsol_all(:,ss_terminal.dmax_pos)=out.ss_terminal.gdsol;
        out.ss_terminal.gss_terminalol_all=zeros(2,ss_terminal.n0.ng);
        out.ss_terminal.gss_terminalol_all(:,ss_terminal.smax_pos)=out.ss_terminal.gss_terminalol;
        
        In=par.out.intervals_out;  %24 intervals on optimization period (1 day)
        T=ss_terminal.c.T/3600;             %time in hours
        int_bounds=[0:T/In:T]; out.ss_terminal.int_bounds=int_bounds;
        [out.ss_terminal.dbase_int]=pts_to_int(ss_terminal.m.xd'/3600,out.ss_terminal.dbase,int_bounds');
        [out.ss_terminal.gsub_int]=pts_to_int(ss_terminal.m.xd'/3600,out.ss_terminal.gsub,int_bounds');
        [out.ss_terminal.gslb_int]=pts_to_int(ss_terminal.m.xd'/3600,out.ss_terminal.gslb,int_bounds');
        [out.ss_terminal.gdub_int]=pts_to_int(ss_terminal.m.xd'/3600,out.ss_terminal.gdub,int_bounds');
        [out.ss_terminal.gdlb_int]=pts_to_int(ss_terminal.m.xd'/3600,out.ss_terminal.gdlb,int_bounds');
        [out.ss_terminal.gdsol_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.gdsol_all,int_bounds');
        [out.ss_terminal.gss_terminalol_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.gss_terminalol_all,int_bounds');
        [out.ss_terminal.dgflows_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.dgflows_all,int_bounds');
        [out.ss_terminal.supp_flow_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.supp_flow,int_bounds');
        [out.ss_terminal.flows_all_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.flows_all',int_bounds');
        [out.ss_terminal.Prslack_int]=pts_to_int(ss_terminal.m.xd'/3600,ss_terminal.m.Prslack',int_bounds');
        [out.ss_terminal.Prs_int]=pts_to_int(ss_terminal.m.xd'/3600,ss_terminal.m.Prs',int_bounds');
        [out.ss_terminal.Prd_int]=pts_to_int(ss_terminal.m.xd'/3600,ss_terminal.m.Prd',int_bounds');
        %------------------
        [out.ss_terminal.ppinopt_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.ppinopt,int_bounds');
        [out.ss_terminal.ppoutopt_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.ppoutopt,int_bounds');
        [out.ss_terminal.qqinopt_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.qqinopt,int_bounds');
        [out.ss_terminal.qqoutopt_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.qqoutopt,int_bounds');
        [out.ss_terminal.ppoptnodal_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.ppoptnodal,int_bounds');
        [out.ss_terminal.dgflows_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.dgflows,int_bounds');
        [out.ss_terminal.supp_flow_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.supp_flow,int_bounds');
        [out.ss_terminal.ccopt_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.ccopt,int_bounds');
        [out.ss_terminal.csetopt_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.csetopt,int_bounds');
        [out.ss_terminal.cpowopt_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.cpowopt,int_bounds');
        [out.ss_terminal.trlmpnodal_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.trlmpnodal,int_bounds');
        [out.ss_terminal.dglmp_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.dglmp,int_bounds');
        [out.ss_terminal.lmpin_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.lmpin,int_bounds');
        [out.ss_terminal.lmpout_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.lmpout,int_bounds');
        [out.ss_terminal.mult0_pmax_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.mult0_pmax,int_bounds');
        [out.ss_terminal.mult0_cmax_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.mult0_cmax,int_bounds');
        [out.ss_terminal.flowbalrel_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.flowbalrel,int_bounds');
        [out.ss_terminal.pipe_mass_int]=pts_to_int(out.ss_terminal.tt0,out.ss_terminal.pipe_mass_0,int_bounds');
        pipe_cols=[1:out.ss_terminal.n0.ne-out.ss_terminal.n0.nc]; comp_cols=[out.ss_terminal.n0.ne-out.ss_terminal.n0.nc+1:out.ss_terminal.n0.ne];
        

        if(par.out.steadystateonly==0)
            %write files
            dlmwrite([mfolder '\output_int_pipe-press_terminalure-in.csv'],double([pipe_cols;out.ss_terminal.ppinopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-press_terminalure-out.ss_terminal.csv'],double([pipe_cols;out.ss_terminal.ppoutopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-press_terminalure-in.csv'],double([1:out.ss_terminal.n0.nc;out.ss_terminal.ppinopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-press_terminalure-out.ss_terminal.csv'],double([1:out.ss_terminal.n0.nc;out.ss_terminal.ppoutopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-flow-in.csv'],double([pipe_cols;out.ss_terminal.qqinopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-flow-out.ss_terminal.csv'],double([pipe_cols;out.ss_terminal.qqoutopt_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-flow-in.csv'],double([1:out.ss_terminal.n0.nc;out.ss_terminal.qqinopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-flow-out.ss_terminal.csv'],double([1:out.ss_terminal.n0.nc;out.ss_terminal.qqoutopt_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_nodal-press_terminalure.csv'],double([[1:out.ss_terminal.n0.nv];out.ss_terminal.ppoptnodal_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_gnode-physical-withdrawals.csv'],double([ss_terminal.m.guniqueind';out.ss_terminal.dgflows_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_gnode-supply-flows.csv'],double([ss_terminal.n0.phys_node';out.ss_terminal.gss_terminalol_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_gnode-demand-flows.csv'],double([ss_terminal.n0.phys_node';out.ss_terminal.gdsol_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_slack-flows.csv'],double([out.ss_terminal.pn';out.ss_terminal.supp_flow_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-ratios.csv'],double([[1:out.ss_terminal.cn];out.ss_terminal.ccopt_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-discharge-press_terminalure.csv'],double([[1:out.ss_terminal.cn];out.ss_terminal.csetopt_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-power.csv'],double([[1:out.ss_terminal.cn];out.ss_terminal.cpowopt_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_lmp-nodal-all.csv'],double([out.ss_terminal.fn';out.ss_terminal.trlmpnodal_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_lmp-bidders.csv'],douint_flow_constble([out.ss_terminal.gunique';out.ss_terminal.dglmp_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-lmp-in.csv'],double([pipe_cols;out.ss_terminal.lmpin_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_pipe-lmp-out.ss_terminal.csv'],double([pipe_cols;out.ss_terminal.lmpout_int(:,pipe_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-lmp-in.csv'],double([1:out.ss_terminal.n0.nc;out.ss_terminal.lmpin_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-lmp-out.ss_terminal.csv'],double([1:out.ss_terminal.n0.nc;out.ss_terminal.lmpout_int(:,comp_cols)]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-pmax-mp.csv'],double([[1:out.ss_terminal.n0.nc];out.ss_terminal.mult0_pmax_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_comp-hpmax-mp.csv'],double([[1:out.ss_terminal.n0.nc];out.ss_terminal.mult0_cmax_int]),'precision',16,'delimiter',',');
            dlmwrite([mfolder '\output_int_flowbalrel.csv'],double([[1:out.ss_terminal.n0.nv];out.ss_terminal.flowbalrel_int]),'precision',16,'delimiter',',');
        end
    end

if(ss_terminal.m.save_state==1)
dlmwrite([mfolder '\output_ss_terminal_state_save.csv'],double(full(ss_terminal.state_save)),'precision',16,'delimiter',',');
end

par.out.ss_terminal=out.ss_terminal;
%par.ss_terminal=ss_terminal;

%out.ss_terminal.mult0_pmax=ss_terminal.mult0_pmax/2*ss_terminal.m.N/(ss_terminal.c.psc/1000000)/mpa_to_psi;    %output press_terminalure marginal prices ($/
%out.ss_terminal.mult0_cmax=ss_terminal.mult0_cmax/2*ss_terminal.m.N*3.6/0.75; %compress_terminalion marginal prices ($/hp)

function [xints]=pts_to_int(tpts,xpts,ibnds)
    xbnds=interp1qr(tpts,xpts,ibnds); In=length(ibnds)-1;
    xints=(xbnds(1:In,:)+xbnds(2:In+1,:))/2;
return;