function [c,cecon_o,ceff_o]=pipe_obj_base_static(x,par)

%index parameters
N1=par.N1; N=par.N; M=par.M; C=par.C; FN=par.FN; NE=par.NE; GN=par.GN;
Ts=par.Ts; w=par.w; D=-par.D*2/Ts; 
clinks=par.comp_pos(:,2); slinks=par.comp_pos(par.spos,2);
%problem parameters
m=par.mpow; xs=par.xs; ew=par.econweight;
%variables
if(par.use_init_state==0)
    q=reshape(x(N1*FN+1:N1*M),NE,N1);
    comps=reshape(x(N1*M+1:N1*M+N1*C),C,N1);
elseif(par.use_init_state==1)
    st_ind=FN;
    %pp_init=par.state_init(st_ind+1:st_ind+FN); st_ind=st_ind+FN;
    qq_init=par.state_init(st_ind+1:st_ind+NE); st_ind=st_ind+NE;
    cc_init=par.state_init(st_ind+1:st_ind+C); st_ind=st_ind+C;
    %p=[pp_init reshape(x(1:N*FN),FN,N)]; 
    q=[qq_init reshape(x(N*FN+1:N*M),NE,N)];
    comps=[cc_init reshape(x(N*M+1:N*M+N*C),C,N)];
end
qcomp=q(clinks,:); fs=q(slinks,:);


%efficiency objective
ceff=sum((diag(xs(clinks))*abs(qcomp)).*((comps).^(m)-1),1)*w;
%economic objective
cecon=sum((diag(xs(slinks))*par.prslack).*fs,1)*w;

%output objective
cecon_o=(Ts/2)*cecon/par.objsc;
ceff_o=(Ts/2)*ceff/par.objsc;
c=(Ts/2)*(ceff*(1-ew)+cecon*ew)/par.objsc;