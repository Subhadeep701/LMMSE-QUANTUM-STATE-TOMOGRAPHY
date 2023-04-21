clc;
clear all;
close all;
Lin=70;
B={};n=6;sample=1000;clnumber=10;
N=sample/clnumber;
stat=rand(1,Lin);stat=stat/sum(stat);


%%%%%%%%%%%%%%%%%%%%%%target state%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i2=1:Lin
    X = (randn(n) + 1i*randn(n))/sqrt(2);
    [Q,R] = qr(X);
    %v=[ones(1,1),zeros(1,n-1)];v=v/sum(v);
    %v=ones(1,n);v=v/sum(v);
    v=rand(1,n); v=v/sum(v);
    %v=[ones(1,3),zeros(1,n-3)];v=v/sum(v);
    D1{1,i2}=Q'*(diag(v))*Q;
    %D1{1,i2}=Q'*(diag(v(i2,:)))*Q;
end
%  llast=measure(stat);
%            Dlast=D1{1,llast};

noise_cov_mtx=zeros(n^2,n^2);
PT_LS=zeros(1,n^2);
%%%%%%%%%%%%%%%%%%%%%% Generalized Pauli and SIC-POVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sigma_operator,tetrahedron_2qb]=sigpov(n);

%%%%%%%%%%%%%%%%H matrix for LINEAR INVERSION%%%%%%%%%%%%%%%%%%%%%%
for d=1:n^2
    for c=1:(n^2)-1
        Anew(d,c)=sum(diag(tetrahedron_2qb(:,:,d)*sigma_operator(:,:,c)))*sqrt(n/(2*(n-1)));
    end
end
W_LS=pinv(Anew);
%U=[0.01 0.1 0.5];

E_thet_cov=zeros(n^2-1,n^2-1);noise_ls_cov=zeros(n^2,n^2);
E_thet=zeros(n^2-1,1);E_Y=zeros(n^2,1);E_thet_noise1=zeros(n^2-1,n^2);E_noise_thet1=zeros(n^2,n^2-1);
R_xy=zeros(n^2-1,n^2);R_yy=zeros(n^2,n^2);

for i=1:clnumber
    PT=zeros(1,n^2);
    C=zeros(1,n^2);
    l=measure(stat);
    D=D1{1,l};
    STORE_D{1,i}=D;
    WD=D;
    probs_2qbe = qmt(WD, tetrahedron_2qb);
    counts_2qbe = histc(rand(N,1), [0; cumsum(probs_2qbe)]);
    counts_2qbe = counts_2qbe(1:end-1);
    PT=counts_2qbe/N;

    PTN{1,i}=PT;
    Y(:,i)= (PT-1/n)*(n/(n-1));
    thet{1,i}=inv(Anew'*Anew)*Anew'*Y(:,i);
    E_Y=Y(:,i)+E_Y;
    E_thet=thet{1,i}+E_thet;
end
E_Y=E_Y/clnumber;
E_thet=E_thet/clnumber;
for i11=1:clnumber
    E_thet_cov=(thet{1,i11}-E_thet)*(thet{1,i11}-E_thet)'+E_thet_cov;
end
E_thet_cov=E_thet_cov/clnumber;
noise_ls_vec=Y-E_Y;
for n1=1:clnumber
    noise_ls_cov1=noise_ls_vec(:,n1)*noise_ls_vec(:,n1)';
    noise_ls_cov=noise_ls_cov1+noise_ls_cov;
    E_thet_noise1=(thet{n1}-E_thet)*noise_ls_vec(:,n1)'+ E_thet_noise1;
    E_noise_thet1=noise_ls_vec(:,n1)*(thet{n1}-E_thet)'+E_noise_thet1;
    R_xy=(thet{n1}-E_thet)*(Y(:,n1)-E_Y)'+R_xy;
    R_yy=(Y(:,n1)-E_Y)*(Y(:,n1)-E_Y)'+R_yy;
end
noise_ls_cov=noise_ls_cov/clnumber;
E_thet_noise1=E_thet_noise1/clnumber;
E_noise_thet1=E_noise_thet1/clnumber;
R_xy=R_xy/clnumber;
R_yy=R_yy/clnumber;
R_yy=R_yy+0.00001*eye(n^2);
W_Lmmse=R_xy*inv(R_yy);
W_Lmmse_init=W_Lmmse;
%%%%%Testing%%%%%
%     for j2=1:Lin
%      for i=1:n^2-1
%      thet_check(i,j2)=sum(diag(D1{1,j2}*sigma_operator(:,:,i)))*sqrt(n/(2*(n-1)));
%      end
%     end
%     E_thet_check=0.1*thet_check(:,1)+0.3*thet_check(:,2)+0.6*thet_check(:,3);
%     thet_check_diff=thet_check-E_thet_check;
%     E_thetcheck_cov=zeros(n^2-1,n^2-1);
%     for j3=1:Lin
%         E_thetcheck_cov=stat(j3)*thet_check_diff(:,j3)*thet_check_diff(:,j3)'+E_thetcheck_cov;
%     end
%     Dln_LS1= eye(n,n)/n;
%     for e=1:n^2-1
%         Dln_LS1=thet_check(e,3)*sigma_operator(:,:,e)*sqrt((n-1)/(2*n))+ Dln_LS1;
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%Adaptive part after initial LMMSE%%%%%%
% R_xy=zeros(n^2-1,n^2);R_yy=zeros(n^2,n^2);
% E_thet=zeros(n^2-1,1);E_Y=zeros(n^2,1);
for mas_i=1:50
   
    Block=50;
    cl_block=300;
    T_ls=zeros(1,Block);
    T_LMMSE=zeros(1,Block);
    T_LMMSE_init=zeros(1,Block);
    for bl_in= 1:Block
        thet_cluster=zeros(n^2-1,cl_block);
        Y_cluster=zeros(n^2,cl_block);
        Block_R_xy=zeros(n^2-1,n^2);
        Block_R_yy=zeros(n^2,n^2);
        %%%% cluster inside block %%%%%
        for cl_in=1:cl_block
            samplecluster=n^2*10;
            lcluster=measure(stat);
            Dcluster=D1{1,lcluster};
            probs_cluster = qmt(Dcluster, tetrahedron_2qb);
            counts_cluster = histc(rand(samplecluster,1), [0; cumsum(probs_cluster)]);
            counts_cluster = counts_cluster(1:end-1);
            PTcluster=counts_cluster/samplecluster;
            Y_cluster(:,cl_in)=(PTcluster-1/n)*(n/(n-1));
            %%%Stroring theta after each cluster%%%
            %thet_cluster(:,cl_in)=W_Lmmse*Y_cluster(:,cl_in);
            thet_cluster(:,cl_in)=W_LS*Y_cluster(:,cl_in);
            %Lmmse_thet_cluster(:,cl_in)=W_Lmmse*(Y_cluster(:,cl_in)-E_Y)+E_thet;

        end
        %%%%Adaptive mean after one cluster %%%%%%
        block_thet_adaptive=sum(thet_cluster,2)/cl_block;
        %%%%%deduction of mean for generation of R_thetY
        cov_vec_cl=thet_cluster-block_thet_adaptive;
        %%%%% adaptive mean Y %%%
        block_Y_adaptive=sum(Y_cluster,2)/cl_block;
        %%%% deduction of mean for generation of covariance vector R_yy and R_thety%%%%
        Y_cov_cl=Y_cluster-block_Y_adaptive;
        %%%%Generation of R_thety and R_yy itertively for one block%%%%
        for i=1:cl_in
            Block_R_xy=cov_vec_cl(:,cl_in)*Y_cov_cl(:,cl_in)'+Block_R_xy;
            Block_R_yy= Y_cov_cl(:,cl_in)*Y_cov_cl(:,cl_in)'+Block_R_yy;
        end
        Block_R_yy=Block_R_yy/cl_in;
        Block_R_xy=Block_R_xy/cl_in;
        temp_ratio=sample/(sample+samplecluster*cl_block);
        sample=(sample+samplecluster*cl_block);
        %temp_ratio=0.9;
        E_Y=temp_ratio*E_Y+(1-temp_ratio)*block_Y_adaptive;
        E_thet=temp_ratio*E_thet+(1-temp_ratio)*block_thet_adaptive;
        R_xy;
        R_yy;
        R_xy=temp_ratio*R_xy+(1-temp_ratio)*Block_R_xy;
        R_yy=temp_ratio*R_yy+(1-temp_ratio)*Block_R_yy;

        W_LMMSE_adapt=(W_Lmmse'+0.5*(R_xy'-R_yy*W_Lmmse'))';
        W_Lmmse=R_xy*inv(R_yy+0.0001*eye(n^2))





        %%%%%New state estimation%%%%%%%
        %repe=1:20
        samplenew=360;
        llast=measure(stat);
        Dlast=D1{1,llast};

        probs_2qb = qmt(Dlast, tetrahedron_2qb);


        counts_2qb = histc(rand(samplenew,1), [0; cumsum(probs_2qb)]);
        counts_2qb = counts_2qb(1:end-1);
        PTnew=counts_2qb/samplenew;
        Y_new=(PTnew-1/n)*(n/(n-1));
        thetnew=pinv(Anew)*Y_new;


        Lmmse_thet=W_Lmmse*(Y_new-E_Y)+E_thet;
        Lmmse_thet_init=W_LMMSE_adapt*(Y_new-E_Y)+E_thet;
        Dln= eye(n,n)/n;
        Dln_init=eye(n,n)/n;
        Dln_LS= eye(n,n)/n;
        for e=1:n^2-1
            Dln=Lmmse_thet(e)*sigma_operator(:,:,e)*sqrt((n-1)/(2*n))+Dln;
            Dln_init=Lmmse_thet_init(e)*sigma_operator(:,:,e)*sqrt((n-1)/(2*n))+Dln_init;
            Dln_LS=thetnew(e)*sigma_operator(:,:,e)*sqrt((n-1)/(2*n))+ Dln_LS;
        end
        Dln=proj_spectrahedron(Dln);
        Dln_init=proj_spectrahedron(Dln_init);
        %Dln_LS=proj_spectrahedron(Dln_LS);
        T_ls(bl_in)=sum(abs(eig(Dlast-Dln_LS)))/2;
        T_LMMSE(bl_in)=sum(abs(eig(Dlast-Dln)))/2;
        T_LMMSE_init(bl_in)=sum(abs(eig(Dlast-Dln_init)))/2;
    end
    MT_ls(mas_i,:)=T_ls;
    MT_lmmse(mas_i)=T_LMMSE;
    MT_Lmmse(mas_i)=T_LMMSE_init;
end
plot([1:3:Block],T_ls(1:3:end),'-ok','LineWidth',1.0,'MarkerSize',10);
hold on;
plot(1:3:Block,T_LMMSE(1:3:end),'-*k','LineWidth',1.0,'MarkerSize',10);
plot(1:3:Block,T_LMMSE_init(1:3:end),'--','LineWidth',1.0,'MarkerSize',10);
hold off;
xlabel('iteration blocks','fontweight','bold','fontsize',10);
ylabel('Sample average trace distance','fontweight','bold','fontsize',10);
legend('LS based estimate','LMMSE based estimate','adaptive LMMSE');